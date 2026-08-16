rm(list=ls())
library(quadprog)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(glmnet)
library(flare)
library(QUIC)
library(huge)
library(matrixcalc)
library(survival)
library(MASS)
library(lpSolve)
library(readxl)
library(dplyr)
library(ggplot2)
library(qqman)

# ---------------------------------------------------------------------------
# Input selection. No paths are hard-coded: each file is chosen at run time.
# In an interactive session a file dialog is shown; otherwise the path is read
# from the console. Set the environment variables listed below to run
# unattended (for example inside a compute capsule).
# ---------------------------------------------------------------------------
pick_file <- function(msg, envvar = NA_character_) {
  if (!is.na(envvar)) {
    v <- Sys.getenv(envvar, unset = "")
    if (nzchar(v) && file.exists(v)) { cat("Using", envvar, "=", v, "\n"); return(v) }
  }
  cat("\n>>> ", msg, "\n", sep = "")
  flush.console()
  path <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() && rstudioapi::hasFun("selectFile")) {
    path <- tryCatch(rstudioapi::selectFile(caption = msg), error = function(e) NULL)
  }
  if ((is.null(path) || !length(path)) && interactive()) {
    path <- tryCatch(file.choose(), error = function(e) NULL)
  }
  while (is.null(path) || !length(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    path <- trimws(readline("    Full path to the file: "))
    path <- gsub('^"|"$', "", gsub("\\\\", "/", path))
    if (!nzchar(path) || !file.exists(path)) { cat("    Not found.\n"); path <- NULL }
  }
  cat("    Selected: ", path, "\n", sep = "")
  path
}


sourceCpp(Sys.getenv("SIM_UNIVLIB", unset = "sim_univLib.cpp"))

#========================================
# Data Preparation
#========================================

df <- read.table(pick_file("Select the aligned PLINK .raw genotype file", "GENO_RAW"), header = TRUE)
summary(df)
df2 <- read.csv(pick_file("Select the cleaned clinical CSV", "CLIN_CSV"), header = TRUE, stringsAsFactors = FALSE)
summary(df2)
snp_list <- readLines("730_snp_list.txt")
head(snp_list)
pattern <- paste0("^", snp_list, collapse = "|")
keep_cols <- c("FID", "IID", grep(pattern, names(df), value = TRUE))
df_snp <- df[, keep_cols]
snp_cols <- setdiff(names(df_snp), c("FID", "IID"))
df_snp[snp_cols] <- lapply(df_snp[snp_cols], as.factor)
df_all <- merge(df_snp, df2, by = "IID")

df_all$SEX <- as.factor(df_all$SEX)
df_all$smoksort <- as.factor(df_all$smoksort)
df_all$early_late <- as.factor(df_all$early_late)
df_all$RADS <- as.factor(df_all$RADS)
df_all$chemotx <- as.factor(df_all$chemotx)
df_all$surgery <- as.factor(df_all$surgery)
df_all <- df_all %>%
  filter(!is.na(OS_month), !is.na(DEAD))

# Missing value imputation function
impute_mode_mean <- function(d) {
  d <- as.data.frame(d)
  
  char_cols <- vapply(d, is.character, logical(1))
  d[char_cols] <- lapply(d[char_cols], as.factor)
  
  fac_cols <- vapply(d, is.factor, logical(1))
  d[fac_cols] <- lapply(d[fac_cols], function(x) {
    tab <- table(x, useNA = "no")
    if (length(tab) == 0) return(x)
    mode_val <- names(tab)[which.max(tab)]
    if (!mode_val %in% levels(x)) x <- factor(x, levels = c(levels(x), mode_val))
    x[is.na(x)] <- mode_val
    x
  })
  
  num_cols <- vapply(d, is.numeric, logical(1))
  d[num_cols] <- lapply(d[num_cols], function(x) {
    m <- mean(x, na.rm = TRUE)
    if (is.nan(m)) return(x)
    x[is.na(x)] <- m
    x
  })
  
  d
}

df_imputed <- impute_mode_mean(df_all)

df_imputed[snp_cols] <- lapply(df_imputed[snp_cols], as.character)
df_imputed[snp_cols] <- lapply(df_imputed[snp_cols], as.numeric)
final_vars <- c(snp_cols, "SEX", "AGE", "smoksort", "pc1", "pc2", "pc3",
                "early_late", "RADS", "chemotx", "surgery")
formula_str <- paste("~", paste(final_vars, collapse = " + "))
X <- model.matrix(as.formula(formula_str), data = df_imputed)
X <- as.matrix(as.data.frame(X))[,-1]

# Save column names of X matrix
var_names <- colnames(X)

n_lambda <- 100
nfold <- 5
tol <- 1.0e-6
maxiter <- 50000

# [Fix 1] Dynamically get p
p <- ncol(X)
cat("Number of predictors (p):", p, "\n")
cat("Number of original variables (final_vars):", length(final_vars), "\n")
cat("Number of SNP columns (snp_cols):", length(snp_cols), "\n\n")

#========================================
# Lasso Estimation
#========================================

cvobj_glmnet <- cv.glmnet(x=X, y=cbind(time=df_imputed$OS_month, status=df_imputed$DEAD), family="cox", 
                          alpha=1, standardize=F, nfolds = nfold, nlambda=n_lambda)

beta_glmnet <- as.vector(coef(glmnet(x=X, y=cbind(time=df_imputed$OS_month, status=df_imputed$DEAD), 
                                     family="cox", alpha=1, lambda=cvobj_glmnet$lambda.min, standardize=F, 
                                     thresh=tol, maxit=maxiter)))

neg_loglik_glmnet <- 0
neg_dloglik_glmnet <- rep(0,p)
neg_ddloglik_glmnet <- matrix(0, nrow=p, ncol=p)
score_sq <- matrix(0, nrow=p, ncol=p)
neg_loglik_functions_cpp_ext(neg_loglik_glmnet, neg_dloglik_glmnet, neg_ddloglik_glmnet, score_sq,
                             X, df_imputed$OS_month, df_imputed$DEAD, beta_glmnet)
r <- eigen(score_sq)
r$values[r$values<=1.0e-14] <- 0

#========================================
# Cross-validation to select tuning parameter
#========================================

n_cv <- 5
n <- dim(df_imputed)[1]
n_multi <- 30
multiplier_seq <- exp(seq(from=log(0.001), to=log(5), length.out=n_multi))
alpha_cv <- 0.1

# [Fix 4] Set random seed
set.seed(12345)

all_cv_idx <- rep(1:n_cv, (n+1)/n_cv)
all_cv_idx <- sample(all_cv_idx, size=n, replace = T)
all_cvpl2 <- array(NA, length(multiplier_seq))

for(jj in 1:length(multiplier_seq)) {
  cat("Processing multiplier", jj, "of", length(multiplier_seq), "\n")
  multiplier <- multiplier_seq[jj]
  
  cv_idx <- NULL
  cvpl2 <- 0
  
  for(k in 1:n_cv) {
    cv_idx <- which(all_cv_idx==k)
    train_x <- X[-c(cv_idx),]
    test_x <- X[cv_idx,]
    train_time <- df_imputed$OS_month[-c(cv_idx)]
    test_time <- df_imputed$OS_month[cv_idx]
    train_delta <- df_imputed$DEAD[-c(cv_idx)]
    test_delta <- df_imputed$DEAD[cv_idx]
    
    cvobj_glmnet_train <- cv.glmnet(x=train_x, y=cbind(time=train_time, status=train_delta), family="cox",
                                    alpha=1, standardize=F, nfolds = nfold, nlambda=n_lambda)
    beta_glmnet_train <- as.vector(coef(glmnet(x=train_x, y=cbind(time=train_time, status=train_delta),
                                               family="cox", alpha=1, lambda=cvobj_glmnet_train$lambda.min, standardize=F,
                                               thresh=tol, maxit=maxiter)))
    neg_loglik_glmnet_train <- 0
    neg_dloglik_glmnet_train <- rep(0,p)
    neg_ddloglik_glmnet_train <- matrix(0, nrow=p, ncol=p)
    score_sq_train <- matrix(0, nrow=p, ncol=p)
    neg_loglik_functions_cpp_ext(neg_loglik_glmnet_train, neg_dloglik_glmnet_train, neg_ddloglik_glmnet_train,
                                 score_sq_train, train_x, train_time, train_delta, beta_glmnet_train)
    r_train <- eigen(score_sq_train)
    r_train$values[r_train$values<=1.0e-14] <- 0
    
    b_hat_new <- rep(NA, p)
    se_new <- rep(NA, p)
    mu <- multiplier*sqrt(log(p)/n)
    my_pos <- which(r_train$values > 0)
    my_rank <- sum(r_train$values > 0)
    Dmat <- diag(r_train$values[my_pos])
    dvec <- rep(0,my_rank)
    Amat <- t(rbind(-r_train$vectors[,my_pos]%*%Dmat, r_train$vectors[,my_pos]%*%Dmat))
    
    for(j in 1:p) {
      Vpos <- r_train$vectors[, my_pos, drop = FALSE]
      e_j <- rep(0, p); e_j[j] <- 1
      ej_proj <- Vpos %*% (t(Vpos) %*% e_j)
      feas_gap <- max(abs(e_j - ej_proj))
      mu_eff <- max(mu, feas_gap + 1e-8)
      bvec <- c(-e_j, e_j) - mu_eff
      res <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec)
      m <- as.vector(r_train$vectors[,my_pos]%*%res$solution) +
        as.vector(as.matrix(r_train$vectors[,-my_pos])%*%rep(0, p-my_rank))
      b_hat_new[j] <- beta_glmnet_train[j] - as.numeric(m%*%neg_dloglik_glmnet_train)
      
      # [Fix 3] Prevent NaN
      if (m[j] > 0) {
        se_new[j] <- sqrt(m[j]/nrow(train_x))
      } else {
        se_new[j] <- NA
      }
    }
    
    pval_new <- 2*pnorm(abs(b_hat_new/se_new), lower.tail = F)
    tmp_beta <- b_hat_new*as.numeric(pval_new < (alpha_cv/p))
    cvpl2 <- cvpl2 + loglik_cpp_ext(X=test_x, time=test_time, delta=test_delta, beta=tmp_beta)
  }
  
  all_cvpl2[jj] <- cvpl2
}

#========================================
# Final Estimation
#========================================

multiplier <- multiplier_seq[which.max(all_cvpl2)]

b_hat_new <- array(NA, p)
se_new <- array(NA, p)
theta_new <- matrix(NA, ncol=p, nrow=p)
mu_new <- multiplier*sqrt(log(p)/n)
my_pos <- which(r$values > 0)
my_rank <- sum(r$values > 0)
Dmat <- diag(r$values[my_pos])
dvec <- rep(0,my_rank)
Amat <- t(rbind(-r$vectors[,my_pos]%*%Dmat, r$vectors[,my_pos]%*%Dmat))

# [Fix 2] Add feasibility check
Vpos <- r$vectors[, my_pos, drop = FALSE]

cat("\nComputing final estimates...\n")
for(j in 1:p) {
  e_j <- rep(0, p); e_j[j] <- 1
  
  ej_proj <- Vpos %*% (t(Vpos) %*% e_j)
  feas_gap <- max(abs(e_j - ej_proj))
  mu_eff <- max(mu_new, feas_gap + 1e-8)
  
  bvec <- (c(-e_j, e_j) - mu_eff*rep(1,2*p))
  res <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec)
  m <- as.vector(r$vectors[,my_pos]%*%res$solution) + 
    as.vector(as.matrix(r$vectors[,-my_pos])%*%rep(0, p-my_rank))
  b_hat_new[j] <- beta_glmnet[j] - as.numeric(m%*%neg_dloglik_glmnet)
  
  # [Fix 3] Prevent NaN
  if (m[j] > 0) {
    se_new[j] <- sqrt(m[j]/n)
  } else {
    se_new[j] <- NA
  }
  
  theta_new[j,] <- m
}

pval_new <- 2*pnorm(abs(b_hat_new/se_new), lower.tail=F)

# Results Summary
cat("\n=== Results Summary ===\n")
cat("Number of predictors:", p, "\n")
cat("Sample size:", n, "\n")
cat("Chosen multiplier:", multiplier, "\n")
cat("Number of significant variables (Bonferroni, alpha=0.05):", sum(pval_new < 0.05/p, na.rm=TRUE), "\n")
cat("Number of NA standard errors:", sum(is.na(se_new)), "\n")

#========================================
# Output results CSV and Plots
#========================================

# 1. Organize result data frame
result <- data.frame(
  Variable = var_names,
  beta = b_hat_new,
  SE = se_new,
  pvalue = pval_new
)

stopifnot(nrow(result) == p)
cat("\nResult data frame dimensions:", nrow(result), "rows x", ncol(result), "cols\n")

# 2. Save results
write.csv(result, "cox_debiased_results_fixed.csv", row.names = FALSE)
cat("Results file saved: cox_debiased_results_fixed.csv\n")

#========================================
# [Fix] SNP variable identification - Multiple matching strategies
#========================================

# Convert snp_cols to format generated by model.matrix
snp_cols_safe <- make.names(snp_cols)

# Multiple matching: Original name | converted by make.names | IDs starting with rs
is_snp <- var_names %in% snp_cols |
  var_names %in% snp_cols_safe |
  grepl("^rs\\d+", var_names, ignore.case = TRUE)

cat("Number of SNP columns identified:", sum(is_snp), "\n")

# If the above methods fail, use positional matching (assuming SNP columns are at the beginning)
if (sum(is_snp) == 0) {
  cat("Warning: Unable to match SNPs by name, trying to use the first", length(snp_cols), "columns as SNPs\n")
  is_snp <- seq_len(p) <= length(snp_cols)
}

# Diagnostic information
cat("Examples of matched SNP variables (first 5):\n")
print(head(var_names[is_snp], 5))

#========================================
# Volcano Plot
#========================================

volcano_data <- result[is_snp, ]

if (nrow(volcano_data) > 0) {
  # Handle extreme values and NA
  volcano_data$pvalue[volcano_data$pvalue == 0] <- 1e-300
  volcano_data$pvalue[is.na(volcano_data$pvalue)] <- 1
  volcano_data$log10p <- -log10(volcano_data$pvalue)
  volcano_data <- volcano_data[!is.na(volcano_data$beta), ]
  
  n_snps <- nrow(volcano_data)
  bonferroni_threshold <- 0.05 / n_snps
  
  p_volcano <- ggplot(volcano_data, aes(x=beta, y=log10p)) +
    geom_point(alpha=0.6, color="steelblue") +
    geom_hline(yintercept = -log10(bonferroni_threshold), 
               linetype="dashed", color="red") +
    annotate("text", x = max(volcano_data$beta, na.rm=TRUE) * 0.7, 
             y = -log10(bonferroni_threshold) + 0.3,
             label = paste0("Bonferroni (p=", format(bonferroni_threshold, digits=2), ")"),
             color = "red", size = 3) +
    xlab("Effect size (beta)") +
    ylab("-log10(p-value)") +
    ggtitle(paste0("Volcano Plot (", n_snps, " SNPs)")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave("volcano_plot.png", p_volcano, width = 8, height = 6, dpi = 150)
  print(p_volcano)
  cat("Volcano plot saved: volcano_plot.png\n")
} else {
  cat("Warning: No SNP data found, skipping volcano plot.\n")
}

#========================================
# [Fix] Manhattan Plot - Support for external Map files
#========================================

if (nrow(volcano_data) > 0) {
  
  # Attempt to find .bim file (PLINK format, containing chromosome and position info)
  bim_files <- list.files(pattern = "\\.bim$", full.names = TRUE)
  
  if (length(bim_files) > 0) {
    # Plan A: Use .bim file
    cat("Found .bim file:", bim_files[1], "\n")
    bim <- read.table(bim_files[1], header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
    
    # Extract rsID from volcano_data$Variable (remove possible _Allele suffix)
    volcano_data$rsID <- gsub("_[ATCG]+$", "", volcano_data$Variable)
    
    manhattan_data <- merge(
      volcano_data[, c("rsID", "Variable", "pvalue")],
      bim[, c("SNP", "CHR", "BP")],
      by.x = "rsID",
      by.y = "SNP",
      all.x = TRUE
    )
    manhattan_data$SNP <- manhattan_data$Variable
    
    # Check merge results
    n_matched <- sum(!is.na(manhattan_data$CHR))
    cat("Number of SNPs successfully matched with chromosome info:", n_matched, "/", nrow(manhattan_data), "\n")
    
    if (n_matched == 0) {
      cat("Warning: SNP IDs in .bim file do not match data, trying other methods.\n")
      manhattan_data <- NULL
    }
    
  } else {
    manhattan_data <- NULL
  }
  
  # Plan B: If no .bim file or match failed, try extracting from file names
  if (is.null(manhattan_data)) {
    cat("No usable .bim file found, attempting to extract chromosome info from SNP names...\n")
    
    manhattan_data <- data.frame(
      SNP = volcano_data$Variable,
      P = volcano_data$pvalue,
      stringsAsFactors = FALSE
    )
    
    # Try multiple naming formats to extract chromosome number
    # Format 1: chr1:12345 or 1:12345
    chr_from_pos <- as.numeric(gsub("^(chr)?(\\d{1,2}):.*", "\\2", manhattan_data$SNP))
    # Format 2: 1_12345_A_G (PLINK2 format)
    chr_from_plink2 <- as.numeric(gsub("^(\\d{1,2})_\\d+_.*", "\\1", manhattan_data$SNP))
    
    # Merge results
    manhattan_data$CHR <- ifelse(!is.na(chr_from_pos) & chr_from_pos <= 26, chr_from_pos,
                                 ifelse(!is.na(chr_from_plink2) & chr_from_plink2 <= 26, chr_from_plink2, NA))
    manhattan_data$BP <- seq_len(nrow(manhattan_data))
  }
  
  # Plan C: If still no chromosome info
  if (all(is.na(manhattan_data$CHR))) {
    cat("\n========================================\n")
    cat("Note: Unable to extract chromosome info from SNP names.\n")
    cat("Common rsID formats (e.g., rs12345) do not contain chromosome numbers.\n")
    cat("To generate a meaningful Manhattan plot, please provide one of the following files:\n")
    cat("  1. PLINK .bim file (in the same directory as data)\n")
    cat("  2. Annotation file containing SNP, CHR, BP columns\n")
    cat("\nWill plot p-value distribution sorted sequentially (non-standard Manhattan plot).\n")
    cat("========================================\n\n")
    
    manhattan_data$CHR <- 1
    manhattan_data$BP <- seq_len(nrow(manhattan_data))
    plot_title <- "P-value Distribution (No CHR info available)"
  } else {
    # Handle partial missing values
    manhattan_data$CHR[is.na(manhattan_data$CHR)] <- 0
    if (any(manhattan_data$CHR == 0)) {
      cat("Note: Some SNPs (", sum(manhattan_data$CHR == 0), ") failed to extract chromosome info, set to CHR 0.\n")
    }
    plot_title <- "Manhattan Plot"
  }
  
  # Remove invalid data
  manhattan_data <- manhattan_data[!is.na(manhattan_data$P) & manhattan_data$P > 0, ]
  
  if (nrow(manhattan_data) > 0) {
    png("manhattan_plot.png", width = 1000, height = 600, res = 100)
    tryCatch({
      manhattan(manhattan_data,
                chr = "CHR", bp = "BP", p = "P", snp = "SNP",
                genomewideline = -log10(5e-8),
                suggestiveline = -log10(1e-5),
                main = plot_title,
                col = c("steelblue", "orange"),
                cex = 0.6,
                chrlabs = if(all(manhattan_data$CHR == 1)) NULL else NULL)
    }, error = function(e) {
      cat("Error generating Manhattan plot:", conditionMessage(e), "\n")
      plot(1, type="n", main="Manhattan Plot Error", xlab="", ylab="")
      text(1, 1, paste("Error:", conditionMessage(e)))
    })
    dev.off()
    cat("Manhattan plot saved: manhattan_plot.png\n")
  }
  
} else {
  cat("Warning: No valid data, skipping Manhattan plot.\n")
}

#========================================
# List of Significant SNPs
#========================================

n_snp_total <- sum(is_snp)
significant_snps <- result[is_snp & !is.na(result$pvalue) & result$pvalue < 0.05/n_snp_total, ]

if (nrow(significant_snps) > 0) {
  cat("\n=== Significant SNPs (Bonferroni corrected, alpha=0.05) ===\n")
  significant_snps <- significant_snps[order(significant_snps$pvalue), ]
  print(significant_snps)
  
  # Save significant SNPs to a separate file
  write.csv(significant_snps, "significant_snps.csv", row.names = FALSE)
  cat("Significant SNP list saved: significant_snps.csv\n")
} else {
  cat("\nNo significant SNPs (Bonferroni corrected, threshold p <", format(0.05/n_snp_total, digits=3), ")\n")
}

cat("\n=== Analysis Complete ===\n")