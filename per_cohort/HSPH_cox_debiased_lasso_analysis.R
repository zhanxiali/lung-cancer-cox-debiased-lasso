rm(list=ls())

#========================================
# HSPH 984 - Cox Debiased Lasso åˆ†æž
# åˆ†æž SNP ä¸Žæ€»ç”Ÿå­˜æœŸçš„å…³è” (ä¸åˆ†æžtype)
# åŸºäºŽ Oncoarray ç‰ˆæœ¬ä¿®æ”¹ no immu 
# åˆ›å»º rsID æ˜ å°„ (ä¿ç•™ Allele åŽç¼€)
#========================================

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
library(dplyr)
library(ggplot2)
library(ggrepel)

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



# åŠ è½½ C++ å‡½æ•°ï¼ˆéœ€è¦ sim_univLib.cpp æ–‡ä»¶ï¼‰
# è¯·ç¡®ä¿æ­¤æ–‡ä»¶åœ¨å·¥ä½œç›®å½•ä¸­
sourceCpp(Sys.getenv("SIM_UNIVLIB", unset = "../core/sim_univLib.cpp"))

#========================================
# è®¾ç½®å·¥ä½œç›®å½•
#========================================


cat("å·¥ä½œç›®å½•:", getwd(), "\n\n")

#========================================
# æ•°æ®å‡†å¤‡
#========================================

cat("=== ç¬¬ä¸€éƒ¨åˆ†ï¼šæ•°æ®å‡†å¤‡ ===\n\n")

# 1. è¯»å– SNP æ•°æ® (.raw æ ¼å¼)
df <- read.table(pick_file("Select the aligned PLINK .raw genotype file (HSPH)", "GENO_RAW"), header = TRUE, stringsAsFactors = FALSE)
cat("SNP æ•°æ®: ", nrow(df), "ä¾‹,", ncol(df), "åˆ—\n")

# 2. è¯»å–ä¸´åºŠæ•°æ®ï¼ˆå« PCAï¼‰
df2 <- read.csv(pick_file("Select the cleaned clinical CSV (HSPH)", "CLIN_CSV"), header = TRUE, stringsAsFactors = FALSE)
cat("ä¸´åºŠæ•°æ®: ", nrow(df2), "ä¾‹,", ncol(df2), "åˆ—\n")

# 3. è¯»å– SNP ä½ç½®æ–‡ä»¶ï¼ˆç”¨äºŽ rsID æ˜ å°„ï¼‰
snp_pos <- read.csv(pick_file("Select the candidate SNP position file", "SNP_POS"), header = TRUE, stringsAsFactors = FALSE)
cat("SNP ä½ç½®ä¿¡æ¯: ", nrow(snp_pos), "æ¡è®°å½•\n")

# åˆ›å»º CHR:POS -> rsID æ˜ å°„
snp_pos$match_key <- paste0(snp_pos$Chr, ":", snp_pos$Pos)
cat("SNP ä½ç½®ç¤ºä¾‹:\n")
print(head(snp_pos))

# 4. æå– SNP åˆ—åï¼ˆç¬¬7åˆ—å¼€å§‹æ˜¯SNPï¼Œå‰6åˆ—æ˜¯ FID, IID, PAT, MAT, SEX, PHENOTYPEï¼‰
snp_cols <- colnames(df)[7:ncol(df)]
cat("\nSNP æ•°é‡: ", length(snp_cols), "\n")
cat("SNP åˆ—åç¤ºä¾‹: ", paste(head(snp_cols, 3), collapse=", "), "\n\n")

# 5. å‡†å¤‡ SNP æ•°æ®æ¡†
df_snp <- df[, c("IID", snp_cols)]

# 6. åˆå¹¶æ•°æ®
df_all <- merge(df_snp, df2, by = "IID")
cat("åˆå¹¶åŽæ•°æ®: ", nrow(df_all), "ä¾‹\n\n")

#========================================
# å˜é‡å¤„ç†
#========================================

cat("=== ç¬¬äºŒéƒ¨åˆ†ï¼šå˜é‡å¤„ç† ===\n\n")

# å°†åˆ†ç±»å˜é‡è½¬ä¸ºå› å­
df_all$SEX <- as.factor(df_all$SEX)
df_all$smoksort <- as.factor(df_all$smoksort)
df_all$early_late <- as.factor(df_all$early_late)
df_all$RADS <- as.factor(df_all$RADS)
df_all$chemotx <- as.factor(df_all$chemotx)
df_all$surgery <- as.factor(df_all$surgery)

# è¿‡æ»¤ç¼ºå¤±ç”Ÿå­˜æ•°æ®çš„æ ·æœ¬
df_all <- df_all %>%
  filter(!is.na(OS_month), !is.na(DEAD))
cat("è¿‡æ»¤åŽæ ·æœ¬æ•°: ", nrow(df_all), "\n")

# ç¼ºå¤±å€¼å¡«è¡¥å‡½æ•°
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

# SNP åˆ—è½¬ä¸ºæ•°å€¼åž‹
df_imputed[snp_cols] <- lapply(df_imputed[snp_cols], as.character)
df_imputed[snp_cols] <- lapply(df_imputed[snp_cols], as.numeric)

#========================================
# æž„å»ºè®¾è®¡çŸ©é˜µ
#========================================

cat("\n=== ç¬¬ä¸‰éƒ¨åˆ†ï¼šæž„å»ºè®¾è®¡çŸ©é˜µ ===\n\n")

# å®šä¹‰æ¨¡åž‹å˜é‡ï¼ˆSNP + åå˜é‡ï¼‰
final_vars <- c(snp_cols, 
                "SEX", "AGE", "smoksort",      # äººå£å­¦
                "pc1", "pc2", "pc3",           # é—ä¼ ä¸»æˆåˆ†ï¼ˆå¿…é¡»ï¼‰
                "early_late",                   # ç–¾ç—…ç‰¹å¾
                "surgery", "chemotx", "RADS")  # æ²»ç–—

formula_str <- paste("~", paste(final_vars, collapse = " + "))
X <- model.matrix(as.formula(formula_str), data = df_imputed)
X <- as.matrix(as.data.frame(X))[,-1]  # åŽ»æŽ‰æˆªè·

# ä¿å­˜ X çŸ©é˜µçš„åˆ—å
var_names <- colnames(X)

# æ¨¡åž‹å‚æ•°
n_lambda <- 100
nfold <- 5
tol <- 1.0e-6
maxiter <- 50000

# åŠ¨æ€èŽ·å– p
p <- ncol(X)
n <- nrow(X)

cat("è®¾è®¡çŸ©é˜µç»´åº¦:\n")
cat("  æ ·æœ¬æ•° (n):", n, "\n")
cat("  å˜é‡æ•° (p):", p, "\n")
cat("  åŽŸå§‹å˜é‡æ•° (final_vars):", length(final_vars), "\n")
cat("  SNP åˆ—æ•° (snp_cols):", length(snp_cols), "\n\n")

#========================================
# Lasso ä¼°è®¡
#========================================

cat("=== ç¬¬å››éƒ¨åˆ†ï¼šLasso åˆå§‹ä¼°è®¡ ===\n\n")

cvobj_glmnet <- cv.glmnet(x=X, y=cbind(time=df_imputed$OS_month, status=df_imputed$DEAD), 
                          family="cox", alpha=1, standardize=F, nfolds = nfold, nlambda=n_lambda)

beta_glmnet <- as.vector(coef(glmnet(x=X, y=cbind(time=df_imputed$OS_month, status=df_imputed$DEAD), 
                                     family="cox", alpha=1, lambda=cvobj_glmnet$lambda.min, standardize=F, 
                                     thresh=tol, maxit=maxiter)))

cat("æœ€ä¼˜ lambda:", cvobj_glmnet$lambda.min, "\n")
cat("éžé›¶ç³»æ•°æ•°é‡:", sum(beta_glmnet != 0), "\n\n")

neg_loglik_glmnet <- 0
neg_dloglik_glmnet <- rep(0,p)
neg_ddloglik_glmnet <- matrix(0, nrow=p, ncol=p)
score_sq <- matrix(0, nrow=p, ncol=p)
neg_loglik_functions_cpp_ext(neg_loglik_glmnet, neg_dloglik_glmnet, neg_ddloglik_glmnet, score_sq,
                             X, df_imputed$OS_month, df_imputed$DEAD, beta_glmnet)
r <- eigen(score_sq)
r$values[r$values<=1.0e-14] <- 0

#========================================
# äº¤å‰éªŒè¯é€‰æ‹©è°ƒä¼˜å‚æ•°
#========================================

cat("=== ç¬¬äº”éƒ¨åˆ†ï¼šäº¤å‰éªŒè¯é€‰æ‹©è°ƒä¼˜å‚æ•° ===\n\n")

n_cv <- 5
n_multi <- 30
multiplier_seq <- exp(seq(from=log(0.001), to=log(5), length.out=n_multi))
alpha_cv <- 0.1

# è®¾ç½®éšæœºç§å­
set.seed(12345)

all_cv_idx <- rep(1:n_cv, (n+1)/n_cv)
all_cv_idx <- sample(all_cv_idx, size=n, replace = T)
all_cvpl2 <- array(NA, length(multiplier_seq))

for(jj in 1:length(multiplier_seq)) {
  cat("Processing multiplier", jj, "of", length(multiplier_seq), "\r")
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
      
      # é˜²æ­¢ NaN
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
# æœ€ç»ˆä¼°è®¡
#========================================

cat("\n\n=== ç¬¬å…­éƒ¨åˆ†ï¼šæœ€ç»ˆä¼°è®¡ ===\n\n")

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

# æ·»åŠ å¯è¡Œæ€§æ£€æŸ¥
Vpos <- r$vectors[, my_pos, drop = FALSE]

cat("Computing final estimates...\n")
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
  
  # é˜²æ­¢ NaN
  if (m[j] > 0) {
    se_new[j] <- sqrt(m[j]/n)
  } else {
    se_new[j] <- NA
  }
  
  theta_new[j,] <- m
}

pval_new <- 2*pnorm(abs(b_hat_new/se_new), lower.tail=F)

# ç»“æžœæ‘˜è¦
cat("\n=== Results Summary ===\n")
cat("Number of predictors:", p, "\n")
cat("Sample size:", n, "\n")
cat("Chosen multiplier:", multiplier, "\n")
cat("Number of significant variables (Bonferroni, alpha=0.05):", sum(pval_new < 0.05/p, na.rm=TRUE), "\n")
cat("Number of NA standard errors:", sum(is.na(se_new)), "\n")

#========================================
# æ•´ç†ç»“æžœ + rsID æ˜ å°„
#========================================

cat("\n=== ç¬¬ä¸ƒéƒ¨åˆ†ï¼šæ•´ç†ç»“æžœ ===\n\n")

result <- data.frame(
  Variable = var_names,
  beta = b_hat_new,
  SE = se_new,
  pvalue = pval_new,
  HR = exp(b_hat_new),
  HR_lower = exp(b_hat_new - 1.96 * se_new),
  HR_upper = exp(b_hat_new + 1.96 * se_new),
  stringsAsFactors = FALSE
)

stopifnot(nrow(result) == p)
cat("Result data frame dimensions:", nrow(result), "rows x", ncol(result), "cols\n")

#========================================
# SNP å˜é‡è¯†åˆ« + rsID è½¬æ¢
#========================================

# å°† snp_cols è½¬æ¢ä¸º model.matrix å¯èƒ½ç”Ÿæˆçš„æ ¼å¼
snp_cols_safe <- make.names(snp_cols)

# å¤šé‡åŒ¹é…ï¼šåŽŸå§‹å | make.namesè½¬æ¢åŽ | Xå¼€å¤´çš„æ ¼å¼
is_snp <- var_names %in% snp_cols |
  var_names %in% snp_cols_safe |
  grepl("^X\\d+\\.\\d+\\.", var_names)

cat("è¯†åˆ«çš„ SNP å˜é‡æ•°é‡:", sum(is_snp), "\n")

# å¦‚æžœä¸Šè¿°æ–¹æ³•éƒ½å¤±è´¥ï¼Œä½¿ç”¨ä½ç½®åŒ¹é…
if (sum(is_snp) == 0) {
  cat("è­¦å‘Šï¼šæ— æ³•é€šè¿‡åç§°åŒ¹é… SNPï¼Œå°è¯•ä½¿ç”¨å‰", length(snp_cols), "åˆ—ä½œä¸º SNP\n")
  is_snp <- seq_len(p) <= length(snp_cols)
}

#========================================
# ==============================================================================
# ä¿®æ”¹ç‰ˆï¼šåˆ›å»º rsID æ˜ å°„ï¼ˆä¿ç•™ç­‰ä½åŸºå› åŽç¼€ï¼‰
# ==============================================================================

cat("\n=== åˆ›å»º rsID æ˜ å°„ (ä¿ç•™ Allele åŽç¼€) ===\n")

# åˆå§‹åŒ–
result$CHR <- NA
result$BP <- NA
result$rsID <- NA

# éåŽ†æ¯ä¸€ä¸ªè¯†åˆ«å‡ºçš„ SNP
for (i in which(is_snp)) {
  var_name <- result$Variable[i]
  
  # 1. å°è¯•æå– CHR å’Œ POS (ç”¨äºŽæŸ¥è¡¨)
  # æ ¼å¼ç¤ºä¾‹: X9.129416317.C.T_T
  chr_match <- gsub("^X(\\d+)\\..*", "\\1", var_name)
  pos_match <- gsub("^X\\d+\\.(\\d+)\\..*", "\\1", var_name)
  
  # ç¡®ä¿æå–æˆåŠŸï¼ˆæ²¡æœ‰å–åˆ°æ•´ä¸ªå­—ç¬¦ä¸²ï¼‰
  if (chr_match != var_name && pos_match != var_name) {
    result$CHR[i] <- as.numeric(chr_match)
    result$BP[i] <- as.numeric(pos_match)
    
    # 2. æŸ¥æ‰¾å¯¹åº”çš„ rsID (ä»Žä½ç½®æ–‡ä»¶ä¸­)
    match_key <- paste0(chr_match, ":", pos_match)
    rs_idx <- which(snp_pos$match_key == match_key)
    
    if (length(rs_idx) > 0) {
      # æ‰¾åˆ°äº†çº¯ rsID (ä¾‹å¦‚ rs12345)
      pure_rsid <- snp_pos$SNP[rs_idx[1]]
      
      # 3. ã€æ ¸å¿ƒä¿®æ”¹ã€‘ä»ŽåŽŸå§‹å˜é‡åä¸­æå–æ•ˆåº”ç­‰ä½åŸºå› åŽç¼€
      # å‡è®¾ var_name ç»“å°¾æ˜¯ "_A", "_T", "_G", "_C"
      if (grepl("_[A-Za-z0-9]+$", var_name)) {
        # æå–æœ€åŽä¸€ä¸ªä¸‹åˆ’çº¿åŽé¢çš„å†…å®¹
        allele_suffix <- sub(".*_([A-Za-z0-9]+)$", "\\1", var_name)
        
        # æ‹¼æŽ¥æˆ rs12345_A çš„æ ¼å¼
        result$rsID[i] <- paste0(pure_rsid, "_", allele_suffix)
      } else {
        # å¦‚æžœåŽŸå§‹åå­—é‡Œæ²¡åŽç¼€ï¼Œå°±åªèƒ½ç”¨çº¯ rsID äº†
        result$rsID[i] <- pure_rsid
      }
    }
  }
}

# æ£€æŸ¥æ˜ å°„ç»“æžœ
cat("æˆåŠŸæ˜ å°„ rsID çš„ SNP æ•°é‡:", sum(!is.na(result$rsID)), "/", sum(is_snp), "\n")
if(sum(!is.na(result$rsID)) > 0) {
  cat("æ˜ å°„åŽçš„ rsID ç¤ºä¾‹ (å‰3ä¸ª):\n")
  print(head(result$rsID[!is.na(result$rsID)], 3))
}

# å¯¹äºŽ SNPï¼Œä½¿ç”¨ æ–°çš„rsID æ›¿æ¢ Variable
result$Variable_original <- result$Variable
result$Variable[!is.na(result$rsID)] <- result$rsID[!is.na(result$rsID)]

#========================================
# æ·»åŠ  FDR æ ¡æ­£
#========================================

# å¯¹ SNP ç»“æžœæ·»åŠ  FDR æ ¡æ­£
snp_pvalues <- result$pvalue[is_snp]
result$pvalue_FDR <- NA
result$pvalue_FDR[is_snp] <- p.adjust(snp_pvalues, method = "BH")

# ä¿å­˜æ‰€æœ‰ç»“æžœ
write.csv(result, "HSPH984_cox_debiased_results_all.csv", row.names = FALSE)
cat("å·²ä¿å­˜: HSPH984_cox_debiased_results_all.csv\n")

# 3. ä¿å­˜ SNP ç»“æžœ
snp_result <- result[is_snp, ]
snp_result <- snp_result[order(snp_result$pvalue), ]
write.csv(snp_result, "HSPH984_cox_debiased_results_SNPs.csv", row.names = FALSE)
cat("å·²ä¿å­˜: HSPH984_cox_debiased_results_SNPs.csv\n")

# 4. ä¿å­˜åå˜é‡ç»“æžœ
covar_result <- result[!is_snp, ]
covar_result <- covar_result[order(covar_result$pvalue), ]
write.csv(covar_result, "HSPH984_cox_debiased_results_covariates.csv", row.names = FALSE)
cat("å·²ä¿å­˜: HSPH984_cox_debiased_results_covariates.csv\n")

#========================================
# ç«å±±å›¾
#========================================

cat("\n=== ç¬¬å…«éƒ¨åˆ†ï¼šå¯è§†åŒ– ===\n\n")

volcano_data <- snp_result

if (nrow(volcano_data) > 0) {
  # å¤„ç†æžç«¯å€¼å’Œ NA
  volcano_data$pvalue[volcano_data$pvalue == 0] <- 1e-300
  volcano_data$pvalue[is.na(volcano_data$pvalue)] <- 1
  volcano_data$log10p <- -log10(volcano_data$pvalue)
  volcano_data <- volcano_data[!is.na(volcano_data$beta), ]
  
  n_snps <- nrow(volcano_data)
  bonferroni_threshold <- 0.05 / n_snps
  
  p_volcano <- ggplot(volcano_data, aes(x=beta, y=log10p)) +
    geom_point(alpha=0.6, color="steelblue", size=2) +
    geom_hline(yintercept = -log10(bonferroni_threshold), 
               linetype="dashed", color="red") +
    geom_hline(yintercept = -log10(0.05), 
               linetype="dotted", color="orange") +
    geom_text_repel(data = volcano_data[volcano_data$pvalue < 0.05, ],
                    aes(label = Variable), hjust = -0.1, vjust = 0, size = 3,
                    max.overlaps = 20) +
    annotate("text", x = max(volcano_data$beta, na.rm=TRUE) * 0.7, 
             y = -log10(bonferroni_threshold) + 0.3,
             label = paste0("Bonferroni (p=", format(bonferroni_threshold, digits=2), ")"),
             color = "red", size = 3) +
    annotate("text", x = max(volcano_data$beta, na.rm=TRUE) * 0.7, 
             y = -log10(0.05) + 0.3,
             label = "p=0.05",
             color = "orange", size = 3) +
    xlab("Effect size (beta)") +
    ylab("-log10(p-value)") +
    ggtitle(paste0("Volcano Plot: ", n_snps, " SNPs - HSPH 984")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave("HSPH984_volcano_plot.png", p_volcano, width = 10, height = 8, dpi = 200)
  print(p_volcano)
  cat("å·²ä¿å­˜ç«å±±å›¾: HSPH984_volcano_plot.png\n")
} else {
  cat("è­¦å‘Šï¼šæœªæ‰¾åˆ° SNP æ•°æ®ï¼Œè·³è¿‡ç«å±±å›¾ã€‚\n")
}

#========================================
# æ£®æž—å›¾ - Top SNPs (ä½¿ç”¨ rsID)
#========================================

if (nrow(snp_result) > 0) {
  top_snps <- head(snp_result, 20)
  top_snps$Variable <- factor(top_snps$Variable, levels = rev(top_snps$Variable))
  
  p_forest <- ggplot(top_snps, aes(x = HR, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
    geom_point(aes(color = pvalue < 0.05), size = 3) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                       labels = c("TRUE" = "P < 0.05", "FALSE" = "P â‰¥ 0.05"),
                       name = "Significance") +
    scale_x_log10() +
    labs(
      title = "Forest Plot: Top 20 SNPs by P-value - HSPH 984",
      subtitle = "Hazard Ratio with 95% CI",
      x = "Hazard Ratio (log scale)",
      y = ""
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave("HSPH984_forest_plot_SNPs.png", p_forest, width = 10, height = 8, dpi = 200)
  print(p_forest)
  cat("å·²ä¿å­˜æ£®æž—å›¾: HSPH984_forest_plot_SNPs.png\n")
}

#========================================
# åå˜é‡æ£®æž—å›¾
#========================================

if (nrow(covar_result) > 0) {
  covar_result$Variable <- factor(covar_result$Variable, levels = rev(covar_result$Variable))
  
  p_covar <- ggplot(covar_result, aes(x = HR, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
    geom_point(aes(color = pvalue < 0.05), size = 3) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                       labels = c("TRUE" = "P < 0.05", "FALSE" = "P â‰¥ 0.05"),
                       name = "Significance") +
    scale_x_log10() +
    labs(
      title = "Forest Plot: Clinical Covariates - HSPH 984",
      subtitle = "Hazard Ratio with 95% CI",
      x = "Hazard Ratio (log scale)",
      y = ""
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave("HSPH984_forest_plot_covariates.png", p_covar, width = 10, height = 6, dpi = 200)
  print(p_covar)
  cat("å·²ä¿å­˜åå˜é‡æ£®æž—å›¾: HSPH984_forest_plot_covariates.png\n")
}

#========================================
# æ›¼å“ˆé¡¿å›¾ï¼ˆggplot2 ç‰ˆæœ¬ - æ›´æ¸…æ™°çš„æ ‡ç­¾ï¼‰
#========================================

cat("\nç»‘åˆ¶æ›¼å“ˆé¡¿å›¾...\n")

manhattan_data <- snp_result[, c("Variable", "pvalue", "CHR", "BP")]
manhattan_data <- manhattan_data[!is.na(manhattan_data$CHR) & 
                                   !is.na(manhattan_data$pvalue) &
                                   manhattan_data$pvalue > 0, ]
manhattan_data$log10P <- -log10(manhattan_data$pvalue)

cat("æœ‰æ•ˆæ›¼å“ˆé¡¿å›¾æ•°æ®:", nrow(manhattan_data), "ä¸ª SNP\n")

if (nrow(manhattan_data) > 0) {
  n_snps <- nrow(manhattan_data)
  bonf_p <- 0.05 / n_snps
  
  p_manhattan <- ggplot(manhattan_data, aes(x = CHR, y = log10P)) +
    geom_point(aes(color = factor(CHR %% 2)), size = 3, alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange", linewidth = 0.8) +
    geom_hline(yintercept = -log10(bonf_p), linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_text_repel(data = manhattan_data[manhattan_data$pvalue < 0.05, ],
                    aes(label = Variable),
                    size = 4, fontface = "bold",
                    box.padding = 0.5,
                    max.overlaps = 20,
                    min.segment.length = 0) +
    scale_color_manual(values = c("steelblue", "darkorange"), guide = "none") +
    scale_x_continuous(breaks = sort(unique(manhattan_data$CHR))) +
    labs(x = "Chromosome", 
         y = expression(-log[10](P)), 
         title = "Manhattan Plot - HSPH 984") +
    annotate("text", x = max(manhattan_data$CHR) - 1, y = -log10(0.05) + 0.2, 
             label = "P = 0.05", color = "orange", size = 3.5, hjust = 0) +
    annotate("text", x = max(manhattan_data$CHR) - 1, y = -log10(bonf_p) + 0.2, 
             label = paste0("Bonferroni (P = ", format(bonf_p, digits = 2), ")"), 
             color = "red", size = 3.5, hjust = 0) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  ggsave("HSPH984_manhattan_plot.png", p_manhattan, width = 12, height = 6, dpi = 200)
  print(p_manhattan)
  cat("å·²ä¿å­˜æ›¼å“ˆé¡¿å›¾: HSPH984_manhattan_plot.png\n")
} else {
  cat("è­¦å‘Šï¼šæ— æœ‰æ•ˆæ•°æ®ï¼Œè·³è¿‡æ›¼å“ˆé¡¿å›¾ã€‚\n")
}

#========================================
# æ˜¾è‘—ç»“æžœæ±‡æ€»
#========================================

cat("\n========================================\n")
cat("=== æœ€ç»ˆç»“æžœæ±‡æ€» ===\n")
cat("========================================\n\n")

n_snp_total <- sum(is_snp)
bonf_threshold <- 0.05 / n_snp_total

# æ˜¾è‘— SNP (P < 0.05)
sig_snps_nominal <- snp_result[snp_result$pvalue < 0.05, ]
cat("P < 0.05 çš„ SNP æ•°é‡:", nrow(sig_snps_nominal), "\n")
if (nrow(sig_snps_nominal) > 0) {
  cat("\n--- P < 0.05 çš„ SNP ---\n")
  print(sig_snps_nominal[, c("Variable", "CHR", "BP", "beta", "HR", "HR_lower", "HR_upper", "pvalue", "pvalue_FDR")])
}

# æ˜¾è‘— SNP (Bonferroni)
sig_snps_bonf <- snp_result[snp_result$pvalue < bonf_threshold, ]
cat("\nBonferroni æ ¡æ­£åŽæ˜¾è‘—çš„ SNP æ•°é‡ (p <", format(bonf_threshold, digits=3), "):", nrow(sig_snps_bonf), "\n")
if (nrow(sig_snps_bonf) > 0) {
  cat("\n--- Bonferroni æ˜¾è‘—çš„ SNP ---\n")
  print(sig_snps_bonf[, c("Variable", "CHR", "BP", "beta", "HR", "HR_lower", "HR_upper", "pvalue", "pvalue_FDR")])
  
  write.csv(sig_snps_bonf, "HSPH984_significant_SNPs_bonferroni.csv", row.names = FALSE)
  cat("\nå·²ä¿å­˜: HSPH984_significant_SNPs_bonferroni.csv\n")
}

# æ˜¾è‘— SNP (FDR)
sig_snps_fdr <- snp_result[!is.na(snp_result$pvalue_FDR) & snp_result$pvalue_FDR < 0.05, ]
cat("\nFDR æ ¡æ­£åŽæ˜¾è‘—çš„ SNP æ•°é‡ (q < 0.05):", nrow(sig_snps_fdr), "\n")
if (nrow(sig_snps_fdr) > 0) {
  cat("\n--- FDR æ˜¾è‘—çš„ SNP ---\n")
  print(sig_snps_fdr[, c("Variable", "CHR", "BP", "beta", "HR", "pvalue", "pvalue_FDR")])
}

# æ˜¾è‘—åå˜é‡
sig_covar <- covar_result[covar_result$pvalue < 0.05, ]
cat("\næ˜¾è‘—åå˜é‡ (P < 0.05):\n")
if (nrow(sig_covar) > 0) {
  print(sig_covar[, c("Variable", "beta", "HR", "HR_lower", "HR_upper", "pvalue")])
} else {
  cat("æ— æ˜¾è‘—åå˜é‡\n")
}

#========================================
# å¤šé‡æ ¡æ­£æ±‡æ€»
#========================================

cat("\n========================================\n")
cat("=== å¤šé‡æ ¡æ­£æ±‡æ€» ===\n")
cat("========================================\n")
cat("SNP æ€»æ•°:", n_snp_total, "\n")
cat("P < 0.05 (æœªæ ¡æ­£):", nrow(sig_snps_nominal), "\n")
cat("Bonferroni æ˜¾è‘— (P <", format(bonf_threshold, digits=3), "):", nrow(sig_snps_bonf), "\n")
cat("FDR æ˜¾è‘— (q < 0.05):", nrow(sig_snps_fdr), "\n")

cat("\n========================================\n")
cat("=== åˆ†æžå®Œæˆ ===\n")
cat("========================================\n")
cat("\nè¾“å‡ºæ–‡ä»¶:\n")
cat("1. HSPH984_cox_debiased_results_all.csv - æ‰€æœ‰å˜é‡ç»“æžœ\n")
cat("2. HSPH984_cox_debiased_results_SNPs.csv - SNP ç»“æžœ (å« rsID å’Œ FDR)\n")
cat("3. HSPH984_cox_debiased_results_covariates.csv - åå˜é‡ç»“æžœ\n")
cat("4. HSPH984_volcano_plot.png - ç«å±±å›¾\n")
cat("5. HSPH984_forest_plot_SNPs.png - SNP æ£®æž—å›¾ (rsID æ ¼å¼)\n")
cat("6. HSPH984_forest_plot_covariates.png - åå˜é‡æ£®æž—å›¾\n")
cat("7. HSPH984_manhattan_plot.png - æ›¼å“ˆé¡¿å›¾ (rsID æ ‡æ³¨)\n")
