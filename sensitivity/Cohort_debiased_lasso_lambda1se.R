rm(list = ls())

# ==============================================================================
# Cox Debiased Lasso under the lambda.1se criterion  (Reviewer 2, Comment 2)
# v2 — robust interactive file selection
#
# HOW TO RUN
#   Open this script in RStudio (or R GUI) and click "Source", or run
#   source("path/to/this/script.R") in the console. Run it THREE times,
#   once per cohort.
#
#   File selection happens FIRST, before any package loading:
#     - A file dialog will pop up for each input. NOTE: on Windows the dialog
#       sometimes opens BEHIND the RStudio window - check the taskbar if
#       nothing seems to happen.
#     - If no dialog can be shown (e.g. non-interactive session), the script
#       falls back to asking you to PASTE the full file path in the console.
#
# WHAT THIS SCRIPT DOES
#   Repeats the primary per-cohort pipeline with the ONLY change being the
#   penalty criterion: lambda.1se instead of lambda.min, applied in BOTH the
#   main penalized fit and the inner cross-validation fits used to tune the
#   debiasing projection multiplier. Covariates, seed (12345), CV folds,
#   multiplier grid, and the debiasing algorithm are identical to the primary
#   analysis scripts.
#
# INPUTS (selected interactively)
#   1. The C++ library file (e.g. path/to/sim_univLib.cpp)
#   2. The cohort's aligned PLINK .raw genotype file
#   3. The cohort's cleaned clinical CSV (with pc1-pc3)
#   4. A cohort label typed at the console (used in output filenames)
#
# OUTPUTS (written to ./lambda1se_results/ under the current directory)
#   <label>_cox_debiased_results_all_1se.csv
#   <label>_1se_run_summary.csv
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Robust file chooser (dialog if possible, paste-path fallback otherwise)
# ------------------------------------------------------------------------------

norm_path <- function(p) {
  if (is.null(p) || length(p) == 0) return(NULL)
  p <- p[1]
  if (is.na(p) || !nzchar(p)) return(NULL)
  p
}

choose_file <- function(msg) {
  cat("\n>>> ", msg, "\n", sep = "")
  flush.console()
  path <- NULL

  # 1) RStudio's own file picker (works in RStudio Desktop AND RStudio Server,
  #    always renders inside RStudio so it cannot hide behind the window)
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() && rstudioapi::hasFun("selectFile")) {
    path <- norm_path(tryCatch(rstudioapi::selectFile(caption = msg),
                               error = function(e) NULL))
  }

  # 2) Native Windows dialog
  if (is.null(path) && .Platform$OS.type == "windows" && interactive()) {
    path <- norm_path(tryCatch(utils::choose.files(caption = msg, multi = FALSE),
                               error = function(e) NULL))
  }

  # 3) Base R dialog
  if (is.null(path) && interactive()) {
    path <- norm_path(tryCatch(file.choose(), error = function(e) NULL))
  }

  # 4) Paste-path fallback
  while (is.null(path) || !file.exists(path)) {
    cat("    (No file selected via dialog.)\n")
    path <- trimws(readline(prompt = "    Paste the FULL file path here and press Enter: "))
    path <- gsub('^"|"$', "", path)          # strip surrounding quotes if pasted
    path <- gsub("\\\\", "/", path)          # allow Windows backslash paths
    if (!nzchar(path) || !file.exists(path)) { cat("    File not found, please try again.\n"); path <- NULL }
  }
  cat("    Selected: ", path, "\n", sep = "")
  path
}

cat("==============================================================\n")
cat(" Lambda.1se rerun - please select the 3 input files.\n")
cat(" (If no dialog appears, check BEHIND the RStudio window /\n")
cat("  in the taskbar; otherwise you will be asked to paste paths.)\n")
cat("==============================================================\n")

repeat {
  cpp_file <- choose_file("Select the C++ library file (sim_univLib.cpp)")
  if (grepl("\\.cpp$", cpp_file, ignore.case = TRUE)) break
  cat("    WRONG FILE - this should be a .cpp file (sim_univLib.cpp). Please select again.\n")
}

repeat {
  raw_file <- choose_file("Select the aligned PLINK .raw genotype file for this cohort")
  hdr <- tryCatch(read.table(raw_file, header = TRUE, nrows = 2, stringsAsFactors = FALSE),
                  error = function(e) NULL)
  if (!is.null(hdr) && "IID" %in% colnames(hdr) && ncol(hdr) > 7) break
  cat("    WRONG FILE - a PLINK .raw file starts with FID/IID/PAT/MAT/SEX/PHENOTYPE",
      "followed by SNP columns. Please select again.\n")
}

repeat {
  clin_file <- choose_file("Select the cleaned clinical CSV for this cohort (with pc1-pc3)")
  hdr <- tryCatch(read.csv(clin_file, header = TRUE, nrows = 2, stringsAsFactors = FALSE),
                  error = function(e) NULL)
  if (!is.null(hdr) && all(c("IID", "OS_month", "DEAD") %in% colnames(hdr))) break
  cat("    WRONG FILE - the clinical CSV must contain IID, OS_month and DEAD columns.",
      "Please select again.\n")
}

default_label <- sub("\\.raw$", "", basename(raw_file))
cohort_label <- trimws(readline(
  prompt = sprintf("Type a cohort label for output filenames (Enter = '%s'): ", default_label)))
if (!nzchar(cohort_label)) cohort_label <- default_label

# ------------------------------------------------------------------------------
# 1. Load packages (auto-install if missing), then compile the C++ library
# ------------------------------------------------------------------------------

pkgs <- c("quadprog", "Rcpp", "RcppArmadillo", "glmnet", "survival", "dplyr")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing missing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("\nCompiling C++ library (this can take a minute)...\n")
sourceCpp(cpp_file)

out_dir <- file.path(getwd(), "lambda1se_results")
if (!dir.exists(out_dir)) dir.create(out_dir)

# ------------------------------------------------------------------------------
# 2. Data preparation (identical to the primary analysis)
# ------------------------------------------------------------------------------

df <- read.table(raw_file, header = TRUE, stringsAsFactors = FALSE)
stopifnot("IID" %in% colnames(df))
snp_cols <- colnames(df)[7:ncol(df)]   # PLINK .raw layout
cat(sprintf("Genotypes: %d samples, %d SNP columns\n", nrow(df), length(snp_cols)))

df2 <- read.csv(clin_file, header = TRUE, stringsAsFactors = FALSE)
stopifnot(all(c("IID", "OS_month", "DEAD") %in% colnames(df2)))

df_all <- merge(df[, c("IID", snp_cols)], df2, by = "IID")
cat(sprintf("Merged: %d samples\n", nrow(df_all)))

vars_to_factor <- c("SEX", "smoksort", "early_late", "RADS", "chemotx", "surgery")
for (v in vars_to_factor) if (v %in% names(df_all)) df_all[[v]] <- as.factor(df_all[[v]])

df_all <- df_all %>% filter(!is.na(OS_month), !is.na(DEAD))

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

# Covariate set: identical across all three primary cohort scripts
final_vars <- c(snp_cols,
                "SEX", "AGE", "smoksort",
                "pc1", "pc2", "pc3",
                "early_late",
                "surgery", "chemotx", "RADS")
final_vars <- final_vars[final_vars %in% colnames(df_imputed)]

formula_str <- paste("~", paste(final_vars, collapse = " + "))
X <- model.matrix(as.formula(formula_str), data = df_imputed)
X <- as.matrix(as.data.frame(X))[, -1]

var_names <- colnames(X)
p <- ncol(X)
n <- nrow(X)
cat(sprintf("Design matrix: n = %d, p = %d, events = %d\n\n", n, p, sum(df_imputed$DEAD)))

# ------------------------------------------------------------------------------
# 3. Penalized estimation — lambda.1se  (THE ONLY CHANGE vs primary analysis)
# ------------------------------------------------------------------------------

n_lambda <- 100
nfold <- 5
tol <- 1.0e-6
maxiter <- 50000

cvobj_glmnet <- cv.glmnet(x = X,
                          y = cbind(time = df_imputed$OS_month, status = df_imputed$DEAD),
                          family = "cox", alpha = 1, standardize = FALSE,
                          nfolds = nfold, nlambda = n_lambda)

beta_1se <- as.vector(coef(glmnet(x = X,
                                  y = cbind(time = df_imputed$OS_month, status = df_imputed$DEAD),
                                  family = "cox", alpha = 1,
                                  lambda = cvobj_glmnet$lambda.1se,   # <<< lambda.1se
                                  standardize = FALSE, thresh = tol, maxit = maxiter)))

# For the reply letter, also record the lambda.min fit's nonzero count from the
# SAME cv object (no extra CV randomness in the comparison).
beta_min <- as.vector(coef(glmnet(x = X,
                                  y = cbind(time = df_imputed$OS_month, status = df_imputed$DEAD),
                                  family = "cox", alpha = 1,
                                  lambda = cvobj_glmnet$lambda.min,
                                  standardize = FALSE, thresh = tol, maxit = maxiter)))

cat(sprintf("lambda.min = %.6f  (nonzero coefficients: %d)\n",
            cvobj_glmnet$lambda.min, sum(beta_min != 0)))
cat(sprintf("lambda.1se = %.6f  (nonzero coefficients: %d)\n\n",
            cvobj_glmnet$lambda.1se, sum(beta_1se != 0)))

beta_glmnet <- beta_1se   # initial estimator for debiasing

neg_loglik_glmnet <- 0
neg_dloglik_glmnet <- rep(0, p)
neg_ddloglik_glmnet <- matrix(0, nrow = p, ncol = p)
score_sq <- matrix(0, nrow = p, ncol = p)
neg_loglik_functions_cpp_ext(neg_loglik_glmnet, neg_dloglik_glmnet, neg_ddloglik_glmnet,
                             score_sq, X, df_imputed$OS_month, df_imputed$DEAD, beta_glmnet)
r <- eigen(score_sq)
r$values[r$values <= 1.0e-14] <- 0

# ------------------------------------------------------------------------------
# 4. Cross-validated selection of the projection multiplier
#    (identical grid, folds, seed; inner fits also use lambda.1se)
# ------------------------------------------------------------------------------

n_cv <- 5
n_multi <- 30
multiplier_seq <- exp(seq(from = log(0.001), to = log(5), length.out = n_multi))
alpha_cv <- 0.1
set.seed(12345)

all_cv_idx <- rep(1:n_cv, (n + 1) / n_cv)
all_cv_idx <- sample(all_cv_idx, size = n, replace = TRUE)
all_cvpl2 <- array(NA, length(multiplier_seq))

for (jj in seq_along(multiplier_seq)) {
  if (jj %% 5 == 0) cat("Processing multiplier", jj, "of", n_multi, "...\r")
  multiplier <- multiplier_seq[jj]
  cvpl2 <- 0

  for (k in 1:n_cv) {
    cv_idx <- which(all_cv_idx == k)
    train_x <- X[-cv_idx, ]
    test_x  <- X[cv_idx, ]
    train_time  <- df_imputed$OS_month[-cv_idx]
    test_time   <- df_imputed$OS_month[cv_idx]
    train_delta <- df_imputed$DEAD[-cv_idx]
    test_delta  <- df_imputed$DEAD[cv_idx]

    cvobj_train <- cv.glmnet(x = train_x, y = cbind(time = train_time, status = train_delta),
                             family = "cox", alpha = 1, standardize = FALSE,
                             nfolds = nfold, nlambda = n_lambda)
    beta_train <- as.vector(coef(glmnet(x = train_x,
                                        y = cbind(time = train_time, status = train_delta),
                                        family = "cox", alpha = 1,
                                        lambda = cvobj_train$lambda.1se,   # <<< lambda.1se
                                        standardize = FALSE, thresh = tol, maxit = maxiter)))

    neg_log_train <- 0
    neg_dlog_train <- rep(0, p)
    neg_ddlog_train <- matrix(0, nrow = p, ncol = p)
    score_sq_train <- matrix(0, nrow = p, ncol = p)
    neg_loglik_functions_cpp_ext(neg_log_train, neg_dlog_train, neg_ddlog_train,
                                 score_sq_train, train_x, train_time, train_delta, beta_train)
    r_train <- eigen(score_sq_train)
    r_train$values[r_train$values <= 1.0e-14] <- 0

    b_hat_new <- rep(NA, p)
    se_new <- rep(NA, p)
    mu <- multiplier * sqrt(log(p) / n)
    my_pos <- which(r_train$values > 0)
    my_rank <- sum(r_train$values > 0)
    Dmat <- diag(r_train$values[my_pos])
    dvec <- rep(0, my_rank)
    Amat <- t(rbind(-r_train$vectors[, my_pos] %*% Dmat, r_train$vectors[, my_pos] %*% Dmat))

    for (j in 1:p) {
      Vpos <- r_train$vectors[, my_pos, drop = FALSE]
      e_j <- rep(0, p); e_j[j] <- 1
      ej_proj <- Vpos %*% (t(Vpos) %*% e_j)
      feas_gap <- max(abs(e_j - ej_proj))
      mu_eff <- max(mu, feas_gap + 1e-8)
      bvec <- c(-e_j, e_j) - mu_eff
      res <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec)
      m <- as.vector(r_train$vectors[, my_pos] %*% res$solution) +
        as.vector(as.matrix(r_train$vectors[, -my_pos]) %*% rep(0, p - my_rank))
      b_hat_new[j] <- beta_train[j] - as.numeric(m %*% neg_dlog_train)
      se_new[j] <- if (m[j] > 0) sqrt(m[j] / nrow(train_x)) else NA
    }
    pval_new <- 2 * pnorm(abs(b_hat_new / se_new), lower.tail = FALSE)
    tmp_beta <- b_hat_new * as.numeric(pval_new < (alpha_cv / p))
    cvpl2 <- cvpl2 + loglik_cpp_ext(X = test_x, time = test_time, delta = test_delta, beta = tmp_beta)
  }
  all_cvpl2[jj] <- cvpl2
}
cat("\nCross-validation complete.\n\n")

# ------------------------------------------------------------------------------
# 5. Final debiased estimates (identical algorithm to the primary analysis)
# ------------------------------------------------------------------------------

multiplier <- multiplier_seq[which.max(all_cvpl2)]
cat("Chosen multiplier:", multiplier, "\n")

b_hat_new <- array(NA, p)
se_new <- array(NA, p)
mu_new <- multiplier * sqrt(log(p) / n)
my_pos <- which(r$values > 0)
my_rank <- sum(r$values > 0)
Dmat <- diag(r$values[my_pos])
dvec <- rep(0, my_rank)
Amat <- t(rbind(-r$vectors[, my_pos] %*% Dmat, r$vectors[, my_pos] %*% Dmat))
Vpos <- r$vectors[, my_pos, drop = FALSE]

for (j in 1:p) {
  e_j <- rep(0, p); e_j[j] <- 1
  ej_proj <- Vpos %*% (t(Vpos) %*% e_j)
  feas_gap <- max(abs(e_j - ej_proj))
  mu_eff <- max(mu_new, feas_gap + 1e-8)
  bvec <- (c(-e_j, e_j) - mu_eff * rep(1, 2 * p))
  res <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec)
  m <- as.vector(r$vectors[, my_pos] %*% res$solution) +
    as.vector(as.matrix(r$vectors[, -my_pos]) %*% rep(0, p - my_rank))
  b_hat_new[j] <- beta_glmnet[j] - as.numeric(m %*% neg_dloglik_glmnet)
  se_new[j] <- if (m[j] > 0) sqrt(m[j] / n) else NA
}
pval_new <- 2 * pnorm(abs(b_hat_new / se_new), lower.tail = FALSE)

# ------------------------------------------------------------------------------
# 6. Assemble, FDR-correct, and save
# ------------------------------------------------------------------------------

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

snp_cols_safe <- make.names(snp_cols)
is_snp <- var_names %in% snp_cols |
  var_names %in% snp_cols_safe |
  grepl("^rs\\d+", var_names, ignore.case = TRUE)
if (sum(is_snp) == 0) is_snp <- seq_len(p) <= length(snp_cols)

result$pvalue_FDR <- NA
result$pvalue_FDR[is_snp] <- p.adjust(result$pvalue[is_snp], method = "BH")

res_file <- file.path(out_dir, paste0(cohort_label, "_cox_debiased_results_all_1se.csv"))
write.csv(result, res_file, row.names = FALSE)

run_summary <- data.frame(
  cohort = cohort_label,
  n = n, p = p, events = sum(df_imputed$DEAD),
  lambda_min = cvobj_glmnet$lambda.min,
  lambda_1se = cvobj_glmnet$lambda.1se,
  nonzero_min = sum(beta_min != 0),
  nonzero_1se = sum(beta_1se != 0),
  chosen_multiplier = multiplier,
  stringsAsFactors = FALSE
)
sum_file <- file.path(out_dir, paste0(cohort_label, "_1se_run_summary.csv"))
write.csv(run_summary, sum_file, row.names = FALSE)

cat("\nSaved:\n  ", res_file, "\n  ", sum_file, "\n")
cat(sprintf("\nFor the reply letter: nonzero coefficients %d (lambda.min) vs %d (lambda.1se)\n",
            sum(beta_min != 0), sum(beta_1se != 0)))
cat("Done. Repeat for the remaining cohorts, then run the comparison script.\n")
