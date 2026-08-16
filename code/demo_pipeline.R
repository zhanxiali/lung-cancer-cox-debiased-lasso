#!/usr/bin/env Rscript
# End-to-end run of the analytic pipeline on the simulated dataset.
#
# The per-cohort scripts in this repository prompt for their inputs interactively.
# This driver takes the same steps without prompting so that the pipeline can run
# unattended inside a compute capsule.
#
# Inputs  : example/example_genotypes.raw, example/example_clinical.csv
# Outputs : results/demo_debiased_results.csv, results/demo_session_info.txt
#
# Arguments (all optional, passed as name=value):
#   n_multi   number of projection multipliers to evaluate   (demo 8, full 30)
#   n_cv      outer cross-validation folds                    (demo 3, full 5)
#   out       output directory                                (default ../results)

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, default) {
  hit <- grep(paste0("^", k, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  val <- sub(paste0("^", k, "="), "", hit[1])
  if (is.numeric(default)) as.numeric(val) else val
}
n_multi <- getarg("n_multi", 8)
n_cv    <- getarg("n_cv", 3)
outdir  <- getarg("out", "../results")

suppressPackageStartupMessages({
  library(quadprog); library(Rcpp); library(RcppArmadillo)
  library(glmnet);   library(survival)
})

here <- function(...) file.path(dirname(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ...)
if (is.na(here())) here <- function(...) file.path(".", ...)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Paths: environment variables take precedence, so the capsule can point at
# /data without editing the script. Otherwise fall back to the repository layout.
env_or <- function(var, fallback) {
  v <- Sys.getenv(var, unset = "")
  if (nzchar(v) && file.exists(v)) v else fallback
}
cpp_path  <- env_or("SIM_UNIVLIB", here("..", "core", "sim_univLib.cpp"))
geno_path <- env_or("GENO_RAW",    here("..", "example", "example_genotypes.raw"))
clin_path <- env_or("CLIN_CSV",    here("..", "example", "example_clinical.csv"))
for (f in c(cpp_path, geno_path, clin_path))
  if (!file.exists(f)) stop("File not found: ", f, call. = FALSE)

message("Compiling C++ routines")
sourceCpp(cpp_path)

message("Reading simulated data")
geno <- read.table(geno_path, header = TRUE, stringsAsFactors = FALSE)
clin <- read.csv(clin_path, stringsAsFactors = FALSE)
snp_cols <- colnames(geno)[7:ncol(geno)]
dat <- merge(geno[, c("IID", snp_cols)], clin, by = "IID")
message(sprintf("  n = %d, SNPs = %d, events = %d", nrow(dat), length(snp_cols), sum(dat$DEAD)))

for (v in c("SEX", "smoksort", "early_late", "RADS", "chemotx", "surgery"))
  if (v %in% names(dat)) dat[[v]] <- as.factor(dat[[v]])
dat[snp_cols] <- lapply(dat[snp_cols], function(x) as.numeric(as.character(x)))

final_vars <- c(snp_cols, "SEX", "AGE", "smoksort", "pc1", "pc2", "pc3",
                "early_late", "surgery", "chemotx", "RADS")
final_vars <- final_vars[final_vars %in% colnames(dat)]
X <- as.matrix(as.data.frame(model.matrix(as.formula(paste("~", paste(final_vars, collapse = " + "))),
                                          data = dat)))[, -1]
p <- ncol(X); n <- nrow(X)
message(sprintf("  design matrix: n = %d, p = %d", n, p))

tol <- 1e-6; maxiter <- 50000

message("Penalized Cox fit")
y <- cbind(time = dat$OS_month, status = dat$DEAD)
cvfit <- cv.glmnet(X, y, family = "cox", alpha = 1, standardize = FALSE,
                   nfolds = 5, nlambda = 100)
beta0 <- as.vector(coef(glmnet(X, y, family = "cox", alpha = 1,
                               lambda = cvfit$lambda.min, standardize = FALSE,
                               thresh = tol, maxit = maxiter)))
message(sprintf("  lambda.min = %.5f, nonzero coefficients = %d",
                cvfit$lambda.min, sum(beta0 != 0)))

nl <- 0; dl <- rep(0, p); ddl <- matrix(0, p, p); ss <- matrix(0, p, p)
neg_loglik_functions_cpp_ext(nl, dl, ddl, ss, X, dat$OS_month, dat$DEAD, beta0)
r <- eigen(ss); r$values[r$values <= 1e-14] <- 0

message(sprintf("Tuning the projection multiplier (%d values, %d folds)", n_multi, n_cv))
multipliers <- exp(seq(log(0.001), log(5), length.out = n_multi))
set.seed(12345)
fold <- sample(rep_len(1:n_cv, n))
cvpl <- numeric(length(multipliers))

debias <- function(Xm, dlm, rm, mult, nrow_used) {
  pos <- which(rm$values > 0)
  Dmat <- diag(rm$values[pos]); dvec <- rep(0, length(pos))
  Amat <- t(rbind(-rm$vectors[, pos] %*% Dmat, rm$vectors[, pos] %*% Dmat))
  V <- rm$vectors[, pos, drop = FALSE]
  mu <- mult * sqrt(log(ncol(Xm)) / nrow_used)
  b <- se <- rep(NA_real_, ncol(Xm))
  for (j in seq_len(ncol(Xm))) {
    e <- rep(0, ncol(Xm)); e[j] <- 1
    gap <- max(abs(e - V %*% (t(V) %*% e)))
    bvec <- c(-e, e) - max(mu, gap + 1e-8)
    sol <- tryCatch(solve.QP(Dmat, dvec, Amat, bvec), error = function(z) NULL)
    if (is.null(sol)) next
    m <- as.vector(rm$vectors[, pos] %*% sol$solution)
    b[j] <- NA; se[j] <- if (m[j] > 0) sqrt(m[j] / nrow_used) else NA
    attr(b, "m") <- NULL
    b[j] <- -as.numeric(m %*% dlm)
  }
  list(shift = b, se = se)
}

for (i in seq_along(multipliers)) {
  total <- 0
  for (k in seq_len(n_cv)) {
    idx <- which(fold == k)
    Xtr <- X[-idx, ]; Xte <- X[idx, ]
    ttr <- dat$OS_month[-idx]; dtr <- dat$DEAD[-idx]
    tte <- dat$OS_month[idx];  dte <- dat$DEAD[idx]
    cvk <- cv.glmnet(Xtr, cbind(time = ttr, status = dtr), family = "cox",
                     alpha = 1, standardize = FALSE, nfolds = 5, nlambda = 100)
    bk <- as.vector(coef(glmnet(Xtr, cbind(time = ttr, status = dtr), family = "cox",
                                alpha = 1, lambda = cvk$lambda.min,
                                standardize = FALSE, thresh = tol, maxit = maxiter)))
    nl2 <- 0; dl2 <- rep(0, p); ddl2 <- matrix(0, p, p); ss2 <- matrix(0, p, p)
    neg_loglik_functions_cpp_ext(nl2, dl2, ddl2, ss2, Xtr, ttr, dtr, bk)
    r2 <- eigen(ss2); r2$values[r2$values <= 1e-14] <- 0
    dd <- debias(Xtr, dl2, r2, multipliers[i], nrow(Xtr))
    bh <- bk + dd$shift
    pv <- 2 * pnorm(abs(bh / dd$se), lower.tail = FALSE)
    bh[is.na(pv) | pv >= 0.1 / p] <- 0
    total <- total + loglik_cpp_ext(X = Xte, time = tte, delta = dte, beta = bh)
  }
  cvpl[i] <- total
  message(sprintf("  multiplier %2d/%d: %.4f -> cvpl %.2f", i, length(multipliers),
                  multipliers[i], total))
}
mult <- multipliers[which.max(cvpl)]
message(sprintf("Selected multiplier: %.4f", mult))

message("Final debiased estimation")
dd <- debias(X, dl, r, mult, n)
b_hat <- beta0 + dd$shift
pval <- 2 * pnorm(abs(b_hat / dd$se), lower.tail = FALSE)

res <- data.frame(Variable = colnames(X), beta = b_hat, SE = dd$se, pvalue = pval,
                  HR = exp(b_hat),
                  HR_lower = exp(b_hat - 1.96 * dd$se),
                  HR_upper = exp(b_hat + 1.96 * dd$se),
                  stringsAsFactors = FALSE)
is_snp <- grepl("^rsSIM", res$Variable)
res$pvalue_FDR <- NA
res$pvalue_FDR[is_snp] <- p.adjust(res$pvalue[is_snp], method = "BH")
res <- res[order(res$pvalue), ]

write.csv(res, file.path(outdir, "demo_debiased_results.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "demo_session_info.txt"))

message("\nTop 10 by P value:")
print(head(res[, c("Variable", "beta", "SE", "pvalue", "HR")], 10), digits = 3)
message(sprintf("\nWritten to %s", normalizePath(outdir, mustWork = FALSE)))
