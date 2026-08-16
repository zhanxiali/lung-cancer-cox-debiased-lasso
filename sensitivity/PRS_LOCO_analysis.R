rm(list = ls())

# ==============================================================================
# Leave-One-Cohort-Out (LOCO) Polygenic Risk Score Analysis for Overall Survival
# Response to CEBP Reviewer 2, Comment 8 (R2-8)
#
# DESIGN
#   For each target cohort, SNP weights are derived EXCLUSIVELY from the
#   inverse-variance-weighted (IVW) fixed-effect meta-analysis of the debiased
#   lasso estimates of the OTHER TWO cohorts. The weighted score is then tested
#   against overall survival in the held-out target cohort with full covariate
#   adjustment. Weights are therefore independent of the evaluation data,
#   implementing the "different dataset" logic requested by the reviewer.
#
# TWO SCORES PER FOLD
#   PRS_nominal (PRIMARY):     SNPs with training-meta P < 0.05
#                              (fold-specific selection; no leakage)
#   PRS_full    (SENSITIVITY): all weighted SNPs, no selection
#
# ALLELE HARMONIZATION (IMPORTANT)
#   The per-cohort results CSVs are on each cohort's own dosage coding, which
#   for a subset of SNPs is OPPOSITE to the harmonized reference used in the
#   published meta-analysis. Therefore:
#     - Weights are taken from the HARMONIZED per-cohort columns
#       (beta_HSPH / beta_Onco / beta_MGH) of the meta-analysis results file,
#       which are all on the reference (HSPH) effect-allele coding.
#     - For scoring, each target-cohort dosage column is flipped (g -> 2 - g)
#       whenever the cohort's own CSV beta and the harmonized beta for that SNP
#       agree in magnitude but differ in sign. This flip map is inferred
#       automatically and reported; ambiguous SNPs are excluded with a warning.
#
# INPUTS (selected interactively via file.choose(); order is prompted)
#    1.   Meta-analysis results file for ALL meta-analyzable SNPs
#         (columns rsID, beta_HSPH, SE_HSPH, beta_Onco, SE_Onco, beta_MGH,
#          SE_MGH), e.g. path/to/Meta_analysis_results_all.csv
#    2-4. Per-cohort debiased results CSVs (Variable, beta, SE, ...)
#    5-7. Per-cohort aligned PLINK .raw genotype files
#    8-10 Per-cohort cleaned clinical CSVs with pc1-pc3
#
# OUTPUTS (written to ./PRS_results/ under the current working directory)
#   PRS_LOCO_per_fold_results.csv        per-cohort per-SD HR, C-index, quartiles
#   PRS_LOCO_meta_summary.csv            IVW meta of the three per-SD log HRs
#   PRS_LOCO_snp_weights_target_<c>.csv  weights + flip flags (audit trail)
#   PRS_LOCO_KM_<c>_<score>.png          quartile Kaplan-Meier curves
#   PRS_LOCO_forest.png                  forest plot of per-SD HRs + pooled
# ==============================================================================

suppressPackageStartupMessages({
  library(survival)
})
have_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

set.seed(12345)

out_dir <- file.path(getwd(), "PRS_results")
if (!dir.exists(out_dir)) dir.create(out_dir)

cohort_names <- c("HSPH", "Oncoarray", "MGH")
# Column stems of the harmonized betas in the meta results file, by cohort
meta_stem <- c(HSPH = "HSPH", Oncoarray = "Onco", MGH = "MGH")

# ------------------------------------------------------------------------------
# 0. Helpers
# ------------------------------------------------------------------------------

strip_allele <- function(x) sub("_[ACGT]+$", "", x)

norm_path <- function(p) {
  if (is.null(p) || length(p) == 0) return(NULL)
  p <- p[1]
  if (is.na(p) || !nzchar(p)) return(NULL)
  p
}

# Robust file chooser: RStudio picker -> Windows dialog -> base dialog -> paste path
choose_file <- function(msg) {
  cat("\n>>> ", msg, "\n", sep = "")
  flush.console()
  path <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() && rstudioapi::hasFun("selectFile")) {
    path <- norm_path(tryCatch(rstudioapi::selectFile(caption = msg),
                               error = function(e) NULL))
  }
  if (is.null(path) && .Platform$OS.type == "windows" && interactive()) {
    path <- norm_path(tryCatch(utils::choose.files(caption = msg, multi = FALSE),
                               error = function(e) NULL))
  }
  if (is.null(path) && interactive()) {
    path <- norm_path(tryCatch(file.choose(), error = function(e) NULL))
  }
  while (is.null(path) || !file.exists(path)) {
    cat("    (No file selected via dialog.)\n")
    path <- trimws(readline(prompt = "    Paste the FULL file path here and press Enter: "))
    path <- gsub('^"|"$', "", path)
    path <- gsub("\\\\", "/", path)
    if (!nzchar(path) || !file.exists(path)) { cat("    File not found, please try again.\n"); path <- NULL }
  }
  cat("    Selected: ", path, "\n", sep = "")
  path
}

ivw_meta <- function(betas, ses) {
  ok <- !is.na(betas) & !is.na(ses) & ses > 0
  betas <- betas[ok]; ses <- ses[ok]
  if (length(betas) == 0) return(c(beta = NA, se = NA, p = NA, k = 0))
  w <- 1 / ses^2
  b <- sum(w * betas) / sum(w)
  se <- sqrt(1 / sum(w))
  c(beta = b, se = se, p = 2 * pnorm(-abs(b / se)), k = length(betas))
}

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

c_index <- function(fit) as.numeric(summary(fit)$concordance[1])

# ------------------------------------------------------------------------------
# 1. Interactive file selection (no hardcoded paths)
# ------------------------------------------------------------------------------

cat("==============================================================\n")
cat(" Please select input files when prompted (10 files total).\n")
cat("==============================================================\n\n")

# Choose a file and verify it looks right (reads only the header rows);
# on a wrong pick, explain and re-prompt instead of aborting.
choose_validated <- function(msg, check, hint) {
  repeat {
    f <- choose_file(msg)
    ok <- isTRUE(tryCatch(check(f), error = function(e) FALSE))
    if (ok) return(f)
    cat("    WRONG FILE - ", hint, "\n    Please select again.\n", sep = "")
  }
}
has_csv_cols <- function(cols) function(f) {
  h <- read.csv(f, header = TRUE, nrows = 2, stringsAsFactors = FALSE)
  all(cols %in% colnames(h))
}
looks_like_raw <- function(f) {
  h <- read.table(f, header = TRUE, nrows = 2, stringsAsFactors = FALSE)
  "IID" %in% colnames(h) && ncol(h) > 7
}

meta_file <- choose_validated(
  "[1/10] Select 'Meta_analysis_results_all.csv' (the META file, NOT a per-cohort file)",
  has_csv_cols(c("rsID", "beta_HSPH", "beta_Onco", "beta_MGH")),
  "This must be Meta_analysis_results_all.csv with beta_HSPH/beta_Onco/beta_MGH columns; a per-cohort results file does not have them.")

res_files <- raw_files <- clin_files <- setNames(vector("list", 3), cohort_names)
for (cn in cohort_names) {
  res_files[[cn]] <- choose_validated(
    sprintf("[%s] Select the DEBIASED RESULTS csv (e.g. %s..._cox_debiased_results_all.csv)", cn, cn),
    has_csv_cols(c("Variable", "beta", "SE")),
    "This must be a per-cohort debiased results file with Variable/beta/SE columns.")
}
for (cn in cohort_names) {
  raw_files[[cn]] <- choose_validated(
    sprintf("[%s] Select the aligned PLINK .raw GENOTYPE file", cn),
    looks_like_raw,
    "A PLINK .raw file starts with FID/IID/PAT/MAT/SEX/PHENOTYPE followed by SNP columns.")
}
for (cn in cohort_names) {
  clin_files[[cn]] <- choose_validated(
    sprintf("[%s] Select the cleaned CLINICAL csv (with pc1-pc3)", cn),
    has_csv_cols(c("IID", "OS_month", "DEAD")),
    "The clinical CSV must contain IID, OS_month and DEAD columns.")
}

# ------------------------------------------------------------------------------
# 2. Load harmonized meta results (weight source) and per-cohort CSVs
# ------------------------------------------------------------------------------

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
need_cols <- c("rsID", paste0("beta_", meta_stem), paste0("SE_", meta_stem))
stopifnot(all(need_cols %in% colnames(meta)))
for (cc in setdiff(need_cols, "rsID")) meta[[cc]] <- suppressWarnings(as.numeric(meta[[cc]]))
cat(sprintf("Meta results loaded: %d SNPs (harmonized to reference coding)\n", nrow(meta)))

res_list <- list()
for (cn in cohort_names) {
  r <- read.csv(res_files[[cn]], stringsAsFactors = FALSE)
  stopifnot(all(c("Variable", "beta") %in% colnames(r)))
  r <- r[grepl("^rs[0-9]+_[ACGT]+$", r$Variable), , drop = FALSE]
  r$rsID <- strip_allele(r$Variable)
  if (!"Variable_original" %in% colnames(r)) r$Variable_original <- NA_character_
  res_list[[cn]] <- r[, c("rsID", "Variable", "Variable_original", "beta")]
  cat(sprintf("[%s] per-cohort results: %d SNPs\n", cn, nrow(r)))
}

# ------------------------------------------------------------------------------
# 3. Infer per-cohort dosage flip map (own coding vs harmonized reference)
#    flip = TRUE  -> score with (2 - dosage)
# ------------------------------------------------------------------------------

tol <- 1e-8
flip_list <- list()
for (cn in cohort_names) {
  r <- res_list[[cn]]
  hb <- meta[[paste0("beta_", meta_stem[cn])]][match(r$rsID, meta$rsID)]
  status <- rep(NA_character_, nrow(r))
  status[!is.na(hb) & abs(r$beta - hb) < tol] <- "same"
  status[!is.na(hb) & abs(r$beta + hb) < tol] <- "flip"
  status[!is.na(hb) & is.na(status)] <- "ambiguous"
  flip_list[[cn]] <- data.frame(rsID = r$rsID, Variable = r$Variable,
                                Variable_original = r$Variable_original,
                                status = status, stringsAsFactors = FALSE)
  cat(sprintf("[%s] flip map: same=%d, flipped=%d, ambiguous=%d, not-in-meta=%d\n",
              cn, sum(status == "same", na.rm = TRUE), sum(status == "flip", na.rm = TRUE),
              sum(status == "ambiguous", na.rm = TRUE), sum(is.na(status))))
  if (any(status == "ambiguous", na.rm = TRUE))
    warning(sprintf("[%s] %d SNP(s) with ambiguous flip status will be EXCLUDED from scoring.",
                    cn, sum(status == "ambiguous", na.rm = TRUE)))
}

# ------------------------------------------------------------------------------
# 4. Load genotypes + clinical data per cohort
# ------------------------------------------------------------------------------

geno_list <- list()
for (cn in cohort_names) {
  g <- read.table(raw_files[[cn]], header = TRUE, stringsAsFactors = FALSE)
  stopifnot("IID" %in% colnames(g))
  snp_cols <- colnames(g)[7:ncol(g)]   # PLINK .raw layout
  cl <- read.csv(clin_files[[cn]], stringsAsFactors = FALSE)
  stopifnot(all(c("IID", "OS_month", "DEAD") %in% colnames(cl)))
  m <- merge(g[, c("IID", snp_cols)], cl, by = "IID")
  cat(sprintf("[%s] merged genotype+clinical: n = %d, SNP columns = %d\n",
              cn, nrow(m), length(snp_cols)))
  geno_list[[cn]] <- list(data = m, snp_cols = snp_cols)
}

# Match one rsID to a genotype column of a cohort:
# try Variable (rsID_allele), then Variable_original (chr.pos.ref.alt_allele),
# then rsID prefix on stripped names.
match_geno_col <- function(rsid, fliprow, snp_cols) {
  cand <- c(fliprow$Variable, fliprow$Variable_original)
  cand <- cand[!is.na(cand)]
  hit <- cand[cand %in% snp_cols]
  if (length(hit) > 0) return(hit[1])
  idx <- which(strip_allele(snp_cols) == rsid)
  if (length(idx) > 0) return(snp_cols[idx[1]])
  NA_character_
}

# ------------------------------------------------------------------------------
# 5. Covariates (aligned with the primary analysis; constants auto-dropped)
# ------------------------------------------------------------------------------

covars_wanted <- c("AGE", "SEX", "smoksort", "early_late",
                   "surgery", "chemotx", "RADS", "immunotherapy",
                   "pc1", "pc2", "pc3")

build_covars <- function(d) {
  keep <- covars_wanted[covars_wanted %in% colnames(d)]
  keep[vapply(keep, function(v) length(unique(na.omit(d[[v]]))) > 1, logical(1))]
}

# ------------------------------------------------------------------------------
# 6. LOCO folds
# ------------------------------------------------------------------------------

fold_results <- list()

for (target in cohort_names) {

  training <- setdiff(cohort_names, target)
  cat(sprintf("\n============ FOLD: target = %s | training = %s + %s ============\n",
              target, training[1], training[2]))

  # ---- 6a. Training-only IVW weights from HARMONIZED betas ----------------
  b1 <- meta[[paste0("beta_", meta_stem[training[1]])]]
  s1 <- meta[[paste0("SE_",   meta_stem[training[1]])]]
  b2 <- meta[[paste0("beta_", meta_stem[training[2]])]]
  s2 <- meta[[paste0("SE_",   meta_stem[training[2]])]]

  wtab <- do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
    m <- ivw_meta(c(b1[i], b2[i]), c(s1[i], s2[i]))
    data.frame(rsID = meta$rsID[i], beta_train = m["beta"], se_train = m["se"],
               p_train = m["p"], k_train = m["k"], stringsAsFactors = FALSE)
  }))
  wtab <- wtab[!is.na(wtab$beta_train), , drop = FALSE]

  # ---- 6b. Restrict to SNPs scoreable in the target ----------------------
  fmap <- flip_list[[target]]
  wtab <- merge(wtab, fmap, by = "rsID")
  wtab <- wtab[wtab$status %in% c("same", "flip"), , drop = FALSE]

  tgt <- geno_list[[target]]
  wtab$geno_col <- vapply(seq_len(nrow(wtab)), function(i)
    match_geno_col(wtab$rsID[i], wtab[i, ], tgt$snp_cols), character(1))
  wtab <- wtab[!is.na(wtab$geno_col), , drop = FALSE]
  cat(sprintf("SNPs with training weights, unambiguous coding, present in target: %d\n",
              nrow(wtab)))

  set_nominal <- wtab[wtab$p_train < 0.05, , drop = FALSE]
  set_full    <- wtab
  cat(sprintf("  PRS_nominal (training P<0.05): %d SNPs [%s]\n",
              nrow(set_nominal), paste(set_nominal$rsID, collapse = ", ")))
  cat(sprintf("  PRS_full    (no selection):    %d SNPs\n", nrow(set_full)))

  write.csv(wtab, file.path(out_dir, sprintf("PRS_LOCO_snp_weights_target_%s.csv", target)),
            row.names = FALSE)

  # ---- 6c. Prepare target data -------------------------------------------
  d <- tgt$data
  gc_needed <- unique(set_full$geno_col)
  d[gc_needed] <- lapply(d[gc_needed], function(x) as.numeric(as.character(x)))
  d <- impute_mode_mean(d)

  covs <- build_covars(d)
  cat(sprintf("  Covariates used: %s\n", paste(covs, collapse = ", ")))

  compute_prs <- function(snpset) {
    G <- as.matrix(d[, snpset$geno_col, drop = FALSE])
    flip <- snpset$status == "flip"
    if (any(flip)) G[, flip] <- 2 - G[, flip]   # align dosages to reference coding
    as.numeric(G %*% snpset$beta_train)
  }

  fit_one <- function(prs_raw, label, n_snps) {
    d$PRS_z <- as.numeric(scale(prs_raw))
    f_base <- as.formula(paste("Surv(OS_month, DEAD) ~", paste(covs, collapse = " + ")))
    f_prs  <- as.formula(paste("Surv(OS_month, DEAD) ~ PRS_z +", paste(covs, collapse = " + ")))
    fit_base <- coxph(f_base, data = d)
    fit_prs  <- coxph(f_prs,  data = d)
    s <- summary(fit_prs)
    hr  <- s$conf.int["PRS_z", "exp(coef)"]
    lcl <- s$conf.int["PRS_z", "lower .95"]
    ucl <- s$conf.int["PRS_z", "upper .95"]
    pv  <- s$coefficients["PRS_z", "Pr(>|z|)"]
    lrt_p <- anova(fit_base, fit_prs)[2, "Pr(>|Chi|)"]

    d$PRS_q <- cut(d$PRS_z, breaks = quantile(d$PRS_z, probs = seq(0, 1, 0.25)),
                   include.lowest = TRUE, labels = paste0("Q", 1:4))
    f_q <- as.formula(paste("Surv(OS_month, DEAD) ~ PRS_q +", paste(covs, collapse = " + ")))
    sq <- summary(coxph(f_q, data = d))
    q4 <- grep("PRS_qQ4", rownames(sq$conf.int))

    km <- survfit(Surv(OS_month, DEAD) ~ PRS_q, data = d)
    png(file.path(out_dir, sprintf("PRS_LOCO_KM_%s_%s.png", target, label)),
        width = 1600, height = 1200, res = 200)
    plot(km, col = c("#1b9e77", "#7570b3", "#d95f02", "#e7298a"), lwd = 2,
         xlab = "Months since diagnosis", ylab = "Overall survival probability",
         main = sprintf("Target: %s | %s (%d SNPs)", target, label, n_snps))
    legend("topright", legend = paste0("PRS Q", 1:4),
           col = c("#1b9e77", "#7570b3", "#d95f02", "#e7298a"), lwd = 2, bty = "n")
    dev.off()

    data.frame(
      target = target, score = label, n = nrow(d), events = sum(d$DEAD),
      n_snps = n_snps,
      log_hr = log(hr), se_log_hr = s$coefficients["PRS_z", "se(coef)"],
      HR_perSD = hr, CI_lower = lcl, CI_upper = ucl, P = pv, LRT_P = lrt_p,
      HR_Q4_vs_Q1 = sq$conf.int[q4, "exp(coef)"],
      Q4_CI_lower = sq$conf.int[q4, "lower .95"],
      Q4_CI_upper = sq$conf.int[q4, "upper .95"],
      Q4_P = sq$coefficients[q4, "Pr(>|z|)"],
      C_base = c_index(fit_base), C_with_PRS = c_index(fit_prs),
      delta_C = c_index(fit_prs) - c_index(fit_base),
      stringsAsFactors = FALSE
    )
  }

  res_rows <- list()
  if (nrow(set_nominal) >= 2) {
    res_rows[["nominal"]] <- fit_one(compute_prs(set_nominal), "PRS_nominal", nrow(set_nominal))
  } else {
    cat("  [warn] <2 SNPs pass training P<0.05; PRS_nominal skipped for this fold\n")
  }
  res_rows[["full"]] <- fit_one(compute_prs(set_full), "PRS_full", nrow(set_full))

  fold_results[[target]] <- do.call(rbind, res_rows)
}

per_fold <- do.call(rbind, fold_results)
rownames(per_fold) <- NULL
write.csv(per_fold, file.path(out_dir, "PRS_LOCO_per_fold_results.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. Pool per-SD log HRs across folds (IVW fixed; DerSimonian-Laird random)
# ------------------------------------------------------------------------------

pool <- function(sub) {
  m <- ivw_meta(sub$log_hr, sub$se_log_hr)
  w <- 1 / sub$se_log_hr^2
  Q <- sum(w * (sub$log_hr - m["beta"])^2)
  k <- nrow(sub)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  wr <- 1 / (sub$se_log_hr^2 + tau2)
  br <- sum(wr * sub$log_hr) / sum(wr)
  ser <- sqrt(1 / sum(wr))
  I2 <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  data.frame(
    score = sub$score[1], k = k,
    HR_fixed = exp(m["beta"]),
    CI_lower_fixed = exp(m["beta"] - 1.96 * m["se"]),
    CI_upper_fixed = exp(m["beta"] + 1.96 * m["se"]),
    P_fixed = m["p"],
    HR_random = exp(br),
    CI_lower_random = exp(br - 1.96 * ser),
    CI_upper_random = exp(br + 1.96 * ser),
    P_random = 2 * pnorm(-abs(br / ser)),
    I2 = I2, Q = Q, Q_p = pchisq(Q, k - 1, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

meta_summary <- do.call(rbind, lapply(split(per_fold, per_fold$score), pool))
rownames(meta_summary) <- NULL
write.csv(meta_summary, file.path(out_dir, "PRS_LOCO_meta_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 8. Forest plot (primary score)
# ------------------------------------------------------------------------------

if (have_ggplot) {
  library(ggplot2)
  prim <- per_fold[per_fold$score == "PRS_nominal", , drop = FALSE]
  if (nrow(prim) == 0) prim <- per_fold[per_fold$score == "PRS_full", , drop = FALSE]
  ms <- meta_summary[meta_summary$score == prim$score[1], ]
  fp <- rbind(
    data.frame(label = prim$target, HR = prim$HR_perSD,
               lo = prim$CI_lower, hi = prim$CI_upper, type = "Cohort"),
    data.frame(label = "Pooled (IVW fixed)", HR = ms$HR_fixed,
               lo = ms$CI_lower_fixed, hi = ms$CI_upper_fixed, type = "Pooled")
  )
  fp$label <- factor(fp$label, levels = rev(fp$label))
  p <- ggplot(fp, aes(x = HR, y = label, color = type)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c(Cohort = "#2166ac", Pooled = "#b2182b")) +
    labs(x = "Hazard ratio per SD of PRS (95% CI)", y = NULL,
         title = sprintf("Leave-one-cohort-out PRS and overall survival (%s)",
                         prim$score[1])) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  ggsave(file.path(out_dir, "PRS_LOCO_forest.png"), p, width = 7, height = 3.2, dpi = 300)
}

# ------------------------------------------------------------------------------
# 9. Console summary + suggested wording (verify numbers before use)
# ------------------------------------------------------------------------------

cat("\n\n==================== SUMMARY ====================\n")
print(per_fold[, c("target", "score", "n", "events", "n_snps",
                   "HR_perSD", "CI_lower", "CI_upper", "P", "delta_C")], digits = 3)
cat("\n---- Pooled across folds ----\n")
print(meta_summary, digits = 3)

prim_row <- if ("PRS_nominal" %in% meta_summary$score)
  meta_summary[meta_summary$score == "PRS_nominal", ] else meta_summary[1, ]

cat("\n---- Suggested Results wording (VERIFY before use) ----\n")
cat(sprintf(
  "In leave-one-cohort-out analyses, a weighted PRS built from the debiased\n"))
cat(sprintf(
  "estimates of the two training cohorts (training-meta P<0.05 variants) was\n"))
cat(sprintf(
  "associated with overall survival in the held-out cohort (pooled HR per SD =\n"))
cat(sprintf(
  "%.3f, 95%% CI %.3f-%.3f, P = %.3g; I2 = %.1f%%), with per-cohort change in\n",
  prim_row$HR_fixed, prim_row$CI_lower_fixed, prim_row$CI_upper_fixed,
  prim_row$P_fixed, prim_row$I2))
cat(sprintf(
  "Harrell's C reported in Supplementary Table Sx.\n"))

cat("\nAll outputs written to:", out_dir, "\n")
cat("Done.\n")
