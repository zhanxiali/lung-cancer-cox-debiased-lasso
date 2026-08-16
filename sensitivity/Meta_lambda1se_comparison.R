rm(list = ls())

# ==============================================================================
# Lambda.1se Meta-Analysis and Robustness Comparison  (Reviewer 2, Comment 2)
# v2 — robust interactive file selection
#
# HOW TO RUN
#   Open in RStudio (or R GUI) and Source. A file dialog pops up for each
#   input. NOTE: on Windows the dialog sometimes opens BEHIND the RStudio
#   window - check the taskbar if nothing seems to happen. If no dialog can
#   be shown, the script asks you to PASTE the full file path instead.
#
# WHAT THIS SCRIPT DOES
#   1. Meta-analyzes the per-cohort lambda.1se debiased estimates using the
#      same inverse-variance-weighted fixed-effect model as the primary
#      analysis (random-effects REML as sensitivity; BH-FDR across SNPs).
#   2. Compares against the published lambda.min meta-analysis and produces
#      the supplementary comparison table plus every number needed for the
#      reply letter.
#
# ALLELE HARMONIZATION (IMPORTANT)
#   Per-cohort CSVs are on each cohort's own .raw dosage coding, which for a
#   subset of SNPs is opposite to the harmonized reference of the published
#   meta-analysis file. Coding is a property of the .raw file, so it is
#   IDENTICAL between the lambda.min and lambda.1se runs of the same cohort.
#   This script infers, per SNP and cohort, whether the cohort coding is
#   flipped relative to the published reference (lambda.min cohort beta vs
#   harmonized beta: equal magnitude, same/opposite sign) and applies the
#   same flip to the lambda.1se betas before meta-analysis. Ambiguous SNPs
#   are excluded with a warning.
#
# INPUTS (selected interactively)
#   1.    Published meta-analysis results file (lambda.min reference;
#         columns rsID, beta_HSPH/Onco/MGH, SE_*)
#   2-4.  Per-cohort LAMBDA.MIN debiased results CSVs (primary analysis)
#   5-7.  Per-cohort LAMBDA.1SE debiased results CSVs (from the rerun script)
#   8-10. (Optional) the three _1se_run_summary.csv files
#
# OUTPUTS (written to ./lambda1se_results/ under the current directory)
#   Meta_analysis_results_all_1se.csv
#   SuppTable_lambda_min_vs_1se.csv   (11 significant SNPs sorted first)
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

  # 1) RStudio's own file picker (Desktop AND Server; cannot hide behind windows)
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
    path <- gsub('^"|"$', "", path)
    path <- gsub("\\\\", "/", path)
    if (!nzchar(path) || !file.exists(path)) { cat("    File not found, please try again.\n"); path <- NULL }
  }
  cat("    Selected: ", path, "\n", sep = "")
  path
}

ok_metafor <- requireNamespace("metafor", quietly = TRUE)

# Choose a CSV and verify it has the required columns; if not, explain what
# was likely mis-selected and re-prompt instead of aborting.
choose_csv_validated <- function(msg, required_cols, hint) {
  repeat {
    f <- choose_file(msg)
    d <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(d)) { cat("    Could not read this file as a CSV. ", hint, "\n", sep = ""); next }
    miss <- setdiff(required_cols, colnames(d))
    if (length(miss) == 0) { attr(d, "source_path") <- f; return(d) }
    cat("    WRONG FILE - missing column(s): ", paste(miss, collapse = ", "),
        "\n    ", hint, "\n    Please select again.\n", sep = "")
  }
}

out_dir <- file.path(getwd(), "lambda1se_results")
if (!dir.exists(out_dir)) dir.create(out_dir)

cohort_names <- c("HSPH", "Oncoarray", "MGH")
meta_stem <- c(HSPH = "HSPH", Oncoarray = "Onco", MGH = "MGH")

strip_allele <- function(x) sub("_[ACGT]+$", "", x)

# Process a LAMBDA.MIN per-cohort table: keep rs-style SNP rows and, when
# present, the Variable_original column (the PLINK-style raw column name,
# e.g. X1.9498113.T.C_C) which is needed to map the lambda.1se output back
# to rsIDs.
read_min_df <- function(r) {
  r <- r[grepl("^rs[0-9]+_[ACGT]+$", r$Variable), , drop = FALSE]
  r$rsID <- strip_allele(r$Variable)
  if (!"Variable_original" %in% colnames(r)) r$Variable_original <- NA_character_
  r[, c("rsID", "Variable", "Variable_original", "beta", "SE", "pvalue")]
}

# Map a LAMBDA.1SE per-cohort table to rsIDs. The rerun script keeps the raw
# design-matrix column names, which are rs-style for cohorts whose .raw file
# uses rsIDs (MGH) but PLINK position-style (X<chr>.<pos>.<ref>.<alt>_<a>)
# for the others. Match against the min table's Variable AND
# Variable_original; fall back to the rs-pattern for directly-named SNPs.
map_se_df <- function(se, minc) {
  idx <- match(se$Variable, minc$Variable)
  idx2 <- match(se$Variable, minc$Variable_original)
  idx[is.na(idx)] <- idx2[is.na(idx)]
  rs <- minc$rsID[idx]
  fb <- strip_allele(se$Variable)
  use_fb <- is.na(rs) & grepl("^rs[0-9]+$", fb)
  rs[use_fb] <- fb[use_fb]
  se$rsID <- rs
  se <- se[!is.na(se$rsID), , drop = FALSE]
  se[!duplicated(se$rsID), c("rsID", "Variable", "beta", "SE", "pvalue")]
}

# ------------------------------------------------------------------------------
# 1. Interactive file selection
# ------------------------------------------------------------------------------

cat("==============================================================\n")
cat(" Lambda.1se meta + comparison - select 7 files (plus 3 optional).\n")
cat("==============================================================\n")

need <- c("rsID", paste0("beta_", meta_stem), paste0("SE_", meta_stem))
meta_min <- choose_csv_validated(
  "Select 'Meta_analysis_results_all.csv' (the PUBLISHED lambda.min META file, NOT a per-cohort file)",
  need,
  "This must be Meta_analysis_results_all.csv - it has beta_HSPH/beta_Onco/beta_MGH columns. A per-cohort results file (e.g. HSPH984_cox_debiased_results_all.csv) does not.")
for (cc in setdiff(need, "rsID")) meta_min[[cc]] <- suppressWarnings(as.numeric(meta_min[[cc]]))

cohort_cols <- c("Variable", "beta", "SE", "pvalue")
min_list <- list()
for (cn in cohort_names) {
  d <- choose_csv_validated(
    sprintf("Select the LAMBDA.MIN per-cohort results csv for %s (e.g. %s..._cox_debiased_results_all.csv)", cn, cn),
    cohort_cols,
    "This must be a per-cohort debiased results file with Variable/beta/SE/pvalue columns.")
  min_list[[cn]] <- read_min_df(d)
  cat(sprintf("    [%s] lambda.min SNP rows: %d\n", cn, nrow(min_list[[cn]])))
}

se_list <- list()
for (cn in cohort_names) {
  d <- choose_csv_validated(
    sprintf("Select the LAMBDA.1SE per-cohort results csv for %s (from ./lambda1se_results/, ends in _1se.csv)", cn),
    cohort_cols,
    "This must be a *_cox_debiased_results_all_1se.csv file produced by the lambda.1se rerun script.")
  se_list[[cn]] <- map_se_df(d, min_list[[cn]])
  cat(sprintf("    [%s] lambda.1se SNP rows mapped to rsIDs: %d\n", cn, nrow(se_list[[cn]])))
  if (nrow(se_list[[cn]]) == 0)
    warning(sprintf("[%s] no lambda.1se rows could be mapped - check that the min and 1se files are from the same cohort.", cn))
}

nz_line <- NULL
ans <- readline(prompt = "\nLoad the three _1se_run_summary.csv files for nonzero counts? (y/n): ")
if (tolower(trimws(ans)) == "y") {
  nz <- do.call(rbind, lapply(seq_len(3), function(i) {
    choose_csv_validated(sprintf("Select run summary file %d of 3 (*_1se_run_summary.csv)", i),
                         c("nonzero_min", "nonzero_1se"),
                         "This must be a *_1se_run_summary.csv file produced by the lambda.1se rerun script.")
  }))
  nz_line <- sprintf("Nonzero lasso coefficients: %d-%d (lambda.min) vs %d-%d (lambda.1se) across cohorts",
                     min(nz$nonzero_min), max(nz$nonzero_min),
                     min(nz$nonzero_1se), max(nz$nonzero_1se))
}

# ------------------------------------------------------------------------------
# 2. Flip inference (cohort coding vs published harmonized reference)
# ------------------------------------------------------------------------------

tol <- 1e-8
flip_list <- list()
for (cn in cohort_names) {
  r <- min_list[[cn]]
  hb <- meta_min[[paste0("beta_", meta_stem[cn])]][match(r$rsID, meta_min$rsID)]
  status <- rep(NA_character_, nrow(r))
  status[!is.na(hb) & abs(r$beta - hb) < tol] <- "same"
  status[!is.na(hb) & abs(r$beta + hb) < tol] <- "flip"
  status[!is.na(hb) & is.na(status)] <- "ambiguous"
  flip_list[[cn]] <- data.frame(rsID = r$rsID, status = status, stringsAsFactors = FALSE)
  cat(sprintf("[%s] flip map: same=%d, flipped=%d, ambiguous=%d, not-in-meta=%d\n",
              cn, sum(status == "same", na.rm = TRUE), sum(status == "flip", na.rm = TRUE),
              sum(status == "ambiguous", na.rm = TRUE), sum(is.na(status))))
  if (any(status == "ambiguous", na.rm = TRUE))
    warning(sprintf("[%s] ambiguous SNPs will be excluded from the 1se meta-analysis.", cn))
}

harm_list <- list()
for (cn in cohort_names) {
  r <- se_list[[cn]]
  st <- flip_list[[cn]]$status[match(r$rsID, flip_list[[cn]]$rsID)]
  keep <- !is.na(st) & st %in% c("same", "flip")
  r <- r[keep, , drop = FALSE]
  st <- st[keep]
  r$beta_h <- ifelse(st == "flip", -r$beta, r$beta)
  harm_list[[cn]] <- r
}

# ------------------------------------------------------------------------------
# 3. Meta-analysis of harmonized lambda.1se estimates
# ------------------------------------------------------------------------------

all_rs <- sort(unique(unlist(lapply(harm_list, function(x) x$rsID))))

meta_one <- function(rs) {
  betas <- ses <- c(); studies <- c()
  for (cn in cohort_names) {
    h <- harm_list[[cn]]
    i <- match(rs, h$rsID)
    if (!is.na(i) && !is.na(h$beta_h[i]) && !is.na(h$SE[i]) && h$SE[i] > 0) {
      betas <- c(betas, h$beta_h[i]); ses <- c(ses, h$SE[i]); studies <- c(studies, cn)
    }
  }
  if (length(betas) < 2) return(NULL)

  w <- 1 / ses^2
  bF <- sum(w * betas) / sum(w)
  seF <- sqrt(1 / sum(w))
  pF <- 2 * pnorm(-abs(bF / seF))
  Q <- sum(w * (betas - bF)^2)
  k <- length(betas)
  I2 <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0

  bR <- seR <- pR <- NA
  if (ok_metafor) {
    fit <- tryCatch(metafor::rma(yi = betas, sei = ses, method = "REML"),
                    error = function(e) NULL)
    if (!is.null(fit)) { bR <- as.numeric(fit$b); seR <- fit$se; pR <- fit$pval }
  }
  if (is.na(bR)) {   # DerSimonian-Laird fallback
    tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
    wr <- 1 / (ses^2 + tau2)
    bR <- sum(wr * betas) / sum(wr); seR <- sqrt(1 / sum(wr))
    pR <- 2 * pnorm(-abs(bR / seR))
  }

  dirs <- vapply(cohort_names, function(cn) {
    if (cn %in% studies) ifelse(betas[studies == cn] > 0, "+", "-") else "."
  }, character(1))

  data.frame(rsID = rs, n_studies_1se = k,
             beta_fixed_1se = bF, SE_fixed_1se = seF, pvalue_fixed_1se = pF,
             HR_fixed_1se = exp(bF),
             HR_lower_1se = exp(bF - 1.96 * seF), HR_upper_1se = exp(bF + 1.96 * seF),
             beta_random_1se = bR, pvalue_random_1se = pR,
             I2_1se = I2, direction_1se = paste0(dirs, collapse = ""),
             stringsAsFactors = FALSE)
}

meta_1se <- do.call(rbind, Filter(Negate(is.null), lapply(all_rs, meta_one)))
if (is.null(meta_1se) || nrow(meta_1se) == 0) {
  stop("No SNP was available in >=2 cohorts after mapping/harmonization.\n",
       "  Most likely cause: a lambda.1se file was paired with the wrong cohort's\n",
       "  lambda.min file, so Variable/Variable_original names did not match.\n",
       "  Check the 'lambda.1se SNP rows mapped to rsIDs' counts printed above:\n",
       "  each cohort should map to roughly its lambda.min SNP count.")
}
meta_1se$pvalue_fixed_FDR_1se <- p.adjust(meta_1se$pvalue_fixed_1se, method = "BH")
cat(sprintf("\nLambda.1se meta-analysis: %d SNPs in >=2 cohorts\n", nrow(meta_1se)))

write.csv(meta_1se, file.path(out_dir, "Meta_analysis_results_all_1se.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Comparison table: lambda.min (published) vs lambda.1se (new)
# ------------------------------------------------------------------------------

keep_min <- c("rsID", "n_studies", "HR_fixed", "HR_lower_fixed", "HR_upper_fixed",
              "pvalue_fixed", "pvalue_fixed_FDR", "I2", "direction")
keep_min <- keep_min[keep_min %in% colnames(meta_min)]
cmp <- merge(meta_min[, keep_min], meta_1se, by = "rsID", all = TRUE)

cmp$beta_fixed_min <- log(cmp$HR_fixed)
cmp$same_direction <- sign(cmp$beta_fixed_min) == sign(cmp$beta_fixed_1se)

sig11 <- c("rs11022690", "rs10987386", "rs6006399", "rs17821105", "rs11227223",
           "rs35684381", "rs17032590", "rs196025", "rs72743477", "rs72811372",
           "rs35797611")
cmp$primary_significant <- cmp$rsID %in% sig11
cmp <- cmp[order(-cmp$primary_significant,
                 ifelse(is.na(cmp$pvalue_fixed), 1, cmp$pvalue_fixed)), ]

write.csv(cmp, file.path(out_dir, "SuppTable_lambda_min_vs_1se.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Reply-letter numbers
# ------------------------------------------------------------------------------

cat("\n==================== REPLY-LETTER NUMBERS ====================\n")
if (!is.null(nz_line)) cat(nz_line, "\n")

shared <- cmp[!is.na(cmp$beta_fixed_min) & !is.na(cmp$beta_fixed_1se), ]
cat(sprintf("SNPs meta-analyzable under both criteria: %d\n", nrow(shared)))
cat(sprintf("Correlation of pooled log-HRs (all shared SNPs): r = %.3f\n",
            cor(shared$beta_fixed_min, shared$beta_fixed_1se)))

s11 <- shared[shared$rsID %in% sig11, ]
if (nrow(s11) > 1) {
  cat(sprintf("Correlation of pooled log-HRs (11 significant SNPs): r = %.3f\n",
              cor(s11$beta_fixed_min, s11$beta_fixed_1se)))
}
cat(sprintf("Direction unchanged: %d of %d significant SNPs (%d of %d overall)\n",
            sum(s11$same_direction), nrow(s11),
            sum(shared$same_direction), nrow(shared)))

r11 <- cmp[cmp$rsID == "rs11022690", ]
if (nrow(r11) == 1 && !is.na(r11$HR_fixed_1se)) {
  cat(sprintf("rs11022690 lambda.min: HR = %.3f (%.3f-%.3f), P = %.3g\n",
              r11$HR_fixed, r11$HR_lower_fixed, r11$HR_upper_fixed, r11$pvalue_fixed))
  cat(sprintf("rs11022690 lambda.1se: HR = %.3f (%.3f-%.3f), P = %.3g, FDR = %.3g\n",
              r11$HR_fixed_1se, r11$HR_lower_1se, r11$HR_upper_1se,
              r11$pvalue_fixed_1se, r11$pvalue_fixed_FDR_1se))
  top_1se <- meta_1se$rsID[which.min(meta_1se$pvalue_fixed_1se)]
  cat(sprintf("Top-ranked SNP under lambda.1se: %s %s\n", top_1se,
              ifelse(top_1se == "rs11022690", "(unchanged)", "(CHANGED - report carefully)")))
}

cat("\nOutputs written to:", out_dir, "\n")
cat("  Meta_analysis_results_all_1se.csv\n")
cat("  SuppTable_lambda_min_vs_1se.csv  (11 significant SNPs sorted first)\n")
cat("Done.\n")
