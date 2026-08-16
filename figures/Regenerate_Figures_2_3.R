rm(list = ls())

# ==============================================================================
# Regenerate Figures 2 and 3 as conventional weighted-square forest plots
# (Reviewer 1, Comment 2; Reviewer 2, Comment 7)
#
# WHAT THIS PRODUCES (in ./figures_revised/ under the current directory)
#   Figure2_Forest_rs11022690.png / .pdf
#     Squares = subgroup-specific debiased HR estimates; square AREA is
#     proportional to the inverse-variance weight; horizontal lines = 95% CIs;
#     diamond = pooled fixed-effect estimate (width = 95% CI).
#   Figure3_Summary_Forest_11SNPs.png / .pdf
#     Squares = pooled fixed-effect HRs for the 11 nominally significant SNPs;
#     square size proportional to the number of contributing subgroups (2 vs 3);
#     color distinguishes Tier 1 (rs11022690) from Tier 2 (the remaining 10
#     nominal signals, per the Methods definition - fixes the R2-7 issue).
#
# INPUT (selected interactively): Meta_analysis_results_all.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- robust file chooser (same as the other revision scripts) ----------------
norm_path <- function(p) {
  if (is.null(p) || length(p) == 0) return(NULL)
  p <- p[1]; if (is.na(p) || !nzchar(p)) return(NULL); p
}
choose_file <- function(msg) {
  cat("\n>>> ", msg, "\n", sep = ""); flush.console()
  path <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() && rstudioapi::hasFun("selectFile")) {
    path <- norm_path(tryCatch(rstudioapi::selectFile(caption = msg), error = function(e) NULL))
  }
  if (is.null(path) && .Platform$OS.type == "windows" && interactive()) {
    path <- norm_path(tryCatch(utils::choose.files(caption = msg, multi = FALSE), error = function(e) NULL))
  }
  if (is.null(path) && interactive()) {
    path <- norm_path(tryCatch(file.choose(), error = function(e) NULL))
  }
  while (is.null(path) || !file.exists(path)) {
    cat("    (No file selected via dialog.)\n")
    path <- trimws(readline(prompt = "    Paste the FULL file path here and press Enter: "))
    path <- gsub('^"|"$', "", path); path <- gsub("\\\\", "/", path)
    if (!nzchar(path) || !file.exists(path)) { cat("    File not found, try again.\n"); path <- NULL }
  }
  cat("    Selected: ", path, "\n", sep = ""); path
}

# ---- load meta results -------------------------------------------------------
repeat {
  meta_file <- choose_file("Select 'Meta_analysis_results_all.csv'")
  meta <- tryCatch(read.csv(meta_file, stringsAsFactors = FALSE), error = function(e) NULL)
  need <- c("rsID", "HR_fixed", "HR_lower_fixed", "HR_upper_fixed", "pvalue_fixed",
            "n_studies", "beta_HSPH", "SE_HSPH", "beta_Onco", "SE_Onco", "beta_MGH", "SE_MGH")
  if (!is.null(meta) && all(need %in% colnames(meta))) break
  cat("    WRONG FILE - this must be Meta_analysis_results_all.csv. Please select again.\n")
}
num_cols <- setdiff(need, "rsID")
for (cc in num_cols) meta[[cc]] <- suppressWarnings(as.numeric(meta[[cc]]))

out_dir <- file.path(getwd(), "figures_revised")
if (!dir.exists(out_dir)) dir.create(out_dir)

# ==============================================================================
# Figure 2: rs11022690 across subgroups + pooled diamond
# ==============================================================================

m <- meta[meta$rsID == "rs11022690", ]
stopifnot(nrow(m) == 1)

sub <- data.frame(
  label = c("HSPH", "OncoArray", "MGH-LPWGS"),
  beta  = c(m$beta_HSPH, m$beta_Onco, m$beta_MGH),
  se    = c(m$SE_HSPH,  m$SE_Onco,  m$SE_MGH),
  stringsAsFactors = FALSE
)
sub <- sub[!is.na(sub$beta), ]
sub$HR    <- exp(sub$beta)
sub$lo    <- exp(sub$beta - 1.96 * sub$se)
sub$hi    <- exp(sub$beta + 1.96 * sub$se)
sub$w     <- 1 / sub$se^2
# square AREA proportional to weight -> side length proportional to sqrt(w)
sub$sq    <- 3 + 5 * sqrt(sub$w / max(sub$w))
sub$y     <- rev(seq_len(nrow(sub))) + 1     # rows from top; leave row 1 for pooled

pooled <- data.frame(HR = m$HR_fixed, lo = m$HR_lower_fixed, hi = m$HR_upper_fixed)
dia_h <- 0.22
diamond <- data.frame(
  x = c(pooled$lo, pooled$HR, pooled$hi, pooled$HR),
  y = c(1, 1 + dia_h, 1, 1 - dia_h)
)

lab_hr <- function(hr, lo, hi) sprintf("%.3f (%.3f\u2013%.3f)", hr, lo, hi)
ann <- rbind(
  data.frame(y = sub$y, label = sub$label, est = lab_hr(sub$HR, sub$lo, sub$hi)),
  data.frame(y = 1, label = "Pooled (fixed effect)", est = lab_hr(pooled$HR, pooled$lo, pooled$hi))
)

xmin <- min(sub$lo, pooled$lo) * 0.97
xmax <- max(sub$hi, pooled$hi) * 1.03

p2 <- ggplot() +
  geom_segment(data = sub, aes(x = lo, xend = hi, y = y, yend = y), linewidth = 0.6) +
  geom_point(data = sub, aes(x = HR, y = y, size = sq), shape = 15, color = "#2166ac") +
  scale_size_identity() +
  geom_polygon(data = diamond, aes(x = x, y = y), fill = "#b2182b", color = "#b2182b") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  scale_x_log10(limits = c(xmin, xmax)) +
  scale_y_continuous(breaks = ann$y, labels = ann$label,
                     sec.axis = sec_axis(~., breaks = ann$y, labels = ann$est),
                     expand = expansion(add = 0.6)) +
  labs(x = "Hazard ratio (95% CI), log scale", y = NULL,
       title = "rs11022690 and overall survival") +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.y.right = element_text(color = "black", size = 10),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggsave(file.path(out_dir, "Figure2_Forest_rs11022690.png"), p2, width = 7.2, height = 3.0, dpi = 600)
ggsave(file.path(out_dir, "Figure2_Forest_rs11022690.pdf"), p2, width = 7.2, height = 3.0)
cat("Figure 2 written.\n")

# ==============================================================================
# Figure 3: pooled estimates for the 11 nominally significant SNPs
# ==============================================================================

sig <- meta[!is.na(meta$pvalue_fixed) & meta$pvalue_fixed < 0.05, ]
sig <- sig[order(sig$pvalue_fixed), ]
cat(sprintf("Nominally significant SNPs found: %d\n", nrow(sig)))

sig$tier <- ifelse(sig$rsID == "rs11022690", "Tier 1", "Tier 2")
sig$y <- rev(seq_len(nrow(sig)))
sig$est_lab <- lab_hr(sig$HR_fixed, sig$HR_lower_fixed, sig$HR_upper_fixed)

p3 <- ggplot(sig) +
  geom_segment(aes(x = HR_lower_fixed, xend = HR_upper_fixed, y = y, yend = y), linewidth = 0.55) +
  geom_point(aes(x = HR_fixed, y = y, size = factor(n_studies), color = tier), shape = 15) +
  scale_size_manual(values = c("2" = 3.2, "3" = 5.2), name = "Contributing subgroups") +
  scale_color_manual(values = c("Tier 1" = "#b2182b", "Tier 2" = "#2166ac"), name = NULL) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  scale_x_log10() +
  scale_y_continuous(breaks = sig$y, labels = sig$rsID,
                     sec.axis = sec_axis(~., breaks = sig$y, labels = sig$est_lab),
                     expand = expansion(add = 0.6)) +
  labs(x = "Pooled fixed-effect hazard ratio (95% CI), log scale", y = NULL,
       title = "Nominally significant SNPs (pooled meta-analysis)") +
  guides(color = guide_legend(order = 1, override.aes = list(size = 4)),
         size = guide_legend(order = 2)) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(color = "black", size = 10.5),
        axis.text.y.right = element_text(color = "black", size = 9.5),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

ggsave(file.path(out_dir, "Figure3_Summary_Forest_11SNPs.png"), p3,
       width = 7.5, height = 0.42 * nrow(sig) + 2.2, dpi = 600)
ggsave(file.path(out_dir, "Figure3_Summary_Forest_11SNPs.pdf"), p3,
       width = 7.5, height = 0.42 * nrow(sig) + 2.2)
cat("Figure 3 written.\n")

# ------------------------------------------------------------------------------
# Captions to paste into the manuscript (consistent with the plots above)
# ------------------------------------------------------------------------------
cat("\n---- Figure 2 caption ----\n")
cat("Figure 2. Association between rs11022690 and overall survival across the three\n")
cat("subgroups. Squares represent subgroup-specific debiased hazard ratio estimates,\n")
cat("with square size proportional to the inverse-variance weight; horizontal lines\n")
cat("indicate 95% confidence intervals; the diamond represents the pooled fixed-effect\n")
cat("estimate, with its width spanning the 95% confidence interval. The dashed line\n")
cat("marks a hazard ratio of 1. Estimates are adjusted for age, sex, smoking, stage,\n")
cat("treatment, and ancestry principal components.\n")
cat("\n---- Figure 3 caption ----\n")
cat("Figure 3. Pooled fixed-effect hazard ratios for the 11 nominally significant SNPs.\n")
cat("Squares show inverse-variance-weighted meta-analytic estimates with 95% confidence\n")
cat("intervals; square size is proportional to the number of contributing subgroups\n")
cat("(2 or 3). rs11022690 (red) met Tier 1 criteria; the remaining 10 SNPs (blue) are\n")
cat("Tier 2 nominal signals. The dashed line marks a hazard ratio of 1.\n")

cat("\nAll figures written to:", out_dir, "\n")
cat("Done.\n")
