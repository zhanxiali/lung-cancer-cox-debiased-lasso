########################################################################
#  Treatment Sensitivity: Meta分析 + 与主分析比较
#  修复版: 加入sign correction（与主分析保持一致的allele编码）
#
#  问题: 各队列per-cohort CSV中的beta使用的是各自平台的allele编码
#        主分析的Meta脚本在汇总时做了sign correction（翻转28个SNP的beta）
#        敏感性分析必须做同样的翻转，否则会出现虚假的方向不一致
#
#  解决: 使用 sign_flip_table.csv（从主分析Meta结果反推的翻转表）
#        对Oncoarray 16个 + MGH 17个SNP进行beta翻转
########################################################################

library(dplyr)
library(meta)
library(metafor)

cat("请选择输出目录...\n")
output_dir <- choose.dir(caption = "选择保存结果的目录")

# ==================================================================
#  Step 1: 读取 sign_flip_table
# ==================================================================

cat("\n=== 读取 Sign-Flip 校正表 ===\n")
cat("请选择 sign_flip_table.csv...\n")
flip_table <- read.csv(file.choose(), stringsAsFactors = FALSE)
cat("  Oncoarray需flip:", sum(flip_table$flip_Onco == 1), "个SNP\n")
cat("  MGH需flip:", sum(flip_table$flip_MGH == 1), "个SNP\n")

# ==================================================================
#  Step 2: 读取敏感性分析结果
# ==================================================================

cat("\n=== 读取敏感性分析结果 ===\n\n")

cat("[HSPH] 请选择 HSPH984_NO_TREATMENT_results_all.csv...\n")
hsph <- read.csv(file.choose(), stringsAsFactors = FALSE)
cat("  HSPH:", nrow(hsph), "个变量\n")

cat("[Oncoarray] 请选择 Oncoarray_NO_TREATMENT_results_all.csv...\n")
onco <- read.csv(file.choose(), stringsAsFactors = FALSE)
cat("  Oncoarray:", nrow(onco), "个变量\n")

cat("[MGH] 请选择 MGH_LPWGS_NO_TREATMENT_results_all.csv...\n")
mgh <- read.csv(file.choose(), stringsAsFactors = FALSE)
cat("  MGH:", nrow(mgh), "个变量\n")

# ==================================================================
#  Step 3: 提取SNP行并标准化rsID
# ==================================================================

cat("\n=== 标准化 SNP 名称 ===\n")

extract_snps <- function(df, label) {
  if ("rsID" %in% names(df)) {
    df$rsID_clean <- gsub("_[ATCG]+$", "", df$rsID)
  } else {
    df$rsID_clean <- gsub("_[ATCG]+$", "", df$Variable)
  }
  snp_df <- df[grepl("^rs", df$rsID_clean), ]
  cat(paste0("  ", label, " SNPs: ", nrow(snp_df), "\n"))
  return(snp_df)
}

hsph_snp <- extract_snps(hsph, "HSPH")
onco_snp <- extract_snps(onco, "Oncoarray")
mgh_snp  <- extract_snps(mgh,  "MGH")

# ==================================================================
#  Step 4: 应用 Sign Correction
#  对Oncoarray和MGH中需要flip的SNP，将beta取反
# ==================================================================

cat("\n=== 应用 Sign Correction ===\n")

# Oncoarray flip
onco_flip_rsids <- flip_table$rsID[flip_table$flip_Onco == 1]
n_flipped_onco <- 0
for (i in seq_len(nrow(onco_snp))) {
  if (onco_snp$rsID_clean[i] %in% onco_flip_rsids) {
    onco_snp$beta[i] <- -onco_snp$beta[i]
    n_flipped_onco <- n_flipped_onco + 1
  }
}
cat("  Oncoarray: 翻转了", n_flipped_onco, "个SNP的beta\n")

# MGH flip
mgh_flip_rsids <- flip_table$rsID[flip_table$flip_MGH == 1]
n_flipped_mgh <- 0
for (i in seq_len(nrow(mgh_snp))) {
  if (mgh_snp$rsID_clean[i] %in% mgh_flip_rsids) {
    mgh_snp$beta[i] <- -mgh_snp$beta[i]
    n_flipped_mgh <- n_flipped_mgh + 1
  }
}
cat("  MGH: 翻转了", n_flipped_mgh, "个SNP的beta\n")

# HSPH不需要翻转（参考队列）

# ==================================================================
#  Step 5: 固定效应Meta分析
# ==================================================================

cat("\n=== Meta分析 ===\n\n")

all_rsids <- unique(c(hsph_snp$rsID_clean, onco_snp$rsID_clean, mgh_snp$rsID_clean))
cat("总unique SNPs:", length(all_rsids), "\n")

meta_results <- data.frame()

for (rs in all_rsids) {
  betas <- c()
  ses <- c()
  labels <- c()
  
  row_h <- hsph_snp[hsph_snp$rsID_clean == rs, ]
  if (nrow(row_h) == 1 && !is.na(row_h$beta) && !is.na(row_h$SE) && row_h$SE > 0) {
    betas <- c(betas, row_h$beta); ses <- c(ses, row_h$SE); labels <- c(labels, "HSPH")
  }
  
  row_o <- onco_snp[onco_snp$rsID_clean == rs, ]
  if (nrow(row_o) == 1 && !is.na(row_o$beta) && !is.na(row_o$SE) && row_o$SE > 0) {
    betas <- c(betas, row_o$beta); ses <- c(ses, row_o$SE); labels <- c(labels, "Onco")
  }
  
  row_m <- mgh_snp[mgh_snp$rsID_clean == rs, ]
  if (nrow(row_m) == 1 && !is.na(row_m$beta) && !is.na(row_m$SE) && row_m$SE > 0) {
    betas <- c(betas, row_m$beta); ses <- c(ses, row_m$SE); labels <- c(labels, "MGH")
  }
  
  if (length(betas) >= 2) {
    w <- 1 / ses^2
    beta_fixed <- sum(w * betas) / sum(w)
    se_fixed <- sqrt(1 / sum(w))
    z <- beta_fixed / se_fixed
    p_fixed <- 2 * pnorm(-abs(z))
    
    Q <- sum(w * (betas - beta_fixed)^2)
    df_q <- length(betas) - 1
    Q_p <- pchisq(Q, df_q, lower.tail = FALSE)
    I2 <- max(0, (Q - df_q) / Q * 100)
    
    dirs <- ifelse(betas < 0, "-", "+")
    direction <- paste0(
      ifelse("HSPH" %in% labels, dirs[labels == "HSPH"], "."),
      ifelse("Onco" %in% labels, dirs[labels == "Onco"], "."),
      ifelse("MGH"  %in% labels, dirs[labels == "MGH"],  ".")
    )
    
    meta_results <- rbind(meta_results, data.frame(
      rsID = rs, n_studies = length(betas),
      beta_fixed_sens = round(beta_fixed, 6),
      SE_fixed_sens = round(se_fixed, 6),
      pvalue_fixed_sens = p_fixed,
      HR_fixed_sens = round(exp(beta_fixed), 4),
      HR_lower_sens = round(exp(beta_fixed - 1.96 * se_fixed), 4),
      HR_upper_sens = round(exp(beta_fixed + 1.96 * se_fixed), 4),
      I2_sens = round(I2, 1),
      direction_sens = direction,
      beta_HSPH = ifelse("HSPH" %in% labels, round(betas[labels=="HSPH"], 6), NA),
      beta_Onco = ifelse("Onco" %in% labels, round(betas[labels=="Onco"], 6), NA),
      beta_MGH  = ifelse("MGH"  %in% labels, round(betas[labels=="MGH"],  6), NA),
      stringsAsFactors = FALSE
    ))
  }
}

meta_results$pvalue_FDR_sens <- p.adjust(meta_results$pvalue_fixed_sens, method = "BH")
meta_results <- meta_results[order(meta_results$pvalue_fixed_sens), ]

cat("Meta分析完成:", nrow(meta_results), "个SNP\n")
cat("  P<0.05:", sum(meta_results$pvalue_fixed_sens < 0.05), "\n")
cat("  FDR<0.05:", sum(meta_results$pvalue_FDR_sens < 0.05, na.rm = TRUE), "\n")

f1 <- file.path(output_dir, "Meta_sensitivity_NO_TREATMENT_all.csv")
write.csv(meta_results, f1, row.names = FALSE)
cat("保存:", f1, "\n")

# ==================================================================
#  Step 6: 与主分析结果比较
# ==================================================================

cat("\n=== 与主分析比较 ===\n\n")
cat("请选择 Meta_analysis_results_all.csv (主分析结果)...\n")
main <- read.csv(file.choose(), stringsAsFactors = FALSE)

comparison <- merge(
  main[, c("rsID", "n_studies", "beta_fixed", "SE_fixed", "pvalue_fixed",
           "HR_fixed", "HR_lower_fixed", "HR_upper_fixed", "I2", "direction")],
  meta_results[, c("rsID", "n_studies", "beta_fixed_sens", "SE_fixed_sens",
                    "pvalue_fixed_sens", "HR_fixed_sens", "HR_lower_sens",
                    "HR_upper_sens", "I2_sens", "direction_sens")],
  by = "rsID", suffixes = c("_main", "_sens")
)

comparison <- comparison %>%
  mutate(
    beta_diff = beta_fixed_sens - beta_fixed,
    HR_diff = HR_fixed_sens - HR_fixed,
    direction_same = sign(beta_fixed) == sign(beta_fixed_sens),
    sens_in_main_CI = HR_fixed_sens >= HR_lower_fixed & HR_fixed_sens <= HR_upper_fixed
  ) %>%
  arrange(pvalue_fixed)

f2 <- file.path(output_dir, "SuppTable_Sy_Treatment_Sensitivity_Comparison.csv")
write.csv(comparison, f2, row.names = FALSE)
cat("比较表保存:", f2, "\n")

# ==================================================================
#  Step 7: 汇总报告
# ==================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("  Treatment Sensitivity Analysis 汇总\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

n_same <- sum(comparison$direction_same, na.rm = TRUE)
n_total <- nrow(comparison)
n_in_ci <- sum(comparison$sens_in_main_CI, na.rm = TRUE)

cat(sprintf("  方向一致性: %d/%d (%.1f%%)\n", n_same, n_total, n_same/n_total*100))
cat(sprintf("  敏感性HR在主分析CI内: %d/%d (%.1f%%)\n", n_in_ci, n_total, n_in_ci/n_total*100))

cat("\n  --- Top 11 SNPs 详细比较 ---\n")
top11 <- c("rs11022690", "rs10987386", "rs6006399", "rs17821105",
           "rs11227223", "rs35684381", "rs17032590", "rs196025",
           "rs72743477", "rs35797611", "rs72811372")

for (rs in top11) {
  row <- comparison[comparison$rsID == rs, ]
  if (nrow(row) == 1) {
    cat(sprintf("  %s: HR %.3f -> %.3f (diff=%.3f), dir=%s, inCI=%s, P: %s -> %s\n",
                row$rsID, row$HR_fixed, row$HR_fixed_sens, row$HR_diff,
                ifelse(row$direction_same, "same", "CHANGED"),
                ifelse(row$sens_in_main_CI, "yes", "NO"),
                format(row$pvalue_fixed, digits = 3, scientific = TRUE),
                format(row$pvalue_fixed_sens, digits = 3, scientific = TRUE)))
  }
}

cat("\n  --- rs11022690 (Tier 1 signal) ---\n")
r11 <- comparison[comparison$rsID == "rs11022690", ]
if (nrow(r11) == 1) {
  cat(sprintf("    主分析:   HR = %.3f (%.3f-%.3f), P = %s\n",
              r11$HR_fixed, r11$HR_lower_fixed, r11$HR_upper_fixed,
              format(r11$pvalue_fixed, scientific = TRUE)))
  cat(sprintf("    敏感性:   HR = %.3f (%.3f-%.3f), P = %s\n",
              r11$HR_fixed_sens, r11$HR_lower_sens, r11$HR_upper_sens,
              format(r11$pvalue_fixed_sens, scientific = TRUE)))
  cat(sprintf("    方向一致: %s\n", ifelse(r11$direction_same, "YES", "NO")))
  cat(sprintf("    在主CI内: %s\n", ifelse(r11$sens_in_main_CI, "YES", "NO")))
}

cat("\n\n  ========== 论文建议措辞 ==========\n\n")
cat(sprintf('  "Sensitivity analysis excluding treatment covariates\n'))
cat(sprintf('   (surgery, chemotherapy, radiation) yielded directionally\n'))
cat(sprintf('   consistent results for %d of %d meta-analyzed SNPs (%.0f%%).\n',
            n_same, n_total, n_same/n_total*100))
if (nrow(r11) == 1) {
  cat(sprintf('   For the top-ranked signal rs11022690, the hazard ratio\n'))
  cat(sprintf('   was %.3f (P = %s) in the sensitivity analysis,\n',
              r11$HR_fixed_sens, format(r11$pvalue_fixed_sens, digits=3, scientific=TRUE)))
}
cat(sprintf('   compared with 0.927 (P = 7.79 x 10^-4) in the primary\n'))
cat(sprintf('   analysis (Supplementary Table Sy), supporting the\n'))
cat(sprintf('   robustness of the primary findings to treatment\n'))
cat(sprintf('   covariate specification."\n'))

cat("\n完成!\n")
