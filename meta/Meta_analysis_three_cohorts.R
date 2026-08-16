# ==============================================================================
# Meta分析: 整合 HSPH984, Oncoarray, MGH_LPWGS 三个队列
# Cox Debiased Lasso 结果的固定效应与随机效应 Meta-Analysis 20260211
# ==============================================================================
#
# 数据集概况:
#   - HSPH984:   68 个 SNP
#   - Oncoarray: 45 个 SNP
#   - MGH_LPWGS: 62 个 SNP
#   - 三者共有:   34 个 SNP
#
# ==============================================================================

rm(list = ls())

cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║     Meta-Analysis: Three-Cohort Cox Debiased Lasso Results      ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ==============================================================================
# 1. 加载R包
# ==============================================================================

cat("【步骤 1/10】加载R包...\n")

packages <- c("meta", "metafor", "dplyr", "ggplot2", "ggrepel", "gridExtra")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)

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


  }
}
cat("✓ R包加载完成\n\n")

# ==============================================================================
# 2. 读取数据
# ==============================================================================

cat("【步骤 2/10】读取数据...\n")

# ⚠️ 请根据实际路径修改
hsph <- read.csv(pick_file("Select the debiased results CSV (HSPH)", "RES_HSPH"), stringsAsFactors = FALSE)
onco <- read.csv(pick_file("Select the debiased results CSV (ONCO)", "RES_ONCO"), stringsAsFactors = FALSE)
mgh <- read.csv(pick_file("Select the debiased results CSV (MGH)", "RES_MGH"), stringsAsFactors = FALSE)

cat("  HSPH984:  ", nrow(hsph), "行\n")
cat("  Oncoarray:", nrow(onco), "行\n")
cat("  MGH_LPWGS:", nrow(mgh), "行\n")
cat("✓ 数据读取完成\n\n")

# ==============================================================================
# 3. 提取纯rsID并筛选SNP
# ==============================================================================

cat("【步骤 3/10】标准化SNP标识符...\n")

# 从 Variable 列 (如 rs12345_A) 提取纯 rsID (如 rs12345)
extract_rsID <- function(x) gsub("_[A-Za-z0-9]+$", "", x)

# 添加纯rsID列
hsph$rsID_clean <- extract_rsID(hsph$Variable)
onco$rsID_clean <- extract_rsID(onco$Variable)
mgh$rsID_clean  <- extract_rsID(mgh$Variable)

# 筛选SNP行 (以rs开头)
hsph_snp <- hsph[grepl("^rs\\d+", hsph$rsID_clean), ]
onco_snp <- onco[grepl("^rs\\d+", onco$rsID_clean), ]
mgh_snp  <- mgh[grepl("^rs\\d+", mgh$rsID_clean), ]

cat("  HSPH984 SNP:  ", nrow(hsph_snp), "\n")
cat("  Oncoarray SNP:", nrow(onco_snp), "\n")
cat("  MGH_LPWGS SNP:", nrow(mgh_snp), "\n")
cat("✓ SNP筛选完成\n\n")

# ==============================================================================
# 4. 计算SNP重叠
# ==============================================================================

cat("【步骤 4/10】分析SNP重叠情况...\n")

rsid_hsph <- unique(hsph_snp$rsID_clean)
rsid_onco <- unique(onco_snp$rsID_clean)
rsid_mgh  <- unique(mgh_snp$rsID_clean)

# 两两交集
common_hsph_onco <- intersect(rsid_hsph, rsid_onco)
common_hsph_mgh  <- intersect(rsid_hsph, rsid_mgh)
common_onco_mgh  <- intersect(rsid_onco, rsid_mgh)

# 三者交集
common_all <- Reduce(intersect, list(rsid_hsph, rsid_onco, rsid_mgh))

# 至少两个有的SNP (用于部分meta)
common_any2 <- unique(c(common_hsph_onco, common_hsph_mgh, common_onco_mgh))

cat("  ┌─────────────────────────────────────┐\n")
cat("  │ HSPH ∩ Oncoarray:  ", sprintf("%3d", length(common_hsph_onco)), " 个 SNP    │\n")
cat("  │ HSPH ∩ MGH:        ", sprintf("%3d", length(common_hsph_mgh)), " 个 SNP    │\n")
cat("  │ Oncoarray ∩ MGH:   ", sprintf("%3d", length(common_onco_mgh)), " 个 SNP    │\n")
cat("  │ 三者共有:          ", sprintf("%3d", length(common_all)), " 个 SNP    │\n")
cat("  │ 至少两个有:        ", sprintf("%3d", length(common_any2)), " 个 SNP    │\n")
cat("  └─────────────────────────────────────┘\n")
cat("✓ 重叠分析完成\n\n")

# ==============================================================================
# 5. 执行Meta分析 (核心部分)
# ==============================================================================

cat("【步骤 5/10】执行Meta分析...\n")

meta_results <- list()
snp_count <- 0

for (snp in common_any2) {
  
  # 获取各数据集中该SNP的数据
  h <- hsph_snp[hsph_snp$rsID_clean == snp, ]
  o <- onco_snp[onco_snp$rsID_clean == snp, ]
  m <- mgh_snp[mgh_snp$rsID_clean == snp, ]
  
  # 收集有效数据
  betas   <- c()
  ses     <- c()
  studies <- c()
  
  if (nrow(h) > 0 && !is.na(h$beta[1]) && !is.na(h$SE[1]) && h$SE[1] > 0) {
    betas   <- c(betas, h$beta[1])
    ses     <- c(ses, h$SE[1])
    studies <- c(studies, "HSPH984")
  }
  
  if (nrow(o) > 0 && !is.na(o$beta[1]) && !is.na(o$SE[1]) && o$SE[1] > 0) {
    betas   <- c(betas, o$beta[1])
    ses     <- c(ses, o$SE[1])
    studies <- c(studies, "Oncoarray")
  }
  
  if (nrow(m) > 0 && !is.na(m$beta[1]) && !is.na(m$SE[1]) && m$SE[1] > 0) {
    betas   <- c(betas, m$beta[1])
    ses     <- c(ses, m$SE[1])
    studies <- c(studies, "MGH_LPWGS")
  }
  
  n_studies <- length(studies)
  
  # 至少2个研究才能做Meta
  if (n_studies >= 2) {
    
    tryCatch({
      # 使用 meta 包的 metagen 函数
      meta_obj <- metagen(
        TE      = betas,
        seTE    = ses,
        studlab = studies,
        sm      = "HR",
        method.tau = "REML",
        hakn    = FALSE
      )
      
      # 构建方向字符串
      dir_h <- ifelse("HSPH984" %in% studies, 
                      ifelse(betas[studies == "HSPH984"] > 0, "+", "-"), ".")
      dir_o <- ifelse("Oncoarray" %in% studies,
                      ifelse(betas[studies == "Oncoarray"] > 0, "+", "-"), ".")
      dir_m <- ifelse("MGH_LPWGS" %in% studies,
                      ifelse(betas[studies == "MGH_LPWGS"] > 0, "+", "-"), ".")
      direction <- paste0(dir_h, dir_o, dir_m)
      
      # 存储结果
      meta_results[[snp]] <- data.frame(
        rsID = snp,
        n_studies = n_studies,
        
        # 固定效应
        beta_fixed = meta_obj$TE.fixed,
        SE_fixed = meta_obj$seTE.fixed,
        pvalue_fixed = meta_obj$pval.fixed,
        HR_fixed = exp(meta_obj$TE.fixed),
        HR_lower_fixed = exp(meta_obj$lower.fixed),
        HR_upper_fixed = exp(meta_obj$upper.fixed),
        
        # 随机效应
        beta_random = meta_obj$TE.random,
        SE_random = meta_obj$seTE.random,
        pvalue_random = meta_obj$pval.random,
        HR_random = exp(meta_obj$TE.random),
        HR_lower_random = exp(meta_obj$lower.random),
        HR_upper_random = exp(meta_obj$upper.random),
        
        # 异质性
        tau2 = meta_obj$tau2,
        I2 = meta_obj$I2 * 100,
        Q = meta_obj$Q,
        Q_pvalue = meta_obj$pval.Q,
        
        # 各研究beta
        beta_HSPH = ifelse("HSPH984" %in% studies, betas[studies == "HSPH984"], NA),
        beta_Onco = ifelse("Oncoarray" %in% studies, betas[studies == "Oncoarray"], NA),
        beta_MGH  = ifelse("MGH_LPWGS" %in% studies, betas[studies == "MGH_LPWGS"], NA),
        
        # 各研究SE
        SE_HSPH = ifelse("HSPH984" %in% studies, ses[studies == "HSPH984"], NA),
        SE_Onco = ifelse("Oncoarray" %in% studies, ses[studies == "Oncoarray"], NA),
        SE_MGH  = ifelse("MGH_LPWGS" %in% studies, ses[studies == "MGH_LPWGS"], NA),
        
        # 方向
        direction = direction,
        
        stringsAsFactors = FALSE
      )
      
      snp_count <- snp_count + 1
      
    }, error = function(e) {
      # 静默处理错误
    })
  }
}

# 合并结果
meta_df <- do.call(rbind, meta_results)
rownames(meta_df) <- NULL

cat("  成功分析SNP数:", nrow(meta_df), "\n")
cat("✓ Meta分析完成\n\n")

# ==============================================================================
# 6. 添加统计校正和注释
# ==============================================================================

cat("【步骤 6/10】统计校正与注释...\n")

# FDR校正
meta_df$pvalue_fixed_FDR  <- p.adjust(meta_df$pvalue_fixed, method = "BH")
meta_df$pvalue_random_FDR <- p.adjust(meta_df$pvalue_random, method = "BH")

# 方向一致性判断 (不含 +- 或 -+ 组合)
meta_df$direction_consistent <- !grepl("\\+[^.]*-|-[^.]*\\+", meta_df$direction)

# 添加CHR和BP (从HSPH获取)
meta_df$CHR <- NA
meta_df$BP  <- NA

for (i in seq_len(nrow(meta_df))) {
  snp <- meta_df$rsID[i]
  match_row <- hsph_snp[hsph_snp$rsID_clean == snp, ]
  if (nrow(match_row) > 0 && "CHR" %in% names(match_row)) {
    meta_df$CHR[i] <- match_row$CHR[1]
    meta_df$BP[i]  <- match_row$BP[1]
  }
}

# 按固定效应p值排序
meta_df <- meta_df[order(meta_df$pvalue_fixed), ]

cat("  方向一致SNP:", sum(meta_df$direction_consistent), "/", nrow(meta_df),
    "(", round(mean(meta_df$direction_consistent) * 100, 1), "%)\n")
cat("✓ 统计校正完成\n\n")

# ==============================================================================
# 7. 保存结果
# ==============================================================================

cat("【步骤 7/10】保存结果...\n")

# 全部结果
write.csv(meta_df, "Meta_analysis_results_all.csv", row.names = FALSE)
cat("  ✓ Meta_analysis_results_all.csv\n")

# 显著结果
sig_df <- meta_df[meta_df$pvalue_fixed < 0.05 | meta_df$pvalue_random < 0.05, ]
write.csv(sig_df, "Meta_analysis_results_significant.csv", row.names = FALSE)
cat("  ✓ Meta_analysis_results_significant.csv (", nrow(sig_df), "个SNP)\n\n")

# ==============================================================================
# 8. 森林图 (Forest Plot)
# ==============================================================================

cat("【步骤 8/10】生成森林图...\n")

# 8a. 汇总森林图 (显著SNP)
if (nrow(sig_df) > 0) {
  
  forest_data <- sig_df[, c("rsID", "HR_random", "HR_lower_random", "HR_upper_random", 
                            "pvalue_random", "I2", "direction")]
  names(forest_data) <- c("SNP", "HR", "HR_lower", "HR_upper", "pvalue", "I2", "direction")
  forest_data$SNP <- factor(forest_data$SNP, levels = rev(forest_data$SNP))
  
  p_forest <- ggplot(forest_data, aes(x = HR, y = SNP)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.3, color = "steelblue") +
    geom_point(aes(color = pvalue < 0.01), size = 3, shape = 18) +
    scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"), guide = "none") +
    scale_x_log10() +
    labs(
      title = "Meta-Analysis Forest Plot (Random Effects)",
      subtitle = paste0("SNPs with P < 0.05 (n = ", nrow(forest_data), ")"),
      x = "Hazard Ratio (log scale)", y = ""
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  
  ggsave("Meta_forest_plot_significant.png", p_forest, 
         width = 10, height = max(5, nrow(forest_data) * 0.5), dpi = 200)
  cat("  ✓ Meta_forest_plot_significant.png\n")
}

# 8b. 单个SNP详细森林图 (Top 5)
top5 <- head(meta_df$rsID, 5)

for (snp in top5) {
  
  h <- hsph_snp[hsph_snp$rsID_clean == snp, ]
  o <- onco_snp[onco_snp$rsID_clean == snp, ]
  m <- mgh_snp[mgh_snp$rsID_clean == snp, ]
  meta_row <- meta_df[meta_df$rsID == snp, ]
  
  plot_df <- data.frame(Study = character(), beta = numeric(), SE = numeric())
  
  if (nrow(h) > 0 && !is.na(h$beta[1])) 
    plot_df <- rbind(plot_df, data.frame(Study = "HSPH984", beta = h$beta[1], SE = h$SE[1]))
  if (nrow(o) > 0 && !is.na(o$beta[1])) 
    plot_df <- rbind(plot_df, data.frame(Study = "Oncoarray", beta = o$beta[1], SE = o$SE[1]))
  if (nrow(m) > 0 && !is.na(m$beta[1])) 
    plot_df <- rbind(plot_df, data.frame(Study = "MGH_LPWGS", beta = m$beta[1], SE = m$SE[1]))
  
  # 添加Meta结果
  plot_df <- rbind(plot_df, data.frame(
    Study = "Meta (Random)", 
    beta = meta_row$beta_random[1], 
    SE = meta_row$SE_random[1]
  ))
  
  plot_df$HR <- exp(plot_df$beta)
  plot_df$HR_lower <- exp(plot_df$beta - 1.96 * plot_df$SE)
  plot_df$HR_upper <- exp(plot_df$beta + 1.96 * plot_df$SE)
  plot_df$Study <- factor(plot_df$Study, 
                          levels = rev(c("HSPH984", "Oncoarray", "MGH_LPWGS", "Meta (Random)")))
  
  p <- ggplot(plot_df, aes(x = HR, y = Study)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
    geom_point(aes(color = Study == "Meta (Random)"), size = 4) +
    scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"), guide = "none") +
    scale_x_log10() +
    labs(
      title = paste0("Forest Plot: ", snp),
      subtitle = paste0("I² = ", round(meta_row$I2[1], 1), 
                        "%, P(heterogeneity) = ", format(meta_row$Q_pvalue[1], digits = 2)),
      x = "Hazard Ratio (log scale)", y = ""
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(paste0("Forest_", snp, ".png"), p, width = 8, height = 4, dpi = 150)
}
cat("  ✓ Forest_[rsID].png (Top 5)\n")

# ==============================================================================
# 9. 漏斗图 (Funnel Plot) - 发表偏倚检测
# ==============================================================================

cat("\n【步骤 9/10】生成漏斗图...\n")

funnel_data <- data.frame()

for (snp in common_all) {
  h <- hsph_snp[hsph_snp$rsID_clean == snp, ]
  o <- onco_snp[onco_snp$rsID_clean == snp, ]
  m <- mgh_snp[mgh_snp$rsID_clean == snp, ]
  
  if (nrow(h) > 0) funnel_data <- rbind(funnel_data, 
    data.frame(beta = h$beta[1], SE = h$SE[1], Study = "HSPH984", SNP = snp))
  if (nrow(o) > 0) funnel_data <- rbind(funnel_data, 
    data.frame(beta = o$beta[1], SE = o$SE[1], Study = "Oncoarray", SNP = snp))
  if (nrow(m) > 0) funnel_data <- rbind(funnel_data, 
    data.frame(beta = m$beta[1], SE = m$SE[1], Study = "MGH_LPWGS", SNP = snp))
}

p_funnel <- ggplot(funnel_data, aes(x = beta, y = SE)) +
  geom_point(aes(color = Study), alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_reverse() +
  scale_color_manual(values = c("HSPH984" = "#4477AA", "Oncoarray" = "#EE6677", "MGH_LPWGS" = "#228833")) +
  labs(
    title = "Funnel Plot - Publication Bias Assessment",
    subtitle = paste0("SNPs present in all 3 studies (n = ", length(common_all), ")"),
    x = "Effect Size (beta)", y = "Standard Error"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

ggsave("Meta_funnel_plot.png", p_funnel, width = 10, height = 8, dpi = 200)
cat("  ✓ Meta_funnel_plot.png\n")

# ==============================================================================
# 10. 曼哈顿图 (Manhattan Plot)
# ==============================================================================

cat("\n【步骤 10/10】生成曼哈顿图...\n")

man_data <- meta_df[!is.na(meta_df$CHR) & meta_df$pvalue_fixed > 0, ]
man_data$log10P <- -log10(man_data$pvalue_fixed)

if (nrow(man_data) > 0) {
  
  bonf <- 0.05 / nrow(man_data)
  
  p_man <- ggplot(man_data, aes(x = CHR, y = log10P)) +
    geom_point(aes(color = factor(CHR %% 2)), size = 3, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange") +
    geom_hline(yintercept = -log10(bonf), linetype = "dashed", color = "red") +
    geom_text_repel(
      data = man_data[man_data$pvalue_fixed < 0.05, ],
      aes(label = rsID), size = 3, fontface = "bold",
      box.padding = 0.5, max.overlaps = 20
    ) +
    scale_color_manual(values = c("#4477AA", "#EE6677"), guide = "none") +
    scale_x_continuous(breaks = sort(unique(man_data$CHR))) +
    labs(
      title = "Meta-Analysis Manhattan Plot (Fixed Effects)",
      x = "Chromosome", y = expression(-log[10](P))
    ) +
    annotate("text", x = max(man_data$CHR), y = -log10(0.05) + 0.2,
             label = "P=0.05", color = "orange", hjust = 1, size = 3) +
    annotate("text", x = max(man_data$CHR), y = -log10(bonf) + 0.2,
             label = "Bonferroni", color = "red", hjust = 1, size = 3) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("Meta_manhattan_plot.png", p_man, width = 12, height = 6, dpi = 200)
  cat("  ✓ Meta_manhattan_plot.png\n")
}

# ==============================================================================
# 结果汇总
# ==============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                        结 果 汇 总                              ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# 异质性分布
cat("【异质性分析 (I²)】\n")
cat("  I² = 0%:        ", sum(meta_df$I2 == 0), "个SNP (无异质性)\n")
cat("  I² < 25%:       ", sum(meta_df$I2 > 0 & meta_df$I2 < 25), "个SNP (低)\n")
cat("  I² 25-50%:      ", sum(meta_df$I2 >= 25 & meta_df$I2 < 50), "个SNP (中)\n")
cat("  I² 50-75%:      ", sum(meta_df$I2 >= 50 & meta_df$I2 < 75), "个SNP (较高)\n")
cat("  I² ≥ 75%:       ", sum(meta_df$I2 >= 75), "个SNP (高)\n\n")

# 显著性检验
bonf <- 0.05 / nrow(meta_df)
cat("【显著性检验】\n")
cat("  P < 0.05 (Fixed):  ", sum(meta_df$pvalue_fixed < 0.05), "个SNP\n")
cat("  P < 0.05 (Random): ", sum(meta_df$pvalue_random < 0.05), "个SNP\n")
cat("  Bonferroni (", format(bonf, digits = 2), "):", sum(meta_df$pvalue_fixed < bonf), "个SNP\n")
cat("  FDR < 0.05:        ", sum(meta_df$pvalue_fixed_FDR < 0.05), "个SNP\n\n")

# Top 10 结果
cat("【Top 10 SNPs (按Fixed P值)】\n")
cat("─────────────────────────────────────────────────────────────────────\n")
top10 <- head(meta_df, 10)
print(top10[, c("rsID", "n_studies", "HR_fixed", "pvalue_fixed", "I2", "direction")])

# 显著SNP详情
if (nrow(sig_df) > 0) {
  cat("\n【显著SNP详情 (P < 0.05)】\n")
  cat("─────────────────────────────────────────────────────────────────────\n")
  print(sig_df[, c("rsID", "CHR", "n_studies", "HR_fixed", "HR_lower_fixed", 
                   "HR_upper_fixed", "pvalue_fixed", "I2", "direction")])
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                     输 出 文 件 列 表                           ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
cat("  1. Meta_analysis_results_all.csv          - 全部Meta分析结果\n")
cat("  2. Meta_analysis_results_significant.csv  - 显著SNP结果\n")
cat("  3. Meta_forest_plot_significant.png       - 显著SNP森林图\n")
cat("  4. Forest_[rsID].png                      - 单SNP森林图 (Top 5)\n")
cat("  5. Meta_funnel_plot.png                   - 漏斗图\n")
cat("  6. Meta_manhattan_plot.png                - 曼哈顿图\n")
cat("\n")
cat("══════════════════════════════════════════════════════════════════\n")
cat("                    Meta分析完成！                                \n")
cat("══════════════════════════════════════════════════════════════════\n")
