rm(list=ls())

# ==============================================================================
# MGH LPWGS - Cox Debiased Lasso 分析 (终极完整版)
# 
# 更新内容:
# 1. 输入: 读取 PLINK 对齐后的 _aligned.raw 数据
# 2. 模型: 协变量中已移除 immunotherapy
# 3. 统计: 添加 FDR (Benjamini-Hochberg) 校正
# 4. 可视化: 全部升级为 ggplot2 风格 (曼哈顿图、森林图、火山图)
# 5. 输出: 打印结果中包含 FDR P值，确保与 Oncoarray 格式统一
# ==============================================================================

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
# library(qqman) # 已弃用旧版画图包

# 检查并自动安装 ggrepel (如果未安装)
if (!require(ggrepel, quietly = TRUE)) {
  install.packages("ggrepel", repos = "https://cloud.r-project.org")
  library(ggrepel)
} else {
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


}

# 加载 C++ 函数（需要 sim_univLib.cpp 文件）
sourceCpp(Sys.getenv("SIM_UNIVLIB", unset = "../core/sim_univLib.cpp"))

# ==============================================================================
# 1. 设置工作目录
# ==============================================================================


cat("当前工作目录:", getwd(), "\n\n")

# ==============================================================================
# 2. 数据准备
# ==============================================================================

cat("=== 第一部分：数据准备 ===\n\n")

# 1. 读取 SNP 数据 (使用对齐后的文件)
raw_file <- pick_file("Select the aligned PLINK .raw genotype file (MGH)", "GENO_RAW")
if(!file.exists(raw_file)) stop("错误：找不到输入文件 ", raw_file)

df <- read.table(raw_file, header = TRUE, stringsAsFactors = FALSE)
cat("SNP 数据 (Aligned): ", nrow(df), "例,", ncol(df), "列\n")

# 2. 读取临床数据
clin_file <- pick_file("Select the cleaned clinical CSV (MGH)", "CLIN_CSV")
df2 <- read.csv(clin_file, header = TRUE, stringsAsFactors = FALSE)
cat("临床数据: ", nrow(df2), "例,", ncol(df2), "列\n")

# 3. 提取 SNP 列名
snp_cols <- colnames(df)[7:ncol(df)]
cat("SNP 数量: ", length(snp_cols), "\n")
cat("SNP 示例: ", paste(head(snp_cols, 3), collapse=", "), "\n\n")

# 4. 合并数据
df_snp <- df[, c("IID", snp_cols)]
df_all <- merge(df_snp, df2, by = "IID")
cat("合并后数据: ", nrow(df_all), "例\n\n")

# ==============================================================================
# 3. 变量处理
# ==============================================================================

cat("=== 第二部分：变量处理 ===\n\n")

# 将分类变量转为因子
vars_to_factor <- c("SEX", "smoksort", "early_late", "RADS", "chemotx", "surgery")
for(v in vars_to_factor) {
  if(v %in% names(df_all)) df_all[[v]] <- as.factor(df_all[[v]])
}

# 过滤缺失生存数据的样本
df_all <- df_all %>% filter(!is.na(OS_month), !is.na(DEAD))
cat("过滤后样本数: ", nrow(df_all), "\n")

# 缺失值填补函数
impute_mode_mean <- function(d) {
  d <- as.data.frame(d)
  
  # 字符转因子
  char_cols <- vapply(d, is.character, logical(1))
  d[char_cols] <- lapply(d[char_cols], as.factor)
  
  # 因子填众数
  fac_cols <- vapply(d, is.factor, logical(1))
  d[fac_cols] <- lapply(d[fac_cols], function(x) {
    tab <- table(x, useNA = "no")
    if (length(tab) == 0) return(x)
    mode_val <- names(tab)[which.max(tab)]
    if (!mode_val %in% levels(x)) x <- factor(x, levels = c(levels(x), mode_val))
    x[is.na(x)] <- mode_val
    x
  })
  
  # 数值填均值
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

# SNP 列确保为数值型
df_imputed[snp_cols] <- lapply(df_imputed[snp_cols], as.character)
df_imputed[snp_cols] <- lapply(df_imputed[snp_cols], as.numeric)

# ==============================================================================
# 4. 构建设计矩阵
# ==============================================================================

cat("\n=== 第三部分：构建设计矩阵 ===\n\n")

# 定义模型变量 (已移除 immunotherapy)
final_vars <- c(snp_cols, 
                "SEX", "AGE", "smoksort",      
                "pc1", "pc2", "pc3",           
                "early_late",           
                "surgery", "chemotx", "RADS")

formula_str <- paste("~", paste(final_vars, collapse = " + "))
X <- model.matrix(as.formula(formula_str), data = df_imputed)
X <- as.matrix(as.data.frame(X))[,-1]  # 去掉截距

var_names <- colnames(X)
p <- ncol(X)
n <- nrow(X)

cat("设计矩阵维度: n=", n, ", p=", p, "\n\n")

# ==============================================================================
# 5. Lasso 估计与交叉验证 (核心算法)
# ==============================================================================

cat("=== 第四-五部分：Lasso 估计与交叉验证 (正在计算，请稍候...) ===\n\n")

# Lasso 初始估计 parameters
n_lambda <- 100
nfold <- 5
tol <- 1.0e-6
maxiter <- 50000

# CV Lasso
cvobj_glmnet <- cv.glmnet(x=X, y=cbind(time=df_imputed$OS_month, status=df_imputed$DEAD), 
                          family="cox", alpha=1, standardize=F, nfolds = nfold, nlambda=n_lambda)

beta_glmnet <- as.vector(coef(glmnet(x=X, y=cbind(time=df_imputed$OS_month, status=df_imputed$DEAD), 
                                     family="cox", alpha=1, lambda=cvobj_glmnet$lambda.min, standardize=F, 
                                     thresh=tol, maxit=maxiter)))

cat("最优 lambda:", cvobj_glmnet$lambda.min, "\n")
cat("非零系数数量:", sum(beta_glmnet != 0), "\n")

# Hessian 矩阵计算
neg_loglik_glmnet <- 0
neg_dloglik_glmnet <- rep(0,p)
neg_ddloglik_glmnet <- matrix(0, nrow=p, ncol=p)
score_sq <- matrix(0, nrow=p, ncol=p)
neg_loglik_functions_cpp_ext(neg_loglik_glmnet, neg_dloglik_glmnet, neg_ddloglik_glmnet, score_sq,
                             X, df_imputed$OS_month, df_imputed$DEAD, beta_glmnet)
r <- eigen(score_sq)
r$values[r$values<=1.0e-14] <- 0

# 交叉验证选择 Multiplier
n_cv <- 5
n_multi <- 30
multiplier_seq <- exp(seq(from=log(0.001), to=log(5), length.out=n_multi))
alpha_cv <- 0.1
set.seed(12345)

all_cv_idx <- rep(1:n_cv, (n+1)/n_cv)
all_cv_idx <- sample(all_cv_idx, size=n, replace = T)
all_cvpl2 <- array(NA, length(multiplier_seq))

# 循环 CV
for(jj in 1:length(multiplier_seq)) {
  # 显示进度
  if(jj %% 5 == 0) cat("Processing multiplier", jj, "of", n_multi, "...\r")
  
  multiplier <- multiplier_seq[jj]
  cvpl2 <- 0
  
  for(k in 1:n_cv) {
    cv_idx <- which(all_cv_idx==k)
    train_x <- X[-c(cv_idx),]
    test_x <- X[cv_idx,]
    train_time <- df_imputed$OS_month[-c(cv_idx)]
    test_time <- df_imputed$OS_month[cv_idx]
    train_delta <- df_imputed$DEAD[-c(cv_idx)]
    test_delta <- df_imputed$DEAD[cv_idx]
    
    cvobj_train <- cv.glmnet(x=train_x, y=cbind(time=train_time, status=train_delta), family="cox",
                             alpha=1, standardize=F, nfolds = nfold, nlambda=n_lambda)
    beta_train <- as.vector(coef(glmnet(x=train_x, y=cbind(time=train_time, status=train_delta),
                                        family="cox", alpha=1, lambda=cvobj_train$lambda.min, standardize=F,
                                        thresh=tol, maxit=maxiter)))
    
    neg_log_train <- 0
    neg_dlog_train <- rep(0,p)
    neg_ddlog_train <- matrix(0, nrow=p, ncol=p)
    score_sq_train <- matrix(0, nrow=p, ncol=p)
    neg_loglik_functions_cpp_ext(neg_log_train, neg_dlog_train, neg_ddlog_train,
                                 score_sq_train, train_x, train_time, train_delta, beta_train)
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
      b_hat_new[j] <- beta_train[j] - as.numeric(m%*%neg_dlog_train)
      if (m[j] > 0) { se_new[j] <- sqrt(m[j]/nrow(train_x)) } else { se_new[j] <- NA }
    }
    pval_new <- 2*pnorm(abs(b_hat_new/se_new), lower.tail = F)
    tmp_beta <- b_hat_new*as.numeric(pval_new < (alpha_cv/p))
    cvpl2 <- cvpl2 + loglik_cpp_ext(X=test_x, time=test_time, delta=test_delta, beta=tmp_beta)
  }
  all_cvpl2[jj] <- cvpl2
}
cat("\n交叉验证完成!\n\n")

# ==============================================================================
# 6. 最终 Debiased 估计
# ==============================================================================

cat("=== 第六部分：最终 Debiased 估计 ===\n\n")

multiplier <- multiplier_seq[which.max(all_cvpl2)]
cat("Chosen multiplier:", multiplier, "\n")

b_hat_new <- array(NA, p)
se_new <- array(NA, p)
mu_new <- multiplier*sqrt(log(p)/n)
my_pos <- which(r$values > 0)
my_rank <- sum(r$values > 0)
Dmat <- diag(r$values[my_pos])
dvec <- rep(0,my_rank)
Amat <- t(rbind(-r$vectors[,my_pos]%*%Dmat, r$vectors[,my_pos]%*%Dmat))
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
  if (m[j] > 0) { se_new[j] <- sqrt(m[j]/n) } else { se_new[j] <- NA }
}

pval_new <- 2*pnorm(abs(b_hat_new/se_new), lower.tail=F)

# ==============================================================================
# 7. 整理结果与 FDR 校正
# ==============================================================================

cat("\n=== 第七部分：整理与保存 ===\n\n")

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

# 识别 SNP 变量
snp_cols_safe <- make.names(snp_cols)
is_snp <- var_names %in% snp_cols |
  var_names %in% snp_cols_safe |
  grepl("^rs\\d+", var_names, ignore.case = TRUE)

# 如果没匹配到，使用位置匹配
if (sum(is_snp) == 0) {
  cat("警告：使用位置匹配前", length(snp_cols), "列为 SNP\n")
  is_snp <- seq_len(p) <= length(snp_cols)
}

# --- 添加 FDR (Benjamini-Hochberg) ---
result$pvalue_FDR <- NA
if(sum(is_snp) > 0) {
  # 只对 SNP 进行校正
  result$pvalue_FDR[is_snp] <- p.adjust(result$pvalue[is_snp], method = "BH")
  cat("已计算 FDR (Benjamini-Hochberg)\n")
}

# 保存文件
write.csv(result, "MGH_LPWGS_cox_debiased_results_all.csv", row.names = FALSE)

# SNP 结果
snp_result <- result[is_snp, ]
snp_result <- snp_result[order(snp_result$pvalue), ]
write.csv(snp_result, "MGH_LPWGS_cox_debiased_results_SNPs.csv", row.names = FALSE)

# 协变量结果
covar_result <- result[!is_snp, ]
covar_result <- covar_result[order(covar_result$pvalue), ]
write.csv(covar_result, "MGH_LPWGS_cox_debiased_results_covariates.csv", row.names = FALSE)

cat("所有 CSV 结果已保存。\n")

# ==============================================================================
# 8. 可视化 (全部升级为 ggplot2 风格)
# ==============================================================================

cat("\n=== 第八部分：可视化 (Full ggplot2 Style) ===\n\n")

# ---------------------------------------------------------
# A. 火山图 (Volcano Plot)
# ---------------------------------------------------------
if (nrow(snp_result) > 0) {
  volcano_data <- snp_result
  # 修复 1e-300 问题
  volcano_data$pvalue[volcano_data$pvalue == 0] <- 1e-300
  
  volcano_data$log10p <- -log10(volcano_data$pvalue)
  bonferroni_threshold <- 0.05 / nrow(volcano_data)
  
  p_volcano <- ggplot(volcano_data, aes(x=beta, y=log10p)) +
    geom_point(alpha=0.6, color="steelblue", size=2) +
    geom_hline(yintercept = -log10(bonferroni_threshold), linetype="dashed", color="red") +
    geom_hline(yintercept = -log10(0.05), linetype="dotted", color="orange") +
    annotate("text", x = max(volcano_data$beta, na.rm=T) * 0.7, y = -log10(bonferroni_threshold) + 0.3,
             label = "Bonferroni", color = "red", size = 3) +
    annotate("text", x = max(volcano_data$beta, na.rm=T) * 0.7, y = -log10(0.05) + 0.3,
             label = "P=0.05", color = "orange", size = 3) +
    xlab("Effect size (beta)") + ylab("-log10(p-value)") +
    ggtitle(paste0("Volcano Plot: ", nrow(volcano_data), " SNPs - MGH LPWGS")) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  ggsave("MGH_LPWGS_volcano_plot.png", p_volcano, width = 10, height = 8, dpi = 150)
  cat("已保存: MGH_LPWGS_volcano_plot.png\n")
}

# ---------------------------------------------------------
# B. 森林图 - Top SNPs (Forest Plot)
# ---------------------------------------------------------
if (nrow(snp_result) > 0) {
  top_snps <- head(snp_result, 20)
  top_snps$Variable <- factor(top_snps$Variable, levels = rev(top_snps$Variable))
  
  p_forest <- ggplot(top_snps, aes(x = HR, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
    geom_point(aes(color = pvalue < 0.05), size = 3) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                       labels = c("TRUE" = "P < 0.05", "FALSE" = "P ≥ 0.05"),
                       name = "Significance") +
    scale_x_log10() +
    labs(title = "Forest Plot: Top 20 SNPs - MGH LPWGS", x = "Hazard Ratio (log scale)", y = "") +
    theme_bw() + theme(legend.position = "bottom")
  
  ggsave("MGH_LPWGS_forest_plot_SNPs.png", p_forest, width = 10, height = 8, dpi = 150)
  cat("已保存: MGH_LPWGS_forest_plot_SNPs.png\n")
}

# ---------------------------------------------------------
# C. 森林图 - 协变量 (Forest Plot Covariates)
# ---------------------------------------------------------
if (nrow(covar_result) > 0) {
  covar_result$Variable <- factor(covar_result$Variable, levels = rev(covar_result$Variable))
  
  p_covar <- ggplot(covar_result, aes(x = HR, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
    geom_point(aes(color = pvalue < 0.05), size = 3) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                       labels = c("TRUE" = "P < 0.05", "FALSE" = "P ≥ 0.05"),
                       name = "Significance") +
    scale_x_log10() +
    labs(title = "Forest Plot: Clinical Covariates - MGH LPWGS", x = "Hazard Ratio (log scale)", y = "") +
    theme_bw() + theme(legend.position = "bottom")
  
  ggsave("MGH_LPWGS_forest_plot_covariates.png", p_covar, width = 10, height = 6, dpi = 150)
  cat("已保存: MGH_LPWGS_forest_plot_covariates.png\n")
}

# ---------------------------------------------------------
# D. 曼哈顿图 (Manhattan Plot) - ggplot2 + ggrepel
# ---------------------------------------------------------
snp_pos_file <- pick_file("Select the candidate SNP position file", "SNP_POS")
if (file.exists(snp_pos_file)) {
  
  snp_pos <- read.csv(snp_pos_file, stringsAsFactors = FALSE)
  plot_data <- snp_result
  # 提取纯 rsID 以匹配位置文件
  plot_data$rsID_clean <- gsub("_[A-Za-z0-9]+$", "", plot_data$Variable)
  
  manhattan_data <- merge(plot_data, snp_pos, by.x = "rsID_clean", by.y = "SNP", all.x = FALSE)
  manhattan_data$log10P <- -log10(manhattan_data$pvalue)
  manhattan_data$CHR <- as.numeric(manhattan_data$Chr)
  manhattan_data <- manhattan_data[!is.na(manhattan_data$CHR), ]
  
  if (nrow(manhattan_data) > 0) {
    n_snps <- nrow(manhattan_data)
    bonf_p <- 0.05 / n_snps
    
    p_manhattan <- ggplot(manhattan_data, aes(x = CHR, y = log10P)) +
      geom_point(aes(color = factor(CHR %% 2)), size = 3, alpha = 0.8) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange", linewidth = 0.8) +
      geom_hline(yintercept = -log10(bonf_p), linetype = "dashed", color = "red", linewidth = 0.8) +
      # 使用 ggrepel 自动避让标签
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
           title = "Manhattan Plot - MGH LPWGS") +
      annotate("text", x = max(manhattan_data$CHR) - 1, y = -log10(0.05) + 0.2, 
               label = "P = 0.05", color = "orange", size = 3.5, hjust = 0) +
      annotate("text", x = max(manhattan_data$CHR) - 1, y = -log10(bonf_p) + 0.2, 
               label = "Bonferroni", color = "red", size = 3.5, hjust = 0) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 12))
    
    ggsave("MGH_LPWGS_manhattan_plot.png", p_manhattan, width = 12, height = 6, dpi = 200)
    cat("已保存: MGH_LPWGS_manhattan_plot.png\n")
  }
} else {
  cat("警告：找不到位置文件 77_SNPs_with_Positions.csv，跳过曼哈顿图。\n")
}

# ==============================================================================
# 9. 最终结果汇总 (Results Summary)
# ==============================================================================

cat("\n========================================\n")
cat("=== 最终结果汇总 (含 FDR P值) ===\n")
cat("========================================\n\n")

# 1. 打印 P < 0.05 的 SNP (包含 pvalue_FDR)
sig_snps_nominal <- snp_result[snp_result$pvalue < 0.05, ]
cat("P < 0.05 的 SNP 数量:", nrow(sig_snps_nominal), "\n")
if (nrow(sig_snps_nominal) > 0) {
  cat("\n--- P < 0.05 的 SNP (Top Results) ---\n")
  # 打印所需列
  cols_to_print <- c("Variable", "beta", "HR", "HR_lower", "HR_upper", "pvalue", "pvalue_FDR")
  print(sig_snps_nominal[, cols_to_print])
}

# 2. 打印并保存 Bonferroni 显著 SNP
n_snp_total <- nrow(snp_result)
bonf_threshold <- 0.05 / n_snp_total
sig_snps_bonf <- snp_result[snp_result$pvalue < bonf_threshold, ]

cat("\nBonferroni 校正后显著 SNP 数量:", nrow(sig_snps_bonf), "\n")
if (nrow(sig_snps_bonf) > 0) {
  cat("\n--- Bonferroni 显著的 SNP ---\n")
  print(sig_snps_bonf[, c("Variable", "beta", "HR", "pvalue")])
  write.csv(sig_snps_bonf, "MGH_LPWGS_significant_SNPs_bonferroni.csv", row.names = FALSE)
  cat("已保存: MGH_LPWGS_significant_SNPs_bonferroni.csv\n")
}

# 3. 打印显著协变量
sig_covar <- covar_result[covar_result$pvalue < 0.05, ]
cat("\n显著协变量 (P < 0.05):\n")
if (nrow(sig_covar) > 0) {
  print(sig_covar[, c("Variable", "beta", "HR", "pvalue")])
} else {
  cat("无显著协变量\n")
}

cat("\n=== 分析圆满完成！(所有图表和数据均已生成) ===\n")