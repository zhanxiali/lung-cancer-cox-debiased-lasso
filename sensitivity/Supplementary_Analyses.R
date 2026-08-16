########################################################################
#  Supplementary analyses for the NSCLC lung-function SNP survival study
#
#  This is the version used to generate the published supplementary output.
#
#  Task 1  Effect-allele frequency concordance across subgroups
#            -> SuppTable_Sx_EAF_Concordance.csv        (Supplementary Table S3)
#  Task 2  Scaled Schoenfeld residuals and Grambsch-Therneau tests
#            -> SuppFig_S2_Schoenfeld_<cohort>.pdf      (Supplementary Figure S1)
#            -> SuppTable_PH_Test_Results.csv           (Supplementary Table S6)
#  Task 3  Treatment-covariate sensitivity analysis (framework)
#            -> SuppTable_Sy_Treatment_Sensitivity.csv
#          The treatment sensitivity results reported in the paper
#          (Supplementary Table S1) were produced by
#          meta/Meta_sensitivity_NO_TREATMENT.R.
#
#  Genotype column naming differs by platform. HSPH and OncoArray .raw files
#  use position-encoded names (X1.22653424.C.G_G); MGH uses rsID names
#  (rs2794359_A). build_rsid_map() handles both: for the first two it matches
#  against the Variable_original column of the per-cohort results CSV, and for
#  MGH it strips the allele suffix directly.
#
#  Inputs are selected at run time; no paths are hard-coded. Set the
#  environment variables below to run without prompts:
#    OUTPUT_DIR
#    RAW_HSPH   RAW_ONCOARRAY   RAW_MGH        aligned PLINK .raw files
#    CLIN_HSPH  CLIN_ONCOARRAY  CLIN_MGH       cleaned clinical CSVs
#    RES_HSPH   RES_ONCOARRAY   RES_MGH        per-cohort debiased results CSVs
#    SIM_UNIVLIB  APPLY_DEBIAS  META_ALL       Task 3 only
#
#  Requires: dplyr, survival (Task 3 additionally requires quadprog, Rcpp,
#  RcppArmadillo, glmnet, MASS, lpSolve, meta, metafor, mvtnorm).
########################################################################

library(dplyr)
library(survival)

# ---------------------------------------------------------------------------
# Input selection. Paths are chosen at run time rather than hard-coded. In an
# interactive session a file dialog is shown; otherwise the path is read from
# the console. Set the environment variable named in the second argument to run
# unattended (for example inside a compute capsule).
# ---------------------------------------------------------------------------
pick_file <- function(msg, envvar = NA_character_) {
  if (!is.na(envvar)) {
    v <- Sys.getenv(envvar, unset = "")
    if (nzchar(v) && file.exists(v)) { cat("Using ", envvar, " = ", v, "\n", sep = ""); return(v) }
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

pick_dir <- function(msg, envvar = "OUTPUT_DIR") {
  v <- Sys.getenv(envvar, unset = "")
  if (nzchar(v)) { dir.create(v, showWarnings = FALSE, recursive = TRUE)
                   cat("Using ", envvar, " = ", v, "\n", sep = ""); return(v) }
  cat("\n>>> ", msg, "\n", sep = "")
  d <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() && rstudioapi::hasFun("selectDirectory")) {
    d <- tryCatch(rstudioapi::selectDirectory(caption = msg), error = function(e) NULL)
  }
  if ((is.null(d) || !length(d)) && interactive() && exists("choose.dir")) {
    d <- tryCatch(utils::choose.dir(caption = msg), error = function(e) NULL)
  }
  while (is.null(d) || !length(d) || is.na(d) || !nzchar(d)) {
    d <- trimws(readline("    Output directory: "))
    d <- gsub('^"|"$', "", gsub("\\\\", "/", d))
    if (nzchar(d)) dir.create(d, showWarnings = FALSE, recursive = TRUE) else d <- NULL
  }
  cat("    Output directory: ", d, "\n", sep = "")
  d
}

# ======================================================================
#  公共工具函数
# ======================================================================

output_dir <- pick_dir("Select the directory for output files", "OUTPUT_DIR")

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

# ------------------------------------------------------------------
# 核心修复: .raw列名 -> rsID 映射
# 使用per-cohort结果CSV中已有的映射信息
# ------------------------------------------------------------------
build_rsid_map <- function(results_csv_path, raw_snp_cols) {
  res <- read.csv(results_csv_path, stringsAsFactors = FALSE)
  res_snp <- res[grepl("^rs", res$Variable), ]
  res_snp$rsID_clean <- gsub("_[ATCG]+$", "", res_snp$Variable)
  
  if ("Variable_original" %in% names(res_snp)) {
    # HSPH/Oncoarray: Variable_original = .raw列名
    mapping <- data.frame(
      rsID = res_snp$rsID_clean,
      Variable_original = res_snp$Variable_original,
      stringsAsFactors = FALSE
    )
    map_result <- data.frame(raw_col = raw_snp_cols, rsID = NA_character_, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(mapping))) {
      vo <- mapping$Variable_original[i]
      idx <- which(raw_snp_cols == vo)
      if (length(idx) == 1) {
        map_result$rsID[idx] <- mapping$rsID[i]
      } else {
        vo_base <- gsub("_[ATCG]+$", "", vo)
        raw_bases <- gsub("_[ATCG]+$", "", raw_snp_cols)
        idx2 <- which(raw_bases == vo_base)
        if (length(idx2) >= 1) map_result$rsID[idx2[1]] <- mapping$rsID[i]
      }
    }
  } else {
    # MGH: Variable本身就是 rsID_allele 格式
    map_result <- data.frame(
      raw_col = raw_snp_cols,
      rsID = gsub("_[ATCG]+$", "", raw_snp_cols),
      stringsAsFactors = FALSE
    )
  }
  
  n_mapped <- sum(!is.na(map_result$rsID))
  cat(paste0("  rsID映射: ", n_mapped, "/", nrow(map_result), " 成功\n"))
  if (n_mapped < nrow(map_result)) {
    unmapped <- map_result$raw_col[is.na(map_result$rsID)]
    cat("  未映射示例:", paste(head(unmapped, 3), collapse=", "), "\n")
  }
  return(map_result)
}


########################################################################
#  TASK 1: EAF Concordance Table
########################################################################
run_task1_eaf <- function() {
  cat("\n============================================================\n")
  cat("  TASK 1: EAF Concordance Table\n")
  cat("============================================================\n\n")
  
  cat("--- 读取三个subgroup的数据 ---\n\n")
  
  # HSPH
  raw_hsph <- read.table(pick_file("Select the aligned PLINK .raw genotype file (HSPH)", "RAW_HSPH"), header=TRUE, stringsAsFactors=FALSE)
  snp_cols_hsph <- colnames(raw_hsph)[7:ncol(raw_hsph)]
  map_hsph <- build_rsid_map(pick_file("Select the debiased results CSV (HSPH)", "RES_HSPH"), snp_cols_hsph)
  
  # Oncoarray
  raw_onco <- read.table(pick_file("Select the aligned PLINK .raw genotype file (OncoArray)", "RAW_ONCO"), header=TRUE, stringsAsFactors=FALSE)
  snp_cols_onco <- colnames(raw_onco)[7:ncol(raw_onco)]
  map_onco <- build_rsid_map(pick_file("Select the debiased results CSV (OncoArray)", "RES_ONCO"), snp_cols_onco)
  
  # MGH
  raw_mgh <- read.table(pick_file("Select the aligned PLINK .raw genotype file (MGH)", "RAW_MGH"), header=TRUE, stringsAsFactors=FALSE)
  snp_cols_mgh <- colnames(raw_mgh)[7:ncol(raw_mgh)]
  map_mgh <- build_rsid_map(pick_file("Select the debiased results CSV (MGH)", "RES_MGH"), snp_cols_mgh)
  
  # 计算EAF
  cat("\n--- 计算EAF ---\n")
  calc_eaf <- function(raw_df, snp_cols, rsid_map) {
    snp_data <- as.data.frame(lapply(raw_df[, snp_cols, drop=FALSE],
                                      function(x) as.numeric(as.character(x))))
    eaf <- colMeans(snp_data, na.rm=TRUE) / 2
    result <- data.frame(raw_col=snp_cols, EAF=eaf, N=colSums(!is.na(snp_data)),
                          stringsAsFactors=FALSE, row.names=NULL)
    result <- merge(result, rsid_map, by="raw_col", all.x=TRUE)
    result <- result[!is.na(result$rsID), ]
    return(result)
  }
  
  eaf_h <- calc_eaf(raw_hsph, snp_cols_hsph, map_hsph)
  eaf_o <- calc_eaf(raw_onco, snp_cols_onco, map_onco)
  eaf_m <- calc_eaf(raw_mgh, snp_cols_mgh, map_mgh)
  cat("  HSPH:", nrow(eaf_h), " Oncoarray:", nrow(eaf_o), " MGH:", nrow(eaf_m), "\n")
  
  # 合并
  cat("\n--- 合并 ---\n")
  eaf_merge <- eaf_h %>% select(rsID, EAF_HSPH=EAF) %>%
    full_join(eaf_o %>% select(rsID, EAF_Onco=EAF), by="rsID") %>%
    full_join(eaf_m %>% select(rsID, EAF_MGH=EAF), by="rsID")
  
  n3 <- sum(!is.na(eaf_merge$EAF_HSPH) & !is.na(eaf_merge$EAF_Onco) & !is.na(eaf_merge$EAF_MGH))
  n2 <- sum(rowSums(!is.na(eaf_merge[,c("EAF_HSPH","EAF_Onco","EAF_MGH")])) >= 2)
  cat("  总:", nrow(eaf_merge), " 三subgroup共有:", n3, " >=2共有:", n2, "\n")
  
  # Pairwise差异
  eaf_merge <- eaf_merge %>% mutate(
    diff_HSPH_Onco = abs(EAF_HSPH - EAF_Onco),
    diff_HSPH_MGH  = abs(EAF_HSPH - EAF_MGH),
    diff_Onco_MGH  = abs(EAF_Onco - EAF_MGH)
  )
  
  cat("\n--- Pairwise |EAF差| ---\n")
  for (col in c("diff_HSPH_Onco","diff_HSPH_MGH","diff_Onco_MGH")) {
    v <- eaf_merge[[col]][!is.na(eaf_merge[[col]])]
    if (length(v) > 0) cat(sprintf("  %s: n=%d, median=%.4f, max=%.4f\n", col, length(v), median(v), max(v)))
    else cat(sprintf("  %s: 无配对数据\n", col))
  }
  
  # Sign-flip
  eaf_merge <- eaf_merge %>% mutate(
    flip_Onco = ifelse(!is.na(EAF_HSPH) & !is.na(EAF_Onco),
                        abs(EAF_HSPH-(1-EAF_Onco)) < abs(EAF_HSPH-EAF_Onco), NA),
    flip_MGH  = ifelse(!is.na(EAF_HSPH) & !is.na(EAF_MGH),
                        abs(EAF_HSPH-(1-EAF_MGH)) < abs(EAF_HSPH-EAF_MGH), NA)
  )
  cat("\n  需flip: Onco=", sum(eaf_merge$flip_Onco==TRUE,na.rm=T),
      " MGH=", sum(eaf_merge$flip_MGH==TRUE,na.rm=T), "\n")
  
  # 异常
  outliers <- eaf_merge %>% filter(diff_HSPH_Onco>0.10 | diff_HSPH_MGH>0.10 | diff_Onco_MGH>0.10)
  cat("\n  |EAF差|>0.10:", nrow(outliers), "个\n")
  if (nrow(outliers)>0) print(outliers %>% select(rsID, EAF_HSPH, EAF_Onco, EAF_MGH, diff_HSPH_Onco, diff_HSPH_MGH, diff_Onco_MGH))
  
  # 保存
  f <- file.path(output_dir, "SuppTable_Sx_EAF_Concordance.csv")
  write.csv(eaf_merge %>% arrange(rsID), f, row.names=FALSE)
  cat("\n保存:", f, "\nTask 1 完成!\n")
  return(invisible(eaf_merge))
}


########################################################################
#  TASK 2: Schoenfeld Residual Plots
########################################################################
run_task2_schoenfeld <- function() {
  cat("\n============================================================\n")
  cat("  TASK 2: Schoenfeld Residuals (cox.zph)\n")
  cat("============================================================\n\n")
  
  top11 <- c("rs11022690","rs10987386","rs6006399","rs17821105",
             "rs11227223","rs35684381","rs17032590","rs196025",
             "rs72743477","rs35797611","rs72811372")
  all_ph <- list()
  
  for (cohort in c("HSPH","Oncoarray","MGH")) {
    cat("\n--- ", cohort, " ---\n")
    
    raw_df <- read.table(pick_file(paste0("Select the aligned PLINK .raw genotype file (", cohort, ")"),
                                   paste0("RAW_", toupper(cohort))), header=TRUE, stringsAsFactors=FALSE)
    clin_df <- read.csv(pick_file(paste0("Select the cleaned clinical CSV (", cohort, ")"),
                                   paste0("CLIN_", toupper(cohort))), stringsAsFactors=FALSE)
    snp_cols <- colnames(raw_df)[7:ncol(raw_df)]
    rsid_map <- build_rsid_map(pick_file(paste0("Select the debiased results CSV (", cohort, ")"),
                                          paste0("RES_", toupper(cohort))), snp_cols)
    
    df_snp <- raw_df[, c("IID", snp_cols)]
    df_snp$IID <- as.character(df_snp$IID)
    clin_df$IID <- as.character(clin_df$IID)
    df_all <- merge(df_snp, clin_df, by="IID")
    for (v in c("SEX","smoksort","early_late","RADS","chemotx","surgery"))
      df_all[[v]] <- as.factor(df_all[[v]])
    df_all <- df_all %>% filter(!is.na(OS_month), !is.na(DEAD))
    df_imp <- impute_mode_mean(df_all)
    df_imp[snp_cols] <- lapply(df_imp[snp_cols], function(x) as.numeric(as.character(x)))
    
    map_top <- rsid_map %>% filter(rsID %in% top11)
    if (nrow(map_top)==0) { cat("  无top SNP，跳过\n"); next }
    cat("  可用top SNPs:", nrow(map_top), "-", paste(map_top$rsID, collapse=", "), "\n")
    
    fvars <- c(map_top$raw_col, "AGE","SEX","smoksort","pc1","pc2","pc3",
               "early_late","surgery","chemotx","RADS")
    fml <- as.formula(paste("Surv(OS_month, DEAD)~", paste(fvars, collapse="+")))
    
    fit <- tryCatch(coxph(fml, data=df_imp), error=function(e) { cat("  错误:",e$message,"\n"); NULL })
    if (is.null(fit)) next
    
    ph <- cox.zph(fit, transform="log")
    ph_tbl <- as.data.frame(ph$table)
    rn <- rownames(ph_tbl)
    for (j in 1:nrow(map_top)) rn <- gsub(map_top$raw_col[j], map_top$rsID[j], rn, fixed=TRUE)
    rownames(ph_tbl) <- rn
    
    cat("\n  PH检验:\n"); print(ph_tbl)
    
    ph_df <- data.frame(Cohort=cohort, Variable=rn, chisq=ph_tbl$chisq,
                          df=ph_tbl$df, p=ph_tbl$p, sig001=ph_tbl$p<0.01,
                          stringsAsFactors=FALSE, row.names=NULL)
    all_ph[[cohort]] <- ph_df
    
    # 图
    pdf_f <- file.path(output_dir, paste0("SuppFig_S2_Schoenfeld_", cohort, ".pdf"))
    pdf(pdf_f, width=10, height=8)
    for (i in 1:nrow(map_top)) {
      vi <- which(names(fit$coefficients)==map_top$raw_col[i])
      if (length(vi)>0) {
        plot(ph, var=vi, main=paste0(cohort,": ",map_top$rsID[i]," (P=",format(ph$table[vi,"p"],digits=3,scientific=T),")"))
        abline(h=0, lty=2, col="red")
      }
    }
    for (vc in c("AGE","early_late1","surgery1")) {
      vi <- which(names(fit$coefficients)==vc)
      if (length(vi)>0) {
        plot(ph, var=vi, main=paste0(cohort,": ",vc," (P=",format(ph$table[vi,"p"],digits=3,scientific=T),")"))
        abline(h=0, lty=2, col="red")
      }
    }
    dev.off()
    cat("  图:", pdf_f, "\n")
  }
  
  if (length(all_ph)>0) {
    pc <- do.call(rbind, all_ph)
    f <- file.path(output_dir, "SuppTable_PH_Test_Results.csv")
    write.csv(pc, f, row.names=FALSE)
    cat("\n保存:", f, "\n\n=== 汇总 ===\n")
    for (co in names(all_ph)) {
      ph <- all_ph[[co]]
      nv <- sum(ph$sig001[ph$Variable!="GLOBAL"], na.rm=T)
      gp <- ph$p[ph$Variable=="GLOBAL"]
      cat(sprintf("  %s: %d变量P<0.01, Global P=%s\n", co, nv, format(gp, digits=3)))
      for (rs in top11) {
        row <- ph[grepl(rs, ph$Variable),]
        if (nrow(row)>0) cat(sprintf("    %s: P=%s %s\n", rs, format(row$p[1],digits=3,scientific=T),
                                       ifelse(row$sig001[1],"VIOLATION","OK")))
      }
    }
  }
  cat("\nTask 2 完成!\n")
}


########################################################################
#  TASK 3: Treatment Sensitivity (框架 — 需要完整pipeline)
########################################################################
run_task3_treatment_sensitivity <- function() {
  cat("\n============================================================\n")
  cat("  TASK 3: Treatment Sensitivity Analysis\n")
  cat("============================================================\n\n")
  
  library(quadprog); library(mvtnorm); library(Rcpp); library(RcppArmadillo)
  library(glmnet); library(MASS); library(lpSolve); library(meta); library(metafor)
  
  sourceCpp(pick_file("Select sim_univLib.cpp", "SIM_UNIVLIB"))
  source(pick_file("Select applydebiaslasso.R", "APPLY_DEBIAS"))
  
  sens_results <- list()
  
  for (cohort in c("HSPH","Oncoarray","MGH")) {
    cat("\n--- ", cohort, " ---\n")
    cat(paste0("  运行", cohort, "? (y/n): "))
    if (tolower(readline()) != "y") { cat("  跳过\n"); next }
    
    raw_df <- read.table(pick_file(paste0("Select the aligned PLINK .raw genotype file (", cohort, ")"),
                                   paste0("RAW_", toupper(cohort))), header=TRUE, stringsAsFactors=FALSE)
    clin_df <- read.csv(pick_file(paste0("Select the cleaned clinical CSV (", cohort, ")"),
                                   paste0("CLIN_", toupper(cohort))), stringsAsFactors=FALSE)
    snp_cols <- colnames(raw_df)[7:ncol(raw_df)]
    rsid_map <- build_rsid_map(pick_file(paste0("Select the debiased results CSV (", cohort, ")"),
                                          paste0("RES_", toupper(cohort))), snp_cols)
    
    df_snp <- raw_df[, c("IID", snp_cols)]
    df_snp$IID <- as.character(df_snp$IID)
    clin_df$IID <- as.character(clin_df$IID)
    df_all <- merge(df_snp, clin_df, by="IID")
    df_all$SEX <- as.factor(df_all$SEX)
    df_all$smoksort <- as.factor(df_all$smoksort)
    df_all$early_late <- as.factor(df_all$early_late)
    df_all <- df_all %>% filter(!is.na(OS_month), !is.na(DEAD))
    df_imp <- impute_mode_mean(df_all)
    df_imp[snp_cols] <- lapply(df_imp[snp_cols], function(x) as.numeric(as.character(x)))
    
    # 去除治疗变量
    final_vars <- c(snp_cols, "SEX","AGE","smoksort","pc1","pc2","pc3","early_late")
    X <- model.matrix(as.formula(paste("~",paste(final_vars,collapse="+"))), data=df_imp)
    X <- as.matrix(as.data.frame(X))[,-1]
    var_names <- colnames(X); p <- ncol(X); n <- nrow(X)
    cat("  n=",n,", p=",p,"(不含treatment)\n")
    
    # Lasso
    cat("  Lasso...\n")
    cvobj <- cv.glmnet(x=X, y=cbind(time=df_imp$OS_month,status=df_imp$DEAD),
                        family="cox",alpha=1,standardize=F,nfolds=5,nlambda=100)
    b_lasso <- as.vector(coef(glmnet(x=X, y=cbind(time=df_imp$OS_month,status=df_imp$DEAD),
                                       family="cox",alpha=1,lambda=cvobj$lambda.min,standardize=F,thresh=1e-6,maxit=50000)))
    cat("  非零:", sum(b_lasso!=0), "\n")
    
    # Score矩阵
    cat("  Score矩阵...\n")
    nll <- 0; ndll <- rep(0,p); nddll <- matrix(0,p,p); ssq <- matrix(0,p,p)
    neg_loglik_functions_cpp_ext(nll, ndll, nddll, ssq, X, df_imp$OS_month, df_imp$DEAD, b_lasso)
    r <- eigen(ssq); r$values[r$values<=1e-14] <- 0
    
    # CV选mu
    cat("  CV选mu...\n")
    n_cv <- 5; n_multi <- 30
    mseq <- exp(seq(log(0.001),log(5),length.out=n_multi))
    set.seed(12345); cv_grp <- sample(rep(1:n_cv,length.out=n))
    cv_nll <- matrix(NA, n_cv, n_multi)
    
    for (ic in 1:n_cv) {
      ite <- which(cv_grp==ic); itr <- which(cv_grp!=ic)
      Xtr <- X[itr,,drop=F]; Xte <- X[ite,,drop=F]
      ttr <- df_imp$OS_month[itr]; tte <- df_imp$OS_month[ite]
      dtr <- df_imp$DEAD[itr]; dte <- df_imp$DEAD[ite]
      ntr <- length(itr)
      
      cv_l <- cv.glmnet(x=Xtr,y=cbind(time=ttr,status=dtr),family="cox",alpha=1,standardize=F,nfolds=5,nlambda=100)
      btr <- as.vector(coef(glmnet(x=Xtr,y=cbind(time=ttr,status=dtr),family="cox",alpha=1,
                                     lambda=cv_l$lambda.min,standardize=F,thresh=1e-6,maxit=50000)))
      
      nll2 <- 0; ndll2 <- rep(0,p); nddll2 <- matrix(0,p,p); ssq2 <- matrix(0,p,p)
      neg_loglik_functions_cpp_ext(nll2,ndll2,nddll2,ssq2,Xtr,ttr,dtr,btr)
      rtr <- eigen(ssq2); rtr$values[rtr$values<=1e-14] <- 0
      
      for (im in 1:n_multi) {
        gn <- mseq[im]*sqrt(log(p)/ntr)
        bd <- rep(NA,p)
        for (j in 1:p) {
          ej <- rep(0,p); ej[j] <- 1
          tryCatch({ mj <- sim_univ_debias_cpp(ej,rtr$vectors,rtr$values,gn)
                      bd[j] <- btr[j]-sum(mj*ndll2)/ntr }, error=function(e) bd[j] <<- NA)
        }
        be <- bd; be[is.na(be)] <- btr[is.na(be)]
        nll3 <- 0; ndll3 <- rep(0,p); nddll3 <- matrix(0,p,p); ssq3 <- matrix(0,p,p)
        neg_loglik_functions_cpp_ext(nll3,ndll3,nddll3,ssq3,Xte,tte,dte,be)
        cv_nll[ic,im] <- nll3/length(ite)
      }
      cat("    CV fold",ic,"完成\n")
    }
    
    bi <- which.min(colMeans(cv_nll,na.rm=T))
    bg <- mseq[bi]*sqrt(log(p)/n)
    cat("  最优gamma:",round(bg,6),"\n")
    
    # Debiasing
    cat("  Debiasing...\n")
    b_new <- se_new <- rep(NA,p)
    for (j in 1:p) {
      ej <- rep(0,p); ej[j] <- 1
      tryCatch({ mj <- sim_univ_debias_cpp(ej,r$vectors,r$values,bg)
                  b_new[j] <- b_lasso[j]-sum(mj*ndll)/n
                  se_new[j] <- sqrt(mj[j]/n) }, error=function(e) { b_new[j] <<- NA; se_new[j] <<- NA })
    }
    
    pval <- 2*pnorm(-abs(b_new/se_new))
    res <- data.frame(Variable=var_names, beta=round(b_new,6), SE=round(se_new,6), pvalue=pval,
                       HR=round(exp(b_new),4), HR_lower=round(exp(b_new-1.96*se_new),4),
                       HR_upper=round(exp(b_new+1.96*se_new),4), stringsAsFactors=FALSE)
    res$raw_col <- var_names
    res <- merge(res, rsid_map[,c("raw_col","rsID")], by="raw_col", all.x=TRUE)
    res$rsID[is.na(res$rsID)] <- res$Variable[is.na(res$rsID)]
    snp_mask <- grepl("^rs",res$rsID)
    res$pvalue_FDR <- NA
    if (sum(snp_mask)>0) res$pvalue_FDR[snp_mask] <- p.adjust(res$pvalue[snp_mask],"BH")
    
    f <- file.path(output_dir, paste0(cohort,"_results_NO_TREATMENT.csv"))
    write.csv(res, f, row.names=FALSE)
    cat("  保存:",f,"\n")
    sens_results[[cohort]] <- res
    
    cat("  Top SNPs:\n")
    print(res %>% filter(grepl("^rs",rsID)) %>% arrange(pvalue) %>% head(5) %>% select(rsID,beta,pvalue,HR))
  }
  
  # Meta
  if (length(sens_results)>=2) {
    cat("\n--- Meta分析 ---\n")
    ssnp <- lapply(sens_results, function(d) d %>% filter(grepl("^rs",rsID)) %>% select(rsID,beta,SE))
    allrs <- unique(unlist(lapply(ssnp, function(x) x$rsID)))
    ms <- data.frame()
    for (rs in allrs) {
      bs <- c(); ss <- c()
      for (co in names(ssnp)) {
        row <- ssnp[[co]] %>% filter(rsID==rs)
        if (nrow(row)==1 && !is.na(row$beta) && !is.na(row$SE) && row$SE>0) { bs <- c(bs,row$beta); ss <- c(ss,row$SE) }
      }
      if (length(bs)>=2) {
        w <- 1/ss^2; bf <- sum(w*bs)/sum(w); sf <- sqrt(1/sum(w)); pf <- 2*pnorm(-abs(bf/sf))
        Q <- sum(w*(bs-bf)^2); I2 <- max(0,(Q-(length(bs)-1))/Q*100)
        ms <- rbind(ms, data.frame(rsID=rs,n=length(bs),beta_sens=round(bf,6),SE_sens=round(sf,6),
                                    p_sens=pf,HR_sens=round(exp(bf),4),I2_sens=round(I2,1),stringsAsFactors=F))
      }
    }
    ms$FDR_sens <- p.adjust(ms$p_sens,"BH")
    
    main <- read.csv(pick_file("Select Meta_analysis_results_all.csv", "META_ALL"), stringsAsFactors=FALSE)
    comp <- merge(main %>% select(rsID,beta_fixed,SE_fixed,pvalue_fixed,HR_fixed), ms, by="rsID") %>%
      mutate(dir_same=sign(beta_fixed)==sign(beta_sens), HR_diff=abs(HR_sens-HR_fixed)) %>% arrange(pvalue_fixed)
    
    f <- file.path(output_dir, "SuppTable_Sy_Treatment_Sensitivity.csv")
    write.csv(comp, f, row.names=FALSE)
    cat("比较表:",f,"\n")
    
    cat("\n  === rs11022690 ===\n")
    r11 <- comp %>% filter(rsID=="rs11022690")
    if (nrow(r11)>0) cat(sprintf("  主分析: HR=%.3f P=%s\n  敏感性: HR=%.3f P=%s\n  方向一致: %s\n",
                                   r11$HR_fixed, format(r11$pvalue_fixed,scientific=T),
                                   r11$HR_sens, format(r11$p_sens,scientific=T),
                                   ifelse(r11$dir_same,"YES","NO")))
    cat(sprintf("  整体方向一致: %d/%d\n", sum(comp$dir_same,na.rm=T), nrow(comp)))
  }
  cat("\nTask 3 完成!\n")
}


########################################################################
cat("\n============================================================\n")
cat("  NSCLC 补充分析 (修复版 v2)\n")
cat("  1=EAF  2=Schoenfeld  3=Treatment  0=全部\n")
cat("============================================================\n")
task <- readline("任务编号: ")
if (task %in% c("0","1")) run_task1_eaf()
if (task %in% c("0","2")) run_task2_schoenfeld()
if (task %in% c("0","3")) run_task3_treatment_sensitivity()
cat("\n完成!\n")
