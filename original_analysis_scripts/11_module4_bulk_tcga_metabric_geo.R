# 11_module4_bulk_tcga_metabric_geo.R
# 模块 4：在本地 TCGA-BRCA (PanCancer Atlas + Xena)、METABRIC 及外部 GEO 队列 (GSE20685 / GSE42568)
# 对全部 8 个纯 cis-MR 阳性基因进行表达差异、亚型分层、生存预后验证与严格方向一致性审计

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(GEOquery)
  library(Biobase)
})

options(timeout = 150)

bulk_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
geo_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/geo_cache")
res_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)

sig_cis <- fread(file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))
target_genes <- unique(as.character(sig_cis$gene_symbol))
cat("Module 4 target genes (8 cis-MR positive genes):", paste(target_genes, collapse = ", "), "\n")

parse_status <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.integer(sub(":.*", "", trimws(x))))
  out
}

# 辅助函数：对单个基因列表跑表达差异 & Cox/KM 生存
run_survival_stats <- function(expr_vec, time_vec, event_vec, cohort_name, endpoint, stratum, gene) {
  ok <- !is.na(expr_vec) & !is.na(time_vec) & !is.na(event_vec) & time_vec > 0
  if (sum(ok) < 30 || sum(event_vec[ok] == 1) < 5 || length(unique(expr_vec[ok])) < 5) return(NULL)
  e <- as.numeric(expr_vec[ok])
  t <- as.numeric(time_vec[ok])
  ev <- as.integer(event_vec[ok])
  
  cox <- tryCatch(summary(coxph(Surv(t, ev) ~ e)), error = function(err) NULL)
  if (is.null(cox)) return(NULL)
  hr  <- cox$conf.int[1, "exp(coef)"]
  lci <- cox$conf.int[1, "lower .95"]
  uci <- cox$conf.int[1, "upper .95"]
  cp  <- cox$coefficients[1, "Pr(>|z|)"]
  
  grp <- ifelse(e >= median(e, na.rm = TRUE), "High", "Low")
  km_p <- NA_real_
  if (length(unique(grp)) == 2) {
    sd <- tryCatch(survdiff(Surv(t, ev) ~ grp), error = function(err) NULL)
    if (!is.null(sd)) km_p <- pchisq(sd$chisq, df = 1, lower.tail = FALSE)
  }
  data.table(
    cohort = cohort_name, stratum = stratum, endpoint = endpoint, gene_symbol = gene,
    n_samples = length(e), n_events = sum(ev == 1),
    HR = hr, HR_LCI95 = lci, HR_UCI95 = uci, cox_pval = cp, km_logrank_pval = km_p
  )
}

# ==============================================================================
# 1. TCGA-BRCA PanCancer Atlas (本地复用 + 内存安全 gc())
# ==============================================================================
cat("\n=== [Cohort 1] Loading TCGA-BRCA PanCancer Atlas ===\n")
tcga_dir <- file.path(bulk_dir, "brca_tcga_pan_can_atlas_2018")

read_rsem_clean <- function(path) {
  dt <- fread(path)
  sym <- trimws(as.character(dt$Hugo_Symbol))
  keep <- !is.na(sym) & nzchar(sym)
  dt <- dt[keep]
  sym <- sym[keep]
  mat <- as.matrix(dt[, -c(1, 2), with = FALSE])
  storage.mode(mat) <- "numeric"
  # 去重基因取均值
  if (any(duplicated(sym))) {
    dt_num <- as.data.table(mat)
    dt_num[, gene := sym]
    dt_agg <- dt_num[, lapply(.SD, mean, na.rm = TRUE), by = gene]
    sym <- dt_agg$gene
    mat <- as.matrix(dt_agg[, -"gene", with = FALSE])
    rm(dt_num, dt_agg)
  }
  rownames(mat) <- sym
  mat[mat < 0] <- NA_real_
  mat <- log2(mat + 1)
  rm(dt); gc()
  mat
}

tcga_tumor  <- read_rsem_clean(file.path(tcga_dir, "data_mrna_seq_v2_rsem.txt"))
tcga_normal <- read_rsem_clean(file.path(tcga_dir, "normals/data_mrna_seq_v2_rsem_normal_samples.txt"))

# 临床表与 Xena ER/PR/HER2 表
clin_pat  <- fread(file.path(tcga_dir, "data_clinical_patient.txt"), skip = 4)
clin_samp <- fread(file.path(tcga_dir, "data_clinical_sample.txt"), skip = 4)
xena_clin <- fread(file.path(bulk_dir, "TCGA_BRCA_Xena_clinicalMatrix.tsv"))

tcga_clin <- merge(clin_samp, clin_pat, by = "PATIENT_ID", all.x = TRUE, suffixes = c(".samp", ".pat"))
tcga_clin <- merge(
  tcga_clin,
  xena_clin[, .(
    SAMPLE_ID = sampleID,
    ER_Xena = breast_carcinoma_estrogen_receptor_status,
    PR_Xena = breast_carcinoma_progesterone_receptor_status,
    HER2_Xena = lab_proc_her2_neu_immunohistochemistry_receptor_status,
    PAM50_Xena = PAM50Call_RNAseq
  )],
  by = "SAMPLE_ID", all.x = TRUE
)
tcga_clin[, OS_event := parse_status(OS_STATUS)]
tcga_clin[, DSS_event := parse_status(DSS_STATUS)]
tcga_clin[, PFS_event := parse_status(PFS_STATUS)]
tcga_clin[, DFS_event := parse_status(DFS_STATUS)]
tcga_clin[, ER_clean := fifelse(ER_Xena %in% c("Positive", "Negative"), ER_Xena,
                        fifelse(SUBTYPE %in% c("BRCA_LumA", "BRCA_LumB"), "Positive_inferred",
                        fifelse(SUBTYPE == "BRCA_Basal", "Negative_inferred", NA_character_)))]

# 保存轻量肿瘤矩阵与临床数据供模块 5/6 直接复用
tcga_tumor_sub <- tcga_tumor[, colnames(tcga_tumor) %in% tcga_clin$SAMPLE_ID, drop = FALSE]
saveRDS(list(tumor_log2 = tcga_tumor_sub, normal_log2 = tcga_normal, clin = tcga_clin),
        file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
cat("TCGA Tumor samples:", ncol(tcga_tumor_sub), "| Normal samples:", ncol(tcga_normal), "\n")

diff_list <- list()
surv_list <- list()
subtype_list <- list()

# 1A. TCGA Tumor vs Normal (全部 + 配对)
pat_t <- substr(colnames(tcga_tumor_sub), 1, 12)
pat_n <- substr(colnames(tcga_normal), 1, 12)
common_pats <- intersect(pat_t, pat_n)

for (g in target_genes) {
  if (!(g %in% rownames(tcga_tumor_sub))) next
  vt <- as.numeric(tcga_tumor_sub[g, ])
  vn <- if (g %in% rownames(tcga_normal)) as.numeric(tcga_normal[g, ]) else numeric(0)
  
  if (length(vn) > 5) {
    wt <- wilcox.test(vt, vn)
    # 配对检验
    paired_p <- NA_real_; paired_fc <- NA_real_
    if (length(common_pats) >= 10) {
      vt_p <- sapply(common_pats, function(p) mean(vt[pat_t == p], na.rm = TRUE))
      vn_p <- sapply(common_pats, function(p) mean(vn[pat_n == p], na.rm = TRUE))
      paired_p  <- wilcox.test(vt_p, vn_p, paired = TRUE)$p.value
      paired_fc <- mean(vt_p - vn_p, na.rm = TRUE)
    }
    diff_list[[paste("TCGA_TN", g)]] <- data.table(
      cohort = "TCGA_BRCA", comparison = "Tumor_vs_Normal", gene_symbol = g,
      n_group1 = sum(!is.na(vt)), n_group2 = sum(!is.na(vn)),
      mean_group1 = round(mean(vt, na.rm = TRUE), 3),
      mean_group2 = round(mean(vn, na.rm = TRUE), 3),
      log2FC = round(mean(vt, na.rm = TRUE) - mean(vn, na.rm = TRUE), 3),
      wilcox_pval = wt$p.value,
      paired_n = length(common_pats), paired_log2FC = round(paired_fc, 3), paired_pval = paired_p
    )
  }
  
  # ER+ vs ER- (用 Xena 精确 IHC 状态)
  idx_erp <- which(colnames(tcga_tumor_sub) %in% tcga_clin[ER_Xena == "Positive"]$SAMPLE_ID)
  idx_ern <- which(colnames(tcga_tumor_sub) %in% tcga_clin[ER_Xena == "Negative"]$SAMPLE_ID)
  if (length(idx_erp) > 10 && length(idx_ern) > 10) {
    wt_er <- wilcox.test(vt[idx_erp], vt[idx_ern])
    diff_list[[paste("TCGA_ER", g)]] <- data.table(
      cohort = "TCGA_BRCA", comparison = "ERpos_vs_ERneg", gene_symbol = g,
      n_group1 = length(idx_erp), n_group2 = length(idx_ern),
      mean_group1 = round(mean(vt[idx_erp], na.rm = TRUE), 3),
      mean_group2 = round(mean(vt[idx_ern], na.rm = TRUE), 3),
      log2FC = round(mean(vt[idx_erp], na.rm = TRUE) - mean(vt[idx_ern], na.rm = TRUE), 3),
      wilcox_pval = wt_er$p.value,
      paired_n = NA_integer_, paired_log2FC = NA_real_, paired_pval = NA_real_
    )
  }
  
  # PAM50 亚型均数
  sub_df <- data.table(SAMPLE_ID = colnames(tcga_tumor_sub), expr = vt)
  sub_df <- merge(sub_df, tcga_clin[, .(SAMPLE_ID, SUBTYPE, ER_Xena)], by = "SAMPLE_ID")
  sub_sm <- sub_df[!is.na(SUBTYPE) & nzchar(SUBTYPE), .(
    cohort = "TCGA_BRCA", gene_symbol = g, n = .N,
    mean_expr = round(mean(expr, na.rm = TRUE), 3),
    sd_expr = round(sd(expr, na.rm = TRUE), 3)
  ), by = SUBTYPE]
  subtype_list[[paste("TCGA", g)]] <- sub_sm
  
  # 生存分析: Overall / ER+ / ER- across OS, DSS, PFS
  cl_ord <- tcga_clin[match(colnames(tcga_tumor_sub), tcga_clin$SAMPLE_ID)]
  strata_idx <- list(
    Overall = seq_len(nrow(cl_ord)),
    ER_pos  = which(cl_ord$ER_Xena == "Positive"),
    ER_neg  = which(cl_ord$ER_Xena == "Negative")
  )
  for (st in names(strata_idx)) {
    ii <- strata_idx[[st]]
    surv_list[[paste("TCGA", g, st, "OS")]]  <- run_survival_stats(vt[ii], cl_ord$OS_MONTHS[ii], cl_ord$OS_event[ii], "TCGA_BRCA", "OS", st, g)
    surv_list[[paste("TCGA", g, st, "DSS")]] <- run_survival_stats(vt[ii], cl_ord$DSS_MONTHS[ii], cl_ord$DSS_event[ii], "TCGA_BRCA", "DSS", st, g)
    surv_list[[paste("TCGA", g, st, "PFS")]] <- run_survival_stats(vt[ii], cl_ord$PFS_MONTHS[ii], cl_ord$PFS_event[ii], "TCGA_BRCA", "PFS", st, g)
  }
}

rm(tcga_tumor, tcga_tumor_sub, tcga_normal); gc()

# ==============================================================================
# 2. METABRIC (n ~ 1,980)
# ==============================================================================
cat("\n=== [Cohort 2] Loading METABRIC ===\n")
met_dir <- file.path(bulk_dir, "brca_metabric")
met_dt  <- fread(file.path(met_dir, "data_mrna_illumina_microarray.txt"))
met_sym <- trimws(as.character(met_dt$Hugo_Symbol))
met_keep <- !is.na(met_sym) & nzchar(met_sym) & (met_sym %in% target_genes)
# 省内存：只保留 8 个目标基因的表达行！
met_sub <- met_dt[met_keep]
met_mat <- as.matrix(met_sub[, -c(1, 2), with = FALSE])
storage.mode(met_mat) <- "numeric"
rownames(met_mat) <- met_sub$Hugo_Symbol
if (any(duplicated(rownames(met_mat)))) {
  tmp_dt <- as.data.table(met_mat)
  tmp_dt[, gene := rownames(met_mat)]
  tmp_agg <- tmp_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene]
  met_mat <- as.matrix(tmp_agg[, -"gene", with = FALSE])
  rownames(met_mat) <- tmp_agg$gene
}
rm(met_dt, met_sub); gc()

met_pat  <- fread(file.path(met_dir, "data_clinical_patient.txt"), skip = 4)
met_samp <- fread(file.path(met_dir, "data_clinical_sample.txt"), skip = 4)
met_clin <- merge(met_samp, met_pat, by = "PATIENT_ID", all.x = TRUE)
met_clin[, OS_event := parse_status(OS_STATUS)]
met_clin[, RFS_event := parse_status(RFS_STATUS)]
met_clin[, ER_clean := fifelse(grepl("pos", ER_STATUS, ignore.case = TRUE) | grepl("pos", ER_IHC, ignore.case = TRUE), "Positive",
                       fifelse(grepl("neg", ER_STATUS, ignore.case = TRUE) | grepl("neg", ER_IHC, ignore.case = TRUE), "Negative", NA_character_))]

# 匹配样本顺序
common_mb <- intersect(colnames(met_mat), met_clin$SAMPLE_ID)
met_mat <- met_mat[, common_mb, drop = FALSE]
met_cl  <- met_clin[match(common_mb, met_clin$SAMPLE_ID)]
saveRDS(list(expr_target = met_mat, clin = met_cl), file.path(bulk_dir, "metabric_target_cache.rds"))
cat("METABRIC matched samples:", length(common_mb), "\n")

for (g in target_genes) {
  if (!(g %in% rownames(met_mat))) next
  vm <- as.numeric(met_mat[g, ])
  
  idx_erp <- which(met_cl$ER_clean == "Positive")
  idx_ern <- which(met_cl$ER_clean == "Negative")
  if (length(idx_erp) > 10 && length(idx_ern) > 10) {
    wt_er <- wilcox.test(vm[idx_erp], vm[idx_ern])
    diff_list[[paste("METABRIC_ER", g)]] <- data.table(
      cohort = "METABRIC", comparison = "ERpos_vs_ERneg", gene_symbol = g,
      n_group1 = length(idx_erp), n_group2 = length(idx_ern),
      mean_group1 = round(mean(vm[idx_erp], na.rm = TRUE), 3),
      mean_group2 = round(mean(vm[idx_ern], na.rm = TRUE), 3),
      log2FC = round(mean(vm[idx_erp], na.rm = TRUE) - mean(vm[idx_ern], na.rm = TRUE), 3),
      wilcox_pval = wt_er$p.value,
      paired_n = NA_integer_, paired_log2FC = NA_real_, paired_pval = NA_real_
    )
  }
  
  strata_mb <- list(
    Overall = seq_len(nrow(met_cl)),
    ER_pos  = idx_erp,
    ER_neg  = idx_ern
  )
  for (st in names(strata_mb)) {
    ii <- strata_mb[[st]]
    surv_list[[paste("METABRIC", g, st, "OS")]]  <- run_survival_stats(vm[ii], met_cl$OS_MONTHS[ii], met_cl$OS_event[ii], "METABRIC", "OS", st, g)
    surv_list[[paste("METABRIC", g, st, "RFS")]] <- run_survival_stats(vm[ii], met_cl$RFS_MONTHS[ii], met_cl$RFS_event[ii], "METABRIC", "RFS", st, g)
  }
}
rm(met_mat); gc()

# ==============================================================================
# 3. 外部 GEO 队列验证 (GSE20685 + GSE42568, 安全 tryCatch 防卡死)
# ==============================================================================
geo_status <- list()

fetch_geo_safe <- function(gse_id) {
  cat("\n=== [GEO] Fetching", gse_id, "===\n")
  cache_rds <- file.path(geo_dir, paste0(gse_id, "_parsed.rds"))
  if (file.exists(cache_rds)) return(readRDS(cache_rds))
  
  gset <- tryCatch({
    res <- GEOquery::getGEO(gse_id, destdir = geo_dir, GSEMatrix = TRUE, getGPL = TRUE)
    if (is.list(res)) res[[1]] else res
  }, error = function(e) {
    message("GEO download/parse failed for ", gse_id, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(gset)) return(NULL)
  
  ex <- Biobase::exprs(gset)
  # 若非 log2 尺度则做 log2(x+1)
  qx <- quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  if ((qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)) {
    ex[ex < 0] <- 0
    ex <- log2(ex + 1)
  }
  fd <- Biobase::fData(gset)
  pd <- Biobase::pData(gset)
  
  # 找 Gene Symbol 列
  sym_col <- intersect(c("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_symbol"), names(fd))[1]
  if (is.na(sym_col)) {
    sym_col <- names(fd)[grepl("symbol", names(fd), ignore.case = TRUE)][1]
  }
  if (!is.na(sym_col)) {
    syms <- trimws(sapply(strsplit(as.character(fd[[sym_col]]), "///|;"), `[`, 1))
    keep <- syms %in% target_genes
    ex_sub <- ex[keep, , drop = FALSE]
    syms_sub <- syms[keep]
    if (nrow(ex_sub) > 0) {
      dt_e <- as.data.table(ex_sub)
      dt_e[, gene := syms_sub]
      dt_agg <- dt_e[, lapply(.SD, mean, na.rm = TRUE), by = gene]
      mat_target <- as.matrix(dt_agg[, -"gene", with = FALSE])
      rownames(mat_target) <- dt_agg$gene
    } else {
      mat_target <- matrix(nrow = 0, ncol = ncol(ex))
    }
  } else {
    mat_target <- matrix(nrow = 0, ncol = ncol(ex))
  }
  out <- list(expr_target = mat_target, pheno = as.data.table(pd))
  saveRDS(out, cache_rds)
  rm(gset, ex); gc()
  out
}

# 3A. GSE42568 (104 Tumor + 17 Normal, 含 OS/RFS)
g42 <- fetch_geo_safe("GSE42568")
if (!is.null(g42) && nrow(g42$expr_target) > 0) {
  pd42 <- g42$pheno
  # 自动识别 Tumor/Normal 与生存列
  char_txt <- apply(pd42[, grepl("title|source|characteristics", names(pd42), ignore.case = TRUE), with = FALSE], 1, paste, collapse = " | ")
  is_norm <- grepl("normal", char_txt, ignore.case = TRUE)
  is_tum  <- !is_norm
  
  # 提取 OS / RFS
  find_num_col <- function(df, pat) {
    cn <- names(df)[grepl(pat, names(df), ignore.case = TRUE)]
    if (length(cn) == 0) return(rep(NA_real_, nrow(df)))
    suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(df[[cn[1]]]))))
  }
  os_t42 <- find_num_col(pd42, "overall survival.*time|os.*time|survival.*time")
  os_e42 <- find_num_col(pd42, "overall survival.*event|os.*event|os_status|status.*overall")
  rfs_t42 <- find_num_col(pd42, "relapse.*time|rfs.*time|recurrence.*time")
  rfs_e42 <- find_num_col(pd42, "relapse.*event|rfs.*event|recurrence.*event")
  
  for (g in rownames(g42$expr_target)) {
    v <- as.numeric(g42$expr_target[g, ])
    if (sum(is_tum) > 10 && sum(is_norm) > 5) {
      wt <- wilcox.test(v[is_tum], v[is_norm])
      diff_list[[paste("GSE42568_TN", g)]] <- data.table(
        cohort = "GSE42568", comparison = "Tumor_vs_Normal", gene_symbol = g,
        n_group1 = sum(is_tum), n_group2 = sum(is_norm),
        mean_group1 = round(mean(v[is_tum], na.rm = TRUE), 3),
        mean_group2 = round(mean(v[is_norm], na.rm = TRUE), 3),
        log2FC = round(mean(v[is_tum], na.rm = TRUE) - mean(v[is_norm], na.rm = TRUE), 3),
        wilcox_pval = wt$p.value,
        paired_n = NA_integer_, paired_log2FC = NA_real_, paired_pval = NA_real_
      )
    }
    surv_list[[paste("GSE42568", g, "OS")]]  <- run_survival_stats(v[is_tum], os_t42[is_tum], os_e42[is_tum], "GSE42568", "OS", "Overall", g)
    surv_list[[paste("GSE42568", g, "RFS")]] <- run_survival_stats(v[is_tum], rfs_t42[is_tum], rfs_e42[is_tum], "GSE42568", "RFS", "Overall", g)
  }
  geo_status[["GSE42568"]] <- data.table(cohort = "GSE42568", status = "SUCCESS", n_samples = ncol(g42$expr_target), n_genes_matched = nrow(g42$expr_target))
} else {
  geo_status[["GSE42568"]] <- data.table(cohort = "GSE42568", status = "FAILED_OR_EMPTY", n_samples = 0, n_genes_matched = 0)
}

# 3B. GSE20685 (n=327 Breast cancer with OS & MFS)
g20 <- fetch_geo_safe("GSE20685")
if (!is.null(g20) && nrow(g20$expr_target) > 0) {
  pd20 <- g20$pheno
  find_num_col <- function(df, pat) {
    cn <- names(df)[grepl(pat, names(df), ignore.case = TRUE)]
    if (length(cn) == 0) return(rep(NA_real_, nrow(df)))
    suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(df[[cn[1]]]))))
  }
  os_t20 <- find_num_col(pd20, "follow_up_time|os.*time|survival.*time")
  os_e20 <- find_num_col(pd20, "event_death|os.*event|death")
  mfs_t20 <- find_num_col(pd20, "metastasis.*time|mfs.*time|dmfs.*time")
  mfs_e20 <- find_num_col(pd20, "event_metastasis|metastasis.*event")
  
  for (g in rownames(g20$expr_target)) {
    v <- as.numeric(g20$expr_target[g, ])
    surv_list[[paste("GSE20685", g, "OS")]]  <- run_survival_stats(v, os_t20, os_e20, "GSE20685", "OS", "Overall", g)
    surv_list[[paste("GSE20685", g, "MFS")]] <- run_survival_stats(v, mfs_t20, mfs_e20, "GSE20685", "MFS", "Overall", g)
  }
  geo_status[["GSE20685"]] <- data.table(cohort = "GSE20685", status = "SUCCESS", n_samples = ncol(g20$expr_target), n_genes_matched = nrow(g20$expr_target))
} else {
  geo_status[["GSE20685"]] <- data.table(cohort = "GSE20685", status = "FAILED_OR_EMPTY", n_samples = 0, n_genes_matched = 0)
}

# ==============================================================================
# 4. 汇总保存 & 自动进行严格方向一致性审计
# ==============================================================================
diff_dt    <- rbindlist(diff_list, fill = TRUE)
surv_dt    <- rbindlist(surv_list, fill = TRUE)
subtype_dt <- rbindlist(subtype_list, fill = TRUE)
geo_st_dt  <- rbindlist(geo_status, fill = TRUE)

fwrite(diff_dt, file.path(res_dir, "module4_differential_expression.csv"))
fwrite(surv_dt, file.path(res_dir, "module4_survival_analysis.csv"))
fwrite(subtype_dt, file.path(res_dir, "module4_tcga_pam50_subtypes.csv"))
fwrite(geo_st_dt, file.path(res_dir, "module4_geo_status.csv"))

# 构建每个基因的多关口严格方向一致性汇总表
coloc_dt <- fread(file.path(res_dir, "module3_cis_positive_genes_coloc.csv"))
audit_rows <- list()
for (g in target_genes) {
  mr_sub <- sig_cis[gene_symbol == g][order(pval)]
  best_mr <- mr_sub[1]
  col_sub <- coloc_dt[gene_symbol == g]
  max_h4 <- if (nrow(col_sub) > 0) max(col_sub$PP.H4, na.rm = TRUE) else NA_real_
  
  tcga_tn <- diff_dt[cohort == "TCGA_BRCA" & comparison == "Tumor_vs_Normal" & gene_symbol == g]
  tcga_er <- diff_dt[cohort == "TCGA_BRCA" & comparison == "ERpos_vs_ERneg" & gene_symbol == g]
  gse_tn  <- diff_dt[cohort == "GSE42568" & comparison == "Tumor_vs_Normal" & gene_symbol == g]
  
  tcga_os  <- surv_dt[cohort == "TCGA_BRCA" & stratum == "Overall" & endpoint == "OS" & gene_symbol == g]
  tcga_dss <- surv_dt[cohort == "TCGA_BRCA" & stratum == "Overall" & endpoint == "DSS" & gene_symbol == g]
  mb_os    <- surv_dt[cohort == "METABRIC" & stratum == "Overall" & endpoint == "OS" & gene_symbol == g]
  mb_rfs   <- surv_dt[cohort == "METABRIC" & stratum == "Overall" & endpoint == "RFS" & gene_symbol == g]
  g20_os   <- surv_dt[cohort == "GSE20685" & endpoint == "OS" & gene_symbol == g]
  
  mr_dir <- ifelse(best_mr$OR > 1, "Risk (OR>1)", "Protective (OR<1)")
  tn_dir <- if (nrow(tcga_tn) > 0) ifelse(tcga_tn$log2FC > 0, "Up in Tumor", "Down in Tumor") else NA_character_
  # MR与肿瘤表达方向是否一致：Risk基因应在肿瘤高表达，Protective基因应在肿瘤低表达
  mr_vs_expr_match <- if (!is.na(tn_dir)) {
    (best_mr$OR > 1 && tcga_tn$log2FC > 0) || (best_mr$OR < 1 && tcga_tn$log2FC < 0)
  } else NA
  
  audit_rows[[g]] <- data.table(
    gene_symbol = g,
    upstream_metabolites = best_mr$upstream_metabolites,
    gutmgene_alteration = best_mr$alterations,
    best_MR_outcome = best_mr$outcome_label,
    MR_OR = round(best_mr$OR, 3),
    MR_pval = signif(best_mr$pval, 3),
    MR_FDR = signif(best_mr$FDR, 3),
    coloc_max_PP.H4 = round(max_h4, 4),
    TCGA_TN_log2FC = if (nrow(tcga_tn) > 0) tcga_tn$log2FC else NA_real_,
    TCGA_TN_pval = if (nrow(tcga_tn) > 0) signif(tcga_tn$wilcox_pval, 3) else NA_real_,
    GSE42568_TN_log2FC = if (nrow(gse_tn) > 0) gse_tn$log2FC else NA_real_,
    GSE42568_TN_pval = if (nrow(gse_tn) > 0) signif(gse_tn$wilcox_pval, 3) else NA_real_,
    MR_vs_TCGA_TN_consistent = mr_vs_expr_match,
    TCGA_ER_log2FC = if (nrow(tcga_er) > 0) tcga_er$log2FC else NA_real_,
    TCGA_OS_HR = if (nrow(tcga_os) > 0) round(tcga_os$HR, 3) else NA_real_,
    TCGA_OS_pval = if (nrow(tcga_os) > 0) signif(tcga_os$cox_pval, 3) else NA_real_,
    METABRIC_OS_HR = if (nrow(mb_os) > 0) round(mb_os$HR, 3) else NA_real_,
    METABRIC_OS_pval = if (nrow(mb_os) > 0) signif(mb_os$cox_pval, 3) else NA_real_,
    METABRIC_RFS_HR = if (nrow(mb_rfs) > 0) round(mb_rfs$HR, 3) else NA_real_,
    METABRIC_RFS_pval = if (nrow(mb_rfs) > 0) signif(mb_rfs$cox_pval, 3) else NA_real_,
    GSE20685_OS_HR = if (nrow(g20_os) > 0) round(g20_os$HR, 3) else NA_real_,
    GSE20685_OS_pval = if (nrow(g20_os) > 0) signif(g20_os$cox_pval, 3) else NA_real_
  )
}

audit_dt <- rbindlist(audit_rows)
fwrite(audit_dt, file.path(res_dir, "module4_direction_audit_summary.csv"))

cat("\n=== 模块 4 多队列验证与方向一致性审计总表 ===\n")
print(audit_dt[, .(gene_symbol, best_MR_outcome, MR_OR, MR_pval, MR_vs_TCGA_TN_consistent,
                   TCGA_TN_log2FC, TCGA_TN_pval, GSE42568_TN_log2FC, GSE42568_TN_pval,
                   TCGA_OS_HR, TCGA_OS_pval, METABRIC_OS_HR, METABRIC_OS_pval, GSE20685_OS_HR, GSE20685_OS_pval)])
cat("\n=== GEO 队列状态 ===\n")
print(geo_st_dt)
