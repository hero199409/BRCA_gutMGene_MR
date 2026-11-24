# 11b_fix_geo_and_cxcl8_audit.R
# 1) 检查 gutMGene 三元组文献摘要 (PPARG, IL1B, FFAR2, CXCL8, TNF)
# 2) 用 Entrez ID 3576 补查 TCGA 中 CXCL8 是否存在
# 3) 补跑 GSE20685 (timeout=360s) 与轻量队列 GSE7390 (n=198 生存验证)

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(GEOquery)
  library(Biobase)
})

options(timeout = 360)

bulk_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
geo_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/geo_cache")
res_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

# 1. 检查 gutMGene 详细文献证据
f5 <- fread(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/gutMGene_v2.0_Human_Microbe_Metabolite_Host_Gene_PubChem_CID_join.csv"))
gs <- c("PPARG", "IL1B", "FFAR2", "CXCL8", "TNF")
ev <- unique(f5[host_gene %in% gs, .(host_gene, gene_table_metabolite, metabolite_gene_alteration,
                                     metabolite_gene_pmids, metabolite_gene_literature_evidence)])
fwrite(ev, file.path(res_dir, "module3_top_genes_gutmgene_literature.csv"))
cat("=== gutMGene 文献证据详情 ===\n")
print(ev)

# 2. 检查 Entrez ID 3576 (CXCL8) 在 TCGA RSEM 中是否存在
hdr_t <- fread(file.path(bulk_dir, "brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt"), select = 1:2)
cat("\nEntrez 3576 matches in TCGA:", nrow(hdr_t[Entrez_Gene_Id == 3576]), "\n")
if (nrow(hdr_t[Entrez_Gene_Id == 3576]) > 0) print(hdr_t[Entrez_Gene_Id == 3576])

# 如果 TCGA 有 Entrez 3576，提取它的表达补跑 TCGA CXCL8
tc <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
if (nrow(hdr_t[Entrez_Gene_Id == 3576]) > 0) {
  row_idx <- which(hdr_t$Entrez_Gene_Id == 3576)[1]
  t_full <- fread(file.path(bulk_dir, "brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt"), skip = row_idx, nrows = 1, header = FALSE)
  cnames_t <- names(fread(file.path(bulk_dir, "brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt"), nrows = 0))
  setnames(t_full, cnames_t)
  
  hdr_n <- fread(file.path(bulk_dir, "brca_tcga_pan_can_atlas_2018/normals/data_mrna_seq_v2_rsem_normal_samples.txt"), select = 1:2)
  row_n <- which(hdr_n$Entrez_Gene_Id == 3576)[1]
  n_full <- fread(file.path(bulk_dir, "brca_tcga_pan_can_atlas_2018/normals/data_mrna_seq_v2_rsem_normal_samples.txt"), skip = row_n, nrows = 1, header = FALSE)
  cnames_n <- names(fread(file.path(bulk_dir, "brca_tcga_pan_can_atlas_2018/normals/data_mrna_seq_v2_rsem_normal_samples.txt"), nrows = 0))
  setnames(n_full, cnames_n)
  
  vt <- log2(as.numeric(t_full[, colnames(tc$tumor_log2), with = FALSE]) + 1)
  vn <- log2(as.numeric(n_full[, colnames(tc$normal_log2), with = FALSE]) + 1)
  cat("TCGA CXCL8 (Entrez 3576) Tumor vs Normal log2FC:", round(mean(vt, na.rm=TRUE) - mean(vn, na.rm=TRUE), 3),
      "p=", signif(wilcox.test(vt, vn)$p.value, 3), "\n")
  
  # 加到 cache 里供模块 5/6 使用
  if (!("CXCL8" %in% rownames(tc$tumor_log2))) {
    tc$tumor_log2  <- rbind(tc$tumor_log2, CXCL8 = vt)
    tc$normal_log2 <- rbind(tc$normal_log2, CXCL8 = vn)
    saveRDS(tc, file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
  }
  # 算 TCGA CXCL8 生存
  cl_ord <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]
  ok <- !is.na(vt) & !is.na(cl_ord$OS_MONTHS) & !is.na(cl_ord$OS_event) & cl_ord$OS_MONTHS > 0
  cx <- summary(coxph(Surv(cl_ord$OS_MONTHS[ok], cl_ord$OS_event[ok]) ~ vt[ok]))
  cat("TCGA CXCL8 Overall OS HR=", round(cx$conf.int[1,"exp(coef)"], 3), "p=", signif(cx$coefficients[1,"Pr(>|z|)"], 3), "\n")
}
rm(tc); gc()

# 3. 补拉取 GSE20685 (n=327) 和 GSE7390 (n=198)
target_genes <- c("IL1B", "FFAR2", "CXCL8", "HDAC11", "CXCR2", "PPARG", "IL12A", "TNF")

fetch_geo_clean <- function(gse_id) {
  cache_rds <- file.path(geo_dir, paste0(gse_id, "_parsed.rds"))
  if (file.exists(cache_rds)) return(readRDS(cache_rds))
  cat("Downloading", gse_id, "...\n")
  gset <- tryCatch({
    res <- GEOquery::getGEO(gse_id, destdir = geo_dir, GSEMatrix = TRUE, getGPL = TRUE)
    if (is.list(res)) res[[1]] else res
  }, error = function(e) { message("  Failed ", gse_id, ": ", conditionMessage(e)); NULL })
  if (is.null(gset)) return(NULL)
  ex <- Biobase::exprs(gset)
  qx <- quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  if ((qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)) {
    ex[ex < 0] <- 0
    ex <- log2(ex + 1)
  }
  fd <- Biobase::fData(gset)
  pd <- Biobase::pData(gset)
  sym_col <- intersect(c("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_symbol"), names(fd))[1]
  if (is.na(sym_col)) sym_col <- names(fd)[grepl("symbol", names(fd), ignore.case = TRUE)][1]
  syms <- trimws(sapply(strsplit(as.character(fd[[sym_col]]), "///|;"), `[`, 1))
  keep <- syms %in% target_genes
  ex_sub <- ex[keep, , drop = FALSE]
  dt_e <- as.data.table(ex_sub); dt_e[, gene := syms[keep]]
  dt_agg <- dt_e[, lapply(.SD, mean, na.rm = TRUE), by = gene]
  mat_t <- as.matrix(dt_agg[, -"gene", with = FALSE]); rownames(mat_t) <- dt_agg$gene
  out <- list(expr_target = mat_t, pheno = as.data.table(pd))
  saveRDS(out, cache_rds)
  rm(gset, ex); gc()
  out
}

g20 <- fetch_geo_clean("GSE20685")
g73 <- fetch_geo_clean("GSE7390")

extra_surv <- list()
run_cox_km <- function(e, t, ev, coh, ep, g) {
  ok <- !is.na(e) & !is.na(t) & !is.na(ev) & t > 0
  if (sum(ok) < 30 || sum(ev[ok] == 1) < 5) return(NULL)
  cx <- summary(coxph(Surv(t[ok], ev[ok]) ~ e[ok]))
  grp <- ifelse(e[ok] >= median(e[ok]), "High", "Low")
  kp <- if (length(unique(grp)) == 2) pchisq(survdiff(Surv(t[ok], ev[ok]) ~ grp)$chisq, df = 1, lower.tail = FALSE) else NA_real_
  data.table(cohort = coh, stratum = "Overall", endpoint = ep, gene_symbol = g,
             n_samples = sum(ok), n_events = sum(ev[ok] == 1),
             HR = cx$conf.int[1,"exp(coef)"], HR_LCI95 = cx$conf.int[1,"lower .95"], HR_UCI95 = cx$conf.int[1,"upper .95"],
             cox_pval = cx$coefficients[1,"Pr(>|z|)"], km_logrank_pval = kp)
}

if (!is.null(g20)) {
  pd <- g20$pheno
  get_num <- function(p) {
    cn <- names(pd)[grepl(p, names(pd), ignore.case = TRUE)]
    if (length(cn) == 0) return(rep(NA_real_, nrow(pd)))
    suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(pd[[cn[1]]]))))
  }
  t_os <- get_num("follow_up_time|os.*time")
  e_os <- get_num("event_death|death")
  for (g in rownames(g20$expr_target)) {
    extra_surv[[paste("GSE20685", g)]] <- run_cox_km(as.numeric(g20$expr_target[g,]), t_os, e_os, "GSE20685", "OS", g)
  }
  cat("GSE20685 parsed successfully! Samples:", ncol(g20$expr_target), "\n")
}

if (!is.null(g73)) {
  pd <- g73$pheno
  get_num <- function(p) {
    cn <- names(pd)[grepl(p, names(pd), ignore.case = TRUE)]
    if (length(cn) == 0) return(rep(NA_real_, nrow(pd)))
    suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(pd[[cn[1]]]))))
  }
  t_os <- get_num("^t\\.os|os\\.time|t\\.dmfs")
  e_os <- get_num("^e\\.os|os\\.event|e\\.dmfs")
  # 打印匹配列名方便核查
  cat("GSE7390 surv cols:", paste(names(pd)[grepl("os|rfs|dmfs", names(pd), ignore.case = TRUE)], collapse = ", "), "\n")
  t_os_col <- names(pd)[grepl("t\\.os", names(pd), ignore.case = TRUE)][1]
  e_os_col <- names(pd)[grepl("e\\.os", names(pd), ignore.case = TRUE)][1]
  t_rfs_col <- names(pd)[grepl("t\\.rfs", names(pd), ignore.case = TRUE)][1]
  e_rfs_col <- names(pd)[grepl("e\\.rfs", names(pd), ignore.case = TRUE)][1]
  if (!is.na(t_os_col) && !is.na(e_os_col)) {
    to <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(pd[[t_os_col]]))))
    eo <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(pd[[e_os_col]]))))
    for (g in rownames(g73$expr_target)) {
      extra_surv[[paste("GSE7390_OS", g)]] <- run_cox_km(as.numeric(g73$expr_target[g,]), to, eo, "GSE7390", "OS", g)
    }
  }
  if (!is.na(t_rfs_col) && !is.na(e_rfs_col)) {
    tr <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(pd[[t_rfs_col]]))))
    er <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(pd[[e_rfs_col]]))))
    for (g in rownames(g73$expr_target)) {
      extra_surv[[paste("GSE7390_RFS", g)]] <- run_cox_km(as.numeric(g73$expr_target[g,]), tr, er, "GSE7390", "RFS", g)
    }
  }
}

extra_dt <- rbindlist(extra_surv, fill = TRUE)
if (nrow(extra_dt) > 0) {
  fwrite(extra_dt, file.path(res_dir, "module4_geo_extra_survival.csv"))
  cat("\n=== 额外 GEO 队列生存验证 (GSE20685 / GSE7390) ===\n")
  print(extra_dt[order(cox_pval), .(cohort, endpoint, gene_symbol, n_samples, n_events, HR = round(HR, 3), cox_pval = signif(cox_pval, 3), km_pval = signif(km_logrank_pval, 3))])
}
