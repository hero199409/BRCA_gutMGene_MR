# 12_module5_tme_gsea_local.R
# 模块 5 本地分析：TCGA-BRCA (1,082 例) 肿瘤免疫微环境 (quanTIseq, MCP-counter, ESTIMATE)、
# 免疫检查点/TMB/MSI 相关性、单基因 Hallmark GSEA 富集，以及 TISCH2/HPA 网页待查清单

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(immunedeconv)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(jsonlite)
})

bulk_dir <- "00_data/bulk_cohorts"
res_dir  <- "02_results"

tc <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
expr_log2 <- tc$tumor_log2
clin      <- tc$clin[match(colnames(expr_log2), tc$clin$SAMPLE_ID)]
rm(tc); gc()

# 确保无 NA / 负值并转为线性尺度供 immunedeconv 使用
expr_clean <- expr_log2
expr_clean[is.na(expr_clean)] <- 0
expr_lin <- 2^expr_clean - 1
expr_lin[expr_lin < 0] <- 0

# 1. 免疫浸润与 ESTIMATE 评分 (带缓存 + tryCatch 防异常)
tme_cache <- file.path(bulk_dir, "tcga_tme_scores.rds")
if (file.exists(tme_cache)) {
  tme_res <- readRDS(tme_cache)
} else {
  cat("=== Running immunedeconv: quantiseq, mcp_counter, estimate ===\n")
  qs  <- tryCatch(immunedeconv::deconvolute(expr_lin, method = "quantiseq", tumor = TRUE, arrays = FALSE, scale_mrna = TRUE),
                  error = function(e) { message("quantiseq error: ", conditionMessage(e)); NULL })
  mcp <- tryCatch(immunedeconv::deconvolute(expr_lin, method = "mcp_counter"),
                  error = function(e) { message("mcp_counter error: ", conditionMessage(e)); NULL })
  est <- tryCatch(immunedeconv::deconvolute(expr_lin, method = "estimate"),
                  error = function(e) { message("estimate error: ", conditionMessage(e)); NULL })
  tme_res <- list(quantiseq = qs, mcp_counter = mcp, estimate = est)
  saveRDS(tme_res, tme_cache)
}
rm(expr_lin); gc()

# 整理为样本 x 特征矩阵
to_sample_mat <- function(df, prefix) {
  if (is.null(df)) return(NULL)
  df <- as.data.frame(df)
  rn <- paste0(prefix, "__", df$cell_type)
  m <- as.matrix(df[, -1, drop = FALSE])
  rownames(m) <- rn
  m
}

m_qs  <- to_sample_mat(tme_res$quantiseq, "quanTIseq")
m_mcp <- to_sample_mat(tme_res$mcp_counter, "MCPcounter")
m_est <- to_sample_mat(tme_res$estimate, "ESTIMATE")
tme_mat <- do.call(rbind, Filter(Negate(is.null), list(m_qs, m_mcp, m_est)))
cat("TME features extracted:", nrow(tme_mat), "across", ncol(tme_mat), "samples\n")

# 2. 核心基因与免疫细胞、免疫检查点、TMB、MSI 的 Spearman 相关性 (Overall / ER+ / ER-)
focus_genes <- c("PPARG", "CXCL8", "IL1B", "FFAR2", "HDAC11", "TNF", "IL12A")
checkpoints <- c("CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD80", "CD86", "IDO1", "FOXP3")
avail_cp <- intersect(checkpoints, rownames(expr_log2))

strata_list <- list(
  Overall = seq_len(ncol(expr_log2)),
  ER_pos  = which(clin$ER_Xena == "Positive"),
  ER_neg  = which(clin$ER_Xena == "Negative")
)

cor_rows <- list()
for (st in names(strata_list)) {
  ii <- strata_list[[st]]
  sids <- colnames(expr_log2)[ii]
  for (g in focus_genes) {
    if (!(g %in% rownames(expr_log2))) next
    gv <- as.numeric(expr_log2[g, ii])
    
    # 免疫细胞与 ESTIMATE
    for (feat in rownames(tme_mat)) {
      fv <- as.numeric(tme_mat[feat, sids])
      ok <- !is.na(gv) & !is.na(fv) & sd(fv[ok = !is.na(gv) & !is.na(fv)]) > 0
      if (sum(ok) >= 20 && sd(gv[ok]) > 0 && sd(fv[ok]) > 0) {
        ct <- suppressWarnings(cor.test(gv[ok], fv[ok], method = "spearman"))
        cor_rows[[paste(st, g, feat)]] <- data.table(
          stratum = st, gene_symbol = g, feature_type = "TME_Cell_Score",
          feature = feat, n = sum(ok), rho = unname(ct$estimate), pval = ct$p.value
        )
      }
    }
    # 免疫检查点
    for (cp in avail_cp) {
      cv <- as.numeric(expr_log2[cp, ii])
      ok <- !is.na(gv) & !is.na(cv)
      if (sum(ok) >= 20) {
        ct <- suppressWarnings(cor.test(gv[ok], cv[ok], method = "spearman"))
        cor_rows[[paste(st, g, "CP", cp)]] <- data.table(
          stratum = st, gene_symbol = g, feature_type = "Immune_Checkpoint",
          feature = cp, n = sum(ok), rho = unname(ct$estimate), pval = ct$p.value
        )
      }
    }
    # TMB 与 MSI
    for (gfeat in c("TMB_NONSYNONYMOUS", "MSI_SENSOR_SCORE")) {
      mv <- suppressWarnings(as.numeric(clin[[gfeat]][ii]))
      ok <- !is.na(gv) & !is.na(mv)
      if (sum(ok) >= 20 && sd(mv[ok]) > 0) {
        ct <- suppressWarnings(cor.test(gv[ok], mv[ok], method = "spearman"))
        cor_rows[[paste(st, g, gfeat)]] <- data.table(
          stratum = st, gene_symbol = g, feature_type = "Genomic_Feature",
          feature = gfeat, n = sum(ok), rho = unname(ct$estimate), pval = ct$p.value
        )
      }
    }
  }
}

cor_dt <- rbindlist(cor_rows)
cor_dt[, FDR := p.adjust(pval, method = "BH"), by = .(stratum, gene_symbol, feature_type)]
fwrite(cor_dt[order(stratum, gene_symbol, feature_type, pval)], file.path(res_dir, "module5_tme_checkpoint_correlations.csv"))

cat("\n=== PPARG (ER+) & CXCL8 (ER-) 与关键免疫浸润相关性 Top 10 ===\n")
print(cor_dt[(gene_symbol == "PPARG" & stratum == "ER_pos") | (gene_symbol == "CXCL8" & stratum == "ER_neg")][
  feature_type == "TME_Cell_Score"][order(gene_symbol, pval)][, .SD[1:min(8, .N)], by = gene_symbol][
  , .(stratum, gene_symbol, feature, rho = round(rho, 3), pval = signif(pval, 3), FDR = signif(FDR, 3))])

# 3. 单基因 Hallmark GSEA (PPARG in Overall & ER+, CXCL8 in Overall & ER-, IL1B/FFAR2 in ER-)
cat("\n=== Running single-gene Hallmark GSEA ===\n")
gmt <- clusterProfiler::read.gmt(file.path(bulk_dir, "hallmark_human_entrez.gmt"))
# 映射表达矩阵基因名到 Entrez ID (预先映射一次省时间)
sym_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = rownames(expr_clean),
                                                  columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL"))
sym_map <- as.data.table(sym_map)[!is.na(ENTREZID) & !duplicated(SYMBOL)]

run_single_gene_gsea <- function(gene, stratum_name, idx_cols) {
  if (!(gene %in% rownames(expr_clean))) return(NULL)
  sub_mat <- expr_clean[, idx_cols, drop = FALSE]
  gv <- as.numeric(sub_mat[gene, ])
  grp <- factor(ifelse(gv >= median(gv, na.rm = TRUE), "High", "Low"), levels = c("Low", "High"))
  design <- model.matrix(~ grp)
  fit <- limma::lmFit(sub_mat, design)
  fit <- limma::eBayes(fit)
  tt  <- limma::topTable(fit, coef = "grpHigh", number = Inf, sort.by = "none")
  tt$SYMBOL <- rownames(tt)
  tt_dt <- merge(as.data.table(tt), sym_map, by = "SYMBOL")
  tt_dt <- tt_dt[!is.na(logFC) & !duplicated(ENTREZID)]
  # 按 logFC 降序 (若相等按 t 统计量打破平局)
  tt_dt <- tt_dt[order(-logFC, -t)]
  gene_list <- setNames(tt_dt$logFC, tt_dt$ENTREZID)
  
  gsea_res <- tryCatch({
    set.seed(2026)
    clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = gmt,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1.0,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
  }, error = function(e) { message("GSEA error for ", gene, ": ", conditionMessage(e)); NULL })
  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) return(NULL)
  res_df <- as.data.table(as.data.frame(gsea_res))
  res_df[, `:=`(target_gene = gene, stratum = stratum_name)]
  res_df
}

gsea_tasks <- list(
  list(g = "PPARG", st = "Overall", idx = strata_list$Overall),
  list(g = "PPARG", st = "ER_pos",  idx = strata_list$ER_pos),
  list(g = "CXCL8", st = "Overall", idx = strata_list$Overall),
  list(g = "CXCL8", st = "ER_neg",  idx = strata_list$ER_neg),
  list(g = "IL1B",  st = "ER_neg",  idx = strata_list$ER_neg),
  list(g = "FFAR2", st = "ER_neg",  idx = strata_list$ER_neg)
)

gsea_all <- list()
for (tk in gsea_tasks) {
  cat("  GSEA:", tk$g, "(", tk$st, ")\n")
  out <- run_single_gene_gsea(tk$g, tk$st, tk$idx)
  if (!is.null(out)) gsea_all[[paste(tk$g, tk$st)]] <- out
}
gsea_dt <- rbindlist(gsea_all, fill = TRUE)
fwrite(gsea_dt[, .(target_gene, stratum, ID, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue)],
       file.path(res_dir, "module5_hallmark_gsea_results.csv"))

cat("\n=== Top Hallmark 通路 (PPARG ER+ & CXCL8 ER-) ===\n")
print(gsea_dt[(target_gene == "PPARG" & stratum == "ER_pos") | (target_gene == "CXCL8" & stratum == "ER_neg")][
  p.adjust < 0.05][order(target_gene, -abs(NES))][, .SD[1:min(6, .N)], by = target_gene][
  , .(target_gene, stratum, ID, NES = round(NES, 2), p.adjust = signif(p.adjust, 3))])

rm(expr_clean, expr_log2); gc()

# 4. 测试 HPA JSON 连通性 & 输出网页端(TISCH2/HPA)操作指引清单 (不卡死)
ens_check <- c(PPARG = "ENSG00000132170", CXCL8 = "ENSG00000169429", IL1B = "ENSG00000125538", FFAR2 = "ENSG00000126262")
hpa_rows <- list()
for (sym in names(ens_check)) {
  eid <- ens_check[[sym]]
  url <- sprintf("https://www.proteinatlas.org/%s.json", eid)
  cmd <- sprintf('curl -4 -s -w "HTTP:%%{http_code}" --connect-timeout 5 --max-time 8 "%s"', url)
  raw <- tryCatch(paste(suppressWarnings(system(cmd, intern = TRUE)), collapse = ""), error = function(e) "ERR")
  ok_http <- grepl("HTTP:200$", raw)
  rna_norm <- NA_character_; rna_brca <- NA_character_; ihc_brca <- NA_character_
  if (ok_http) {
    js_txt <- sub("HTTP:200$", "", raw)
    js <- tryCatch(jsonlite::fromJSON(js_txt), error = function(e) NULL)
    if (!is.null(js)) {
      rna_brca <- paste(js[["Pathology prognostics - Breast cancer"]], collapse = "; ")
      ihc_brca <- paste(js[["Antibody"]] %||% "", collapse = "")
    }
  }
  hpa_rows[[sym]] <- data.table(
    gene_symbol = sym,
    ensembl_id = eid,
    hpa_api_status = ifelse(ok_http, "HTTP_200_OK", "SKIPPED_TIMEOUT_OR_BLOCKED"),
    hpa_gene_url = sprintf("https://www.proteinatlas.org/%s-%s/pathology/breast+cancer", eid, sym),
    tisch2_url = "http://tisch.comp-genomics.org/gallery/?cancer=BRCA",
    tisch2_recommended_cohorts = "BRCA_GSE176078; BRCA_EMTAB8107; BRCA_GSE114727_inDrop",
    hpa_pathology_note = ifelse(nzchar(rna_brca), rna_brca, "Check web manually")
  )
}
web_dt <- rbindlist(hpa_rows)
fwrite(web_dt, file.path(res_dir, "module5_web_lookup_checklist.csv"))
cat("\n=== HPA & TISCH2 网页查证状态 ===\n")
print(web_dt[, .(gene_symbol, hpa_api_status, hpa_pathology_note)])
