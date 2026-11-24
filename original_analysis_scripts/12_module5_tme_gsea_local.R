# 12_module5_tme_gsea_local.R
# 模块 5 本地分析：TCGA-BRCA 肿瘤免疫微环境 (quanTIseq, MCP-counter, ESTIMATE)、
# 免疫检查点/TMB/MSI 相关性、单基因 Hallmark GSEA 富集，以及 TISCH2/HPA 网页待查清单

# 此脚本始终离线，仅使用本地缓存和已装 R 包。
suppressPackageStartupMessages({
  library(data.table)
  library(immunedeconv)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

bulk_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
res_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

tc <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
expr_log2 <- tc$tumor_log2
sids_tcga <- colnames(expr_log2)
stopifnot(!anyDuplicated(sids_tcga), !anyDuplicated(tc$clin$SAMPLE_ID),
          !anyDuplicated(rownames(expr_log2)), all(sids_tcga %in% tc$clin$SAMPLE_ID))
clin <- tc$clin[match(sids_tcga, tc$clin$SAMPLE_ID)]
stopifnot(identical(sids_tcga, as.character(clin$SAMPLE_ID)))
rm(tc); gc()

# 直接用原始本地 Xena 临床表按表达样本 ID 复核缓存 ER 状态。
xena <- fread(file.path(bulk_dir, "TCGA_BRCA_Xena_clinicalMatrix.tsv"),
              select = c("sampleID", "breast_carcinoma_estrogen_receptor_status"))
stopifnot(!anyDuplicated(xena$sampleID))
xena_er <- xena$breast_carcinoma_estrogen_receptor_status[match(sids_tcga, xena$sampleID)]
normalize_er <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- "Missing"
  x
}
stopifnot(identical(normalize_er(clin$ER_Xena), normalize_er(xena_er)))
er_verification <- data.table(er_status = normalize_er(xena_er))[, .(n_samples = .N), by = er_status]
er_verification <- rbind(data.table(er_status = "Overall", n_samples = length(sids_tcga)), er_verification)
fwrite(er_verification, file.path(res_dir, "module5_tcga_er_sample_counts.csv"))
print(er_verification)

# 还原线性 RSEM。TPM-scale 输入只将每个样本的列和缩放至 1e6，未重按基因长度计算。
expr_clean <- expr_log2
expr_clean[is.na(expr_clean)] <- 0
expr_lin <- pmax(2^expr_clean - 1, 0)
col_sums <- colSums(expr_lin)
stopifnot(all(is.finite(col_sums)), all(col_sums > 0))
expr_tpm_scale <- sweep(expr_lin, 2L, 1e6 / col_sums, "*")
stopifnot(max(abs(colSums(expr_tpm_scale) - 1e6)) < 1e-6)
cat("RSEM input column sum range:", range(col_sums), "| TPM-scale sum range:",
    range(colSums(expr_tpm_scale)), "\n")

# 按 quantiseqr 本地 HGNC 表解析基因别名。先合并少量映射重复行，避免
# mapGenes 对整个宽矩阵逐列、逐基因组执行 aggregate(median)。
quantiseq_mapping <- function(symbols) {
  ref <- read.csv(system.file("extdata", "HGNC_genenames_20170418.txt.gz",
                              package = "quantiseqr", mustWork = TRUE), header = TRUE, sep = "\t")
  # 本地包的重注释代码读取 PreviousSymbols 与 Synonyms；有些包版本只包含前者。
  source_aliases <- Filter(Negate(is.null), list(ref$PreviousSymbols, ref$Synonyms))
  stopifnot(all(lengths(c(list(ref$ApprovedSymbol, ref$Status, ref$ApprovedName),
                          source_aliases)) == nrow(ref)))
  direct <- match(symbols, ref$ApprovedSymbol)
  approved <- which(!is.na(direct) & ref$Status[direct] == "Approved")
  withdrawn <- which(!is.na(direct) & !ref$Status[direct] %in% c("Approved", "EntryWithdrawn"))
  mapped <- alias <- rep(NA_character_, length(symbols))
  mapped[approved] <- symbols[approved]
  alias[withdrawn] <- gsub("symbolwithdrawn,see", "", ref$ApprovedName[direct[withdrawn]])
  unmatched <- which(is.na(direct))
  for (column in source_aliases) {
    aliases <- strsplit(as.character(column), ",", fixed = TRUE)
    keys <- unlist(aliases, use.names = FALSE)
    rows <- rep.int(seq_len(nrow(ref)), lengths(aliases))
    last <- !duplicated(keys, fromLast = TRUE)
    hits <- rows[last][match(symbols[unmatched], keys[last])]
    found <- !is.na(hits)
    alias[unmatched[found]] <- ref$ApprovedSymbol[hits[found]]
  }
  alias[alias %in% mapped[!is.na(mapped)]] <- NA_character_
  mapped[is.na(mapped)] <- alias[is.na(mapped)]
  mapped
}
merge_quantiseq_rows <- function(mat, mapped) {
  mat <- mat[!is.na(mapped), , drop = FALSE]
  mapped <- mapped[!is.na(mapped)]
  duplicate_names <- unique(mapped[duplicated(mapped)])
  for (symbol in duplicate_names) {
    ii <- which(mapped == symbol)
    mat[ii[1L], ] <- apply(mat[ii, , drop = FALSE], 2L, median)
  }
  # 保留每组代表行的原始符号，随后由 quanTIseq 完成一次正式重注释。
  mat[!duplicated(mapped), , drop = FALSE]
}
qmap <- quantiseq_mapping(rownames(expr_lin))
# 用原包在三个真实对齐样本上逐行验证映射和重复项聚合。
probe <- expr_lin[, unique(round(seq(1L, ncol(expr_lin), length.out = 3L))), drop = FALSE]
probe_collapsed <- as.matrix(quantiseqr:::mapGenes(merge_quantiseq_rows(probe, qmap)))
probe_pkg <- as.matrix(quantiseqr:::mapGenes(probe))
stopifnot(setequal(rownames(probe_collapsed), rownames(probe_pkg)),
          identical(probe_collapsed, probe_pkg[rownames(probe_collapsed), , drop = FALSE]))
cat("quanTIseq HGNC mapping verified across", nrow(probe_pkg), "genes and", ncol(probe), "samples\n")

# 对这两种输入重新运行 quanTIseq，参数和样本完全相同。
run_quantiseq <- function(mat) {
  immunedeconv::deconvolute(merge_quantiseq_rows(mat, qmap), method = "quantiseq",
                          tumor = TRUE, arrays = FALSE, scale_mrna = TRUE)
}
cat("=== Running quanTIseq on linear RSEM ===\n")
qs_rsem <- run_quantiseq(expr_lin)
cat("=== Running quanTIseq on column-rescaled RSEM (TPM-scale) ===\n")
qs_tpm_scale <- run_quantiseq(expr_tpm_scale)
rm(expr_tpm_scale); gc()

# MCP-counter 和 ESTIMATE 保持原有本地缓存路径；quanTIseq 结果以上方重算为准。
tme_cache <- file.path(bulk_dir, "tcga_tme_scores.rds")
if (file.exists(tme_cache)) {
  tme_res <- readRDS(tme_cache)
} else {
  tme_res <- list(mcp_counter = immunedeconv::deconvolute(expr_lin, method = "mcp_counter"),
                  estimate = immunedeconv::deconvolute(expr_lin, method = "estimate"))
  saveRDS(c(list(quantiseq = qs_rsem), tme_res), tme_cache)
}
rm(expr_lin); gc()

# 整理为特征 × 样本矩阵并按表达样本顺序索引。
to_sample_mat <- function(df, prefix) {
  stopifnot(!is.null(df), !anyDuplicated(df$cell_type))
  m <- as.matrix(as.data.frame(df)[, -1, drop = FALSE])
  rownames(m) <- paste0(prefix, "__", df$cell_type)
  stopifnot(!anyDuplicated(colnames(m)), setequal(colnames(m), sids_tcga))
  m[, sids_tcga, drop = FALSE]
}

m_qs  <- to_sample_mat(qs_rsem, "quanTIseq")
m_mcp <- to_sample_mat(tme_res$mcp_counter, "MCPcounter")
m_est <- to_sample_mat(tme_res$estimate, "ESTIMATE")
tme_mat <- rbind(m_qs, m_mcp, m_est)
cat("TME features extracted:", nrow(tme_mat), "across", ncol(tme_mat), "samples\n")

# 2. 核心基因与免疫细胞、免疫检查点、TMB、MSI 的 Spearman 相关性 (Overall / ER+ / ER-)
focus_genes <- c("PPARG", "CXCL8", "IL1B", "FFAR2", "HDAC11", "TNF", "IL12A")
checkpoints <- c("CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD80", "CD86", "IDO1", "FOXP3")
avail_cp <- intersect(checkpoints, rownames(expr_log2))

strata_list <- list(
  Overall = seq_len(ncol(expr_log2)),
  ER_pos  = which(normalize_er(xena_er) == "Positive"),
  ER_neg  = which(normalize_er(xena_er) == "Negative")
)

testable <- function(x) length(x) >= 20L && all(is.finite(x)) && sd(x) > 0
spearman <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman"))

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
      ok <- is.finite(gv) & is.finite(fv)
      if (testable(gv[ok]) && testable(fv[ok])) {
        ct <- spearman(gv[ok], fv[ok])
        cor_rows[[paste(st, g, feat)]] <- data.table(
          stratum = st, gene_symbol = g, feature_type = "TME_Cell_Score",
          feature = feat, n = sum(ok), rho = unname(ct$estimate), pval = ct$p.value
        )
      }
    }
    # 免疫检查点
    for (cp in avail_cp) {
      cv <- as.numeric(expr_log2[cp, ii])
      ok <- is.finite(gv) & is.finite(cv)
      if (testable(gv[ok]) && testable(cv[ok])) {
        ct <- spearman(gv[ok], cv[ok])
        cor_rows[[paste(st, g, "CP", cp)]] <- data.table(
          stratum = st, gene_symbol = g, feature_type = "Immune_Checkpoint",
          feature = cp, n = sum(ok), rho = unname(ct$estimate), pval = ct$p.value
        )
      }
    }
    # TMB 与 MSI
    for (gfeat in c("TMB_NONSYNONYMOUS", "MSI_SENSOR_SCORE")) {
      mv <- suppressWarnings(as.numeric(clin[[gfeat]][ii]))
      ok <- is.finite(gv) & is.finite(mv)
      if (testable(gv[ok]) && testable(mv[ok])) {
        ct <- spearman(gv[ok], mv[ok])
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
cat("Exact stratum sample counts in TCGA TME: Overall =", length(strata_list$Overall),
    "| ER_pos =", length(strata_list$ER_pos), "| ER_neg =", length(strata_list$ER_neg), "\n")

# 2b. 同一批样本和同一个 quanTIseq 检验族用于 RSEM 与 TPM-scale 的比较。
# 各 stratum × gene 中，仅以两种输入均可检验的、相同的 quanTIseq 特征作 BH 校正。
# 主表的 FDR 仍覆盖相应 stratum × gene × feature_type 的全部特征。
m_qs_tpm_scale <- to_sample_mat(qs_tpm_scale, "quanTIseq")
stopifnot(setequal(rownames(m_qs), rownames(m_qs_tpm_scale)))
m_qs_tpm_scale <- m_qs_tpm_scale[rownames(m_qs), , drop = FALSE]
stopifnot(all(is.finite(m_qs)), all(is.finite(m_qs_tpm_scale)))
score_agreement <- rbindlist(lapply(names(strata_list), function(st) {
  ii <- strata_list[[st]]
  rbindlist(lapply(rownames(m_qs), function(ft) {
    raw <- m_qs[ft, ii]
    scaled <- m_qs_tpm_scale[ft, ii]
    data.table(stratum = st, feature = ft, n = length(ii),
      score_spearman = if (testable(raw) && testable(scaled)) cor(raw, scaled, method = "spearman") else NA_real_,
      median_abs_fraction_delta = median(abs(scaled - raw)), max_abs_fraction_delta = max(abs(scaled - raw)))
  }))
}))
fwrite(score_agreement, file.path(res_dir, "module5_quantiseq_score_agreement.csv"))

sens_rows <- list()
for (st in names(strata_list)) {
  ii <- strata_list[[st]]
  for (g in intersect(focus_genes, rownames(expr_log2))) {
    gv <- as.numeric(expr_log2[g, ii])
    for (ft in rownames(m_qs)) {
      raw <- m_qs[ft, ii]
      scaled <- m_qs_tpm_scale[ft, ii]
      ok <- is.finite(gv) & is.finite(raw) & is.finite(scaled)
      if (!all(vapply(list(gv[ok], raw[ok], scaled[ok]), testable, logical(1)))) next
      a <- spearman(gv[ok], raw[ok])
      b <- spearman(gv[ok], scaled[ok])
      sens_rows[[length(sens_rows) + 1L]] <- data.table(
        stratum = st, gene_symbol = g, feature = ft, n = sum(ok),
        rho_rsem = unname(a$estimate), pval_rsem = a$p.value,
        rho_tpm_scale = unname(b$estimate), pval_tpm_scale = b$p.value)
    }
  }
}
tpm_cmp <- rbindlist(sens_rows)
stopifnot(nrow(tpm_cmp) > 0L, !anyDuplicated(tpm_cmp[, .(stratum, gene_symbol, feature)]))
tpm_cmp[, `:=`(FDR_rsem = p.adjust(pval_rsem, "BH"), FDR_tpm_scale = p.adjust(pval_tpm_scale, "BH"),
               n_tests_in_family = .N, delta_rho = rho_tpm_scale - rho_rsem), by = .(stratum, gene_symbol)]
primary_qs <- cor_dt[feature_type == "TME_Cell_Score" & feature %in% rownames(m_qs)]
checked <- primary_qs[tpm_cmp, on = .(stratum, gene_symbol, feature)]
stopifnot(all(checked$n == checked$i.n),
          isTRUE(all.equal(checked$pval, checked$pval_rsem)), isTRUE(all.equal(checked$rho, checked$rho_rsem)))
fwrite(tpm_cmp[order(stratum, gene_symbol, pval_tpm_scale)], file.path(res_dir, "module5_tme_tpm_sensitivity.csv"))
cat("quanTIseq RSEM vs TPM-scale Pearson r across paired gene-feature Spearman coefficients:",
    cor(tpm_cmp$rho_rsem, tpm_cmp$rho_tpm_scale), "| max absolute delta rho:", max(abs(tpm_cmp$delta_rho)), "\n")
print(tpm_cmp[(gene_symbol == "PPARG" & stratum == "ER_pos") | (gene_symbol == "CXCL8" & stratum == "ER_neg")][
  order(gene_symbol, pval_tpm_scale), .(stratum, gene_symbol, feature, n, rho_rsem, rho_tpm_scale, FDR_rsem, FDR_tpm_scale)])

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

# 4. HPA/TISCH2 本地待核清单。保留参考链接供手工核查。
ens_check <- c(PPARG = "ENSG00000132170", CXCL8 = "ENSG00000169429", IL1B = "ENSG00000125538", FFAR2 = "ENSG00000126262")
web_dt <- data.table(
  gene_symbol = names(ens_check), ensembl_id = unname(ens_check),
  hpa_api_status = "SKIPPED_OFFLINE",
  hpa_gene_url = sprintf("https://www.proteinatlas.org/%s-%s/pathology/breast+cancer", ens_check, names(ens_check)),
  tisch2_url = "http://tisch.comp-genomics.org/gallery/?cancer=BRCA",
  tisch2_recommended_cohorts = "BRCA_GSE176078; BRCA_EMTAB8107; BRCA_GSE114727_inDrop",
  hpa_pathology_note = "Check web manually")
fwrite(web_dt, file.path(res_dir, "module5_web_lookup_checklist.csv"))
