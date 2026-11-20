# 09_run_eqtlgen_all49_mr.R
# 模块 3A：对全部 61 个候选宿主基因中匹配到 eQTLGen 的 49 个基因跑全量 2SMR（严格保存完整分母与阴性结果）
# 严格阈值：p1 = 5e-8, r2 = 0.001, kb = 10000, F > 10，并做 cis (±1Mb 基因位置) 标注与 Steiger 方向检验

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

cache_dir <- "00_data/mr_cache"
res_dir   <- "02_results"

univ <- fread(file.path(res_dir, "module3_candidate_host_genes_universe.csv"))
targets <- univ[!is.na(eqtl_id) & nzchar(eqtl_id)]
cat("eQTLGen matched candidate genes to test (full denominator):", nrow(targets), "\n")

# 1. 提取全部 49 个基因的强工具变量 (p < 5e-8, r2 = 0.001, kb = 10000, F > 10)
safe_extract_eqtl_iv <- function(gid, max_try = 4) {
  rds_file <- file.path(cache_dir, paste0("iv_eqtl_", gid, ".rds"))
  if (file.exists(rds_file)) return(readRDS(rds_file))
  for (att in seq_len(max_try)) {
    res <- tryCatch({
      Sys.sleep(0.35)
      TwoSampleMR::extract_instruments(outcomes = gid, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
    }, error = function(e) {
      msg <- conditionMessage(e)
      message("  Attempt ", att, " error for ", gid, ": ", msg)
      Sys.sleep(4 * att)
      NULL
    })
    if (!is.null(res) || att == max_try) {
      if (is.null(res)) res <- data.frame()
      saveRDS(res, rds_file)
      return(res)
    }
  }
}

iv_list <- list()
status_dt <- copy(univ)
status_dt[, iv_count_raw := 0L]
status_dt[, iv_count_F10 := 0L]

for (i in seq_len(nrow(targets))) {
  gid <- targets$eqtl_id[i]
  sym <- targets$gene_symbol[i]
  if (i %% 10 == 1 || i == nrow(targets)) {
    cat(sprintf("[%d/%d] Extracting eQTLGen IVs for %s (%s)\n", i, nrow(targets), sym, gid))
  }
  iv <- safe_extract_eqtl_iv(gid)
  if (!is.null(iv) && nrow(iv) > 0) {
    status_dt[gene_symbol == sym, iv_count_raw := nrow(iv)]
    iv$F_stat <- (iv$beta.exposure / iv$se.exposure)^2
    iv <- iv[iv$F_stat > 10, ]
    if (nrow(iv) > 0) {
      iv$gene_symbol <- sym
      status_dt[gene_symbol == sym, iv_count_F10 := nrow(iv)]
      iv_list[[gid]] <- iv
    }
  }
}

exp_eqtl <- rbindlist(iv_list, fill = TRUE)
saveRDS(exp_eqtl, file.path(cache_dir, "eqtlgen_all49_exp_ivs.rds"))
fwrite(status_dt, file.path(res_dir, "module3_candidate_genes_iv_audit.csv"))
cat("Genes with >=1 strong IV (F>10):", length(iv_list), "/", nrow(targets), "| Total SNPs:", nrow(exp_eqtl), "\n")

if (nrow(exp_eqtl) == 0) stop("No eQTLGen IVs extracted!")

# 2. 批量提取 BCAC 结局数据 (Overall, ER+, ER-)
out_cache <- file.path(cache_dir, "outcome_bcac_eqtlgen_all49.rds")
outcomes_bcac <- c("ieu-a-1126", "ieu-a-1127", "ieu-a-1128")
unique_snps <- unique(exp_eqtl$SNP)
cat("Unique eQTLGen SNPs to query in BCAC:", length(unique_snps), "\n")

if (file.exists(out_cache)) {
  out_all <- readRDS(out_cache)
} else {
  snp_chunks <- split(unique_snps, ceiling(seq_along(unique_snps) / 250))
  out_chunks <- list()
  for (o_id in outcomes_bcac) {
    cat("Fetching BCAC outcome:", o_id, "\n")
    for (ci in seq_along(snp_chunks)) {
      Sys.sleep(0.45)
      r <- tryCatch(
        TwoSampleMR::extract_outcome_data(snps = snp_chunks[[ci]], outcomes = o_id, proxies = TRUE, rsq = 0.8),
        error = function(e) { Sys.sleep(5); tryCatch(TwoSampleMR::extract_outcome_data(snps = snp_chunks[[ci]], outcomes = o_id, proxies = TRUE, rsq = 0.8), error = function(e2) NULL) }
      )
      if (!is.null(r) && nrow(r) > 0) out_chunks[[paste(o_id, ci, sep = "_")]] <- r
    }
  }
  out_all <- rbindlist(out_chunks, fill = TRUE)
  saveRDS(out_all, out_cache)
}

# 3. Harmonise & 跑全部候选基因的 MR + 异质性 + 多效性 + Steiger 方向检验
dat_harm <- TwoSampleMR::harmonise_data(as.data.frame(exp_eqtl), as.data.frame(out_all), action = 2)
saveRDS(dat_harm, file.path(cache_dir, "eqtlgen_bcac_harmonised.rds"))

mr_res  <- as.data.table(TwoSampleMR::mr(dat_harm, method_list = c("mr_ivw", "mr_wald_ratio", "mr_egger_regression", "mr_weighted_median")))
mr_het  <- as.data.table(TwoSampleMR::mr_heterogeneity(dat_harm))
mr_plei <- as.data.table(TwoSampleMR::mr_pleiotropy_test(dat_harm))
steiger <- tryCatch(as.data.table(TwoSampleMR::directionality_test(dat_harm)), error = function(e) data.table())

mr_res[, OR := exp(b)]
mr_res[, OR_LCI95 := exp(b - 1.96 * se)]
mr_res[, OR_UCI95 := exp(b + 1.96 * se)]

het_ivw <- mr_het[method == "Inverse variance weighted", .(id.exposure, id.outcome, Q, Q_pval)]
plei_dt <- mr_plei[, .(id.exposure, id.outcome, egger_intercept, pleiotropy_pval = pval)]
dir_dt  <- mr_res[, .(all_dir_consistent = (all(b > 0) | all(b < 0))), by = .(id.exposure, id.outcome)]

# 合并上游 gutMGene 元信息（保留完整分母记录）
mr_full <- merge(mr_res, targets[, .(id.exposure = eqtl_id, gene_symbol, upstream_microbes, upstream_metabolites, alterations)],
                 by = "id.exposure", all.x = TRUE)
mr_full <- merge(mr_full, het_ivw, by = c("id.exposure", "id.outcome"), all.x = TRUE)
mr_full <- merge(mr_full, plei_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)
mr_full <- merge(mr_full, dir_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)
if (nrow(steiger) > 0) {
  mr_full <- merge(mr_full, steiger[, .(id.exposure, id.outcome, correct_causal_direction, steiger_pval)],
                   by = c("id.exposure", "id.outcome"), all.x = TRUE)
}
mr_full[, outcome_label := fifelse(id.outcome == "ieu-a-1126", "BRCA_Overall",
                            fifelse(id.outcome == "ieu-a-1127", "BRCA_ER_pos", "BRCA_ER_neg"))]

prim_all <- mr_full[method %in% c("Inverse variance weighted", "Wald ratio")]
# 按结局在完整分母上计算 FDR
prim_all[, FDR := p.adjust(pval, method = "BH"), by = id.outcome]

fwrite(mr_full, file.path(res_dir, "module3_eqtlgen_all49_mr_all_methods.csv"))
fwrite(prim_all[order(outcome_label, pval)], file.path(res_dir, "module3_eqtlgen_all49_mr_primary_full_denominator.csv"))

sig_eqtl <- prim_all[pval < 0.05][order(outcome_label, pval)]
fwrite(sig_eqtl, file.path(res_dir, "module3_eqtlgen_mr_positive_p005.csv"))

cat("\n=== 模块 3A eQTLGen -> BCAC 全量 49 基因 MR 结果汇总 ===\n")
cat("Tested pairs (primary):", nrow(prim_all), "| Unique genes tested:", length(unique(prim_all$gene_symbol)), "\n")
cat("Nominal significant pairs (P < 0.05):", nrow(sig_eqtl), "| Unique positive genes:", length(unique(sig_eqtl$gene_symbol)), "\n")
cat("FDR < 0.05 pairs:", sum(prim_all$FDR < 0.05, na.rm = TRUE), "| FDR < 0.10 pairs:", sum(prim_all$FDR < 0.10, na.rm = TRUE), "\n\n")
if (nrow(sig_eqtl) > 0) {
  print(sig_eqtl[, .(outcome_label, gene_symbol, method, nsnp, OR = round(OR, 3), pval = signif(pval, 3),
                     FDR = signif(FDR, 3), all_dir_consistent, Q_pval = signif(Q_pval, 2),
                     pleiotropy_pval = signif(pleiotropy_pval, 2), correct_causal_direction,
                     upstream_metabolites, alterations)])
}
