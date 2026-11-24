# 05_mibiogen_finngen_and_reverse_mr.R
# 1) 用 FinnGen 独立队列 (finn-b-C3_BREAST, finn-b-C3_BREAST_EXALLC, HER+/HER-) 验证 7 个 MiBioGen 阳性菌群并做 Meta 分析
# 2) 做反向 MR (BCAC 总乳腺癌/ER+/ER- -> 7 个阳性菌群)，排除因果倒置

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(meta)
})

cache_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/mr_cache")
res_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

sig_bcac <- fread(file.path(res_dir, "mibiogen_bcac_mr_significant_p005.csv"))
pos_ids  <- unique(sig_bcac$id.exposure)
pos_meta <- unique(sig_bcac[, .(id.exposure, gwas_rank, gwas_name)])
cat("Positive MiBioGen taxa to validate:", length(pos_ids), "\n")
print(pos_meta)

# 1. 读取已缓存的 7 个阳性菌群 IVs（0 API 消耗）
exp_all <- readRDS(file.path(cache_dir, "mibiogen_aligned_exp_ivs.rds"))
exp_pos <- as.data.frame(exp_all[exp_all$id.exposure %in% pos_ids, ])
unique_snps <- unique(exp_pos$SNP)
cat("Positive taxa SNPs:", length(unique_snps), "\n")

# 2. 提取 FinnGen 结局验证数据（带缓存）
finn_cache <- file.path(cache_dir, "outcome_finngen_pos_mibiogen.rds")
finngen_ids <- c("finn-b-C3_BREAST_EXALLC", "finn-b-C3_BREAST", "finn-b-C3_BREAST_HERPLUS_EXALLC", "finn-b-C3_BREAST_HERNEG_EXALLC")

if (file.exists(finn_cache)) {
  out_finn <- readRDS(finn_cache)
} else {
  out_chunks <- list()
  for (fid in finngen_ids) {
    cat("Fetching FinnGen outcome:", fid, "\n")
    Sys.sleep(0.5)
    res <- tryCatch(
      TwoSampleMR::extract_outcome_data(snps = unique_snps, outcomes = fid, proxies = TRUE, rsq = 0.8),
      error = function(e) { message("  Error: ", conditionMessage(e)); NULL }
    )
    if (!is.null(res) && nrow(res) > 0) out_chunks[[fid]] <- res
  }
  out_finn <- as.data.frame(rbindlist(out_chunks, fill = TRUE))
  saveRDS(out_finn, finn_cache)
}

dat_finn <- TwoSampleMR::harmonise_data(exp_pos, out_finn, action = 2)
mr_finn  <- as.data.table(TwoSampleMR::mr(dat_finn, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median", "mr_egger_regression")))
mr_finn[, OR := exp(b)]
mr_finn[, OR_LCI95 := exp(b - 1.96 * se)]
mr_finn[, OR_UCI95 := exp(b + 1.96 * se)]
mr_finn <- merge(mr_finn, pos_meta, by = "id.exposure", all.x = TRUE)
fwrite(mr_finn, file.path(res_dir, "mibiogen_finngen_validation_mr.csv"))

finn_prim <- mr_finn[method %in% c("Inverse variance weighted", "Wald ratio")][order(pval)]
cat("\n=== FinnGen 独立验证结果 (按 P 值排序) ===\n")
print(finn_prim[, .(id.outcome, gwas_rank, gwas_name, nsnp, OR = round(OR, 3), pval = signif(pval, 3))])

# 3. 合并 BCAC_Overall + FinnGen (finn-b-C3_BREAST_EXALLC) 做 MR Meta 分析
bcac_prim <- fread(file.path(res_dir, "mibiogen_bcac_mr_primary.csv"))
bcac_ov   <- bcac_prim[id.exposure %in% pos_ids & id.outcome == "ieu-a-1126"]
finn_ov   <- finn_prim[id.outcome == "finn-b-C3_BREAST_EXALLC"]

meta_rows <- list()
for (tid in pos_ids) {
  b1 <- bcac_ov[id.exposure == tid]
  b2 <- finn_ov[id.exposure == tid]
  if (nrow(b1) == 1 && nrow(b2) == 1) {
    m <- meta::metagen(
      TE = c(b1$b, b2$b),
      seTE = c(b1$se, b2$se),
      studlab = c("BCAC", "FinnGen"),
      sm = "OR"
    )
    use_random <- (m$pval.Q < 0.05)
    te_use <- ifelse(use_random, m$TE.random, m$TE.fixed)
    se_use <- ifelse(use_random, m$seTE.random, m$seTE.fixed)
    p_use  <- ifelse(use_random, m$pval.random, m$pval.fixed)
    meta_rows[[tid]] <- data.table(
      id.exposure = tid,
      gwas_name = b1$gwas_name,
      BCAC_OR = round(exp(b1$b), 3), BCAC_P = signif(b1$pval, 3),
      FinnGen_OR = round(exp(b2$b), 3), FinnGen_P = signif(b2$pval, 3),
      same_direction = (sign(b1$b) == sign(b2$b)),
      Meta_model = ifelse(use_random, "Random", "Fixed"),
      Meta_OR = round(exp(te_use), 3),
      Meta_LCI95 = round(exp(te_use - 1.96 * se_use), 3),
      Meta_UCI95 = round(exp(te_use + 1.96 * se_use), 3),
      Meta_P = signif(p_use, 3),
      I2 = round(m$I2 * 100, 1),
      Het_P = signif(m$pval.Q, 3)
    )
  }
}
meta_dt <- rbindlist(meta_rows)[order(Meta_P)]
fwrite(meta_dt, file.path(res_dir, "mibiogen_bcac_finngen_meta.csv"))
cat("\n=== BCAC (Overall) + FinnGen Meta 分析汇总 ===\n")
print(meta_dt)

# 4. 反向 MR：BCAC (Overall / ER+ / ER-) -> 7 个阳性菌群
cat("\n=== Running Reverse MR: BCAC -> 7 Positive Taxa ===\n")
rev_exp_cache <- file.path(cache_dir, "reverse_bcac_ivs.rds")
bcac_ids <- c("ieu-a-1126", "ieu-a-1127", "ieu-a-1128")

if (file.exists(rev_exp_cache)) {
  rev_exp <- readRDS(rev_exp_cache)
} else {
  rev_list <- list()
  for (bid in bcac_ids) {
    cat("Extracting BRCA IVs (p1=5e-8):", bid, "\n")
    Sys.sleep(0.5)
    iv <- TwoSampleMR::extract_instruments(outcomes = bid, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
    if (!is.null(iv) && nrow(iv) > 0) {
      iv$F_stat <- (iv$beta.exposure / iv$se.exposure)^2
      rev_list[[bid]] <- iv[iv$F_stat > 10, ]
    }
  }
  rev_exp <- as.data.frame(rbindlist(rev_list, fill = TRUE))
  saveRDS(rev_exp, rev_exp_cache)
}
cat("Reverse MR BRCA IV count:", nrow(rev_exp), "unique SNPs:", length(unique(rev_exp$SNP)), "\n")

rev_out_cache <- file.path(cache_dir, "reverse_mibiogen_outcome.rds")
if (file.exists(rev_out_cache)) {
  rev_out <- readRDS(rev_out_cache)
} else {
  rev_snps <- unique(rev_exp$SNP)
  snp_chunks <- split(rev_snps, ceiling(seq_along(rev_snps) / 250))
  ro_list <- list()
  for (tid in pos_ids) {
    cat("Querying reverse outcome taxon:", tid, "\n")
    for (ci in seq_along(snp_chunks)) {
      Sys.sleep(0.4)
      r <- tryCatch(
        TwoSampleMR::extract_outcome_data(snps = snp_chunks[[ci]], outcomes = tid, proxies = TRUE, rsq = 0.8),
        error = function(e) NULL
      )
      if (!is.null(r) && nrow(r) > 0) ro_list[[paste(tid, ci, sep = "_")]] <- r
    }
  }
  rev_out <- as.data.frame(rbindlist(ro_list, fill = TRUE))
  saveRDS(rev_out, rev_out_cache)
}

dat_rev <- TwoSampleMR::harmonise_data(rev_exp, rev_out, action = 2)
mr_rev  <- as.data.table(TwoSampleMR::mr(dat_rev, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
mr_rev <- merge(mr_rev, pos_meta[, .(id.outcome = id.exposure, gwas_name)], by = "id.outcome", all.x = TRUE)
mr_rev[, exp_label := fifelse(id.exposure == "ieu-a-1126", "BRCA_Overall",
                       fifelse(id.exposure == "ieu-a-1127", "BRCA_ER_pos", "BRCA_ER_neg"))]
fwrite(mr_rev, file.path(res_dir, "mibiogen_reverse_mr.csv"))

rev_ivw <- mr_rev[method == "Inverse variance weighted"][order(pval)]
cat("\n=== 反向 MR (IVW) 结果：检验是否存在因果倒置 (期望 P > 0.05) ===\n")
print(rev_ivw[, .(exp_label, target_taxon = gwas_name, nsnp, b = round(b, 4), pval = signif(pval, 3))])
