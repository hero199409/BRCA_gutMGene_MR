# 03_run_mibiogen_brca_mr.R
# 与 gutMGene 对齐的 MiBioGen 肠道菌群 -> BCAC 乳腺癌（总体 / ER+ / ER-）两样本 MR。
# 对主方法 P < 0.10 的配对计算 Mode 法与 MR-PRESSO。

project_root <- Sys.getenv("MR_PROJECT_ROOT", unset = Sys.getenv("BRCA_PROJECT_ROOT", unset = "."))
source(file.path(project_root, "01_scripts", "00_mr_extract_helpers.R"))
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

cache_dir <- Sys.getenv("MR_CACHE_DIR", unset = file.path(project_root, "00_data", "mr_cache"))
input_dir <- file.path(project_root, "02_results")
res_dir <- Sys.getenv("MR_OUTPUT_DIR", unset = input_dir)

aligned <- fread(file.path(input_dir, "aligned_microbiome_gwas.csv"))
mib <- unique(aligned[cohort == "MiBioGen_2021"])
stopifnot(nrow(mib) > 0L, !anyDuplicated(mib$id))
bcac_ids <- names(mr_bcac_labels)

# 优先读协调后的 RDS。若缺失，先找原始暴露/结局 RDS，再按需提取。
harm_cache <- file.path(cache_dir, "mibiogen_bcac_harmonised.rds")
dat_harm <- as.data.frame(mr_cache_rds(harm_cache, function() {
  exp_cache <- file.path(cache_dir, "mibiogen_aligned_exp_ivs.rds")
  exp_all <- mr_cache_rds(exp_cache, function() {
    as.data.frame(rbindlist(lapply(mib$id, mr_fetch_iv, p_threshold = 1e-5), fill = TRUE))
  })
  exp_all <- as.data.frame(exp_all)
  exp_all <- mr_filter_iv(exp_all[exp_all$id.exposure %in% mib$id, , drop = FALSE], 1e-5)
  if (!nrow(exp_all)) stop("No eligible MiBioGen instruments.")
  out_cache <- file.path(cache_dir, "outcome_bcac_mibiogen_aligned.rds")
  out_all <- mr_cache_rds(out_cache, function() mr_fetch_outcomes(exp_all$SNP, bcac_ids, proxies = TRUE))
  out_all <- mr_select_outcomes(out_all, bcac_ids, exp_all$SNP, basename(out_cache))
  mr_harmonise(exp_all, out_all)
}))
dat_harm <- dat_harm[dat_harm$id.exposure %in% mib$id & dat_harm$id.outcome %in% bcac_ids, , drop = FALSE]
if (!all(bcac_ids %in% dat_harm$id.outcome)) stop("Incomplete BCAC harmonised cache.")
print(as.data.table(dat_harm)[, .(rows = .N, kept = sum(mr_keep %in% TRUE),
                                eligible_taxa = uniqueN(id.exposure[mr_keep %in% TRUE])), by = id.outcome])
cat("Aligned MiBioGen taxa:", nrow(mib), "\n")
set.seed(20261122) # 保证 bootstrap 和 MR-PRESSO 的 Monte Carlo 结果可复跑


# 2. 先跑全部菌群的主方法 + Egger + Weighted Median（几秒钟完成）
cat("Running fast MR (IVW / Wald / Egger / Weighted Median) on all pairs...\n")
mr_fast <- TwoSampleMR::mr(
  dat_harm,
  method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_wald_ratio")
)
mr_het  <- TwoSampleMR::mr_heterogeneity(dat_harm)
mr_plei <- TwoSampleMR::mr_pleiotropy_test(dat_harm)

mr_fast_dt <- as.data.table(mr_fast)
primary_dt <- mr_fast_dt[method %in% c("Inverse variance weighted", "Wald ratio")]
primary_dt[, FDR := p.adjust(pval, method = "BH"), by = id.outcome]

# 选出 P_primary < 0.10 的候选配对，补全 Simple/Weighted Mode 和 MR-PRESSO (NbDistribution=1000)
cand_pairs <- unique(primary_dt[pval < 0.10, .(id.exposure, id.outcome)])
cat("Candidate pairs with P_primary < 0.10:", nrow(cand_pairs), "\n")

mr_mode_dt <- data.table()
presso_dt <- data.table()

if (nrow(cand_pairs) > 0) {
  dat_cand <- merge(as.data.table(dat_harm), cand_pairs, by = c("id.exposure", "id.outcome"))
  mr_mode <- TwoSampleMR::mr(
    as.data.frame(dat_cand),
    method_list = c("mr_simple_mode", "mr_weighted_mode")
  )
  mr_mode_dt <- as.data.table(mr_mode)

  dat_keep <- dat_cand[mr_keep == TRUE]
  for (k in seq_len(nrow(cand_pairs))) {
    sub_d <- as.data.frame(dat_keep[id.exposure == cand_pairs$id.exposure[k] &
                                      id.outcome == cand_pairs$id.outcome[k]])
    if (nrow(sub_d) >= 4) {
      pr <- tryCatch({
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome", SdExposure = "se.exposure",
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sub_d,
          NbDistribution = 1000, SignifThreshold = 0.05
        )
      }, error = function(e) NULL)
      if (!is.null(pr)) {
        gp <- pr$`MR-PRESSO results`$`Global Test`$Pvalue
        presso_dt <- rbind(presso_dt, data.table(
          id.exposure = cand_pairs$id.exposure[k],
          id.outcome  = cand_pairs$id.outcome[k],
          presso_global_p = as.character(gp)
        ))
      }
    }
  }
}

# 3. 合并所有方法与敏感性检验
mr_dt <- rbindlist(list(mr_fast_dt, mr_mode_dt), fill = TRUE)
mr_dt[, OR := exp(b)]
mr_dt[, OR_LCI95 := exp(b - 1.96 * se)]
mr_dt[, OR_UCI95 := exp(b + 1.96 * se)]

het_ivw <- as.data.table(mr_het)[method == "Inverse variance weighted", .(id.exposure, id.outcome, Q, Q_pval)]
plei_dt <- as.data.table(mr_plei)[, .(id.exposure, id.outcome, egger_intercept, pleiotropy_pval = pval)]

dir_dt <- mr_dt[, .(
  n_methods = .N,
  all_dir_consistent = (all(b > 0) | all(b < 0))
), by = .(id.exposure, id.outcome)]

mr_dt <- merge(mr_dt, mib[, .(id.exposure = id, gwas_rank, gwas_name, has_direct_or_triplet_gene, has_metabolite)],
               by = "id.exposure", all.x = TRUE)
mr_dt <- merge(mr_dt, het_ivw, by = c("id.exposure", "id.outcome"), all.x = TRUE)
mr_dt <- merge(mr_dt, plei_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)
mr_dt <- merge(mr_dt, dir_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)
if (nrow(presso_dt) > 0) {
  mr_dt <- merge(mr_dt, presso_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)
} else {
  mr_dt[, presso_global_p := NA_character_]
}

mr_dt[, outcome_label := unname(mr_bcac_labels[id.outcome])]

primary_out <- mr_dt[method %in% mr_primary_methods]
primary_out[, FDR := p.adjust(pval, method = "BH"), by = id.outcome]

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(mr_dt, file.path(res_dir, "mibiogen_bcac_mr_all_methods.csv"))
fwrite(primary_out[order(pval)], file.path(res_dir, "mibiogen_bcac_mr_primary.csv"))

sig_dt <- primary_out[pval < 0.05][order(outcome_label, pval)]
fwrite(sig_dt, file.path(res_dir, "mibiogen_bcac_mr_significant_p005.csv"))

cat("\n=== MiBioGen -> BCAC 名义显著结果 (P_IVW/Wald < 0.05) 汇总 ===\n")
cat("Total primary MR pairs tested:", nrow(primary_out), "\n")
cat("Nominal significant (P < 0.05):", nrow(sig_dt), "\n")
cat("FDR < 0.05:", sum(primary_out$FDR < 0.05, na.rm = TRUE), " | FDR < 0.10:", sum(primary_out$FDR < 0.10, na.rm = TRUE), "\n\n")
if (nrow(sig_dt) > 0) {
  print(sig_dt[, .(outcome_label, gwas_rank, gwas_name, nsnp, OR = round(OR, 3), pval = signif(pval, 3), FDR = signif(FDR, 3),
                   all_dir_consistent, Q_pval = signif(Q_pval, 2), pleiotropy_pval = signif(pleiotropy_pval, 2),
                   presso_global_p, has_direct_or_triplet_gene)])
}
