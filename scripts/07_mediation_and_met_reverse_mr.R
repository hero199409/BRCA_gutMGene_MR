# 07_mediation_and_met_reverse_mr.R (快速修复版：避免 met-a-* 上百个缺失 SNP 的代理检索卡死)
source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
})

cache_dir <- "00_data/mr_cache"
res_dir   <- "02_results"
sig_met_ids <- c("met-d-Acetate", "met-c-893", "met-a-462", "met-a-385", "met-a-446")

# --- Part 2: 代谢物反向 MR (BCAC -> 5 个阳性代谢物, proxies=FALSE 防 Shin HapMap2 卡死) ---
rev_exp <- readRDS(file.path(cache_dir, "reverse_bcac_ivs.rds"))
rev_snps <- unique(rev_exp$SNP)

rev_met_cache <- file.path(cache_dir, "reverse_met_outcome.rds")
if (file.exists(rev_met_cache)) {
  rev_met_out <- readRDS(rev_met_cache)
} else {
  snp_chunks <- split(rev_snps, ceiling(seq_along(rev_snps) / 250))
  rm_list <- list()
  for (mid in sig_met_ids) {
    cat("Querying reverse outcome metabolite (no slow proxy):", mid, "\n")
    use_proxy <- grepl("^met-d-|^met-c-", mid)
    for (ci in seq_along(snp_chunks)) {
      Sys.sleep(0.35)
      r <- tryCatch(
        TwoSampleMR::extract_outcome_data(snps = snp_chunks[[ci]], outcomes = mid, proxies = use_proxy, rsq = 0.8),
        error = function(e) NULL
      )
      if (!is.null(r) && nrow(r) > 0) rm_list[[paste(mid, ci, sep = "_")]] <- r
    }
  }
  rev_met_out <- as.data.frame(rbindlist(rm_list, fill = TRUE))
  saveRDS(rev_met_out, rev_met_cache)
}

dat_rev_met <- TwoSampleMR::harmonise_data(rev_exp, rev_met_out, action = 2)
mr_rev_met  <- as.data.table(TwoSampleMR::mr(dat_rev_met, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression")))
fwrite(mr_rev_met, file.path(res_dir, "metabolite_reverse_mr.csv"))

cat("\n=== 代谢物反向 MR (BCAC -> 5 个阳性代谢物, 期望 P > 0.05) ===\n")
print(mr_rev_met[method == "Inverse variance weighted"][order(pval), .(id.exposure, id.outcome, nsnp, b = round(b, 4), pval = signif(pval, 3))])

# --- Part 3: 5 个阳性代谢物 -> FinnGen 乳腺癌独立验证 ---
exp_met <- readRDS(file.path(cache_dir, "met_aligned_exp_ivs.rds"))
exp_met_sig <- as.data.frame(exp_met[exp_met$id.exposure %in% sig_met_ids, ])
met_snps <- unique(exp_met_sig$SNP)

finn_met_cache <- file.path(cache_dir, "outcome_finngen_sig_met.rds")
if (file.exists(finn_met_cache)) {
  out_finn_met <- readRDS(finn_met_cache)
} else {
  cat("Fetching FinnGen outcome for 5 significant metabolites...\n")
  Sys.sleep(0.4)
  out_finn_met <- TwoSampleMR::extract_outcome_data(
    snps = met_snps,
    outcomes = c("finn-b-C3_BREAST_EXALLC", "finn-b-C3_BREAST"),
    proxies = TRUE, rsq = 0.8
  )
  saveRDS(out_finn_met, finn_met_cache)
}

dat_finn_met <- TwoSampleMR::harmonise_data(exp_met_sig, out_finn_met, action = 2)
mr_finn_met  <- as.data.table(TwoSampleMR::mr(dat_finn_met, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression")))
mr_finn_met[, OR := exp(b)]
fwrite(mr_finn_met, file.path(res_dir, "metabolite_finngen_validation_mr.csv"))

cat("\n=== 5 个阳性代谢物 -> FinnGen 验证结果 ===\n")
print(mr_finn_met[method == "Inverse variance weighted"][order(pval), .(id.outcome, id.exposure, exposure, nsnp, OR = round(OR, 3), pval = signif(pval, 3))])
