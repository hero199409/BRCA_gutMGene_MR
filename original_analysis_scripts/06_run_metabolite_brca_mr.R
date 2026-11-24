# 06_run_metabolite_brca_mr.R
# 血浆代谢物 (与 gutMGene 匹配的 47 个 met-* GWAS) -> BCAC 乳腺癌 (Overall / ER+ / ER-) 两样本 MR

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
})

cache_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/mr_cache")
res_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")
ao <- readRDS(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/opengwas_available_outcomes.rds"))

syns <- c("hydroxybutyrate","bohbutyrate","indolepropionate","indoleacetate","kynurenate",
          "lithocholate","deoxycholate","nicotinate","propionate","butyrate","succinate",
          "serotonin","pyruvate","acetate","cresol","phenylacetyl","tmao","trimethylamine",
          "equol","urolithin","caffeate","mannose","valerate","lactate","linoleic")

met_list <- ao[grepl("^met-", id) & grepl(paste(syns, collapse = "|"), trait, ignore.case = TRUE),
               .(id, trait, author, year, sample_size)]
# 排除比率表型仅保留浓度表型？保留全部 47 个并在结果中标注
cat("Target gutMGene-aligned metabolite GWAS count:", nrow(met_list), "\n")
fwrite(met_list, file.path(res_dir, "target_47_metabolite_gwas.csv"))

safe_extract_met_iv <- function(gid, author_name, max_try = 3) {
  rds_file <- file.path(cache_dir, paste0("iv_met_", gid, ".rds"))
  if (file.exists(rds_file)) return(readRDS(rds_file))
  
  # Borges (n~115k) 优先用 5e-8，不足 3 个则回退 1e-5；Kettunen/Shin 用 1e-5
  p_thresh <- if (author_name == "Borges CM") 5e-8 else 1e-5
  
  for (att in seq_len(max_try)) {
    res <- tryCatch({
      Sys.sleep(0.35)
      iv <- TwoSampleMR::extract_instruments(gid, p1 = p_thresh, clump = TRUE, r2 = 0.001, kb = 10000)
      if ((is.null(iv) || nrow(iv) < 3) && p_thresh == 5e-8) {
        Sys.sleep(0.35)
        iv <- TwoSampleMR::extract_instruments(gid, p1 = 1e-5, clump = TRUE, r2 = 0.001, kb = 10000)
      }
      iv
    }, error = function(e) {
      message("  Attempt ", att, " error for ", gid, ": ", conditionMessage(e))
      Sys.sleep(3 * att)
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
for (i in seq_len(nrow(met_list))) {
  gid <- met_list$id[i]
  if (i %% 10 == 1 || i == nrow(met_list)) {
    cat(sprintf("[%d/%d] Extracting IVs for %s (%s)\n", i, nrow(met_list), gid, met_list$trait[i]))
  }
  iv <- safe_extract_met_iv(gid, met_list$author[i])
  if (!is.null(iv) && nrow(iv) > 0) {
    iv$F_stat <- (iv$beta.exposure / iv$se.exposure)^2
    iv <- iv[iv$F_stat > 10, ]
    if (nrow(iv) > 0) iv_list[[gid]] <- iv
  }
}

exp_met <- rbindlist(iv_list, fill = TRUE)
cat("Metabolites with strong IVs:", length(iv_list), "| Total SNPs:", nrow(exp_met), "\n")
saveRDS(exp_met, file.path(cache_dir, "met_aligned_exp_ivs.rds"))

out_cache <- file.path(cache_dir, "outcome_bcac_met_aligned.rds")
outcomes_bcac <- c("ieu-a-1126", "ieu-a-1127", "ieu-a-1128")
unique_snps <- unique(exp_met$SNP)
cat("Unique metabolite SNPs to query in BCAC:", length(unique_snps), "\n")

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
        error = function(e) { Sys.sleep(4); tryCatch(TwoSampleMR::extract_outcome_data(snps = snp_chunks[[ci]], outcomes = o_id, proxies = TRUE, rsq = 0.8), error = function(e2) NULL) }
      )
      if (!is.null(r) && nrow(r) > 0) out_chunks[[paste(o_id, ci, sep = "_")]] <- r
    }
  }
  out_all <- rbindlist(out_chunks, fill = TRUE)
  saveRDS(out_all, out_cache)
}

dat_harm <- TwoSampleMR::harmonise_data(as.data.frame(exp_met), as.data.frame(out_all), action = 2)
saveRDS(dat_harm, file.path(cache_dir, "met_bcac_harmonised.rds"))

mr_fast <- as.data.table(TwoSampleMR::mr(dat_harm, method_list = c("mr_ivw", "mr_wald_ratio", "mr_egger_regression", "mr_weighted_median")))
mr_het  <- as.data.table(TwoSampleMR::mr_heterogeneity(dat_harm))
mr_plei <- as.data.table(TwoSampleMR::mr_pleiotropy_test(dat_harm))

mr_fast[, OR := exp(b)]
mr_fast[, OR_LCI95 := exp(b - 1.96 * se)]
mr_fast[, OR_UCI95 := exp(b + 1.96 * se)]

het_ivw <- mr_het[method == "Inverse variance weighted", .(id.exposure, id.outcome, Q_pval)]
plei_dt <- mr_plei[, .(id.exposure, id.outcome, pleiotropy_pval = pval)]
dir_dt  <- mr_fast[, .(all_dir_consistent = (all(b > 0) | all(b < 0))), by = .(id.exposure, id.outcome)]

mr_fast <- merge(mr_fast, met_list[, .(id.exposure = id, trait_name = trait, author, sample_size)], by = "id.exposure", all.x = TRUE)
mr_fast <- merge(mr_fast, het_ivw, by = c("id.exposure", "id.outcome"), all.x = TRUE)
mr_fast <- merge(mr_fast, plei_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)
mr_fast <- merge(mr_fast, dir_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)
mr_fast[, outcome_label := fifelse(id.outcome == "ieu-a-1126", "BRCA_Overall",
                            fifelse(id.outcome == "ieu-a-1127", "BRCA_ER_pos", "BRCA_ER_neg"))]

prim <- mr_fast[method %in% c("Inverse variance weighted", "Wald ratio")]
prim[, FDR := p.adjust(pval, method = "BH"), by = id.outcome]

# 仅对 P < 0.05 候选跑 MR-PRESSO (最多十几个，很快)
sig_pairs <- unique(prim[pval < 0.05, .(id.exposure, id.outcome)])
cat("Significant metabolite-BRCA pairs (P < 0.05):", nrow(sig_pairs), "\n")
presso_dt <- data.table()
if (nrow(sig_pairs) > 0) {
  dat_keep <- as.data.table(dat_harm)[mr_keep == TRUE]
  for (k in seq_len(nrow(sig_pairs))) {
    sub_d <- as.data.frame(dat_keep[id.exposure == sig_pairs$id.exposure[k] & id.outcome == sig_pairs$id.outcome[k]])
    if (nrow(sub_d) >= 4) {
      pr <- tryCatch(
        MRPRESSO::mr_presso(BetaOutcome="beta.outcome", BetaExposure="beta.exposure",
                            SdOutcome="se.outcome", SdExposure="se.exposure",
                            OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data=sub_d,
                            NbDistribution=500, SignifThreshold=0.05),
        error = function(e) NULL
      )
      if (!is.null(pr)) {
        presso_dt <- rbind(presso_dt, data.table(
          id.exposure = sig_pairs$id.exposure[k],
          id.outcome  = sig_pairs$id.outcome[k],
          presso_global_p = as.character(pr$`MR-PRESSO results`$`Global Test`$Pvalue)
        ))
      }
    }
  }
}

if (nrow(presso_dt) > 0) {
  prim <- merge(prim, presso_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)
}

fwrite(mr_fast, file.path(res_dir, "metabolite_bcac_mr_all_methods.csv"))
fwrite(prim[order(pval)], file.path(res_dir, "metabolite_bcac_mr_primary.csv"))
sig_met <- prim[pval < 0.05][order(outcome_label, pval)]
fwrite(sig_met, file.path(res_dir, "metabolite_bcac_mr_significant_p005.csv"))

cat("\n=== 血浆代谢物 -> BCAC 名义显著结果 (P < 0.05) ===\n")
if (nrow(sig_met) > 0) {
  print(sig_met[, .(outcome_label, id.exposure, trait_name, author, nsnp, OR = round(OR, 3),
                    pval = signif(pval, 3), FDR = signif(FDR, 3), all_dir_consistent,
                    Q_pval = signif(Q_pval, 2), pleiotropy_pval = signif(pleiotropy_pval, 2))])
} else {
  cat("No metabolite reached P < 0.05.\n")
}
