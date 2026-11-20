# 19_module3_dmp_fr02_and_pqtl.R
# 模块 3 加餐：
# 1) DMP_2022 (84个) + FR02_2022 (36个) 独立菌群队列全量 MR 复现与跨队列 Meta
# 2) OpenGWAS 血浆 IL-8 / CXCL8 pQTL 蛋白层两样本 MR + 共定位 (及 PPARG pQTL 检索审计)

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(coloc)
  library(meta)
})

cache_dir <- "00_data/mr_cache"
res_dir   <- "02_results"
ao        <- as.data.table(readRDS("00_data/opengwas_available_outcomes.rds"))

bcac_ids   <- c("ieu-a-1126", "ieu-a-1127", "ieu-a-1128")
finngen_ids <- c("finn-b-C3_BREAST", "finn-b-C3_BREAST_EXALLC")
out_labels <- c(
  "ieu-a-1126" = "BRCA_Overall",
  "ieu-a-1127" = "BRCA_ER_pos",
  "ieu-a-1128" = "BRCA_ER_neg",
  "finn-b-C3_BREAST" = "FinnGen_BREAST",
  "finn-b-C3_BREAST_EXALLC" = "FinnGen_BREAST_EXALLC"
)

safe_extract_ivs <- function(gid, p_thresh = 1e-5, r2_val = 0.001, kb_val = 10000, prefix = "iv_extra_") {
  rds_fp <- file.path(cache_dir, sprintf("%s%s_%g.rds", prefix, gid, p_thresh))
  if (file.exists(rds_fp)) return(readRDS(rds_fp))
  for (att in 1:3) {
    res <- tryCatch({
      Sys.sleep(0.35)
      TwoSampleMR::extract_instruments(gid, p1 = p_thresh, clump = TRUE, r2 = r2_val, kb = kb_val)
    }, error = function(e) {
      message("  Attempt ", att, " error for ", gid, ": ", conditionMessage(e))
      Sys.sleep(2.5 * att)
      NULL
    })
    if (!is.null(res) || att == 3) {
      if (is.null(res)) res <- data.frame()
      saveRDS(res, rds_fp)
      return(res)
    }
  }
}

safe_extract_outcomes <- function(snps, outcomes, cache_fp) {
  if (file.exists(cache_fp)) return(readRDS(cache_fp))
  snps <- unique(snps[!is.na(snps) & nzchar(snps)])
  if (length(snps) == 0) return(data.frame())
  chunks <- split(snps, ceiling(seq_along(snps) / 250))
  out_list <- list()
  for (ci in seq_along(chunks)) {
    for (oid in outcomes) {
      Sys.sleep(0.35)
      r <- tryCatch(
        TwoSampleMR::extract_outcome_data(snps = chunks[[ci]], outcomes = oid, proxies = TRUE, rsq = 0.8),
        error = function(e) { message("  Outcome chunk err ", oid, ": ", conditionMessage(e)); NULL }
      )
      if (!is.null(r) && nrow(r) > 0) out_list[[paste(oid, ci, sep = "_")]] <- r
    }
  }
  res <- if (length(out_list) > 0) as.data.frame(rbindlist(out_list, fill = TRUE)) else data.frame()
  saveRDS(res, cache_fp)
  res
}

# ==============================================================================
# Part 1: DMP_2022 (84) + FR02_2022 (36) 全量提取与 BCAC MR 检验
# ==============================================================================
cat("=== [1/2] DMP_2022 (84) + FR02_2022 (36) Independent Microbiome Cohorts ===\n")
aligned <- fread(file.path(res_dir, "aligned_microbiome_gwas.csv"))
extra_micro <- aligned[cohort %in% c("DMP_2022", "FR02_2022")]
cat("Aligned taxa: DMP_2022 =", sum(extra_micro$cohort == "DMP_2022"),
    "| FR02_2022 =", sum(extra_micro$cohort == "FR02_2022"), "\n")

iv_rows <- list()
for (i in seq_len(nrow(extra_micro))) {
  gid <- extra_micro$id[i]
  if (i %% 20 == 1 || i == nrow(extra_micro)) {
    cat(sprintf("  [%d/%d] IV extract %s (%s: %s)\n", i, nrow(extra_micro), gid, extra_micro$cohort[i], extra_micro$gwas_name[i]))
  }
  iv <- safe_extract_ivs(gid, p_thresh = 1e-5, prefix = "iv_dmpfr_")
  if (!is.null(iv) && nrow(iv) > 0) iv_rows[[gid]] <- iv
}
iv_extra_all <- as.data.table(rbindlist(iv_rows, fill = TRUE))
cat("Raw IVs extracted (P < 1e-5):", nrow(iv_extra_all), "across", length(unique(iv_extra_all$id.exposure)), "taxa\n")

if (nrow(iv_extra_all) > 0) {
  iv_extra_all[, F_stat := (beta.exposure / se.exposure)^2]
  iv_extra_strong <- iv_extra_all[F_stat > 10]
  cat("Strong IVs (F > 10):", nrow(iv_extra_strong), "across", length(unique(iv_extra_strong$id.exposure)), "taxa\n")

  out_bcac_extra <- safe_extract_outcomes(
    unique(iv_extra_strong$SNP), bcac_ids,
    file.path(cache_dir, "outcome_bcac_dmp_fr02.rds")
  )
  harm_extra <- TwoSampleMR::harmonise_data(as.data.frame(iv_extra_strong), out_bcac_extra, action = 2)
  saveRDS(harm_extra, file.path(cache_dir, "dmp_fr02_bcac_harmonised.rds"))

  mr_ex <- as.data.table(TwoSampleMR::mr(harm_extra, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median", "mr_egger_regression")))
  het_ex <- as.data.table(TwoSampleMR::mr_heterogeneity(harm_extra))[method == "Inverse variance weighted", .(id.exposure, id.outcome, Q_pval)]
  pl_ex  <- as.data.table(TwoSampleMR::mr_pleiotropy_test(harm_extra))[, .(id.exposure, id.outcome, pleiotropy_pval = pval)]

  mr_ex <- merge(mr_ex, extra_micro[, .(id.exposure = id, cohort, gwas_rank, gwas_name, matched_exact, matched_genus)], by = "id.exposure", all.x = TRUE)
  mr_ex <- merge(mr_ex, het_ex, by = c("id.exposure", "id.outcome"), all.x = TRUE)
  mr_ex <- merge(mr_ex, pl_ex,  by = c("id.exposure", "id.outcome"), all.x = TRUE)
  mr_ex[, `:=`(
    outcome_label = out_labels[id.outcome],
    OR = exp(b), OR_LCI95 = exp(b - 1.96 * se), OR_UCI95 = exp(b + 1.96 * se)
  )]
  fwrite(mr_ex, file.path(res_dir, "module3_dmp_fr02_bcac_mr_all_methods.csv"))

  prim_ex <- mr_ex[method %in% c("Inverse variance weighted", "Wald ratio")]
  prim_ex[, FDR := p.adjust(pval, method = "BH"), by = .(cohort, id.outcome)]
  fwrite(prim_ex, file.path(res_dir, "module3_dmp_fr02_bcac_mr_full_denominator.csv"))

  # Check replication specifically for our 6 MiBioGen causal taxa!
  core6 <- c("Adlercreutzia", "Desulfovibrio", "Parabacteroides", "Ruminococcaceae", "Dorea", "Blautia")
  rep_hits <- prim_ex[grepl(paste(core6, collapse = "|"), paste(gwas_name, trait_name = exposure), ignore.case = TRUE)]
  fwrite(rep_hits, file.path(res_dir, "module3_dmp_fr02_core6_replication.csv"))
  cat("\n=== Core 6 Causal Microbes tested in DMP_2022 / FR02_2022 ===\n")
  print(rep_hits[order(pval), .(cohort, gwas_rank, gwas_name, outcome_label, nsnp, OR = round(OR, 3), pval = signif(pval, 3), FDR = signif(FDR, 3))])

  cat("\n=== Any Nominal P < 0.05 in DMP_2022 / FR02_2022 ===\n")
  print(prim_ex[pval < 0.05][order(pval), .(cohort, gwas_rank, gwas_name, outcome_label, nsnp, OR = round(OR, 3), pval = signif(pval, 3), FDR = signif(FDR, 3))][1:min(25, .N)])

  # Cross-cohort Meta for overlapping taxa matching the 6 MiBioGen taxa + BCAC outcomes
  mib_prim <- fread(file.path(res_dir, "mibiogen_bcac_mr_primary.csv"))
  mib_c6 <- mib_prim[gwas_name %in% core6, .(cohort = "MiBioGen_2021", gwas_name, outcome_label, id.exposure, b, se, pval, nsnp)]
  rep_c6 <- rep_hits[, .(cohort, gwas_name, outcome_label, id.exposure, b, se, pval, nsnp)]
  # Standardize genus matching name
  rep_c6[, core_match := sapply(gwas_name, function(x) { m <- core6[sapply(core6, function(c6) grepl(c6, x, ignore.case = TRUE))]; if (length(m) > 0) m[1] else x })]
  mib_c6[, core_match := gwas_name]
  comb_c6 <- rbind(mib_c6, rep_c6, fill = TRUE)
  meta_rows <- list()
  for (cm in unique(rep_c6$core_match)) {
    for (ol in unique(comb_c6$outcome_label)) {
      sub <- comb_c6[core_match == cm & outcome_label == ol & !is.na(b) & !is.na(se) & se > 0]
      # Average or pick top species within same cohort if multiple DMP species, or keep unique cohorts:
      if (length(unique(sub$cohort)) >= 2) {
        # For clean 1-per-cohort meta, pick exact genus match or smallest SE (largest power) per cohort
        sub1 <- sub[order(se)][, .SD[1], by = cohort]
        mm <- meta::metagen(TE = sub1$b, seTE = sub1$se, studlab = paste0(sub1$cohort, ":", sub1$gwas_name), sm = "OR", common = TRUE, random = TRUE)
        meta_rows[[paste(cm, ol)]] <- data.table(
          core_taxon = cm, outcome_label = ol, k_cohorts = nrow(sub1),
          studies = paste(paste0(sub1$cohort, "(", sub1$gwas_name, ")"), collapse = "; "),
          FE_OR = round(exp(mm$TE.common), 4), FE_LCI = round(exp(mm$lower.common), 4), FE_UCI = round(exp(mm$upper.common), 4), FE_pval = signif(mm$pval.common, 4),
          RE_OR = round(exp(mm$TE.random), 4), RE_pval = signif(mm$pval.random, 4), I2_pct = round(mm$I2 * 100, 1)
        )
      }
    }
  }
  meta_micro <- if (length(meta_rows) > 0) rbindlist(meta_rows) else data.table()
  fwrite(meta_micro, file.path(res_dir, "module3_microbiome_cross_cohort_meta.csv"))
  cat("\n=== Cross-Cohort Microbiome Meta (MiBioGen + DMP + FR02) ===\n")
  print(meta_micro[order(FE_pval)][1:min(15, .N)])
}

# ==============================================================================
# Part 2: CXCL8 / IL-8 Plasma pQTL MR & Coloc + PPARG pQTL Catalog Audit
# ==============================================================================
cat("\n=== [2/2] Plasma pQTL Search & CXCL8/IL-8 Protein-Level MR ===\n")
pparg_pqtl <- ao[grepl("PPARG|NR1C3|peroxisome proliferator-activated receptor gamma", trait, ignore.case = TRUE) & !grepl("^eqtl-", id)]
cat("OpenGWAS PPARG plasma pQTL matches:", nrow(pparg_pqtl), "\n")
if (nrow(pparg_pqtl) > 0) print(pparg_pqtl[, .(id, trait, author, year, sample_size)])

il8_ids <- c("ebi-a-GCST90011994", "prot-a-749", "ebi-a-GCST004445", "prot-c-3447_64_2", "ebi-a-GCST90000480")
pqtl_iv_list <- list()
for (pid in il8_ids) {
  # Try 5e-8 first, fallback 1e-5 if < 2 SNPs
  iv_5e8 <- safe_extract_ivs(pid, p_thresh = 5e-8, prefix = "iv_pqtl_")
  if (nrow(iv_5e8) >= 1) {
    iv_5e8$thresh_used <- "5e-8"
    pqtl_iv_list[[paste0(pid, "_5e8")]] <- iv_5e8
  }
  iv_1e5 <- safe_extract_ivs(pid, p_thresh = 1e-5, prefix = "iv_pqtl_")
  if (nrow(iv_1e5) >= 1) {
    iv_1e5$thresh_used <- "1e-5"
    pqtl_iv_list[[paste0(pid, "_1e5")]] <- iv_1e5
  }
}

# Combine genome-wide (5e-8 primary; plus 1e-5 & strict cis chr4:72.7-74.7Mb)
pqtl_rows <- list()
if (length(pqtl_iv_list) > 0) {
  all_pqtl_snps <- unique(unlist(lapply(pqtl_iv_list, function(x) x$SNP)))
  out_pqtl <- safe_extract_outcomes(
    all_pqtl_snps, c(bcac_ids, finngen_ids),
    file.path(cache_dir, "outcome_bcac_finngen_il8_pqtl.rds")
  )
  for (nm in names(pqtl_iv_list)) {
    ivd <- as.data.table(pqtl_iv_list[[nm]])
    ivd[, F_stat := (beta.exposure / se.exposure)^2]
    ivd <- ivd[F_stat > 10]
    if (nrow(ivd) == 0) next
    # Run BOTH Full (genome-wide) and cis-only (chr4, 72.74-74.74 Mb around CXCL8)
    subsets <- list(
      All_IVs = ivd,
      Cis_1Mb_CXCL8 = ivd[as.character(chr.exposure) == "4" & as.numeric(pos.exposure) >= 72740568 & as.numeric(pos.exposure) <= 74743716]
    )
    for (sname in names(subsets)) {
      sd <- subsets[[sname]]
      if (nrow(sd) == 0) next
      hm <- TwoSampleMR::harmonise_data(as.data.frame(sd), out_pqtl, action = 2)
      mr_p <- as.data.table(TwoSampleMR::mr(hm, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median", "mr_egger_regression")))
      stg_p <- tryCatch(as.data.table(TwoSampleMR::directionality_test(hm)), error = function(e) data.table())
      mr_p[, `:=`(
        pqtl_set = nm,
        iv_scope = sname,
        thresh = sd$thresh_used[1],
        outcome_label = out_labels[id.outcome],
        OR = exp(b), OR_LCI95 = exp(b - 1.96 * se), OR_UCI95 = exp(b + 1.96 * se)
      )]
      if (nrow(stg_p) > 0) mr_p <- merge(mr_p, stg_p[, .(id.exposure, id.outcome, steiger_correct = correct_causal_direction)], by = c("id.exposure", "id.outcome"), all.x = TRUE)
      pqtl_rows[[paste(nm, sname)]] <- mr_p
    }
  }
}
pqtl_dt <- if (length(pqtl_rows) > 0) rbindlist(pqtl_rows, fill = TRUE) else data.table(note = "No IL-8 pQTL IVs matched")
pqtl_dt <- merge(pqtl_dt, ao[, .(id.exposure = id, trait, author, year, pQTL_n = sample_size)], by = "id.exposure", all.x = TRUE)
fwrite(pqtl_dt, file.path(res_dir, "planA_pqtl_il8_mr.csv"))
cat("\n=== CXCL8 / IL-8 Plasma pQTL MR Results (Primary Methods) ===\n")
if ("method" %in% names(pqtl_dt)) {
  print(pqtl_dt[method %in% c("Inverse variance weighted", "Wald ratio") & outcome_label %like% "BRCA",
                .(id.exposure, author, pQTL_n, thresh, iv_scope, outcome_label, nsnp, OR = round(OR, 3), pval = signif(pval, 3))])
}

# Colocalization for largest IL-8 pQTL (ebi-a-GCST90011994, n=21758) cis-CXCL8 region vs BCAC ER-
cat("\n=== Attempting cis-pQTL Coloc for SCALLOP IL-8 (ebi-a-GCST90011994) vs BCAC ER- ===\n")
reg_il8_e <- tryCatch(as.data.table(ieugwasr::associations(variants = "4:73240568-74243716", id = "ebi-a-GCST90011994", proxies = 0)), error = function(e) NULL)
reg_il8_o <- tryCatch(as.data.table(ieugwasr::associations(variants = "4:73240568-74243716", id = "ieu-a-1128", proxies = 0)), error = function(e) NULL)
if (!is.null(reg_il8_e) && nrow(reg_il8_e) > 20 && !is.null(reg_il8_o) && nrow(reg_il8_o) > 20) {
  e2 <- reg_il8_e[!is.na(rsid) & !is.na(beta) & se > 0][!duplicated(rsid)]
  o2 <- reg_il8_o[!is.na(rsid) & !is.na(beta) & se > 0][!duplicated(rsid)]
  mg <- merge(e2, o2, by = "rsid", suffixes = c(".e", ".o"))
  if (nrow(mg) >= 20) {
    maf <- ifelse(!is.na(mg$eaf.e) & mg$eaf.e > 0 & mg$eaf.e < 1, pmin(mg$eaf.e, 1 - mg$eaf.e), 0.2)
    maf[is.na(maf) | maf <= 0 | maf >= 0.5] <- 0.2
    cres <- coloc::coloc.abf(
      dataset1 = list(snp = mg$rsid, beta = mg$beta.e, varbeta = mg$se.e^2, N = 21758, MAF = maf, type = "quant"),
      dataset2 = list(snp = mg$rsid, beta = mg$beta.o, varbeta = mg$se.o^2, s = 21468 / (21468 + 100594), type = "cc")
    )
    print(cres$summary)
    fwrite(as.data.table(as.list(cres$summary))[, `:=`(n_snps = nrow(mg), pair = "SCALLOP_IL8_pQTL_vs_BCAC_ERneg")],
           file.path(res_dir, "planA_pqtl_il8_coloc.csv"))
  } else {
    cat("Too few merged cis SNPs in pQTL region:", nrow(mg), "\n")
  }
} else {
  cat("Region query returned limited/no SNPs in IL-8 cis window; checking lead SNP.\n")
}
cat("Module 3 extra analyses complete!\n")
