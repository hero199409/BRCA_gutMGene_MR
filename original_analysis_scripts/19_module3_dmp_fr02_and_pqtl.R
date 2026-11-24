# 19_module3_dmp_fr02_and_pqtl.R
# 模块 3 加餐：
# 1) DMP_2022 (84个) + FR02_2022 (36个) 独立菌群队列全量 MR 复现与跨队列 Meta
# 2) OpenGWAS 血浆 IL-8 / CXCL8 pQTL 蛋白层两样本 MR + 共定位 (及 PPARG pQTL 检索审计)

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(coloc)
  library(meta)
})

cache_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/mr_cache")
res_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")
ao        <- as.data.table(readRDS(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/opengwas_available_outcomes.rds")))

set.seed(2026)
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
  if (mr_offline) stop("Missing instrument cache: ", rds_fp)
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
  if (mr_offline) stop("Missing outcome cache: ", cache_fp)
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

# Compare genome-wide and cis sets at each available threshold.
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
    ivd[, is_cis_GRCh37 := sub("^chr", "", as.character(chr.exposure)) == "4" & between(as.numeric(pos.exposure), 73606223, 75609433)]
    pqtl_iv_list[[nm]] <- ivd

    if (nrow(ivd) == 0) next
    # Run BOTH Full (genome-wide) and cis-only (GRCh37 CXCL8 chr4:74,606,223-74,609,433 ±1Mb = 73.61-75.61 Mb)
    subsets <- list(
      All_IVs = ivd,
      Cis_1Mb_CXCL8 = ivd[is_cis_GRCh37 == TRUE]
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
pqtl_dt[, primary_estimator := method %in% c("Inverse variance weighted", "Wald ratio")]
pqtl_dt[primary_estimator == TRUE, FDR_within_scope_threshold := p.adjust(pval, "BH"), by = .(iv_scope, thresh)]
fwrite(pqtl_dt, file.path(res_dir, "planA_pqtl_il8_mr.csv"))
fwrite(rbindlist(pqtl_iv_list, idcol = "pqtl_set", fill = TRUE), file.path(res_dir, "planA_pqtl_il8_iv_audit.csv"))
cat("\n=== CXCL8 / IL-8 Plasma pQTL MR Results (Primary Methods) ===\n")
if ("method" %in% names(pqtl_dt)) {
  print(pqtl_dt[method %in% c("Inverse variance weighted", "Wald ratio") & outcome_label %like% "BRCA",
                .(id.exposure, author, pQTL_n, thresh, iv_scope, outcome_label, nsnp, OR = round(OR, 3), pval = signif(pval, 3))])
}

# Same-rank primary meta-analysis. Species substitution is reported separately.
mib_full <- fread(file.path(res_dir, "mibiogen_bcac_mr_primary.csv"))[, cohort := "MiBioGen_2021"]
dmp_full <- fread(file.path(res_dir, "module3_dmp_fr02_bcac_mr_full_denominator.csv"))
core_taxa <- c("Ruminococcaceae", "Adlercreutzia", "Parabacteroides", "Blautia", "Dorea", "Desulfovibrio")
make_meta <- function(taxon, scope) {
  rank <- if (taxon == "Ruminococcaceae") "family" else "genus"
  d <- rbindlist(list(mib_full, dmp_full), fill = TRUE)
  d[, same_rank := gwas_rank == rank & gwas_name == taxon]
  d[, species := gwas_rank == "species" & startsWith(gwas_name, paste0(taxon, " "))]
  d <- d[(same_rank | (scope == "species_substitution" & species)) & is.finite(b) & is.finite(se) & se > 0]
  if (!nrow(d)) return(NULL)
  # Select a matched genus/family first. For an unmatched cohort, use accession order
  # for the species sensitivity, so selection does not depend on an outcome statistic.
  setorderv(d, c("same_rank", "id.exposure"), c(-1, 1))
  d <- unique(d, by = c("cohort", "outcome_label"))
  d[, {
    sources <- paste(paste0(cohort, " (", gwas_name, ", ", id.exposure, ")"), collapse = " | ")
    if (.N < 2L) {
      list(core_taxon = taxon, analysis_scope = scope, k_cohorts = .N, studies = sources,
           FE_OR = NA_real_, FE_LCI = NA_real_, FE_UCI = NA_real_, FE_pval = NA_real_,
           RE_OR = NA_real_, RE_LCI = NA_real_, RE_UCI = NA_real_, RE_pval = NA_real_, I2_pct = NA_real_)
    } else {
      fit <- meta::metagen(b, se, studlab = cohort, sm = "OR", common = TRUE, random = TRUE, method.tau = "REML")
      list(core_taxon = taxon, analysis_scope = scope, k_cohorts = .N, studies = sources,
           FE_OR = exp(fit$TE.common), FE_LCI = exp(fit$lower.common), FE_UCI = exp(fit$upper.common), FE_pval = fit$pval.common,
           RE_OR = exp(fit$TE.random), RE_LCI = exp(fit$lower.random), RE_UCI = exp(fit$upper.random), RE_pval = fit$pval.random,
           I2_pct = 100 * fit$I2)
    }
  }, by = outcome_label]
}
strict <- rbindlist(lapply(core_taxa, make_meta, scope = "same_rank"))
species <- rbindlist(lapply(core_taxa, make_meta, scope = "species_substitution"))
meta_audit <- rbindlist(list(strict, species))
meta_audit[, RE_FDR := p.adjust(RE_pval, "BH"), by = analysis_scope]
fwrite(meta_audit, file.path(res_dir, "module3_microbiome_cross_cohort_meta_strict.csv"))
fwrite(meta_audit[analysis_scope == "same_rank"], file.path(res_dir, "module3_microbiome_cross_cohort_meta.csv"))

pq_eid <- "ebi-a-GCST90011994"
pq_oid <- "ieu-a-1128"
region <- "4:74106223-75109433"
reg_cache <- file.path(cache_dir, "coloc_scallop_il8_bcac_erneg_grch37.rds")
ld_cache <- file.path(cache_dir, "ld_scallop_interval_cxcl8_EUR.rds")
if (file.exists(reg_cache)) {
  reg <- readRDS(reg_cache)
} else if (!mr_offline) {
  reg <- list(e = ieugwasr::associations(region, pq_eid, proxies = 0),
              o = ieugwasr::associations(region, pq_oid, proxies = 0))
  if (!is.data.frame(reg$e) || !is.data.frame(reg$o) || min(nrow(reg$e), nrow(reg$o)) < 20L) {
    stop("Regional response lacks sufficient summary data.")
  }
  saveRDS(reg, reg_cache)
} else {
  reg <- NULL
}
reg_status <- "GRCh37 summary statistics not cached"
if (!is.null(reg)) {
  e <- as.data.table(reg$e)[is.finite(beta) & is.finite(se) & se > 0 & is.finite(eaf) & eaf > 0 & eaf < 1]
  o <- as.data.table(reg$o)[is.finite(beta) & is.finite(se) & se > 0]
  m <- merge(e[!duplicated(rsid)], o[!duplicated(rsid)], by = "rsid", suffixes = c(".e", ".o"))
  if (nrow(m) >= 20L) {
    bm <- ao[id == pq_oid]
    ne <- ao[id == pq_eid, sample_size][1]
    pq_coloc <- rbindlist(lapply(c(1e-6, 1e-5, 5e-5), function(prior) {
      fit <- coloc::coloc.abf(
        dataset1 = list(snp = m$rsid, beta = m$beta.e, varbeta = m$se.e^2, N = ne,
                        MAF = pmin(m$eaf.e, 1 - m$eaf.e), type = "quant"),
        dataset2 = list(snp = m$rsid, beta = m$beta.o, varbeta = m$se.o^2,
                        s = bm$ncase / (bm$ncase + bm$ncontrol), type = "cc"),
        p1 = 1e-4, p2 = 1e-4, p12 = prior)
      as.data.table(as.list(fit$summary))[, `:=`(pair = "SCALLOP_IL8_pQTL_vs_BCAC_ERneg", n_snps = nrow(m),
        region_hg19 = region, p12_prior = prior,
        rs7655660_present = "rs7655660" %in% m$rsid)]
    }))
    fwrite(pq_coloc, file.path(res_dir, "planA_pqtl_il8_coloc_sensitivity.csv"))
    fwrite(pq_coloc[p12_prior == 1e-5], file.path(res_dir, "planA_pqtl_il8_coloc.csv"))
    reg_status <- "Completed with GRCh37 regional summary statistics"
  } else reg_status <- "Fewer than 20 eligible shared SNPs"
}
ld_snps <- c("rs7655660", "rs74327761", "rs1951703", "rs10010970")
if (file.exists(ld_cache)) {
  ld <- readRDS(ld_cache)
} else if (!mr_offline) {
  ld <- ieugwasr::ld_matrix(variants = ld_snps, pop = "EUR", with_alleles = FALSE)
  if (is.matrix(ld) && length(ld)) saveRDS(ld, ld_cache)
} else ld <- NULL
if (is.matrix(ld) && length(ld)) {
  ld_tab <- as.data.table(as.table(ld^2))
  setnames(ld_tab, c("snp1", "snp2", "r2"))
  fwrite(ld_tab, file.path(res_dir, "planA_pqtl_il8_ld_check.csv"))
}
fwrite(data.table(analysis = c("SCALLOP coloc", "EUR LD"), region = c(region, NA_character_),
                  status = c(reg_status, if (is.null(ld)) "EUR LD matrix not cached" else "See LD table")),
       file.path(res_dir, "planA_pqtl_il8_availability.csv"))
cat("Module 3 extra analyses complete. ", reg_status, "\n")
