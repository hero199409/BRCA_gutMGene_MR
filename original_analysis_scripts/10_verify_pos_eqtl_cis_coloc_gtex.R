# 10_verify_pos_eqtl_cis_coloc_gtex.R
# 严格基于纯 cis-IV (±1Mb) 的 eQTLGen MR 结果，对全部 cis-MR 阳性基因跑：
# 1) cis-IV 敏感性（异质性/多效性/Steiger 方向检验）
# 2) coloc 贝叶斯共定位 (eQTLGen vs BCAC, ±500kb 窗口)
# 3) FinnGen 独立队列验证
# 4) GTEx v8 Breast_Mammary_Tissue API v2 乳腺组织 eQTL 查证

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
  library(coloc)
  library(jsonlite)
})

cache_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/mr_cache")
res_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

ao <- as.data.table(readRDS(file.path(project_dir, "00_data/opengwas_available_outcomes.rds")))
set.seed(2026)
univ <- fread(file.path(res_dir, "module3_candidate_host_genes_universe.csv"))
ivs_all <- as.data.table(readRDS(file.path(cache_dir, "eqtlgen_all49_exp_ivs.rds")))
all_tested_genes <- unique(as.character(ivs_all$gene_symbol))

# 1. 严格标注全量 41 个有 IV 基因的 cis/trans 状态
# 使用权威 Ensembl GRCh37 (hg19) 坐标以与 eQTLGen/BCAC (GRCh37) 严格对齐
grch37_fp <- file.path(res_dir, "module3_grch37_gene_coords_lookup.csv")
ivs_all[, ensembl_id := sub("^eqtl-a-", "", id.exposure)]
if (file.exists(grch37_fp)) {
  g37 <- fread(grch37_fp)
} else {
  g37 <- data.table(gene_symbol = character(), ensembl_id = character(), chr_hg19 = character(), start_hg19 = numeric(), end_hg19 = numeric())
}
missing_ids <- setdiff(unique(ivs_all$ensembl_id), g37$ensembl_id)
if (!mr_offline && length(missing_ids) > 0 && requireNamespace("jsonlite", quietly = TRUE)) {
  extra_g37 <- rbindlist(lapply(missing_ids, function(eid) {
    sym <- ivs_all[ensembl_id == eid, gene_symbol][1]
    cmd <- sprintf('curl.exe -4 -s %s --connect-timeout 6 --max-time 12 "https://grch37.rest.ensembl.org/lookup/id/%s?content-type=application/json"', curl_proxy_flag, eid)
    js <- tryCatch(jsonlite::fromJSON(paste(suppressWarnings(system(cmd, intern = TRUE)), collapse = "")), error = function(e) NULL)
    if (!is.null(js) && identical(js$assembly_name, "GRCh37") && !is.null(js$seq_region_name)) {
      data.table(gene_symbol = sym, ensembl_id = eid, ensembl_display_name = js$display_name,
                 chr_hg19 = as.character(js$seq_region_name), start_hg19 = as.numeric(js$start),
                 end_hg19 = as.numeric(js$end), strand_hg19 = as.numeric(js$strand), assembly = "GRCh37")
    } else NULL
  }), fill = TRUE)
  if (nrow(extra_g37) > 0) {
    g37 <- rbindlist(list(g37, extra_g37), fill = TRUE)
    fwrite(g37, grch37_fp)
  }
}
stopifnot(all(g37$assembly == "GRCh37"), !anyDuplicated(g37$ensembl_id),
          all(ivs_all$ensembl_id %in% g37$ensembl_id))
gene_coords_hg19 <- unique(g37[, .(ensembl_id, gene_chr = as.character(chr_hg19),
                                   gene_start = as.numeric(start_hg19), gene_end = as.numeric(end_hg19))])
ivs_annot <- merge(ivs_all, gene_coords_hg19, by = "ensembl_id", all.x = TRUE)
ivs_annot[, chr_clean := sub("^chr", "", as.character(chr.exposure))]
ivs_annot[, pos_num := as.numeric(pos.exposure)]
ivs_annot[, is_cis_1mb := (chr_clean == gene_chr) &
            (pos_num >= (gene_start - 1e6)) &
            (pos_num <= (gene_end + 1e6))]
ivs_annot[, dist_to_gene_body_bp := fifelse(chr_clean != gene_chr, NA_real_,
                                   fifelse(pos_num < gene_start, gene_start - pos_num,
                                   fifelse(pos_num > gene_end, pos_num - gene_end, 0)))]
fwrite(ivs_annot[, .(gene_symbol, id.exposure, ensembl_id, SNP, chr_clean, pos_num, gene_chr, gene_start, gene_end,
                     dist_to_gene_body_bp, is_cis_1mb,
                     effect_allele.exposure, other_allele.exposure, eaf.exposure,
                     beta.exposure, se.exposure, pval.exposure, F_stat)],
       file.path(res_dir, "module3_all41_genes_iv_cis_audit.csv"))
fwrite(ivs_annot[order(gene_symbol, SNP),
                 .(gene_symbol, id.exposure, ensembl_id, SNP,
                   chr_snp_hg19 = chr_clean, pos_snp_hg19 = pos_num,
                   gene_chr_hg19 = gene_chr, gene_start_hg19 = gene_start, gene_end_hg19 = gene_end,
                   dist_to_gene_body_bp, is_cis_1mb_hg19 = is_cis_1mb,
                   effect_allele.exposure, other_allele.exposure, eaf.exposure,
                   beta.exposure, se.exposure, pval.exposure, F_stat)],
       file.path(res_dir, "module3_all41_genes_157_iv_hg19_audit.csv"))

# 2. 用纯 cis-IV 跑全量分母 MR + 异质性 + 多效性 + Steiger
cis_ivs <- ivs_annot[is_cis_1mb == TRUE]
out_bcac <- readRDS(file.path(cache_dir, "outcome_bcac_eqtlgen_all49.rds"))
dat_cis <- TwoSampleMR::harmonise_data(as.data.frame(cis_ivs), as.data.frame(out_bcac), action = 2)
saveRDS(dat_cis, file.path(cache_dir, "eqtlgen_cis_only_bcac_harmonised.rds"))

mr_cis  <- as.data.table(TwoSampleMR::mr(dat_cis, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median", "mr_egger_regression")))
het_cis <- as.data.table(TwoSampleMR::mr_heterogeneity(dat_cis))
plei_cis <- as.data.table(TwoSampleMR::mr_pleiotropy_test(dat_cis))
steiger_cis <- data.table() # Binary-outcome liability-scale directionality is computed below.

mr_cis[, OR := exp(b)]
mr_cis[, OR_LCI95 := exp(b - 1.96 * se)]
mr_cis[, OR_UCI95 := exp(b + 1.96 * se)]

het_ivw <- het_cis[method == "Inverse variance weighted", .(id.exposure, id.outcome, Q, Q_pval)]
plei_dt <- plei_cis[, .(id.exposure, id.outcome, egger_intercept, pleiotropy_pval = pval)]
dir_dt  <- mr_cis[, .(all_dir_consistent = all(b > 0) || all(b < 0)), by = .(id.exposure, id.outcome)]

mr_cis <- merge(mr_cis, univ[, .(id.exposure = eqtl_id, gene_symbol, upstream_microbes, upstream_metabolites, alterations)],
                by = "id.exposure", all.x = TRUE)
mr_cis <- merge(mr_cis, het_ivw, by = c("id.exposure", "id.outcome"), all.x = TRUE)
mr_cis <- merge(mr_cis, plei_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)
mr_cis <- merge(mr_cis, dir_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)
if (nrow(steiger_cis) > 0) {
  mr_cis <- merge(mr_cis, steiger_cis[, .(id.exposure, id.outcome, correct_causal_direction, steiger_pval)],
                  by = c("id.exposure", "id.outcome"), all.x = TRUE)
}
mr_cis[, outcome_label := fifelse(id.outcome == "ieu-a-1126", "BRCA_Overall",
                           fifelse(id.outcome == "ieu-a-1127", "BRCA_ER_pos", "BRCA_ER_neg"))]

prim_cis <- mr_cis[method %in% c("Inverse variance weighted", "Wald ratio")]
prim_cis[, FDR := p.adjust(pval, method = "BH"), by = id.outcome]

fwrite(mr_cis, file.path(res_dir, "module3_eqtlgen_cis_only_mr_all_methods.csv"))
for (name in c("module3_eqtlgen_cis_only_mr_primary_full_denominator.csv", "module3_eqtlgen_cis_only_mr_primary.csv")) {
  fwrite(prim_cis[order(outcome_label, pval)], file.path(res_dir, name))
}

sig_cis <- prim_cis[pval < 0.05][order(outcome_label, pval)]
fwrite(sig_cis, file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))

pos_genes <- unique(as.character(sig_cis$gene_symbol))
cat("=== 纯 cis-IV eQTLGen MR 汇总 ===\n")
cat("Total cis-IV genes tested:", length(unique(prim_cis$gene_symbol)),
    "| Pairs tested:", nrow(prim_cis),
    "| Nominal P<0.05 pairs:", nrow(sig_cis),
    "| Unique positive cis-genes:", length(pos_genes), "(", paste(pos_genes, collapse=", "), ")\n")
if (!"correct_causal_direction" %in% names(sig_cis)) sig_cis[, correct_causal_direction := NA]
print(sig_cis[, .(outcome_label, gene_symbol, method, nsnp, OR=round(OR,3), pval=signif(pval,3), FDR=signif(FDR,3),
                  all_dir_consistent, Q_pval=signif(Q_pval,2), pleiotropy_pval=signif(pleiotropy_pval,2),
                  correct_causal_direction, upstream_metabolites, alterations)])

# 3. 仅对纯 cis-MR 阳性基因跑 coloc 共定位 (±500kb 围绕基因顺式 lead SNP)
cat("\n=== Running coloc for cis-MR positive pairs ===\n")
bcac_meta <- ao[id %in% c("ieu-a-1126", "ieu-a-1127", "ieu-a-1128"), .(id, ncase, ncontrol, sample_size)]
coloc_rows <- list()
coloc_sensitivity <- list()

for (i in seq_len(nrow(sig_cis))) {
  eid <- sig_cis$id.exposure[i]
  oid <- sig_cis$id.outcome[i]
  sym <- sig_cis$gene_symbol[i]
  olb <- sig_cis$outcome_label[i]
  
  sub_iv <- cis_ivs[gene_symbol == sym][order(pval.exposure)][1]
  lead_chr <- sub_iv$chr_clean
  lead_pos <- sub_iv$pos_num
  win_start <- max(1, floor(lead_pos - 500000))
  win_end   <- ceiling(lead_pos + 500000)
  region_str <- sprintf("%s:%d-%d", lead_chr, win_start, win_end)
  
  coloc_cache <- file.path(cache_dir, paste0("coloc_", eid, "_", oid, "_", lead_chr, "_", win_start, ".rds"))
  cat(sprintf("[%d/%d] Coloc %s (%s) vs %s at %s\n", i, nrow(sig_cis), sym, eid, olb, region_str))
  
  reg_data <- NULL
  if (file.exists(coloc_cache)) {
    reg_data <- readRDS(coloc_cache)
  } else if (!mr_offline) {
    Sys.sleep(0.45)
    reg_data <- tryCatch({
      e_reg <- ieugwasr::associations(variants = region_str, id = eid, proxies = 0)
      Sys.sleep(0.45)
      o_reg <- ieugwasr::associations(variants = region_str, id = oid, proxies = 0)
      list(e_reg = as.data.table(e_reg), o_reg = as.data.table(o_reg))
    }, error = function(e) {
      message("  Coloc query failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(reg_data)) saveRDS(reg_data, coloc_cache)
  }
  
  if (!is.null(reg_data) && nrow(reg_data$e_reg) > 10 && nrow(reg_data$o_reg) > 10) {
    e_dt <- reg_data$e_reg[!is.na(rsid) & !is.na(beta) & !is.na(se) & se > 0][!duplicated(rsid)]
    o_dt <- reg_data$o_reg[!is.na(rsid) & !is.na(beta) & !is.na(se) & se > 0][!duplicated(rsid)]
    merged <- merge(e_dt, o_dt, by = "rsid", suffixes = c(".e", ".o"))
    if (nrow(merged) >= 20) {
      merged <- merged[is.finite(eaf.e) & eaf.e > 0 & eaf.e < 1]
      if (nrow(merged) < 20L) next
      n_e <- ao[id == eid, sample_size][1]
      bm <- bcac_meta[id == oid]
      s_prop <- bm$ncase[1] / (bm$ncase[1] + bm$ncontrol[1])
      row_key <- paste(eid, oid, sep = "_")
      priors <- c(1e-6, 1e-5, 5e-5)
      all_priors <- rbindlist(lapply(priors, function(prior) {
        fit <- coloc::coloc.abf(
          dataset1 = list(snp = merged$rsid, beta = merged$beta.e, varbeta = merged$se.e^2,
                          N = n_e, MAF = pmin(merged$eaf.e, 1 - merged$eaf.e), type = "quant"),
          dataset2 = list(snp = merged$rsid, beta = merged$beta.o, varbeta = merged$se.o^2,
                          s = s_prop, type = "cc"), p1 = 1e-4, p2 = 1e-4, p12 = prior
        )
        sm <- fit$summary
        data.table(gene_symbol = sym, id.exposure = eid, outcome_label = olb,
                   region = region_str, n_shared_snps = nrow(merged), p12_prior = prior,
                   PP.H0 = sm[["PP.H0.abf"]], PP.H1 = sm[["PP.H1.abf"]],
                   PP.H2 = sm[["PP.H2.abf"]], PP.H3 = sm[["PP.H3.abf"]], PP.H4 = sm[["PP.H4.abf"]])
      }))
      coloc_sensitivity[[row_key]] <- all_priors
      coloc_rows[[row_key]] <- all_priors[p12_prior == 1e-5]
    }
  }
}

coloc_dt <- if (length(coloc_rows) > 0) rbindlist(coloc_rows) else data.table()
fwrite(coloc_dt, file.path(res_dir, "module3_cis_positive_genes_coloc.csv"))
cat("\n=== 纯 cis-MR 阳性基因 coloc 共定位结果 ===\n")
print(coloc_dt)
fwrite(rbindlist(coloc_sensitivity), file.path(res_dir, "module3_cis_coloc_prior_sensitivity.csv"))

# 4. FinnGen 独立队列验证纯 cis-MR 阳性基因
cat("\n=== Running FinnGen validation for cis-MR positive genes ===\n")
pos_exp_ivs <- as.data.frame(cis_ivs[gene_symbol %in% pos_genes])
finn_cis_cache <- file.path(cache_dir, "outcome_finngen_cis_pos_eqtl.rds")
if (mr_offline && !file.exists(finn_cis_cache)) stop("FinnGen summary cache is required in offline mode.")
if (file.exists(finn_cis_cache)) {
  out_finn_eqtl <- readRDS(finn_cis_cache)
} else {
  Sys.sleep(0.5)
  out_finn_eqtl <- TwoSampleMR::extract_outcome_data(
    snps = unique(pos_exp_ivs$SNP),
    outcomes = c("finn-b-C3_BREAST_EXALLC", "finn-b-C3_BREAST"),
    proxies = TRUE, rsq = 0.8
  )
  saveRDS(out_finn_eqtl, finn_cis_cache)
}
dat_fe <- TwoSampleMR::harmonise_data(pos_exp_ivs, out_finn_eqtl, action = 2)
mr_fe  <- as.data.table(TwoSampleMR::mr(dat_fe, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median")))
mr_fe[, OR := exp(b)]
mr_fe <- merge(mr_fe, unique(cis_ivs[, .(id.exposure, gene_symbol)]), by = "id.exposure", all.x = TRUE)
fwrite(mr_fe, file.path(res_dir, "module3_cis_positive_genes_finngen_validation.csv"))
cat("FinnGen validation (primary):\n")
print(mr_fe[method %in% c("Inverse variance weighted", "Wald ratio")][order(pval), .(id.outcome, gene_symbol, nsnp, OR = round(OR, 3), pval = signif(pval, 3))])

# 5. GTEx v8 Breast_Mammary_Tissue API v2 查证 (仅查纯 cis-MR 阳性基因)
cat("\n=== Querying GTEx v8 Breast_Mammary_Tissue eQTLs via API v2 ===\n")
gtex_rows <- list()
if (!mr_offline) {
  get_gtex <- function(url) {
    response <- httr::GET(url, httr::timeout(20))
    httr::stop_for_status(response)
    jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))$data
  }
  for (sym in pos_genes) {
    ref <- get_gtex(paste0("https://gtexportal.org/api/v2/reference/gene?gencodeVersion=v26&geneId=", sym))
    gid <- ref$gencodeId[1]
    if (length(gid) != 1L || is.na(gid)) next
    eq <- get_gtex(paste0("https://gtexportal.org/api/v2/association/independentEqtl?datasetId=gtex_v8&tissueSiteDetailId=Breast_Mammary_Tissue&gencodeId=", gid))
    if (is.data.frame(eq) && nrow(eq)) gtex_rows[[sym]] <- as.data.table(eq)[, gene_symbol := sym]
  }
  if (length(gtex_rows)) fwrite(rbindlist(gtex_rows, fill = TRUE), file.path(res_dir, "module3_cis_positive_genes_gtex_breast_eqtl.csv"))
} else {
  message("Offline mode: existing GTEx results retained.")
}

# ==============================================================================
# [5/5] 补充稳健性审计：
#   (a) 40->36 cis-IV 协调丢失审计
#   (b) ER+ vs ER- 亚型效应异质性检验 (Cochran Q P_het)
#   (c) Coloc 先验灵敏度 (p12 = 1e-6, 1e-5, 5e-5)
#   (d) GWAS 元数据、队列重叠与 Steiger 患病率敏感性 (1%, 5%, 12%)
# ==============================================================================
cat("\n=== [5/5] Cis Drop Audit, ER+/ER- P_het, Coloc Prior & Steiger Sensitivity ===\n")

# (a) 40 cis-IV genes -> 36 harmonised genes audit
ivs_cis_df <- cis_ivs
dat_harm_cis_dt <- as.data.table(dat_cis)
prim_cis_dt <- prim_cis
sig_cis_dt <- sig_cis
out_all_dt <- as.data.table(out_bcac)

cis_genes_40 <- unique(ivs_cis_df[, .(gene_symbol, id.exposure)])
drop_rows <- lapply(seq_len(nrow(cis_genes_40)), function(k) {
  sym <- cis_genes_40$gene_symbol[k]
  eid <- cis_genes_40$id.exposure[k]
  raw_snps <- unique(ivs_cis_df[id.exposure == eid]$SNP)
  in_out   <- intersect(raw_snps, unique(out_all_dt$SNP))
  h_sub    <- dat_harm_cis_dt[id.exposure == eid]
  h_keep   <- h_sub[mr_keep == TRUE]
  reason <- if (length(in_out) == 0) {
    "Absent from BCAC outcome summary statistics (all cis-IVs/proxies missing)"
  } else if (uniqueN(h_keep$SNP) == 0) {
    "All available instruments excluded by harmonisation"
  } else {
    "Retained"
  }
  data.table(
    gene_symbol = sym, id.exposure = eid,
    n_cis_ivs_raw = length(raw_snps), raw_snps = paste(raw_snps, collapse = ";"),
    n_snps_in_bcac = length(in_out),
    n_unique_snps_post_harm = uniqueN(h_keep$SNP),
    retained_in_cis_mr = uniqueN(h_keep$SNP) > 0,
    drop_reason = reason
  )
})
drop_dt <- rbindlist(drop_rows)[order(-retained_in_cis_mr, gene_symbol)]
fwrite(drop_dt, file.path(res_dir, "module3_cis_harmonization_drop_audit.csv"))
cat("Cis harmonization audit: raw genes =", nrow(drop_dt), "| retained =", sum(drop_dt$retained_in_cis_mr),
    "| dropped =", sum(!drop_dt$retained_in_cis_mr), "\n")
print(drop_dt[retained_in_cis_mr == FALSE])

# (b) ER+ (ieu-a-1127) vs ER- (ieu-a-1128) heterogeneity test across 36 genes
p_erpos <- prim_cis_dt[id.outcome == "ieu-a-1127", .(gene_symbol, id.exposure, nsnp_erpos = nsnp, b_erpos = b, se_erpos = se, OR_erpos = OR, pval_erpos = pval, FDR_erpos = FDR)]
p_erneg <- prim_cis_dt[id.outcome == "ieu-a-1128", .(gene_symbol, id.exposure, nsnp_erneg = nsnp, b_erneg = b, se_erneg = se, OR_erneg = OR, pval_erneg = pval, FDR_erneg = FDR)]
het_er  <- merge(p_erpos, p_erneg, by = c("gene_symbol", "id.exposure"), all = TRUE)
# BCAC subtype summary estimates share controls. Their covariance is unavailable.
# Evaluate a range of correlations; zero is a reference calculation, not an independent-cohort test.
retained <- dat_harm_cis_dt[mr_keep == TRUE & id.outcome %in% c("ieu-a-1127", "ieu-a-1128"),
                           .(instruments = paste(sort(unique(SNP)), collapse = ";")), by = .(id.exposure, id.outcome)]
retained <- dcast(retained, id.exposure ~ id.outcome, value.var = "instruments")
het_er <- merge(het_er, retained, by = "id.exposure", all.x = TRUE)
het_er[, same_instruments := `ieu-a-1127` == `ieu-a-1128`]
het_grid <- rbindlist(lapply(c(0, 0.25, 0.5, 0.75), function(rho) {
  d <- copy(het_er)
  d[, `:=`(assumed_estimate_correlation = rho, contrast_beta = b_erpos - b_erneg)]
  d[, contrast_se := sqrt(se_erpos^2 + se_erneg^2 - 2 * rho * se_erpos * se_erneg)]
  d[, P_het_erpos_vs_erneg := fifelse(same_instruments & contrast_se > 0,
                                      2 * pnorm(-abs(contrast_beta / contrast_se)), NA_real_)]
  d[, FDR_het := p.adjust(P_het_erpos_vs_erneg, method = "BH")]
  d
}))
fwrite(het_grid, file.path(res_dir, "module3_cis_er_subtype_heterogeneity.csv"))

# (d) GWAS metadata, cohort overlap & Steiger prevalence sensitivity (1%, 5%, 12%)
ao_dt <- as.data.table(readRDS(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/opengwas_available_outcomes.rds")))

core_gwas_ids <- c(
  "ieu-a-1126", "ieu-a-1127", "ieu-a-1128",
  "finn-b-C3_BREAST", "finn-b-C3_BREAST_EXALLC",
  "ebi-a-GCST90016960", "ebi-a-GCST90016987", "ebi-a-GCST90017039",
  "ebi-a-GCST90016951", "ebi-a-GCST90016972", "ebi-a-GCST90016989",
  "met-d-Acetate", "met-a-462", "met-a-385", "met-c-893", "met-a-616",
  "ebi-a-GCST90011994", "prot-a-749", "ebi-a-GCST004445", "ebi-a-GCST90000480", "prot-c-3447_64_2",
  "eqtl-a-ENSG00000125538", "eqtl-a-ENSG00000126262", "eqtl-a-ENSG00000169429", "eqtl-a-ENSG00000132170"
)
meta_sub <- ao_dt[id %in% core_gwas_ids, .(
  opengwas_id = id, trait, consortium = as.character(consortium),
  author, year, pmid, population, sample_size, ncase, ncontrol, unit, build, note
)]
fwrite(meta_sub, file.path(res_dir, "module3_gwas_study_metadata.csv"))
cat("Saved GWAS study metadata rows:", nrow(meta_sub), "\n")
fwrite(data.table(
  comparison = c("BCAC ER positive vs ER negative", "BCAC overall vs ER subtypes", "FinnGen C3 BREAST endpoints", "Microbiome, metabolite, eQTL and pQTL exposures vs BCAC"),
  overlap = c("Catalog lists the same 105974 controls for each subtype", "Nested outcomes in the same consortium", "Alternative endpoint/control definitions within FinnGen", "Participant-level overlap cannot be quantified from cached study metadata")),
  file.path(res_dir, "module3_gwas_overlap_audit.csv"))

prev_grid <- c(0.01, 0.05, 0.12)
steiger_sets <- list()

run_steiger_grid <- function(hm_df, pairs_dt, domain_lbl) {
  hdt <- as.data.table(hm_df)
  out_rows <- list()
  for (r in seq_len(nrow(pairs_dt))) {
    eid <- pairs_dt$id.exposure[r]
    oid <- pairs_dt$id.outcome[r]
    lbl <- pairs_dt$label[r]
    sub <- copy(hdt[id.exposure == eid & id.outcome == oid & mr_keep == TRUE])
    if (nrow(sub) == 0) next
    ao_o <- ao_dt[id == oid]
    ao_e <- ao_dt[id == eid]
    if (!nrow(ao_o) || !nrow(ao_e) || anyNA(c(ao_o$ncase[1], ao_o$ncontrol[1]))) next
    if (!"samplesize.exposure" %in% names(sub)) sub[, samplesize.exposure := NA_real_]
    sub[is.na(samplesize.exposure), samplesize.exposure := ao_e$sample_size[1]]
    sub[, samplesize.outcome := ao_o$ncase[1] + ao_o$ncontrol[1]]
    initial_n <- nrow(sub)
    sub <- sub[is.finite(eaf.outcome) & eaf.outcome > 0 & eaf.outcome < 1 &
                 samplesize.exposure > 2 & is.finite(beta.exposure) & se.exposure > 0]
    if (!nrow(sub)) next
    sub[, r.exposure := abs(beta.exposure / se.exposure) /
          sqrt((beta.exposure / se.exposure)^2 + samplesize.exposure - 2)]
    for (pv in prev_grid) {
      sub_p <- copy(sub)
      sub_p[, r.outcome := TwoSampleMR::get_r_from_lor(beta.outcome, eaf.outcome,
                            ao_o$ncase[1], ao_o$ncontrol[1], pv)]
      stg <- as.data.table(TwoSampleMR::directionality_test(as.data.frame(sub_p)))
      if (!is.null(stg) && nrow(stg) > 0) {
        out_rows[[length(out_rows) + 1]] <- data.table(
          domain = domain_lbl, pair_label = lbl, id.exposure = eid, id.outcome = oid,
          prevalence_assumed = pv, n_snps = nrow(sub_p), n_snps_before_metadata_filter = initial_n,
          snp_r2_exposure = stg$snp_r2.exposure[1], snp_r2_outcome = stg$snp_r2.outcome[1],
          correct_causal_direction = stg$correct_causal_direction[1], steiger_pval = stg$steiger_pval[1]
        )
      }
    }
  }
  rbindlist(out_rows, fill = TRUE)
}

pos_pairs <- unique(sig_cis_dt[, .(id.exposure, id.outcome, label = paste(gene_symbol, outcome_label, sep = "->"))])
steiger_sets[["cis_eQTL"]] <- run_steiger_grid(dat_harm_cis_dt, pos_pairs, "cis-eQTL -> BCAC")

hm_mib_fp <- file.path(cache_dir, "mibiogen_bcac_harmonised.rds")
mib_sig_fp <- file.path(res_dir, "mibiogen_bcac_mr_significant_p005.csv")
if (file.exists(hm_mib_fp) && file.exists(mib_sig_fp)) {
  mib_s <- fread(mib_sig_fp)[, .(id.exposure, id.outcome, label = paste(gwas_name, outcome_label, sep = "->"))]
  steiger_sets[["Microbiota"]] <- run_steiger_grid(readRDS(hm_mib_fp), unique(mib_s), "MiBioGen -> BCAC")
}
hm_met_fp <- file.path(cache_dir, "met_bcac_harmonised.rds")
met_sig_fp <- file.path(res_dir, "metabolite_bcac_mr_significant_p005.csv")
if (file.exists(hm_met_fp) && file.exists(met_sig_fp)) {
  met_s <- fread(met_sig_fp)[, .(id.exposure, id.outcome, label = paste(trait_name, outcome_label, sep = "->"))]
  steiger_sets[["Metabolite"]] <- run_steiger_grid(readRDS(hm_met_fp), unique(met_s), "Metabolite -> BCAC")
}
steiger_all <- rbindlist(steiger_sets, fill = TRUE)
fwrite(steiger_all, file.path(res_dir, "module3_gwas_metadata_and_steiger_sensitivity.csv"))
cat("Saved Steiger prevalence sensitivity (1%, 5%, 12%) rows:", nrow(steiger_all), "\n")
if (nrow(steiger_all)) print(steiger_all[domain == "cis-eQTL -> BCAC" & pair_label %in% c("IL1B->BRCA_ER_neg","FFAR2->BRCA_ER_neg","CXCL8->BRCA_ER_neg","PPARG->BRCA_ER_pos")])


