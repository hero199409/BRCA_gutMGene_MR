# 10_verify_pos_eqtl_cis_coloc_gtex.R
# 严格基于纯 cis-IV (±1Mb) 的 eQTLGen MR 结果，对全部 cis-MR 阳性基因跑：
# 1) cis-IV 敏感性（异质性/多效性/Steiger 方向检验）
# 2) coloc 贝叶斯共定位 (eQTLGen vs BCAC, ±500kb 窗口)
# 3) FinnGen 独立队列验证
# 4) GTEx v8 Breast_Mammary_Tissue API v2 乳腺组织 eQTL 查证

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
  library(coloc)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(jsonlite)
})

cache_dir <- "00_data/mr_cache"
res_dir   <- "02_results"

ao <- readRDS("00_data/opengwas_available_outcomes.rds")
univ <- fread(file.path(res_dir, "module3_candidate_host_genes_universe.csv"))
ivs_all <- as.data.table(readRDS(file.path(cache_dir, "eqtlgen_all49_exp_ivs.rds")))
all_tested_genes <- unique(as.character(ivs_all$gene_symbol))

# 1. 严格标注全量 41 个有 IV 基因的 cis/trans 状态
gene_coords <- as.data.table(suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_tested_genes,
  columns = c("SYMBOL", "ENSEMBL", "CHR", "CHRLOC", "CHRLOCEND"),
  keytype = "SYMBOL"
)))
gene_coords[, start_bp := abs(CHRLOC)]
gene_coords[, end_bp := abs(CHRLOCEND)]
gene_coords <- gene_coords[!is.na(CHR) & !is.na(start_bp), .(
  gene_chr = as.character(CHR[1]),
  gene_start = min(start_bp),
  gene_end = max(end_bp),
  ensembl = na.omit(ENSEMBL)[1]
), by = .(gene_symbol = SYMBOL)]

ivs_annot <- merge(ivs_all, gene_coords, by = "gene_symbol", all.x = TRUE)
ivs_annot[, chr_clean := sub("^chr", "", as.character(chr.exposure))]
ivs_annot[, pos_num := as.numeric(pos.exposure)]
ivs_annot[, is_cis_1mb := (chr_clean == gene_chr) &
            (pos_num >= (gene_start - 1e6)) &
            (pos_num <= (gene_end + 1e6))]
fwrite(ivs_annot[, .(gene_symbol, id.exposure, SNP, chr_clean, pos_num, gene_chr, gene_start, gene_end, is_cis_1mb,
                     beta.exposure, se.exposure, pval.exposure, F_stat)],
       file.path(res_dir, "module3_all41_genes_iv_cis_audit.csv"))

# 2. 用纯 cis-IV 跑全量分母 MR + 异质性 + 多效性 + Steiger
cis_ivs <- ivs_annot[is_cis_1mb == TRUE]
out_bcac <- readRDS(file.path(cache_dir, "outcome_bcac_eqtlgen_all49.rds"))
dat_cis <- TwoSampleMR::harmonise_data(as.data.frame(cis_ivs), as.data.frame(out_bcac), action = 2)
saveRDS(dat_cis, file.path(cache_dir, "eqtlgen_cis_only_bcac_harmonised.rds"))

mr_cis  <- as.data.table(TwoSampleMR::mr(dat_cis, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median", "mr_egger_regression")))
het_cis <- as.data.table(TwoSampleMR::mr_heterogeneity(dat_cis))
plei_cis <- as.data.table(TwoSampleMR::mr_pleiotropy_test(dat_cis))
steiger_cis <- tryCatch(as.data.table(TwoSampleMR::directionality_test(dat_cis)), error = function(e) data.table())

mr_cis[, OR := exp(b)]
mr_cis[, OR_LCI95 := exp(b - 1.96 * se)]
mr_cis[, OR_UCI95 := exp(b + 1.96 * se)]

het_ivw <- het_cis[method == "Inverse variance weighted", .(id.exposure, id.outcome, Q, Q_pval)]
plei_dt <- plei_cis[, .(id.exposure, id.outcome, egger_intercept, pleiotropy_pval = pval)]
dir_dt  <- mr_cis[, .(all_dir_consistent = (all(b > 0) | all(b < 0))), by = .(id.exposure, id.outcome)]

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
fwrite(prim_cis[order(outcome_label, pval)], file.path(res_dir, "module3_eqtlgen_cis_only_mr_primary_full_denominator.csv"))

sig_cis <- prim_cis[pval < 0.05][order(outcome_label, pval)]
fwrite(sig_cis, file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))

pos_genes <- unique(as.character(sig_cis$gene_symbol))
cat("=== 纯 cis-IV eQTLGen MR 汇总 ===\n")
cat("Total cis-IV genes tested:", length(unique(prim_cis$gene_symbol)),
    "| Pairs tested:", nrow(prim_cis),
    "| Nominal P<0.05 pairs:", nrow(sig_cis),
    "| Unique positive cis-genes:", length(pos_genes), "(", paste(pos_genes, collapse=", "), ")\n")
print(sig_cis[, .(outcome_label, gene_symbol, method, nsnp, OR=round(OR,3), pval=signif(pval,3), FDR=signif(FDR,3),
                  all_dir_consistent, Q_pval=signif(Q_pval,2), pleiotropy_pval=signif(pleiotropy_pval,2),
                  correct_causal_direction, upstream_metabolites, alterations)])

# 3. 仅对纯 cis-MR 阳性基因跑 coloc 共定位 (±500kb 围绕基因顺式 lead SNP)
cat("\n=== Running coloc for cis-MR positive pairs ===\n")
bcac_meta <- ao[id %in% c("ieu-a-1126", "ieu-a-1127", "ieu-a-1128"), .(id, ncase, ncontrol, sample_size)]
coloc_rows <- list()

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
  } else {
    Sys.sleep(0.45)
    reg_data <- tryCatch({
      e_reg <- ieugwasr::associations(variants = region_str, id = eid)
      Sys.sleep(0.45)
      o_reg <- ieugwasr::associations(variants = region_str, id = oid)
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
      maf_e <- ifelse(!is.na(merged$eaf.e) & merged$eaf.e > 0 & merged$eaf.e < 1,
                      pmin(merged$eaf.e, 1 - merged$eaf.e), NA_real_)
      if (all(is.na(maf_e))) {
        maf_e <- ifelse(!is.na(merged$eaf.o) & merged$eaf.o > 0 & merged$eaf.o < 1,
                        pmin(merged$eaf.o, 1 - merged$eaf.o), 0.1)
      }
      maf_e[is.na(maf_e) | maf_e <= 0 | maf_e >= 0.5] <- 0.1
      
      n_e <- as.numeric(ao[id == eid]$sample_size[1])
      if (length(n_e) == 0 || is.na(n_e)) n_e <- 31684
      
      bm <- bcac_meta[id == oid]
      s_prop <- bm$ncase[1] / (bm$ncase[1] + bm$ncontrol[1])
      
      res_coloc <- tryCatch({
        coloc::coloc.abf(
          dataset1 = list(snp = merged$rsid, beta = merged$beta.e, varbeta = merged$se.e^2,
                          N = n_e, MAF = maf_e, type = "quant"),
          dataset2 = list(snp = merged$rsid, beta = merged$beta.o, varbeta = merged$se.o^2,
                          s = s_prop, type = "cc")
        )
      }, error = function(e) NULL)
      
      if (!is.null(res_coloc)) {
        sm <- res_coloc$summary
        coloc_rows[[paste(eid, oid, sep = "_")]] <- data.table(
          gene_symbol = sym, id.exposure = eid, outcome_label = olb,
          region = region_str, n_shared_snps = nrow(merged),
          PP.H0 = round(sm["PP.H0.abf"], 4),
          PP.H1 = round(sm["PP.H1.abf"], 4),
          PP.H2 = round(sm["PP.H2.abf"], 4),
          PP.H3 = round(sm["PP.H3.abf"], 4),
          PP.H4 = round(sm["PP.H4.abf"], 4)
        )
      }
    }
  }
}

coloc_dt <- if (length(coloc_rows) > 0) rbindlist(coloc_rows) else data.table()
fwrite(coloc_dt, file.path(res_dir, "module3_cis_positive_genes_coloc.csv"))
cat("\n=== 纯 cis-MR 阳性基因 coloc 共定位结果 ===\n")
print(coloc_dt)

# 4. FinnGen 独立队列验证纯 cis-MR 阳性基因
cat("\n=== Running FinnGen validation for cis-MR positive genes ===\n")
pos_exp_ivs <- as.data.frame(cis_ivs[gene_symbol %in% pos_genes])
finn_cis_cache <- file.path(cache_dir, "outcome_finngen_cis_pos_eqtl.rds")
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
for (sym in pos_genes) {
  ref_cmd <- sprintf('curl -4 -s --connect-timeout 8 --max-time 12 "https://gtexportal.org/api/v2/reference/gene?geneId=%s&gencodeVersion=v26"', sym)
  ref_txt <- paste(suppressWarnings(system(ref_cmd, intern = TRUE)), collapse = "")
  ref_js  <- tryCatch(jsonlite::fromJSON(ref_txt), error = function(e) NULL)
  gencode_id <- NA_character_
  if (!is.null(ref_js) && length(ref_js$data) > 0 && "gencodeId" %in% names(ref_js$data)) {
    gencode_id <- ref_js$data$gencodeId[1]
  }
  if (!is.na(gencode_id)) {
    eq_cmd <- sprintf('curl -4 -s --connect-timeout 8 --max-time 12 "https://gtexportal.org/api/v2/association/independentEqtl?tissueSiteDetailId=Breast_Mammary_Tissue&gencodeId=%s&datasetId=gtex_v8"', gencode_id)
    eq_txt <- paste(suppressWarnings(system(eq_cmd, intern = TRUE)), collapse = "")
    eq_js  <- tryCatch(jsonlite::fromJSON(eq_txt), error = function(e) NULL)
    if (!is.null(eq_js) && length(eq_js$data) > 0 && is.data.frame(eq_js$data)) {
      d <- as.data.table(eq_js$data)
      d[, gene_symbol := sym]
      gtex_rows[[sym]] <- d
      cat("  GTEx Breast_Mammary_Tissue hit for", sym, "(", gencode_id, "): n_indep_eQTL =", nrow(d),
          "| lead SNP =", d$snpId[1], "| p =", d$pValue[1], "| nes =", round(d$nes[1], 3), "\n")
    } else {
      cat("  GTEx Breast_Mammary_Tissue no independent eQTL for", sym, "(", gencode_id, ")\n")
    }
  } else {
    cat("  GTEx gene lookup failed/empty for", sym, "\n")
  }
}
if (length(gtex_rows) > 0) {
  gtex_dt <- rbindlist(gtex_rows, fill = TRUE)
  fwrite(gtex_dt, file.path(res_dir, "module3_cis_positive_genes_gtex_breast_eqtl.csv"))
} else {
  fwrite(data.table(note = "No GTEx Breast_Mammary_Tissue independent eQTL returned or API unavailable"),
         file.path(res_dir, "module3_cis_positive_genes_gtex_breast_eqtl.csv"))
}
