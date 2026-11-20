# 15_planA_deep_validation_and_hpa.R
# 方案 A（CXCL8 ER- 风险轴 + PPARG ER+ 保护轴）发表级深度验证脚本：
# 1) MR 敏感性（单 SNP Wald、Leave-One-Out、F/PVE、Steiger）与局部 coloc 区域 SNP 提取
# 2) 官方 HPA XML 蛋白免疫组化（IHC 患者计数）、正常/癌对照与 HPA 单细胞 nCPM 解析（并保存代表性 IHC 图）
# 3) TCGA + METABRIC PAM50 五亚型统计学比较（Kruskal-Wallis + 成对 Wilcoxon）
# 4) 5 临床队列（TCGA / METABRIC / GSE20685 / GSE7390 / GSE42568）Meta-HR 合并与异质性 I2 计算
# 5) GDSC2 x DepMap 47 株细胞系亚型内分层（TNBC n=21 / ER+ n=14）及亚型校正偏回归药敏去混杂验证

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(xml2)
  library(meta)
  library(readxl)
})

cache_dir <- "00_data/mr_cache"
bulk_dir  <- "00_data/bulk_cohorts"
res_dir   <- "02_results"
tmp_dir   <- "tmp"
ihc_dir   <- "03_figures/HPA_IHC_Images"
dir.create(ihc_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. MR 敏感性深度审计：单 SNP Wald、F 统计量/PVE、Leave-One-Out、多方法汇总
# ==============================================================================
cat("=== [1/5] MR Sensitivity & Single-SNP Analysis ===\n")

hm_mib  <- as.data.table(readRDS(file.path(cache_dir, "mibiogen_bcac_harmonised.rds")))
hm_met  <- as.data.table(readRDS(file.path(cache_dir, "met_bcac_harmonised.rds")))
hm_cis  <- as.data.table(readRDS(file.path(cache_dir, "eqtlgen_cis_only_bcac_harmonised.rds")))

out_label_map <- c(
  "ieu-a-1126" = "BRCA_Overall",
  "ieu-a-1127" = "BRCA_ER_pos",
  "ieu-a-1128" = "BRCA_ER_neg"
)

# 筛选方案 A 相关与对照展示的暴露-结局组合
key_mib_pairs <- data.table(
  id.exposure = c("ebi-a-GCST90016960","ebi-a-GCST90016960","ebi-a-GCST90016987","ebi-a-GCST90016987",
                  "ebi-a-GCST90017039","ebi-a-GCST90017039","ebi-a-GCST90017039",
                  "ebi-a-GCST90016951","ebi-a-GCST90016951","ebi-a-GCST90016951",
                  "ebi-a-GCST90016989","ebi-a-GCST90016972"),
  id.outcome  = c("ieu-a-1127","ieu-a-1126","ieu-a-1128","ieu-a-1126",
                  "ieu-a-1126","ieu-a-1127","ieu-a-1128",
                  "ieu-a-1127","ieu-a-1126","ieu-a-1128",
                  "ieu-a-1128","ieu-a-1128"),
  domain = "Microbiota"
)
key_met_pairs <- data.table(
  id.exposure = c("met-d-Acetate","met-d-Acetate","met-d-Acetate",
                  "met-a-462","met-a-462","met-a-385","met-a-385","met-c-893","met-a-446"),
  id.outcome  = c("ieu-a-1126","ieu-a-1127","ieu-a-1128",
                  "ieu-a-1126","ieu-a-1127","ieu-a-1126","ieu-a-1128","ieu-a-1128","ieu-a-1126"),
  domain = "Metabolite"
)

cis_pos <- fread(file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))
key_cis_pairs <- unique(cis_pos[, .(id.exposure, id.outcome, domain = "cis_eQTL")])

extract_pair_details <- function(hm_dt, pairs_dt, gene_map = NULL) {
  single_rows <- list()
  loo_rows    <- list()
  sens_rows   <- list()

  for (i in seq_len(nrow(pairs_dt))) {
    eid <- pairs_dt$id.exposure[i]
    oid <- pairs_dt$id.outcome[i]
    dom <- pairs_dt$domain[i]
    sub_df <- as.data.frame(hm_dt[id.exposure == eid & id.outcome == oid & mr_keep == TRUE])
    if (nrow(sub_df) == 0) next

    # Clean exposure display name
    exp_lbl <- if (!is.null(gene_map) && eid %in% names(gene_map)) gene_map[[eid]] else sub(" \\|\\| id:.*$", "", sub_df$exposure[1])
    out_lbl <- out_label_map[[oid]] %||% oid

    # Calculate single-SNP Wald ratio & F-stat manually + via TwoSampleMR
    sub_dt <- as.data.table(sub_df)
    sub_dt[, F_stat := (beta.exposure / se.exposure)^2]
    sub_dt[, R2_approx := {
      maf <- ifelse(!is.na(eaf.exposure) & eaf.exposure > 0 & eaf.exposure < 1, eaf.exposure, NA_real_)
      2 * (beta.exposure^2) * maf * (1 - maf) / (2 * (beta.exposure^2) * maf * (1 - maf) + (se.exposure^2) * 2 * 31684 * maf * (1 - maf))
    }]
    sub_dt[, wald_b := beta.outcome / beta.exposure]
    sub_dt[, wald_se := abs(se.outcome / beta.exposure)]
    sub_dt[, wald_pval := 2 * pnorm(-abs(wald_b / wald_se))]
    sub_dt[, `:=`(
      domain = dom,
      exposure_clean = exp_lbl,
      outcome_label = out_lbl,
      wald_OR = exp(wald_b),
      wald_LCI95 = exp(wald_b - 1.96 * wald_se),
      wald_UCI95 = exp(wald_b + 1.96 * wald_se)
    )]
    keep_cols <- intersect(c("domain","id.exposure","exposure_clean","id.outcome","outcome_label","SNP",
                             "effect_allele.exposure","other_allele.exposure","eaf.exposure",
                             "beta.exposure","se.exposure","pval.exposure","F_stat",
                             "beta.outcome","se.outcome","pval.outcome",
                             "wald_b","wald_se","wald_pval","wald_OR","wald_LCI95","wald_UCI95"), names(sub_dt))
    single_rows[[paste(eid, oid, sep = "_")]] <- sub_dt[, ..keep_cols]

    # Sensitivity summary: Heterogeneity, Egger intercept, Steiger
    het <- tryCatch(TwoSampleMR::mr_heterogeneity(sub_df), error = function(e) NULL)
    plt <- tryCatch(TwoSampleMR::mr_pleiotropy_test(sub_df), error = function(e) NULL)
    stg <- tryCatch(TwoSampleMR::directionality_test(sub_df), error = function(e) NULL)
    mrs <- tryCatch(TwoSampleMR::mr(sub_df, method_list = c("mr_ivw","mr_wald_ratio","mr_weighted_median","mr_egger_regression","mr_weighted_mode")), error = function(e) NULL)

    if (!is.null(mrs) && nrow(mrs) > 0) {
      mr_dt <- as.data.table(mrs)
      mr_dt[, `:=`(
        domain = dom,
        exposure_clean = exp_lbl,
        outcome_label = out_lbl,
        OR = exp(b),
        OR_LCI95 = exp(b - 1.96 * se),
        OR_UCI95 = exp(b + 1.96 * se),
        Q_pval_IVW = if (!is.null(het) && any(het$method == "Inverse variance weighted")) het$Q_pval[het$method == "Inverse variance weighted"][1] else NA_real_,
        egger_intercept = if (!is.null(plt)) plt$egger_intercept[1] else NA_real_,
        egger_intercept_pval = if (!is.null(plt)) plt$pval[1] else NA_real_,
        steiger_correct = if (!is.null(stg)) stg$correct_causal_direction[1] else NA,
        steiger_pval = if (!is.null(stg)) stg$steiger_pval[1] else NA_real_
      )]
      sens_rows[[paste(eid, oid, sep = "_")]] <- mr_dt
    }

    if (nrow(sub_df) >= 3) {
      loo <- tryCatch(TwoSampleMR::mr_leaveoneout(sub_df), error = function(e) NULL)
      if (!is.null(loo) && nrow(loo) > 0) {
        ldt <- as.data.table(loo)
        ldt[, `:=`(
          domain = dom,
          exposure_clean = exp_lbl,
          outcome_label = out_lbl,
          OR = exp(b),
          OR_LCI95 = exp(b - 1.96 * se),
          OR_UCI95 = exp(b + 1.96 * se)
        )]
        loo_rows[[paste(eid, oid, sep = "_")]] <- ldt
      }
    }
  }
  list(
    single = rbindlist(single_rows, fill = TRUE),
    sens   = rbindlist(sens_rows, fill = TRUE),
    loo    = rbindlist(loo_rows, fill = TRUE)
  )
}

gene_map <- setNames(cis_pos$gene_symbol, cis_pos$id.exposure)
res_mib <- extract_pair_details(hm_mib, key_mib_pairs)
res_met <- extract_pair_details(hm_met, key_met_pairs)
res_cis <- extract_pair_details(hm_cis, key_cis_pairs, gene_map = gene_map)

single_all <- rbindlist(list(res_mib$single, res_met$single, res_cis$single), fill = TRUE)
sens_all   <- rbindlist(list(res_mib$sens,   res_met$sens,   res_cis$sens),   fill = TRUE)
loo_all    <- rbindlist(list(res_mib$loo,    res_met$loo,    res_cis$loo),    fill = TRUE)

fwrite(single_all, file.path(res_dir, "planA_mr_single_snp_wald.csv"))
fwrite(sens_all,   file.path(res_dir, "planA_mr_sensitivity_summary.csv"))
fwrite(loo_all,    file.path(res_dir, "planA_mr_leaveoneout.csv"))

cat("Saved Plan A MR single-SNP rows:", nrow(single_all), "| Sensitivity rows:", nrow(sens_all), "| LOO rows:", nrow(loo_all), "\n")
print(single_all[exposure_clean %in% c("CXCL8", "PPARG") & ((exposure_clean == "CXCL8" & outcome_label == "BRCA_ER_neg") | (exposure_clean == "PPARG" & outcome_label == "BRCA_ER_pos")),
                 .(exposure_clean, outcome_label, SNP, beta.exposure = round(beta.exposure, 3), F_stat = round(F_stat, 1), wald_OR = round(wald_OR, 3), wald_pval = signif(wald_pval, 3))])

# ==============================================================================
# 2. Coloc 区域级 SNP 数据整理（用于 LocusCompare 绘图与 lead SNP 审计）
# ==============================================================================
cat("\n=== [2/5] Coloc Regional SNP Extraction ===\n")
coloc_files <- c(
  "IL1B_ER_neg"   = "coloc_eqtl-a-ENSG00000125538_ieu-a-1128_2_112964709.rds",
  "FFAR2_ER_neg"  = "coloc_eqtl-a-ENSG00000126262_ieu-a-1128_19_35431073.rds",
  "CXCL8_ER_neg"  = "coloc_eqtl-a-ENSG00000169429_ieu-a-1128_4_74072039.rds",
  "PPARG_ER_pos"  = "coloc_eqtl-a-ENSG00000132170_ieu-a-1127_3_12004290.rds",
  "CXCR2_Overall" = "coloc_eqtl-a-ENSG00000180871_ieu-a-1126_2_218491005.rds"
)

coloc_plot_list <- list()
coloc_lead_rows <- list()
for (nm in names(coloc_files)) {
  fp <- file.path(cache_dir, coloc_files[[nm]])
  if (!file.exists(fp)) next
  d <- readRDS(fp)
  e_dt <- as.data.table(d$e_reg)[!is.na(rsid) & !is.na(p) & p > 0][!duplicated(rsid)]
  o_dt <- as.data.table(d$o_reg)[!is.na(rsid) & !is.na(p) & p > 0][!duplicated(rsid)]
  m <- merge(e_dt[, .(rsid, chr, pos = as.numeric(position), beta_e = beta, p_e = p)],
             o_dt[, .(rsid, beta_o = beta, p_o = p)], by = "rsid")
  m[, `:=`(
    pair = nm,
    log10p_e = -log10(pmax(p_e, 1e-300)),
    log10p_o = -log10(pmax(p_o, 1e-300))
  )]
  coloc_plot_list[[nm]] <- m
  lead_e <- m[which.min(p_e)]
  lead_o <- m[which.min(p_o)]
  coloc_lead_rows[[nm]] <- data.table(
    pair = nm, n_shared_snps = nrow(m),
    lead_eqtl_snp = lead_e$rsid[1], lead_eqtl_pos = lead_e$pos[1], lead_eqtl_p = lead_e$p_e[1],
    lead_brca_snp = lead_o$rsid[1], lead_brca_pos = lead_o$pos[1], lead_brca_p = lead_o$p_o[1],
    bp_distance = abs(lead_e$pos[1] - lead_o$pos[1])
  )
}
coloc_plot_dt <- rbindlist(coloc_plot_list)
saveRDS(coloc_plot_dt, file.path(cache_dir, "planA_coloc_plot_snps.rds"))
fwrite(rbindlist(coloc_lead_rows), file.path(res_dir, "planA_coloc_lead_snps_summary.csv"))
print(rbindlist(coloc_lead_rows))

# ==============================================================================
# 3. HPA 官方 XML 蛋白 IHC + 组织 RNA + 单细胞 nCPM 解析（并保存 IHC 微观图片）
# ==============================================================================
cat("\n=== [3/5] HPA Official XML IHC, Bulk RNA & Single-Cell Extraction ===\n")
hpa_genes <- c(PPARG = "ENSG00000132170", CXCL8 = "ENSG00000169429", IL1B = "ENSG00000125538", FFAR2 = "ENSG00000126262")
hpa_sum_rows <- list()
hpa_sc_rows  <- list()

for (sym in names(hpa_genes)) {
  eid <- hpa_genes[[sym]]
  xml_fp <- file.path(tmp_dir, paste0(tolower(sym), "_hpa.xml"))
  if (!file.exists(xml_fp)) {
    cmd <- sprintf('curl -s --retry 2 --connect-timeout 10 --max-time 25 -o "%s" "https://www.proteinatlas.org/%s.xml"', xml_fp, eid)
    system(cmd)
  }
  if (!file.exists(xml_fp) || file.info(xml_fp)$size < 1000) next
  doc <- read_xml(xml_fp)

  # RNA in Normal Breast
  rna_cons <- xml_find_first(doc, "//rnaExpression[@assayType='consensusTissue']/data[tissue='Breast']/level[@type='normalizedRNAExpression']")
  rna_hpa  <- xml_find_first(doc, "//rnaExpression[@assayType='tissue']/data[tissue='Breast']/level[@type='normalizedRNAExpression']")
  cons_ntpm <- if (!is.na(rna_cons)) as.numeric(xml_attr(rna_cons, "expRNA")) else NA_real_
  hpa_ntpm  <- if (!is.na(rna_hpa))  as.numeric(xml_attr(rna_hpa,  "expRNA")) else NA_real_

  # Antibody IHC in Breast Cancer & Normal Breast
  ab_nodes <- xml_find_all(doc, "//antibody")
  ab_used <- NA_character_; norm_ihc <- NA_character_; canc_ihc <- NA_character_
  n_high <- 0L; n_med <- 0L; n_low <- 0L; n_nd <- 0L; canc_total <- 0L
  norm_img_url <- NA_character_; canc_img_url <- NA_character_

  for (ab in ab_nodes) {
    ab_id <- xml_attr(ab, "id")
    c_bc <- xml_find_first(ab, ".//tissueExpression[@assayType='cancer']/data[tissue='Breast cancer']")
    n_bc <- xml_find_first(ab, ".//tissueExpression[@assayType='tissue']/data[tissue='Breast']")
    if (!is.na(c_bc) || !is.na(n_bc)) {
      ab_used <- paste(na.omit(c(ab_used, ab_id)), collapse = "; ")
      if (!is.na(n_bc) && is.na(norm_ihc)) {
        tc <- xml_find_all(n_bc, "tissueCell")
        norm_ihc <- paste(xml_text(xml_find_first(tc, "cellType")), xml_text(xml_find_first(tc, "level[@type='staining']")), sep = "=", collapse = "; ")
        ni <- xml_find_first(n_bc, ".//assayImage/image/imageUrl")
        if (!is.na(ni)) norm_img_url <- xml_text(ni)
      }
      if (!is.na(c_bc) && canc_total == 0L) {
        lv <- xml_find_all(c_bc, "tissueCell/level")
        if (length(lv) > 0) {
          cts <- as.integer(xml_attr(lv, "count"))
          nms <- tolower(trimws(xml_text(lv)))
          for (k in seq_along(nms)) {
            if (nms[k] == "high") n_high <- cts[k]
            if (nms[k] == "medium") n_med <- cts[k]
            if (nms[k] == "low") n_low <- cts[k]
            if (nms[k] == "not detected") n_nd <- cts[k]
          }
          canc_total <- n_high + n_med + n_low + n_nd
          canc_ihc <- paste(paste(nms, cts, sep = ":"), collapse = "; ")
        }
        ci <- xml_find_first(c_bc, ".//assayImage/image/imageUrl")
        if (!is.na(ci)) canc_img_url <- xml_text(ci)
      }
    }
  }

  # For PPARG & CXCL8 download 1 normal breast & 1 breast cancer representative IHC micrograph
  if (sym %in% c("PPARG", "CXCL8")) {
    if (!is.na(norm_img_url) && nzchar(norm_img_url)) {
      dst_n <- file.path(ihc_dir, sprintf("HPA_%s_NormalBreast_IHC.jpg", sym))
      if (!file.exists(dst_n)) system(sprintf('curl -s --max-time 15 -o "%s" "%s"', dst_n, norm_img_url))
    }
    if (!is.na(canc_img_url) && nzchar(canc_img_url)) {
      dst_c <- file.path(ihc_dir, sprintf("HPA_%s_BreastCancer_IHC.jpg", sym))
      if (!file.exists(dst_c)) system(sprintf('curl -s --max-time 15 -o "%s" "%s"', dst_c, canc_img_url))
    }
  }

  # Prognostics in BRCA from XML cancer RNA data
  prog_node <- xml_find_first(doc, "//rnaExpression[@assayType='cancer']/data[contains(tissue, 'Breast Invasive Carcinoma (TCGA)')]/survivalAnalysis")
  if (is.na(prog_node)) {
    prog_node <- xml_find_first(doc, "//data[contains(tissue, 'Breast Invasive Carcinoma (TCGA)')]/survivalAnalysis")
  }
  prog_type <- if (!is.na(prog_node)) xml_attr(prog_node, "prognosticType") else NA_character_
  prog_cat  <- if (!is.na(prog_node)) xml_attr(prog_node, "prognostic") else NA_character_
  prog_p    <- if (!is.na(prog_node)) as.numeric(xml_attr(prog_node, "pValue")) else NA_real_

  hpa_sum_rows[[sym]] <- data.table(
    gene_symbol = sym, ensembl_id = eid, antibodies = ab_used,
    normal_breast_consensus_nTPM = cons_ntpm, normal_breast_hpa_nTPM = hpa_ntpm,
    normal_breast_IHC_cells = norm_ihc,
    cancer_BRCA_IHC_summary = canc_ihc,
    IHC_High = n_high, IHC_Medium = n_med, IHC_Low = n_low, IHC_NotDetected = n_nd, IHC_TotalPatients = canc_total,
    HPA_TCGA_prognostic = prog_cat, HPA_TCGA_prog_type = prog_type, HPA_TCGA_best_pval = prog_p,
    normal_breast_IHC_url = norm_img_url, cancer_BRCA_IHC_url = canc_img_url
  )

  # Single-cell cell types in HPA
  sc_nodes <- xml_find_all(doc, "//singleCellTypeExpression")
  if (length(sc_nodes) > 0) {
    sc_dt <- data.table(
      gene_symbol = sym,
      cell_type = xml_attr(sc_nodes, "name"),
      nCPM = as.numeric(xml_attr(sc_nodes, "expRNA"))
    )
    hpa_sc_rows[[sym]] <- sc_dt
  }
}

hpa_sum_dt <- rbindlist(hpa_sum_rows, fill = TRUE)
hpa_sc_dt  <- rbindlist(hpa_sc_rows, fill = TRUE)
fwrite(hpa_sum_dt, file.path(res_dir, "planA_hpa_protein_ihc_and_rna.csv"))
fwrite(hpa_sc_dt,  file.path(res_dir, "planA_hpa_single_cell_ncpm.csv"))
print(hpa_sum_dt[, .(gene_symbol, normal_breast_consensus_nTPM, IHC_High, IHC_Medium, IHC_Low, IHC_NotDetected, IHC_TotalPatients, HPA_TCGA_best_pval)])

# ==============================================================================
# 4. PAM50 分子亚型检验（TCGA & METABRIC）+ 5 队列 Meta-HR 汇总
# ==============================================================================
cat("\n=== [4/5] PAM50 Subtype Tests & 5-Cohort Meta-HR ===\n")
tc <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
mb <- readRDS(file.path(bulk_dir, "metabric_target_cache.rds"))

cl_tc <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]
pam_genes <- c("PPARG", "CXCL8", "IL1B", "FFAR2")

pam_rows <- list()
for (g in pam_genes) {
  # TCGA PAM50
  dt_tc <- data.table(
    subtype = factor(gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), levels = c("Normal","LumA","LumB","Her2","Basal")),
    expr = as.numeric(tc$tumor_log2[g, ])
  )[!is.na(subtype) & !is.na(expr)]
  kw_tc <- kruskal.test(expr ~ subtype, data = dt_tc)$p.value
  w_bl_a_tc <- wilcox.test(dt_tc[subtype == "Basal"]$expr, dt_tc[subtype == "LumA"]$expr)$p.value
  w_bl_b_tc <- wilcox.test(dt_tc[subtype == "Basal"]$expr, dt_tc[subtype == "LumB"]$expr)$p.value
  agg_tc <- dt_tc[, .(cohort = "TCGA_BRCA", gene_symbol = g, n = .N,
                      mean = round(mean(expr), 3), sd = round(sd(expr), 3),
                      median = round(median(expr), 3), iqr = round(IQR(expr), 3),
                      kw_pval = kw_tc, wilcox_Basal_vs_LumA = w_bl_a_tc, wilcox_Basal_vs_LumB = w_bl_b_tc), by = subtype]
  pam_rows[[paste("TCGA", g)]] <- agg_tc

  # METABRIC CLAUDIN_SUBTYPE
  dt_mb <- data.table(
    subtype = factor(as.character(mb$clin$CLAUDIN_SUBTYPE), levels = c("Normal","LumA","LumB","Her2","Basal","claudin-low")),
    expr = as.numeric(mb$expr_target[g, ])
  )[!is.na(subtype) & !is.na(expr)]
  kw_mb <- kruskal.test(expr ~ subtype, data = dt_mb)$p.value
  w_bl_a_mb <- wilcox.test(dt_mb[subtype == "Basal"]$expr, dt_mb[subtype == "LumA"]$expr)$p.value
  w_bl_b_mb <- wilcox.test(dt_mb[subtype == "Basal"]$expr, dt_mb[subtype == "LumB"]$expr)$p.value
  agg_mb <- dt_mb[, .(cohort = "METABRIC", gene_symbol = g, n = .N,
                      mean = round(mean(expr), 3), sd = round(sd(expr), 3),
                      median = round(median(expr), 3), iqr = round(IQR(expr), 3),
                      kw_pval = kw_mb, wilcox_Basal_vs_LumA = w_bl_a_mb, wilcox_Basal_vs_LumB = w_bl_b_mb), by = subtype]
  pam_rows[[paste("MB", g)]] <- agg_mb
}
pam_dt <- rbindlist(pam_rows)
fwrite(pam_dt, file.path(res_dir, "planA_pam50_subtypes_stats.csv"))
cat("PAM50 KW tests (PPARG & CXCL8):\n")
print(unique(pam_dt[gene_symbol %in% c("PPARG", "CXCL8"), .(cohort, gene_symbol, kw_pval = signif(kw_pval, 3), wilcox_Basal_vs_LumA = signif(wilcox_Basal_vs_LumA, 3))]))

rm(tc, mb); gc()

# 5-Cohort Meta-HR using exact Cox logHR and SE = (log(UCI) - log(LCI)) / (2 * 1.96)
surv_all <- fread(file.path(res_dir, "module4_survival_analysis.csv"))
run_meta_hr <- function(dt, label) {
  dt <- copy(dt)[!is.na(HR) & !is.na(HR_LCI95) & !is.na(HR_UCI95) & HR_LCI95 > 0]
  if (nrow(dt) < 2) return(NULL)
  dt[, TE := log(HR)]
  dt[, seTE := (log(HR_UCI95) - log(HR_LCI95)) / (2 * 1.96)]
  m <- meta::metagen(TE = TE, seTE = seTE, studlab = paste(cohort, stratum, endpoint), sm = "HR", data = dt, common = TRUE, random = TRUE)
  data.table(
    analysis = label,
    k_cohorts = nrow(dt),
    total_n = sum(dt$n_samples),
    total_events = sum(dt$n_events),
    FE_HR = round(exp(m$TE.common), 4),
    FE_LCI = round(exp(m$lower.common), 4),
    FE_UCI = round(exp(m$upper.common), 4),
    FE_pval = signif(m$pval.common, 4),
    RE_HR = round(exp(m$TE.random), 4),
    RE_LCI = round(exp(m$lower.random), 4),
    RE_UCI = round(exp(m$upper.random), 4),
    RE_pval = signif(m$pval.random, 4),
    I2_pct = round(m$I2 * 100, 1),
    Q_pval = signif(m$pval.Q, 4)
  )
}

meta_list <- list(
  run_meta_hr(surv_all[gene_symbol == "PPARG" & stratum == "Overall" & endpoint == "OS"], "PPARG_5Cohort_Overall_OS"),
  run_meta_hr(surv_all[gene_symbol == "PPARG" & stratum == "ER_pos"  & endpoint == "OS"], "PPARG_2Cohort_ERpos_OS"),
  run_meta_hr(surv_all[gene_symbol == "CXCL8" & stratum == "Overall" & endpoint == "OS"], "CXCL8_5Cohort_Overall_OS"),
  run_meta_hr(surv_all[gene_symbol == "CXCL8" & stratum == "Overall" & endpoint %in% c("RFS","MFS","PFS")], "CXCL8_5Cohort_Progression_RFS_MFS_PFS"),
  run_meta_hr(surv_all[gene_symbol == "CXCL8" & cohort %in% c("METABRIC","GSE20685","GSE7390","GSE42568") & stratum == "Overall" & endpoint %in% c("RFS","MFS")], "CXCL8_4Ext_RFS_MFS")
)
meta_hr_dt <- rbindlist(meta_list)
fwrite(meta_hr_dt, file.path(res_dir, "planA_multicohort_meta_survival.csv"))
print(meta_hr_dt)

# ==============================================================================
# 5. GDSC2 x DepMap 47 株乳腺癌细胞系亚型内去混杂检验（TNBC / ER+ 分层 + OLS 偏回归）
# ==============================================================================
cat("\n=== [5/5] GDSC2 Subtype-Stratified & Subtype-Adjusted Sensitivity Check ===\n")
mod <- fread(file.path(bulk_dir, "DepMap24Q4_Model.csv"))
brca_expr <- fread(file.path(bulk_dir, "DepMap24Q4_BRCA_expr.csv"))
setnames(brca_expr, 1, "ModelID")
raw_cols <- names(brca_expr)[-1]
clean_syms <- sub(" \\(.*$", "", raw_cols)
keep_idx <- which(clean_syms %in% c("PPARG", "CXCL8", "IL1B", "FFAR2"))
expr_sub <- brca_expr[, c(1, keep_idx + 1), with = FALSE]
setnames(expr_sub, c("ModelID", clean_syms[keep_idx]))
rm(brca_expr); gc()

mod_sub <- mod[ModelID %in% expr_sub$ModelID, .(ModelID, CellLineName, StrippedCellLineName, ModelSubtypeFeatures, SangerModelID, COSMICID)]
mod_sub[, COSMIC_STR := ifelse(!is.na(COSMICID), as.character(as.integer(COSMICID)), NA_character_)]
expr_meta <- merge(expr_sub, mod_sub, by = "ModelID")

# 统一细胞系大类亚型 (TNBC, ER_pos, HER2_pos) — 注意避开 HER2+ 误匹配 ER
expr_meta[, Subtype_Group := fifelse(
  grepl("TNBC", ModelSubtypeFeatures, ignore.case = TRUE), "TNBC",
  fifelse(grepl("ER\\+|ER,", ModelSubtypeFeatures), "ER_pos", "HER2_pos")
)]
cat("All 70 DepMap BRCA cell line groups:\n")
print(table(expr_meta$Subtype_Group))

gdsc2 <- as.data.table(readxl::read_excel(file.path(bulk_dir, "GDSC2_fitted_dose_response_27Oct23.xlsx")))
gdsc2[, COSMIC_STR := ifelse(!is.na(COSMIC_ID), as.character(as.integer(COSMIC_ID)), NA_character_)]
gdsc2[, SANGER_STR := trimws(as.character(SANGER_MODEL_ID))]
gdsc2[, CELL_CLEAN := toupper(gsub("[^A-Za-z0-9]", "", as.character(CELL_LINE_NAME)))]
expr_meta[, STRIP_CLEAN := toupper(gsub("[^A-Za-z0-9]", "", as.character(StrippedCellLineName)))]

m1 <- merge(gdsc2, expr_meta[!is.na(SangerModelID) & nzchar(SangerModelID)], by.x = "SANGER_STR", by.y = "SangerModelID", allow.cartesian = TRUE)
rem1 <- setdiff(expr_meta$ModelID, unique(m1$ModelID))
m2 <- merge(gdsc2, expr_meta[ModelID %in% rem1 & !is.na(COSMIC_STR)], by = "COSMIC_STR", allow.cartesian = TRUE)
rem2 <- setdiff(expr_meta$ModelID, unique(c(m1$ModelID, m2$ModelID)))
m3 <- merge(gdsc2, expr_meta[ModelID %in% rem2], by.x = "CELL_CLEAN", by.y = "STRIP_CLEAN", allow.cartesian = TRUE)

gb <- rbindlist(list(m1, m2, m3), use.names = TRUE, fill = TRUE)
gb <- gb[, .(
  LN_IC50 = mean(LN_IC50, na.rm = TRUE),
  AUC = mean(AUC, na.rm = TRUE),
  PUTATIVE_TARGET = PUTATIVE_TARGET[1],
  PATHWAY_NAME = PATHWAY_NAME[1]
), by = .(ModelID, CellLineName, ModelSubtypeFeatures, Subtype_Group, DRUG_ID, DRUG_NAME, PPARG, CXCL8, IL1B, FFAR2)]
rm(gdsc2, m1, m2, m3); gc()
cat("47 GDSC2-matched BRCA cell lines by subtype:\n")
print(table(unique(gb[, .(ModelID, Subtype_Group)])$Subtype_Group))

# Focus on all P < 0.05 drugs + key clinical BRCA drugs for PPARG and CXCL8
all_d <- fread(file.path(res_dir, "module6_gdsc2_all_drugs_correlation.csv"))
key_d_ids <- unique(c(
  all_d[gene_symbol %in% c("CXCL8", "PPARG") & pval_IC50 < 0.05]$DRUG_ID,
  all_d[grepl("Fulvestrant|Tamoxifen|Alpelisib|Palbociclib|Ribociclib|Olaparib|Talazoparib|Paclitaxel|Docetaxel|Cisplatin|Doxorubicin|Gemcitabine|Dasatinib|Trametinib|PD0325901|ERK_2440|AZD7762|Vorinostat|Pevonedistat|Lapatinib", DRUG_NAME, ignore.case = TRUE)]$DRUG_ID
))

sub_adj_rows <- list()
for (g in c("CXCL8", "PPARG")) {
  sub_group_target <- ifelse(g == "CXCL8", "TNBC", "ER_pos")
  for (did in key_d_ids) {
    dd <- gb[DRUG_ID == did]
    if (nrow(dd) < 15) next
    x_all <- dd[[g]]
    y_all <- dd$LN_IC50
    ok_all <- is.finite(x_all) & is.finite(y_all)
    if (sum(ok_all) < 15) next
    dd <- dd[ok_all]

    ct_all <- suppressWarnings(cor.test(dd[[g]], dd$LN_IC50, method = "spearman"))
    # Subtype-adjusted OLS: LN_IC50 ~ expr + factor(Subtype_Group)
    fit <- tryCatch(lm(as.formula(paste("LN_IC50 ~", g, "+ factor(Subtype_Group)")), data = dd), error = function(e) NULL)
    sm_f <- if (!is.null(fit)) summary(fit)$coefficients else NULL
    p_adj <- if (!is.null(sm_f) && g %in% rownames(sm_f)) sm_f[g, "Pr(>|t|)"] else NA_real_
    b_adj <- if (!is.null(sm_f) && g %in% rownames(sm_f)) sm_f[g, "Estimate"] else NA_real_

    # Within target subtype (TNBC for CXCL8, ER_pos for PPARG)
    dd_st <- dd[Subtype_Group == sub_group_target]
    n_st <- nrow(dd_st)
    rho_st <- NA_real_; p_st <- NA_real_
    if (n_st >= 6 && sd(dd_st[[g]]) > 0) {
      ct_st <- suppressWarnings(cor.test(dd_st[[g]], dd_st$LN_IC50, method = "spearman"))
      rho_st <- as.numeric(ct_st$estimate)
      p_st   <- ct_st$p.value
    }

    sub_adj_rows[[paste(g, did, sep = "_")]] <- data.table(
      gene_symbol = g, DRUG_ID = did, DRUG_NAME = dd$DRUG_NAME[1],
      PUTATIVE_TARGET = dd$PUTATIVE_TARGET[1], PATHWAY_NAME = dd$PATHWAY_NAME[1],
      n_all = nrow(dd),
      all_spearman_rho = as.numeric(ct_all$estimate), all_spearman_pval = ct_all$p.value,
      subtype_adjusted_beta = b_adj, subtype_adjusted_pval = p_adj,
      target_subtype = sub_group_target, n_target_subtype = n_st,
      within_subtype_rho = rho_st, within_subtype_pval = p_st
    )
  }
}

sub_adj_dt <- rbindlist(sub_adj_rows)
sub_adj_dt[, same_direction_after_adj := sign(all_spearman_rho) == sign(subtype_adjusted_beta)]
sub_adj_dt[, same_direction_within_subtype := sign(all_spearman_rho) == sign(within_subtype_rho)]
fwrite(sub_adj_dt, file.path(res_dir, "planA_gdsc2_subtype_adjusted_drugs.csv"))
saveRDS(gb[DRUG_ID %in% key_d_ids], file.path(bulk_dir, "planA_gdsc2_key_cell_drug_matrix.rds"))

cat("\nTop subtype-robust drugs for CXCL8 (TNBC + Subtype-Adjusted):\n")
print(sub_adj_dt[gene_symbol == "CXCL8" & (all_spearman_pval < 0.01 | subtype_adjusted_pval < 0.05)][order(subtype_adjusted_pval)][1:min(12, .N),
      .(DRUG_NAME, all_rho = round(all_spearman_rho, 3), all_p = signif(all_spearman_pval, 3),
        adj_beta = round(subtype_adjusted_beta, 3), adj_p = signif(subtype_adjusted_pval, 3),
        tnbc_n = n_target_subtype, tnbc_rho = round(within_subtype_rho, 3), tnbc_p = signif(within_subtype_pval, 3))])

cat("\nTop subtype-robust drugs for PPARG (ER+ + Subtype-Adjusted):\n")
print(sub_adj_dt[gene_symbol == "PPARG" & (all_spearman_pval < 0.02 | subtype_adjusted_pval < 0.05)][order(subtype_adjusted_pval)][1:min(12, .N),
      .(DRUG_NAME, all_rho = round(all_spearman_rho, 3), all_p = signif(all_spearman_pval, 3),
        adj_beta = round(subtype_adjusted_beta, 3), adj_p = signif(subtype_adjusted_pval, 3),
        erpos_n = n_target_subtype, erpos_rho = round(within_subtype_rho, 3), erpos_p = signif(within_subtype_pval, 3))])

cat("\nAll Plan A deep validations complete!\n")
