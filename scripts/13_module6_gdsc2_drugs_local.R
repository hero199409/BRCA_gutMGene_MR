# 13_module6_gdsc2_drugs_local.R
# 模块 6 本地药敏分析：基于 DepMap24Q4 乳腺癌细胞系转录组 + GDSC2 实测药物反应 (LN_IC50 & AUC)
# 分析核心基因 (PPARG, CXCL8, IL1B, FFAR2, HDAC11, TNF, IL12A) 与全量 ~290 种抗肿瘤药物及乳腺癌临床关键药物的敏感性关联

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

bulk_dir <- "00_data/bulk_cohorts"
res_dir  <- "02_results"

# 1. 读取 DepMap 乳腺癌细胞系表达与模型信息
mod <- fread(file.path(bulk_dir, "DepMap24Q4_Model.csv"))
brca_expr_dt <- fread(file.path(bulk_dir, "DepMap24Q4_BRCA_expr.csv"))
setnames(brca_expr_dt, 1, "ModelID")

# 清理列名 "GENE (ENTREZ)" -> "GENE"
raw_cols <- names(brca_expr_dt)[-1]
clean_syms <- sub(" \\(.*$", "", raw_cols)

focus_genes <- c("PPARG", "CXCL8", "IL1B", "FFAR2", "HDAC11", "TNF", "IL12A")
keep_idx <- which(clean_syms %in% focus_genes)
cat("Matched focus genes in DepMap24Q4_BRCA_expr:", paste(clean_syms[keep_idx], collapse = ", "), "\n")

expr_focus <- brca_expr_dt[, c(1, keep_idx + 1), with = FALSE]
setnames(expr_focus, c("ModelID", clean_syms[keep_idx]))
rm(brca_expr_dt); gc()

mod_brca <- mod[ModelID %in% expr_focus$ModelID, .(
  ModelID, CellLineName, StrippedCellLineName, OncotreeSubtype, ModelSubtypeFeatures, SangerModelID, COSMICID
)]
mod_brca[, COSMIC_STR := ifelse(!is.na(COSMICID), as.character(as.integer(COSMICID)), NA_character_)]
expr_meta <- merge(expr_focus, mod_brca, by = "ModelID")

# 2. 读取 GDSC2 实测药敏数据 (仅保留乳腺癌细胞系匹配条目)
cat("Reading GDSC2_fitted_dose_response_27Oct23.xlsx...\n")
gdsc2 <- as.data.table(readxl::read_excel(file.path(bulk_dir, "GDSC2_fitted_dose_response_27Oct23.xlsx")))
gdsc2[, COSMIC_STR := ifelse(!is.na(COSMIC_ID), as.character(as.integer(COSMIC_ID)), NA_character_)]
gdsc2[, SANGER_STR := trimws(as.character(SANGER_MODEL_ID))]
gdsc2[, CELL_CLEAN := toupper(gsub("[^A-Za-z0-9]", "", as.character(CELL_LINE_NAME)))]
expr_meta[, STRIP_CLEAN := toupper(gsub("[^A-Za-z0-9]", "", as.character(StrippedCellLineName)))]

# 多重 ID 精确匹配: 1) SangerModelID 2) COSMICID 3) StrippedCellLineName
m1 <- merge(gdsc2, expr_meta[!is.na(SangerModelID) & nzchar(SangerModelID)],
            by.x = "SANGER_STR", by.y = "SangerModelID", allow.cartesian = TRUE)
rem_ids <- setdiff(expr_meta$ModelID, unique(m1$ModelID))
m2 <- merge(gdsc2, expr_meta[ModelID %in% rem_ids & !is.na(COSMIC_STR)],
            by = "COSMIC_STR", allow.cartesian = TRUE)
rem_ids2 <- setdiff(expr_meta$ModelID, unique(c(m1$ModelID, m2$ModelID)))
m3 <- merge(gdsc2, expr_meta[ModelID %in% rem_ids2],
            by.x = "CELL_CLEAN", by.y = "STRIP_CLEAN", allow.cartesian = TRUE)

gdsc_brca <- rbindlist(list(m1, m2, m3), use.names = TRUE, fill = TRUE)
# 去重同细胞系-同药物记录取均值 (LN_IC50, AUC)
gdsc_brca <- gdsc_brca[, .(
  LN_IC50 = mean(LN_IC50, na.rm = TRUE),
  AUC     = mean(AUC, na.rm = TRUE),
  PUTATIVE_TARGET = PUTATIVE_TARGET[1],
  PATHWAY_NAME    = PATHWAY_NAME[1]
), by = .(ModelID, CellLineName, ModelSubtypeFeatures, DRUG_ID, DRUG_NAME, PPARG, CXCL8, IL1B, FFAR2, HDAC11, TNF, IL12A)]

rm(gdsc2, m1, m2, m3); gc()
cat("Matched GDSC2 Breast Cancer cell lines:", length(unique(gdsc_brca$ModelID)),
    "| Unique drugs:", length(unique(gdsc_brca$DRUG_NAME)), "\n")

# 保存匹配的乳腺癌细胞系药敏清单
fwrite(unique(gdsc_brca[, .(ModelID, CellLineName, ModelSubtypeFeatures)]),
       file.path(res_dir, "module6_gdsc2_matched_brca_cell_lines.csv"))

# 3. 对每个基因 x 每种药物跑 Spearman 相关性 + 高低表达组 Wilcoxon LN_IC50 检验
drug_rows <- list()
drugs_u <- unique(gdsc_brca[, .(DRUG_ID, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME)])

for (g in focus_genes) {
  if (!(g %in% names(gdsc_brca))) next
  for (i in seq_len(nrow(drugs_u))) {
    did <- drugs_u$DRUG_ID[i]
    dnm <- drugs_u$DRUG_NAME[i]
    sub <- gdsc_brca[DRUG_ID == did]
    gv  <- as.numeric(sub[[g]])
    ic  <- as.numeric(sub$LN_IC50)
    au  <- as.numeric(sub$AUC)
    ok  <- !is.na(gv) & !is.na(ic)
    if (sum(ok) < 15 || sd(gv[ok]) == 0) next
    
    ct_ic <- suppressWarnings(cor.test(gv[ok], ic[ok], method = "spearman"))
    ct_au <- suppressWarnings(cor.test(gv[ok], au[ok], method = "spearman"))
    
    grp <- ifelse(gv[ok] >= median(gv[ok]), "High", "Low")
    wp <- NA_real_; diff_ic <- NA_real_
    if (length(unique(grp)) == 2) {
      wt <- suppressWarnings(wilcox.test(ic[ok] ~ grp))
      wp <- wt$p.value
      diff_ic <- mean(ic[ok][grp == "High"]) - mean(ic[ok][grp == "Low"])
    }
    drug_rows[[paste(g, did)]] <- data.table(
      gene_symbol = g, DRUG_ID = did, DRUG_NAME = dnm,
      PUTATIVE_TARGET = drugs_u$PUTATIVE_TARGET[i],
      PATHWAY_NAME = drugs_u$PATHWAY_NAME[i],
      n_cell_lines = sum(ok),
      spearman_rho_IC50 = round(unname(ct_ic$estimate), 4),
      pval_IC50 = ct_ic$p.value,
      spearman_rho_AUC = round(unname(ct_au$estimate), 4),
      pval_AUC = ct_au$p.value,
      high_vs_low_diff_LN_IC50 = round(diff_ic, 4),
      wilcox_pval = wp
    )
  }
}

drug_dt <- rbindlist(drug_rows)
drug_dt[, FDR_IC50 := p.adjust(pval_IC50, method = "BH"), by = gene_symbol]
drug_dt[, sensitivity_direction := fifelse(spearman_rho_IC50 < 0, "Higher expr -> Lower IC50 (More Sensitive)",
                                    "Higher expr -> Higher IC50 (More Resistant)")]

fwrite(drug_dt[order(gene_symbol, pval_IC50)], file.path(res_dir, "module6_gdsc2_all_drugs_correlation.csv"))

# 筛选名义显著结果 & 乳腺癌临床常用关键药物
brca_key_drugs <- c("Paclitaxel", "Docetaxel", "Doxorubicin", "Epirubicin", "Cisplatin",
                    "Tamoxifen", "Fulvestrant", "Palbociclib", "Ribociclib", "Olaparib",
                    "Talazoparib", "Lapatinib", "Alpelisib", "Vorinostat", "Entinostat",
                    "5-Fluorouracil", "Gemcitabine", "Methotrexate", "Erlotinib", "Trametinib")

key_dt <- drug_dt[DRUG_NAME %in% brca_key_drugs][order(gene_symbol, pval_IC50)]
fwrite(key_dt, file.path(res_dir, "module6_gdsc2_key_brca_drugs.csv"))

cat("\n=== 模块 6 PPARG & CXCL8 显著关联药物 Top 12 (P < 0.05) ===\n")
print(drug_dt[gene_symbol %in% c("PPARG", "CXCL8") & pval_IC50 < 0.05][order(gene_symbol, pval_IC50)][
  , .SD[1:min(8, .N)], by = gene_symbol][
  , .(gene_symbol, DRUG_NAME, PATHWAY_NAME, n_cell_lines, rho = spearman_rho_IC50,
      pval = signif(pval_IC50, 3), FDR = signif(FDR_IC50, 3), sensitivity_direction)])

cat("\n=== 乳腺癌临床常用药物中 P < 0.05 的组合 ===\n")
print(key_dt[pval_IC50 < 0.05, .(gene_symbol, DRUG_NAME, n_cell_lines, rho = spearman_rho_IC50,
                                  pval = signif(pval_IC50, 3), sensitivity_direction)])
