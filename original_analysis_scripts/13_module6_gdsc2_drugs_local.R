# 13_module6_gdsc2_drugs_local.R
# 本地 DepMap 24Q4 × GDSC2 药敏分析。此模块无需外部配置、凭据或网络。

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

gdsc2_genes <- c("PPARG", "CXCL8", "IL1B", "FFAR2", "HDAC11", "TNF", "IL12A")
gdsc2_key_brca_drugs <- c(
  "Paclitaxel", "Docetaxel", "Doxorubicin", "Epirubicin", "Cisplatin",
  "Tamoxifen", "Fulvestrant", "Palbociclib", "Ribociclib", "Olaparib",
  "Talazoparib", "Lapatinib", "Alpelisib", "Vorinostat", "Entinostat",
  "5-Fluorouracil", "Gemcitabine", "Methotrexate", "Erlotinib", "Trametinib"
)
gdsc2_subtypes <- c("TNBC", "ER_pos", "HER2_pos")

gdsc2_match_local <- function(bulk_dir) {
  model <- fread(file.path(bulk_dir, "DepMap24Q4_Model.csv"), select = c(
    "ModelID", "CellLineName", "StrippedCellLineName", "ModelSubtypeFeatures", "SangerModelID", "COSMICID"
  ))
  expr_path <- file.path(bulk_dir, "DepMap24Q4_BRCA_expr.csv")
  expr_names <- names(fread(expr_path, nrows = 0L))
  symbols <- sub(" \\(.*$", "", expr_names[-1L])
  keep <- which(symbols %in% gdsc2_genes)
  stopifnot(setequal(symbols[keep], gdsc2_genes), !anyDuplicated(symbols[keep]))
  expr <- fread(expr_path, select = c(1L, keep + 1L))
  setnames(expr, c("ModelID", symbols[keep]))
  stopifnot(!anyDuplicated(expr$ModelID), !anyDuplicated(model$ModelID))
  expr_meta <- merge(expr, model, by = "ModelID")
  stopifnot(nrow(expr_meta) == nrow(expr))

  # 沿用 Plan A 标签优先级 TNBC > ER+ > HER2+，双阳性归 ER_pos。
  # 仅有明确 HER2+ 标签时才进入 HER2_pos，缺失或无法识别的标签保留为 NA。
  flags <- trimws(expr_meta$ModelSubtypeFeatures)
  expr_meta[, Subtype_Group := fcase(
    grepl("TNBC", flags, ignore.case = TRUE), "TNBC",
    grepl("\\bER\\s*(\\+|,)", flags, ignore.case = TRUE, perl = TRUE), "ER_pos",
    grepl("\\bHER2\\s*\\+", flags, ignore.case = TRUE, perl = TRUE), "HER2_pos"
  )]
  clean_name <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
  expr_meta[, `:=`(COSMIC_STR = as.character(as.integer(COSMICID)),
                    SANGER_STR = trimws(as.character(SangerModelID)),
                    CELL_CLEAN = clean_name(StrippedCellLineName))]
  gdsc <- as.data.table(read_excel(file.path(bulk_dir, "GDSC2_fitted_dose_response_27Oct23.xlsx")))
  gdsc[, `:=`(COSMIC_STR = as.character(as.integer(COSMIC_ID)),
              SANGER_STR = trimws(as.character(SANGER_MODEL_ID)),
              CELL_CLEAN = clean_name(CELL_LINE_NAME))]

  # 各行按 Sanger、COSMIC、标准化名称的顺序匹配，每株模型只使用最高优先级的命中。
  match_rows <- function(left, right, id_col, label) {
    out <- merge(left, right[!is.na(get(id_col)) & nzchar(get(id_col))],
                 by = id_col, allow.cartesian = TRUE)
    out[, match_method := label]
    out
  }
  m1 <- match_rows(gdsc, expr_meta, "SANGER_STR", "SangerModelID")
  m2 <- match_rows(gdsc, expr_meta[!ModelID %in% m1$ModelID], "COSMIC_STR", "COSMICID")
  m3 <- match_rows(gdsc, expr_meta[!ModelID %in% c(m1$ModelID, m2$ModelID)],
                   "CELL_CLEAN", "CellLineName")
  matched <- rbindlist(list(m1, m2, m3), use.names = TRUE, fill = TRUE)
  # 若同一来源身份映射到多个 DepMap 模型则停止，避免重复计入细胞系。
  id_map <- unique(matched[!is.na(SANGER_MODEL_ID) & nzchar(SANGER_MODEL_ID), .(SANGER_MODEL_ID, ModelID)])
  stopifnot(all(id_map[, .N, by = SANGER_MODEL_ID]$N == 1L), !anyNA(matched$DRUG_ID))

  drug_meta <- unique(matched[, .(DRUG_ID, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME)])
  stopifnot(!anyDuplicated(drug_meta$DRUG_ID))
  match_cols <- c("ModelID", "CellLineName", "ModelSubtypeFeatures", "Subtype_Group", "match_method", gdsc2_genes)
  model_meta <- unique(matched[, ..match_cols])
  stopifnot(!anyDuplicated(model_meta$ModelID))
  gb <- matched[, .(LN_IC50 = mean(LN_IC50[is.finite(LN_IC50)]),
                    AUC = mean(AUC[is.finite(AUC)])), by = .(ModelID, DRUG_ID)]
  gb <- merge(merge(gb, model_meta, by = "ModelID"), drug_meta, by = "DRUG_ID")
  setorder(gb, DRUG_ID, ModelID)
  stopifnot(!anyDuplicated(gb, by = c("ModelID", "DRUG_ID")))
  cat(sprintf("Matched %d cell lines, %d DRUG_ID (%d unique names), %d ModelID × DRUG_ID pairs.\n",
              nrow(model_meta), nrow(drug_meta), uniqueN(drug_meta$DRUG_NAME), nrow(gb)))
  print(model_meta[, .N, by = Subtype_Group][order(Subtype_Group)])
  print(model_meta[, .N, by = match_method][order(match_method)])
  list(matrix = gb, drugs = drug_meta, models = model_meta)
}

gdsc2_eligible_family <- function(inputs, min_n = 15L) {
  # 以药物层面的有限 LN_IC50 样本数界定家族，独立于基因表达和所有 P 值。
  coverage <- inputs$matrix[, .(n_LN_IC50 = sum(is.finite(LN_IC50))), by = DRUG_ID]
  eligible_ids <- coverage[n_LN_IC50 >= min_n, DRUG_ID]
  excluded <- coverage[n_LN_IC50 < min_n]
  cat(sprintf("Eligible family: %d DRUG_ID (>= %d finite LN_IC50); %d excluded.\n",
              length(eligible_ids), min_n, nrow(excluded)))
  if (nrow(excluded)) print(merge(excluded, inputs$drugs[, .(DRUG_ID, DRUG_NAME)], by = "DRUG_ID"))
  inputs$matrix <- inputs$matrix[DRUG_ID %in% eligible_ids]
  inputs$drugs <- inputs$drugs[DRUG_ID %in% eligible_ids]
  inputs
}

gdsc2_spearman <- function(x, y, min_n) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_n || uniqueN(x[ok]) < 2L || uniqueN(y[ok]) < 2L) {
    return(list(rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  list(rho = unname(ct$estimate), p = ct$p.value)
}

gdsc2_raw_screen <- function(inputs, min_n = 15L) {
  raw <- rbindlist(lapply(gdsc2_genes, function(gene) {
    inputs$matrix[, {
      expr <- get(gene)
      ok <- is.finite(expr) & is.finite(LN_IC50)
      x <- expr[ok]
      y <- LN_IC50[ok]
      n <- length(x)
      raw_ic <- gdsc2_spearman(x, y, min_n)
      raw_auc <- gdsc2_spearman(x, AUC[ok], min_n)
      group <- if (n >= min_n) x >= median(x) else logical()
      diff_ic <- NA_real_
      wp <- NA_real_
      if (any(group) && any(!group)) {
        diff_ic <- mean(y[group]) - mean(y[!group])
        wp <- suppressWarnings(wilcox.test(y[group], y[!group]))$p.value
      }
      .(gene_symbol = gene, n_cell_lines = n,
        spearman_rho_IC50 = raw_ic$rho, pval_IC50 = raw_ic$p,
        spearman_rho_AUC = raw_auc$rho, pval_AUC = raw_auc$p,
        high_vs_low_diff_LN_IC50 = diff_ic, wilcox_pval = wp)
    }, by = DRUG_ID]
  }))
  raw <- merge(raw, inputs$drugs, by = "DRUG_ID")
  raw[, FDR_IC50 := p.adjust(pval_IC50, method = "BH", n = .N), by = gene_symbol]
  raw[, sensitivity_direction := fcase(
    spearman_rho_IC50 < 0, "Higher expr -> Lower IC50 (More Sensitive)",
    spearman_rho_IC50 > 0, "Higher expr -> Higher IC50 (More Resistant)",
    spearman_rho_IC50 == 0, "No rank correlation"
  )]
  stopifnot(nrow(raw) == nrow(inputs$drugs) * length(gdsc2_genes),
            !anyDuplicated(raw, by = c("gene_symbol", "DRUG_ID")))
  setcolorder(raw, c("gene_symbol", "DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME",
                     "n_cell_lines", "spearman_rho_IC50", "pval_IC50", "spearman_rho_AUC", "pval_AUC",
                     "high_vs_low_diff_LN_IC50", "wilcox_pval", "FDR_IC50", "sensitivity_direction"))
  setorder(raw, gene_symbol, pval_IC50, DRUG_ID, na.last = TRUE)
  raw
}

run_module6_gdsc2 <- function(project_dir = Sys.getenv("BRCA_PROJECT_ROOT", unset = ".")) {
  inputs <- gdsc2_match_local(file.path(project_dir, "00_data", "bulk_cohorts"))
  res_dir <- file.path(project_dir, "02_results")
  raw <- gdsc2_raw_screen(gdsc2_eligible_family(inputs))
  # 保留现有原始 CSV 精度，后续分层结果使用未舍入值。
  rounded <- c("spearman_rho_IC50", "spearman_rho_AUC", "high_vs_low_diff_LN_IC50")
  raw[, (rounded) := lapply(.SD, round, digits = 4L), .SDcols = rounded]
  fwrite(raw, file.path(res_dir, "module6_gdsc2_all_drugs_correlation.csv"))
  fwrite(raw[tolower(DRUG_NAME) %in% tolower(gdsc2_key_brca_drugs)],
         file.path(res_dir, "module6_gdsc2_key_brca_drugs.csv"))
  fwrite(inputs$models[order(ModelID), .(ModelID, CellLineName, ModelSubtypeFeatures,
                                         Subtype_Group, match_method)],
         file.path(res_dir, "module6_gdsc2_matched_brca_cell_lines.csv"))
  cat(sprintf("Wrote %d raw gene–drug rows.\n", nrow(raw)))
  invisible(raw)
}

# 脚本 15 第 5 节复用函数并自行写入调整结果。独立运行本脚本时生成原始结果。
if (!exists(".gdsc2_functions_only", inherits = FALSE)) run_module6_gdsc2()
