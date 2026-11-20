# 18_module2_therapy_response_cohorts.R
# 升级模块 2: 真实临床治疗队列验证 (取代纯 GDSC2 体外细胞系单药敏局限)
# 1) ER+ 内分泌治疗 (Tamoxifen / Hormone Therapy) 队列:
#    - METABRIC ER+ Hormone Therapy == YES (n~1220+)
#    - GSE9195 (GPL570, n=77, ER+ Tamoxifen monotherapy, RFS & DMFS)
#    - GSE17705 (GPL96, n=298, ER+ Tamoxifen 5-yr, Distant Relapse / DRFS)
# 2) 新辅助化疗 & 免疫治疗 (pCR & DRFS, Overall / ER- / TNBC) 队列:
#    - GSE25066 (GPL96, n=508, Taxane-Anthracycline NACT, pCR & DRFS)
#    - GSE32646 (GPL570, n=115, Paclitaxel-FEC NACT, pCR)
#    - GSE173839 (GPL20078, n=105, I-SPY2 Durvalumab+Olaparib vs Control HER2-, pCR)

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(Biobase)
  library(survival)
})

geo_dir  <- "00_data/geo_cache"
bulk_dir <- "00_data/bulk_cohorts"
res_dir  <- "02_results"

target_genes <- c("PPARG", "CXCL8", "IL1B", "FFAR2", "ADIPOQ", "PLIN1", "ESR1")

# ------------------------------------------------------------------------------
# Helper: Parse GPL platform symbol map and extract target genes
# ------------------------------------------------------------------------------
load_gpl_map <- function(gpl_id) {
  if (gpl_id == "GPL20078") {
    tbl <- fread(file.path(geo_dir, "GPL20078_platform_table.tsv"))
    sym_col <- if ("ORF" %in% names(tbl)) "ORF" else names(tbl)[grepl("sym|gene", names(tbl), ignore.case = TRUE)][1]
    return(list(tbl = tbl, id_col = "ID", sym_col = sym_col))
  }
  rds_p <- file.path(geo_dir, paste0(gpl_id, "_symmap.rds"))
  if (file.exists(rds_p)) return(readRDS(rds_p))
  g <- getGEO(filename = file.path(geo_dir, paste0(gpl_id, ".soft.gz")))
  tbl <- as.data.table(Table(g))
  sym_col <- intersect(c("Gene Symbol", "GENE_SYMBOL", "Symbol"), names(tbl))[1]
  out <- list(tbl = tbl[, .(ID = as.character(ID), SYM = as.character(get(sym_col)))])
  saveRDS(out, rds_p)
  out
}

extract_geo_targets <- function(gse_id, gpl_id) {
  cache_rds <- file.path(geo_dir, paste0(gse_id, "_planA_therapy_parsed.rds"))
  if (file.exists(cache_rds)) return(readRDS(cache_rds))

  gset <- getGEO(filename = file.path(geo_dir, paste0(gse_id, "_series_matrix.txt.gz")), getGPL = FALSE)
  ex <- Biobase::exprs(gset)
  pd <- as.data.table(Biobase::pData(gset))

  # Auto-log2 transform only if expression looks un-logged
  qx <- quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  if ((qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)) {
    cat(sprintf("  [%s] Auto log2(x+1) applied (q99=%.1f)\n", gse_id, qx[5]))
    ex[ex < 0] <- 0
    ex <- log2(ex + 1)
  } else {
    cat(sprintf("  [%s] Already log-scale (range %.2f to %.2f)\n", gse_id, qx[1], qx[6]))
  }

  gm <- load_gpl_map(gpl_id)
  if (gpl_id == "GPL20078") {
    dt_map <- data.table(ID = as.character(gm$tbl[[gm$id_col]]), SYM = as.character(gm$tbl[[gm$sym_col]]))
  } else {
    dt_map <- gm$tbl
  }
  dt_map[, primary_sym := trimws(sapply(strsplit(SYM, "///|;| "), `[`, 1))]
  dt_map <- dt_map[primary_sym %in% target_genes & ID %in% rownames(ex)]

  ex_sub <- ex[dt_map$ID, , drop = FALSE]
  dt_e <- as.data.table(ex_sub)
  dt_e[, gene := dt_map$primary_sym]
  dt_agg <- dt_e[, lapply(.SD, mean, na.rm = TRUE), by = gene]
  mat <- as.matrix(dt_agg[, -"gene", with = FALSE])
  rownames(mat) <- dt_agg$gene
  colnames(mat) <- colnames(ex)

  out <- list(expr = mat, pheno = pd)
  saveRDS(out, cache_rds)
  rm(gset, ex); gc()
  out
}

# Helper: Cox & KM stats
run_surv <- function(cohort, arm, stratum, endpoint, gene, expr, time_v, event_v, extra_adj = NULL) {
  ok <- !is.na(expr) & !is.na(time_v) & !is.na(event_v) & time_v > 0
  if (sum(ok) < 25 || sum(event_v[ok] == 1) < 5) return(NULL)
  e <- as.numeric(expr[ok])
  t <- as.numeric(time_v[ok])
  ev <- as.integer(event_v[ok])
  df <- data.frame(t = t, ev = ev, e = e, ez = as.numeric(scale(e)))
  cx <- summary(coxph(Surv(t, ev) ~ e, data = df))
  cxz <- summary(coxph(Surv(t, ev) ~ ez, data = df))
  grp <- ifelse(e >= median(e), "High", "Low")
  km_p <- if (length(unique(grp)) == 2) pchisq(survdiff(Surv(t, ev) ~ grp)$chisq, df = 1, lower.tail = FALSE) else NA_real_

  data.table(
    domain = "Survival_Therapy", cohort = cohort, treatment_arm = arm, stratum = stratum,
    endpoint = endpoint, gene = gene, n_samples = nrow(df), n_events = sum(ev == 1),
    effect_metric = "HR_per_unit", estimate = cx$conf.int[1, "exp(coef)"],
    lci95 = cx$conf.int[1, "lower .95"], uci95 = cx$conf.int[1, "upper .95"],
    pval = cx$coefficients[1, "Pr(>|z|)"], km_or_wilcox_pval = km_p,
    std_estimate = cxz$conf.int[1, "exp(coef)"], std_lci = cxz$conf.int[1, "lower .95"], std_uci = cxz$conf.int[1, "upper .95"],
    auc = NA_real_, resp_n = NA_integer_, nonresp_n = NA_integer_
  )
}

# Helper: Logistic OR + Wilcoxon + exact MWU AUC for binary pCR response
calc_mwu_auc <- function(pred, y) {
  ok <- !is.na(pred) & !is.na(y)
  pred <- pred[ok]; y <- y[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(pred)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

run_pcr <- function(cohort, arm, stratum, gene, expr, y_bin) {
  ok <- !is.na(expr) & !is.na(y_bin) & (y_bin %in% c(0L, 1L))
  e <- as.numeric(expr[ok])
  y <- as.integer(y_bin[ok])
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 < 4 || n0 < 5) return(NULL)
  ez <- as.numeric(scale(e))
  fit <- glm(y ~ e, family = binomial())
  s_fit <- summary(fit)$coefficients
  ci <- tryCatch(exp(confint.default(fit)["e", ]), error = function(er) c(NA_real_, NA_real_))
  fitz <- glm(y ~ ez, family = binomial())
  ciz <- tryCatch(exp(confint.default(fitz)["ez", ]), error = function(er) c(NA_real_, NA_real_))
  wp <- wilcox.test(e[y == 1], e[y == 0])$p.value
  auc_val <- calc_mwu_auc(e, y)
  if (auc_val < 0.5) auc_val <- 1 - auc_val

  data.table(
    domain = "Neoadjuvant_pCR", cohort = cohort, treatment_arm = arm, stratum = stratum,
    endpoint = "pCR_vs_RD", gene = gene, n_samples = length(y), n_events = n1,
    effect_metric = "OR_per_unit", estimate = exp(s_fit["e", "Estimate"]),
    lci95 = ci[1], uci95 = ci[2],
    pval = s_fit["e", "Pr(>|z|)"], km_or_wilcox_pval = wp,
    std_estimate = exp(coef(fitz)["ez"]), std_lci = ciz[1], std_uci = ciz[2],
    auc = auc_val, resp_n = n1, nonresp_n = n0
  )
}

res_list <- list()
plot_cache <- list()

# ==============================================================================
# PART 1: METABRIC Therapy-Stratified Outcomes (Endocrine & Chemo)
# ==============================================================================
cat("\n=== [1/6] METABRIC Therapy Stratification ===\n")
mb <- readRDS(file.path(bulk_dir, "metabric_planA_adipose_clin.rds"))
cat("METABRIC HORMONE_THERAPY x ER_clean:\n")
print(table(ER = mb$ER_clean, HT = mb$HORMONE_THERAPY, useNA = "ifany"))
cat("METABRIC CHEMOTHERAPY x ER_clean:\n")
print(table(ER = mb$ER_clean, CT = mb$CHEMOTHERAPY, useNA = "ifany"))

for (g in c("PPARG", "CXCL8")) {
  # ER+ Endocrine Treated vs Untreated
  ht_yes <- mb[ER_clean == "Positive" & toupper(trimws(HORMONE_THERAPY)) == "YES"]
  ht_no  <- mb[ER_clean == "Positive" & toupper(trimws(HORMONE_THERAPY)) == "NO"]
  res_list[[length(res_list) + 1]] <- run_surv("METABRIC", "Hormone_Therapy_YES", "ER_pos", "RFS", g, ht_yes[[g]], ht_yes$RFS_MONTHS, ht_yes$RFS_event)
  res_list[[length(res_list) + 1]] <- run_surv("METABRIC", "Hormone_Therapy_YES", "ER_pos", "OS",  g, ht_yes[[g]], ht_yes$OS_MONTHS,  ht_yes$OS_event)
  res_list[[length(res_list) + 1]] <- run_surv("METABRIC", "Hormone_Therapy_NO",  "ER_pos", "RFS", g, ht_no[[g]],  ht_no$RFS_MONTHS,  ht_no$RFS_event)
  res_list[[length(res_list) + 1]] <- run_surv("METABRIC", "Hormone_Therapy_NO",  "ER_pos", "OS",  g, ht_no[[g]],  ht_no$OS_MONTHS,   ht_no$OS_event)

  # ER- Chemo Treated vs Untreated
  ct_yes <- mb[ER_clean == "Negative" & toupper(trimws(CHEMOTHERAPY)) == "YES"]
  ct_no  <- mb[ER_clean == "Negative" & toupper(trimws(CHEMOTHERAPY)) == "NO"]
  res_list[[length(res_list) + 1]] <- run_surv("METABRIC", "Chemotherapy_YES", "ER_neg", "RFS", g, ct_yes[[g]], ct_yes$RFS_MONTHS, ct_yes$RFS_event)
  res_list[[length(res_list) + 1]] <- run_surv("METABRIC", "Chemotherapy_YES", "ER_neg", "OS",  g, ct_yes[[g]], ct_yes$OS_MONTHS,  ct_yes$OS_event)
  res_list[[length(res_list) + 1]] <- run_surv("METABRIC", "Chemotherapy_NO",  "ER_neg", "RFS", g, ct_no[[g]],  ct_no$RFS_MONTHS,  ct_no$RFS_event)
  res_list[[length(res_list) + 1]] <- run_surv("METABRIC", "Chemotherapy_NO",  "ER_neg", "OS",  g, ct_no[[g]],  ct_no$OS_MONTHS,   ct_no$OS_event)
}
# Store METABRIC ER+ HT RFS/OS data for plotting
plot_cache$metabric_ht_erpos <- mb[ER_clean == "Positive" & toupper(trimws(HORMONE_THERAPY)) == "YES",
                                   .(SAMPLE_ID, PPARG, CXCL8, ADIPOQ, PLIN1, RFS_MONTHS, RFS_event, OS_MONTHS, OS_event)]

# ==============================================================================
# PART 2: GSE9195 (n=77, ER+ Tamoxifen Monotherapy, GPL570)
# ==============================================================================
cat("\n=== [2/6] GSE9195 (Tamoxifen ER+) ===\n")
g91 <- extract_geo_targets("GSE9195", "GPL570")
pd91 <- g91$pheno
t_rfs91 <- suppressWarnings(as.numeric(as.character(pd91[["t.rfs:ch1"]]))) / 365.25
e_rfs91 <- suppressWarnings(as.integer(as.character(pd91[["e.rfs:ch1"]])))
t_dmf91 <- suppressWarnings(as.numeric(as.character(pd91[["t.dmfs:ch1"]]))) / 365.25
e_dmf91 <- suppressWarnings(as.integer(as.character(pd91[["e.dmfs:ch1"]])))
for (g in intersect(c("PPARG", "CXCL8"), rownames(g91$expr))) {
  res_list[[length(res_list) + 1]] <- run_surv("GSE9195", "Tamoxifen_5yr", "ER_pos", "RFS",  g, g91$expr[g, ], t_rfs91, e_rfs91)
  res_list[[length(res_list) + 1]] <- run_surv("GSE9195", "Tamoxifen_5yr", "ER_pos", "DMFS", g, g91$expr[g, ], t_dmf91, e_dmf91)
}
plot_cache$gse9195 <- data.table(PPARG = as.numeric(g91$expr["PPARG", ]), CXCL8 = as.numeric(g91$expr["CXCL8", ]),
                                 t_rfs = t_rfs91, e_rfs = e_rfs91, t_dmfs = t_dmf91, e_dmfs = e_dmf91)

# ==============================================================================
# PART 3: GSE17705 (n=298, ER+ Tamoxifen 5-yr, GPL96)
# ==============================================================================
cat("\n=== [3/6] GSE17705 (Tamoxifen ER+, n=298) ===\n")
g177 <- extract_geo_targets("GSE17705", "GPL96")
pd177 <- g177$pheno
t_177 <- suppressWarnings(as.numeric(as.character(pd177[["event time (years):ch1"]])))
rel_col177 <- names(pd177)[grepl("relapse", names(pd177), ignore.case = TRUE)]
cat("  GSE17705 relapse cols:", paste(rel_col177, collapse = ", "), "\n")
raw_rel177 <- as.character(pd177[[rel_col177[1]]])
e_177 <- ifelse(grepl(": ?1$|1=dr", raw_rel177), 1L,
         ifelse(grepl(": ?0$|0 censored", raw_rel177), 0L, suppressWarnings(as.integer(gsub(".*:", "", raw_rel177)))))
cat("  GSE17705 Distant Relapse table:\n"); print(table(e_177, useNA = "ifany"))
for (g in intersect(c("PPARG", "CXCL8"), rownames(g177$expr))) {
  res_list[[length(res_list) + 1]] <- run_surv("GSE17705", "Tamoxifen_5yr", "ER_pos", "DMFS", g, g177$expr[g, ], t_177, e_177)
}
# Also check PPARG adjusted for ADIPOQ + PLIN1 in Tamoxifen cohorts
if (all(c("PPARG", "ADIPOQ", "PLIN1") %in% rownames(g177$expr))) {
  df177 <- data.table(t = t_177, ev = e_177, PPARG = as.numeric(g177$expr["PPARG", ]),
                      ADIPOQ = as.numeric(g177$expr["ADIPOQ", ]), PLIN1 = as.numeric(g177$expr["PLIN1", ]))
  df177 <- df177[!is.na(t) & t > 0 & !is.na(ev)]
  cx177_adj <- summary(coxph(Surv(t, ev) ~ PPARG + ADIPOQ + PLIN1, data = df177))
  cat(sprintf("  GSE17705 PPARG DMFS adj ADIPOQ+PLIN1: HR=%.3f (%.3f-%.3f), P=%.4g\n",
              cx177_adj$conf.int["PPARG","exp(coef)"], cx177_adj$conf.int["PPARG","lower .95"],
              cx177_adj$conf.int["PPARG","upper .95"], cx177_adj$coefficients["PPARG","Pr(>|z|)"]))
  plot_cache$gse17705 <- df177
}

# ==============================================================================
# PART 4: GSE25066 (n=508, Neoadjuvant Taxane-Anthracycline, GPL96)
# ==============================================================================
cat("\n=== [4/6] GSE25066 (Neoadjuvant Taxane-Anthracycline, n=508) ===\n")
g25 <- extract_geo_targets("GSE25066", "GPL96")
pd25 <- g25$pheno
pcr_raw25 <- trimws(as.character(pd25[["pathologic_response_pcr_rd:ch1"]]))
y_pcr25 <- ifelse(pcr_raw25 == "pCR", 1L, ifelse(pcr_raw25 == "RD", 0L, NA_integer_))
er25 <- trimws(as.character(pd25[["er_status_ihc:ch1"]]))
pr_col25 <- names(pd25)[grepl("^pr_status", names(pd25), ignore.case = TRUE)][1]
pr25 <- if (!is.na(pr_col25)) trimws(as.character(pd25[[pr_col25]])) else ifelse(grepl("pr_status_ihc: P", as.character(pd25$characteristics_ch1.4)), "P",
         ifelse(grepl("pr_status_ihc: N", as.character(pd25$characteristics_ch1.4)), "N", NA_character_))
her25 <- trimws(as.character(pd25[["her2_status:ch1"]]))
erbb2_25 <- trimws(as.character(pd25[["erbb2_status:ch1"]]))
tnbc25 <- (er25 == "N") & (is.na(pr25) | pr25 != "P") & ((her25 == "N") | (is.na(her25) & erbb2_25 == "N"))

t_drfs25 <- suppressWarnings(as.numeric(as.character(pd25[["drfs_even_time_years:ch1"]])))
e_drfs25 <- suppressWarnings(as.integer(as.character(pd25[["drfs_1_event_0_censored:ch1"]])))
cat("  GSE25066 pCR breakdown (ER_pos / ER_neg / TNBC):\n")
print(table(ER = er25, pCR = y_pcr25, useNA = "ifany"))
print(table(TNBC = tnbc25, pCR = y_pcr25, useNA = "ifany"))

strata25 <- list(
  Overall = rep(TRUE, nrow(pd25)),
  ER_neg  = !is.na(er25) & er25 == "N",
  ER_pos  = !is.na(er25) & er25 == "P",
  TNBC    = !is.na(tnbc25) & tnbc25
)
for (st in names(strata25)) {
  idx <- strata25[[st]]
  for (g in intersect(c("CXCL8", "PPARG"), rownames(g25$expr))) {
    res_list[[length(res_list) + 1]] <- run_pcr("GSE25066", "Taxane_Anthracycline_NACT", st, g, g25$expr[g, idx], y_pcr25[idx])
    res_list[[length(res_list) + 1]] <- run_surv("GSE25066", "Taxane_Anthracycline_NACT", st, "DRFS", g, g25$expr[g, idx], t_drfs25[idx], e_drfs25[idx])
  }
}
plot_cache$gse25066 <- data.table(
  PPARG = as.numeric(g25$expr["PPARG", ]), CXCL8 = as.numeric(g25$expr["CXCL8", ]),
  ER = er25, TNBC = tnbc25, pCR = y_pcr25, t_drfs = t_drfs25, e_drfs = e_drfs25
)

# ==============================================================================
# PART 5: GSE32646 (n=115, Neoadjuvant Paclitaxel + FEC, GPL570)
# ==============================================================================
cat("\n=== [5/6] GSE32646 (Neoadjuvant PTX+FEC, n=115) ===\n")
g32 <- extract_geo_targets("GSE32646", "GPL570")
pd32 <- g32$pheno
pcr_raw32 <- trimws(as.character(pd32[["pathologic response pcr ncr:ch1"]]))
y_pcr32 <- ifelse(pcr_raw32 == "pCR", 1L, ifelse(pcr_raw32 == "nCR", 0L, NA_integer_))
er32  <- tolower(trimws(as.character(pd32[["er status ihc:ch1"]])))
pr32  <- tolower(trimws(as.character(pd32[["pr status ihc:ch1"]])))
her32 <- tolower(trimws(as.character(pd32[["her2 status fish:ch1"]])))
tnbc32 <- (er32 == "negative") & (pr32 == "negative") & (her32 == "negative")
cat("  GSE32646 pCR breakdown:\n")
print(table(ER = er32, pCR = y_pcr32, useNA = "ifany"))
print(table(TNBC = tnbc32, pCR = y_pcr32, useNA = "ifany"))

strata32 <- list(Overall = rep(TRUE, nrow(pd32)), ER_neg = er32 == "negative", ER_pos = er32 == "positive", TNBC = tnbc32)
for (st in names(strata32)) {
  idx <- strata32[[st]]
  for (g in intersect(c("CXCL8", "PPARG"), rownames(g32$expr))) {
    res_list[[length(res_list) + 1]] <- run_pcr("GSE32646", "Paclitaxel_FEC_NACT", st, g, g32$expr[g, idx], y_pcr32[idx])
  }
}
plot_cache$gse32646 <- data.table(
  PPARG = as.numeric(g32$expr["PPARG", ]), CXCL8 = as.numeric(g32$expr["CXCL8", ]),
  ER = er32, TNBC = tnbc32, pCR = y_pcr32
)

# ==============================================================================
# PART 6: GSE173839 (n=105, I-SPY2 Durvalumab/Olaparib vs Control HER2-, GPL20078)
# ==============================================================================
cat("\n=== [6/6] GSE173839 (I-SPY2 Durvalumab+Olaparib Immunotherapy/PARPi) ===\n")
g173 <- extract_geo_targets("GSE173839", "GPL20078")
pd173 <- g173$pheno
arm173 <- trimws(as.character(pd173[["arm:ch1"]]))
pcr_raw173 <- trimws(as.character(pd173[["pcr:ch1"]]))
y_pcr173 <- ifelse(pcr_raw173 == "1", 1L, ifelse(pcr_raw173 == "0", 0L, NA_integer_))
hr_raw173 <- trimws(as.character(pd173[["characteristics_ch1.2"]]))
if (!any(grepl("hr:", hr_raw173))) {
  hr_col <- names(pd173)[grepl("^hr", names(pd173), ignore.case = TRUE)][1]
  hr_raw173 <- as.character(pd173[[hr_col]])
}
hr_neg173 <- grepl("0$", trimws(hr_raw173))
cat("  GSE173839 Arm x HR- x pCR (valid 0/1):\n")
print(table(Arm = arm173, HR_neg = hr_neg173, pCR = y_pcr173, useNA = "ifany"))

arm_list173 <- list(
  list(arm = "Durvalumab_Olaparib_Arm", idx = arm173 == "durvalumab/olaparib", st = "HER2neg_All"),
  list(arm = "Durvalumab_Olaparib_Arm", idx = (arm173 == "durvalumab/olaparib") & hr_neg173, st = "TNBC_HRneg"),
  list(arm = "Control_Paclitaxel_Arm",  idx = arm173 == "control", st = "HER2neg_All"),
  list(arm = "Control_Paclitaxel_Arm",  idx = (arm173 == "control") & hr_neg173, st = "TNBC_HRneg"),
  list(arm = "All_Arms_Combined",       idx = rep(TRUE, nrow(pd173)), st = "HER2neg_All"),
  list(arm = "All_Arms_Combined",       idx = hr_neg173, st = "TNBC_HRneg")
)
for (al in arm_list173) {
  for (g in intersect(c("CXCL8", "PPARG"), rownames(g173$expr))) {
    res_list[[length(res_list) + 1]] <- run_pcr("GSE173839", al$arm, al$st, g, g173$expr[g, al$idx], y_pcr173[al$idx])
  }
}
plot_cache$gse173839 <- data.table(
  PPARG = as.numeric(g173$expr["PPARG", ]), CXCL8 = as.numeric(g173$expr["CXCL8", ]),
  Arm = arm173, HR_neg = hr_neg173, pCR = y_pcr173
)

# ------------------------------------------------------------------------------
# Save unified table & print summary
# ------------------------------------------------------------------------------
res_dt <- rbindlist(res_list, fill = TRUE)
fwrite(res_dt, file.path(res_dir, "planA_therapy_response_summary.csv"))
saveRDS(plot_cache, file.path(bulk_dir, "planA_therapy_plot_cache.rds"))

cat("\n============================================================\n")
cat("MODULE 2 THERAPY COHORT RESULTS SUMMARY\n")
cat("============================================================\n")
print(res_dt[, .(domain, cohort, treatment_arm, stratum, endpoint, gene, n_samples, n_events,
                 est = round(estimate, 3), std_est = round(std_estimate, 3),
                 pval = signif(pval, 3), nonpara_p = signif(km_or_wilcox_pval, 3), auc = round(auc, 3))])
