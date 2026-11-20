# 20_module4_os_meta_and_nomogram.R
# 模块 4 加餐：
# 1) CXCL8 5 队列 OS Meta 异质性 (I2 = 78.3%) 来源排查 + 按 ER- / 同质平台分层敏感性 Meta
# 2) PPARG 在 TCGA & METABRIC "ER+ 亚组内部" 校正 Age/Stage/Grade 后的 C-index、Time-ROC (3/5/10年) 增益与 Nomogram 校准验证

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(meta)
  library(timeROC)
  library(rms)
})

bulk_dir <- "00_data/bulk_cohorts"
res_dir  <- "02_results"

# ==============================================================================
# Part 1: CXCL8 5-Cohort OS Heterogeneity (I2 = 78.3%) Deep Decomposition
# ==============================================================================
cat("=== [1/2] CXCL8 OS Meta Heterogeneity Investigation ===\n")
surv_all <- fread(file.path(res_dir, "module4_survival_analysis.csv"))
cx_os <- surv_all[gene_symbol == "CXCL8" & stratum == "Overall" & endpoint == "OS" & !is.na(HR) & HR_LCI95 > 0]
cx_os[, TE := log(HR)]
cx_os[, seTE := (log(HR_UCI95) - log(HR_LCI95)) / (2 * 1.96)]
print(cx_os[, .(cohort, n_samples, n_events, HR = round(HR, 3), LCI = round(HR_LCI95, 3), UCI = round(HR_UCI95, 3), cox_pval = signif(cox_pval, 3))])

meta_summary_row <- function(dt, label) {
  dt <- copy(dt)[!is.na(HR) & !is.na(HR_LCI95) & !is.na(HR_UCI95) & HR_LCI95 > 0]
  if (nrow(dt) < 2) return(NULL)
  dt[, TE := log(HR)]
  dt[, seTE := (log(HR_UCI95) - log(HR_LCI95)) / (2 * 1.96)]
  m <- meta::metagen(TE = TE, seTE = seTE, studlab = paste(cohort, stratum, endpoint), sm = "HR", data = dt, common = TRUE, random = TRUE)
  data.table(
    analysis = label, k = nrow(dt), cohorts = paste(dt$cohort, collapse = "; "),
    total_n = sum(dt$n_samples), total_events = sum(dt$n_events),
    FE_HR = round(exp(m$TE.common), 4), FE_LCI = round(exp(m$lower.common), 4), FE_UCI = round(exp(m$upper.common), 4), FE_pval = signif(m$pval.common, 4),
    RE_HR = round(exp(m$TE.random), 4), RE_LCI = round(exp(m$lower.random), 4), RE_UCI = round(exp(m$upper.random), 4), RE_pval = signif(m$pval.random, 4),
    I2_pct = round(m$I2 * 100, 1), Q_pval = signif(m$pval.Q, 4)
  )
}

het_rows <- list(meta_summary_row(cx_os, "All_5_Cohorts_Overall_OS"))
# Leave-one-out cohort sensitivity
for (drop_c in cx_os$cohort) {
  het_rows[[paste0("LOO_Drop_", drop_c)]] <- meta_summary_row(cx_os[cohort != drop_c], paste0("LOO_Excluding_", drop_c))
}
# Stratified subgroups:
# 1) External microarray cohorts (excl. TCGA RNA-seq)
het_rows[["Microarray_4Cohorts"]] <- meta_summary_row(cx_os[cohort != "TCGA_BRCA"], "Microarray_4Cohorts_Overall_OS")
# 2) Pure Affymetrix GPL96/GPL570 GEO cohorts (excl. TCGA & METABRIC)
het_rows[["GEO_3Cohorts_Affy"]] <- meta_summary_row(cx_os[cohort %in% c("GSE42568", "GSE20685", "GSE7390")], "GEO_3Cohorts_Affy_OS")
# 3) Untreated / node-negative or clean homogeneous Affy (GSE20685 + GSE7390)
het_rows[["GEO_GSE20685_GSE7390"]] <- meta_summary_row(cx_os[cohort %in% c("GSE20685", "GSE7390")], "GEO_Large_GSE20685_GSE7390_OS")
# 4) ER_neg OS across cohorts that have ER_neg stratum in module4_survival_analysis.csv
cx_ern <- surv_all[gene_symbol == "CXCL8" & stratum == "ER_neg" & endpoint == "OS"]
print(cx_ern[, .(cohort, n_samples, n_events, HR = round(HR, 3), LCI = round(HR_LCI95, 3), UCI = round(HR_UCI95, 3), cox_pval = signif(cox_pval, 3), km_pval = signif(km_logrank_pval, 3))])
het_rows[["ER_neg_Stratum_OS"]] <- meta_summary_row(cx_ern, "ER_neg_Stratum_TCGA_METABRIC_OS")

# Also check METABRIC ER_neg RFS & TCGA ER_neg PFS & Multivariable OS in METABRIC ER_neg!
mb_cl <- readRDS(file.path(bulk_dir, "metabric_planA_adipose_clin.rds"))
tc_cl <- readRDS(file.path(bulk_dir, "tcga_planA_adipose_clin.rds"))
mb_ern_d <- mb_cl[ER_clean == "Negative" & !is.na(OS_MONTHS) & OS_MONTHS > 0 & !is.na(OS_event) & !is.na(AGE) & !is.na(STAGE)]
f_mb_ern <- coxph(Surv(OS_MONTHS, OS_event) ~ CXCL8 + AGE + STAGE, data = mb_ern_d)
cat("METABRIC ER- Multivar OS (CXCL8 + Age + Stage, n=", nrow(mb_ern_d), "): HR=",
    round(exp(coef(f_mb_ern)["CXCL8"]), 3), "P=", signif(summary(f_mb_ern)$coefficients["CXCL8", "Pr(>|z|)"], 3), "\n")

het_dt <- rbindlist(het_rows)
fwrite(het_dt, file.path(res_dir, "planA_cxcl8_meta_heterogeneity_audit.csv"))
print(het_dt[, .(analysis, k, total_n, FE_HR, FE_pval, RE_HR, RE_pval, I2_pct, Q_pval)])

# ==============================================================================
# Part 2: PPARG ER+ Internal ROC (3/5/10 yr), C-index Gain & Nomogram
# ==============================================================================
cat("\n=== [2/2] PPARG ER+ ROC Gain & Nomogram Calibration ===\n")

eval_erpos_model <- function(df, time_col, event_col, cohort_nm, has_grade = FALSE) {
  cols_need <- c(time_col, event_col, "PPARG", "AGE", "Stage_num", "ADIPOQ", "PLIN1")
  if (has_grade) cols_need <- c(cols_need, "GRADE")
  if ("StromalScore" %in% names(df)) cols_need <- c(cols_need, "StromalScore")
  d <- copy(df)[complete.cases(df[, ..cols_need])]
  d <- d[get(time_col) > 0]
  d[, time_yr := get(time_col) / 12.0]
  d[, ev := get(event_col)]

  m_specs <- list(
    "M1_PPARG_Alone"    = "PPARG",
    "M2_Clinical_Base"  = if (has_grade) "AGE + Stage_num + GRADE" else "AGE + Stage_num",
    "M3_Clin_Plus_PPARG" = if (has_grade) "AGE + Stage_num + GRADE + PPARG" else "AGE + Stage_num + PPARG"
  )
  if ("StromalScore" %in% names(d)) {
    m_specs[["M4_Clin_Adipose_Base"]] <- "AGE + Stage_num + ADIPOQ + PLIN1 + StromalScore"
    m_specs[["M5_Clin_Adipose_Plus_PPARG"]] <- "AGE + Stage_num + ADIPOQ + PLIN1 + StromalScore + PPARG"
  }

  out_rows <- list()
  lp_list <- list()
  for (mn in names(m_specs)) {
    f <- coxph(as.formula(paste0("Surv(time_yr, ev) ~ ", m_specs[[mn]])), data = d, x = TRUE, y = TRUE)
    lp <- predict(f, type = "lp")
    lp_list[[mn]] <- lp
    cs <- summary(f)$concordance
    tr <- tryCatch(
      timeROC::timeROC(T = d$time_yr, delta = d$ev, marker = lp, cause = 1, times = c(3, 5, 10), iid = FALSE),
      error = function(e) NULL
    )
    auc3 <- if (!is.null(tr) && length(tr$AUC) >= 1) tr$AUC[1] else NA_real_
    auc5 <- if (!is.null(tr) && length(tr$AUC) >= 2) tr$AUC[2] else NA_real_
    auc10 <- if (!is.null(tr) && length(tr$AUC) >= 3) tr$AUC[3] else NA_real_
    lr_p <- if ("PPARG" %in% names(coef(f))) summary(f)$coefficients["PPARG", "Pr(>|z|)"] else NA_real_
    out_rows[[mn]] <- data.table(
      cohort = cohort_nm, model = mn, n = nrow(d), events = sum(d$ev),
      pparg_pval = signif(lr_p, 4),
      C_index = round(as.numeric(cs[1]), 4), C_SE = round(as.numeric(cs[2]), 4),
      AUC_3yr = round(as.numeric(auc3), 4), AUC_5yr = round(as.numeric(auc5), 4), AUC_10yr = round(as.numeric(auc10), 4),
      AIC = round(AIC(f), 2)
    )
  }
  # Likelihood ratio test comparing Clinical Base vs Clinical + PPARG
  f_base <- coxph(as.formula(paste0("Surv(time_yr, ev) ~ ", m_specs[["M2_Clinical_Base"]])), data = d)
  f_full <- coxph(as.formula(paste0("Surv(time_yr, ev) ~ ", m_specs[["M3_Clin_Plus_PPARG"]])), data = d)
  lrt <- anova(f_base, f_full)
  cat(cohort_nm, "LRT P (Clinical Base vs Clinical + PPARG):", signif(lrt[["Pr(>|Chi|)"]][2], 4), "\n")
  if ("M4_Clin_Adipose_Base" %in% names(m_specs)) {
    f_ab <- coxph(as.formula(paste0("Surv(time_yr, ev) ~ ", m_specs[["M4_Clin_Adipose_Base"]])), data = d)
    f_af <- coxph(as.formula(paste0("Surv(time_yr, ev) ~ ", m_specs[["M5_Clin_Adipose_Plus_PPARG"]])), data = d)
    lrt2 <- anova(f_ab, f_af)
    cat(cohort_nm, "LRT P (Clin+Adipose vs Clin+Adipose+PPARG):", signif(lrt2[["Pr(>|Chi|)"]][2], 4), "\n")
  }
  rbindlist(out_rows)
}

tc_erp <- tc_cl[ER == "Positive"]
mb_erp <- mb_cl[ER_clean == "Positive"][, Stage_num := STAGE]

roc_tc <- eval_erpos_model(tc_erp, "OS_MONTHS", "OS_event", "TCGA_ERpos", has_grade = FALSE)
roc_mb <- eval_erpos_model(mb_erp, "OS_MONTHS", "OS_event", "METABRIC_ERpos", has_grade = TRUE)
roc_all <- rbind(roc_tc, roc_mb)
fwrite(roc_all, file.path(res_dir, "planA_erpos_roc_cindex_gain.csv"))
print(roc_all)

# Build and calibrate rms::cph Nomogram on TCGA ER+ (OS in years)
tc_nom <- copy(tc_erp)[!is.na(OS_MONTHS) & OS_MONTHS > 0 & !is.na(OS_event) & !is.na(AGE) & !is.na(Stage_num) & !is.na(PPARG) & !is.na(ADIPOQ) & !is.na(PLIN1) & !is.na(StromalScore)]
tc_nom[, time_yr := OS_MONTHS / 12.0]
tc_nom[, Stage_f := factor(STAGE, levels = c("I", "II", "III", "IV"))]
dd_nom <- datadist(tc_nom[, .(time_yr, OS_event, PPARG, AGE, Stage_f, ADIPOQ, PLIN1, StromalScore)])
options(datadist = "dd_nom")
f_cph <- rms::cph(Surv(time_yr, OS_event) ~ AGE + Stage_f + PPARG + ADIPOQ + StromalScore,
                  data = tc_nom, x = TRUE, y = TRUE, surv = TRUE, time.inc = 5)
cal_5yr <- tryCatch(rms::calibrate(f_cph, cmethod = "KM", method = "boot", u = 5, m = 150, B = 200), error = function(e) NULL)
saveRDS(list(data = tc_nom, cph_fit = f_cph, cal_5yr = cal_5yr), file.path(bulk_dir, "planA_tcga_erpos_nomogram.rds"))
if (!is.null(cal_5yr)) {
  cal_df <- as.data.table(unclass(cal_5yr))
  fwrite(cal_df, file.path(res_dir, "planA_tcga_erpos_nomogram_5yr_cal.csv"))
}
cat("Module 4 statistical upgrades complete!\n")
