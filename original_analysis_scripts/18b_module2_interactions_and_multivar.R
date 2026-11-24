# 18b_module2_interactions_and_multivar.R
# Therapy interactions and multivariable adjustments using local analysis caches.
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

bulk_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
res_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

mb <- readRDS(file.path(bulk_dir, "metabric_planA_adipose_clin.rds"))
pc <- readRDS(file.path(bulk_dir, "planA_therapy_plot_cache.rds"))

rows <- list()
add_res <- function(cohort, test_name, term, est_type, est, lci, uci, pval, n, events,
                    fdr_family, stratum, model_formula,
                    n_control = NA_integer_, events_control = NA_integer_,
                    n_durva = NA_integer_, events_durva = NA_integer_) {
  rows[[length(rows) + 1L]] <<- data.table(
    cohort = cohort, test_name = test_name, term = term, est_type = est_type,
    n = n, events = events, estimate = est, lci95 = lci, uci95 = uci, pval = pval,
    stratum = stratum, model_formula = paste(deparse(model_formula), collapse = " "),
    fdr_family = fdr_family, n_control = n_control, events_control = events_control,
    n_durva = n_durva, events_durva = events_durva
  )
}

add_cox <- function(cohort, test_name, term, stratum, fdr_family, formula, data, est_type = "HR") {
  fit <- summary(coxph(formula, data = data))
  add_res(cohort, test_name, term, est_type,
          fit$conf.int[term, "exp(coef)"], fit$conf.int[term, "lower .95"],
          fit$conf.int[term, "upper .95"], fit$coefficients[term, "Pr(>|z|)"],
          fit$n, fit$nevent, fdr_family, stratum, formula)
}

add_logistic <- function(cohort, test_name, term, stratum, fdr_family, formula, data,
                         est_type = "OR") {
  model_vars <- all.vars(formula)
  model_data <- data[complete.cases(data[, ..model_vars])]
  fit <- glm(formula, data = model_data, family = binomial(), na.action = na.fail)
  if (!fit$converged || fit$boundary || anyNA(coef(fit))) {
    stop("Logistic fit failed for ", cohort, " / ", test_name, " / ", term, call. = FALSE)
  }
  coefs <- summary(fit)$coefficients
  ci <- exp(confint.default(fit)[term, ])
  n_control <- events_control <- n_durva <- events_durva <- NA_integer_
  if ("Arm_Durva" %in% names(model_data)) {
    n_control <- model_data[Arm_Durva == 0L, .N]
    events_control <- model_data[Arm_Durva == 0L, sum(pCR == 1L)]
    n_durva <- model_data[Arm_Durva == 1L, .N]
    events_durva <- model_data[Arm_Durva == 1L, sum(pCR == 1L)]
  }
  add_res(cohort, test_name, term, est_type,
          exp(coefs[term, "Estimate"]), ci[1], ci[2], coefs[term, "Pr(>|z|)"],
          nobs(fit), sum(model.response(model.frame(fit)) == 1L),
          fdr_family, stratum, formula, n_control, events_control, n_durva, events_durva)
}

# 1. METABRIC ER+ hormone therapy interactions and adipose adjustment.
mb_erp <- mb[ER_clean == "Positive" & toupper(trimws(HORMONE_THERAPY)) %in% c("YES", "NO")]
mb_erp[, HT_bin := as.integer(toupper(trimws(HORMONE_THERAPY)) == "YES")]
f_os <- mb_erp[!is.na(OS_MONTHS) & OS_MONTHS > 0 & !is.na(OS_event)]
f_rfs <- mb_erp[!is.na(RFS_MONTHS) & RFS_MONTHS > 0 & !is.na(RFS_event)]
add_cox("METABRIC_ERpos", "OS_Interaction_PPARG_x_HT", "PPARG:HT_bin", "ER_pos",
        "METABRIC_HT_interactions", Surv(OS_MONTHS, OS_event) ~ PPARG * HT_bin, f_os, "Ratio_HR")
add_cox("METABRIC_ERpos", "RFS_Interaction_PPARG_x_HT", "PPARG:HT_bin", "ER_pos",
        "METABRIC_HT_interactions", Surv(RFS_MONTHS, RFS_event) ~ PPARG * HT_bin, f_rfs, "Ratio_HR")
add_cox("METABRIC_ERpos_HT_YES", "OS_Adj_ADIPOQ_PLIN1", "PPARG", "ER_pos",
        "METABRIC_HT_adipose_adjusted", Surv(OS_MONTHS, OS_event) ~ PPARG + ADIPOQ + PLIN1,
        f_os[HT_bin == 1L])
add_cox("METABRIC_ERpos_HT_YES", "RFS_Adj_ADIPOQ_PLIN1", "PPARG", "ER_pos",
        "METABRIC_HT_adipose_adjusted", Surv(RFS_MONTHS, RFS_event) ~ PPARG + ADIPOQ + PLIN1,
        f_rfs[HT_bin == 1L])

# 2. NACT analyses adjusted for ER status and (for DRFS) pCR.
g25 <- copy(pc$gse25066)
g25[, ER_neg := fifelse(ER == "N", 1L, fifelse(ER == "P", 0L, NA_integer_))]
add_logistic("GSE25066", "pCR_Multivar_ER_Adj", "CXCL8", "Overall",
             "NACT_pCR_ER_adjusted", pCR ~ CXCL8 + ER_neg, g25)
add_cox("GSE25066", "DRFS_Multivar_ER_and_pCR_Adj", "CXCL8", "Overall",
        "NACT_DRFS_ER_pCR_adjusted", Surv(t_drfs, e_drfs) ~ CXCL8 + ER_neg + pCR,
        g25[!is.na(t_drfs) & t_drfs > 0 & !is.na(e_drfs)])

g32 <- copy(pc$gse32646)
g32[, ER_neg := fifelse(ER == "negative", 1L, fifelse(ER == "positive", 0L, NA_integer_))]
add_logistic("GSE32646", "pCR_Multivar_ER_Adj", "CXCL8", "Overall",
             "NACT_pCR_ER_adjusted", pCR ~ CXCL8 + ER_neg, g32)

# 3. GSE173839 formal interactions. Control is the arm reference. MP is an
# additive covariate with High 1 as reference, verified by script 18.
required173 <- c("Arm", "pCR", "pCR_evaluable", "HER2_neg", "HR_neg", "MP_code", "MammaPrint", "CXCL8", "PPARG")
if (!all(required173 %in% names(pc$gse173839))) stop("Run script 18 to regenerate the audited I-SPY2 cache.")
g173 <- copy(pc$gse173839)[pCR_evaluable & HER2_neg]
stopifnot(all(g173$pCR %in% c(0L, 1L)), all(g173$Arm %in% c("control", "durvalumab/olaparib")),
          !anyNA(g173$HR_neg), identical(levels(g173$MammaPrint), c("High_1", "UltraHigh_2")),
          identical(as.character(g173$MammaPrint), c("High_1", "UltraHigh_2")[g173$MP_code + 1L]))
g173[, `:=`(Arm_Durva = as.integer(Arm == "durvalumab/olaparib"), HRneg = as.integer(HR_neg))]
g173_hrneg <- g173[HR_neg == TRUE]

add_ispy2_interaction <- function(gene, suffix, stratum, family, adjust, data) {
  formula <- reformulate(c(paste0(gene, " * Arm_Durva"), adjust), response = "pCR")
  test_name <- paste0("pCR_Interaction_", gene, "_x_DurvaArm", suffix)
  add_logistic("GSE173839_ISPY2", test_name, paste0(gene, ":Arm_Durva"), stratum,
               family, formula, data, est_type = "Ratio_OR")
}
for (gene in c("CXCL8", "PPARG")) {
  add_ispy2_interaction(gene, "(Adj_HR)", "HER2neg_All", "ISPY2_interactions_HER2neg",
                        "HRneg", g173)
  add_ispy2_interaction(gene, "(HRneg_Subset)", "TNBC_HRneg", "ISPY2_interactions_HRneg_exploratory",
                        NULL, g173_hrneg)
  add_ispy2_interaction(gene, "(Adj_HR_MP)", "HER2neg_All", "ISPY2_interactions_MP_sensitivity",
                        c("HRneg", "MammaPrint"), g173)
  add_ispy2_interaction(gene, "(HRneg_Subset_Adj_MP)", "TNBC_HRneg", "ISPY2_interactions_MP_sensitivity",
                        "MammaPrint", g173_hrneg)
  add_logistic("GSE173839_DurvaArm", "pCR_Multivar_HR_Adj", gene, "HER2neg_All",
               "ISPY2_Durva_arm_HR_adjusted", reformulate(c(gene, "HRneg"), response = "pCR"),
               g173[Arm_Durva == 1L])
}

out_dt <- rbindlist(rows)
# A family comprises biomarkers, endpoints, or NACT cohorts addressing the
# same question. MP sensitivity is a family across both strata. Pooled I-SPY2
# adjustment covers all interaction tests for comparisons across models.
out_dt[, `:=`(fdr_family_n = sum(!is.na(pval)), fdr_bh = p.adjust(pval, method = "BH")),
       by = fdr_family]
out_dt[cohort == "GSE173839_ISPY2" & est_type == "Ratio_OR",
       `:=`(fdr_bh_all_ispy2_interactions = p.adjust(pval, method = "BH"),
            fdr_n_all_ispy2_interactions = sum(!is.na(pval)))]
fwrite(out_dt, file.path(res_dir, "planA_therapy_interactions_and_multivar.csv"))
cat("\n=== MODULE 2 INTERACTIONS & MULTIVARIABLE ADJUSTMENTS ===\n")
print(out_dt[, .(cohort, test_name, term, est_type, stratum, n, events,
                 est = round(estimate, 3), lci = round(lci95, 3), uci = round(uci95, 3),
                 pval = signif(pval, 3), fdr_bh = signif(fdr_bh, 3),
                 fdr_bh_all_ispy2_interactions = signif(fdr_bh_all_ispy2_interactions, 3))])
