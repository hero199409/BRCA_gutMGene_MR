# 18b_module2_interactions_and_multivar.R
# 补充模块 2 严谨多因素校正与治疗交互作用检验 (P_interaction)
source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

bulk_dir <- "00_data/bulk_cohorts"
res_dir  <- "02_results"

mb <- readRDS(file.path(bulk_dir, "metabric_planA_adipose_clin.rds"))
pc <- readRDS(file.path(bulk_dir, "planA_therapy_plot_cache.rds"))

rows <- list()
add_res <- function(cohort, test_name, term, est_type, est, lci, uci, pval, n, events) {
  rows[[length(rows) + 1]] <<- data.table(
    cohort = cohort, test_name = test_name, term = term, est_type = est_type,
    n = n, events = events, estimate = est, lci95 = lci, uci95 = uci, pval = pval
  )
}

# ------------------------------------------------------------------------------
# 1. METABRIC ER+: PPARG x HORMONE_THERAPY interaction + Adipose adj inside HT=YES
# ------------------------------------------------------------------------------
mb_erp <- mb[ER_clean == "Positive" & toupper(trimws(HORMONE_THERAPY)) %in% c("YES", "NO")]
mb_erp[, HT_bin := ifelse(toupper(trimws(HORMONE_THERAPY)) == "YES", 1L, 0L)]

# Interaction OS
f_os <- mb_erp[!is.na(OS_MONTHS) & OS_MONTHS > 0 & !is.na(OS_event)]
m_os_int <- summary(coxph(Surv(OS_MONTHS, OS_event) ~ PPARG * HT_bin, data = f_os))
add_res("METABRIC_ERpos", "OS_Interaction_PPARG_x_HT", "PPARG:HT_bin", "Ratio_HR",
        m_os_int$conf.int["PPARG:HT_bin", "exp(coef)"], m_os_int$conf.int["PPARG:HT_bin", "lower .95"],
        m_os_int$conf.int["PPARG:HT_bin", "upper .95"], m_os_int$coefficients["PPARG:HT_bin", "Pr(>|z|)"],
        nrow(f_os), sum(f_os$OS_event == 1))

# Interaction RFS
f_rfs <- mb_erp[!is.na(RFS_MONTHS) & RFS_MONTHS > 0 & !is.na(RFS_event)]
m_rfs_int <- summary(coxph(Surv(RFS_MONTHS, RFS_event) ~ PPARG * HT_bin, data = f_rfs))
add_res("METABRIC_ERpos", "RFS_Interaction_PPARG_x_HT", "PPARG:HT_bin", "Ratio_HR",
        m_rfs_int$conf.int["PPARG:HT_bin", "exp(coef)"], m_rfs_int$conf.int["PPARG:HT_bin", "lower .95"],
        m_rfs_int$conf.int["PPARG:HT_bin", "upper .95"], m_rfs_int$coefficients["PPARG:HT_bin", "Pr(>|z|)"],
        nrow(f_rfs), sum(f_rfs$RFS_event == 1))

# HT=YES Adipose-adjusted
ht_yes <- f_os[HT_bin == 1L]
m_ht_ad <- summary(coxph(Surv(OS_MONTHS, OS_event) ~ PPARG + ADIPOQ + PLIN1, data = ht_yes))
add_res("METABRIC_ERpos_HT_YES", "OS_Adj_ADIPOQ_PLIN1", "PPARG", "HR",
        m_ht_ad$conf.int["PPARG", "exp(coef)"], m_ht_ad$conf.int["PPARG", "lower .95"],
        m_ht_ad$conf.int["PPARG", "upper .95"], m_ht_ad$coefficients["PPARG", "Pr(>|z|)"],
        nrow(ht_yes), sum(ht_yes$OS_event == 1))

ht_yes_rfs <- f_rfs[HT_bin == 1L]
m_ht_rfs_ad <- summary(coxph(Surv(RFS_MONTHS, RFS_event) ~ PPARG + ADIPOQ + PLIN1, data = ht_yes_rfs))
add_res("METABRIC_ERpos_HT_YES", "RFS_Adj_ADIPOQ_PLIN1", "PPARG", "HR",
        m_ht_rfs_ad$conf.int["PPARG", "exp(coef)"], m_ht_rfs_ad$conf.int["PPARG", "lower .95"],
        m_ht_rfs_ad$conf.int["PPARG", "upper .95"], m_ht_rfs_ad$coefficients["PPARG", "Pr(>|z|)"],
        nrow(ht_yes_rfs), sum(ht_yes_rfs$RFS_event == 1))

# ------------------------------------------------------------------------------
# 2. GSE25066 & GSE32646: Multivar pCR & DRFS adjusted for ER_status (+ pCR)
# ------------------------------------------------------------------------------
g25 <- copy(pc$gse25066)
g25_pcr <- g25[!is.na(pCR) & ER %in% c("P", "N")]
g25_pcr[, ER_neg := ifelse(ER == "N", 1L, 0L)]
m25_pcr <- glm(pCR ~ CXCL8 + ER_neg, data = g25_pcr, family = binomial())
c25 <- summary(m25_pcr)$coefficients; ci25 <- exp(confint.default(m25_pcr)["CXCL8", ])
add_res("GSE25066", "pCR_Multivar_ER_Adj", "CXCL8", "OR", exp(c25["CXCL8","Estimate"]), ci25[1], ci25[2], c25["CXCL8","Pr(>|z|)"],
        nrow(g25_pcr), sum(g25_pcr$pCR == 1))

g25_drfs <- g25[!is.na(t_drfs) & t_drfs > 0 & !is.na(e_drfs) & ER %in% c("P", "N") & !is.na(pCR)]
g25_drfs[, ER_neg := ifelse(ER == "N", 1L, 0L)]
m25_dr <- summary(coxph(Surv(t_drfs, e_drfs) ~ CXCL8 + ER_neg + pCR, data = g25_drfs))
add_res("GSE25066", "DRFS_Multivar_ER_and_pCR_Adj", "CXCL8", "HR",
        m25_dr$conf.int["CXCL8","exp(coef)"], m25_dr$conf.int["CXCL8","lower .95"],
        m25_dr$conf.int["CXCL8","upper .95"], m25_dr$coefficients["CXCL8","Pr(>|z|)"],
        nrow(g25_drfs), sum(g25_drfs$e_drfs == 1))

g32 <- copy(pc$gse32646)
g32_pcr <- g32[!is.na(pCR) & ER %in% c("positive", "negative")]
g32_pcr[, ER_neg := ifelse(ER == "negative", 1L, 0L)]
m32_pcr <- glm(pCR ~ CXCL8 + ER_neg, data = g32_pcr, family = binomial())
c32 <- summary(m32_pcr)$coefficients; ci32 <- exp(confint.default(m32_pcr)["CXCL8", ])
add_res("GSE32646", "pCR_Multivar_ER_Adj", "CXCL8", "OR", exp(c32["CXCL8","Estimate"]), ci32[1], ci32[2], c32["CXCL8","Pr(>|z|)"],
        nrow(g32_pcr), sum(g32_pcr$pCR == 1))

# ------------------------------------------------------------------------------
# 3. GSE173839 (I-SPY2 Durvalumab/Olaparib vs Control): Formal Arm Interaction
# ------------------------------------------------------------------------------
g173 <- copy(pc$gse173839)[!is.na(pCR)]
g173[, Arm_Durva := ifelse(Arm == "durvalumab/olaparib", 1L, 0L)]
g173[, HRneg := as.integer(HR_neg)]

m173_c <- glm(pCR ~ CXCL8 * Arm_Durva + HRneg, data = g173, family = binomial())
sc173 <- summary(m173_c)$coefficients; ci_173c <- exp(confint.default(m173_c)["CXCL8:Arm_Durva", ])
add_res("GSE173839_ISPY2", "pCR_Interaction_CXCL8_x_DurvaArm(Adj_HR)", "CXCL8:Arm_Durva", "Ratio_OR",
        exp(sc173["CXCL8:Arm_Durva","Estimate"]), ci_173c[1], ci_173c[2], sc173["CXCL8:Arm_Durva","Pr(>|z|)"],
        nrow(g173), sum(g173$pCR == 1))

m173_p <- glm(pCR ~ PPARG * Arm_Durva + HRneg, data = g173, family = binomial())
sp173 <- summary(m173_p)$coefficients; ci_173p <- exp(confint.default(m173_p)["PPARG:Arm_Durva", ])
add_res("GSE173839_ISPY2", "pCR_Interaction_PPARG_x_DurvaArm(Adj_HR)", "PPARG:Arm_Durva", "Ratio_OR",
        exp(sp173["PPARG:Arm_Durva","Estimate"]), ci_173p[1], ci_173p[2], sp173["PPARG:Arm_Durva","Pr(>|z|)"],
        nrow(g173), sum(g173$pCR == 1))

# Also within Durvalumab arm adjusted for HR status
g173_d <- g173[Arm_Durva == 1L]
md_c <- glm(pCR ~ CXCL8 + HRneg, data = g173_d, family = binomial())
sd_c <- summary(md_c)$coefficients; cid_c <- exp(confint.default(md_c)["CXCL8", ])
add_res("GSE173839_DurvaArm", "pCR_Multivar_HR_Adj", "CXCL8", "OR",
        exp(sd_c["CXCL8","Estimate"]), cid_c[1], cid_c[2], sd_c["CXCL8","Pr(>|z|)"],
        nrow(g173_d), sum(g173_d$pCR == 1))

md_p <- glm(pCR ~ PPARG + HRneg, data = g173_d, family = binomial())
sd_p <- summary(md_p)$coefficients; cid_p <- exp(confint.default(md_p)["PPARG", ])
add_res("GSE173839_DurvaArm", "pCR_Multivar_HR_Adj", "PPARG", "OR",
        exp(sd_p["PPARG","Estimate"]), cid_p[1], cid_p[2], sd_p["PPARG","Pr(>|z|)"],
        nrow(g173_d), sum(g173_d$pCR == 1))

out_dt <- rbindlist(rows)
fwrite(out_dt, file.path(res_dir, "planA_therapy_interactions_and_multivar.csv"))
cat("\n=== MODULE 2 INTERACTIONS & MULTIVARIABLE ADJUSTMENTS ===\n")
print(out_dt[, .(cohort, test_name, term, est_type, n, events,
                 est = round(estimate, 3), lci = round(lci95, 3), uci = round(uci95, 3), pval = signif(pval, 3))])
