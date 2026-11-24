# Extract adipose & survival markers from TCGA + raw METABRIC cleanly
# 输入为固定的本地缓存，不需要 OpenGWAS 配置或凭据。
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

bulk_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
res_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

# --- Part 1: TCGA ER+ Adipose Marker & Stromal Score Survival Control ---
tc  <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
tme <- readRDS(file.path(bulk_dir, "tcga_tme_scores.rds"))
cl_tc <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]
stopifnot(identical(cl_tc$SAMPLE_ID, colnames(tc$tumor_log2)),
          !anyDuplicated(cl_tc$PATIENT_ID))

ad_genes <- intersect(c("PPARG", "CXCL8", "ADIPOQ", "PLIN1", "FABP4", "LEP"), rownames(tc$tumor_log2))
cat("TCGA matched genes:", paste(ad_genes, collapse = ", "), "\n")

# Normal vs ER+ Tumor correlations of PPARG with Adipocytes (ADIPOQ/PLIN1)
df_tc_norm <- data.table(
  PPARG  = as.numeric(tc$normal_log2["PPARG", ]),
  ADIPOQ = as.numeric(tc$normal_log2["ADIPOQ", ]),
  PLIN1  = as.numeric(tc$normal_log2["PLIN1", ])
)
est_df <- as.data.frame(tme$estimate)
rownames(est_df) <- est_df$cell_type
stroma_v <- as.numeric(est_df["stroma score", colnames(tc$tumor_log2)])
purity_v <- as.numeric(est_df["tumor purity", colnames(tc$tumor_log2)])

df_tc_all <- data.table(
  SAMPLE_ID = cl_tc$SAMPLE_ID,
  ER = cl_tc$ER_Xena,
  AGE = suppressWarnings(as.numeric(cl_tc$AGE)),
  STAGE = gsub("[ABCDabcd]$", "", gsub("^STAGE ", "", toupper(trimws(as.character(cl_tc$AJCC_PATHOLOGIC_TUMOR_STAGE))))),
  OS_MONTHS = suppressWarnings(as.numeric(cl_tc$OS_MONTHS)),
  OS_event  = suppressWarnings(as.integer(cl_tc$OS_event)),
  PFS_MONTHS = suppressWarnings(as.numeric(cl_tc$PFS_MONTHS)),
  PFS_event  = suppressWarnings(as.integer(cl_tc$PFS_event)),
  PPARG  = as.numeric(tc$tumor_log2["PPARG", ]),
  CXCL8  = as.numeric(tc$tumor_log2["CXCL8", ]),
  ADIPOQ = as.numeric(tc$tumor_log2["ADIPOQ", ]),
  PLIN1  = as.numeric(tc$tumor_log2["PLIN1", ]),
  FABP4  = as.numeric(tc$tumor_log2["FABP4", ]),
  StromalScore = stroma_v,
  TumorPurity  = purity_v
)
df_tc_all[!STAGE %in% c("I","II","III","IV"), STAGE := NA_character_]
df_tc_all[, Stage_num := match(STAGE, c("I","II","III","IV"))]
saveRDS(df_tc_all, file.path(bulk_dir, "tcga_planA_adipose_clin.rds"))

valid_endpoint <- function(d, time, event) {
  d[is.finite(get(time)) & get(time) > 0 & get(event) %in% 0:1]
}
complete_numeric <- function(d, variables) {
  d[Reduce(`&`, lapply(d[, ..variables], is.finite))]
}
tc_erp <- valid_endpoint(df_tc_all[ER == "Positive"], "OS_MONTHS", "OS_event")
tc_cc <- complete_numeric(tc_erp, c("PPARG", "AGE", "Stage_num", "ADIPOQ", "PLIN1", "StromalScore"))
tc_pfs <- valid_endpoint(df_tc_all[ER == "Positive"], "PFS_MONTHS", "PFS_event")
tc_pfs_cc <- complete_numeric(tc_pfs, c("PPARG", "AGE", "Stage_num", "ADIPOQ", "PLIN1", "StromalScore"))
cat("TCGA ER+ survival n =", nrow(tc_erp), "| events =", sum(tc_erp$OS_event), "\n")
cat("TCGA Normal r(PPARG, ADIPOQ) =", round(cor(df_tc_norm$PPARG, df_tc_norm$ADIPOQ, method="spearman"), 3),
    "| ER+ Tumor r(PPARG, ADIPOQ) =", round(cor(tc_erp$PPARG, tc_erp$ADIPOQ, method="spearman"), 3),
    "| ER+ Tumor r(PPARG, PLIN1) =", round(cor(tc_erp$PPARG, tc_erp$PLIN1, method="spearman"), 3), "\n")

run_cox_row <- function(dt, cohort_nm, model_nm, formula_str, target_var = "PPARG", sample_set = "ERpos_OS_common_complete_cases") {
  f <- as.formula(formula_str)
  stopifnot(all(complete.cases(dt[, all.vars(f), with = FALSE])))
  fit <- survival::coxph(f, data = dt, na.action = na.fail)
  sm <- summary(fit)
  if (!(target_var %in% rownames(sm$coefficients))) return(NULL)
  data.table(
    cohort = cohort_nm,
    model = model_nm,
    gene = target_var,
    endpoint = if (all.vars(f)[1] == "OS_MONTHS") "OS" else "PFS",
    sample_set = sample_set,
    expression_scale = "per native expression unit",
    n = fit$n,
    events = fit$nevent,
    HR = sm$conf.int[target_var, "exp(coef)"],
    HR_LCI95 = sm$conf.int[target_var, "lower .95"],
    HR_UCI95 = sm$conf.int[target_var, "upper .95"],
    pval = sm$coefficients[target_var, "Pr(>|z|)"]
  )
}

cox_rows <- list(
  run_cox_row(tc_cc, "TCGA_ERpos", "M1_Univariate", "Surv(OS_MONTHS, OS_event) ~ PPARG"),
  run_cox_row(tc_cc, "TCGA_ERpos", "M2_Adj_Adipose(ADIPOQ+PLIN1)", "Surv(OS_MONTHS, OS_event) ~ PPARG + ADIPOQ + PLIN1"),
  run_cox_row(tc_cc, "TCGA_ERpos", "M3_Adj_Adipose+Stromal", "Surv(OS_MONTHS, OS_event) ~ PPARG + ADIPOQ + PLIN1 + StromalScore"),
  run_cox_row(tc_cc, "TCGA_ERpos", "M4_Adj_Age+Stage", "Surv(OS_MONTHS, OS_event) ~ PPARG + AGE + Stage_num"),
  run_cox_row(tc_cc, "TCGA_ERpos", "M5_Adj_Full(Age+Stage+Adipose+Stromal)", "Surv(OS_MONTHS, OS_event) ~ PPARG + AGE + Stage_num + ADIPOQ + PLIN1 + StromalScore"),
  run_cox_row(tc_cc[ADIPOQ < median(ADIPOQ)], "TCGA_ERpos_LowADIPOQ", "M6_LowAdipose_Subgroup(ADIPOQ<med)", "Surv(OS_MONTHS, OS_event) ~ PPARG", sample_set = "common_OS_cases_below_ADIPOQ_median"),
  run_cox_row(tc_cc[PLIN1 < median(PLIN1)],   "TCGA_ERpos_LowPLIN1",  "M7_LowAdipose_Subgroup(PLIN1<med)",  "Surv(OS_MONTHS, OS_event) ~ PPARG", sample_set = "common_OS_cases_below_PLIN1_median"),
  run_cox_row(tc_pfs_cc, "TCGA_ERpos_PFS", "M8_Adj_Full_PFS(Age+Stage+Adipose+Stromal)", "Surv(PFS_MONTHS, PFS_event) ~ PPARG + AGE + Stage_num + ADIPOQ + PLIN1 + StromalScore", sample_set = "ERpos_PFS_complete_cases")
)
rm(tc, tme); gc()

# --- Part 2: Extract ADIPOQ / PLIN1 / FABP4 from METABRIC microarray cleanly ---
mb <- readRDS(file.path(bulk_dir, "metabric_target_cache.rds"))
mb_txt <- file.path(bulk_dir, "brca_metabric/data_mrna_illumina_microarray.txt")
# Read header and scan line-by-line or fread first two cols to find target row indices cleanly
mb_id_col <- fread(mb_txt, select = c(1, 2), data.table = TRUE)
want_g <- c("ADIPOQ", "PLIN1", "FABP4", "LEP", "PPARG", "CXCL8")
hit_idx <- which(mb_id_col[[1]] %in% want_g)
cat("METABRIC raw hits:", paste(mb_id_col[[1]][hit_idx], collapse = ", "), "\n")
# Read only matching lines using grep-like fread on subset or skip/nrows? Fast streaming in R:
con <- file(mb_txt, "r")
hdr <- strsplit(readLines(con, n = 1), "\t")[[1]]
sub_mat_rows <- list()
repeat {
  lns <- readLines(con, n = 2000)
  if (length(lns) == 0) break
  g_nm <- sub("\t.*$", "", lns)
  m_idx <- which(g_nm %in% want_g)
  for (mi in m_idx) {
    vals <- strsplit(lns[mi], "\t")[[1]]
    sub_mat_rows[[vals[1]]] <- suppressWarnings(as.numeric(vals[-c(1, 2)]))
  }
}
close(con)
mb_sids <- hdr[-c(1, 2)]
mb_cl <- copy(mb$clin)
midx <- match(mb_cl$SAMPLE_ID, mb_sids)
expr_idx <- match(mb_cl$SAMPLE_ID, colnames(mb$expr_target))
stopifnot(!anyNA(midx), !anyNA(expr_idx), !anyDuplicated(mb_cl$PATIENT_ID))
for (gn in names(sub_mat_rows)) {
  if (!(gn %in% names(mb_cl))) mb_cl[, (gn) := sub_mat_rows[[gn]][midx]]
}
mb_cl[, PPARG := as.numeric(mb$expr_target["PPARG", expr_idx])]
mb_cl[, CXCL8 := as.numeric(mb$expr_target["CXCL8", expr_idx])]
mb_cl[, AGE := suppressWarnings(as.numeric(AGE_AT_DIAGNOSIS))]
mb_cl[, STAGE := suppressWarnings(as.numeric(TUMOR_STAGE))]
mb_cl[, GRADE := suppressWarnings(as.numeric(GRADE))]
mb_cl[, OS_MONTHS := suppressWarnings(as.numeric(OS_MONTHS))]
mb_cl[, RFS_MONTHS := suppressWarnings(as.numeric(RFS_MONTHS))]

saveRDS(mb_cl, file.path(bulk_dir, "metabric_planA_adipose_clin.rds"))

mb_erp <- valid_endpoint(mb_cl[ER_clean == "Positive"], "OS_MONTHS", "OS_event")
mb_cc <- complete_numeric(mb_erp, c("PPARG", "AGE", "STAGE", "GRADE", "ADIPOQ", "PLIN1"))
cat("METABRIC ER+ survival n =", nrow(mb_erp), "| events =", sum(mb_erp$OS_event),
    "| r(PPARG, ADIPOQ) =", round(cor(mb_erp$PPARG, mb_erp$ADIPOQ, use="complete.obs", method="spearman"), 3),
    "| r(PPARG, PLIN1) =", round(cor(mb_erp$PPARG, mb_erp$PLIN1, use="complete.obs", method="spearman"), 3), "\n")

cox_rows <- c(cox_rows, list(
  run_cox_row(mb_cc, "METABRIC_ERpos", "M1_Univariate", "Surv(OS_MONTHS, OS_event) ~ PPARG"),
  run_cox_row(mb_cc, "METABRIC_ERpos", "M2_Adj_Adipose(ADIPOQ+PLIN1)", "Surv(OS_MONTHS, OS_event) ~ PPARG + ADIPOQ + PLIN1"),
  run_cox_row(mb_cc, "METABRIC_ERpos", "M3_Adj_Age+Stage", "Surv(OS_MONTHS, OS_event) ~ PPARG + AGE + STAGE"),
  run_cox_row(mb_cc, "METABRIC_ERpos", "M4_Adj_Age+Stage+Grade", "Surv(OS_MONTHS, OS_event) ~ PPARG + AGE + STAGE + GRADE"),
  run_cox_row(mb_cc, "METABRIC_ERpos", "M5_Adj_Full(Age+Stage+Grade+Adipose)", "Surv(OS_MONTHS, OS_event) ~ PPARG + AGE + STAGE + GRADE + ADIPOQ + PLIN1"),
  run_cox_row(mb_cc[ADIPOQ < median(ADIPOQ, na.rm=TRUE)], "METABRIC_ERpos_LowADIPOQ", "M6_LowAdipose_Subgroup(ADIPOQ<med)", "Surv(OS_MONTHS, OS_event) ~ PPARG", sample_set = "common_OS_cases_below_ADIPOQ_median"),
  run_cox_row(mb_cc[PLIN1 < median(PLIN1, na.rm=TRUE)],   "METABRIC_ERpos_LowPLIN1",  "M7_LowAdipose_Subgroup(PLIN1<med)",  "Surv(OS_MONTHS, OS_event) ~ PPARG", sample_set = "common_OS_cases_below_PLIN1_median")
))

res_ad <- rbindlist(cox_rows)
fwrite(res_ad, file.path(res_dir, "planA_pparg_adipose_and_multivar_cox.csv"))
print(res_ad[, .(cohort, model, n, events, HR = round(HR, 3), LCI = round(HR_LCI95, 3), UCI = round(HR_UCI95, 3), pval = signif(pval, 3))])
rm(mb, mb_id_col, mb_cl); gc()
