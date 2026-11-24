# 生存尺度、跨队列 OS 及 ER+ OS 区分度审计。仅用本地缓存和 11c 的尺度审计结果。
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(meta)
  library(timeROC)
})

bulk_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
res_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

scale_file <- file.path(res_dir, "module4_survival_scale_audit.csv")
if (!file.exists(scale_file)) stop("Run 11c_update_module4_complete.R --survival-scale-audit first.")
cohort_results <- fread(scale_file)[sensitivity == "positive_times"]
stopifnot(!anyDuplicated(cohort_results[, .(cohort, stratum, endpoint, gene_symbol)]),
          all(cohort_results$expression_scale == "per within-cohort/stratum endpoint SD"))

meta_summary_row <- function(d, label) {
  stopifnot(nrow(d) >= 1L, !anyDuplicated(d$cohort), uniqueN(d$endpoint) == 1L,
            uniqueN(d$stratum) == 1L, uniqueN(d$gene_symbol) == 1L,
            all(is.finite(d$log_HR) & is.finite(d$SE_log_HR) & d$SE_log_HR > 0))
  common <- c(d$HR, d$HR_LCI95, d$HR_UCI95, d$cox_pval)
  random <- rep(NA_real_, 4L)
  hk <- rep(NA_real_, 3L)
  i2 <- q_p <- tau2 <- NA_real_
  method <- "single_cohort_Cox"
  if (nrow(d) > 1L) {
    m <- metagen(TE = log_HR, seTE = SE_log_HR, studlab = cohort, data = d,
                 sm = "HR", common = TRUE, random = TRUE, method.tau = "REML",
                 method.random.ci = c("classic", "HK"), adhoc.hakn.ci = "se")
    common <- c(exp(c(m$TE.common, m$lower.common, m$upper.common)), m$pval.common)
    random <- c(exp(c(m$TE.random[1], m$lower.random[1], m$upper.random[1])), m$pval.random[1])
    hk <- c(exp(c(m$lower.random[2], m$upper.random[2])), m$pval.random[2])
    i2 <- 100 * m$I2
    q_p <- m$pval.Q
    tau2 <- m$tau2
    method <- "REML_random_effects_Wald; modified_Hartung_Knapp_sensitivity"
  }
  data.table(analysis = label, gene_symbol = d$gene_symbol[1], stratum = d$stratum[1],
             endpoint = d$endpoint[1], method, k = nrow(d), cohorts = paste(d$cohort, collapse = "; "),
             total_n = sum(d$n_samples), total_events = sum(d$n_events),
             FE_HR = if (nrow(d) > 1L) common[1] else NA_real_,
             FE_LCI = if (nrow(d) > 1L) common[2] else NA_real_,
             FE_UCI = if (nrow(d) > 1L) common[3] else NA_real_,
             FE_pval = if (nrow(d) > 1L) common[4] else NA_real_,
             RE_HR = random[1], RE_LCI = random[2], RE_UCI = random[3], RE_pval = random[4],
             RE_HK_LCI = hk[1], RE_HK_UCI = hk[2], RE_HK_pval = hk[3],
             single_HR = if (nrow(d) == 1L) common[1] else NA_real_,
             single_LCI = if (nrow(d) == 1L) common[2] else NA_real_,
             single_UCI = if (nrow(d) == 1L) common[3] else NA_real_,
             single_pval = if (nrow(d) == 1L) common[4] else NA_real_,
             I2_pct = i2, Q_pval = q_p, tau2, expression_scale = d$expression_scale[1])
}

os_rows <- cohort_results[endpoint == "OS"]
stopifnot(setequal(os_rows[stratum == "Overall", unique(cohort)],
                  c("TCGA_BRCA", "METABRIC", "GSE42568", "GSE20685", "GSE7390")))
os_meta <- os_rows[, meta_summary_row(.SD, paste(.BY$gene_symbol, .BY$stratum, "OS", sep = "_")),
                   by = .(gene_symbol, stratum), .SDcols = names(os_rows)]
os_meta <- os_meta[, !duplicated(names(os_meta)), with = FALSE]
fwrite(os_meta, file.path(res_dir, "planA_standardized_os_meta.csv"))

endpoint_meta <- rbindlist(lapply(c("RFS", "MFS", "PFS"), function(ep) {
  d <- cohort_results[endpoint == ep & stratum == "Overall"]
  rbindlist(lapply(unique(d$gene_symbol), function(g) {
    meta_summary_row(d[gene_symbol == g], paste(g, "Overall", ep, sep = "_"))
  }))
}))
fwrite(endpoint_meta, file.path(res_dir, "planA_standardized_non_os_endpoints.csv"))

# 复用历史敏感性分析 ID，标准化 HR 用独立文件，便于区分旧图的原始单位结果。
cx_os <- os_rows[gene_symbol == "CXCL8" & stratum == "Overall"]
het_rows <- list(meta_summary_row(cx_os, "All_5_Cohorts_Overall_OS"))
for (drop_c in cx_os$cohort) {
  het_rows[[length(het_rows) + 1L]] <- meta_summary_row(cx_os[cohort != drop_c], paste0("LOO_Excluding_", drop_c))
}
het_rows <- c(het_rows, list(
  meta_summary_row(cx_os[cohort != "TCGA_BRCA"], "Microarray_4Cohorts_Overall_OS"),
  meta_summary_row(cx_os[cohort %in% c("GSE42568", "GSE20685", "GSE7390")], "GEO_3Cohorts_Affy_OS"),
  meta_summary_row(cx_os[cohort %in% c("GSE20685", "GSE7390")], "GEO_Large_GSE20685_GSE7390_OS"),
  meta_summary_row(os_rows[gene_symbol == "CXCL8" & stratum == "ER_neg"], "ER_neg_Stratum_TCGA_METABRIC_OS")
))
het_dt <- rbindlist(het_rows)
fwrite(het_dt, file.path(res_dir, "planA_cxcl8_standardized_os_heterogeneity.csv"))
print(os_meta[, .(analysis, total_n, total_events, RE_HR, RE_LCI, RE_UCI, RE_pval, I2_pct, Q_pval)])
print(het_dt[, .(analysis, k, RE_HR, RE_LCI, RE_UCI, RE_pval, I2_pct, Q_pval)])
print(endpoint_meta[, .(analysis, k, total_n, total_events, RE_HR, RE_LCI, RE_UCI,
                       RE_pval, single_HR, single_LCI, single_UCI, single_pval, I2_pct)])

bootstrap_delta_c <- function(d, f_base, f_plus, base_c, plus_c, bootstrap_reps, seed) {
  # 每次对患者成对重采样，重新拟合两模型。同一组重采样同时估计表观增量不确定性及乐观偏倚。
  delta <- plus_c - base_c
  estimate_c <- function(fit, d_test) {
    score <- predict(fit, newdata = d_test, type = "lp", reference = "zero")
    concordance(Surv(time_yr, ev) ~ score, data = d_test, reverse = TRUE)$concordance
  }
  set.seed(seed)
  reps <- replicate(bootstrap_reps, {
    train <- d[sample.int(nrow(d), replace = TRUE)]
    tryCatch({
      fb <- coxph(f_base, data = train, na.action = na.fail, model = TRUE)
      fp <- coxph(f_plus, data = train, na.action = na.fail, model = TRUE)
      stopifnot(all(is.finite(coef(fb))), all(is.finite(coef(fp))))
      delta_train <- unname(fp$concordance["concordance"] - fb$concordance["concordance"])
      delta_test <- estimate_c(fp, d) - estimate_c(fb, d)
      c(apparent = delta_train, optimism = delta_train - delta_test)
    }, error = function(e) c(apparent = NA_real_, optimism = NA_real_),
       warning = function(w) c(apparent = NA_real_, optimism = NA_real_))
  })
  ok <- colSums(is.finite(reps)) == nrow(reps)
  if (sum(ok) < 0.95 * bootstrap_reps) stop("Fewer than 95% of the paired bootstrap fits succeeded.")
  interval <- quantile(reps["apparent", ok], c(0.025, 0.975), names = FALSE)
  optimism <- mean(reps["optimism", ok])
  data.table(delta_C = delta, delta_C_LCI95 = interval[1], delta_C_UCI95 = interval[2],
             optimism_delta_C = optimism, optimism_corrected_delta_C = delta - optimism,
             bootstrap_reps, bootstrap_valid = sum(ok), seed,
             interval_method = "paired patient bootstrap percentile CI for apparent delta C")
}

model_rows <- function(df, cohort_nm, clinical, adipose, bootstrap_seed) {
  columns <- c("OS_MONTHS", "OS_event", "PPARG", clinical, adipose)
  stopifnot(!anyDuplicated(df$SAMPLE_ID), all(columns %in% names(df)))
  d <- df[is.finite(OS_MONTHS) & OS_MONTHS > 0 & OS_event %in% 0:1]
  d <- d[Reduce(`&`, lapply(d[, ..columns], is.finite))]
  d[, `:=`(time_yr = OS_MONTHS / 12, ev = OS_event)]
  specs <- list(M1_PPARG_Alone = "PPARG", M2_Clinical_Base = clinical,
                M3_Clin_Plus_PPARG = c(clinical, "PPARG"),
                M4_Clin_Adipose_Base = c(clinical, adipose),
                M5_Clin_Adipose_Plus_PPARG = c(clinical, adipose, "PPARG"))
  forms <- lapply(specs, function(p) reformulate(p, response = "Surv(time_yr, ev)"))
  fits <- lapply(forms, coxph, data = d, x = TRUE, y = TRUE, model = TRUE, na.action = na.fail)
  stopifnot(all(vapply(fits, function(f) f$n == nrow(d) && f$nevent == sum(d$ev), logical(1))))
  times <- c(3, 5, 10)
  disc <- rbindlist(lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    s <- summary(f)
    tr <- timeROC(T = d$time_yr, delta = d$ev, marker = predict(f, type = "lp"),
                  cause = 1, times = times, iid = FALSE)
    stopifnot(all(times %in% tr$times))
    auc <- tr$AUC[match(times, tr$times)]
    data.table(cohort = cohort_nm, endpoint = "OS", model = nm, n = f$n, events = f$nevent,
               covariates = paste(specs[[nm]], collapse = " + "),
               pparg_pval = if ("PPARG" %in% specs[[nm]]) s$coefficients["PPARG", "Pr(>|z|)"] else NA_real_,
               C_index = unname(s$concordance[1]), C_SE = unname(s$concordance[2]),
               AUC_3yr = auc[1], AUC_5yr = auc[2], AUC_10yr = auc[3], AIC = AIC(f),
               cases_before_3yr = sum(d$ev == 1 & d$time_yr < 3), observed_beyond_3yr = sum(d$time_yr > 3),
               cases_before_5yr = sum(d$ev == 1 & d$time_yr < 5), observed_beyond_5yr = sum(d$time_yr > 5),
               cases_before_10yr = sum(d$ev == 1 & d$time_yr < 10), observed_beyond_10yr = sum(d$time_yr > 10))
  }))
  comparisons <- list(Clinical = c("M2_Clinical_Base", "M3_Clin_Plus_PPARG"),
                      Clinical_Adipose = c("M4_Clin_Adipose_Base", "M5_Clin_Adipose_Plus_PPARG"))
  paired <- rbindlist(lapply(names(comparisons), function(nm) {
    b <- comparisons[[nm]][1]
    p <- comparisons[[nm]][2]
    c_b <- disc[model == b, C_index]
    c_p <- disc[model == p, C_index]
    lrt <- anova(fits[[b]], fits[[p]])
    cbind(data.table(cohort = cohort_nm, endpoint = "OS", comparison = nm, n = nrow(d), events = sum(d$ev),
                     base_model = b, plus_model = p, C_index_base = c_b, C_index_plus = c_p,
                     LRT_pval = lrt[["Pr(>|Chi|)"]][2]),
          bootstrap_delta_c(d, forms[[b]], forms[[p]], c_b, c_p, 2000L, bootstrap_seed))
  }))
  ph <- rbindlist(lapply(names(fits), function(nm) {
    z <- cox.zph(fits[[nm]], transform = "km")$table
    data.table(cohort = cohort_nm, endpoint = "OS", model = nm, n = nrow(d), events = sum(d$ev),
               transform = "km", term = rownames(z), chisq = z[, "chisq"], df = z[, "df"], pval = z[, "p"])
  }))
  list(disc = disc, paired = paired, ph = ph)
}

tc <- readRDS(file.path(bulk_dir, "tcga_planA_adipose_clin.rds"))
mb <- readRDS(file.path(bulk_dir, "metabric_planA_adipose_clin.rds"))
mb[, Stage_num := STAGE]
results <- list(
  model_rows(tc[ER == "Positive"], "TCGA_ERpos", c("AGE", "Stage_num"), c("ADIPOQ", "PLIN1", "StromalScore"), 20261122L),
  model_rows(mb[ER_clean == "Positive"], "METABRIC_ERpos", c("AGE", "Stage_num", "GRADE"), c("ADIPOQ", "PLIN1"), 20261123L)
)
disc <- rbindlist(lapply(results, `[[`, "disc"))
paired <- rbindlist(lapply(results, `[[`, "paired"))
ph <- rbindlist(lapply(results, `[[`, "ph"))
fitness <- fread(file.path(res_dir, "planA_pparg_adipose_and_multivar_cox.csv"))
fit_counts <- unique(fitness[sample_set == "ERpos_OS_common_complete_cases", .(cohort, n, events)])
stopifnot(isTRUE(all.equal(unique(disc[, .(cohort, n, events)]), fit_counts)),
          all(abs(paired$delta_C - (paired$C_index_plus - paired$C_index_base)) < 1e-12))
fwrite(disc, file.path(res_dir, "planA_erpos_roc_cindex_gain.csv"))
fwrite(paired, file.path(res_dir, "planA_erpos_delta_c_bootstrap.csv"))
fwrite(ph, file.path(res_dir, "planA_erpos_os_ph_diagnostics.csv"))
print(disc[, .(cohort, model, n, events, C_index, AUC_3yr, AUC_5yr, AUC_10yr, pparg_pval)])
print(paired[, .(cohort, comparison, n, events, delta_C, delta_C_LCI95, delta_C_UCI95,
                 optimism_corrected_delta_C, LRT_pval, bootstrap_valid)])
print(ph[cohort == "METABRIC_ERpos"])
# 临床列线图及校准退出投稿图。历史缓存留给旧版绘图与补充表的追溯。
