# 11c_update_module4_complete.R
# 补全 TCGA CXCL8 (Entrez 3576)、GSE20685 (n=327 OS/MFS)、GSE7390 (n=198 OS/RFS) 以及多因素 Cox (年龄+分期校正)

# 独立的本地生存尺度审计入口。不改写旧的模块 4、差异表达及方向审计结果。
if ("--survival-scale-audit" %in% commandArgs(trailingOnly = TRUE)) {
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
  })
  local({
    root <- Sys.getenv("BRCA_PROJECT_ROOT", unset = ".")
    bulk <- file.path(root, "00_data/bulk_cohorts")
    geo <- file.path(root, "00_data/geo_cache")
    genes <- c("PPARG", "CXCL8")
    num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    aligned <- function(expr, pheno, id) {
      idx <- match(colnames(expr), pheno[[id]])
      stopifnot(!anyDuplicated(colnames(expr)), !anyDuplicated(pheno[[id]]), !anyNA(idx))
      as.data.table(pheno)[idx]
    }
    estimate <- function(expr, when, status, cohort, endpoint, unit, scale_to_year,
                         time_source, event_source, strata = list(Overall = seq_along(when)),
                         sensitivity = "positive_times") {
      stopifnot(all(genes %in% rownames(expr)), length(when) == ncol(expr),
                length(status) == ncol(expr), all(is.na(status) | status %in% 0:1))
      rbindlist(lapply(names(strata), function(stratum) rbindlist(lapply(genes, function(gene) {
        i <- strata[[stratum]]
        d <- data.table(x = as.numeric(expr[gene, i]), time = when[i], event = status[i])
        zeros <- d[time == 0 & event == 1, .N]
        d <- d[is.finite(time) & event %in% 0:1 & is.finite(x) &
                 (time > 0 | (sensitivity == "include_recorded_zero_events" & time == 0 & event == 1))]
        stopifnot(nrow(d) > 1L, sum(d$event) > 0L, sd(d$x) > 0)
        expr_sd <- sd(d$x)
        d[, `:=`(z = (x - mean(x)) / expr_sd, time_yr = time * scale_to_year)]
        fit <- coxph(Surv(time_yr, event) ~ z, data = d, na.action = na.fail)
        native <- coxph(Surv(time, event) ~ x, data = d, na.action = na.fail)
        years <- coxph(Surv(time_yr, event) ~ x, data = d, na.action = na.fail)
        # 验证时间单位换算与表达尺度转换对拟合的预期作用。
        stopifnot(isTRUE(all.equal(unname(coef(native)), unname(coef(years)), tolerance = 1e-7)),
                  isTRUE(all.equal(unname(coef(native)) * expr_sd, unname(coef(fit)), tolerance = 1e-7)),
                  isTRUE(all.equal(unname(sqrt(vcov(native)[1, 1])) * expr_sd,
                                   unname(sqrt(vcov(fit)[1, 1])), tolerance = 1e-7)))
        s <- summary(fit)
        data.table(cohort, stratum, endpoint, gene_symbol = gene, sensitivity,
                   n_samples = fit$n, n_events = fit$nevent, zero_time_events = zeros,
                   time_source, event_source, time_unit_source = unit,
                   time_unit_analysis = "years", median_observed_years = median(d$time_yr),
                   max_observed_years = max(d$time_yr), expression_scale = "per within-cohort/stratum endpoint SD",
                   expression_sd_native = expr_sd, log_HR = unname(coef(fit)),
                   SE_log_HR = sqrt(vcov(fit)[1, 1]), HR_native = unname(exp(coef(native))),
                   HR = s$conf.int[1, "exp(coef)"], HR_LCI95 = s$conf.int[1, "lower .95"],
                   HR_UCI95 = s$conf.int[1, "upper .95"], cox_pval = s$coefficients[1, "Pr(>|z|)"])
      }))))
    }
    strata_for <- function(er) list(Overall = seq_along(er), ER_pos = which(er == "Positive"),
                                     ER_neg = which(er == "Negative"))
    tc <- readRDS(file.path(bulk, "tcga_brca_clean_cache.rds"))
    cl <- aligned(tc$tumor_log2, tc$clin, "SAMPLE_ID")
    tc_rows <- rbind(
      estimate(tc$tumor_log2, num(cl$OS_MONTHS), num(cl$OS_event), "TCGA_BRCA", "OS",
               "months", 1/12, "OS_MONTHS", "OS_event", strata_for(cl$ER_Xena)),
      estimate(tc$tumor_log2, num(cl$PFS_MONTHS), num(cl$PFS_event), "TCGA_BRCA", "PFS",
               "months", 1/12, "PFS_MONTHS", "PFS_event", strata_for(cl$ER_Xena)))
    rm(tc, cl)
    mb <- readRDS(file.path(bulk, "metabric_target_cache.rds"))
    cl <- aligned(mb$expr_target, mb$clin, "SAMPLE_ID")
    mb_rows <- rbind(
      estimate(mb$expr_target, num(cl$OS_MONTHS), num(cl$OS_event), "METABRIC", "OS",
               "months", 1/12, "OS_MONTHS", "OS_event", strata_for(cl$ER_clean)),
      estimate(mb$expr_target, num(cl$RFS_MONTHS), num(cl$RFS_event), "METABRIC", "RFS",
               "months", 1/12, "RFS_MONTHS", "RFS_event", strata_for(cl$ER_clean)))
    open_geo <- function(id) {
      g <- readRDS(file.path(geo, paste0(id, "_parsed.rds")))
      g$pheno <- aligned(g$expr_target, g$pheno, "geo_accession")
      g
    }
    geo_endpoint <- function(g, cohort, endpoint, time_col, event_col, unit, scale_to_year, ...) {
      stopifnot(all(c(time_col, event_col) %in% names(g$pheno)))
      estimate(g$expr_target, num(g$pheno[[time_col]]), num(g$pheno[[event_col]]),
               cohort, endpoint, unit, scale_to_year, time_col, event_col, ...)
    }
    g42 <- open_geo("GSE42568")
    g42_rows <- rbind(
      geo_endpoint(g42, "GSE42568", "OS", "overall survival time_days:ch1", "overall survival event:ch1", "days", 1/365.25),
      geo_endpoint(g42, "GSE42568", "RFS", "relapse free survival time_days:ch1", "relapse free survival event:ch1", "days", 1/365.25))
    g73 <- open_geo("GSE7390")
    # GSE7390 本地字段未注明单位。按数值量级以天换算，并校验 HR 对单位换算不变。
    g73_rows <- rbind(
      geo_endpoint(g73, "GSE7390", "OS", "t.os:ch1", "e.os:ch1", "days (inferred from cached values)", 1/365.25),
      geo_endpoint(g73, "GSE7390", "RFS", "t.rfs:ch1", "e.rfs:ch1", "days (inferred from cached values)", 1/365.25))
    g20 <- open_geo("GSE20685")
    os20 <- geo_endpoint(g20, "GSE20685", "OS", "follow_up_duration (years):ch1", "event_death:ch1", "years", 1)
    e_met <- num(g20$pheno[["event_metastasis:ch1"]])
    t_fu <- num(g20$pheno[["follow_up_duration (years):ch1"]])
    t_met <- num(g20$pheno[["time_to_metastasis (years):ch1"]])
    # 只对无转移病例使用随访时长作删失时间。零时点事件另做敏感性分析。
    t_mfs <- fifelse(e_met == 1, t_met, t_fu)
    mfs20 <- rbindlist(lapply(c("positive_times", "include_recorded_zero_events"), function(rule) {
      estimate(g20$expr_target, t_mfs, e_met, "GSE20685", "MFS", "years", 1,
               "time_to_metastasis for events; follow_up_duration for censoring",
               "event_metastasis:ch1", sensitivity = rule)
    }))
    out <- rbindlist(list(tc_rows, mb_rows, g42_rows, g73_rows, os20, mfs20))
    fwrite(out, file.path(root, "02_results/module4_survival_scale_audit.csv"))
    print(out[, .(cohort, stratum, endpoint, gene_symbol, sensitivity, n_samples, n_events,
                  zero_time_events, HR, HR_LCI95, HR_UCI95, cox_pval)])
  })
  quit(save = "no", status = 0L)
}

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

bulk_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
geo_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/geo_cache")
res_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

sig_cis <- fread(file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))
target_genes <- unique(as.character(sig_cis$gene_symbol))

run_surv <- function(e, t, ev, coh, st, ep, g) {
  ok <- !is.na(e) & !is.na(t) & !is.na(ev) & t > 0
  if (sum(ok) < 30 || sum(ev[ok] == 1) < 5 || length(unique(e[ok])) < 5) return(NULL)
  e <- as.numeric(e[ok]); t <- as.numeric(t[ok]); ev <- as.integer(ev[ok])
  cx <- tryCatch(summary(coxph(Surv(t, ev) ~ e)), error = function(err) NULL)
  if (is.null(cx)) return(NULL)
  grp <- ifelse(e >= median(e), "High", "Low")
  kp <- if (length(unique(grp)) == 2) tryCatch(pchisq(survdiff(Surv(t, ev) ~ grp)$chisq, df = 1, lower.tail = FALSE), error = function(err) NA_real_) else NA_real_
  data.table(cohort = coh, stratum = st, endpoint = ep, gene_symbol = g,
             n_samples = length(e), n_events = sum(ev == 1),
             HR = cx$conf.int[1, "exp(coef)"], HR_LCI95 = cx$conf.int[1, "lower .95"], HR_UCI95 = cx$conf.int[1, "upper .95"],
             cox_pval = cx$coefficients[1, "Pr(>|z|)"], km_logrank_pval = kp)
}

# 1. 补全 TCGA CXCL8 差异/亚型/生存 + 全 8 基因多因素 Cox
tc <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
cl <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]
cl[, stage_num := fifelse(grepl("Stage IV", AJCC_PATHOLOGIC_TUMOR_STAGE), 4L,
                  fifelse(grepl("Stage III", AJCC_PATHOLOGIC_TUMOR_STAGE), 3L,
                  fifelse(grepl("Stage II", AJCC_PATHOLOGIC_TUMOR_STAGE), 2L,
                  fifelse(grepl("Stage I", AJCC_PATHOLOGIC_TUMOR_STAGE), 1L, NA_integer_))))]

diff_old <- fread(file.path(res_dir, "module4_differential_expression.csv"))
surv_old <- fread(file.path(res_dir, "module4_survival_analysis.csv"))

# 补 CXCL8 TCGA
if ("CXCL8" %in% rownames(tc$tumor_log2)) {
  vt <- as.numeric(tc$tumor_log2["CXCL8", ])
  vn <- as.numeric(tc$normal_log2["CXCL8", ])
  pat_t <- substr(colnames(tc$tumor_log2), 1, 12)
  pat_n <- substr(colnames(tc$normal_log2), 1, 12)
  cp <- intersect(pat_t, pat_n)
  vt_p <- sapply(cp, function(p) mean(vt[pat_t == p], na.rm = TRUE))
  vn_p <- sapply(cp, function(p) mean(vn[pat_n == p], na.rm = TRUE))
  
  new_diffs <- list()
  new_diffs[[1]] <- data.table(
    cohort = "TCGA_BRCA", comparison = "Tumor_vs_Normal", gene_symbol = "CXCL8",
    n_group1 = sum(!is.na(vt)), n_group2 = sum(!is.na(vn)),
    mean_group1 = round(mean(vt, na.rm=TRUE), 3), mean_group2 = round(mean(vn, na.rm=TRUE), 3),
    log2FC = round(mean(vt, na.rm=TRUE) - mean(vn, na.rm=TRUE), 3),
    wilcox_pval = wilcox.test(vt, vn)$p.value,
    paired_n = length(cp), paired_log2FC = round(mean(vt_p - vn_p, na.rm=TRUE), 3),
    paired_pval = wilcox.test(vt_p, vn_p, paired = TRUE)$p.value
  )
  idx_erp <- which(cl$ER_Xena == "Positive")
  idx_ern <- which(cl$ER_Xena == "Negative")
  new_diffs[[2]] <- data.table(
    cohort = "TCGA_BRCA", comparison = "ERpos_vs_ERneg", gene_symbol = "CXCL8",
    n_group1 = length(idx_erp), n_group2 = length(idx_ern),
    mean_group1 = round(mean(vt[idx_erp], na.rm=TRUE), 3), mean_group2 = round(mean(vt[idx_ern], na.rm=TRUE), 3),
    log2FC = round(mean(vt[idx_erp], na.rm=TRUE) - mean(vt[idx_ern], na.rm=TRUE), 3),
    wilcox_pval = wilcox.test(vt[idx_erp], vt[idx_ern])$p.value,
    paired_n = NA_integer_, paired_log2FC = NA_real_, paired_pval = NA_real_
  )
  diff_old <- diff_old[!(cohort == "TCGA_BRCA" & gene_symbol == "CXCL8")]
  diff_old <- rbindlist(list(diff_old, rbindlist(new_diffs)), fill = TRUE)
  
  cxcl8_surv <- list()
  strata_idx <- list(Overall = seq_len(nrow(cl)), ER_pos = idx_erp, ER_neg = idx_ern)
  for (st in names(strata_idx)) {
    ii <- strata_idx[[st]]
    cxcl8_surv[[paste(st,"OS")]]  <- run_surv(vt[ii], cl$OS_MONTHS[ii], cl$OS_event[ii], "TCGA_BRCA", st, "OS", "CXCL8")
    cxcl8_surv[[paste(st,"DSS")]] <- run_surv(vt[ii], cl$DSS_MONTHS[ii], cl$DSS_event[ii], "TCGA_BRCA", st, "DSS", "CXCL8")
    cxcl8_surv[[paste(st,"PFS")]] <- run_surv(vt[ii], cl$PFS_MONTHS[ii], cl$PFS_event[ii], "TCGA_BRCA", st, "PFS", "CXCL8")
  }
  surv_old <- surv_old[!(cohort == "TCGA_BRCA" & gene_symbol == "CXCL8")]
  surv_old <- rbindlist(list(surv_old, rbindlist(cxcl8_surv)), fill = TRUE)
}

# 多因素 Cox (TCGA: expr + AGE + stage_num)
mv_rows <- list()
for (g in target_genes) {
  if (!(g %in% rownames(tc$tumor_log2))) next
  e <- as.numeric(tc$tumor_log2[g, ])
  for (st in c("Overall", "ER_pos")) {
    ii <- if (st == "Overall") seq_len(nrow(cl)) else which(cl$ER_Xena == "Positive")
    ok <- ii[!is.na(e[ii]) & !is.na(cl$OS_MONTHS[ii]) & !is.na(cl$OS_event[ii]) & cl$OS_MONTHS[ii] > 0 &
               !is.na(cl$AGE[ii]) & !is.na(cl$stage_num[ii])]
    if (length(ok) >= 50 && sum(cl$OS_event[ok] == 1) >= 10) {
      cx <- tryCatch(summary(coxph(Surv(cl$OS_MONTHS[ok], cl$OS_event[ok]) ~ e[ok] + cl$AGE[ok] + cl$stage_num[ok])), error = function(err) NULL)
      if (!is.null(cx)) {
        mv_rows[[paste("TCGA", st, g)]] <- data.table(
          cohort = "TCGA_BRCA", stratum = st, endpoint = "OS_Multivar(Age+Stage)", gene_symbol = g,
          n = length(ok), events = sum(cl$OS_event[ok] == 1),
          HR = round(cx$conf.int[1, "exp(coef)"], 3),
          HR_LCI95 = round(cx$conf.int[1, "lower .95"], 3),
          HR_UCI95 = round(cx$conf.int[1, "upper .95"], 3),
          pval = signif(cx$coefficients[1, "Pr(>|z|)"], 3)
        )
      }
    }
  }
}
rm(tc); gc()

# METABRIC 多因素 Cox (expr + AGE_AT_DIAGNOSIS + TUMOR_STAGE)
mb <- readRDS(file.path(bulk_dir, "metabric_target_cache.rds"))
mcl <- mb$clin
mcl[, stage_num := suppressWarnings(as.integer(TUMOR_STAGE))]
for (g in target_genes) {
  if (!(g %in% rownames(mb$expr_target))) next
  e <- as.numeric(mb$expr_target[g, ])
  for (st in c("Overall", "ER_pos")) {
    ii <- if (st == "Overall") seq_len(nrow(mcl)) else which(mcl$ER_clean == "Positive")
    ok <- ii[!is.na(e[ii]) & !is.na(mcl$OS_MONTHS[ii]) & !is.na(mcl$OS_event[ii]) & mcl$OS_MONTHS[ii] > 0 &
               !is.na(mcl$AGE_AT_DIAGNOSIS[ii]) & !is.na(mcl$stage_num[ii])]
    if (length(ok) >= 50 && sum(mcl$OS_event[ok] == 1) >= 10) {
      cx <- tryCatch(summary(coxph(Surv(mcl$OS_MONTHS[ok], mcl$OS_event[ok]) ~ e[ok] + mcl$AGE_AT_DIAGNOSIS[ok] + mcl$stage_num[ok])), error = function(err) NULL)
      if (!is.null(cx)) {
        mv_rows[[paste("METABRIC", st, g)]] <- data.table(
          cohort = "METABRIC", stratum = st, endpoint = "OS_Multivar(Age+Stage)", gene_symbol = g,
          n = length(ok), events = sum(mcl$OS_event[ok] == 1),
          HR = round(cx$conf.int[1, "exp(coef)"], 3),
          HR_LCI95 = round(cx$conf.int[1, "lower .95"], 3),
          HR_UCI95 = round(cx$conf.int[1, "upper .95"], 3),
          pval = signif(cx$coefficients[1, "Pr(>|z|)"], 3)
        )
      }
    }
  }
}
rm(mb); gc()

mv_dt <- rbindlist(mv_rows)
fwrite(mv_dt, file.path(res_dir, "module4_multivariable_cox.csv"))
cat("=== 多因素 Cox (校正年龄+分期) 显著结果 (P < 0.05) ===\n")
print(mv_dt[pval < 0.05][order(pval)])

# 2. 解析 GSE20685 (follow_up_duration (years):ch1, event_death:ch1, event_metastasis:ch1)
g20 <- readRDS(file.path(geo_dir, "GSE20685_parsed.rds"))
pd20 <- as.data.table(g20$pheno)
t_fu  <- suppressWarnings(as.numeric(as.character(pd20[["follow_up_duration (years):ch1"]])))
e_os  <- suppressWarnings(as.integer(as.character(pd20[["event_death:ch1"]])))
t_met <- suppressWarnings(as.numeric(as.character(pd20[["time_to_metastasis (years):ch1"]])))
e_met <- suppressWarnings(as.integer(as.character(pd20[["event_metastasis:ch1"]])))
t_mfs <- ifelse(!is.na(t_met) & t_met > 0, t_met, t_fu)

g20_surv <- list()
for (g in rownames(g20$expr_target)) {
  v <- as.numeric(g20$expr_target[g, ])
  g20_surv[[paste("OS", g)]]  <- run_surv(v, t_fu, e_os, "GSE20685", "Overall", "OS", g)
  g20_surv[[paste("MFS", g)]] <- run_surv(v, t_mfs, e_met, "GSE20685", "Overall", "MFS", g)
}
g20_dt <- rbindlist(g20_surv)

g73_dt <- if (file.exists(file.path(res_dir, "module4_geo_extra_survival.csv"))) fread(file.path(res_dir, "module4_geo_extra_survival.csv")) else data.table()

surv_all <- surv_old[!(cohort %in% c("GSE20685", "GSE7390"))]
surv_all <- rbindlist(list(surv_all, g20_dt, g73_dt), fill = TRUE)

fwrite(diff_old, file.path(res_dir, "module4_differential_expression.csv"))
fwrite(surv_all, file.path(res_dir, "module4_survival_analysis.csv"))

# 更新 GEO 状态表
geo_st <- data.table(
  cohort = c("GSE42568", "GSE20685", "GSE7390"),
  status = "SUCCESS",
  n_samples = c(121L, ncol(g20$expr_target), 198L),
  notes = c("104 Tumor + 17 Normal, OS/RFS", "327 Tumor, OS/MFS", "198 Tumor (node-negative), OS/RFS")
)
fwrite(geo_st, file.path(res_dir, "module4_geo_status.csv"))

cat("\n=== GSE20685 (n=327) 生存验证结果 ===\n")
print(g20_dt[order(cox_pval), .(endpoint, gene_symbol, n_samples, n_events, HR = round(HR, 3), cox_pval = signif(cox_pval, 3), km_pval = signif(km_logrank_pval, 3))])

# 3. 更新最终模块 4 多队列方向审计总表
coloc_dt <- fread(file.path(res_dir, "module3_cis_positive_genes_coloc.csv"))
audit_rows <- list()
for (g in target_genes) {
  best_mr <- sig_cis[gene_symbol == g][order(pval)][1]
  col_sub <- coloc_dt[gene_symbol == g]
  max_h4 <- if (nrow(col_sub) > 0) max(col_sub$PP.H4, na.rm = TRUE) else NA_real_
  
  tcga_tn <- diff_old[cohort == "TCGA_BRCA" & comparison == "Tumor_vs_Normal" & gene_symbol == g]
  tcga_er <- diff_old[cohort == "TCGA_BRCA" & comparison == "ERpos_vs_ERneg" & gene_symbol == g]
  gse_tn  <- diff_old[cohort == "GSE42568" & comparison == "Tumor_vs_Normal" & gene_symbol == g]
  
  tcga_os  <- surv_all[cohort == "TCGA_BRCA" & stratum == "Overall" & endpoint == "OS" & gene_symbol == g]
  tcga_erp <- surv_all[cohort == "TCGA_BRCA" & stratum == "ER_pos" & endpoint == "OS" & gene_symbol == g]
  mb_os    <- surv_all[cohort == "METABRIC" & stratum == "Overall" & endpoint == "OS" & gene_symbol == g]
  mb_erp   <- surv_all[cohort == "METABRIC" & stratum == "ER_pos" & endpoint == "OS" & gene_symbol == g]
  g20_os   <- surv_all[cohort == "GSE20685" & endpoint == "OS" & gene_symbol == g]
  g73_os   <- surv_all[cohort == "GSE7390" & endpoint == "OS" & gene_symbol == g]
  
  mr_vs_tn <- if (nrow(tcga_tn) > 0) ((best_mr$OR > 1 && tcga_tn$log2FC > 0) || (best_mr$OR < 1 && tcga_tn$log2FC < 0)) else NA
  mr_vs_os <- if (nrow(mb_os) > 0) ((best_mr$OR > 1 && mb_os$HR > 1) || (best_mr$OR < 1 && mb_os$HR < 1)) else NA
  
  audit_rows[[g]] <- data.table(
    gene_symbol = g,
    upstream_metabolites = best_mr$upstream_metabolites,
    gutmgene_alteration = best_mr$alterations,
    best_MR_outcome = best_mr$outcome_label,
    MR_OR = round(best_mr$OR, 3),
    MR_pval = signif(best_mr$pval, 3),
    MR_FDR = signif(best_mr$FDR, 3),
    coloc_max_PP.H4 = round(max_h4, 4),
    TCGA_TN_log2FC = if (nrow(tcga_tn) > 0) tcga_tn$log2FC else NA_real_,
    TCGA_TN_pval = if (nrow(tcga_tn) > 0) signif(tcga_tn$wilcox_pval, 3) else NA_real_,
    GSE42568_TN_log2FC = if (nrow(gse_tn) > 0) gse_tn$log2FC else NA_real_,
    GSE42568_TN_pval = if (nrow(gse_tn) > 0) signif(gse_tn$wilcox_pval, 3) else NA_real_,
    TCGA_OS_HR = if (nrow(tcga_os) > 0) round(tcga_os$HR, 3) else NA_real_,
    TCGA_OS_pval = if (nrow(tcga_os) > 0) signif(tcga_os$cox_pval, 3) else NA_real_,
    TCGA_ERpos_OS_HR = if (nrow(tcga_erp) > 0) round(tcga_erp$HR, 3) else NA_real_,
    TCGA_ERpos_OS_pval = if (nrow(tcga_erp) > 0) signif(tcga_erp$cox_pval, 3) else NA_real_,
    METABRIC_OS_HR = if (nrow(mb_os) > 0) round(mb_os$HR, 3) else NA_real_,
    METABRIC_OS_pval = if (nrow(mb_os) > 0) signif(mb_os$cox_pval, 3) else NA_real_,
    METABRIC_ERpos_OS_HR = if (nrow(mb_erp) > 0) round(mb_erp$HR, 3) else NA_real_,
    METABRIC_ERpos_OS_pval = if (nrow(mb_erp) > 0) signif(mb_erp$cox_pval, 3) else NA_real_,
    GSE20685_OS_HR = if (nrow(g20_os) > 0) round(g20_os$HR, 3) else NA_real_,
    GSE20685_OS_pval = if (nrow(g20_os) > 0) signif(g20_os$cox_pval, 3) else NA_real_,
    GSE7390_OS_HR = if (nrow(g73_os) > 0) round(g73_os$HR, 3) else NA_real_,
    GSE7390_OS_pval = if (nrow(g73_os) > 0) signif(g73_os$cox_pval, 3) else NA_real_,
    MR_vs_TCGA_TN_consistent = mr_vs_tn,
    MR_vs_METABRIC_OS_consistent = mr_vs_os
  )
}
audit_dt <- rbindlist(audit_rows)
fwrite(audit_dt, file.path(res_dir, "module4_direction_audit_summary.csv"))

cat("\n=== 最终模块 4 五队列方向审计总表 ===\n")
print(audit_dt[, .(gene_symbol, best_MR_outcome, MR_OR, MR_pval, MR_vs_TCGA_TN_consistent, MR_vs_METABRIC_OS_consistent,
                   TCGA_ERpos_OS_HR, TCGA_ERpos_OS_pval, METABRIC_ERpos_OS_HR, METABRIC_ERpos_OS_pval,
                   GSE20685_OS_HR, GSE20685_OS_pval, GSE7390_OS_HR, GSE7390_OS_pval)])
