# 05_ispy2_missing_pcr_32scenarios.R
# Exhaustive 2^5 = 32 sensitivity analysis for the 5 control-arm I-SPY2 (GSE173839)
# FFPE patients with missing pathological response (pCR coded as -1).

suppressPackageStartupMessages({
  library(data.table)
})

proj_root <- Sys.getenv("BRCA_PROJECT_ROOT", unset = "")
if (!nzchar(proj_root)) {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  proj_root <- if (length(file_arg) > 0) normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..", ".."), winslash = "/", mustWork = FALSE) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
cache_rds <- file.path(proj_root, "00_data", "bulk_cohorts", "planA_therapy_plot_cache.rds")
out_csv   <- file.path(proj_root, "05_revision_v2", "ispy2_missing_pcr_32scenarios_sensitivity.csv")

th_pc <- readRDS(cache_rds)
df <- copy(as.data.table(th_pc$gse173839))

cat("Total rows in gse173839 cache:", nrow(df), "\n")
cat("pCR table by Arm:\n")
print(table(df$Arm, df$pCR, useNA = "ifany"))

miss_idx <- which(!(df$pCR %in% c(0L, 1L)))
cat("Missing pCR indices count:", length(miss_idx), "\n")
stopifnot(length(miss_idx) == 5L)

# Check column names for HR status and MammaPrint
cat("Columns:", paste(names(df), collapse = ", "), "\n")
hr_col <- intersect(c("HR_status", "HR", "HR_pos", "HR_neg", "hr_status"), names(df))[1]
mp_col <- intersect(c("MP", "MammaPrint", "MP_high", "MP2", "mp_cat"), names(df))[1]
cat("Using HR col:", hr_col, "| MP col:", mp_col, "\n")

df[, Arm_bin := as.integer(Arm == "durvalumab/olaparib")]
df[, HR_bin  := as.integer(as.logical(get(hr_col)))]
if (!is.na(mp_col)) {
  df[, MP_bin := as.integer(as.factor(get(mp_col))) - 1L]
}

grid32 <- as.matrix(expand.grid(rep(list(c(0L, 1L)), 5L)))
res_list <- vector("list", nrow(grid32))

for (i in seq_len(nrow(grid32))) {
  d_i <- copy(df)
  d_i[miss_idx, pCR := grid32[i, ]]
  n_ctrl_pcr <- sum(d_i[Arm_bin == 0L]$pCR == 1L)

  m_hr <- glm(pCR ~ CXCL8 * Arm_bin + HR_bin, data = d_i, family = binomial())
  s_hr <- summary(m_hr)$coefficients["CXCL8:Arm_bin", ]
  or_hr  <- exp(s_hr["Estimate"])
  lci_hr <- exp(s_hr["Estimate"] - 1.96 * s_hr["Std. Error"])
  uci_hr <- exp(s_hr["Estimate"] + 1.96 * s_hr["Std. Error"])
  p_hr   <- s_hr["Pr(>|z|)"]

  or_mp <- lci_mp <- uci_mp <- p_mp <- NA_real_
  if (!is.na(mp_col)) {
    m_mp <- glm(pCR ~ CXCL8 * Arm_bin + HR_bin + MP_bin, data = d_i, family = binomial())
    s_mp <- summary(m_mp)$coefficients["CXCL8:Arm_bin", ]
    or_mp  <- exp(s_mp["Estimate"])
    lci_mp <- exp(s_mp["Estimate"] - 1.96 * s_mp["Std. Error"])
    uci_mp <- exp(s_mp["Estimate"] + 1.96 * s_mp["Std. Error"])
    p_mp   <- s_mp["Pr(>|z|)"]
  }

  res_list[[i]] <- data.table(
    scenario = i,
    imputed_pattern = paste(grid32[i, ], collapse = ""),
    imputed_pCR_count = sum(grid32[i, ]),
    total_ctrl_pCR = n_ctrl_pcr,
    OR_HR = or_hr, LCI_HR = lci_hr, UCI_HR = uci_hr, P_HR = p_hr,
    OR_HR_MP = or_mp, LCI_HR_MP = lci_mp, UCI_HR_MP = uci_mp, P_HR_MP = p_mp
  )
}

res <- rbindlist(res_list)
fwrite(res, out_csv)

cat("\n=== 32-SCENARIO SUMMARY (CXCL8 x Arm) ===\n")
cat(sprintf("HR-adjusted:    OR range = [%.3f, %.3f], P range = [%.5f, %.5f], n(P<0.05) = %d/32\n",
            min(res$OR_HR), max(res$OR_HR), min(res$P_HR), max(res$P_HR), sum(res$P_HR < 0.05)))
if (!all(is.na(res$P_HR_MP))) {
  cat(sprintf("HR+MP-adjusted: OR range = [%.3f, %.3f], P range = [%.5f, %.5f], n(P<0.05) = %d/32\n",
              min(res$OR_HR_MP), max(res$OR_HR_MP), min(res$P_HR_MP), max(res$P_HR_MP), sum(res$P_HR_MP < 0.05)))
}
