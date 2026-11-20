# 22_export_supplementary_tables_upgrade.R
# 整理并导出抗审稿 5 大加餐模块的新增附表 (Table_S21 - Table_S28)

suppressPackageStartupMessages(library(data.table))

res_dir <- "02_results"
sup_dir <- file.path(res_dir, "Supplementary_Tables_PlanA")
dir.create(sup_dir, recursive = TRUE, showWarnings = FALSE)

copy_if_exists <- function(src_name, dst_name) {
  src <- file.path(res_dir, src_name)
  dst <- file.path(sup_dir, dst_name)
  if (file.exists(src)) {
    dt <- fread(src)
    fwrite(dt, dst)
    cat(sprintf("Exported %-58s (%d rows)\n", dst_name, nrow(dt)))
  } else {
    cat("WARNING missing:", src_name, "\n")
  }
}

copy_if_exists("module3_dmp_fr02_bcac_mr_full_denominator.csv", "Table_S21_DMP2022_FR022022_360pairs_FullDenominator.csv")
copy_if_exists("module3_dmp_fr02_core6_replication.csv",        "Table_S22a_Core6_DMP_FR02_Replication.csv")
copy_if_exists("module3_microbiome_cross_cohort_meta.csv",       "Table_S22b_Microbiome_CrossCohort_Meta.csv")
copy_if_exists("planA_pqtl_il8_mr.csv",                          "Table_S23a_CXCL8_IL8_Plasma_pQTL_MR.csv")
copy_if_exists("planA_pqtl_il8_coloc.csv",                       "Table_S23b_SCALLOP_IL8_pQTL_Coloc.csv")
copy_if_exists("planA_pparg_adipose_and_multivar_cox.csv",       "Table_S24_PPARG_Adipose_Confounding_and_MultivarCox.csv")
copy_if_exists("planA_scrna_celltype_summary.csv",               "Table_S25a_Tumor_scRNA_TISCH2_CellType_Summary.csv")
copy_if_exists("planA_scrna_subtype_comparison.csv",             "Table_S25b_Tumor_scRNA_TISCH2_Subtype_Comparison.csv")
copy_if_exists("planA_cxcl8_meta_heterogeneity_audit.csv",       "Table_S26_CXCL8_OS_Meta_Heterogeneity_Decomposition.csv")
copy_if_exists("planA_erpos_roc_cindex_gain.csv",                "Table_S27a_ERpos_Nomogram_ROC_Cindex_Gain.csv")
copy_if_exists("planA_tcga_erpos_nomogram_5yr_cal.csv",          "Table_S27b_TCGA_ERpos_5yr_Bootstrap_Calibration.csv")
copy_if_exists("planA_therapy_response_summary.csv",             "Table_S28a_Clinical_Therapy_Response_Cohorts_58tests.csv")
copy_if_exists("planA_therapy_interactions_and_multivar.csv",    "Table_S28b_Therapy_Interaction_and_Multivar_Tests.csv")
