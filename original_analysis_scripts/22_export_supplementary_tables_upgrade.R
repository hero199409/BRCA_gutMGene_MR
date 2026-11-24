# Export the current tables to the journal bundle and derive the inventory from disk.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("BRCA_PROJECT_ROOT", Sys.getenv("BRCA_PROJECT_ROOT", unset = "."))
res <- file.path(root, "02_results")
out <- file.path(res, "Supplementary_Tables_PlanA")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

spec <- fread(text = "id|label|source
S1|MiBioGen BCAC primary MR|mibiogen_bcac_mr_primary.csv
S2|MiBioGen reverse MR|mibiogen_reverse_mr.csv
S3|MiBioGen BCAC FinnGen meta-analysis|mibiogen_bcac_finngen_meta.csv
S4|Metabolite BCAC primary MR|metabolite_bcac_mr_primary.csv
S5a|Metabolite reverse MR|metabolite_reverse_mr.csv
S5b|Metabolite FinnGen MR|metabolite_finngen_validation_mr.csv
S6|Taxon to metabolite 30 primary MR tests|mediation_step1_taxa_to_metabolite_primary.csv
S7|Curated host gene candidates|module3_candidate_host_genes_universe.csv
S8|Blood cis MR 108 primary tests|module3_eqtlgen_cis_only_mr_primary_full_denominator.csv
S9|GRCh37 audit of 157 instrumental variants|module3_all41_genes_157_iv_hg19_audit.csv
S10|Selected MR sensitivity models|planA_mr_sensitivity_summary.csv
S11|Single variant Wald estimates|planA_mr_single_snp_wald.csv
S12a|Blood eQTL colocalization across priors|module3_cis_coloc_prior_sensitivity.csv
S12b|GTEx v8 independent breast eQTL|module3_cis_positive_genes_gtex_breast_eqtl.csv
S12c|GTEx v8 aligned dynamic breast eQTL|module3_cis_positive_genes_gtex_breast_dyneqtl.csv
S13|PAM50 subtype comparisons|planA_pam50_subtypes_stats.csv
S14|Univariable Cox and Kaplan Meier tests|module4_survival_analysis.csv
S15|Standardized OS meta-analysis|planA_standardized_os_meta.csv
S16|HPA protein and RNA summary|planA_hpa_protein_ihc_and_rna.csv
S17|HPA single-cell reference expression|planA_hpa_single_cell_ncpm.csv
S18|Hallmark GSEA|module5_hallmark_gsea_results.csv
S19|Full GDSC2 subtype-adjusted results|planA_gdsc2_subtype_adjusted_drugs.csv
S20|Full unadjusted GDSC2 tests|module6_gdsc2_all_drugs_correlation.csv
S21|DMP and FR02 primary MR tests|module3_dmp_fr02_bcac_mr_full_denominator.csv
S22a|Selected taxa in DMP and FR02|module3_dmp_fr02_core6_replication.csv
S22b|Same-rank and species substitution microbiome meta-analysis|module3_microbiome_cross_cohort_meta_strict.csv
S23a|Plasma IL-8 pQTL MR|planA_pqtl_il8_mr.csv
S23b|SCALLOP regional colocalization and LD availability|planA_pqtl_il8_availability.csv
S24|ER-positive Cox models|planA_pparg_adipose_and_multivar_cox.csv
S25a|Single-cell lineage summaries|planA_scrna_celltype_summary.csv
S25b|Wu donor-level subtype tests|planA_scrna_wu_donor_subtype_tests.csv
S26|CXCL8 OS heterogeneity and leave-one-out|planA_cxcl8_meta_heterogeneity_audit.csv
S27a|ER-positive apparent model performance|planA_erpos_roc_cindex_gain.csv
S27b|Paired bootstrap C-index differences|planA_erpos_delta_c_bootstrap.csv
S28a|Therapy cohort association models and directional AUC|planA_therapy_response_summary.csv
S28b|Treatment interactions and multivariable models|planA_therapy_interactions_and_multivar.csv
S29|Study metadata|module3_gwas_study_metadata.csv
S30|GWAS sample overlap audit|module3_gwas_overlap_audit.csv
S31|Liability-scale Steiger sensitivity|module3_gwas_metadata_and_steiger_sensitivity.csv
S32|ER subtype contrast covariance sensitivity|module3_cis_er_subtype_heterogeneity.csv
S33|Cis IV harmonization audit|module3_cis_harmonization_drop_audit.csv
S34|Wu donor lineage expression|planA_scrna_wu_donor_subtype_summary.csv
S35|Wu and TISCH2 myeloid composition|planA_scrna_wu_tisch2_myeloid_composition.csv
S36|GSE148673 DCIS exclusion|planA_scrna_gse148673_ecis05t_sensitivity.csv
S37|I-SPY2 pCR assessment flow|planA_gse173839_pcr_flow.csv
S38|Standardized endpoint-specific survival|planA_standardized_non_os_endpoints.csv
S39|Proportional hazards diagnostics|planA_erpos_os_ph_diagnostics.csv
S40|Plasma IL-8 instrument audit|planA_pqtl_il8_iv_audit.csv
S41|Tumor immune score correlations|module5_tme_checkpoint_correlations.csv
S42|quanTIseq input column-sum sensitivity|module5_tme_tpm_sensitivity.csv
S43|Blood cis MR FinnGen overall disease check|module3_cis_positive_genes_finngen_validation.csv", sep = "|")

valid <- file.exists(file.path(res, spec$source))
if (!all(valid)) message("Unavailable optional exports: ", paste(spec$source[!valid], collapse = ", "))
spec <- spec[valid]
# Write one current file for each table ID, retaining previous exports in a separate folder.
archive <- file.path(res, "archived_pre_revision_supplementary_tables")
dir.create(archive, showWarnings = FALSE)
for (i in seq_len(nrow(spec))) {
  entry <- spec[i]
  name <- paste0("Table_", entry$id, "_", gsub("[^A-Za-z0-9]+", "_", entry$label), ".csv")
  existing <- list.files(out, pattern = paste0("^Table_", entry$id, "_.*[.]csv$"), full.names = TRUE)
  case_match <- basename(existing)[tolower(basename(existing)) == tolower(name)]
  if (length(case_match) == 1L) name <- case_match
  prior <- existing[basename(existing) != name]
  d <- fread(file.path(res, entry$source), na.strings = c("", "NA"))
  stopifnot(ncol(d) > 0L)
  fwrite(d, file.path(out, name), na = "")
  if (length(prior)) {
    backups <- file.path(archive, basename(prior))
    if (any(file.exists(backups))) stop("Archived table name already exists: ", paste(backups, collapse = ", "))
    stopifnot(all(file.rename(prior, backups)))
  }
  spec[i, `:=`(file = name, rows = nrow(d), columns = ncol(d))]
}
fwrite(spec, file.path(res, "supplementary_table_index.csv"), na = "")
cat("Supplementary CSV exports:", nrow(spec), "\n")
