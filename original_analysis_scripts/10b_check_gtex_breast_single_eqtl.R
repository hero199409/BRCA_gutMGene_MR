# Exploratory BCAC lookup for the FFAR2 normal-breast independent lead eQTL.
# Script 10c handles the eQTLGen-targeted allelic comparisons.
source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages(library(data.table))
res_dir <- file.path(project_dir, "02_results")
cache_dir <- file.path(project_dir, "00_data", "mr_cache")
lead <- fread(file.path(res_dir, "module3_cis_positive_genes_gtex_breast_eqtl.csv"))
lead <- lead[gene_symbol == "FFAR2" & snpId == "rs73927976" & tissueSiteDetailId == "Breast_Mammary_Tissue"]
stopifnot(nrow(lead) == 1L, is.finite(lead$nes), is.finite(lead$pValue), lead$pValue > 0, lead$pValue < 1)

outcomes <- c("ieu-a-1126", "ieu-a-1127", "ieu-a-1128")
cache <- file.path(cache_dir, "outcome_bcac_gtex_ffar2_rs73927976.rds")
if (file.exists(cache)) {
  outcome <- readRDS(cache)
} else if (!mr_offline) {
  outcome <- TwoSampleMR::extract_outcome_data(snps = lead$snpId, outcomes = outcomes, proxies = FALSE)
  if (is.data.frame(outcome) && nrow(outcome) > 0L) saveRDS(outcome, cache)
} else {
  outcome <- NULL
}
status <- data.table(gene = "FFAR2", gtex_rsid = lead$snpId, nes = lead$nes, pval = lead$pValue,
                     bcac_status = "Direct marker summary statistics unavailable in this run")
alleles <- strsplit(lead$variantId, "_", fixed = TRUE)[[1L]]
stopifnot(length(alleles) == 5L, alleles[5L] == "b38", all(alleles[3:4] %chin% c("A", "C", "G", "T")))
if (is.data.frame(outcome) && nrow(outcome) > 0L) {
  stopifnot(all(outcome$SNP == lead$snpId), all(outcome$id.outcome %chin% outcomes))
  exposure <- data.frame(SNP = lead$snpId, id.exposure = "GTEx_Breast_FFAR2", exposure = "FFAR2 GTEx breast NES",
    beta.exposure = lead$nes, se.exposure = abs(lead$nes) / qnorm(lead$pValue / 2, lower.tail = FALSE),
    effect_allele.exposure = alleles[4L], other_allele.exposure = alleles[3L],
    eaf.exposure = NA_real_, pval.exposure = lead$pValue)
  harmonized <- TwoSampleMR::harmonise_data(exposure, outcome, action = 2)
  if (any(harmonized$mr_keep %in% TRUE)) {
    result <- as.data.table(TwoSampleMR::mr(harmonized, method_list = "mr_wald_ratio"))
    result[, `:=`(OR = exp(b), OR_LCI95 = exp(b - 1.96 * se), OR_UCI95 = exp(b + 1.96 * se),
                  exposure_se_method = "Approximate normal inversion of the GTEx p value")]
    fwrite(result, file.path(res_dir, "module3_ffar2_gtex_breast_mr.csv"))
    status[, bcac_status := "Direct marker Wald estimates available"]
  }
}
fwrite(status, file.path(res_dir, "module3_ffar2_gtex_breast_mr_availability.csv"))
print(status)
