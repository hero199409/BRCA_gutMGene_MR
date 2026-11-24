# Audit eQTLGen/GTEx directions on the GTEx alternate allele.
source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages(library(data.table))
res_dir <- file.path(project_dir, "02_results")
sig <- fread(file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))
ivs <- fread(file.path(res_dir, "module3_all41_genes_157_iv_hg19_audit.csv"))
targets <- ivs[gene_symbol %in% sig$gene_symbol & is_cis_1mb_hg19 == TRUE,
  .(gene_symbol, ensembl_id, rsid = SNP, eqtlgen_effect_allele = effect_allele.exposure,
    eqtlgen_other_allele = other_allele.exposure, eqtlgen_beta = beta.exposure, eqtlgen_pval = pval.exposure)]
stopifnot(!anyDuplicated(targets, by = c("gene_symbol", "rsid")))
fp <- file.path(res_dir, "module3_cis_positive_genes_gtex_breast_dyneqtl.csv")
raw_cols <- c("gene_symbol", "rsid", "gencodeId", "variantId_b38", "gtex_alt", "gtex_ref", "gtex_nes", "gtex_se", "gtex_pval", "gtex_t")
if (mr_offline) {
  if (!file.exists(fp)) stop("The GTEx breast eQTL snapshot is required in offline mode.")
  ref <- fread(fp)[, ..raw_cols]
} else {
  get_gtex <- function(endpoint, parameters) {
    resp <- httr::GET(paste0("https://gtexportal.org/api/v2/", endpoint), query = parameters, httr::timeout(20))
    httr::stop_for_status(resp)
    jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
  }
  gene_ids <- setNames(lapply(unique(targets$gene_symbol), function(symbol) {
    get_gtex("reference/gene", list(geneId = symbol, gencodeVersion = "v26"))$data$gencodeId[1]
  }), unique(targets$gene_symbol))
  ref <- rbindlist(lapply(seq_len(nrow(targets)), function(idx) {
    symbol <- targets$gene_symbol[idx]
    gene_id <- gene_ids[[symbol]]
    rsid <- targets$rsid[idx]
    if (length(gene_id) != 1L || is.na(gene_id)) return(NULL)
    variant <- get_gtex("dataset/variant", list(snpId = rsid, datasetId = "gtex_v8"))$data
    if (!is.data.frame(variant) || nrow(variant) != 1L) return(NULL)
    fit <- get_gtex("association/dyneqtl", list(tissueSiteDetailId = "Breast_Mammary_Tissue",
      gencodeId = gene_id, variantId = variant$variantId[1]))
    scalar <- function(x) if (length(x) == 1L) as.numeric(x) else NA_real_
    data.table(gene_symbol = symbol, rsid = rsid, gencodeId = gene_id,
      variantId_b38 = variant$variantId[1], gtex_alt = variant$alt[1], gtex_ref = variant$ref[1],
      gtex_nes = scalar(fit$nes), gtex_se = scalar(fit$error),
      gtex_pval = scalar(fit$pValue), gtex_t = scalar(fit$tStatistic))
  }))
  if (!nrow(ref)) stop("GTEx returned no usable records.")
}
stopifnot(!anyDuplicated(ref, by = c("gene_symbol", "rsid")))
aud <- merge(targets, ref, by = c("gene_symbol", "rsid"), all.x = TRUE)
complement <- function(x) chartr("ACGT", "TGCA", x)
for (col in c("eqtlgen_effect_allele", "eqtlgen_other_allele", "gtex_alt", "gtex_ref")) set(aud, j = col, value = toupper(aud[[col]]))
aud[, allele_alignment := fcase(
  is.na(gtex_alt) | is.na(gtex_ref) | is.na(eqtlgen_effect_allele) | is.na(eqtlgen_other_allele), "missing",
  !grepl("^[ACGT]$", gtex_alt) | !grepl("^[ACGT]$", gtex_ref), "non-SNV unresolved",
  complement(gtex_alt) == gtex_ref, "palindromic unresolved",
  eqtlgen_effect_allele == gtex_alt & eqtlgen_other_allele == gtex_ref, "same alleles",
  eqtlgen_effect_allele == gtex_ref & eqtlgen_other_allele == gtex_alt, "swapped alleles",
  complement(eqtlgen_effect_allele) == gtex_alt & complement(eqtlgen_other_allele) == gtex_ref, "same alleles after strand flip",
  complement(eqtlgen_effect_allele) == gtex_ref & complement(eqtlgen_other_allele) == gtex_alt, "swapped alleles after strand flip",
  default = "incompatible")]
aud[, eqtlgen_beta_aligned_to_gtex_alt := fcase(startsWith(allele_alignment, "same"), eqtlgen_beta,
  startsWith(allele_alignment, "swapped"), -eqtlgen_beta, default = NA_real_)]
aud[, eqtlgen_vs_gtex_concordant := sign(gtex_nes) == sign(eqtlgen_beta_aligned_to_gtex_alt)]
fwrite(aud, fp)
print(aud[, .N, by = .(allele_alignment, eqtlgen_vs_gtex_concordant)])
