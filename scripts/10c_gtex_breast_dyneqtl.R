# 10c_gtex_breast_dyneqtl.R
# 用 GTEx API v2 /association/dyneqtl 对 8 个纯 cis-MR 阳性基因的 13 个 cis-IV SNP 在正常乳腺组织 (Breast_Mammary_Tissue) 做完整查证

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

res_dir <- "02_results"
sig_cis <- fread(file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))
pos_genes <- unique(sig_cis$gene_symbol)
cis_audit <- fread(file.path(res_dir, "module3_all41_genes_iv_cis_audit.csv"))
pos_snps <- cis_audit[gene_symbol %in% pos_genes & is_cis_1mb == TRUE]

# 缓存每个基因的 gencodeId
gmap <- list()
for (sym in pos_genes) {
  ref_cmd <- sprintf('curl -4 -s --connect-timeout 6 --max-time 10 "https://gtexportal.org/api/v2/reference/gene?geneId=%s&gencodeVersion=v26"', sym)
  ref_js  <- tryCatch(jsonlite::fromJSON(paste(suppressWarnings(system(ref_cmd, intern = TRUE)), collapse = "")), error = function(e) NULL)
  if (!is.null(ref_js) && length(ref_js$data) > 0 && "gencodeId" %in% names(ref_js$data)) {
    gmap[[sym]] <- ref_js$data$gencodeId[1]
  }
}

rows <- list()
for (i in seq_len(nrow(pos_snps))) {
  sym <- pos_snps$gene_symbol[i]
  rs  <- pos_snps$SNP[i]
  gid <- gmap[[sym]]
  if (is.null(gid)) next
  
  var_cmd <- sprintf('curl -4 -s --connect-timeout 6 --max-time 10 "https://gtexportal.org/api/v2/dataset/variant?snpId=%s&datasetId=gtex_v8"', rs)
  var_js  <- tryCatch(jsonlite::fromJSON(paste(suppressWarnings(system(var_cmd, intern = TRUE)), collapse = "")), error = function(e) NULL)
  if (is.null(var_js) || length(var_js$data) == 0 || !("variantId" %in% names(var_js$data))) next
  vid <- var_js$data$variantId[1]
  alt <- var_js$data$alt[1]
  ref <- var_js$data$ref[1]
  
  dyn_cmd <- sprintf('curl -4 -s --connect-timeout 6 --max-time 10 "https://gtexportal.org/api/v2/association/dyneqtl?tissueSiteDetailId=Breast_Mammary_Tissue&gencodeId=%s&variantId=%s"', gid, vid)
  dyn_js  <- tryCatch(jsonlite::fromJSON(paste(suppressWarnings(system(dyn_cmd, intern = TRUE)), collapse = "")), error = function(e) NULL)
  if (!is.null(dyn_js) && !is.null(dyn_js$nes)) {
    rows[[paste(sym, rs)]] <- data.table(
      gene_symbol = sym,
      gencodeId = gid,
      rsid = rs,
      variantId_b38 = vid,
      gtex_alt = alt,
      gtex_ref = ref,
      gtex_nes = as.numeric(dyn_js$nes),
      gtex_se = as.numeric(dyn_js$error),
      gtex_pval = as.numeric(dyn_js$pValue),
      gtex_t = as.numeric(dyn_js$tStatistic),
      eqtlgen_beta = pos_snps$beta.exposure[i],
      eqtlgen_pval = pos_snps$pval.exposure[i]
    )
  } else {
    rows[[paste(sym, rs)]] <- data.table(
      gene_symbol = sym, gencodeId = gid, rsid = rs, variantId_b38 = vid,
      gtex_alt = alt, gtex_ref = ref, gtex_nes = NA_real_, gtex_se = NA_real_,
      gtex_pval = NA_real_, gtex_t = NA_real_,
      eqtlgen_beta = pos_snps$beta.exposure[i], eqtlgen_pval = pos_snps$pval.exposure[i]
    )
  }
}

dt <- rbindlist(rows, fill = TRUE)
fwrite(dt, file.path(res_dir, "module3_cis_positive_genes_gtex_breast_dyneqtl.csv"))
cat("=== GTEx v8 Breast_Mammary_Tissue dyneqtl 核查结果 ===\n")
print(dt[, .(gene_symbol, rsid, gtex_nes = round(gtex_nes, 3), gtex_pval = signif(gtex_pval, 3),
             eqtlgen_beta = round(eqtlgen_beta, 3), eqtlgen_pval = signif(eqtlgen_pval, 3))])
