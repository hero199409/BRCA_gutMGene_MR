# 10b_check_gtex_breast_single_eqtl.R
# 对 8 个纯 cis-MR 阳性基因的精准 eQTLGen SNP + GTEx Breast_Mammary_Tissue 独立 lead SNP 查证
source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(TwoSampleMR)
})

res_dir   <- "02_results"
cache_dir <- "00_data/mr_cache"

sig_cis <- fread(file.path(res_dir, "module3_eqtlgen_cis_only_positive_p005.csv"))
pos_genes <- unique(sig_cis$gene_symbol)
cis_audit <- fread(file.path(res_dir, "module3_all41_genes_iv_cis_audit.csv"))
pos_snps <- cis_audit[gene_symbol %in% pos_genes & is_cis_1mb == TRUE]

# 1. 查询每个基因在 GTEx Breast_Mammary_Tissue 中针对这些具体 SNP 的表达关联 (singleTissueEqtl)
gtex_snp_rows <- list()
for (i in seq_len(nrow(pos_snps))) {
  sym <- pos_snps$gene_symbol[i]
  rs  <- pos_snps$SNP[i]
  ref_cmd <- sprintf('curl -4 -s --connect-timeout 6 --max-time 10 "https://gtexportal.org/api/v2/reference/gene?geneId=%s&gencodeVersion=v26"', sym)
  ref_js  <- tryCatch(jsonlite::fromJSON(paste(suppressWarnings(system(ref_cmd, intern = TRUE)), collapse = "")), error = function(e) NULL)
  gid <- if (!is.null(ref_js) && length(ref_js$data) > 0 && "gencodeId" %in% names(ref_js$data)) ref_js$data$gencodeId[1] else NA
  if (is.na(gid)) next
  
  # 查 variantId 映射 via rsid
  var_cmd <- sprintf('curl -4 -s --connect-timeout 6 --max-time 10 "https://gtexportal.org/api/v2/dataset/variant?snpId=%s&datasetId=gtex_v8"', rs)
  var_js  <- tryCatch(jsonlite::fromJSON(paste(suppressWarnings(system(var_cmd, intern = TRUE)), collapse = "")), error = function(e) NULL)
  vid <- if (!is.null(var_js) && length(var_js$data) > 0 && "variantId" %in% names(var_js$data)) var_js$data$variantId[1] else NA
  if (is.na(vid)) next
  
  sq_cmd <- sprintf('curl -4 -s --connect-timeout 6 --max-time 10 "https://gtexportal.org/api/v2/association/singleTissueEqtl?tissueSiteDetailId=Breast_Mammary_Tissue&gencodeId=%s&variantId=%s&datasetId=gtex_v8"', gid, vid)
  sq_js  <- tryCatch(jsonlite::fromJSON(paste(suppressWarnings(system(sq_cmd, intern = TRUE)), collapse = "")), error = function(e) NULL)
  if (!is.null(sq_js) && length(sq_js$data) > 0 && is.data.frame(sq_js$data)) {
    d <- as.data.table(sq_js$data)
    d[, gene_symbol := sym]
    d[, rsid := rs]
    d[, eqtlgen_beta := pos_snps$beta.exposure[i]]
    gtex_snp_rows[[paste(sym, rs)]] <- d
    cat(sprintf("GTEx Breast %s (%s): nes=%.3f, p=%.3g | eQTLGen beta=%.3f\n", sym, rs, d$nes[1], d$pValue[1], pos_snps$beta.exposure[i]))
  } else {
    cat(sprintf("GTEx Breast %s (%s): no singleTissueEqtl record\n", sym, rs))
  }
}

if (length(gtex_snp_rows) > 0) {
  out_dt <- rbindlist(gtex_snp_rows, fill = TRUE)
  fwrite(out_dt, file.path(res_dir, "module3_cis_positive_genes_gtex_breast_single_snp.csv"))
}

# 2. 检验 FFAR2 在 GTEx Breast_Mammary_Tissue 的独立 lead eQTL (rs73927976) 本身在 BCAC 中的 MR 效应
gtex_ffar2 <- fread(file.path(res_dir, "module3_cis_positive_genes_gtex_breast_eqtl.csv"))
if ("snpId" %in% names(gtex_ffar2) && nrow(gtex_ffar2) > 0) {
  # 获取 rs73927976 的等位基因信息
  var_cmd <- 'curl -4 -s --connect-timeout 6 --max-time 10 "https://gtexportal.org/api/v2/dataset/variant?snpId=rs73927976&datasetId=gtex_v8"'
  vjs <- jsonlite::fromJSON(paste(suppressWarnings(system(var_cmd, intern = TRUE)), collapse = ""))
  vdf <- as.data.table(vjs$data)
  if (nrow(vdf) > 0) {
    exp_ffar2 <- data.frame(
      SNP = "rs73927976",
      beta.exposure = gtex_ffar2$nes[1],
      se.exposure = abs(gtex_ffar2$nes[1]) / qnorm(1 - gtex_ffar2$pValue[1] / 2),
      effect_allele.exposure = vdf$alt[1],
      other_allele.exposure = vdf$ref[1],
      eaf.exposure = vdf$maf[1],
      pval.exposure = gtex_ffar2$pValue[1],
      exposure = "FFAR2_GTEx_Breast",
      id.exposure = "GTEx_Breast_FFAR2"
    )
    out_bcac_f <- TwoSampleMR::extract_outcome_data(snps = "rs73927976", outcomes = c("ieu-a-1126","ieu-a-1127","ieu-a-1128"), proxies = TRUE)
    if (!is.null(out_bcac_f) && nrow(out_bcac_f) > 0) {
      dh <- TwoSampleMR::harmonise_data(exp_ffar2, out_bcac_f, action = 2)
      mr_f <- as.data.table(TwoSampleMR::mr(dh, method_list = "mr_wald_ratio"))
      mr_f[, OR := exp(b)]
      fwrite(mr_f, file.path(res_dir, "module3_ffar2_gtex_breast_mr.csv"))
      cat("\n=== FFAR2 GTEx Breast_Mammary_Tissue lead eQTL (rs73927976) -> BCAC MR ===\n")
      print(mr_f[, .(id.outcome, nsnp, OR = round(OR, 3), pval = signif(pval, 3))])
    }
  }
}
