# 25_1_prepare_cns_plot_cache.R
# 为 CNS 级图形升级准备轻量缓存（仅约数百 KB，不修改任何已有统计结果）：
# 1) Fig 3D 区域位点 GRCh37 蛋白编码基因/外显子轨 + 1000G EUR LD (r2) 信息
# 2) Fig 7B 核心通路（CXCL8 ER- & PPARG ER+）原始 clusterProfiler GSEA 运行富集曲线对象

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ieugwasr)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

cache_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/mr_cache")
bulk_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
res_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

curl_json <- function(url, max_try = 3L, timeout_s = 25L) {
  for (k in seq_len(max_try)) {
    cmd <- sprintf('curl.exe -4 -s --connect-timeout 8 --max-time %d "%s"', timeout_s, url)
    raw <- tryCatch(paste(suppressWarnings(system(cmd, intern = TRUE)), collapse = ""), error = function(e) "")
    if (nzchar(raw) && substr(trimws(raw), 1, 1) %in% c("{", "[")) {
      parsed <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = TRUE), error = function(e) NULL)
      if (!is.null(parsed) && is.null(parsed$error)) return(parsed)
    }
    Sys.sleep(1.0)
  }
  NULL
}

# ------------------------------------------------------------------------------
# Part 1: Fig 3D 4 个位点 GRCh37 基因/外显子注释 + 1000G EUR LD (r2)
# ------------------------------------------------------------------------------
cat("=== [1/2] Fetching GRCh37 gene/exon tracks & 1000G EUR LD for 4 loci ===\n")
reg_snps <- readRDS(file.path(cache_dir, "planA_coloc_plot_snps.rds"))
lead_sum <- fread(file.path(res_dir, "planA_coloc_lead_snps_summary.csv"))

loci_spec <- list(
  PPARG_ER_pos  = list(chr = "3", start = 12004290L,  end = 13004290L,  target = "PPARG"),
  IL1B_ER_neg   = list(chr = "2", start = 112964709L, end = 113964709L, target = "IL1B"),
  CXCL8_ER_neg  = list(chr = "4", start = 74072039L,  end = 75072039L,  target = "CXCL8"),
  CXCR2_Overall = list(chr = "2", start = 218491005L, end = 219491005L, target = "CXCR2")
)

genes_list <- list()
exons_list <- list()
ld_list    <- list()

for (pk in names(loci_spec)) {
  sp <- loci_spec[[pk]]
  reg_str <- sprintf("%s:%d-%d", sp$chr, sp$start, sp$end)
  cat("  Processing locus:", pk, "(", reg_str, ")...\n")

  # 1a. Protein-coding genes in GRCh37 region
  g_url <- sprintf("https://grch37.rest.ensembl.org/overlap/region/human/%s?feature=gene;content-type=application/json", reg_str)
  gdf <- curl_json(g_url)
  if (is.data.frame(gdf) && nrow(gdf) > 0) {
    gdt <- as.data.table(gdf)[biotype == "protein_coding"]
    if (!("external_name" %in% names(gdt))) gdt[, external_name := id]
    gdt[external_name == "IL8", external_name := "CXCL8"]
    gdt[, `:=`(pair = pk, is_target = (external_name == sp$target))]
    gdt <- gdt[!is.na(external_name) & nzchar(external_name)]
    genes_list[[pk]] <- gdt[, .(pair, id, external_name, chr = as.character(seq_region_name), start, end, strand, is_target)]

    # Expand canonical exons for target gene + up to 5 closest protein-coding genes
    mid_pos <- (sp$start + sp$end) / 2
    gdt[, dist_mid := abs((start + end) / 2 - mid_pos)]
    g_pick <- unique(c(gdt[is_target == TRUE]$id, gdt[order(dist_mid)][1:min(6, .N)]$id))
    for (ensg in g_pick) {
      lk <- curl_json(sprintf("https://grch37.rest.ensembl.org/lookup/id/%s?expand=1;content-type=application/json", ensg))
      if (!is.null(lk) && is.data.frame(lk$Transcript) && nrow(lk$Transcript) > 0) {
        tr_idx <- which(lk$Transcript$is_canonical == 1)[1]
        if (is.na(tr_idx)) tr_idx <- 1L
        ex <- lk$Transcript$Exon[[tr_idx]]
        if (is.data.frame(ex) && nrow(ex) > 0) {
          ex_dt <- as.data.table(ex)[, .(pair = pk, gene_id = ensg, gene_name = ifelse(identical(lk$display_name, "IL8"), "CXCL8", lk$display_name),
                                         start = start, end = end, strand = strand)]
          exons_list[[paste(pk, ensg)]] <- ex_dt
        }
      }
    }
  }

  # 1b. 1000G EUR LD (r2) relative to lead_eqtl_snp & lead_brca_snp
  ls_row <- lead_sum[pair == pk]
  sub_snps <- reg_snps[pair == pk]
  lead_e <- ls_row$lead_eqtl_snp[1]
  lead_o <- ls_row$lead_brca_snp[1]

  fetch_ld_vs_lead <- function(lead_rs) {
    # Fast Ensembl 1000G Phase 3 EUR LD endpoint (500 kb window)
    ld_u <- sprintf("https://grch37.rest.ensembl.org/ld/human/%s/1000GENOMES:phase_3:EUR?window_size=500;r2=0.05;content-type=application/json", lead_rs)
    ld_df <- curl_json(ld_u)
    out_dt <- data.table(rsid = lead_rs, r2 = 1.0)
    if (is.data.frame(ld_df) && nrow(ld_df) > 0) {
      tmp <- as.data.table(ld_df)[, .(rsid = as.character(variation2), r2 = as.numeric(r2))]
      out_dt <- rbind(out_dt, tmp)
    }
    # Supplement top 180 SNPs by min(p_e, p_o) via ieugwasr::ld_matrix if missing
    top_miss <- setdiff(sub_snps[order(pmin(p_e, p_o))][1:min(180, .N)]$rsid, out_dt$rsid)
    if (length(top_miss) > 0) {
      q_vec <- unique(c(lead_rs, top_miss))
      mat <- tryCatch(suppressWarnings(ieugwasr::ld_matrix(q_vec, pop = "EUR")), error = function(e) NULL)
      if (is.matrix(mat) && nrow(mat) > 1) {
        rn_clean <- sub("_.*$", "", rownames(mat))
        hit_i <- which(rn_clean == lead_rs)[1]
        if (!is.na(hit_i)) {
          r2_vals <- (mat[hit_i, ])^2
          add_dt <- data.table(rsid = rn_clean, r2 = as.numeric(r2_vals))
          out_dt <- rbind(out_dt, add_dt)
        }
      }
    }
    out_dt[, .(r2 = max(r2, na.rm = TRUE)), by = rsid]
  }

  ld_e <- fetch_ld_vs_lead(lead_e)
  setnames(ld_e, "r2", "r2_eqtl_lead")
  ld_o <- fetch_ld_vs_lead(lead_o)
  setnames(ld_o, "r2", "r2_brca_lead")
  ld_mg <- merge(ld_e, ld_o, by = "rsid", all = TRUE)
  ld_mg[, pair := pk]
  ld_list[[pk]] <- ld_mg
  cat("    LD hits vs eQTL lead:", sum(!is.na(ld_mg$r2_eqtl_lead)), "| vs BRCA lead:", sum(!is.na(ld_mg$r2_brca_lead)), "\n")
}

loci_cache <- list(
  genes = rbindlist(genes_list, fill = TRUE),
  exons = rbindlist(exons_list, fill = TRUE),
  ld    = rbindlist(ld_list, fill = TRUE)
)
saveRDS(loci_cache, file.path(cache_dir, "planA_cns_loci_anno_cache.rds"))
cat("Saved:", file.path(cache_dir, "planA_cns_loci_anno_cache.rds"), "\n")

# ------------------------------------------------------------------------------
# Part 2: 抽取 Fig 7B 经典 GSEA Running Score 曲线对象 (CXCL8 ER- & PPARG ER+)
# ------------------------------------------------------------------------------
cat("\n=== [2/2] Extracting GSEA gseaResult objects for Fig 7B ===\n")
tc <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
expr_clean <- tc$tumor_log2
expr_clean[is.na(expr_clean)] <- 0
clin <- tc$clin[match(colnames(expr_clean), tc$clin$SAMPLE_ID)]
rm(tc); gc()

gmt <- clusterProfiler::read.gmt(file.path(bulk_dir, "hallmark_human_entrez.gmt"))
sym_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = rownames(expr_clean),
                                                  columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL"))
sym_map <- as.data.table(sym_map)[!is.na(ENTREZID) & !duplicated(SYMBOL)]

extract_gsea_obj <- function(gene, idx_cols) {
  sub_mat <- expr_clean[, idx_cols, drop = FALSE]
  gv <- as.numeric(sub_mat[gene, ])
  grp <- factor(ifelse(gv >= median(gv, na.rm = TRUE), "High", "Low"), levels = c("Low", "High"))
  design <- model.matrix(~ grp)
  fit <- limma::eBayes(limma::lmFit(sub_mat, design))
  tt  <- limma::topTable(fit, coef = "grpHigh", number = Inf, sort.by = "none")
  tt$SYMBOL <- rownames(tt)
  tt_dt <- merge(as.data.table(tt), sym_map, by = "SYMBOL")[!is.na(logFC) & !duplicated(ENTREZID)][order(-logFC, -t)]
  gene_list <- setNames(tt_dt$logFC, tt_dt$ENTREZID)
  set.seed(2026)
  clusterProfiler::GSEA(
    geneList = gene_list, TERM2GENE = gmt,
    minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1.0, pAdjustMethod = "BH", verbose = FALSE
  )
}

gsea_objs <- list(
  CXCL8_ERneg = extract_gsea_obj("CXCL8", which(clin$ER_Xena == "Negative")),
  PPARG_ERpos = extract_gsea_obj("PPARG", which(clin$ER_Xena == "Positive"))
)
saveRDS(gsea_objs, file.path(bulk_dir, "planA_cns_gsea_objects.rds"))
cat("Saved:", file.path(bulk_dir, "planA_cns_gsea_objects.rds"), "\n")
cat("PREP CACHE COMPLETE!\n")
