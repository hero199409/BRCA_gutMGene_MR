# 17b_tisch2_tumor_scrna_analysis.R
# 升级模块 1A: 真实乳腺癌肿瘤单细胞转录组 (TISCH2: BRCA_EMTAB8107 [33k cells] + BRCA_GSE176078 [89k cells])
# 验证 PPARG / CXCL8 / IL1B / FFAR2 细胞定位与脂肪标记物 (ADIPOQ / PLIN1 / FABP4 / LEP) 在肿瘤/基质细胞中的缺失与特异性

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
  library(Matrix)
})

sc_dir  <- "00_data/scRNA_cache"
res_dir <- "02_results"

target_genes <- c(
  "PPARG", "CXCL8", "IL1B", "FFAR2",
  "ADIPOQ", "PLIN1", "FABP4", "LEP",
  "EPCAM", "KRT19", "PTPRC", "CD68", "S100A8", "CSF3R", "VWF", "COL1A1"
)

# Helper: Read only target gene rows from TISCH2 H5 cleanly under 16GB RAM limit
read_tisch_h5_targets <- function(h5_path, genes) {
  ls_df <- rhdf5::h5ls(h5_path)
  cat("  H5 datasets in", basename(h5_path), ":\n")
  print(ls_df[ls_df$otype == "H5I_DATASET", c("group", "name", "dim")])

  find_path <- function(patterns) {
    fullp <- ifelse(ls_df$group == "/", paste0("/", ls_df$name), paste0(ls_df$group, "/", ls_df$name))
    fullp <- gsub("//+", "/", fullp)
    ds <- fullp[ls_df$otype == "H5I_DATASET"]
    for (pt in patterns) {
      hit <- grep(pt, ds, ignore.case = TRUE, value = TRUE)
      if (length(hit) > 0) return(hit[1])
    }
    NA_character_
  }

  p_barcodes <- find_path(c("/barcodes$"))
  p_genes    <- find_path(c("/features/name$", "/gene_names$", "/genes$"))
  p_data     <- find_path(c("/data$"))
  p_indices  <- find_path(c("/indices$"))
  p_indptr   <- find_path(c("/indptr$"))
  p_shape    <- find_path(c("/shape$"))

  barcodes <- as.character(rhdf5::h5read(h5_path, p_barcodes))
  all_genes <- as.character(rhdf5::h5read(h5_path, p_genes))
  cat(sprintf("  barcodes=%d | genes=%d\n", length(barcodes), length(all_genes)))

  idx0 <- match(genes, all_genes) - 1L
  valid <- !is.na(idx0)
  found_genes <- genes[valid]
  idx0_valid <- idx0[valid]
  cat("  Matched target genes:", paste(found_genes, collapse = ", "), "\n")
  if (!all(valid)) cat("  Missing genes:", paste(genes[!valid], collapse = ", "), "\n")

  indptr  <- as.numeric(rhdf5::h5read(h5_path, p_indptr))
  indices <- as.integer(rhdf5::h5read(h5_path, p_indices))
  vals    <- as.numeric(rhdf5::h5read(h5_path, p_data))
  cat(sprintf("  Sparse entries=%d | expr range=[%.3f, %.3f]\n", length(vals), min(vals), max(vals)))

  # Build dense (targets x cells) matrix selectively
  mat <- matrix(0, nrow = length(found_genes), ncol = length(barcodes),
                dimnames = list(found_genes, barcodes))

  if (length(indptr) == length(barcodes) + 1L) {
    # CSC: columns = barcodes, indices = 0-based gene indices
    keep <- which(indices %in% idx0_valid)
    if (length(keep) > 0) {
      r_idx <- match(indices[keep], idx0_valid)
      c_idx <- findInterval(keep - 1, indptr)
      ok <- r_idx >= 1L & r_idx <= nrow(mat) & c_idx >= 1L & c_idx <= ncol(mat)
      mat[cbind(r_idx[ok], c_idx[ok])] <- vals[keep[ok]]
    }
  } else if (length(indptr) == length(all_genes) + 1L) {
    # CSC transposed / CSR: outer dimension = genes, indices = 0-based cell indices
    for (gi in seq_along(idx0_valid)) {
      g1 <- idx0_valid[gi] + 1L
      s <- indptr[g1] + 1
      e <- indptr[g1 + 1]
      if (e >= s) {
        c_idx <- indices[s:e] + 1L
        mat[gi, c_idx] <- vals[s:e]
      }
    }
  } else {
    stop("Unexpected indptr length: ", length(indptr))
  }

  rm(indptr, indices, vals); gc()
  rhdf5::h5closeAll()
  mat
}

# Harmonized Wu et al. 2021 (GSE176078) loader using pre-extracted 16-gene sparse TSV + full 100,064-cell metadata + TISCH2 UMAP
process_gse176078_wu <- function() {
  ds_name   <- "BRCA_GSE176078"
  cnt_file  <- file.path(sc_dir, "GSE176078_target16_counts.tsv")
  wu_meta_f <- file.path(sc_dir, "GSE176078_raw", "Wu_etal_2021_BRCA_scRNASeq", "metadata.csv")
  ti_meta_f <- file.path(sc_dir, "BRCA_GSE176078_CellMetainfo_table.tsv")
  if (!file.exists(cnt_file) || !file.exists(wu_meta_f)) {
    cat("Skipping missing GSE176078 files.\n")
    return(NULL)
  }
  cat("\n====================================================\n")
  cat("Processing scRNA dataset (Wu et al. 2021 full): ", ds_name, "\n")
  cat("====================================================\n")

  cnt   <- fread(cnt_file)
  m_wu  <- fread(wu_meta_f)
  setnames(m_wu, "V1", "Cell")

  # Attach TISCH2 UMAP coordinates where available (88,978 overlapping cells)
  if (file.exists(ti_meta_f)) {
    m_ti <- fread(ti_meta_f, select = c("Cell", "UMAP_1", "UMAP_2"))
    m_wu <- merge(m_wu, m_ti, by = "Cell", all.x = TRUE, sort = FALSE)
  } else {
    m_wu[, `:=`(UMAP_1 = NA_real_, UMAP_2 = NA_real_)]
  }

  # Harmonize Wu et al. author-curated celltype_major to consistent lineage nomenclature
  lineage_map <- c(
    "Cancer Epithelial" = "Malignant",
    "Normal Epithelial" = "Epithelial",
    "Myeloid"           = "Mono/Macro",
    "Endothelial"       = "Endothelial",
    "CAFs"              = "Fibroblasts",
    "PVL"               = "PVL",
    "T-cells"           = "T_cells",
    "B-cells"           = "B",
    "Plasmablasts"      = "Plasma"
  )
  m_wu[, Celltype_major := unname(lineage_map[celltype_major])]
  m_wu[, Malignancy := fifelse(Celltype_major == "Malignant", "Malignant cells",
                       fifelse(Celltype_major %in% c("Mono/Macro", "T_cells", "B", "Plasma"), "Immune cells", "Stromal/Normal"))]

  cat(sprintf("  Wu full QC cells=%d | with TISCH2 UMAP=%d\n",
              nrow(m_wu), sum(!is.na(m_wu$UMAP_1) & !is.na(m_wu$UMAP_2))))

  dt_plot <- m_wu[, .(
    Dataset        = ds_name,
    Cell           = Cell,
    UMAP_1         = as.numeric(UMAP_1),
    UMAP_2         = as.numeric(UMAP_2),
    Malignancy     = Malignancy,
    Celltype_major = Celltype_major,
    Celltype_minor = as.character(celltype_minor),
    Patient        = as.character(orig.ident),
    Subtype        = as.character(subtype)
  )]

  norm_fac <- 1e4 / as.numeric(m_wu$nCount_RNA)
  for (g in target_genes) {
    sub_g <- cnt[gene == g]
    idx <- match(dt_plot$Cell, sub_g$barcode)
    raw_c <- ifelse(is.na(idx), 0, as.numeric(sub_g$count[idx]))
    dt_plot[[g]] <- log1p(raw_c * norm_fac)
  }

  sum_rows <- list()
  for (g in target_genes) {
    s <- dt_plot[, .(
      dataset = ds_name,
      gene = g,
      n_cells = .N,
      n_expr = sum(get(g) > 0),
      pct_expr = round(100 * mean(get(g) > 0), 2),
      mean_expr = round(mean(get(g)), 4),
      pos_mean_expr = round(ifelse(sum(get(g) > 0) > 0, mean(get(g)[get(g) > 0]), 0), 4),
      sd_expr = round(sd(get(g)), 4)
    ), by = .(celltype_major = Celltype_major)]
    sum_rows[[length(sum_rows) + 1]] <- s
  }
  sum_dt <- rbindlist(sum_rows)

  cat("\n--- Cell-type summary (key genes) in", ds_name, "---\n")
  print(sum_dt[gene %in% c("PPARG", "CXCL8", "ADIPOQ", "PLIN1", "IL1B", "FFAR2")][order(gene, -mean_expr)])

  sub_rows <- list()
  sub_dt <- dt_plot[Subtype %in% c("ER+", "HER2+", "TNBC")]
  for (ct in c("ALL_CELLS", "Malignant", "Mono/Macro", "Endothelial", "Fibroblasts", "Epithelial", "PVL", "T_cells")) {
    df_ct <- if (ct == "ALL_CELLS") sub_dt else sub_dt[Celltype_major == ct]
    if (nrow(df_ct) < 50) next
    for (g in c("PPARG", "CXCL8", "IL1B", "FFAR2", "ADIPOQ", "PLIN1")) {
      if (!(g %in% names(df_ct))) next
      kw_p <- tryCatch(kruskal.test(get(g) ~ factor(Subtype), data = df_ct)$p.value, error = function(e) NA_real_)
      er_v <- df_ct[Subtype == "ER+", get(g)]
      tn_v <- df_ct[Subtype == "TNBC", get(g)]
      w_p  <- tryCatch(wilcox.test(er_v, tn_v)$p.value, error = function(e) NA_real_)
      by_st <- df_ct[, .(n = .N, pct = 100 * mean(get(g) > 0), mean = mean(get(g))), by = Subtype]
      for (k in seq_len(nrow(by_st))) {
        sub_rows[[length(sub_rows) + 1]] <- data.table(
          dataset = ds_name, celltype = ct, gene = g, subtype = by_st$Subtype[k],
          n_cells = by_st$n[k], pct_expr = round(by_st$pct[k], 2),
          mean_expr = round(by_st$mean[k], 4), kw_pval = signif(kw_p, 4),
          er_vs_tnbc_wilcox_p = signif(w_p, 4)
        )
      }
    }
  }
  sub_res <- rbindlist(sub_rows)
  cat("\n--- Subtype comparison (ER+ vs HER2+ vs TNBC) in", ds_name, "---\n")
  print(sub_res[celltype %in% c("ALL_CELLS", "Malignant", "Mono/Macro", "Endothelial") & gene %in% c("PPARG", "CXCL8")])

  rm(cnt, m_wu); gc()
  list(plot_dt = dt_plot, celltype_sum = sum_dt, subtype_sum = sub_res)
}

process_dataset <- function(ds_name) {
  if (ds_name == "BRCA_GSE176078" && file.exists(file.path(sc_dir, "GSE176078_target16_counts.tsv"))) {
    return(process_gse176078_wu())
  }
  h5_file   <- file.path(sc_dir, paste0(ds_name, "_expression.h5"))
  meta_file <- file.path(sc_dir, paste0(ds_name, "_CellMetainfo_table.tsv"))
  if (!file.exists(h5_file) || !file.exists(meta_file)) {
    cat("Skipping missing dataset:", ds_name, "\n")
    return(NULL)
  }
  cat("\n====================================================\n")
  cat("Processing scRNA dataset:", ds_name, "\n")
  cat("====================================================\n")
  meta <- fread(meta_file)
  mat  <- read_tisch_h5_targets(h5_file, target_genes)

  common <- intersect(meta$Cell, colnames(mat))
  cat("  Matched cells meta vs H5:", length(common), "/", nrow(meta), "\n")
  if (length(common) == 0) {
    # Fallback strip suffix after - if needed
    colnames(mat) <- sub("-[0-9]+$", "", colnames(mat))
    meta[, Cell_clean := sub("-[0-9]+$", "", Cell)]
    common <- intersect(meta$Cell_clean, colnames(mat))
    meta <- meta[match(common, Cell_clean)]
  } else {
    meta <- meta[match(common, Cell)]
  }
  mat <- mat[, common, drop = FALSE]

  # Build long-ready / merged data.table for our key genes
  dt_plot <- meta[, .(
    Dataset = ds_name,
    Cell = Cell,
    UMAP_1 = as.numeric(UMAP_1),
    UMAP_2 = as.numeric(UMAP_2),
    Malignancy = `Celltype (malignancy)`,
    Celltype_major = `Celltype (major-lineage)`,
    Celltype_minor = `Celltype (minor-lineage)`,
    Patient = Patient,
    Subtype = if ("Subtype" %in% names(meta)) as.character(Subtype) else NA_character_
  )]
  for (g in rownames(mat)) {
    dt_plot[[g]] <- as.numeric(mat[g, ])
  }

  # Summarize across Celltype (major-lineage)
  sum_rows <- list()
  for (g in rownames(mat)) {
    s <- dt_plot[, .(
      dataset = ds_name,
      gene = g,
      n_cells = .N,
      n_expr = sum(get(g) > 0),
      pct_expr = round(100 * mean(get(g) > 0), 2),
      mean_expr = round(mean(get(g)), 4),
      pos_mean_expr = round(ifelse(sum(get(g) > 0) > 0, mean(get(g)[get(g) > 0]), 0), 4),
      sd_expr = round(sd(get(g)), 4)
    ), by = .(celltype_major = Celltype_major)]
    sum_rows[[length(sum_rows) + 1]] <- s
  }
  sum_dt <- rbindlist(sum_rows)

  # Print highlights for PPARG, CXCL8, IL1B, FFAR2, ADIPOQ, PLIN1
  cat("\n--- Cell-type summary (key genes) in", ds_name, "---\n")
  print(sum_dt[gene %in% c("PPARG", "CXCL8", "ADIPOQ", "PLIN1", "IL1B", "FFAR2")][order(gene, -mean_expr)])

  # Subtype analysis if available (GSE176078 has ER+ / HER2+ / TNBC)
  sub_rows <- list()
  if (any(!is.na(dt_plot$Subtype) & dt_plot$Subtype %in% c("ER+", "HER2+", "TNBC"))) {
    sub_dt <- dt_plot[Subtype %in% c("ER+", "HER2+", "TNBC")]
    for (ct in c("ALL_CELLS", "Malignant", "Mono/Macro", "Endothelial", "Fibroblasts", "CD8Tex")) {
      df_ct <- if (ct == "ALL_CELLS") sub_dt else sub_dt[Celltype_major == ct]
      if (nrow(df_ct) < 50) next
      for (g in c("PPARG", "CXCL8", "IL1B", "FFAR2", "ADIPOQ", "PLIN1")) {
        if (!(g %in% names(df_ct))) next
        kw_p <- tryCatch(kruskal.test(get(g) ~ factor(Subtype), data = df_ct)$p.value, error = function(e) NA_real_)
        er_v <- df_ct[Subtype == "ER+", get(g)]
        tn_v <- df_ct[Subtype == "TNBC", get(g)]
        w_p  <- tryCatch(wilcox.test(er_v, tn_v)$p.value, error = function(e) NA_real_)
        by_st <- df_ct[, .(n = .N, pct = 100 * mean(get(g) > 0), mean = mean(get(g))), by = Subtype]
        for (k in seq_len(nrow(by_st))) {
          sub_rows[[length(sub_rows) + 1]] <- data.table(
            dataset = ds_name, celltype = ct, gene = g, subtype = by_st$Subtype[k],
            n_cells = by_st$n[k], pct_expr = round(by_st$pct[k], 2),
            mean_expr = round(by_st$mean[k], 4), kw_pval = signif(kw_p, 4),
            er_vs_tnbc_wilcox_p = signif(w_p, 4)
          )
        }
      }
    }
  }
  sub_res <- if (length(sub_rows) > 0) rbindlist(sub_rows) else data.table()
  if (nrow(sub_res) > 0) {
    cat("\n--- Subtype comparison (ER+ vs HER2+ vs TNBC) in", ds_name, "---\n")
    print(sub_res[celltype %in% c("Malignant", "Mono/Macro", "Endothelial") & gene %in% c("PPARG", "CXCL8")])
  }

  rm(mat, meta); gc()
  list(plot_dt = dt_plot, celltype_sum = sum_dt, subtype_sum = sub_res)
}

datasets <- c("BRCA_EMTAB8107", "BRCA_GSE148673", "BRCA_GSE176078")
all_ct <- list()
all_st <- list()
plot_caches <- list()

for (ds in datasets) {
  out <- process_dataset(ds)
  if (!is.null(out)) {
    all_ct[[ds]] <- out$celltype_sum
    if (nrow(out$subtype_sum) > 0) all_st[[ds]] <- out$subtype_sum
    plot_caches[[ds]] <- out$plot_dt
  }
}

if (length(all_ct) > 0) {
  ct_dt <- rbindlist(all_ct)
  fwrite(ct_dt, file.path(res_dir, "planA_scrna_celltype_summary.csv"))
}
if (length(all_st) > 0) {
  st_dt <- rbindlist(all_st)
  fwrite(st_dt, file.path(res_dir, "planA_scrna_subtype_comparison.csv"))
}
if (length(plot_caches) > 0) {
  saveRDS(plot_caches, file.path(sc_dir, "planA_scrna_plot_cache.rds"))
  cat("\nSaved scRNA summary CSVs and plot cache cleanly!\n")
}
