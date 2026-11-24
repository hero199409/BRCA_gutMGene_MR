# 17b_tisch2_tumor_scrna_analysis.R
# 乳腺癌单细胞转录组：Wu 原始注释、TISCH2 注释和本地敏感性分析
# 验证 PPARG / CXCL8 / IL1B / FFAR2 细胞定位与脂肪标记物 (ADIPOQ / PLIN1 / FABP4 / LEP) 在肿瘤/基质细胞中的缺失与特异性

# 仅使用本地数据。
suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
})

sc_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/scRNA_cache")
res_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

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

# Wu 表达以供体为推断单位。每个供体内用所有该细胞群细胞（含零值），
# 先求 log1p(CP10K) 的均值和非零检出百分比，然后比较供体，供体等权。
wu_subtypes <- c("ER+", "HER2+", "TNBC")
subtype_genes <- c("PPARG", "CXCL8", "IL1B", "FFAR2", "ADIPOQ", "PLIN1")

wu_donor_subtype <- function(cells) {
  stopifnot(all(cells$Subtype %in% wu_subtypes), all(nzchar(cells$Patient)))
  donors <- unique(cells[, .(donor = Patient, subtype = Subtype)])
  stopifnot(!anyDuplicated(donors$donor))
  panels <- rbindlist(list(
    cells[, c("Patient", "Subtype", "Celltype_major", subtype_genes), with = FALSE],
    cells[, c("Patient", "Subtype", subtype_genes), with = FALSE][, Celltype_major := "ALL_CELLS"]
  ), use.names = TRUE)
  panels[Celltype_major == "Mono/Macro", Celltype_major := "Myeloid"]
  n_cells <- panels[, .(n_cells = .N), by = .(donor = Patient, subtype = Subtype, celltype = Celltype_major)]
  lineage_grid <- merge(donors[, .(lineage_key = 1L, donor, subtype)],
                        data.table(lineage_key = 1L, celltype = sort(unique(panels$Celltype_major))),
                        by = "lineage_key", allow.cartesian = TRUE)[, lineage_key := NULL]
  n_cells <- merge(lineage_grid, n_cells, by = c("donor", "subtype", "celltype"), all.x = TRUE)
  n_cells[is.na(n_cells), n_cells := 0L]
  fwrite(n_cells, file.path(res_dir, "planA_scrna_wu_donor_cell_counts.csv"))

  donor_expr <- melt(panels, id.vars = c("Patient", "Subtype", "Celltype_major"),
                     measure.vars = subtype_genes, variable.name = "gene", value.name = "expr",
                     variable.factor = FALSE)[,
    .(n_cells = .N, mean_expr = mean(expr), pct_detected = 100 * mean(expr > 0)),
    by = .(celltype = Celltype_major, donor = Patient, subtype = Subtype, gene)]
  fwrite(donor_expr, file.path(res_dir, "planA_scrna_wu_donor_expression.csv"))

  summarize_groups <- function(threshold) {
    n_cells[, {
      kept <- n_cells >= threshold
      nc <- n_cells[kept]
      .(min_cells = threshold, n_donors_total = .N, n_donors = sum(kept), n_cells_retained = sum(nc),
        min_cells_per_donor = if (length(nc)) min(nc) else NA_integer_,
        median_cells_per_donor = if (length(nc)) as.numeric(median(nc)) else NA_real_,
        max_cells_per_donor = if (length(nc)) max(nc) else NA_integer_)
    }, by = .(celltype, subtype)]
  }
  group_counts <- rbindlist(lapply(c(1L, 10L), summarize_groups))
  fwrite(group_counts, file.path(res_dir, "planA_scrna_wu_donor_groups.csv"))

  donor_long <- melt(donor_expr, measure.vars = c("mean_expr", "pct_detected"),
                     variable.name = "endpoint", value.name = "value", variable.factor = FALSE)
  tests <- rbindlist(lapply(c(1L, 10L), function(threshold) {
    donor_long[, {
      keep <- n_cells >= threshold
      x <- value[keep]
      groups <- factor(subtype[keep], levels = wu_subtypes)
      n <- as.integer(table(groups))
      er <- x[groups == "ER+"]
      tn <- x[groups == "TNBC"]
      kw_ok <- all(n >= 2L) && uniqueN(x) > 1L
      pair_ok <- all(n[c(1L, 3L)] >= 2L) && uniqueN(c(er, tn)) > 1L
      p_kw <- if (kw_ok) kruskal.test(x, groups)$p.value else NA_real_
      p_pair <- if (pair_ok) wilcox.test(er, tn, exact = FALSE)$p.value else NA_real_
      .(min_cells = threshold, n_ER_pos = n[1L], n_HER2_pos = n[2L], n_TNBC = n[3L],
        kw_pval = p_kw, er_vs_tnbc_wilcox_p = p_pair,
        er_minus_tnbc_median = if (length(er) && length(tn)) median(er) - median(tn) else NA_real_)
    }, by = .(celltype, gene, endpoint)]
  }))
  # 独立检验族：每个阈值 × 结局 × 检验类型，覆盖全部可检验的基因和细胞群。
  tests[, `:=`(kw_fdr = p.adjust(kw_pval, "BH"),
               er_vs_tnbc_fdr = p.adjust(er_vs_tnbc_wilcox_p, "BH"),
               n_kw_tests = sum(!is.na(kw_pval)), n_pair_tests = sum(!is.na(er_vs_tnbc_wilcox_p))),
        by = .(min_cells, endpoint)]
  fwrite(tests, file.path(res_dir, "planA_scrna_wu_donor_subtype_tests.csv"))

  group_expr <- rbindlist(lapply(c(1L, 10L), function(threshold) {
    donor_long[n_cells >= threshold, .(min_cells = threshold, n_donors = .N,
      donor_mean = mean(value), donor_median = median(value),
      donor_q1 = unname(quantile(value, 0.25)), donor_q3 = unname(quantile(value, 0.75))),
      by = .(celltype, gene, endpoint, subtype)]
  }))
  fwrite(group_expr, file.path(res_dir, "planA_scrna_wu_donor_subtype_summary.csv"))
  cat("\nWu donor representation (ER+, HER2+, TNBC) by cell type and threshold:\n")
  print(dcast(group_counts, celltype + min_cells ~ subtype, value.var = "n_donors"), nrows = Inf)
  list(tests = tests, panels = panels)
}

write_myeloid_composition <- function(wu, tisch) {
  # TISCH2 标准注释列带括号。无括号列保留 Wu 标签，不能用来代替 TISCH2 重注释。
  stopifnot(!anyDuplicated(wu$Cell), !anyDuplicated(tisch$Cell))
  wu_labels <- wu[, .(Cell, wu_major = celltype_major, wu_minor = celltype_minor)]
  ti_labels <- tisch[, .(Cell, tisch_major = `Celltype (major-lineage)`,
                        tisch_minor = `Celltype (minor-lineage)`)]
  x <- merge(wu_labels, ti_labels, by = "Cell", all = TRUE)
  wu_myeloid <- x$wu_major %in% "Myeloid"
  ti_myeloid <- x$tisch_major %in% c("Mono/Macro", "DC", "Neutrophil", "Mast")
  comp <- rbindlist(list(
    x[wu_myeloid, .(n_cells = .N, n_shared_cells = sum(!is.na(tisch_major))),
      by = .(major = wu_major, minor = wu_minor)][, annotation := "Wu"],
    x[ti_myeloid, .(n_cells = .N, n_shared_cells = sum(!is.na(wu_major))),
      by = .(major = tisch_major, minor = tisch_minor)][, annotation := "TISCH2"]
  ))
  comp[, n_annotation_cells := ifelse(annotation == "Wu", nrow(wu), nrow(tisch))]
  comp[, pct_within_major := 100 * n_cells / sum(n_cells), by = .(annotation, major)]
  fwrite(comp, file.path(res_dir, "planA_scrna_wu_tisch2_myeloid_composition.csv"))
  cross <- x[wu_myeloid | ti_myeloid, .(n_cells = .N),
             by = .(wu_major, wu_minor, tisch_major, tisch_minor)]
  fwrite(cross[order(wu_major, wu_minor, -n_cells)], file.path(res_dir, "planA_scrna_wu_tisch2_myeloid_overlap.csv"))
  cat("\nMyeloid annotation composition. Counts from full metadata per source:\n")
  print(comp)
}

# Wu GSE176078: 使用目标基因稀疏计数、原始元数据和 TISCH2 UMAP。
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

  if (!file.exists(ti_meta_f)) stop("Missing TISCH2 GSE176078 metadata: ", ti_meta_f)
  m_ti <- fread(ti_meta_f, select = c("Cell", "UMAP_1", "UMAP_2",
                                     "Celltype (major-lineage)", "Celltype (minor-lineage)"))
  write_myeloid_composition(m_wu, m_ti)
  m_wu <- merge(m_wu, m_ti[, .(Cell, UMAP_1, UMAP_2)], by = "Cell", all.x = TRUE, sort = FALSE)

  # Harmonize Wu et al. author-curated celltype_major to consistent lineage nomenclature
  lineage_map <- c(
    "Cancer Epithelial" = "Malignant",
    "Normal Epithelial" = "Epithelial",
    "Myeloid"           = "Mono/Macro", # 既有绘图标签；Wu 实际包含 DCs 和 Cycling_Myeloid。
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

  stopifnot(!anyNA(m_wu$Celltype_major), all(is.finite(m_wu$nCount_RNA)), all(m_wu$nCount_RNA > 0),
            !anyDuplicated(cnt[, .(gene, barcode)]), all(cnt$barcode %in% dt_plot$Cell),
            all(target_genes %in% cnt$gene), all(is.finite(cnt$count)), all(cnt$count > 0))
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

  donor_results <- wu_donor_subtype(dt_plot)
  # 保留旧图所需的按细胞加权描述统计及字段名。p 值改为供体均值检验。
  legacy_groups <- c("ALL_CELLS", "Malignant", "Myeloid", "Endothelial", "Fibroblasts",
                     "Epithelial", "PVL", "T_cells")
  sub_res <- melt(donor_results$panels[Celltype_major %in% legacy_groups],
                  measure.vars = subtype_genes, variable.name = "gene", value.name = "expr",
                  variable.factor = FALSE)[,
    .(n_cells = .N, pct_expr = round(100 * mean(expr > 0), 2), mean_expr = round(mean(expr), 4)),
    by = .(celltype = Celltype_major, gene, subtype = Subtype)]
  sub_res <- merge(sub_res, donor_results$tests[min_cells == 1L & endpoint == "mean_expr",
    .(celltype, gene, kw_pval, er_vs_tnbc_wilcox_p, kw_fdr, er_vs_tnbc_fdr)],
    by = c("celltype", "gene"))
  sub_res[celltype == "Myeloid", celltype := "Mono/Macro"]
  sub_res[, `:=`(dataset = ds_name, descriptive_unit = "pooled_cells", test_unit = "donor_mean_expr")]

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
  stopifnot(length(common) > 0L, !anyDuplicated(common), nrow(meta) == length(common))
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
    Sample = Sample,
    Tissue = Tissue,
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

  if (ds_name == "BRCA_GSE148673") write_gse148673_sensitivity(dt_plot, rownames(mat))
  list(plot_dt = dt_plot, celltype_sum = sum_dt, subtype_sum = data.table())
}

write_gse148673_sensitivity <- function(cells, genes) {
  # ECIS05T 对应唯一的 DCIS 样本，比较全队列和去除此样本的描述性统计。
  inventory <- cells[, .(n_cells = .N), by = .(Patient, Sample, Tissue)]
  inventory[, excluded := Patient == "ECIS05T"]
  stopifnot(sum(inventory$excluded) == 1L, grepl("^DCIS", inventory[excluded == TRUE, Sample]),
            all(grepl("^TNBC", inventory[excluded == FALSE, Sample])))
  fwrite(inventory, file.path(res_dir, "planA_scrna_gse148673_sample_inventory.csv"))
  summaries <- function(x) {
    y <- rbindlist(list(x, copy(x)[, Celltype_major := "ALL_CELLS"]))
    melt(y[, c("Patient", "Celltype_major", genes), with = FALSE], measure.vars = genes,
         variable.name = "gene", value.name = "expr", variable.factor = FALSE)[,
      .(n_donors = uniqueN(Patient), n_cells = .N, n_expr = sum(expr > 0),
        pct_expr = 100 * mean(expr > 0), mean_expr = mean(expr)),
      by = .(celltype = Celltype_major, gene)]
  }
  both <- merge(summaries(cells), summaries(cells[Patient != "ECIS05T"]),
                by = c("celltype", "gene"), all.x = TRUE, suffixes = c("_all", "_exclude_ECIS05T"))
  both[is.na(n_cells_exclude_ECIS05T), `:=`(n_cells_exclude_ECIS05T = 0L,
                                             n_donors_exclude_ECIS05T = 0L, n_expr_exclude_ECIS05T = 0L)]
  both[, `:=`(delta_pct_expr = pct_expr_exclude_ECIS05T - pct_expr_all,
              delta_mean_expr = mean_expr_exclude_ECIS05T - mean_expr_all)]
  fwrite(both, file.path(res_dir, "planA_scrna_gse148673_ecis05t_sensitivity.csv"))
  cat("\nGSE148673 retained sample inventory:\n")
  print(inventory)
  cat("\nGSE148673 ECIS05T sensitivity (overall and key lineages):\n")
  print(both[celltype %in% c("ALL_CELLS", "Malignant", "Mono/Macro", "Endothelial") &
               gene %in% c("PPARG", "CXCL8", "IL1B", "FFAR2"),
             .(celltype, gene, n_cells_all, n_cells_exclude_ECIS05T,
               pct_expr_all, pct_expr_exclude_ECIS05T, delta_mean_expr)])
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
