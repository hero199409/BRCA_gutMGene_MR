# 26a_generate_cns_figs6to7.R
# Frontiers 修订：图 6C 使用供体层表达与检验。图 7A 从结果读取分层总数与有效样本数。
# 两张图以 180 mm 成品宽度直接绘制，沿用原有输出基础文件名。

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/25_0_cns_plot_helpers.R"))
suppressPackageStartupMessages({
  library(enrichplot)
})

sc_cache <- readRDS(file.path(sc_dir, "planA_scrna_plot_cache.rds"))
sc_sum   <- fread(file.path(res_dir, "planA_scrna_celltype_summary.csv"))

# 极简角标 UMAP 绘图函数（含左下角微型 L 型坐标箭头 + 簇重心直接标签）
plot_umap_cns_ct <- function(umap_dt, ct_pal, ttl, pt_sz = 0.22, label_clusters = TRUE) {
  d <- copy(umap_dt)
  d[, Lineage_clean := clean_lineage_str(Celltype_major)]
  d[Celltype_major == "Mono/Macro", Lineage_clean := "Myeloid"]
  xr <- range(d$UMAP_1, na.rm = TRUE); yr <- range(d$UMAP_2, na.rm = TRUE)
  dx <- diff(xr); dy <- diff(yr)
  ax_x0 <- xr[1] - 0.04 * dx; ax_y0 <- yr[1] - 0.04 * dy
  ax_x1 <- ax_x0 + 0.16 * dx; ax_y1 <- ax_y0 + 0.16 * dy

  p <- ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = Lineage_clean)) +
    ggrastr::rasterise(geom_point(size = pt_sz, alpha = 0.65), dpi = 320) +
    annotate("segment", x = ax_x0, xend = ax_x1, y = ax_y0, yend = ax_y0,
             color = "#0F172A", linewidth = 0.55, arrow = arrow(length = unit(0.11, "cm"), type = "closed")) +
    annotate("segment", x = ax_x0, xend = ax_x0, y = ax_y0, yend = ax_y1,
             color = "#0F172A", linewidth = 0.55, arrow = arrow(length = unit(0.11, "cm"), type = "closed")) +
    annotate("text", x = (ax_x0 + ax_x1)/2, y = ax_y0 - 0.035 * dy, label = "UMAP 1", size = 2.3, fontface = "bold", color = "#334155") +
    annotate("text", x = ax_x0 - 0.035 * dx, y = (ax_y0 + ax_y1)/2, label = "UMAP 2", size = 2.3, fontface = "bold", color = "#334155", angle = 90) +
    scale_color_manual(values = ct_pal) +
    coord_cartesian(xlim = c(xr[1] - 0.08 * dx, xr[2] + 0.04 * dx),
                    ylim = c(yr[1] - 0.08 * dy, yr[2] + 0.04 * dy), clip = "off") +
    labs(title = ttl) +
    theme_void(base_family = "sans") +
    theme(plot.title = element_text(face = "bold", size = 9.8, color = "#0F172A", margin = margin(b = 2)),
          legend.position = "none", plot.margin = margin(4, 4, 4, 4))

  if (label_clusters) {
    cent <- d[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), N = .N), by = Lineage_clean]
    p <- p + ggrepel::geom_label_repel(
      data = cent, aes(x = UMAP_1, y = UMAP_2, label = Lineage_clean, color = Lineage_clean),
      fill = alpha("white", 0.90), size = 2.65, fontface = "bold",
      label.padding = unit(0.14, "lines"), box.padding = 0.25, max.overlaps = 20,
      segment.size = 0.35, min.segment.length = 0.2, seed = 101
    )
  }
  p
}

plot_umap_cns_gene <- function(umap_dt, gene, ttl, pal_low = "#E2E8F0", pal_mid = NULL, pal_hi = "#DC2626", pt_sz = 0.22) {
  d <- copy(umap_dt[, .(UMAP_1, UMAP_2, val = get(gene))])[order(val)]
  xr <- range(d$UMAP_1, na.rm = TRUE); yr <- range(d$UMAP_2, na.rm = TRUE)
  dx <- diff(xr); dy <- diff(yr)
  ax_x0 <- xr[1] - 0.04 * dx; ax_y0 <- yr[1] - 0.04 * dy
  ax_x1 <- ax_x0 + 0.16 * dx; ax_y1 <- ax_y0 + 0.16 * dy

  col_sc <- if (!is.null(pal_mid)) {
    scale_color_gradientn(colors = c(pal_low, pal_mid, pal_hi), name = NULL)
  } else {
    scale_color_gradient(low = pal_low, high = pal_hi, name = NULL)
  }

  ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = val)) +
    ggrastr::rasterise(geom_point(size = pt_sz, alpha = ifelse(d$val > 0, 0.88, 0.25)), dpi = 320) +
    annotate("segment", x = ax_x0, xend = ax_x1, y = ax_y0, yend = ax_y0,
             color = "#0F172A", linewidth = 0.50, arrow = arrow(length = unit(0.10, "cm"), type = "closed")) +
    annotate("segment", x = ax_x0, xend = ax_x0, y = ax_y0, yend = ax_y1,
             color = "#0F172A", linewidth = 0.50, arrow = arrow(length = unit(0.10, "cm"), type = "closed")) +
    col_sc +
    coord_cartesian(xlim = c(xr[1] - 0.06 * dx, xr[2] + 0.04 * dx),
                    ylim = c(yr[1] - 0.06 * dy, yr[2] + 0.04 * dy), clip = "off") +
    labs(title = ttl) +
    theme_void(base_family = "sans") +
    theme(plot.title = element_text(face = "bold", size = 9.8, color = "#0F172A", margin = margin(b = 2)),
          legend.position = "inside", legend.position.inside = c(0.95, 0.22),
          legend.key.height = unit(0.45, "cm"), legend.key.width = unit(0.22, "cm"),
          legend.text = element_text(size = 7.2, color = "#334155"),
          plot.margin = margin(4, 4, 4, 4))
}

# ==============================================================================
# FIG 6: 肿瘤单细胞定位、供体层亚型表达和 HPA 病理图
# ==============================================================================
cat("=== Building Fig 6 (UMAP, donor expression, HPA) ===\n")

sc_df    <- copy(sc_cache[["BRCA_GSE176078"]])
sc_umap  <- sc_df[!is.na(UMAP_1) & !is.na(UMAP_2)]
sc_sum_p <- copy(sc_sum[dataset == "BRCA_GSE176078"])

ct_raw6   <- sc_umap[, .N, by = Celltype_major][order(-N)]$Celltype_major
ct_clean6 <- ifelse(ct_raw6 == "Mono/Macro", "Myeloid", clean_lineage_str(ct_raw6))
ct_pal6   <- setNames(colorRampPalette(c(pal_cns[["blue"]], pal_cns[["red"]], pal_cns[["green"]], pal_cns[["cyan"]],
                                         pal_cns[["coral"]], pal_cns[["slate"]], pal_cns[["amber"]], "#7C3AED", "#EC4899"))(length(ct_clean6)), ct_clean6)

p6a_ct <- plot_umap_cns_ct(sc_umap, ct_pal6, "A1  Cell lineages (GSE176078)")
p6a_cx <- plot_umap_cns_gene(sc_umap, "CXCL8", "A2  CXCL8 expression", pal_mid = "#F59E0B", pal_hi = "#9F1239")
p6a_pp <- plot_umap_cns_gene(sc_umap, "PPARG", "A3  PPARG expression", pal_mid = "#60A5FA", pal_hi = "#1E3A8A")
p6a_pl <- plot_umap_cns_gene(sc_umap, "PLIN1", "A4  Adipocyte marker PLIN1", pal_hi = "#334155")

dot_genes6 <- c("CXCL8","IL1B","FFAR2","PPARG","FABP4","ADIPOQ","PLIN1","LEP","EPCAM","VWF","CD68")
dot_dt6 <- copy(sc_sum_p[gene %in% dot_genes6])
dot_dt6[, gene := factor(gene, levels = dot_genes6)]
ct_ord6 <- clean_lineage_str(c("Malignant","Epithelial","Endothelial","PVL","Fibroblasts","Mono/Macro","T_cells","B","Plasma"))
ct_ord6[ct_ord6 == "Mono/Macrophages"] <- "Myeloid"
dot_dt6[, Lineage := factor(ifelse(celltype_major == "Mono/Macro", "Myeloid",
                                  clean_lineage_str(celltype_major)), levels = rev(ct_ord6))]

p6b <- ggplot(dot_dt6, aes(x = gene, y = Lineage, size = pct_expr, color = mean_expr)) +
  geom_point() +
  scale_size_area(max_size = 4.3, name = "% expressing") +
  scale_color_gradientn(colors = c("#F1F5F9", "#38BDF8", "#F97316", "#BE123C"), name = "Mean expression") +
  labs(title = "B  Lineage expression (GSE176078)", x = NULL, y = NULL) +
  theme_cns_grid(9.2) + theme(axis.text.x = element_text(face = "italic"),
                             legend.position = "right")

# 17b 的主分析阈值为每供体每谱系 >= 1 个细胞。图中只展示 mean_expr 端点。
lineages6c <- c("ALL_CELLS", "Endothelial", "Malignant", "Myeloid")
genes6c <- c("PPARG", "CXCL8")
subtypes6c <- c("ER+", "HER2+", "TNBC")
min_cells6c <- 1L
keys6c <- c("celltype", "gene", "subtype")
files6c <- file.path(res_dir, paste0("planA_scrna_wu_donor_", c("expression", "subtype_summary", "subtype_tests"), ".csv"))
preview6c <- !all(file.exists(files6c))
if (!preview6c) {
  donor6c <- fread(files6c[1])[celltype %chin% lineages6c & gene %chin% genes6c & n_cells >= min_cells6c]
  group6c <- fread(files6c[2])[celltype %chin% lineages6c & gene %chin% genes6c &
                               endpoint == "mean_expr" & min_cells == min_cells6c]
  tests6c <- fread(files6c[3])[celltype %chin% lineages6c & gene %chin% genes6c &
                               endpoint == "mean_expr" & min_cells == min_cells6c]
  counts6c <- donor6c[, .(plotted = .N), by = keys6c]
  check6c <- merge(counts6c, group6c, by = keys6c, all = TRUE)
  tested_counts6c <- melt(tests6c, id.vars = c("celltype", "gene"),
                          measure.vars = c("n_ER_pos", "n_HER2_pos", "n_TNBC"),
                          variable.name = "group_id", value.name = "tested")
  tested_counts6c[, subtype := subtypes6c[as.integer(group_id)]]
  check6c <- merge(check6c, tested_counts6c, by = keys6c, all = TRUE)
  stopifnot(nrow(tests6c) == length(lineages6c) * length(genes6c),
            !anyDuplicated(tests6c, by = c("celltype", "gene")),
            !anyDuplicated(check6c, by = keys6c),
            all(fcoalesce(check6c$plotted, 0L) == fcoalesce(check6c$n_donors, 0L)),
            all(fcoalesce(check6c$plotted, 0L) == check6c$tested))
  cat("Fig 6C: using 17b donor results (mean_expr, min_cells =", min_cells6c, ")\n")
} else {
  # 仅作预览：缓存按供体聚合，不绘制检验值。
  stopifnot(all(c("Patient", "Subtype", "Celltype_major", genes6c) %in% names(sc_df)))
  cells6c <- melt(sc_df[!is.na(Patient) & nzchar(Patient) & Subtype %chin% subtypes6c],
                   id.vars = c("Patient", "Subtype", "Celltype_major"), measure.vars = genes6c,
                   variable.name = "gene", value.name = "expr", variable.factor = FALSE)
  cells6c[Celltype_major == "Mono/Macro", Celltype_major := "Myeloid"]
  cells6c <- rbindlist(list(cells6c[Celltype_major %chin% lineages6c],
                            cells6c[, .(Patient, Subtype, Celltype_major = "ALL_CELLS", gene, expr)]))
  donor6c <- cells6c[is.finite(expr), .(mean_expr = mean(expr), n_cells = .N),
                    by = .(celltype = Celltype_major, gene, subtype = Subtype, donor = Patient)][n_cells >= min_cells6c]
  group6c <- donor6c[, .(n_donors = .N, donor_median = median(mean_expr),
                         donor_q1 = quantile(mean_expr, 0.25), donor_q3 = quantile(mean_expr, 0.75)), by = keys6c]
  tests6c <- CJ(celltype = lineages6c, gene = genes6c)[, `:=`(kw_pval = NA_real_, kw_fdr = NA_real_)]
  message("Fig 6C: PROVISIONAL, donor summaries from the plot cache. Missing: ",
          paste(basename(files6c[!file.exists(files6c)]), collapse = ", "))
}
stopifnot(nrow(donor6c) > 0L, all(donor6c$subtype %chin% subtypes6c),
          all(is.finite(donor6c$mean_expr)), !anyDuplicated(donor6c, by = c(keys6c, "donor")),
          !anyDuplicated(unique(donor6c[, .(donor, subtype)])$donor))
format_6c_facet <- function(dt) {
  dt[, `:=`(celltype_lbl = factor(celltype, levels = lineages6c,
                                 labels = c("All cells", "Endothelial", "Malignant", "Myeloid")),
            gene = factor(gene, levels = genes6c))]
}
counts6c <- merge(CJ(celltype = lineages6c, gene = genes6c, subtype = subtypes6c),
                  group6c[, c(keys6c, "n_donors"), with = FALSE], by = keys6c, all.x = TRUE)
counts6c[is.na(n_donors), n_donors := 0L]
row_span6c <- donor6c[, .(span = max(0.01, max(mean_expr))), by = gene]
counts6c <- merge(counts6c, row_span6c, by = "gene")
tests6c <- merge(tests6c, row_span6c, by = "gene")
tests6c[, annotation := ifelse(is.na(kw_pval), if (preview6c) "" else "KW not tested",
                               sprintf("P = %.2g  |  q = %.2g", kw_pval, kw_fdr))]
for (dt6c in list(donor6c, group6c, tests6c, counts6c)) format_6c_facet(dt6c)
for (dt6c in list(donor6c, group6c, counts6c)) dt6c[, subtype := factor(subtype, levels = subtypes6c)]

p6c <- ggplot(donor6c, aes(subtype, mean_expr)) +
  geom_hline(yintercept = 0, color = "#CBD5E1", linewidth = 0.3) +
  geom_point(aes(color = subtype), position = position_jitter(width = 0.12, height = 0, seed = 26), size = 1.45, alpha = 0.72) +
  geom_errorbar(data = group6c, aes(x = subtype, ymin = donor_q1, ymax = donor_q3),
                inherit.aes = FALSE, width = 0.13, linewidth = 0.5, color = "#0F172A") +
  geom_point(data = group6c, aes(x = subtype, y = donor_median), inherit.aes = FALSE,
             shape = 23, size = 1.9, stroke = 0.55, fill = "white", color = "#0F172A") +
  geom_text(data = counts6c, aes(x = subtype, y = -0.12 * span, label = paste0("n=", n_donors)),
            inherit.aes = FALSE, size = 2.65, color = "#334155") +
  geom_text(data = tests6c, aes(x = 2, y = 1.12 * span, label = annotation),
            inherit.aes = FALSE, size = 2.65, color = "#334155") +
  facet_grid(gene ~ celltype_lbl, scales = "free_y") +
  scale_x_discrete(limits = subtypes6c, drop = FALSE) +
  scale_y_continuous(breaks = function(limits) { b <- scales::breaks_pretty(n = 3)(limits); b[b >= 0] },
                     expand = expansion(mult = c(0.04, 0.05))) +
  scale_color_manual(values = c("ER+" = pal_cns[["blue"]], "HER2+" = pal_cns[["amber"]], "TNBC" = pal_cns[["red"]])) +
  labs(title = "C  Donor expression by tumor subtype", x = NULL, y = "Donor mean\nnormalized expression",
       subtitle = if (preview6c) "PROVISIONAL: donor means from plot cache" else
         "Dots: donors. Diamonds: median. Bars: IQR. P: Kruskal–Wallis. q: BH-adjusted P.") +
  theme_cns_grid(9.2) + theme(legend.position = "none", panel.grid.major.x = element_blank(),
                             plot.subtitle = element_text(size = 8.3, margin = margin(b = 4)))

# 将局部放大标注缩短到 180 mm 成品布局可容纳的宽度。
make_hpa_panel26a <- function(...) {
  panel <- make_hpa_inset_panel(...)
  inset_label <- which(vapply(panel$layers, function(layer) inherits(layer$geom, "GeomLabel"), logical(1)))
  stopifnot(length(inset_label) == 1L)
  panel$layers[[inset_label]]$aes_params$label <- "Inset"
  panel
}

# 6D1-D3: 原图 3000x3000 HPA 切片 + 局部高倍放大镜插图 (Zoom-in Inset)
p6d1 <- make_hpa_panel26a(
  file.path(ihc_dir, "HPA_PPARG_NormalBreast_IHC.jpg"),
  crop_box = c(0.54, 0.80, 0.68, 0.94),
  title_lbl = "D1  PPARG | Normal breast",
  caption_txt = "Low in adipocytes\nNot detected in ducts",
  max_px = 820,
  accent_col = pal_cns[["blue"]]
)
p6d2 <- make_hpa_panel26a(
  file.path(ihc_dir, "HPA_PPARG_BreastCancer_IHC.jpg"),
  crop_box = c(0.38, 0.64, 0.48, 0.74),
  title_lbl = "D2  PPARG | Carcinoma",
  caption_txt = "4 medium, 6 low\n2 not detected",
  max_px = 820,
  accent_col = pal_cns[["blue"]]
)
p6d3 <- make_hpa_panel26a(
  file.path(ihc_dir, "HPA_CXCL8_BreastCancer_IHC.jpg"),
  crop_box = c(0.45, 0.71, 0.61, 0.87),
  title_lbl = "D3  CXCL8 | Carcinoma",
  caption_txt = "Not detected (12/12)",
  max_px = 820,
  accent_col = pal_cns[["red"]]
)

fig6_a <- patchwork::wrap_elements(full = ((p6a_ct | p6a_cx) / (p6a_pp | p6a_pl)))
fig6_d <- patchwork::wrap_elements(full = (p6d1 | p6d2 | p6d3))
fig6_all <- wrap_plots(fig6_a, wrap_elements(full = p6b), wrap_elements(full = p6c), fig6_d,
                       ncol = 1, heights = c(87, 42, 63, 40))
save_cns(fig6_all, "Fig6_Tumor_scRNA_and_HPA_Pathology", w = 180 / 25.4, h = 232 / 25.4)
rm(fig6_all, fig6_a, fig6_d, p6a_ct, p6a_cx, p6a_pp, p6a_pl, p6b, p6c, p6d1, p6d2, p6d3, sc_df, sc_umap)
invisible(gc())

# ==============================================================================
# FIG 7: 斑马纹 TME 相关性 (7A) + 双轴分色 GSEA & Running Curve 峰图 (7B) + 跨队列单细胞 (7C)
# ==============================================================================
cat("=== Building Fig 7 (TME + Dual-Axis GSEA Running Curves + Cross-Cohort scRNA) ===\n")

tme_map <- c(
  "quanTIseq__B cell"                  = "B Cells (quanTIseq)",
  "quanTIseq__Macrophage M1"           = "M1 Macrophages (quanTIseq)",
  "MCPcounter__Neutrophil"             = "Neutrophils (MCPcounter)",
  "quanTIseq__NK cell"                 = "NK Cells (quanTIseq)",
  "TMB_NONSYNONYMOUS"                  = "Nonsynonymous TMB",
  "MCPcounter__T cell CD8+"            = "CD8+ T Cells (MCPcounter)",
  "quanTIseq__Neutrophil"              = "Neutrophils (quanTIseq)",
  "quanTIseq__T cell CD8+"             = "CD8+ T Cells (quanTIseq)",
  "MCPcounter__T cell"                 = "T Cells (MCPcounter)",
  "MCPcounter__B cell"                 = "B Cells (MCPcounter)",
  "MCPcounter__Endothelial cell"       = "Endothelial Cells (MCPcounter)",
  "ESTIMATE__stroma score"             = "Stromal Score (ESTIMATE)",
  "ESTIMATE__tumor purity"             = "Tumor Purity (ESTIMATE)",
  "ESTIMATE__estimate score"           = "ESTIMATE Score",
  "MCPcounter__cytotoxicity score"     = "Cytotoxicity Score (MCPcounter)",
  "quanTIseq__uncharacterized cell"    = "Uncharacterized Cells (quanTIseq)",
  "MCPcounter__Myeloid dendritic cell" = "Myeloid Dendritic Cells (MCPcounter)"
)

cor5 <- fread(file.path(res_dir, "module5_tme_checkpoint_correlations.csv"))
strata7a <- fread(file.path(res_dir, "module5_tcga_er_sample_counts.csv"))[er_status %chin% c("Negative", "Positive")]
stopifnot(nrow(strata7a) == 2L, !anyDuplicated(strata7a$er_status), all(strata7a$n_samples > 0L))
strata7a[, stratum := ifelse(er_status == "Negative", "ER_neg", "ER_pos")]
cor5 <- merge(cor5, strata7a[, .(stratum, n_stratum = n_samples)], by = "stratum")
tme_top <- cor5[((gene_symbol == "CXCL8" & stratum == "ER_neg") | (gene_symbol == "PPARG" & stratum == "ER_pos")) &
                  feature_type %chin% c("TME_Cell_Score", "Checkpoint", "Immune_Checkpoint", "Genomic_Feature") & FDR < 0.05][
                    order(gene_symbol, -abs(rho))][, head(.SD, 10L), by = gene_symbol]
stopifnot(uniqueN(tme_top$gene_symbol) == 2L, all(tme_top$n > 0L & tme_top$n <= tme_top$n_stratum))
tme_top[, feat_clean := ifelse(feature %chin% names(tme_map), unname(tme_map[feature]), gsub("_", " ", feature))]
tme_top[, facet_lbl := sprintf("%s (N = %d)",
                              ifelse(gene_symbol == "CXCL8", "A1  CXCL8 | ER−", "A2  PPARG | ER+"), n_stratum)]
tme_top[, dir := ifelse(rho > 0, "Positive", "Negative")]
tme_top[, feat_ord := paste0(gene_symbol, "___", feat_clean)]
tme_top <- tme_top[order(facet_lbl, rho)]
tme_top[, feat_ord := factor(feat_ord, levels = unique(feat_ord))]
tme_top[, y_local := seq_len(.N), by = facet_lbl]
tme_zebra <- tme_top[y_local %% 2 == 1, .(facet_lbl, y_local)]
cat("Fig 7A strata and plotted complete-pair counts:\n")
print(tme_top[, .(stratum_N = unique(n_stratum), min_n = min(n), max_n = max(n)), by = facet_lbl])
label_7a <- function(x) {
  x <- sub("^[^_]+___", "", x)
  sub(" \\((quanTIseq|MCPcounter|ESTIMATE)\\)$", "\n(\\1)", x)
}

p7a <- ggplot(tme_top, aes(x = rho, y = feat_ord, color = dir)) +
  geom_rect(data = tme_zebra, aes(ymin = y_local - 0.47, ymax = y_local + 0.47, xmin = -Inf, xmax = Inf),
            inherit.aes = FALSE, fill = "#F8FAFC", color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = rho, yend = feat_ord), linewidth = 0.72) +
  geom_point(size = 2.1) +
  geom_text(aes(x = 0.89, label = paste0("n=", n)), hjust = 1, size = 2.65, color = "#1E293B") +
  facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = label_7a) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  scale_color_manual(values = c("Positive" = pal_cns[["red"]], "Negative" = pal_cns[["blue"]])) +
  coord_cartesian(xlim = c(min(-0.55, min(tme_top$rho) - 0.05), 0.96)) +
  labs(title = "A  Tumor microenvironment correlations (TCGA-BRCA)", x = "Spearman ρ", y = NULL,
       subtitle = "N denotes ER stratum total. n denotes complete pairs for each feature.") +
  theme_cns(9.2) + theme(legend.position = "none", axis.text.y = element_text(lineheight = 0.88),
                         plot.subtitle = element_text(size = 8.3, margin = margin(b = 4)))

# 7B: 双轴严密分色 GSEA Lollipop + 真实 GSEA Running Enrichment Score 峰图
gsea_map <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"           = "TNF\u03b1\u2013NF-\u03baB Signaling",
  "HALLMARK_MTORC1_SIGNALING"                  = "mTORC1 Signaling",
  "HALLMARK_G2M_CHECKPOINT"                    = "G2/M Checkpoint",
  "HALLMARK_HYPOXIA"                           = "Hypoxia",
  "HALLMARK_GLYCOLYSIS"                        = "Glycolysis",
  "HALLMARK_INFLAMMATORY_RESPONSE"             = "Inflammatory Response",
  "HALLMARK_E2F_TARGETS"                       = "E2F Targets",
  "HALLMARK_ANGIOGENESIS"                      = "Angiogenesis",
  "HALLMARK_ALLOGRAFT_REJECTION"               = "Allograft Rejection",
  "HALLMARK_KRAS_SIGNALING_UP"                 = "KRAS Signaling Up",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"           = "IL-6/JAK/STAT3 Signaling",
  "HALLMARK_COAGULATION"                       = "Coagulation",
  "HALLMARK_MYOGENESIS"                        = "Myogenesis"
)

gsea5 <- fread(file.path(res_dir, "module5_hallmark_gsea_results.csv"))
gs_top <- gsea5[((target_gene=="CXCL8" & stratum=="ER_neg") | (target_gene=="PPARG" & stratum=="ER_pos")) & p.adjust < 0.05][
  order(target_gene, -abs(NES))][, .SD[1:min(8, .N)], by = target_gene]
gs_top[, clean_id := ifelse(ID %in% names(gsea_map), unname(gsea_map[ID]), gsub("_", " ", sub("^HALLMARK_", "", ID)))]
gs_top[, facet_lbl := ifelse(target_gene=="CXCL8", "B1  CXCL8 | ER\u2212", "B2  PPARG | ER+")]
gs_top[, axis_col := ifelse(target_gene=="CXCL8", "ER\u2212 CXCL8 Axis", "ER+ PPARG Axis")]
gs_top[, id_ord := paste0(target_gene, "___", clean_id)]
gs_top <- gs_top[order(facet_lbl, NES)]; gs_top[, id_ord := factor(id_ord, levels = id_ord)]
gs_top[, y_local := seq_len(.N), by = facet_lbl]
gs_zebra <- gs_top[y_local %% 2 == 1, .(facet_lbl, y_local)]

p7b_bar <- ggplot(gs_top, aes(x = NES, y = id_ord, color = axis_col)) +
  geom_rect(data = gs_zebra, aes(ymin = y_local - 0.45, ymax = y_local + 0.45, xmin = -Inf, xmax = Inf),
            inherit.aes = FALSE, fill = "#F8FAFC", color = NA) +
  geom_vline(xintercept = 0, color = "#94A3B8", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = NES, y = id_ord, yend = id_ord), linewidth = 1.05, show.legend = FALSE) +
  geom_point(aes(size = -log10(p.adjust))) +
  geom_text(aes(x = NES + 0.19, label = sprintf("%.2f", NES)), hjust = 0, size = 2.65, fontface = "bold", color = "#1E293B", show.legend = FALSE) +
  facet_wrap(~ facet_lbl, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^_]+___", "", x)) +
  scale_color_manual(values = c("ER\u2212 CXCL8 Axis" = pal_cns[["red"]], "ER+ PPARG Axis" = pal_cns[["blue"]]), guide = "none") +
  scale_size(range = c(2.4, 5.0), name = expression(bold(-log[10](FDR)))) +
  coord_cartesian(xlim = c(min(0, min(gs_top$NES) - 0.15), max(gs_top$NES) + 0.75)) +
  labs(title = "B  Hallmark pathway enrichment", x = "Normalized enrichment score (NES)", y = NULL) +
  theme_cns(9.2) + theme(legend.position = "bottom", axis.text.y = element_text(size = 8.2))

# 真实 GSEA Running Score 山峰曲线 (CXCL8 Top 3 + PPARG Top 2)
gsea_objs <- readRDS(file.path(bulk_dir, "planA_cns_gsea_objects.rds"))
curve_specs <- list(
  list(obj = gsea_objs$CXCL8_ERneg, id = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", lbl = "CXCL8: EMT",         col = "#DC2626"),
  list(obj = gsea_objs$CXCL8_ERneg, id = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",           lbl = "CXCL8: TNF\u03b1–NF-\u03baB", col = "#F97316"),
  list(obj = gsea_objs$PPARG_ERpos, id = "HALLMARK_ALLOGRAFT_REJECTION",               lbl = "PPARG: Allograft rejection", col = "#1E3A8A"),
  list(obj = gsea_objs$PPARG_ERpos, id = "HALLMARK_ANGIOGENESIS",                      lbl = "PPARG: Angiogenesis",  col = "#00A087")
)
curve_dt_list <- vector("list", length(curve_specs))
tick_dt_list  <- vector("list", length(curve_specs))
for (k in seq_along(curve_specs)) {
  sp <- curve_specs[[k]]
  gi <- as.data.table(enrichplot:::gsInfo(sp$obj, geneSetID = sp$id))
  # 降采样折线到 ~400 点保持矢量轻量化
  idx_keep <- unique(c(seq(1, nrow(gi), by = max(1L, floor(nrow(gi)/420))), which.max(abs(gi$runningScore)), nrow(gi)))
  curve_dt_list[[k]] <- gi[idx_keep, .(x, runningScore, pathway = sp$lbl)]
  tick_dt_list[[k]]  <- gi[position == 1, .(x, pathway = sp$lbl, y0 = -0.05 - (k - 1) * 0.045, y1 = -0.085 - (k - 1) * 0.045)]
}
curve_dt <- rbindlist(curve_dt_list)
tick_dt  <- rbindlist(tick_dt_list)
pw_pal <- setNames(sapply(curve_specs, `[[`, "col"), sapply(curve_specs, `[[`, "lbl"))

p7b_run <- ggplot(curve_dt, aes(x = x, y = runningScore, color = pathway)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.45) +
  geom_line(linewidth = 0.95) +
  geom_segment(data = tick_dt, aes(x = x, xend = x, y = y0, yend = y1, color = pathway), linewidth = 0.25, alpha = 0.75) +
  scale_color_manual(values = pw_pal, name = NULL) +
  labs(title = "B3  Representative running enrichment curves", x = "Gene rank", y = "Running ES") +
  theme_cns(9.2) + theme(legend.position = "right", legend.text = element_text(size = 8.2))

# 7C: 跨队列单细胞一致性气泡图
shared_ct <- c("Endothelial", "Mono/Macro", "Fibroblasts", "Malignant")
rep_genes <- c("PPARG", "FABP4", "CXCL8", "IL1B", "PLIN1", "ADIPOQ")
grid_all  <- CJ(dataset = c("BRCA_GSE176078", "BRCA_EMTAB8107", "BRCA_GSE148673"),
                celltype_major = shared_ct, gene = rep_genes)
rep_dt <- merge(grid_all, sc_sum[, .(dataset, celltype_major, gene, pct_expr, mean_expr)],
                by = c("dataset", "celltype_major", "gene"), all.x = TRUE)
cohorts7c <- c(BRCA_GSE176078 = "GSE176078", BRCA_EMTAB8107 = "E-MTAB-8107", BRCA_GSE148673 = "GSE148673")
counts7c <- unique(sc_sum[dataset %chin% names(cohorts7c), .(dataset, celltype_major, n_cells)])
stopifnot(!anyDuplicated(counts7c, by = c("dataset", "celltype_major")),
          setequal(counts7c$dataset, names(cohorts7c)))
counts7c <- counts7c[, .(cells = sum(n_cells)), by = dataset]
counts7c[, ds_lbl := paste0(cohorts7c[dataset], " (", format(cells, big.mark = ",", trim = TRUE), ")")]
label_order7c <- counts7c$ds_lbl[match(names(cohorts7c), counts7c$dataset)]
rep_dt <- merge(rep_dt, counts7c[, .(dataset, ds_lbl)], by = "dataset")
rep_dt[, Lineage := factor(celltype_major, levels = shared_ct,
                           labels = c("Endo.", "Myel.*", "Fibro.", "Malig."))]
rep_dt[, gene := factor(gene, levels = rep_genes)]
rep_dt[, is_na := is.na(pct_expr)]

p7c <- ggplot(rep_dt, aes(x = Lineage, y = ds_lbl)) +
  geom_tile(data = rep_dt[is_na == TRUE], fill = "#F1F5F9", color = "#E2E8F0") +
  geom_point(data = rep_dt[is_na == FALSE], aes(size = pct_expr, color = mean_expr)) +
  geom_text(data = rep_dt[is_na == TRUE], label = "NA", size = 2.65, color = "#64748B") +
  facet_wrap(~ gene, ncol = 3) +
  scale_size_area(max_size = 4.9, name = "% expressing") +
  scale_color_gradientn(colors = c("#E0F2FE", "#38BDF8", "#F97316", "#BE123C"), name = "Mean expression") +
  scale_y_discrete(limits = rev(label_order7c)) +
  labs(title = "C  Cross-cohort lineage expression", x = NULL, y = NULL,
       subtitle = sprintf("%s cells. *Myeloid (Wu GSE176078) or Mono/Macro (other cohorts).\nNA indicates unavailable data.",
                          format(sum(counts7c$cells), big.mark = ","))) +
  theme_cns_grid(9.2) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1), legend.position = "right",
        strip.text = element_text(face = "bold.italic", size = 8.7),
        plot.subtitle = element_text(size = 8.3, margin = margin(b = 4)),
        panel.spacing.x = unit(0.6, "lines"))

fig7_rows <- lapply(list(p7a, p7b_bar, p7b_run, p7c), function(panel) wrap_elements(full = panel))
fig7_all <- wrap_plots(fig7_rows, ncol = 1, heights = c(76, 57, 36, 60))
save_cns(fig7_all, "Fig7_TME_GSEA_and_CrossCohort_scRNA", w = 180 / 25.4, h = 229 / 25.4)
cat("Part 2A (Fig 6 - Fig 7) Complete!\n")
