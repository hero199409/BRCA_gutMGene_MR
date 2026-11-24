# 02_fig2_scrna_hpa.R
# Figure 2 CNS-Grade Polish (Nature Medicine / Nature Immunology / Cancer Cell aesthetic)
# Strict data & statistical parity with 05_revision_v2/01_scripts/02_fig2_scrna_hpa.R

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "06_figures_cns_polish/01_scripts/00_cns_theme.R"))

cat("=== Generating CNS-Polished Fig 2 ===\n")

sc_cache <- readRDS(file.path(sc_dir, "planA_scrna_plot_cache.rds"))
sc_sum   <- fread(file.path(res_dir, "planA_scrna_celltype_summary.csv"))
g14_sens <- fread(file.path(res_dir, "planA_scrna_gse148673_ecis05t_sensitivity.csv"))

# ------------------------------------------------------------------------------
# Panel A: Cross-atlas single-cell localization
# ------------------------------------------------------------------------------
shared_ct    <- c("Mono/Macro", "Endothelial", "Fibroblasts", "Malignant")
lineage_lvls <- c("Myeloid", "Endothelial", "Fibroblasts", "Malignant")
genes_2a     <- c("CXCL8", "IL1B", "PPARG", "FABP4", "PLIN1")
cohorts_2a   <- c("Wu (GSE176078)", "E-MTAB-8107", "GSE148673*")

wu_em <- sc_sum[dataset %chin% c("BRCA_GSE176078", "BRCA_EMTAB8107") &
                  celltype_major %chin% shared_ct & gene %chin% genes_2a,
                .(atlas = fifelse(dataset == "BRCA_GSE176078", "Wu (GSE176078)", "E-MTAB-8107"),
                  celltype = celltype_major, gene, pct_expr, mean_expr)]

g14_clean <- g14_sens[celltype %chin% shared_ct & gene %chin% genes_2a,
                      .(atlas = "GSE148673*", celltype, gene,
                        pct_expr = pct_expr_exclude_ECIS05T,
                        mean_expr = mean_expr_exclude_ECIS05T)]

grid_2a <- CJ(atlas = cohorts_2a, celltype = shared_ct, gene = genes_2a)
dt_2a <- merge(grid_2a, rbind(wu_em, g14_clean), by = c("atlas", "celltype", "gene"), all.x = TRUE)
dt_2a[, atlas := factor(atlas, levels = cohorts_2a)]
dt_2a[, gene_md := factor(gene, levels = genes_2a, labels = sprintf("<i>%s</i>", genes_2a))]
dt_2a[, Lineage := factor(celltype, levels = rev(shared_ct), labels = rev(lineage_lvls))]
dt_2a[, missing_gene := is.na(pct_expr)]
dt_2a[, pct_lbl := fifelse(!missing_gene & pct_expr >= 12, sprintf("%.0f%%", pct_expr), "")]
dt_2a[, is_big := !missing_gene & pct_expr >= 45]

hl_2a <- data.table(
  gene_md = factor(c("<i>CXCL8</i>", "<i>PPARG</i>"), levels = sprintf("<i>%s</i>", genes_2a)),
  ymin = c(3.55, 2.55), ymax = c(4.45, 3.45),
  fill_col = c(pal_cns$cxcl8_soft, pal_cns$pparg_soft)
)

# Subtle top hairline under each borderless facet strip
strip_line_2a <- data.table(
  gene_md = factor(sprintf("<i>%s</i>", genes_2a), levels = sprintf("<i>%s</i>", genes_2a))
)

p2a <- ggplot(dt_2a, aes(x = atlas, y = Lineage)) +
  geom_rect(data = hl_2a, aes(ymin = ymin, ymax = ymax, fill = fill_col),
            xmin = -Inf, xmax = Inf, inherit.aes = FALSE, alpha = 0.75) +
  scale_fill_identity() +
  geom_hline(yintercept = Inf, color = pal_cns$ref_line, linewidth = 0.40) +
  geom_tile(data = dt_2a[missing_gene == TRUE], width = 0.80, height = 0.74,
            fill = pal_cns$zebra, color = "#E2E8F0", linewidth = 0.30) +
  geom_text(data = dt_2a[missing_gene == TRUE], label = "NP", size = 2.15,
            color = pal_cns$muted) +
  geom_point(data = dt_2a[missing_gene == FALSE],
             aes(size = pct_expr, color = mean_expr)) +
  geom_text(data = dt_2a[nzchar(pct_lbl) & !is_big],
            aes(label = pct_lbl), nudge_y = -0.37, size = 2.0,
            color = pal_cns$ink, show.legend = FALSE) +
  geom_text(data = dt_2a[nzchar(pct_lbl) & is_big],
            aes(label = pct_lbl), size = 2.0, fontface = "bold",
            color = "white", show.legend = FALSE) +
  facet_wrap(~ gene_md, nrow = 1) +
  scale_x_discrete(limits = cohorts_2a) +
  scale_y_discrete(limits = rev(lineage_lvls), expand = expansion(add = c(0.65, 0.55))) +
  scale_size_area(max_size = 6.2, breaks = c(10, 25, 45), name = "% cells") +
  scale_color_gradientn(
    colors = c("#E2E8F0", "#F39B7F", "#E64B35", "#991B1B"),
    name = "Mean expr."
  ) +
  labs(title = cns_panel_title("A", "Cross-atlas single-cell localization"), x = NULL, y = NULL) +
  theme_cns_grid(8.2) +
  theme(
    axis.line.x       = element_blank(),
    axis.line.y       = element_line(color = pal_cns$ink, linewidth = 0.36),
    axis.ticks.x      = element_line(color = pal_cns$ref_line, linewidth = 0.30),
    axis.text.x       = element_text(angle = 22, hjust = 1, vjust = 1, size = 7.1, color = pal_cns$subink),
    axis.text.y       = element_text(size = 7.6, color = pal_cns$ink),
    strip.text        = ggtext::element_markdown(size = 8.4, color = pal_cns$ink, margin = margin(t = 1, b = 3)),
    panel.spacing.x   = unit(3.2, "mm"),
    legend.position   = "right",
    legend.box        = "vertical",
    legend.key.height = unit(3.2, "mm"),
    legend.key.width  = unit(2.6, "mm"),
    legend.margin     = margin(0, 0, 0, 2)
  )

# ------------------------------------------------------------------------------
# Panel B: Cell-level Wu UMAPs (Lineages, CXCL8, PPARG)
# ------------------------------------------------------------------------------
sc_df <- copy(sc_cache[["BRCA_GSE176078"]])
sc_umap <- sc_df[!is.na(UMAP_1) & !is.na(UMAP_2)]
sc_umap[, Lineage := clean_lineage_rev(Celltype_major)]

ct_levels <- sc_umap[, .N, by = Lineage][order(-N)]$Lineage
ct_pal <- c(
  "T Cells"          = "#3C5488",
  "Malignant"        = "#E64B35",
  "Myeloid"          = "#00A087",
  "Endothelial"      = "#4DBBD5",
  "Fibroblasts"      = "#F39B7F",
  "PVL"              = "#8491B4",
  "Norm. Epithelial" = "#E69F00",
  "Plasma"           = "#B09C85",
  "B Cells"          = "#7E6148"
)
ct_pal <- ct_pal[intersect(ct_levels, names(ct_pal))]

xr <- range(sc_umap$UMAP_1); yr <- range(sc_umap$UMAP_2)
dx <- diff(xr); dy <- diff(yr)
ax_x0 <- xr[1] - 0.03 * dx; ax_y0 <- yr[1] - 0.03 * dy
ax_x1 <- ax_x0 + 0.16 * dx; ax_y1 <- ax_y0 + 0.16 * dy

cent <- sc_umap[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = Lineage]

p2b1 <- ggplot(sc_umap, aes(UMAP_1, UMAP_2, color = Lineage)) +
  ggrastr::rasterise(geom_point(size = 0.18, alpha = 0.62), dpi = 320) +
  annotate("segment", x = ax_x0, xend = ax_x1, y = ax_y0, yend = ax_y0,
           color = pal_cns$subink, linewidth = 0.36, arrow = arrow(length = unit(0.08, "cm"), type = "closed")) +
  annotate("segment", x = ax_x0, xend = ax_x0, y = ax_y0, yend = ax_y1,
           color = pal_cns$subink, linewidth = 0.36, arrow = arrow(length = unit(0.08, "cm"), type = "closed")) +
  annotate("text", x = (ax_x0 + ax_x1) / 2, y = ax_y0 - 0.042 * dy, label = "UMAP 1",
           size = 2.0, color = pal_cns$subink) +
  annotate("text", x = ax_x0 - 0.042 * dx, y = (ax_y0 + ax_y1) / 2, label = "UMAP 2",
           size = 2.0, color = pal_cns$subink, angle = 90) +
  ggrepel::geom_label_repel(
    data = cent, aes(label = Lineage, color = Lineage),
    fill = alpha("white", 0.82), label.size = 0,
    size = 2.0, fontface = "bold",
    label.padding = unit(0.08, "lines"), box.padding = 0.24, point.padding = 0.10,
    force = 6, force_pull = 0.3, max.overlaps = Inf,
    segment.size = 0.25, segment.color = pal_cns$muted, min.segment.length = 0.1, seed = 42
  ) +
  scale_color_manual(values = ct_pal, guide = "none") +
  coord_cartesian(xlim = c(xr[1] - 0.08 * dx, xr[2] + 0.04 * dx),
                  ylim = c(yr[1] - 0.08 * dy, yr[2] + 0.05 * dy), clip = "off") +
  labs(title = "Cell lineages (n = 100,064)") +
  theme_void(base_family = "sans") +
  theme(plot.title = element_text(face = "plain", size = 8.0, color = pal_cns$subink, hjust = 0.5, margin = margin(b = 2)),
        plot.margin = margin(2, 5, 3, 3))

make_umap_gene_cns <- function(gene, ttl_md, pal_cols) {
  d <- copy(sc_umap[, .(UMAP_1, UMAP_2, val = get(gene))])[order(val)]
  ggplot(d, aes(UMAP_1, UMAP_2, color = val)) +
    ggrastr::rasterise(geom_point(size = 0.18, alpha = ifelse(d$val > 0, 0.88, 0.20)), dpi = 320) +
    annotate("segment", x = ax_x0, xend = ax_x1, y = ax_y0, yend = ax_y0,
             color = pal_cns$subink, linewidth = 0.36, arrow = arrow(length = unit(0.08, "cm"), type = "closed")) +
    annotate("segment", x = ax_x0, xend = ax_x0, y = ax_y0, yend = ax_y1,
             color = pal_cns$subink, linewidth = 0.36, arrow = arrow(length = unit(0.08, "cm"), type = "closed")) +
    scale_color_gradientn(colors = pal_cols, name = NULL) +
    coord_cartesian(xlim = c(xr[1] - 0.06 * dx, xr[2] + 0.05 * dx),
                    ylim = c(yr[1] - 0.07 * dy, yr[2] + 0.04 * dy)) +
    labs(title = ttl_md) +
    theme_void(base_family = "sans") +
    theme(
      plot.title = ggtext::element_markdown(size = 8.2, color = pal_cns$ink, hjust = 0.5, margin = margin(b = 2)),
      legend.position = "right",
      legend.key.height = unit(0.40, "cm"),
      legend.key.width  = unit(0.16, "cm"),
      legend.text = element_text(size = 6.8, color = pal_cns$subink),
      legend.margin = margin(0, 3, 0, -2),
      plot.margin = margin(2, 4, 3, 4)
    )
}

p2b2 <- make_umap_gene_cns("CXCL8", "<i>CXCL8</i>", c("#E2E8F0", "#F39B7F", "#E64B35", "#991B1B"))
p2b3 <- make_umap_gene_cns("PPARG", "<i>PPARG</i>", c("#E2E8F0", "#4DBBD5", "#3C5488", "#1E2D52"))

fig2_b <- wrap_elements(full = (p2b1 | p2b2 | p2b3) +
  plot_annotation(
    title = cns_panel_title("B", "Cell-level UMAP (Wu et al., GSE176078)"),
    theme = theme(
      plot.title = ggtext::element_markdown(margin = margin(b = 3)),
      plot.margin = margin(2, 2, 0, 2)
    )
  ))

# ------------------------------------------------------------------------------
# Panel C: Donor-level subtype test (Wu, n = 26 donors)
# ------------------------------------------------------------------------------
lineages2c  <- c("Myeloid", "Endothelial", "Malignant", "ALL_CELLS")
genes2c     <- c("CXCL8", "PPARG")
subtypes2c  <- c("ER+", "HER2+", "TNBC")
min_cells2c <- 1L
keys2c      <- c("celltype", "gene", "subtype")

donor2c <- fread(file.path(res_dir, "planA_scrna_wu_donor_expression.csv"))[
  celltype %chin% lineages2c & gene %chin% genes2c & n_cells >= min_cells2c]
group2c <- fread(file.path(res_dir, "planA_scrna_wu_donor_subtype_summary.csv"))[
  celltype %chin% lineages2c & gene %chin% genes2c & endpoint == "mean_expr" & min_cells == min_cells2c]
tests2c <- fread(file.path(res_dir, "planA_scrna_wu_donor_subtype_tests.csv"))[
  celltype %chin% lineages2c & gene %chin% genes2c & endpoint == "mean_expr" & min_cells == min_cells2c]

counts2c <- merge(CJ(celltype = lineages2c, gene = genes2c, subtype = subtypes2c),
                  group2c[, c(keys2c, "n_donors"), with = FALSE], by = keys2c, all.x = TRUE)
counts2c[is.na(n_donors), n_donors := 0L]

row_bounds2c <- donor2c[, .(ymax = max(0.05, max(mean_expr))), by = gene]
counts2c <- merge(counts2c, row_bounds2c, by = "gene")
tests2c  <- merge(tests2c, row_bounds2c, by = "gene")
tests2c[, annotation := sprintf("italic(P)*' = %.3g  ('*italic(q)*' = %.3g)'", kw_pval, kw_fdr)]

fmt2c <- function(dt) {
  dt[, `:=`(
    celltype_lbl = factor(celltype, levels = lineages2c,
                          labels = c("Myeloid", "Endothelial", "Malignant", "All cells")),
    gene_md = factor(gene, levels = genes2c, labels = sprintf("<i>%s</i>", genes2c))
  )]
}
for (dt2c in list(donor2c, group2c, tests2c, counts2c)) fmt2c(dt2c)
for (dt2c in list(donor2c, group2c, counts2c)) dt2c[, subtype := factor(subtype, levels = subtypes2c)]

pad_2c <- rbind(
  row_bounds2c[, .(gene_md = factor(sprintf("<i>%s</i>", gene), levels = sprintf("<i>%s</i>", genes2c)),
                   celltype_lbl = factor("Myeloid", levels = c("Myeloid", "Endothelial", "Malignant", "All cells")),
                   subtype = factor("ER+", levels = subtypes2c), y_lim = -0.27 * ymax)],
  row_bounds2c[, .(gene_md = factor(sprintf("<i>%s</i>", gene), levels = sprintf("<i>%s</i>", genes2c)),
                   celltype_lbl = factor("Myeloid", levels = c("Myeloid", "Endothelial", "Malignant", "All cells")),
                   subtype = factor("ER+", levels = subtypes2c), y_lim = 1.24 * ymax)]
)

p2c <- ggplot(donor2c, aes(subtype, mean_expr)) +
  geom_blank(data = pad_2c, aes(x = subtype, y = y_lim), inherit.aes = FALSE) +
  geom_hline(yintercept = Inf, color = pal_cns$ref_line, linewidth = 0.35) +
  geom_hline(yintercept = 0, color = pal_cns$ref_line, linewidth = 0.30, linetype = "dashed") +
  geom_point(aes(color = subtype), position = position_jitter(width = 0.13, height = 0, seed = 26),
             size = 1.45, alpha = 0.82) +
  geom_errorbar(data = group2c, aes(x = subtype, ymin = donor_q1, ymax = donor_q3),
                inherit.aes = FALSE, width = 0.14, linewidth = 0.42, color = pal_cns$ink) +
  geom_point(data = group2c, aes(x = subtype, y = donor_median), inherit.aes = FALSE,
             shape = 23, size = 1.85, stroke = 0.45, fill = "white", color = pal_cns$ink) +
  geom_text(data = counts2c, aes(x = subtype, y = -0.14 * ymax, label = paste0("n=", n_donors)),
            inherit.aes = FALSE, size = 2.05, color = pal_cns$subink) +
  geom_text(data = tests2c, aes(x = 2, y = 1.15 * ymax, label = annotation),
            inherit.aes = FALSE, parse = TRUE, size = 2.15, color = pal_cns$ink) +
  facet_grid(gene_md ~ celltype_lbl, scales = "free_y") +
  scale_x_discrete(limits = subtypes2c, drop = FALSE) +
  scale_y_continuous(breaks = function(limits) { b <- scales::breaks_pretty(n = 3)(c(0, limits[2])); b[b >= 0] },
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_color_manual(values = c("ER+" = pal_cns$pparg, "HER2+" = pal_cns$teal, "TNBC" = pal_cns$cxcl8)) +
  labs(title = cns_panel_title("C", "Donor-level expression by subtype (Wu et al., n = 26 donors)"),
       x = NULL, y = "Donor mean\nnormalized expr.") +
  theme_cns_grid(8.2) +
  theme(
    legend.position = "none",
    axis.line       = element_line(color = pal_cns$ink, linewidth = 0.36),
    strip.text.x    = ggtext::element_markdown(size = 8.0, color = pal_cns$ink, margin = margin(t = 1, b = 3)),
    strip.text.y    = ggtext::element_markdown(size = 8.4, color = pal_cns$ink, margin = margin(l = 3, r = 2)),
    axis.text.x     = element_text(size = 7.4, color = pal_cns$ink),
    panel.spacing.x = unit(3.5, "mm"),
    panel.spacing.y = unit(3.0, "mm")
  )

# ------------------------------------------------------------------------------
# Panel D: Clean HPA IHC panels (italic gene symbols in titles)
# ------------------------------------------------------------------------------
make_hpa_panel_cns <- function(img_path, crop_box, title_md, short_sub,
                               core_px_width, accent_col = pal_cns$pparg) {
  mg_full <- magick::image_read(img_path)
  info <- magick::image_info(mg_full)
  W <- info$width; H <- info$height

  mg_main <- magick::image_scale(mg_full, "960x960")
  tmp1 <- tempfile(fileext = ".jpg")
  magick::image_write(mg_main, tmp1, format = "jpeg", quality = 92)
  arr_main <- jpeg::readJPEG(tmp1); unlink(tmp1)

  x0 <- round(crop_box[1] * W); cw <- round((crop_box[2] - crop_box[1]) * W)
  y0 <- round(crop_box[3] * H); ch <- round((crop_box[4] - crop_box[3]) * H)
  mg_zoom <- magick::image_crop(mg_full, sprintf("%dx%d+%d+%d", cw, ch, x0, y0))
  mg_zoom <- magick::image_scale(mg_zoom, "720x720")
  tmp2 <- tempfile(fileext = ".jpg")
  magick::image_write(mg_zoom, tmp2, format = "jpeg", quality = 94)
  arr_zoom <- jpeg::readJPEG(tmp2); unlink(tmp2)

  g_main <- grid::rasterGrob(arr_main, interpolate = TRUE, width = unit(1, "npc"), height = unit(1, "npc"))
  g_zoom <- grid::rasterGrob(arr_zoom, interpolate = TRUE, width = unit(1, "npc"), height = unit(1, "npc"))

  bx1 <- crop_box[1]; bx2 <- crop_box[2]
  by1 <- 1 - crop_box[4]; by2 <- 1 - crop_box[3]
  gx1 <- 0.01; gx2 <- 0.96; gy1 <- 0.04; gy2 <- 0.99
  map_x <- function(x) gx1 + x * (gx2 - gx1)
  map_y <- function(y) gy1 + y * (gy2 - gy1)
  zx1 <- 1.02; zx2 <- 1.64; zy1 <- 0.20; zy2 <- 0.82

  sb_frac_img <- (0.20 * core_px_width) / W
  sb_len_plot <- sb_frac_img * (gx2 - gx1)
  sb_x0 <- gx1 + 0.045; sb_x1 <- sb_x0 + sb_len_plot; sb_y <- gy1 + 0.065

  ggplot() +
    annotation_custom(g_main, xmin = gx1, xmax = gx2, ymin = gy1, ymax = gy2) +
    annotation_custom(g_zoom, xmin = zx1, xmax = zx2, ymin = zy1, ymax = zy2) +
    annotate("rect", xmin = map_x(bx1), xmax = map_x(bx2), ymin = map_y(by1), ymax = map_y(by2),
             fill = NA, color = accent_col, linewidth = 0.55) +
    annotate("segment", x = map_x(bx2), xend = zx1, y = map_y(by2), yend = zy2,
             color = accent_col, linewidth = 0.38, linetype = "dashed") +
    annotate("segment", x = map_x(bx2), xend = zx1, y = map_y(by1), yend = zy1,
             color = accent_col, linewidth = 0.38, linetype = "dashed") +
    annotate("rect", xmin = zx1, xmax = zx2, ymin = zy1, ymax = zy2,
             fill = NA, color = accent_col, linewidth = 0.65) +
    annotate("rect", xmin = sb_x0 - 0.02, xmax = sb_x1 + 0.02, ymin = sb_y - 0.022, ymax = sb_y + 0.082,
             fill = alpha("white", 0.88), color = NA) +
    annotate("segment", x = sb_x0, xend = sb_x1, y = sb_y, yend = sb_y,
             color = pal_cns$ink, linewidth = 1.1) +
    annotate("text", x = (sb_x0 + sb_x1) / 2, y = sb_y + 0.044, label = "~200 \u00b5m",
             size = 1.95, fontface = "bold", color = pal_cns$ink) +
    coord_fixed(ratio = 1, xlim = c(0, 1.66), ylim = c(0.02, 1.01), expand = FALSE) +
    labs(title = title_md, caption = short_sub) +
    theme_void(base_family = "sans") +
    theme(
      plot.title   = ggtext::element_markdown(size = 7.9, color = pal_cns$ink, hjust = 0, margin = margin(b = 2)),
      plot.caption = element_text(size = 6.9, color = pal_cns$subink, hjust = 0, margin = margin(t = 2)),
      plot.margin  = margin(2, 4, 2, 4)
    )
}

p2d1 <- make_hpa_panel_cns(
  file.path(ihc_dir, "HPA_PPARG_NormalBreast_IHC.jpg"),
  crop_box = c(0.54, 0.80, 0.68, 0.94),
  title_md = "<i>PPARG</i> \u2014 Normal breast",
  short_sub = "CAB004282 \u00b7 Adipocytes low",
  core_px_width = 2351, accent_col = pal_cns$pparg
)
p2d2 <- make_hpa_panel_cns(
  file.path(ihc_dir, "HPA_PPARG_BreastCancer_IHC.jpg"),
  crop_box = c(0.38, 0.64, 0.48, 0.74),
  title_md = "<i>PPARG</i> \u2014 Duct carcinoma",
  short_sub = "CAB004282 \u00b7 4 med, 6 low, 2 ND",
  core_px_width = 2556, accent_col = pal_cns$pparg
)
p2d3 <- make_hpa_panel_cns(
  file.path(ihc_dir, "HPA_CXCL8_BreastCancer_IHC.jpg"),
  crop_box = c(0.45, 0.71, 0.61, 0.87),
  title_md = "<i>CXCL8</i> \u2014 Duct carcinoma",
  short_sub = "HPA057179 \u00b7 12/12 not detected",
  core_px_width = 2405, accent_col = pal_cns$cxcl8
)

fig2_d <- wrap_elements(full = (p2d1 | p2d2 | p2d3) +
  plot_annotation(
    title = cns_panel_title("D", "Illustrative HPA breast immunohistochemistry"),
    theme = theme(
      plot.title = ggtext::element_markdown(margin = margin(b = 3)),
      plot.margin = margin(2, 2, 1, 2)
    )
  ))

fig2_cns <- wrap_plots(
  wrap_elements(full = p2a),
  fig2_b,
  wrap_elements(full = p2c),
  fig2_d,
  ncol = 1,
  heights = c(50, 56, 66, 52)
)

save_cns_figure(fig2_cns, "Fig2_Tumor_scRNA_and_HPA_Pathology", width_mm = 180, height_mm = 224, dpi_png = 350)
cat("CNS-Polished Fig 2 complete.\n")
