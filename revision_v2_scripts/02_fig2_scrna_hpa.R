# 01a_fig6_new.R
# Figure 6 Redesign (Clean Top-Journal Minimalist Typography):
# Zero prose subtitles or multi-line metadata footers inside the graphic;
# technical provenance, crop coordinates, and CC BY-SA 3.0 details belong in Figure Legend.

local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "05_revision_v2", "01_scripts")
  helper_path <- file.path(script_dir, "00_rev2_helpers.R")
  if (!file.exists(helper_path)) helper_path <- file.path("05_revision_v2", "01_scripts", "00_rev2_helpers.R")
  source(helper_path, local = FALSE)
})

cat("=== Generating Clean Minimalist Fig 6 ===\n")

sc_cache <- readRDS(file.path(sc_dir, "planA_scrna_plot_cache.rds"))
sc_sum   <- fread(file.path(res_dir, "planA_scrna_celltype_summary.csv"))
g14_sens <- fread(file.path(res_dir, "planA_scrna_gse148673_ecis05t_sensitivity.csv"))

# ------------------------------------------------------------------------------
# Panel A (First Sightline): Cross-atlas single-cell lineage localization
# ------------------------------------------------------------------------------
shared_ct    <- c("Mono/Macro", "Endothelial", "Fibroblasts", "Malignant")
lineage_lvls <- c("Myeloid", "Endothelial", "Fibroblasts", "Malignant")
genes_6a     <- c("CXCL8", "IL1B", "PPARG", "FABP4", "PLIN1")
cohorts_6a   <- c("Wu (GSE176078)", "E-MTAB-8107", "GSE148673*")

wu_em <- sc_sum[dataset %chin% c("BRCA_GSE176078", "BRCA_EMTAB8107") &
                  celltype_major %chin% shared_ct & gene %chin% genes_6a,
                .(atlas = fifelse(dataset == "BRCA_GSE176078", "Wu (GSE176078)", "E-MTAB-8107"),
                  celltype = celltype_major, gene, pct_expr, mean_expr)]

g14_clean <- g14_sens[celltype %chin% shared_ct & gene %chin% genes_6a,
                      .(atlas = "GSE148673*", celltype, gene,
                        pct_expr = pct_expr_exclude_ECIS05T,
                        mean_expr = mean_expr_exclude_ECIS05T)]

grid_6a <- CJ(atlas = cohorts_6a, celltype = shared_ct, gene = genes_6a)
dt_6a <- merge(grid_6a, rbind(wu_em, g14_clean), by = c("atlas", "celltype", "gene"), all.x = TRUE)
dt_6a[, atlas := factor(atlas, levels = cohorts_6a)]
dt_6a[, gene := factor(gene, levels = genes_6a)]
dt_6a[, Lineage := factor(celltype, levels = rev(shared_ct), labels = rev(lineage_lvls))]
dt_6a[, missing_gene := is.na(pct_expr)]
dt_6a[, pct_lbl := fifelse(!missing_gene & pct_expr >= 12, sprintf("%.0f%%", pct_expr), "")]
dt_6a[, is_big := !missing_gene & pct_expr >= 45]

hl_6a <- data.table(
  gene = factor(c("CXCL8", "PPARG"), levels = genes_6a),
  ymin = c(3.55, 2.55), ymax = c(4.45, 3.45),
  fill_col = c(pal_rev$cxcl8_soft, pal_rev$pparg_soft)
)

p6a <- ggplot(dt_6a, aes(x = atlas, y = Lineage)) +
  geom_rect(data = hl_6a, aes(ymin = ymin, ymax = ymax, fill = fill_col),
            xmin = -Inf, xmax = Inf, inherit.aes = FALSE, alpha = 0.55) +
  scale_fill_identity() +
  geom_tile(data = dt_6a[missing_gene == TRUE], width = 0.84, height = 0.78,
            fill = "#F1F5F9", color = pal_rev$border, linewidth = 0.35) +
  geom_text(data = dt_6a[missing_gene == TRUE], label = "NP", size = 2.35,
            fontface = "bold", color = pal_rev$slate) +
  geom_point(data = dt_6a[missing_gene == FALSE],
             aes(size = pct_expr, color = mean_expr)) +
  # Put % cleanly below normal dots, or white-centered inside very large dark dots (>=45%)
  geom_text(data = dt_6a[nzchar(pct_lbl) & !is_big],
            aes(label = pct_lbl), nudge_y = -0.37, size = 2.05, fontface = "bold",
            color = pal_rev$ink, show.legend = FALSE) +
  geom_text(data = dt_6a[nzchar(pct_lbl) & is_big],
            aes(label = pct_lbl), size = 2.05, fontface = "bold",
            color = "white", show.legend = FALSE) +
  facet_wrap(~ gene, nrow = 1) +
  scale_x_discrete(limits = cohorts_6a) +
  scale_y_discrete(limits = rev(lineage_lvls), expand = expansion(add = c(0.65, 0.55))) +
  scale_size_area(max_size = 6.4, breaks = c(10, 25, 45), name = "% cells") +
  scale_color_gradientn(colors = c("#E2E8F0", "#60A5FA", "#F59E0B", "#B91C1C"),
                        name = "Mean expr.") +
  labs(title = "A   Cross-atlas single-cell localization", x = NULL, y = NULL) +
  theme_rev_grid(8.6) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 7.3),
    axis.text.y = element_text(size = 7.8, face = "bold"),
    strip.text  = element_text(face = "bold.italic", size = 8.6),
    panel.spacing.x = unit(2.5, "mm"),
    legend.position = "right",
    legend.box = "vertical",
    legend.key.height = unit(3.4, "mm"),
    legend.key.width  = unit(2.8, "mm"),
    legend.margin = margin(0, 0, 0, 2)
  )

# ------------------------------------------------------------------------------
# Panel B: Cell-level Wu UMAPs (Lineages, CXCL8, PPARG)
# ------------------------------------------------------------------------------
sc_df <- copy(sc_cache[["BRCA_GSE176078"]])
sc_umap <- sc_df[!is.na(UMAP_1) & !is.na(UMAP_2)]
sc_umap[, Lineage := clean_lineage_rev(Celltype_major)]

ct_levels <- sc_umap[, .N, by = Lineage][order(-N)]$Lineage
ct_pal <- c(
  "T Cells"          = "#3B5B8C",
  "Malignant"        = "#DC2626",
  "Myeloid"          = "#0F766E",
  "Endothelial"      = "#0284C7",
  "Fibroblasts"      = "#D97706",
  "PVL"              = "#64748B",
  "Norm. Epithelial" = "#EA580C",
  "Plasma"           = "#DB2777",
  "B Cells"          = "#7C3AED"
)
ct_pal <- ct_pal[intersect(ct_levels, names(ct_pal))]

xr <- range(sc_umap$UMAP_1); yr <- range(sc_umap$UMAP_2)
dx <- diff(xr); dy <- diff(yr)
ax_x0 <- xr[1] - 0.03 * dx; ax_y0 <- yr[1] - 0.03 * dy
ax_x1 <- ax_x0 + 0.16 * dx; ax_y1 <- ax_y0 + 0.16 * dy

cent <- sc_umap[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = Lineage]

p6b1 <- ggplot(sc_umap, aes(UMAP_1, UMAP_2, color = Lineage)) +
  ggrastr::rasterise(geom_point(size = 0.20, alpha = 0.62), dpi = 320) +
  annotate("segment", x = ax_x0, xend = ax_x1, y = ax_y0, yend = ax_y0,
           color = pal_rev$ink, linewidth = 0.45, arrow = arrow(length = unit(0.09, "cm"), type = "closed")) +
  annotate("segment", x = ax_x0, xend = ax_x0, y = ax_y0, yend = ax_y1,
           color = pal_rev$ink, linewidth = 0.45, arrow = arrow(length = unit(0.09, "cm"), type = "closed")) +
  annotate("text", x = (ax_x0 + ax_x1) / 2, y = ax_y0 - 0.04 * dy, label = "UMAP 1", size = 2.1, fontface = "bold", color = pal_rev$subink) +
  annotate("text", x = ax_x0 - 0.04 * dx, y = (ax_y0 + ax_y1) / 2, label = "UMAP 2", size = 2.1, fontface = "bold", color = pal_rev$subink, angle = 90) +
  ggrepel::geom_label_repel(
    data = cent, aes(label = Lineage, color = Lineage),
    fill = alpha("white", 0.93), size = 2.0, fontface = "bold",
    label.padding = unit(0.11, "lines"), box.padding = 0.28, point.padding = 0.12,
    force = 6, force_pull = 0.3, max.overlaps = Inf,
    segment.size = 0.30, min.segment.length = 0.1, seed = 42
  ) +
  scale_color_manual(values = ct_pal, guide = "none") +
  coord_cartesian(xlim = c(xr[1] - 0.08 * dx, xr[2] + 0.04 * dx),
                  ylim = c(yr[1] - 0.08 * dy, yr[2] + 0.05 * dy), clip = "off") +
  labs(title = "Cell lineages (n = 100,064)") +
  theme_void(base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 8.2, color = pal_rev$subink, hjust = 0.5, margin = margin(b = 2)),
        plot.margin = margin(2, 5, 3, 3))

make_umap_gene <- function(gene, ttl, pal_mid, pal_hi) {
  d <- copy(sc_umap[, .(UMAP_1, UMAP_2, val = get(gene))])[order(val)]
  ggplot(d, aes(UMAP_1, UMAP_2, color = val)) +
    ggrastr::rasterise(geom_point(size = 0.20, alpha = ifelse(d$val > 0, 0.88, 0.22)), dpi = 320) +
    annotate("segment", x = ax_x0, xend = ax_x1, y = ax_y0, yend = ax_y0,
             color = pal_rev$ink, linewidth = 0.45, arrow = arrow(length = unit(0.09, "cm"), type = "closed")) +
    annotate("segment", x = ax_x0, xend = ax_x0, y = ax_y0, yend = ax_y1,
             color = pal_rev$ink, linewidth = 0.45, arrow = arrow(length = unit(0.09, "cm"), type = "closed")) +
    scale_color_gradientn(colors = c("#E2E8F0", pal_mid, pal_hi), name = NULL) +
    coord_cartesian(xlim = c(xr[1] - 0.06 * dx, xr[2] + 0.05 * dx),
                    ylim = c(yr[1] - 0.07 * dy, yr[2] + 0.04 * dy)) +
    labs(title = ttl) +
    theme_void(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold.italic", size = 8.4, color = pal_rev$subink, hjust = 0.5, margin = margin(b = 2)),
      legend.position = "right",
      legend.key.height = unit(0.42, "cm"),
      legend.key.width  = unit(0.18, "cm"),
      legend.text = element_text(size = 7.0, color = pal_rev$subink),
      legend.margin = margin(0, 3, 0, -2),
      plot.margin = margin(2, 4, 3, 4)
    )
}

p6b2 <- make_umap_gene("CXCL8", "CXCL8", "#F59E0B", "#9F1239")
p6b3 <- make_umap_gene("PPARG", "PPARG", "#60A5FA", "#1E3A8A")

fig6_b <- wrap_elements(full = (p6b1 | p6b2 | p6b3) +
  plot_annotation(
    title = "B   Cell-level UMAP (Wu et al., GSE176078)",
    theme = theme(plot.title = element_text(face = "bold", size = 9.6, color = pal_rev$ink, margin = margin(b = 2)),
                  plot.margin = margin(2, 2, 0, 2))
  ))

# ------------------------------------------------------------------------------
# Panel C: Donor-level subtype test (Wu, n = 26 donors)
# ------------------------------------------------------------------------------
lineages6c  <- c("Myeloid", "Endothelial", "Malignant", "ALL_CELLS")
genes6c     <- c("CXCL8", "PPARG")
subtypes6c  <- c("ER+", "HER2+", "TNBC")
min_cells6c <- 1L
keys6c      <- c("celltype", "gene", "subtype")

donor6c <- fread(file.path(res_dir, "planA_scrna_wu_donor_expression.csv"))[
  celltype %chin% lineages6c & gene %chin% genes6c & n_cells >= min_cells6c]
group6c <- fread(file.path(res_dir, "planA_scrna_wu_donor_subtype_summary.csv"))[
  celltype %chin% lineages6c & gene %chin% genes6c & endpoint == "mean_expr" & min_cells == min_cells6c]
tests6c <- fread(file.path(res_dir, "planA_scrna_wu_donor_subtype_tests.csv"))[
  celltype %chin% lineages6c & gene %chin% genes6c & endpoint == "mean_expr" & min_cells == min_cells6c]

counts6c <- merge(CJ(celltype = lineages6c, gene = genes6c, subtype = subtypes6c),
                  group6c[, c(keys6c, "n_donors"), with = FALSE], by = keys6c, all.x = TRUE)
counts6c[is.na(n_donors), n_donors := 0L]

row_bounds6c <- donor6c[, .(ymax = max(0.05, max(mean_expr))), by = gene]
counts6c <- merge(counts6c, row_bounds6c, by = "gene")
tests6c  <- merge(tests6c, row_bounds6c, by = "gene")
tests6c[, annotation := sprintf("P = %.3g  (q = %.3g)", kw_pval, kw_fdr)]

fmt6c <- function(dt) {
  dt[, `:=`(
    celltype_lbl = factor(celltype, levels = lineages6c,
                          labels = c("Myeloid", "Endothelial", "Malignant", "All cells")),
    gene = factor(gene, levels = genes6c)
  )]
}
for (dt6c in list(donor6c, group6c, tests6c, counts6c)) fmt6c(dt6c)
for (dt6c in list(donor6c, group6c, counts6c)) dt6c[, subtype := factor(subtype, levels = subtypes6c)]

pad_6c <- rbind(
  row_bounds6c[, .(gene, celltype_lbl = factor("Myeloid", levels = c("Myeloid", "Endothelial", "Malignant", "All cells")),
                   subtype = factor("ER+", levels = subtypes6c), y_lim = -0.27 * ymax)],
  row_bounds6c[, .(gene, celltype_lbl = factor("Myeloid", levels = c("Myeloid", "Endothelial", "Malignant", "All cells")),
                   subtype = factor("ER+", levels = subtypes6c), y_lim = 1.24 * ymax)]
)

p6c <- ggplot(donor6c, aes(subtype, mean_expr)) +
  geom_blank(data = pad_6c, aes(x = subtype, y = y_lim), inherit.aes = FALSE) +
  geom_hline(yintercept = 0, color = "#CBD5E1", linewidth = 0.32) +
  geom_point(aes(color = subtype), position = position_jitter(width = 0.13, height = 0, seed = 26),
             size = 1.5, alpha = 0.80) +
  geom_errorbar(data = group6c, aes(x = subtype, ymin = donor_q1, ymax = donor_q3),
                inherit.aes = FALSE, width = 0.15, linewidth = 0.50, color = pal_rev$ink) +
  geom_point(data = group6c, aes(x = subtype, y = donor_median), inherit.aes = FALSE,
             shape = 23, size = 1.95, stroke = 0.52, fill = "white", color = pal_rev$ink) +
  geom_text(data = counts6c, aes(x = subtype, y = -0.14 * ymax, label = paste0("n=", n_donors)),
            inherit.aes = FALSE, size = 2.2, color = pal_rev$subink) +
  geom_text(data = tests6c, aes(x = 2, y = 1.15 * ymax, label = annotation),
            inherit.aes = FALSE, size = 2.25, fontface = "bold", color = pal_rev$ink) +
  facet_grid(gene ~ celltype_lbl, scales = "free_y") +
  scale_x_discrete(limits = subtypes6c, drop = FALSE) +
  scale_y_continuous(breaks = function(limits) { b <- scales::breaks_pretty(n = 3)(c(0, limits[2])); b[b >= 0] },
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_color_manual(values = c("ER+" = pal_rev$pparg, "HER2+" = pal_rev$amber, "TNBC" = pal_rev$cxcl8)) +
  labs(title = "C   Donor-level expression by subtype (Wu et al., n = 26 donors)",
       x = NULL, y = "Donor mean\nnormalized expr.") +
  theme_rev_grid(8.6) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    strip.text.y = element_text(face = "bold.italic", size = 8.6),
    axis.text.x = element_text(size = 7.6, face = "bold"),
    panel.spacing = unit(2.5, "mm")
  )

# ------------------------------------------------------------------------------
# Panel D: Clean HPA IHC panels (1-line minimal caption; full provenance in Legend)
# ------------------------------------------------------------------------------
make_hpa_panel_clean <- function(img_path, crop_box, title_lbl, short_sub,
                                 core_px_width, accent_col = "#2B5B84") {
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
             fill = NA, color = accent_col, linewidth = 0.78) +
    annotate("segment", x = map_x(bx2), xend = zx1, y = map_y(by2), yend = zy2,
             color = accent_col, linewidth = 0.45, linetype = "dashed") +
    annotate("segment", x = map_x(bx2), xend = zx1, y = map_y(by1), yend = zy1,
             color = accent_col, linewidth = 0.45, linetype = "dashed") +
    annotate("rect", xmin = zx1, xmax = zx2, ymin = zy1, ymax = zy2,
             fill = NA, color = accent_col, linewidth = 0.95) +
    annotate("rect", xmin = sb_x0 - 0.02, xmax = sb_x1 + 0.02, ymin = sb_y - 0.022, ymax = sb_y + 0.082,
             fill = alpha("white", 0.90), color = NA) +
    annotate("segment", x = sb_x0, xend = sb_x1, y = sb_y, yend = sb_y,
             color = pal_rev$ink, linewidth = 1.25) +
    annotate("text", x = (sb_x0 + sb_x1) / 2, y = sb_y + 0.044, label = "~200 µm",
             size = 2.0, fontface = "bold", color = pal_rev$ink) +
    coord_fixed(ratio = 1, xlim = c(0, 1.66), ylim = c(0.02, 1.01), expand = FALSE) +
    labs(title = title_lbl, caption = short_sub) +
    theme_void(base_family = "sans") +
    theme(
      plot.title   = element_text(face = "bold", size = 8.1, color = pal_rev$ink, hjust = 0, margin = margin(b = 2)),
      plot.caption = element_text(size = 7.2, color = pal_rev$subink, hjust = 0, margin = margin(t = 2)),
      plot.margin  = margin(2, 4, 2, 4)
    )
}

p6d1 <- make_hpa_panel_clean(
  file.path(ihc_dir, "HPA_PPARG_NormalBreast_IHC.jpg"),
  crop_box = c(0.54, 0.80, 0.68, 0.94),
  title_lbl = "PPARG — Normal breast",
  short_sub = "CAB004282 · Adipocytes low",
  core_px_width = 2351, accent_col = pal_rev$pparg
)
p6d2 <- make_hpa_panel_clean(
  file.path(ihc_dir, "HPA_PPARG_BreastCancer_IHC.jpg"),
  crop_box = c(0.38, 0.64, 0.48, 0.74),
  title_lbl = "PPARG — Duct carcinoma",
  short_sub = "CAB004282 · 4 med, 6 low, 2 ND",
  core_px_width = 2556, accent_col = pal_rev$pparg
)
p6d3 <- make_hpa_panel_clean(
  file.path(ihc_dir, "HPA_CXCL8_BreastCancer_IHC.jpg"),
  crop_box = c(0.45, 0.71, 0.61, 0.87),
  title_lbl = "CXCL8 — Duct carcinoma",
  short_sub = "HPA057179 · 12/12 not detected",
  core_px_width = 2405, accent_col = pal_rev$cxcl8
)

fig6_d <- wrap_elements(full = (p6d1 | p6d2 | p6d3) +
  plot_annotation(
    title = "D   Illustrative HPA breast immunohistochemistry",
    theme = theme(
      plot.title = element_text(face = "bold", size = 9.6, color = pal_rev$ink, margin = margin(b = 2)),
      plot.margin = margin(2, 2, 1, 2)
    )
  ))

fig6_new <- wrap_plots(
  wrap_elements(full = p6a),
  fig6_b,
  wrap_elements(full = p6c),
  fig6_d,
  ncol = 1,
  heights = c(50, 56, 66, 52)
)
save_rev_figure(fig6_new, "Fig2_Tumor_scRNA_and_HPA_Pathology", width_mm = 180, height_mm = 224, dpi_png = 350)
cat("Clean Minimalist Fig 6 complete.\n")
