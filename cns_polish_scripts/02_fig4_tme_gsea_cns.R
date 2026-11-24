# 02_fig4_tme_gsea_cns.R
# Redesigned True CNS-Grade Figure 4 (Nature Medicine / Nature Immunology / Cancer Cell layout):
#   A (Top, Full-Width)   : Horizontal 4-row x 17-column TME Correlation Dot-Heatmap with 3-group top lineage brackets
#   B (Middle, Left/Right): MSigDB Hallmark Enrichment (*CXCL8* ER- vs *PPARG* ER+) with single clean shared legend
#   C (Bottom, Full-Width): Crisp Multi-Curve GSEA Running Enrichment Plot (no muddy fill) + docked 4-track barcode rug

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggtext)
  library(ragg)
  library(cowplot)
  library(grid)
  library(scales)
  library(enrichplot)
})

root_dir    <- Sys.getenv("BRCA_PROJECT_ROOT", unset = ".")
bulk_dir    <- file.path(root_dir, "00_data/bulk_cohorts")
res_dir     <- file.path(root_dir, "02_results")
cns_fig_dir <- file.path(root_dir, "06_figures_cns_polish/03_figures")
dir.create(cns_fig_dir, recursive = TRUE, showWarnings = FALSE)

pal_cns <- list(
  cxcl8      = "#E64B35",
  cxcl8_soft = "#FDE8E5",
  pparg      = "#3C5488",
  pparg_soft = "#E6F4F8",
  teal       = "#00A087",
  ink        = "#0F172A",
  subink     = "#334155"
)

cns_panel_title <- function(letter, subtitle_md) {
  sprintf("<span style='font-size:11.5pt; color:#0F172A'><b>%s</b></span>\u2003<span style='font-size:8.6pt; color:#1E293B'>%s</span>",
          letter, subtitle_md)
}

theme_cns_fig4 <- function(base_size = 8.0, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line           = element_line(color = pal_cns$ink, linewidth = 0.36),
      axis.ticks          = element_line(color = pal_cns$subink, linewidth = 0.32),
      axis.ticks.length   = unit(2.0, "pt"),
      axis.title          = element_text(face = "plain", size = base_size - 0.1, color = pal_cns$ink),
      axis.text           = element_text(size = base_size - 0.8, color = pal_cns$subink),
      panel.grid          = element_blank(),
      panel.background    = element_rect(fill = "white", color = NA),
      plot.background     = element_rect(fill = "white", color = NA),
      plot.title          = element_markdown(size = base_size + 0.6, color = pal_cns$ink, hjust = 0, margin = margin(b = 4)),
      plot.title.position = "plot",
      legend.background   = element_blank(),
      legend.key          = element_blank()
    )
}

save_cns_fig4 <- function(p, name, width_mm = 185, height_mm = 210, dpi_png = 350) {
  png_path <- file.path(cns_fig_dir, paste0(name, ".png"))
  pdf_path <- file.path(cns_fig_dir, paste0(name, ".pdf"))
  ggsave(png_path, p, width = width_mm / 25.4, height = height_mm / 25.4, dpi = dpi_png, device = ragg::agg_png, bg = "white")
  ggsave(pdf_path, p, width = width_mm / 25.4, height = height_mm / 25.4, device = cairo_pdf, bg = "white")
  cat(sprintf("Saved %s | %.0fx%.0f mm | PNG=%.2f MB | PDF=%.2f MB\n",
              name, width_mm, height_mm,
              file.info(png_path)$size / 1048576,
              file.info(pdf_path)$size / 1048576))
}

cat("=== Generating Redesigned CNS-Grade Fig 4 ===\n")

# ---------------------------------------------------------------------------
# 1. Panel A: Horizontal 4 x 17 TME Correlation Dot-Heatmap (TCGA-BRCA)
# ---------------------------------------------------------------------------
tme_map <- c(
  "quanTIseq__B cell"                  = "B cells (quanTIseq)",
  "quanTIseq__Macrophage M1"           = "M1 macrophages (quanTIseq)",
  "MCPcounter__Neutrophil"             = "Neutrophils (MCP)",
  "quanTIseq__NK cell"                 = "NK cells (quanTIseq)",
  "TMB_NONSYNONYMOUS"                  = "Nonsynonymous TMB",
  "MCPcounter__T cell CD8+"            = "CD8+ T cells (MCP)",
  "quanTIseq__Neutrophil"              = "Neutrophils (quanTIseq)",
  "quanTIseq__T cell CD8+"             = "CD8+ T cells (quanTIseq)",
  "MCPcounter__T cell"                 = "T cells (MCP)",
  "MCPcounter__B cell"                 = "B cells (MCP)",
  "MCPcounter__Endothelial cell"       = "Endothelial cells (MCP)",
  "ESTIMATE__stroma score"             = "Stromal score (ESTIMATE)",
  "ESTIMATE__tumor purity"             = "Tumor purity (ESTIMATE)",
  "ESTIMATE__estimate score"           = "ESTIMATE score",
  "MCPcounter__cytotoxicity score"     = "Cytotoxicity (MCP)",
  "quanTIseq__uncharacterized cell"    = "Uncharacterized (quanTIseq)",
  "MCPcounter__Myeloid dendritic cell" = "Myeloid DCs (MCP)"
)

lineage_map <- c(
  "quanTIseq__Macrophage M1"           = "Myeloid / Innate",
  "MCPcounter__Neutrophil"             = "Myeloid / Innate",
  "quanTIseq__Neutrophil"              = "Myeloid / Innate",
  "MCPcounter__Myeloid dendritic cell" = "Myeloid / Innate",
  "quanTIseq__NK cell"                 = "Myeloid / Innate",
  "MCPcounter__cytotoxicity score"     = "Lymphoid",
  "MCPcounter__T cell"                 = "Lymphoid",
  "MCPcounter__T cell CD8+"            = "Lymphoid",
  "quanTIseq__T cell CD8+"             = "Lymphoid",
  "MCPcounter__B cell"                 = "Lymphoid",
  "quanTIseq__B cell"                  = "Lymphoid",
  "MCPcounter__Endothelial cell"       = "Stromal / Genomic",
  "ESTIMATE__stroma score"             = "Stromal / Genomic",
  "ESTIMATE__estimate score"           = "Stromal / Genomic",
  "ESTIMATE__tumor purity"             = "Stromal / Genomic",
  "quanTIseq__uncharacterized cell"    = "Stromal / Genomic",
  "TMB_NONSYNONYMOUS"                  = "Stromal / Genomic"
)

cor5     <- fread(file.path(res_dir, "module5_tme_checkpoint_correlations.csv"))
strata7a <- fread(file.path(res_dir, "module5_tcga_er_sample_counts.csv"))[er_status %chin% c("Negative", "Positive")]
stopifnot(nrow(strata7a) == 2L, !anyDuplicated(strata7a$er_status))
strata7a[, stratum := ifelse(er_status == "Negative", "ER_neg", "ER_pos")]
cor5 <- merge(cor5, strata7a[, .(stratum, n_stratum = n_samples)], by = "stratum")

# Exact top-10 significant features for focal strata (100% consistent with v2)
top_feats <- cor5[((gene_symbol == "CXCL8" & stratum == "ER_neg") | (gene_symbol == "PPARG" & stratum == "ER_pos")) &
                    feature_type %chin% c("TME_Cell_Score", "Checkpoint", "Immune_Checkpoint", "Genomic_Feature") &
                    FDR < 0.05][order(gene_symbol, -abs(rho))][, head(.SD, 10L), by = gene_symbol]
feat_union <- unique(top_feats$feature)

hm_dt <- cor5[gene_symbol %chin% c("CXCL8", "PPARG") & stratum %chin% c("ER_neg", "ER_pos") & feature %chin% feat_union]
hm_dt[, feat_clean := ifelse(feature %chin% names(tme_map), unname(tme_map[feature]), gsub("_", " ", feature))]
hm_dt[, lineage := factor(unname(lineage_map[feature]), levels = c("Myeloid / Innate", "Lymphoid", "Stromal / Genomic"))]

# Order features cleanly within each lineage
ord_helper <- dcast(hm_dt, feature + feat_clean + lineage ~ gene_symbol + stratum, value.var = "rho")
ord_helper[, sort_key := ifelse(lineage == "Stromal / Genomic", PPARG_ER_pos, CXCL8_ER_neg)]
ord_helper <- ord_helper[order(lineage, -sort_key)]
feat_levels <- ord_helper$feat_clean

hm_dt[, feat_f := factor(feat_clean, levels = feat_levels)]
hm_dt[, x_num  := as.numeric(feat_f)]

# Rows (top to bottom): CXCL8 ER-, CXCL8 ER+, PPARG ER-, PPARG ER+
row_levels <- c("PPARG_ER_pos", "PPARG_ER_neg", "CXCL8_ER_pos", "CXCL8_ER_neg")
row_labels <- c(
  "CXCL8_ER_neg" = "<b>*CXCL8*</b> | ER\u2212 (n=237)",
  "CXCL8_ER_pos" = "*CXCL8* | ER+ (n=794)",
  "PPARG_ER_neg" = "*PPARG* | ER\u2212 (n=237)",
  "PPARG_ER_pos" = "<b>*PPARG*</b> | ER+ (n=794)"
)
hm_dt[, row_key := paste0(gene_symbol, "_", stratum)]
hm_dt[, row_f   := factor(row_key, levels = row_levels)]
hm_dt[, y_num   := as.numeric(row_f)]
hm_dt[, sig     := fifelse(FDR < 0.001, "***", fifelse(FDR < 0.01, "**", fifelse(FDR < 0.05, "*", "")))]
hm_dt[, is_focal := (row_key == "CXCL8_ER_neg") | (row_key == "PPARG_ER_pos")]
hm_dt[, txt_col  := ifelse(abs(rho) > 0.34, "white", "#0F172A")]

# Top lineage annotation blocks (exact x positions)
lin_blocks <- ord_helper[, .(xmin = min(match(feat_clean, feat_levels)) - 0.46,
                             xmax = max(match(feat_clean, feat_levels)) + 0.46,
                             xmid = mean(match(feat_clean, feat_levels))), by = lineage]
lin_pal <- c("Myeloid / Innate" = "#E64B35", "Lymphoid" = "#4DBBD5", "Stromal / Genomic" = "#00A087")
lin_blocks[, fill_col := unname(lin_pal[as.character(lineage)])]

# Vertical group separators
x_div1 <- lin_blocks[lineage == "Myeloid / Innate", xmax] + 0.04
x_div2 <- lin_blocks[lineage == "Lymphoid", xmax] + 0.04

# Subtle focal row frames
focal_rows <- data.table(
  ymin = c(3.52, 0.52), ymax = c(4.48, 1.48),
  xmin = 0.52, xmax = length(feat_levels) + 0.48,
  col  = c(pal_cns$cxcl8, pal_cns$pparg)
)

p4a <- ggplot(hm_dt) +
  # Clean crisp matrix tiles
  geom_tile(aes(x = x_num, y = y_num), fill = "white", color = "#E2E8F0", linewidth = 0.32, width = 1, height = 1) +
  geom_tile(data = hm_dt[is_focal == TRUE], aes(x = x_num, y = y_num), fill = "#F8FAFC", color = "#CBD5E1", linewidth = 0.32, width = 1, height = 1) +
  # Group vertical dividers
  annotate("segment", x = c(x_div1, x_div2), xend = c(x_div1, x_div2), y = 0.5, yend = 4.5,
           color = "#94A3B8", linewidth = 0.45) +
  # Horizontal dashed divider between CXCL8 and PPARG
  annotate("segment", x = 0.5, xend = length(feat_levels) + 0.5, y = 2.5, yend = 2.5,
           color = "#64748B", linewidth = 0.42, linetype = "dashed") +
  # Top lineage color ribbon + crisp category titles
  geom_rect(data = lin_blocks, aes(xmin = xmin, xmax = xmax, ymin = 4.66, ymax = 4.86),
            inherit.aes = FALSE, fill = lin_blocks$fill_col, color = NA) +
  geom_text(data = lin_blocks, aes(x = xmid, y = 5.06, label = lineage),
            inherit.aes = FALSE, size = 2.55, fontface = "bold", color = "#1E293B") +
  # Circles scaled by |rho| and colored by Spearman rho
  geom_point(aes(x = x_num, y = y_num, size = abs(rho), fill = rho),
             shape = 21, color = "white", stroke = 0.35) +
  # Significance stars inside circles
  geom_text(aes(x = x_num, y = y_num, label = sig, color = txt_col),
            size = 2.25, vjust = 0.72, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() +
  # Subtle highlight border around focal rows
  geom_rect(data = focal_rows, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = NA, color = focal_rows$col, linewidth = 0.58) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#C83E2B",
    midpoint = 0, limits = c(-0.65, 0.65), oob = scales::squish,
    breaks = c(-0.5, 0, 0.5),
    name = "Spearman \u03c1"
  ) +
  scale_size_continuous(range = c(1.6, 5.8), limits = c(0, 0.65), breaks = c(0.2, 0.4, 0.6), name = "| \u03c1 |") +
  scale_x_continuous(breaks = seq_along(feat_levels), labels = feat_levels,
                     limits = c(0.45, length(feat_levels) + 0.55), expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:4, labels = unname(row_labels[row_levels]),
                     limits = c(0.42, 5.28), expand = c(0, 0)) +
  labs(
    title = cns_panel_title("A", "Tumor microenvironment correlations across TCGA-BRCA strata"),
    x = NULL, y = NULL
  ) +
  theme_cns_fig4(8.2) +
  theme(
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    axis.text.y      = element_markdown(size = 7.4, color = pal_cns$ink, margin = margin(r = 4)),
    axis.text.x      = element_text(size = 6.9, color = pal_cns$subink, angle = 38, hjust = 1, vjust = 1, margin = margin(t = 2)),
    plot.title       = element_markdown(size = 8.6, hjust = 0, margin = margin(b = 4)),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(2.5, "mm"),
    legend.title     = element_text(face = "bold", size = 7.0, color = pal_cns$ink),
    legend.text      = element_text(size = 6.8, color = pal_cns$subink),
    plot.margin      = margin(4, 4, 2, 4)
  ) +
  guides(
    fill = guide_colorbar(barwidth = unit(2.8, "mm"), barheight = unit(16, "mm"),
                          frame.colour = "#CBD5E1", ticks.colour = pal_cns$subink, order = 1),
    size = guide_legend(override.aes = list(fill = "#94A3B8", color = "white"), order = 2)
  )

# ---------------------------------------------------------------------------
# 2. Panel B: MSigDB Hallmark Enrichment (Clean Side-by-Side Dual Panel)
# ---------------------------------------------------------------------------
gsea_map <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"           = "TNF\u03b1\u2013NF-\u03baB signaling",
  "HALLMARK_MTORC1_SIGNALING"                  = "mTORC1 signaling",
  "HALLMARK_G2M_CHECKPOINT"                    = "G2/M checkpoint",
  "HALLMARK_HYPOXIA"                           = "Hypoxia",
  "HALLMARK_GLYCOLYSIS"                        = "Glycolysis",
  "HALLMARK_INFLAMMATORY_RESPONSE"             = "Inflammatory response",
  "HALLMARK_E2F_TARGETS"                       = "E2F targets",
  "HALLMARK_ANGIOGENESIS"                      = "Angiogenesis",
  "HALLMARK_ALLOGRAFT_REJECTION"               = "Allograft rejection",
  "HALLMARK_KRAS_SIGNALING_UP"                 = "KRAS signaling up",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"           = "IL-6/JAK/STAT3 signaling",
  "HALLMARK_COAGULATION"                       = "Coagulation",
  "HALLMARK_MYOGENESIS"                        = "Myogenesis"
)

gsea5  <- fread(file.path(res_dir, "module5_hallmark_gsea_results.csv"))
gs_top <- gsea5[((target_gene == "CXCL8" & stratum == "ER_neg") | (target_gene == "PPARG" & stratum == "ER_pos")) & p.adjust < 0.05][
  order(target_gene, -abs(NES))][, .SD[1:min(8, .N)], by = target_gene]
gs_top[, clean_id := ifelse(ID %in% names(gsea_map), unname(gsea_map[ID]), gsub("_", " ", sub("^HALLMARK_", "", ID)))]
gs_top[, neglog10_fdr := pmin(-log10(p.adjust), 10)]

make_hallmark_panel <- function(dt_sub, accent_col, soft_col, sub_hdr, show_title = FALSE) {
  d <- copy(dt_sub)[order(NES)]
  d[, id_f := factor(clean_id, levels = clean_id)]
  ggplot(d, aes(x = NES, y = id_f)) +
    geom_vline(xintercept = 0, color = pal_cns$ink, linewidth = 0.36) +
    geom_segment(aes(x = 0, xend = NES, yend = id_f), color = soft_col, linewidth = 3.6, lineend = "round") +
    geom_segment(aes(x = 0, xend = NES, yend = id_f), color = accent_col, linewidth = 0.72) +
    geom_point(aes(size = neglog10_fdr), shape = 21, fill = accent_col, color = "white", stroke = 0.45) +
    geom_text(aes(x = NES + 0.14, label = sprintf("%.2f", NES)), hjust = 0, size = 2.35, color = pal_cns$ink) +
    scale_size_continuous(range = c(2.4, 4.8), limits = c(1.3, 10), breaks = c(2, 5, 8, 10),
                          labels = c("2", "5", "8", "\u226510"), name = expression(-log[10](FDR))) +
    scale_x_continuous(breaks = c(0, 1, 2, 3), expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(0, 3.45)) +
    labs(
      title    = if (show_title) cns_panel_title("B", "MSigDB Hallmark enrichment") else " ",
      subtitle = sub_hdr,
      x        = "Normalized enrichment score (NES)",
      y        = NULL
    ) +
    theme_cns_fig4(8.0) +
    theme(
      axis.line.y   = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.text.y   = element_text(size = 7.2, color = pal_cns$ink),
      plot.title    = element_markdown(size = 8.6, hjust = 0, margin = margin(b = 3)),
      plot.subtitle = element_markdown(size = 7.7, hjust = 0, margin = margin(b = 4)),
      plot.margin   = margin(4, 6, 2, 4)
    ) +
    guides(size = guide_legend(override.aes = list(fill = "#64748B", color = "white")))
}

p4b_left <- make_hallmark_panel(
  gs_top[target_gene == "CXCL8"],
  accent_col = pal_cns$cxcl8, soft_col = pal_cns$cxcl8_soft,
  sub_hdr    = "<span style='color:#E64B35;'>&#9644;</span> &nbsp;<b>*CXCL8*</b> | ER\u2212 (n = 237)",
  show_title = TRUE
)
p4b_right <- make_hallmark_panel(
  gs_top[target_gene == "PPARG"],
  accent_col = pal_cns$pparg, soft_col = pal_cns$pparg_soft,
  sub_hdr    = "<span style='color:#3C5488;'>&#9644;</span> &nbsp;<b>*PPARG*</b> | ER+ (n = 794)",
  show_title = FALSE
)

# Extract ONE single neutral shared legend for Panel B
leg_b <- cowplot::get_legend(
  p4b_left + theme(legend.position = "bottom", legend.margin = margin(0, 0, 2, 0),
                   legend.title = element_text(face = "bold", size = 7.2),
                   legend.text  = element_text(size = 7.0))
)
p4b_row <- cowplot::plot_grid(
  p4b_left  + theme(legend.position = "none"),
  p4b_right + theme(legend.position = "none"),
  ncol = 2, align = "h", rel_widths = c(1, 1)
)
p4b_full <- cowplot::plot_grid(p4b_row, leg_b, ncol = 1, rel_heights = c(1, 0.13))

# ---------------------------------------------------------------------------
# 3. Panel C: Crisp Multi-Curve GSEA Running ES (No muddy fill) + Barcode Rug
# ---------------------------------------------------------------------------
gsea_objs <- readRDS(file.path(bulk_dir, "planA_cns_gsea_objects.rds"))
curve_specs <- list(
  list(obj = gsea_objs$CXCL8_ERneg, id = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
       lbl = "*CXCL8* (ER\u2212): EMT (NES = 2.87)", col = "#E64B35"),
  list(obj = gsea_objs$CXCL8_ERneg, id = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
       lbl = "*CXCL8* (ER\u2212): TNF\u03b1\u2013NF-\u03baB (NES = 2.77)", col = "#E08E2B"),
  list(obj = gsea_objs$PPARG_ERpos, id = "HALLMARK_ALLOGRAFT_REJECTION",
       lbl = "*PPARG* (ER+): Allograft rejection (NES = 1.58)", col = "#3C5488"),
  list(obj = gsea_objs$PPARG_ERpos, id = "HALLMARK_ANGIOGENESIS",
       lbl = "*PPARG* (ER+): Angiogenesis (NES = 1.43)", col = "#00A087")
)

curve_dt_list <- vector("list", length(curve_specs))
tick_dt_list  <- vector("list", length(curve_specs))
for (k in seq_along(curve_specs)) {
  sp <- curve_specs[[k]]
  gi <- as.data.table(enrichplot:::gsInfo(sp$obj, geneSetID = sp$id))
  idx_keep <- unique(c(seq(1, nrow(gi), by = max(1L, floor(nrow(gi) / 480))), which.max(abs(gi$runningScore)), nrow(gi)))
  curve_dt_list[[k]] <- gi[idx_keep, .(x, runningScore, pathway = sp$lbl)]
  tick_dt_list[[k]]  <- gi[position == 1, .(x, pathway = sp$lbl, track = 5 - k)]
}
curve_dt <- rbindlist(curve_dt_list)
tick_dt  <- rbindlist(tick_dt_list)
pw_lvls  <- vapply(curve_specs, `[[`, character(1), "lbl")
pw_pal   <- setNames(vapply(curve_specs, `[[`, character(1), "col"), pw_lvls)
curve_dt[, pathway := factor(pathway, levels = pw_lvls)]
tick_dt[,  pathway := factor(pathway, levels = pw_lvls)]
max_rank <- max(curve_dt$x, na.rm = TRUE)

p4c_curve <- ggplot(curve_dt, aes(x = x, y = runningScore, color = pathway)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#CBD5E1", linewidth = 0.40) +
  geom_line(linewidth = 0.95) +
  scale_color_manual(values = pw_pal, name = NULL) +
  scale_x_continuous(limits = c(0, max_rank), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75), expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    title = cns_panel_title("C", "Representative GSEA running enrichment curves"),
    x     = NULL,
    y     = "Running ES"
  ) +
  theme_cns_fig4(8.0) +
  theme(
    axis.line.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.x       = element_blank(),
    plot.title        = element_markdown(size = 8.6, hjust = 0, margin = margin(b = 4)),
    legend.position   = c(0.76, 0.72),
    legend.key.width  = unit(6.0, "mm"),
    legend.key.height = unit(4.0, "mm"),
    legend.text       = element_markdown(size = 7.2, color = pal_cns$ink),
    plot.margin       = margin(4, 8, 0, 6)
  )

p4c_rug <- ggplot(tick_dt, aes(x = x, xend = x, y = track - 0.36, yend = track + 0.36, color = pathway)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.4, ymax = 4.6, fill = "#F8FAFC", color = NA) +
  geom_segment(linewidth = 0.28, alpha = 0.88, show.legend = FALSE) +
  scale_color_manual(values = pw_pal) +
  scale_x_continuous(limits = c(0, max_rank), breaks = c(0, 5000, 10000, 15000), labels = scales::comma, expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0.4, 4.6), expand = c(0, 0)) +
  labs(x = expression(paste("Gene rank (limma ", log[2], "FC, high to low)")), y = NULL) +
  theme_cns_fig4(8.0) +
  theme(
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    plot.margin  = margin(0, 8, 4, 6)
  )

p4c_full <- cowplot::plot_grid(p4c_curve, p4c_rug, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.76, 0.24))

# ---------------------------------------------------------------------------
# 4. Assemble Full CNS Figure 4 with Generous Breathing Room
# ---------------------------------------------------------------------------
fig4_cns <- cowplot::plot_grid(
  p4a,
  p4b_full,
  p4c_full,
  ncol = 1,
  rel_heights = c(0.36, 0.34, 0.30)
)

save_cns_fig4(fig4_cns, "Fig4_TME_Correlations_and_Hallmark_GSEA_CNS", width_mm = 185, height_mm = 210, dpi_png = 350)
cat("Redesigned CNS Fig 4 complete.\n")
