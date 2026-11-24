# 02_fig7_upstream_mr.R
# Figure 7 CNS-Grade Redesign (100-pt Nature Medicine / Nature Microbiology standard)
#   - Unified color semantics across all sub-panels:
#       Steel blue (#3C5488) = protective (OR < 1); Coral red (#E64B35) = risk (OR > 1)
#       Sage green halo (#00A087) + star = FDR < 0.05
#   - Subtypes encoded cleanly by shape & fixed vertical track (Circle=Overall, Square=ER+, Diamond=ER-)
#   - Complete continuous colorbar (OR) + bubble size legend (-log10 P) + shape/FDR key
#   - Panels C & D: Split Bipartite Bezier Evidence Networks (exact gutMGene v2.0 PMIDs, no cross-table chaining)
# Strict data & statistical parity with 05_revision_v2/01_scripts/02_fig7_upstream_mr.R

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "06_figures_cns_polish/01_scripts/00_cns_theme.R"))

cns_panel_title <- function(letter, subtitle_md) {
  sprintf(
    "<span style='font-size:11.5pt; color:#111827'><b>%s</b></span>\u2003<span style='font-size:8.6pt; color:#111827'>%s</span>",
    letter, subtitle_md
  )
}

cat("=== Generating 100-pt Multi-Modal CNS Fig 7 ===\n")

read_res <- function(name) fread(file.path(res_dir, paste0(name, ".csv")))

nominal_mib <- read_res("mibiogen_bcac_mr_significant_p005")[gwas_name != "Marvinbryantia"]
nominal_met <- read_res("metabolite_bcac_mr_significant_p005")[
  id.exposure %chin% c("met-d-Acetate", "met-a-462", "met-a-385", "met-c-893", "met-a-446")]

nominal_met[, short_name := fcase(
  id.exposure == "met-a-385", "3-Phenylpropionate",
  id.exposure == "met-c-893", "Linoleic acid (18:2)",
  id.exposure == "met-a-462", "3-(4-hydroxyphenyl)<br>lactate",
  default = trait_name
)]

format_taxon_md <- function(x) {
  fifelse(grepl("aceae$", x), x, sprintf("<i>%s</i>", x))
}

sub_levels <- c("Overall", "ER+", "ER\u2212")
sub_short  <- c("Overall" = "Ov", "ER+" = "ER+", "ER\u2212" = "ER\u2212")
sub_shapes <- c("Overall" = 21L, "ER+" = 22L, "ER\u2212" = 23L)
sub_offset <- c("Overall" = 0.25, "ER+" = 0.00, "ER\u2212" = -0.25)

# ------------------------------------------------------------------------------
# Composite Builder: Bubble Matrix (Left) | Track Forest (Mid) | Exact Stats (Right)
# ------------------------------------------------------------------------------
build_dual_mr_panel <- function(dt_in, name_col, exp_order, panel_letter, panel_sub,
                                is_taxon = FALSE, xlim_vals = c(0.50, 1.98)) {
  d <- copy(dt_in)
  d[, exp_name := get(name_col)]
  d[, outcome  := factor(clean_subtype_cns(outcome_label), levels = sub_levels)]
  d[, log_or   := log(OR)]
  d[, neglogp  := -log10(pval)]
  d[, is_fdr05 := FDR < 0.05]
  d[, is_fdr10 := FDR >= 0.05 & FDR < 0.10]
  d[, is_2snp  := nsnp == 2L]
  d[, sig_lbl  := fcase(is_fdr05, "\u2605", is_fdr10, "*", default = "")]

  # Unified effect-direction color (Blue = OR < 1, Coral red = OR > 1)
  d[, dir_col  := fifelse(OR >= 1, pal_cns$cxcl8, pal_cns$pparg)]
  d[, pt_fill  := fifelse(is_fdr05 | is_fdr10, dir_col, "white")]
  d[, pt_shape := unname(sub_shapes[as.character(outcome)])]
  d[, oc_tag   := unname(sub_short[as.character(outcome)])]

  exp_meta <- unique(d[, .(exp_name, nsnp)])
  exp_meta <- exp_meta[match(exp_order, exp_name)]
  exp_meta[, row := rev(seq_len(.N))]
  raw_nm <- if (is_taxon) format_taxon_md(exp_meta$exp_name) else exp_meta$exp_name
  snp_tag <- fifelse(exp_meta$nsnp == 2L, "2\u2020", as.character(exp_meta$nsnp))
  exp_meta[, y_lbl_md := sprintf("%s <span style='color:#9CA3AF; font-size:6.3pt'>(%s)</span>", raw_nm, snp_tag)]

  d <- merge(d, exp_meta[, .(exp_name, row)], by = "exp_name")
  d[, y_track := row + unname(sub_offset[as.character(outcome)])]

  grid_all <- CJ(exp_name = exp_order, outcome = factor(sub_levels, levels = sub_levels))
  grid_all <- merge(grid_all, exp_meta[, .(exp_name, row)], by = "exp_name")
  grid_all <- merge(grid_all, d[, .(exp_name, outcome, OR, log_or, neglogp, pval, FDR, is_fdr05, is_fdr10, sig_lbl)],
                    by = c("exp_name", "outcome"), all.x = TRUE)
  grid_all[, has_sig := !is.na(OR)]

  n_exp <- nrow(exp_meta)
  zebra <- exp_meta[row %% 2L == 1L]
  y_lim <- c(0.45, n_exp + 0.55)

  # --- Left: Cross-subtype Effect Bubble Matrix ---
  p_bub <- ggplot(grid_all, aes(x = outcome, y = row)) +
    geom_rect(data = zebra, aes(ymin = row - 0.48, ymax = row + 0.48),
              xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = pal_cns$zebra) +
    geom_hline(yintercept = n_exp + 0.52, color = pal_cns$ink, linewidth = 0.45) +
    geom_point(data = grid_all[has_sig == FALSE], shape = 16, size = 1.05, color = "#E2E8F0") +
    geom_point(data = grid_all[has_sig == TRUE],
               aes(size = neglogp, fill = log_or),
               shape = 21, color = pal_cns$ink, stroke = 0.35) +
    geom_point(data = grid_all[is_fdr05 == TRUE],
               aes(size = neglogp), shape = 21, fill = NA,
               color = pal_cns$teal, stroke = 1.08) +
    geom_text(data = grid_all[nzchar(sig_lbl)], aes(label = sig_lbl),
              color = "white", size = 2.3, fontface = "bold", vjust = 0.52) +
    scale_fill_gradient2(
      low = pal_cns$pparg, mid = "#F8FAFC", high = pal_cns$cxcl8,
      midpoint = 0, limits = c(-0.38, 0.38), oob = scales::squish, guide = "none"
    ) +
    scale_size_continuous(range = c(3.2, 6.2), limits = c(1.2, 3.6), guide = "none") +
    scale_x_discrete(position = "top", expand = expansion(add = 0.55)) +
    scale_y_continuous(breaks = exp_meta$row, labels = exp_meta$y_lbl_md, limits = y_lim, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_cns(7.8) +
    theme(
      axis.line   = element_blank(),
      axis.ticks  = element_blank(),
      axis.text.x = ggtext::element_markdown(face = "bold", size = 7.2, color = pal_cns$subink, margin = margin(b = 3)),
      axis.text.y = ggtext::element_markdown(size = 7.3, color = pal_cns$ink, hjust = 1, margin = margin(r = 4)),
      plot.margin = margin(1, 2, 1, 2)
    )

  # --- Middle: Multi-Subtype Parallel Track Forest (Color = Direction, Shape = Subtype) ---
  p_for <- ggplot(d, aes(x = OR, y = y_track)) +
    annotate("rect", xmin = xlim_vals[1], xmax = 1, ymin = y_lim[1], ymax = n_exp + 0.52,
             fill = "#F4F7FB", alpha = 0.65) +
    annotate("rect", xmin = 1, xmax = xlim_vals[2], ymin = y_lim[1], ymax = n_exp + 0.52,
             fill = "#FFF5F4", alpha = 0.65) +
    geom_hline(yintercept = seq(1.5, n_exp - 0.5, by = 1), color = "#E2E8F0", linewidth = 0.30) +
    geom_hline(yintercept = n_exp + 0.52, color = pal_cns$ink, linewidth = 0.45) +
    geom_vline(xintercept = 1, linetype = "dashed", color = pal_cns$muted, linewidth = 0.38) +
    geom_segment(aes(x = OR_LCI95, xend = OR_UCI95, y = y_track, yend = y_track, color = dir_col), linewidth = 0.68) +
    geom_point(aes(color = dir_col, fill = pt_fill, shape = pt_shape), size = 2.05, stroke = 0.68) +
    geom_point(data = d[is_fdr05 == TRUE], aes(shape = pt_shape), size = 3.25, fill = NA,
               color = pal_cns$teal, stroke = 0.72) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_shape_identity() +
    scale_x_log10(breaks = c(0.6, 0.8, 1.0, 1.4, 1.8), labels = c("0.6", "0.8", "1.0", "1.4", "1.8")) +
    coord_cartesian(xlim = xlim_vals, ylim = y_lim, expand = FALSE) +
    labs(title = "<b>Odds ratio (95% CI, log scale)</b>", x = NULL, y = NULL) +
    theme_cns(7.8) +
    theme(
      axis.line.y  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x  = ggtext::element_markdown(size = 6.8, color = pal_cns$subink),
      plot.title   = ggtext::element_markdown(size = 7.2, color = pal_cns$subink, hjust = 0.5, margin = margin(b = 2, l = 0)),
      plot.margin  = margin(1, 3, 1, 3)
    )

  # --- Right: Track-aligned exact OR (95% CI) and P values (direction-colored + subtype tag) ---
  d[, stat_md := sprintf(
    "<span style='color:#6B7280; font-size:5.7pt'>%s</span>\u2002<span style='color:%s'>%s</span>",
    oc_tag, dir_col,
    fifelse(is_fdr05,
            sprintf("<b>%.2f (%.2f\u2013%.2f) <i>P</i>=%.2g</b>", OR, OR_LCI95, OR_UCI95, pval),
            sprintf("%.2f (%.2f\u2013%.2f) <i>P</i>=%.2g", OR, OR_LCI95, OR_UCI95, pval))
  )]

  p_txt <- ggplot(d, aes(y = y_track)) +
    geom_rect(data = zebra, aes(ymin = row - 0.48, ymax = row + 0.48),
              xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = pal_cns$zebra) +
    geom_hline(yintercept = seq(1.5, n_exp - 0.5, by = 1), color = "#E2E8F0", linewidth = 0.30) +
    geom_hline(yintercept = n_exp + 0.52, color = pal_cns$ink, linewidth = 0.45) +
    ggtext::geom_richtext(
      aes(x = 0.01, label = stat_md),
      hjust = 0, vjust = 0.5, size = 1.98,
      fill = NA, label.color = NA, label.padding = unit(0, "pt")
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    labs(title = "<b>OR (95% CI) &amp; <i>P</i></b>") +
    theme_void(base_family = "sans") +
    theme(
      plot.title  = ggtext::element_markdown(size = 7.2, color = pal_cns$subink, hjust = 0, margin = margin(b = 2, l = 2)),
      plot.margin = margin(1, 2, 1, 1)
    )

  combo <- (p_bub | p_for | p_txt) + plot_layout(widths = c(1.62, 1.78, 1.75))
  wrap_elements(full = combo +
    plot_annotation(
      title = cns_panel_title(panel_letter, panel_sub),
      theme = theme(
        plot.title  = ggtext::element_markdown(margin = margin(b = 2, l = 0)),
        plot.margin = margin(2, 2, 1, 2)
      )
    ))
}

mib_order <- c("Adlercreutzia", "Parabacteroides", "Desulfovibrio", "Dorea", "Ruminococcaceae", "Blautia")
met_order <- c("Acetate", "3-(4-hydroxyphenyl)<br>lactate", "3-Phenylpropionate", "Indoleacetate", "Linoleic acid (18:2)")

fig7_a <- build_dual_mr_panel(nominal_mib, "gwas_name",  mib_order, "A",
                              "Gut microbiota MR across BCAC subtypes (MiBioGen; SNPs in parentheses, \u20202 SNPs)", is_taxon = TRUE)
fig7_b <- build_dual_mr_panel(nominal_met, "short_name", met_order, "B",
                              "Plasma metabolite MR across BCAC subtypes (OpenGWAS; SNPs in parentheses)", is_taxon = FALSE)

# ------------------------------------------------------------------------------
# Complete Continuous + Discrete Legend Strip (Zero Reviewer Gaps)
# ------------------------------------------------------------------------------
# 1. Continuous OR gradient bar
grad_dt <- data.table(x = seq(0.08, 0.28, length.out = 60),
                      val = seq(-0.38, 0.38, length.out = 60))
# 2. Bubble size reference points (-log10 P = 1.5, 2.5, 3.5)
sz_dt <- data.table(x = c(0.415, 0.455, 0.500),
                    val = c(1.5, 2.5, 3.5),
                    lbl = c("1.5", "2.5", "3.5"))
# 3. Shape & FDR key on the right
key_right <- paste0(
  "<span style='color:#374151'><b>\u25cf</b> Overall\u2002<b>\u25a0</b> ER+\u2002<b>\u25c6</b> ER\u2212</span>\u2003|\u2003",
  "<span style='color:#00A087'><b>\u25ce \u2605</b> FDR&lt;0.05</span>\u2002",
  "<span style='color:#374151'><b>\u25cf *</b> FDR&lt;0.10\u2002\u25cb <i>P</i>&lt;0.05</span>"
)

p_leg_bar <- ggplot() +
  # OR colorbar
  annotate("text", x = 0.068, y = 0.62, label = "OR", hjust = 1, size = 2.25, fontface = "bold", color = pal_cns$ink) +
  geom_tile(data = grad_dt, aes(x = x, y = 0.62, fill = val), width = 0.0045, height = 0.28) +
  annotate("rect", xmin = 0.078, xmax = 0.282, ymin = 0.48, ymax = 0.76, fill = NA, color = pal_cns$muted, linewidth = 0.28) +
  annotate("text", x = 0.08, y = 0.26, label = "0.7", size = 1.95, color = pal_cns$pparg, fontface = "bold") +
  annotate("text", x = 0.18, y = 0.26, label = "1.0", size = 1.95, color = pal_cns$subink) +
  annotate("text", x = 0.28, y = 0.26, label = "1.4", size = 1.95, color = pal_cns$cxcl8, fontface = "bold") +
  scale_fill_gradient2(low = pal_cns$pparg, mid = "#F8FAFC", high = pal_cns$cxcl8, midpoint = 0, guide = "none") +
  # -log10 P size legend
  ggtext::geom_richtext(aes(x = 0.39, y = 0.62, label = "<b>\u2212log<sub>10</sub> <i>P</i></b>"),
                        hjust = 1, size = 2.15, color = pal_cns$ink, fill = NA, label.color = NA) +
  geom_point(data = sz_dt, aes(x = x, y = 0.62, size = val), shape = 21, fill = "#CBD5E1", color = pal_cns$ink, stroke = 0.32) +
  geom_text(data = sz_dt, aes(x = x, y = 0.22, label = lbl), size = 1.95, color = pal_cns$subink) +
  scale_size_continuous(range = c(2.6, 5.2), limits = c(1.2, 3.6), guide = "none") +
  # Divider + Shape/FDR key
  annotate("segment", x = 0.545, xend = 0.545, y = 0.25, yend = 0.85, color = "#CBD5E1", linewidth = 0.35) +
  ggtext::geom_richtext(aes(x = 0.565, y = 0.55, label = key_right),
                        hjust = 0, vjust = 0.5, size = 2.12, fill = NA, label.color = NA) +
  coord_cartesian(xlim = c(0.02, 0.99), ylim = c(0.05, 0.95), expand = FALSE) +
  theme_void(base_family = "sans") +
  theme(plot.margin = margin(1, 0, 2, 0))

fig7_key <- wrap_elements(full = p_leg_bar)

# ------------------------------------------------------------------------------
# Panels C & D: Split Bipartite Bezier Evidence Networks (gutMGene v2.0 PMIDs)
# ------------------------------------------------------------------------------
gut_dir <- dirname(cache_dir)
microbe <- fread(file.path(gut_dir, "gutMGene_v2.0_Human_Gut_Microbe_Microbial_Metabolite_official_browse.csv"))
gene_tb <- fread(file.path(gut_dir, "gutMGene_v2.0_Human_Microbial_Metabolite_Host_Gene_official_browse.csv"))

bezier_curve_df <- function(x0, y0, x1, y1, edge_id, col, n = 50) {
  t <- seq(0, 1, length.out = n)
  cx1 <- x0 + 0.48 * (x1 - x0); cy1 <- y0
  cx2 <- x0 + 0.52 * (x1 - x0); cy2 <- y1
  bx <- (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * cx1 + 3 * (1 - t) * t^2 * cx2 + t^3 * x1
  by <- (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * cy1 + 3 * (1 - t) * t^2 * cy2 + t^3 * y1
  data.table(x = bx, y = by, edge_id = edge_id, col = col)
}

bezier_pt <- function(x0, y0, x1, y1, t) {
  cx1 <- x0 + 0.48 * (x1 - x0); cy1 <- y0
  cx2 <- x0 + 0.52 * (x1 - x0); cy2 <- y1
  bx <- (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * cx1 + 3 * (1 - t) * t^2 * cx2 + t^3 * x1
  by <- (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * cy1 + 3 * (1 - t) * t^2 * cy2 + t^3 * y1
  c(bx, by)
}

# --- Panel C: Taxon -> Metabolite Bipartite Network ---
ev_c <- microbe[gut_microbiota %chin% c("Dorea", "Blautia") &
                  metabolite %chin% c("Acetate", "Butyrate", "Succinate") & !is.na(pmid),
                .(pmid = paste(sort(unique(pmid)), collapse = "/")),
                by = .(left = gut_microbiota, right = metabolite)]

nodes_c_left <- data.table(
  name = c("Dorea", "Blautia"),
  lbl  = c("<i>Dorea</i>", "<i>Blautia</i>"),
  x = 0.0, y = c(2.65, 1.35),
  fill = c(pal_cns$pparg_soft, pal_cns$cxcl8_soft),
  edge = c(pal_cns$pparg, pal_cns$cxcl8),
  txt  = c(pal_cns$pparg_dark, pal_cns$cxcl8_dark),
  lty  = "solid"
)
nodes_c_right <- data.table(
  name = c("Acetate", "Butyrate", "Succinate"),
  lbl  = c("Acetate", "Butyrate", "Succinate"),
  x = 1.0, y = c(3.2, 2.0, 0.8),
  fill = c(pal_cns$teal_soft, pal_cns$teal_soft, "#F8FAFC"),
  edge = c(pal_cns$teal, pal_cns$teal, "#CBD5E1"),
  txt  = c(pal_cns$teal_dark, pal_cns$teal_dark, "#9CA3AF"),
  lty  = c("solid", "solid", "dashed")
)

ev_c <- merge(ev_c, nodes_c_left[, .(left = name, x0 = x, y0 = y)], by = "left")
ev_c <- merge(ev_c, nodes_c_right[, .(right = name, x1 = x, y1 = y)], by = "right")
ev_c[, edge_id := paste(left, right, sep = "_")]
ev_c[, t_pos := fcase(
  left == "Dorea"   & right == "Acetate",  0.30,
  left == "Dorea"   & right == "Butyrate", 0.72,
  left == "Blautia" & right == "Acetate",  0.72,
  left == "Blautia" & right == "Butyrate", 0.42
)]
ev_c[, ecol := fifelse(left == "Dorea", pal_cns$pparg, pal_cns$cxcl8)]

curves_c <- rbindlist(lapply(seq_len(nrow(ev_c)), function(i) {
  bezier_curve_df(ev_c$x0[i] + 0.17, ev_c$y0[i], ev_c$x1[i] - 0.17, ev_c$y1[i], ev_c$edge_id[i], ev_c$ecol[i])
}))
lbl_c <- rbindlist(lapply(seq_len(nrow(ev_c)), function(i) {
  pt <- bezier_pt(ev_c$x0[i] + 0.17, ev_c$y0[i], ev_c$x1[i] - 0.17, ev_c$y1[i], ev_c$t_pos[i])
  data.table(x = pt[1], y = pt[2], pmid = ev_c$pmid[i], col = ev_c$ecol[i])
}))
nodes_c <- rbind(nodes_c_left, nodes_c_right)

p7c <- ggplot() +
  geom_path(data = curves_c, aes(x = x, y = y, group = edge_id, color = col),
            linewidth = 1.05, alpha = 0.52, lineend = "round") +
  geom_tile(data = nodes_c, aes(x = x, y = y, fill = fill, color = edge, linetype = lty),
            width = 0.34, height = 0.58, linewidth = 0.55) +
  ggtext::geom_richtext(data = nodes_c, aes(x = x, y = y, label = lbl, color = txt),
                        size = 2.45, fontface = "bold", fill = NA, label.color = NA) +
  geom_label(data = lbl_c, aes(x = x, y = y, label = pmid, color = col),
             fill = "white", size = 1.95, fontface = "bold",
             label.padding = unit(0.10, "lines"), label.r = unit(0.14, "lines"), linewidth = 0.25) +
  annotate("text", x = 0.0, y = 3.82, label = "Gut taxon", size = 2.35, fontface = "bold", color = pal_cns$slate) +
  annotate("text", x = 1.0, y = 3.82, label = "Metabolite", size = 2.35, fontface = "bold", color = pal_cns$slate) +
  scale_color_identity() + scale_fill_identity() + scale_linetype_identity() +
  coord_cartesian(xlim = c(-0.22, 1.22), ylim = c(0.28, 4.02), expand = FALSE) +
  labs(title = cns_panel_title("C", "gutMGene taxon\u2013metabolite PMIDs")) +
  theme_void(base_family = "sans") +
  theme(
    plot.title  = ggtext::element_markdown(margin = margin(b = 3, l = 0)),
    plot.margin = margin(3, 6, 2, 2)
  )

# --- Panel D: Metabolite -> Host Gene Bipartite Network ---
ev_d <- gene_tb[gene %chin% c("CXCL8", "IL1B", "FFAR2", "PPARG") &
                  metabolite %chin% c("Acetate", "Butyrate", "Succinate") & !is.na(pmid),
                .(pmid = paste(sort(unique(pmid)), collapse = "/")),
                by = .(left = metabolite, right = gene)]

nodes_d_left <- data.table(
  name = c("Acetate", "Butyrate", "Succinate"),
  lbl  = c("Acetate", "Butyrate", "Succinate"),
  x = 0.0, y = c(3.2, 2.0, 0.8),
  fill = pal_cns$teal_soft, edge = pal_cns$teal, txt = pal_cns$teal_dark, lty = "solid"
)
nodes_d_right <- data.table(
  name = c("CXCL8", "IL1B", "FFAR2", "PPARG"),
  lbl  = c("<i>CXCL8</i>", "<i>IL1B</i>", "<i>FFAR2</i>", "<i>PPARG</i>"),
  x = 1.0, y = c(3.35, 2.45, 1.55, 0.65),
  fill = c(pal_cns$cxcl8_soft, "#F8FAFC", "#F8FAFC", pal_cns$pparg_soft),
  edge = c(pal_cns$cxcl8, "#94A3B8", "#94A3B8", pal_cns$pparg),
  txt  = c(pal_cns$cxcl8_dark, pal_cns$subink, pal_cns$subink, pal_cns$pparg_dark),
  lty  = "solid"
)

ev_d <- merge(ev_d, nodes_d_left[, .(left = name, x0 = x, y0 = y)], by = "left")
ev_d <- merge(ev_d, nodes_d_right[, .(right = name, x1 = x, y1 = y)], by = "right")
ev_d[, edge_id := paste(left, right, sep = "_")]
ev_d[, ecol := fcase(
  right == "CXCL8", pal_cns$cxcl8,
  right == "PPARG", pal_cns$pparg,
  default = pal_cns$teal
)]
ev_d[, t_pos := fcase(
  left == "Acetate"   & right == "CXCL8", 0.36,
  left == "Butyrate"  & right == "CXCL8", 0.68,
  left == "Succinate" & right == "CXCL8", 0.22,
  left == "Butyrate"  & right == "IL1B",  0.30,
  left == "Succinate" & right == "IL1B",  0.76,
  left == "Acetate"   & right == "FFAR2", 0.24,
  left == "Butyrate"  & right == "FFAR2", 0.70,
  left == "Butyrate"  & right == "PPARG", 0.56,
  default = 0.50
)]

curves_d <- rbindlist(lapply(seq_len(nrow(ev_d)), function(i) {
  bezier_curve_df(ev_d$x0[i] + 0.17, ev_d$y0[i], ev_d$x1[i] - 0.17, ev_d$y1[i], ev_d$edge_id[i], ev_d$ecol[i])
}))
lbl_d <- rbindlist(lapply(seq_len(nrow(ev_d)), function(i) {
  pt <- bezier_pt(ev_d$x0[i] + 0.17, ev_d$y0[i], ev_d$x1[i] - 0.17, ev_d$y1[i], ev_d$t_pos[i])
  data.table(x = pt[1], y = pt[2], pmid = ev_d$pmid[i], col = ev_d$ecol[i])
}))
nodes_d <- rbind(nodes_d_left, nodes_d_right)

p7d <- ggplot() +
  geom_path(data = curves_d, aes(x = x, y = y, group = edge_id, color = col),
            linewidth = 1.0, alpha = 0.50, lineend = "round") +
  geom_tile(data = nodes_d, aes(x = x, y = y, fill = fill, color = edge, linetype = lty),
            width = 0.34, height = 0.55, linewidth = 0.55) +
  ggtext::geom_richtext(data = nodes_d, aes(x = x, y = y, label = lbl, color = txt),
                        size = 2.45, fontface = "bold", fill = NA, label.color = NA) +
  geom_label(data = lbl_d, aes(x = x, y = y, label = pmid, color = col),
             fill = "white", size = 1.85, fontface = "bold",
             label.padding = unit(0.09, "lines"), label.r = unit(0.14, "lines"), linewidth = 0.25) +
  annotate("text", x = 0.0, y = 3.82, label = "Metabolite", size = 2.35, fontface = "bold", color = pal_cns$slate) +
  annotate("text", x = 1.0, y = 3.82, label = "Host gene", size = 2.35, fontface = "bold", color = pal_cns$slate) +
  scale_color_identity() + scale_fill_identity() + scale_linetype_identity() +
  coord_cartesian(xlim = c(-0.22, 1.22), ylim = c(0.28, 4.02), expand = FALSE) +
  labs(title = cns_panel_title("D", "gutMGene metabolite\u2013gene PMIDs")) +
  theme_void(base_family = "sans") +
  theme(
    plot.title  = ggtext::element_markdown(margin = margin(b = 3, l = 0)),
    plot.margin = margin(3, 2, 2, 6)
  )

fig7_cd <- wrap_elements(full = (p7c | p7d) + plot_layout(widths = c(1, 1.08)))

fig7_cns <- wrap_plots(
  fig7_a,
  fig7_b,
  fig7_key,
  fig7_cd,
  ncol = 1,
  heights = c(72, 62, 14, 62)
)

save_cns_figure(fig7_cns, "Fig7_Upstream_MR_and_gutMGene_Records", width_mm = 180, height_mm = 212, dpi_png = 350)
cat("100-pt Multi-Modal CNS Fig 7 complete.\n")
