# 02_fig3_ispy2.R
# Figure 3 CNS-Grade Polish v3 (Nature Medicine / Cancer Cell aesthetic)
# Strict data & statistical parity with 05_revision_v2/01_scripts/02_fig3_ispy2.R

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "06_figures_cns_polish/01_scripts/00_cns_theme.R"))

cat("=== Generating CNS-Polished Fig 3 (v3 High-Impact) ===\n")

th_dt  <- fread(file.path(res_dir, "planA_therapy_response_summary.csv"))
th_int <- fread(file.path(res_dir, "planA_therapy_interactions_and_multivar.csv"))
th_pc  <- readRDS(file.path(bulk_dir, "planA_therapy_plot_cache.rds"))

g173 <- copy(as.data.table(th_pc$gse173839))[pCR %in% c(0L, 1L)]
stopifnot(nrow(g173) == 100L, sum(g173$pCR) == 36L)

# ------------------------------------------------------------------------------
# Panel A: NEJM / Nature Medicine style aligned table-forest plot
#          with subtle directional background tint & significance-encoded markers
# ------------------------------------------------------------------------------
int_all <- copy(th_int[cohort == "GSE173839_ISPY2" & est_type == "Ratio_OR"])
stopifnot(nrow(int_all) == 8L)

int_all[, gene := fifelse(grepl("CXCL8", test_name), "CXCL8", "PPARG")]
int_all[, is_primary := test_name == "pCR_Interaction_CXCL8_x_DurvaArm(Adj_HR)"]
int_all[, is_sig := fdr_bh_all_ispy2_interactions < 0.05]

wrap_bold <- function(txt, is_b) fifelse(is_b, sprintf("<b>%s</b>", txt), txt)

int_all[, model_raw := fcase(
  test_name == "pCR_Interaction_CXCL8_x_DurvaArm(Adj_HR)",              "<i>CXCL8</i> \u00d7 Arm (adj. HR)",
  test_name == "pCR_Interaction_PPARG_x_DurvaArm(Adj_HR)",              "<i>PPARG</i> \u00d7 Arm (adj. HR)",
  test_name == "pCR_Interaction_CXCL8_x_DurvaArm(Adj_HR_MP)",           "<i>CXCL8</i> \u00d7 Arm (adj. HR + MP)",
  test_name == "pCR_Interaction_PPARG_x_DurvaArm(Adj_HR_MP)",           "<i>PPARG</i> \u00d7 Arm (adj. HR + MP)",
  test_name == "pCR_Interaction_CXCL8_x_DurvaArm(HRneg_Subset)",        "<i>CXCL8</i> \u00d7 Arm (HR\u2212, unadj.)",
  test_name == "pCR_Interaction_CXCL8_x_DurvaArm(HRneg_Subset_Adj_MP)", "<i>CXCL8</i> \u00d7 Arm (HR\u2212 + MP)",
  test_name == "pCR_Interaction_PPARG_x_DurvaArm(HRneg_Subset)",        "<i>PPARG</i> \u00d7 Arm (HR\u2212, unadj.)",
  test_name == "pCR_Interaction_PPARG_x_DurvaArm(HRneg_Subset_Adj_MP)", "<i>PPARG</i> \u00d7 Arm (HR\u2212 + MP)"
)]
int_all[, model_md := wrap_bold(model_raw, is_primary)]
int_all[, arms_md  := wrap_bold(sprintf("%d/%d vs %d/%d", events_control, n_control, events_durva, n_durva), is_primary)]
int_all[, or_raw   := fifelse(
  estimate < 0.05,
  sprintf("%.4f (%.4f\u2013%.3f)", estimate, lci95, uci95),
  sprintf("%.3f (%.3f\u2013%.3f)", estimate, lci95, uci95)
)]
int_all[, or_md := wrap_bold(or_raw, is_primary)]
int_all[, p_md  := wrap_bold(sprintf("%.3g", pval), is_primary)]
int_all[, q_md  := wrap_bold(sprintf("%.3g", fdr_bh_all_ispy2_interactions), is_primary)]

build_int_block_cns <- function(dt_sub, order_tests, xlim_vals, x_breaks, x_labels,
                                show_headers = TRUE, x_lab = NULL, block_title = NULL,
                                highlight_primary = FALSE) {
  d <- copy(dt_sub)[match(order_tests, test_name)]
  d[, row := rev(seq_len(.N))]
  d[, fill_col := fifelse(is_sig, fifelse(gene == "CXCL8", pal_cns$cxcl8, pal_cns$pparg), "white")]
  d[, edge_col := fifelse(gene == "CXCL8", pal_cns$cxcl8, pal_cns$pparg)]
  d[, pt_shape := fifelse(grepl("HRneg", test_name), 23L, 21L)]
  zebra <- d[row %% 2L == 1L]
  y_sc <- scale_y_continuous(limits = c(0.42, nrow(d) + 0.58), expand = c(0, 0))

  bg_layer <- if (highlight_primary) {
    list(
      geom_rect(data = zebra, aes(ymin = row - 0.48, ymax = row + 0.48),
                xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = pal_cns$zebra),
      geom_rect(data = d[is_primary == TRUE],
                aes(ymin = row - 0.48, ymax = row + 0.48),
                xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = pal_cns$cxcl8_soft)
    )
  } else {
    list(geom_rect(data = zebra, aes(ymin = row - 0.48, ymax = row + 0.48),
                   xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = pal_cns$zebra))
  }

  top_rule <- if (show_headers) {
    geom_hline(yintercept = nrow(d) + 0.55, color = pal_cns$ink, linewidth = 0.45)
  } else {
    geom_hline(yintercept = nrow(d) + 0.55, color = pal_cns$ref_line, linewidth = 0.32)
  }

  col_text_md <- function(col, header_md, align = "left") {
    x_pos <- if (align == "right") 0.96 else 0.02
    hj    <- if (align == "right") 1 else 0
    ggplot(d, aes(y = row)) + bg_layer + top_rule +
      ggtext::geom_richtext(
        aes(x = x_pos, label = .data[[col]]),
        hjust = hj, vjust = 0.5, size = 2.25,
        fill = NA, label.color = NA, label.padding = unit(0, "pt"),
        color = pal_cns$ink
      ) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + y_sc +
      labs(title = if (show_headers) header_md else NULL) +
      theme_void(base_family = "sans") +
      theme(
        plot.title = ggtext::element_markdown(size = 7.3, color = pal_cns$subink, hjust = hj, margin = margin(b = 2)),
        plot.margin = margin(1, 1.5, 1, 1.5)
      )
  }

  p_mid <- ggplot(d, aes(x = estimate, y = row)) +
    # Subtle directional background zones (Ratio OR < 1 vs > 1)
    annotate("rect", xmin = 0, xmax = 1, ymin = 0.42, ymax = nrow(d) + 0.55,
             fill = "#F4F7FB", alpha = 0.65) +
    annotate("rect", xmin = 1, xmax = Inf, ymin = 0.42, ymax = nrow(d) + 0.55,
             fill = "#FFF5F4", alpha = 0.65) +
    top_rule +
    geom_vline(xintercept = 1, linetype = "dashed", color = pal_cns$muted, linewidth = 0.40) +
    geom_segment(aes(x = lci95, xend = uci95, y = row, yend = row, color = edge_col), linewidth = 0.80) +
    geom_point(aes(fill = fill_col, color = edge_col, shape = pt_shape), size = 2.6, stroke = 0.78) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_shape_identity() +
    scale_x_log10(breaks = x_breaks, labels = x_labels) +
    coord_cartesian(xlim = xlim_vals) +
    y_sc +
    labs(title = if (show_headers) "<b>Interaction Ratio of ORs</b>" else NULL, x = x_lab, y = NULL) +
    theme_cns(8.0) +
    theme(
      axis.line.y  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x  = element_text(size = 6.9, color = pal_cns$subink),
      axis.title.x = ggtext::element_markdown(size = 7.3, color = pal_cns$ink, margin = margin(t = 2)),
      plot.title   = ggtext::element_markdown(size = 7.3, color = pal_cns$subink, hjust = 0.5, margin = margin(b = 2)),
      plot.margin  = margin(1, 4, 1, 4)
    )

  row_plot <- (
    col_text_md("model_md", "<b>Model &amp; Stratum</b>", "left") |
    col_text_md("arms_md",  "<b>Ctrl vs Durva+Ola</b>", "left") |
    p_mid |
    col_text_md("or_md",    "<b>Ratio OR (95% CI)</b>", "right") |
    col_text_md("p_md",     "<b><i>P</i></b>", "right") |
    col_text_md("q_md",     "<b><i>q</i> (m=8)</b>", "right")
  ) + plot_layout(widths = c(1.52, 0.80, 1.24, 1.15, 0.44, 0.46))

  if (!is.null(block_title)) {
    row_plot <- row_plot + plot_annotation(
      subtitle = block_title,
      theme = theme(
        plot.subtitle = element_text(size = 7.4, face = "bold", color = pal_cns$ink, margin = margin(t = 1, b = 1)),
        plot.margin = margin(0.5, 0, 0.5, 0)
      )
    )
  }
  wrap_elements(full = row_plot)
}

tests_full <- c(
  "pCR_Interaction_CXCL8_x_DurvaArm(Adj_HR)",
  "pCR_Interaction_PPARG_x_DurvaArm(Adj_HR)",
  "pCR_Interaction_CXCL8_x_DurvaArm(Adj_HR_MP)",
  "pCR_Interaction_PPARG_x_DurvaArm(Adj_HR_MP)"
)
tests_sub <- c(
  "pCR_Interaction_CXCL8_x_DurvaArm(HRneg_Subset)",
  "pCR_Interaction_CXCL8_x_DurvaArm(HRneg_Subset_Adj_MP)",
  "pCR_Interaction_PPARG_x_DurvaArm(HRneg_Subset)",
  "pCR_Interaction_PPARG_x_DurvaArm(HRneg_Subset_Adj_MP)"
)

p3a1 <- build_int_block_cns(
  int_all, tests_full, xlim_vals = c(0.05, 8.5),
  x_breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
  x_labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8"),
  show_headers = TRUE,
  block_title = "Full HER2\u2212 cohort (n = 100; 36 pCR)",
  highlight_primary = TRUE
)
p3a2 <- build_int_block_cns(
  int_all, tests_sub, xlim_vals = c(0.0004, 25),
  x_breaks = c(0.001, 0.01, 0.1, 1, 10),
  x_labels = c("0.001", "0.01", "0.1", "1", "10"),
  show_headers = FALSE,
  x_lab = "Ratio of treatment ORs per log<sub>2</sub> expression unit (log scale)",
  block_title = "Exploratory HR\u2212 subgroup (n = 38; 15 pCR)"
)

fig3_a <- wrap_elements(full = (p3a1 / p3a2) +
  plot_layout(heights = c(1.00, 1.12)) +
  plot_annotation(
    title = cns_panel_title("A", "I-SPY2 treatment interaction on pCR (GSE173839)"),
    theme = theme(
      plot.title = ggtext::element_markdown(margin = margin(b = 2)),
      plot.margin = margin(2, 2, 1, 2)
    )
  ))

# ------------------------------------------------------------------------------
# Panel B: Polished Raincloud Plots (clean y-axis dedup, luminous white-rimmed dots)
# ------------------------------------------------------------------------------
make_within_arm_panel_cns <- function(gene_name, population, sub_title_md, accent_col, show_ylab = TRUE) {
  d <- copy(g173)[is.finite(get(gene_name))]
  if (population == "TNBC_HRneg") d <- d[HR_neg == TRUE]

  arm_tbl <- data.table(
    Arm = c("control", "durvalumab/olaparib"),
    treatment_arm = c("Control_Paclitaxel_Arm", "Durvalumab_Olaparib_Arm"),
    arm_short = c("Control", "Durva+Ola")
  )
  tests <- th_dt[cohort == "GSE173839" & domain == "Neoadjuvant_pCR" &
                   stratum == population & gene == gene_name & treatment_arm %chin% arm_tbl$treatment_arm]
  tests <- tests[match(arm_tbl$treatment_arm, treatment_arm)]
  counts <- d[, .(n = .N, pcr = sum(pCR)), by = Arm][match(arm_tbl$Arm, Arm)]
  stopifnot(nrow(tests) == 2L, all(counts$n == tests$n_samples), all(counts$pcr == tests$n_events))

  facet_labs <- sprintf("%s (%d/%d)", arm_tbl$arm_short, counts$pcr, counts$n)
  d[, Arm_lbl := factor(Arm, levels = arm_tbl$Arm, labels = facet_labs)]
  d[, Resp := factor(pCR, levels = c(0L, 1L), labels = c("RD", "pCR"))]
  d[, expr := get(gene_name)]

  vio <- build_half_violin_dt(d, "Resp", "expr", facet_col = "Arm_lbl", width = 0.28, side = "right")
  yr <- range(c(d$expr, vio$vy))
  ys <- diff(yr)
  or_str <- ifelse(tests$estimate < 0.1, sprintf("OR = %.3f", tests$estimate), sprintf("OR = %.2f", tests$estimate))
  brk <- data.table(
    Arm_lbl = factor(facet_labs, levels = facet_labs),
    y = yr[2] + 0.08 * ys, tip = yr[2] + 0.03 * ys, y_txt = yr[2] + 0.14 * ys,
    label_md = sprintf("%s<br><i>P</i> = %.3g", or_str, tests$pval)
  )

  # Subtle background highlight on Durva+Ola facet
  durva_bg <- data.table(Arm_lbl = factor(facet_labs[2], levels = facet_labs))

  ggplot(d, aes(y = expr, fill = Resp, color = Resp)) +
    geom_rect(data = durva_bg, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
              inherit.aes = FALSE, fill = "#FAFBFC") +
    geom_hline(yintercept = Inf, color = pal_cns$ref_line, linewidth = 0.35) +
    geom_polygon(data = vio, aes(x = vx, y = vy, group = grp_id, fill = Resp, color = Resp),
                 inherit.aes = FALSE, alpha = 0.22, linewidth = 0.32) +
    ggbeeswarm::geom_quasirandom(aes(x = as.numeric(Resp) - 0.13, group = Resp, fill = Resp),
                                 shape = 21, color = "white", stroke = 0.22,
                                 width = 0.070, size = 1.50, alpha = 0.88) +
    geom_boxplot(aes(x = as.numeric(Resp) - 0.13, group = Resp), width = 0.13,
                 outlier.shape = NA, fill = "white", alpha = 0.80, linewidth = 0.36) +
    geom_segment(data = brk, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.32, color = pal_cns$ink) +
    geom_segment(data = brk, aes(x = 1, xend = 1, y = tip, yend = y), inherit.aes = FALSE, linewidth = 0.32, color = pal_cns$ink) +
    geom_segment(data = brk, aes(x = 2, xend = 2, y = tip, yend = y), inherit.aes = FALSE, linewidth = 0.32, color = pal_cns$ink) +
    ggtext::geom_richtext(data = brk, aes(x = 1.5, y = y_txt, label = label_md), inherit.aes = FALSE,
                          size = 2.02, lineheight = 1.12, vjust = 0, color = pal_cns$ink,
                          fill = NA, label.color = NA, label.padding = unit(0, "pt")) +
    facet_wrap(~ Arm_lbl) +
    scale_x_continuous(breaks = c(1, 2), labels = c("RD", "pCR"), limits = c(0.66, 2.36)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.24))) +
    scale_fill_manual(values = c("RD" = "#94A3B8", "pCR" = accent_col)) +
    scale_color_manual(values = c("RD" = "#64748B", "pCR" = accent_col)) +
    labs(
      title = sub_title_md,
      x = NULL,
      y = if (show_ylab) sprintf("<i>%s</i> expr. (log<sub>2</sub>)", gene_name) else NULL
    ) +
    theme_cns_grid(8.0) +
    theme(
      legend.position = "none",
      axis.title.y    = ggtext::element_markdown(size = 7.7, color = pal_cns$ink),
      strip.text      = ggtext::element_markdown(size = 7.0, color = pal_cns$subink, margin = margin(t = 1, b = 3)),
      plot.title      = ggtext::element_markdown(size = 8.1, color = pal_cns$ink, hjust = 0.5, margin = margin(b = 2)),
      panel.spacing   = unit(2.6, "mm"),
      plot.margin     = margin(2, 3, 2, 3)
    )
}

p3b1 <- make_within_arm_panel_cns("CXCL8", "HER2neg_All", "<i>CXCL8</i> (HER2\u2212 All, n = 100)", pal_cns$cxcl8, show_ylab = TRUE)
p3b2 <- make_within_arm_panel_cns("PPARG", "HER2neg_All", "<i>PPARG</i> (HER2\u2212 All, n = 100)", pal_cns$pparg, show_ylab = TRUE)
p3b3 <- make_within_arm_panel_cns("PPARG", "TNBC_HRneg",  "<i>PPARG</i> (HR\u2212 subgroup, n = 38)", pal_cns$pparg, show_ylab = FALSE)

fig3_b <- wrap_elements(full = (p3b1 | p3b2 | p3b3) +
  plot_layout(widths = c(1.04, 1.04, 0.94)) +
  plot_annotation(
    title = cns_panel_title("B", "Pretreatment tumor expression by treatment arm and pCR status"),
    theme = theme(
      plot.title = ggtext::element_markdown(margin = margin(b = 2)),
      plot.margin = margin(2, 2, 0, 2)
    )
  ))

# ------------------------------------------------------------------------------
# Panel C: Nature Medicine style ROC — diagonal-deviation veil + clean direct labels
#          ZERO cluttered header legends; 95% CI kept crisp in subtle second line
# ------------------------------------------------------------------------------
make_roc_panel_cns <- function(gene_name, population, sub_title_md, exp_col,
                               pos_durva, pos_ctrl, show_ylab = TRUE) {
  d <- copy(g173)[is.finite(get(gene_name))]
  if (population == "TNBC_HRneg") d <- d[HR_neg == TRUE]

  arm_tbl <- data.table(
    Arm = c("control", "durvalumab/olaparib"),
    treatment_arm = c("Control_Paclitaxel_Arm", "Durvalumab_Olaparib_Arm"),
    arm_lbl = c("Control", "Durva+Ola")
  )
  tests <- th_dt[cohort == "GSE173839" & domain == "Neoadjuvant_pCR" &
                   stratum == population & gene == gene_name & treatment_arm %chin% arm_tbl$treatment_arm]
  tests <- tests[match(arm_tbl$treatment_arm, treatment_arm)]

  roc_list <- lapply(seq_len(nrow(arm_tbl)), function(i) {
    sub <- d[Arm == arm_tbl$Arm[i]]
    rc <- compute_directional_roc(sub[[gene_name]], sub$pCR)
    rc[, arm_lbl := arm_tbl$arm_lbl[i]]
    rc
  })
  roc_dt <- rbindlist(roc_list)
  roc_dt[, arm_lbl := factor(arm_lbl, levels = arm_tbl$arm_lbl)]

  # Build polygon between Durva+Ola step curve and diagonal y = x (visualizes directional departure)
  durva_step <- roc_dt[arm_lbl == "Durva+Ola"]
  n_s <- nrow(durva_step)
  hv_x <- numeric(2 * n_s - 1); hv_y <- numeric(2 * n_s - 1)
  hv_x[1] <- durva_step$fpr[1]; hv_y[1] <- durva_step$tpr[1]
  if (n_s > 1) {
    for (j in 2:n_s) {
      hv_x[2*j - 2] <- durva_step$fpr[j]; hv_y[2*j - 2] <- durva_step$tpr[j - 1]
      hv_x[2*j - 1] <- durva_step$fpr[j]; hv_y[2*j - 1] <- durva_step$tpr[j]
    }
  }
  veil_poly <- data.table(x = c(hv_x, rev(hv_x)), y = c(hv_y, rev(hv_x)))

  ctrl_col <- "#94A3B8"
  col_map  <- c("Durva+Ola" = exp_col, "Control" = ctrl_col)
  lwd_map  <- c("Durva+Ola" = 1.15,    "Control" = 0.78)

  lbl_durva <- sprintf(
    "<span style='color:%s'><b>Durva+Ola</b><br>AUC %.3f</span><br><span style='color:#64748B; font-size:5.6pt'>(%.2f\u2013%.2f)</span>",
    exp_col, tests$auc[2], tests$auc_lci95[2], tests$auc_uci95[2]
  )
  lbl_ctrl <- sprintf(
    "<span style='color:#475569'>Control<br>AUC %.3f</span><br><span style='color:#94A3B8; font-size:5.6pt'>(%.2f\u2013%.2f)</span>",
    tests$auc[1], tests$auc_lci95[1], tests$auc_uci95[1]
  )

  ggplot(roc_dt, aes(x = fpr, y = tpr, color = arm_lbl, linewidth = arm_lbl)) +
    geom_polygon(data = veil_poly, aes(x = x, y = y), inherit.aes = FALSE,
                 fill = exp_col, alpha = 0.13) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#CBD5E1", linewidth = 0.42) +
    geom_step(direction = "hv", lineend = "round") +
    ggtext::geom_richtext(
      x = pos_durva[1], y = pos_durva[2], hjust = pos_durva[3], vjust = pos_durva[4],
      label = lbl_durva, size = 2.0, lineheight = 1.05,
      fill = NA, label.color = NA, label.padding = unit(0, "pt"), inherit.aes = FALSE
    ) +
    ggtext::geom_richtext(
      x = pos_ctrl[1], y = pos_ctrl[2], hjust = pos_ctrl[3], vjust = pos_ctrl[4],
      label = lbl_ctrl, size = 1.9, lineheight = 1.05,
      fill = NA, label.color = NA, label.padding = unit(0, "pt"), inherit.aes = FALSE
    ) +
    scale_color_manual(values = col_map, guide = "none") +
    scale_linewidth_manual(values = lwd_map, guide = "none") +
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0.015, 0.015)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0.015, 0.015)) +
    coord_fixed(ratio = 1) +
    labs(
      title = sub_title_md,
      x     = "1 \u2212 Specificity",
      y     = if (show_ylab) "Sensitivity" else NULL
    ) +
    theme_cns(8.0) +
    theme(
      panel.grid.major = element_line(color = "#F1F5F9", linewidth = 0.26),
      panel.border     = element_rect(color = "#E2E8F0", fill = NA, linewidth = 0.35),
      plot.title       = ggtext::element_markdown(size = 8.1, color = pal_cns$ink, hjust = 0.5, margin = margin(b = 3)),
      plot.margin      = margin(2, 5, 2, 5)
    )
}

# Hand-tuned open-space coordinates (x, y, hjust, vjust) — zero curve collision
p3c1 <- make_roc_panel_cns(
  "CXCL8", "HER2neg_All", "<i>CXCL8</i> (HER2\u2212 All)", pal_cns$cxcl8,
  pos_durva = c(0.03, 0.96, 0, 1),
  pos_ctrl  = c(0.97, 0.04, 1, 0),
  show_ylab = TRUE
)
p3c2 <- make_roc_panel_cns(
  "PPARG", "HER2neg_All", "<i>PPARG</i> (HER2\u2212 All)", pal_cns$pparg,
  pos_durva = c(0.97, 0.04, 1, 0),
  pos_ctrl  = c(0.03, 0.96, 0, 1),
  show_ylab = FALSE
)
p3c3 <- make_roc_panel_cns(
  "PPARG", "TNBC_HRneg", "<i>PPARG</i> (HR\u2212 subgroup)", pal_cns$pparg,
  pos_durva = c(0.97, 0.04, 1, 0),
  pos_ctrl  = c(0.03, 0.96, 0, 1),
  show_ylab = FALSE
)

fig3_c <- wrap_elements(full = (p3c1 | p3c2 | p3c3) +
  plot_layout(widths = c(1.06, 1, 1)) +
  plot_annotation(
    title = cns_panel_title("C", "Prespecified directional ROC curves for pCR"),
    theme = theme(
      plot.title = ggtext::element_markdown(margin = margin(b = 3)),
      plot.margin = margin(2, 2, 1, 2)
    )
  ))

fig3_cns <- wrap_plots(fig3_a, fig3_b, fig3_c, ncol = 1, heights = c(84, 64, 62))

save_cns_figure(fig3_cns, "Fig3_ISPY2_Treatment_Interaction_and_ROC", width_mm = 180, height_mm = 210, dpi_png = 350)
cat("CNS-Polished Fig 3 v3 complete.\n")
