# 01b_fig8_new.R
# Figure 8 Redesign (Clean Top-Journal Minimalist Typography):
# Zero prose subtitles or footnote paragraphs inside the graphic; shared legend in C.

local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "05_revision_v2", "01_scripts")
  helper_path <- file.path(script_dir, "00_rev2_helpers.R")
  if (!file.exists(helper_path)) helper_path <- file.path("05_revision_v2", "01_scripts", "00_rev2_helpers.R")
  source(helper_path, local = FALSE)
})

cat("=== Generating Clean Minimalist Fig 8 ===\n")

th_dt  <- fread(file.path(res_dir, "planA_therapy_response_summary.csv"))
th_int <- fread(file.path(res_dir, "planA_therapy_interactions_and_multivar.csv"))
th_pc  <- readRDS(file.path(bulk_dir, "planA_therapy_plot_cache.rds"))

g173 <- copy(as.data.table(th_pc$gse173839))[pCR %in% c(0L, 1L)]
stopifnot(nrow(g173) == 100L, sum(g173$pCR) == 36L)

# ------------------------------------------------------------------------------
# Panel A: Hero Interaction Forest Plot (Split into Full HER2- vs Exploratory HR-)
# ------------------------------------------------------------------------------
int_all <- copy(th_int[cohort == "GSE173839_ISPY2" & est_type == "Ratio_OR"])
stopifnot(nrow(int_all) == 8L)

int_all[, gene := fifelse(grepl("CXCL8", test_name), "CXCL8", "PPARG")]
int_all[, tier := fcase(
  grepl("\\(Adj_HR\\)", test_name), "Primary (adj. HR)",
  grepl("\\(Adj_HR_MP\\)", test_name), "Sensitivity (adj. HR + MP)",
  grepl("\\(HRneg_Subset\\)", test_name), "Exploratory HR\u2212 (unadj.)",
  grepl("\\(HRneg_Subset_Adj_MP\\)", test_name), "Exploratory HR\u2212 (+ MP)"
)]
int_all[, model_lbl := fcase(
  test_name == "pCR_Interaction_CXCL8_x_DurvaArm(Adj_HR)",             "Primary: CXCL8 \u00d7 Arm (adj. HR)",
  test_name == "pCR_Interaction_PPARG_x_DurvaArm(Adj_HR)",             "Primary: PPARG \u00d7 Arm (adj. HR)",
  test_name == "pCR_Interaction_CXCL8_x_DurvaArm(Adj_HR_MP)",          "Sensitivity: CXCL8 \u00d7 Arm (adj. HR + MP)",
  test_name == "pCR_Interaction_PPARG_x_DurvaArm(Adj_HR_MP)",          "Sensitivity: PPARG \u00d7 Arm (adj. HR + MP)",
  test_name == "pCR_Interaction_CXCL8_x_DurvaArm(HRneg_Subset)",       "Subgroup: CXCL8 \u00d7 Arm (HR\u2212)",
  test_name == "pCR_Interaction_CXCL8_x_DurvaArm(HRneg_Subset_Adj_MP)", "Subgroup: CXCL8 \u00d7 Arm (HR\u2212 + MP)",
  test_name == "pCR_Interaction_PPARG_x_DurvaArm(HRneg_Subset)",       "Subgroup: PPARG \u00d7 Arm (HR\u2212)",
  test_name == "pCR_Interaction_PPARG_x_DurvaArm(HRneg_Subset_Adj_MP)", "Subgroup: PPARG \u00d7 Arm (HR\u2212 + MP)"
)]
int_all[, arms_txt := sprintf("%d/%d | %d/%d", events_control, n_control, events_durva, n_durva)]
int_all[, stats_txt := fifelse(
  estimate < 0.05,
  sprintf("%.4f [%.4f, %.3f] | %.3g | %.3g", estimate, lci95, uci95, pval, fdr_bh_all_ispy2_interactions),
  sprintf("%.3f [%.3f, %.3f] | %.3g | %.3g", estimate, lci95, uci95, pval, fdr_bh_all_ispy2_interactions)
)]

build_int_block <- function(dt_sub, order_tests, xlim_vals, x_breaks, x_labels,
                            show_headers = TRUE, x_lab = NULL, block_title = NULL,
                            highlight_primary = FALSE) {
  d <- copy(dt_sub)[match(order_tests, test_name)]
  d[, row := rev(seq_len(.N))]
  zebra <- d[row %% 2L == 1L]
  y_sc <- scale_y_continuous(limits = c(0.42, nrow(d) + 0.58), expand = c(0, 0))

  bg_layer <- if (highlight_primary) {
    list(
      geom_rect(data = zebra, aes(ymin = row - 0.48, ymax = row + 0.48),
                xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = "#F8FAFC"),
      geom_rect(data = d[test_name == "pCR_Interaction_CXCL8_x_DurvaArm(Adj_HR)"],
                aes(ymin = row - 0.48, ymax = row + 0.48),
                xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = "#FEF2F2")
    )
  } else {
    list(geom_rect(data = zebra, aes(ymin = row - 0.48, ymax = row + 0.48),
                   xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = "#F8FAFC"))
  }

  col_text <- function(col, header) {
    ggplot(d, aes(y = row)) + bg_layer +
      geom_text(aes(x = 0.01, label = .data[[col]]), hjust = 0, size = 2.35,
                fontface = ifelse(d$test_name == "pCR_Interaction_CXCL8_x_DurvaArm(Adj_HR)", "bold", "plain"),
                color = pal_rev$ink) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + y_sc +
      labs(title = if (show_headers) header else NULL) +
      theme_void(base_family = "sans") +
      theme(plot.title = element_text(size = 7.5, face = "bold", color = pal_rev$subink, margin = margin(b = 2)),
            plot.margin = margin(1, 2, 1, 2))
  }

  tier_shapes <- c(
    "Primary (adj. HR)" = 16,
    "Sensitivity (adj. HR + MP)" = 15,
    "Exploratory HR\u2212 (unadj.)" = 18,
    "Exploratory HR\u2212 (+ MP)" = 5
  )

  p_mid <- ggplot(d, aes(x = estimate, y = row, color = gene, shape = tier)) +
    bg_layer +
    geom_vline(xintercept = 1, linetype = "dashed", color = pal_rev$muted, linewidth = 0.45) +
    geom_errorbar(aes(xmin = lci95, xmax = uci95), orientation = "y", width = 0.24, linewidth = 0.68) +
    geom_point(size = 2.4, stroke = 0.78) +
    scale_x_log10(breaks = x_breaks, labels = x_labels) +
    coord_cartesian(xlim = xlim_vals) +
    y_sc +
    scale_color_manual(values = c("CXCL8" = pal_rev$cxcl8, "PPARG" = pal_rev$pparg), guide = "none") +
    scale_shape_manual(values = tier_shapes, guide = "none") +
    labs(title = if (show_headers) "Interaction Ratio of ORs" else NULL, x = x_lab, y = NULL) +
    theme_rev(8.2) +
    theme(
      axis.line.y  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x  = element_text(size = 7.2),
      axis.title.x = element_text(size = 7.7, margin = margin(t = 2)),
      plot.title   = element_text(size = 7.5, face = "bold", color = pal_rev$subink, hjust = 0.5, margin = margin(b = 2)),
      plot.margin  = margin(1, 4, 1, 4)
    )

  row_plot <- (col_text("model_lbl", "Model & Stratum") |
                 col_text("arms_txt", "Ctrl | Durva+Ola") |
                 p_mid |
                 col_text("stats_txt", "Ratio OR [95% CI] | P | q (m=8)")) +
    plot_layout(widths = c(1.58, 0.70, 1.25, 1.65))

  if (!is.null(block_title)) {
    row_plot <- row_plot + plot_annotation(
      subtitle = block_title,
      theme = theme(plot.subtitle = element_text(size = 7.6, face = "bold", color = pal_rev$ink, margin = margin(t = 2, b = 1)),
                    plot.margin = margin(1, 0, 1, 0))
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

p8a1 <- build_int_block(
  int_all, tests_full, xlim_vals = c(0.05, 8.5),
  x_breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
  x_labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8"),
  show_headers = TRUE,
  block_title = "Full HER2\u2212 cohort (n = 100; 36 pCR)",
  highlight_primary = TRUE
)
p8a2 <- build_int_block(
  int_all, tests_sub, xlim_vals = c(0.0004, 25),
  x_breaks = c(0.001, 0.01, 0.1, 1, 10),
  x_labels = c("0.001", "0.01", "0.1", "1", "10"),
  show_headers = FALSE,
  x_lab = "Ratio of treatment ORs per log2 expression unit (log scale)",
  block_title = "Exploratory HR\u2212 subgroup (n = 38; 15 pCR)"
)

fig8_a <- wrap_elements(full = (p8a1 / p8a2) +
  plot_layout(heights = c(1.02, 1.12)) +
  plot_annotation(
    title = "A   I-SPY2 treatment interaction on pCR (GSE173839)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 9.6, color = pal_rev$ink, margin = margin(b = 2)),
      plot.margin = margin(2, 2, 1, 2)
    )
  ))

# ------------------------------------------------------------------------------
# Panel B: Within-arm pretreatment bulk expression (B1, B2, B3)
# ------------------------------------------------------------------------------
make_within_arm_panel <- function(gene_name, population, panel_title, accent_col) {
  d <- copy(g173)[is.finite(get(gene_name))]
  if (population == "TNBC_HRneg") d <- d[HR_neg == TRUE]

  arm_tbl <- data.table(
    Arm = c("control", "durvalumab/olaparib"),
    treatment_arm = c("Control_Paclitaxel_Arm", "Durvalumab_Olaparib_Arm"),
    arm_short = c("Control (PTX)", "Durva+Ola+PTX")
  )
  tests <- th_dt[cohort == "GSE173839" & domain == "Neoadjuvant_pCR" &
                   stratum == population & gene == gene_name & treatment_arm %chin% arm_tbl$treatment_arm]
  tests <- tests[match(arm_tbl$treatment_arm, treatment_arm)]
  counts <- d[, .(n = .N, pcr = sum(pCR)), by = Arm][match(arm_tbl$Arm, Arm)]
  stopifnot(nrow(tests) == 2L, all(counts$n == tests$n_samples), all(counts$pcr == tests$n_events))

  facet_labs <- sprintf("%s\n(%d/%d pCR)", arm_tbl$arm_short, counts$pcr, counts$n)
  d[, Arm_lbl := factor(Arm, levels = arm_tbl$Arm, labels = facet_labs)]
  d[, Resp := factor(pCR, levels = c(0L, 1L), labels = c("RD", "pCR"))]
  d[, expr := get(gene_name)]

  vio <- build_half_violin_dt(d, "Resp", "expr", facet_col = "Arm_lbl", width = 0.26, side = "right")
  yr <- range(c(d$expr, vio$vy))
  ys <- diff(yr)
  or_str <- ifelse(tests$estimate < 0.1, sprintf("OR=%.3f", tests$estimate), sprintf("OR=%.2f", tests$estimate))
  brk <- data.table(
    Arm_lbl = factor(facet_labs, levels = facet_labs),
    y = yr[2] + 0.08 * ys, tip = yr[2] + 0.03 * ys, y_txt = yr[2] + 0.16 * ys,
    label = sprintf("%s\nP=%.3g", or_str, tests$pval)
  )

  ggplot(d, aes(y = expr, fill = Resp, color = Resp)) +
    geom_polygon(data = vio, aes(x = vx, y = vy, group = grp_id, fill = Resp),
                 inherit.aes = FALSE, alpha = 0.30, color = pal_rev$slate, linewidth = 0.28) +
    ggbeeswarm::geom_quasirandom(aes(x = as.numeric(Resp) - 0.12, group = Resp),
                                 width = 0.075, size = 1.0, alpha = 0.74) +
    geom_boxplot(aes(x = as.numeric(Resp) - 0.12, group = Resp), width = 0.15,
                 outlier.shape = NA, alpha = 0.65, color = pal_rev$ink, linewidth = 0.38) +
    geom_segment(data = brk, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.38, color = pal_rev$ink) +
    geom_segment(data = brk, aes(x = 1, xend = 1, y = tip, yend = y), inherit.aes = FALSE, linewidth = 0.38, color = pal_rev$ink) +
    geom_segment(data = brk, aes(x = 2, xend = 2, y = tip, yend = y), inherit.aes = FALSE, linewidth = 0.38, color = pal_rev$ink) +
    geom_text(data = brk, aes(x = 1.5, y = y_txt, label = label), inherit.aes = FALSE,
              size = 2.15, lineheight = 0.95, vjust = 0, color = pal_rev$ink) +
    facet_wrap(~ Arm_lbl) +
    scale_x_continuous(breaks = c(1, 2), labels = c("RD", "pCR"), limits = c(0.66, 2.36)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.22))) +
    scale_fill_manual(values = c("RD" = pal_rev$muted, "pCR" = accent_col)) +
    scale_color_manual(values = c("RD" = pal_rev$slate, "pCR" = accent_col)) +
    labs(title = panel_title, x = NULL, y = sprintf("%s expr. (log2)", gene_name)) +
    theme_rev(8.2) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 7.0, lineheight = 1.02),
      plot.title = element_text(size = 8.4),
      panel.spacing = unit(1.8, "mm"),
      plot.margin = margin(3, 3, 2, 3)
    )
}

p8b1 <- make_within_arm_panel("CXCL8", "HER2neg_All", "B1  CXCL8 (HER2\u2212 All, n=100)", pal_rev$cxcl8)
p8b2 <- make_within_arm_panel("PPARG", "HER2neg_All", "B2  PPARG (HER2\u2212 All, n=100)", pal_rev$pparg)
p8b3 <- make_within_arm_panel("PPARG", "TNBC_HRneg",  "B3  PPARG (HR\u2212 subgroup, n=38)", pal_rev$pparg)

# ------------------------------------------------------------------------------
# Panel C: Prespecified directional ROC curves (shared compact legend at bottom)
# ------------------------------------------------------------------------------
make_roc_panel <- function(gene_name, population, panel_title, exp_col) {
  d <- copy(g173)[is.finite(get(gene_name))]
  if (population == "TNBC_HRneg") d <- d[HR_neg == TRUE]

  arm_tbl <- data.table(
    Arm = c("durvalumab/olaparib", "control"),
    treatment_arm = c("Durvalumab_Olaparib_Arm", "Control_Paclitaxel_Arm"),
    arm_lbl = c("Durva + Ola + PTX", "Control (PTX)")
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

  lbl_txt <- sprintf(
    "Durva+Ola (%d/%d): %.3f (%.2f\u2013%.2f)\nControl (%d/%d): %.3f (%.2f\u2013%.2f)",
    tests$resp_n[1], tests$n_samples[1], tests$auc[1], tests$auc_lci95[1], tests$auc_uci95[1],
    tests$resp_n[2], tests$n_samples[2], tests$auc[2], tests$auc_lci95[2], tests$auc_uci95[2]
  )

  box_above <- tests$auc[1] < 0.45
  y_box <- if (box_above) 0.97 else 0.03
  v_box <- if (box_above) 1 else 0
  x_box <- if (box_above) 0.03 else 0.97
  h_box <- if (box_above) 0 else 1

  col_map <- c("Durva + Ola + PTX" = exp_col, "Control (PTX)" = pal_rev$slate)
  lty_map <- c("Durva + Ola + PTX" = "solid", "Control (PTX)" = "dashed")

  ggplot(roc_dt, aes(x = fpr, y = tpr, color = arm_lbl, linetype = arm_lbl)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = pal_rev$muted, linewidth = 0.45) +
    geom_step(linewidth = 0.90, direction = "hv") +
    annotate("label", x = x_box, y = y_box, hjust = h_box, vjust = v_box,
             label = lbl_txt, size = 1.92, lineheight = 1.08,
             fill = alpha("white", 0.92), color = pal_rev$ink, linewidth = 0.25,
             label.padding = unit(0.18, "lines")) +
    scale_color_manual(values = col_map, name = NULL) +
    scale_linetype_manual(values = lty_map, name = NULL) +
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0.01, 0.01)) +
    coord_fixed(ratio = 1) +
    labs(title = panel_title,
         subtitle = sprintf("Higher %s \u2192 pCR", gene_name),
         x = "1 \u2212 Specificity", y = "Sensitivity") +
    theme_rev_grid(8.2) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8.4),
      plot.subtitle = element_text(size = 7.1, margin = margin(b = 2)),
      plot.margin = margin(3, 4, 2, 4)
    )
}

p8c1 <- make_roc_panel("CXCL8", "HER2neg_All", "C1  CXCL8 ROC (HER2\u2212 All)", pal_rev$cxcl8)
p8c2 <- make_roc_panel("PPARG", "HER2neg_All", "C2  PPARG ROC (HER2\u2212 All)", pal_rev$pparg)
p8c3 <- make_roc_panel("PPARG", "TNBC_HRneg",  "C3  PPARG ROC (HR\u2212 subgroup)", pal_rev$pparg)

# Single crisp shared bottom legend for C1-C3
leg_dt <- data.table(
  x = c(1, 2, 3), y = 1,
  lbl = factor(c("Durva + Ola + PTX (CXCL8)", "Durva + Ola + PTX (PPARG)", "Control (PTX)"),
               levels = c("Durva + Ola + PTX (CXCL8)", "Durva + Ola + PTX (PPARG)", "Control (PTX)"))
)
p_leg <- ggplot(leg_dt, aes(x, y, color = lbl, linetype = lbl)) +
  geom_line(linewidth = 0.95) +
  scale_color_manual(values = c("Durva + Ola + PTX (CXCL8)" = pal_rev$cxcl8,
                                "Durva + Ola + PTX (PPARG)" = pal_rev$pparg,
                                "Control (PTX)" = pal_rev$slate), name = NULL) +
  scale_linetype_manual(values = c("Durva + Ola + PTX (CXCL8)" = "solid",
                                   "Durva + Ola + PTX (PPARG)" = "solid",
                                   "Control (PTX)" = "dashed"), name = NULL) +
  theme_void(base_family = "sans") +
  theme(legend.position = "bottom", legend.text = element_text(size = 7.6, color = pal_rev$ink),
        legend.key.width = unit(6, "mm"))
shared_leg <- cowplot::get_legend(p_leg)

fig8_b <- wrap_elements(full = (p8b1 | p8b2 | p8b3))
fig8_c <- wrap_elements(full = ((p8c1 | p8c2 | p8c3) / wrap_elements(full = shared_leg) +
                                  plot_layout(heights = c(1, 0.11))))

fig8_new <- wrap_plots(fig8_a, fig8_b, fig8_c, ncol = 1, heights = c(86, 62, 68))

save_rev_figure(fig8_new, "Fig8_Therapy_Response_ISPY2_and_GDSC2", width_mm = 180, height_mm = 216, dpi_png = 350)
cat("Clean Minimalist Fig 8 complete.\n")
