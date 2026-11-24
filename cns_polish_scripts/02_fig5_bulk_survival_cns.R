# 02_fig5_bulk_survival_cns.R
# Figure 5 CNS Polish: Bulk tumor expression (Tumor vs Normal, PAM50),
# representative KM survival curves, and standardized OS random-effects meta-analysis.
# 100% identical statistical inputs, sample sizes, P values, and CIs to 05_revision_v2.

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "06_figures_cns_polish/01_scripts/00_cns_theme.R"))

cat("=== Generating CNS-Polished Fig 5 (Bulk Expression & Survival Meta) ===\n")

# ------------------------------------------------------------------------------
# Local self-contained helpers (built on shared 00_cns_theme.R palette & theme)
# ------------------------------------------------------------------------------
cns_panel_title <- function(letter, title_md, letter_pt = 11.2, title_pt = 8.5) {
  sprintf(
    "<span style='font-size:%gpt; color:%s;'><b>%s</b></span><span style='font-size:%gpt; color:%s;'>\u2002%s</span>",
    letter_pt, pal_cns$ink, letter, title_pt, pal_cns$ink, title_md
  )
}

format_p_fig5 <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-16 && p == 0) return("< 2.2e-16")
  if (p < 0.001) return(formatC(p, format = "e", digits = 1))
  sprintf("%.3g", p)
}

build_half_violin_dt <- function(dt, group_col, val_col, facet_col = NULL, width = 0.30, side = "right") {
  grp_cols <- if (is.null(facet_col)) group_col else c(facet_col, group_col)
  lvls <- levels(dt[[group_col]])
  if (is.null(lvls)) lvls <- unique(as.character(dt[[group_col]]))
  res <- dt[, {
    v <- get(val_col)
    v <- v[is.finite(v)]
    if (length(v) >= 5 && sd(v) > 0) {
      dns <- density(v, n = 128)
      sc  <- width / max(dns$y)
      x0  <- match(as.character(get(group_col)[1]), lvls)
      vx  <- if (side == "right") x0 + dns$y * sc else x0 - dns$y * sc
      data.table(vx = c(x0, vx, x0), vy = c(min(dns$x), dns$x, max(dns$x)))
    } else NULL
  }, by = grp_cols]
  res[, grp_id := do.call(paste, c(.SD, sep = "___")), .SDcols = grp_cols]
  res
}

survfit_to_step_cns <- function(fit, strata_levels) {
  s_names <- names(fit$strata)
  if (is.null(s_names)) s_names <- paste0(".g=", strata_levels[1])
  idx_start <- 1L
  out_list <- vector("list", length(s_names))
  cens_list <- vector("list", length(s_names))
  for (k in seq_along(s_names)) {
    n_k <- fit$strata[k]
    idx <- idx_start:(idx_start + n_k - 1L)
    idx_start <- idx_start + n_k
    g_lbl <- sub("^\\.g=", "", s_names[k])
    tt <- c(0, fit$time[idx]); ss <- c(1, fit$surv[idx])
    lo <- c(1, ifelse(is.na(fit$lower[idx]), fit$surv[idx], fit$lower[idx]))
    hi <- c(1, ifelse(is.na(fit$upper[idx]), fit$surv[idx], fit$upper[idx]))
    m <- length(tt)
    sx <- numeric(2 * m - 1); sy <- sx; slo <- sx; shi <- sx
    sx[1] <- tt[1]; sy[1] <- ss[1]; slo[1] <- lo[1]; shi[1] <- hi[1]
    if (m > 1) {
      for (j in 2:m) {
        p1 <- 2 * j - 2; p2 <- 2 * j - 1
        sx[p1] <- tt[j]; sx[p2] <- tt[j]
        sy[p1] <- ss[j - 1]; sy[p2] <- ss[j]
        slo[p1] <- lo[j - 1]; slo[p2] <- lo[j]
        shi[p1] <- hi[j - 1]; shi[p2] <- hi[j]
      }
    }
    out_list[[k]] <- data.table(time = sx, surv = sy, lower = slo, upper = shi,
                                strata = factor(g_lbl, levels = strata_levels))
    cm <- fit$n.censor[idx] > 0
    if (any(cm)) {
      cens_list[[k]] <- data.table(time = fit$time[idx][cm], surv = fit$surv[idx][cm],
                                   strata = factor(g_lbl, levels = strata_levels))
    }
  }
  list(steps = rbindlist(out_list), cens = rbindlist(cens_list))
}

plot_km_cns <- function(df, time_col, status_col, group_col,
                        panel_letter, title_md, xlab_txt, ylab_txt,
                        legend_title_md, legend_labs = c("Low", "High"),
                        pal_colors, break_by, max_x = NULL, base_size = 8.0) {
  d <- copy(df)[!is.na(get(time_col)) & get(time_col) > 0 & !is.na(get(status_col))]
  d[, .t := as.numeric(get(time_col))]
  d[, .e := as.integer(get(status_col))]
  d[, .g := factor(as.character(get(group_col)), levels = legend_labs)]
  if (is.null(max_x)) max_x <- ceiling(max(d$.t, na.rm = TRUE) / break_by) * break_by
  ticks <- seq(0, max_x, by = break_by)

  fit  <- survival::survfit(survival::Surv(.t, .e) ~ .g, data = as.data.frame(d))
  sd   <- survival::survdiff(survival::Surv(.t, .e) ~ .g, data = as.data.frame(d))
  km_p <- stats::pchisq(sd$chisq, df = 1, lower.tail = FALSE)
  cx   <- summary(survival::coxph(survival::Surv(.t, .e) ~ .g, data = as.data.frame(d)))
  hr   <- cx$conf.int[1, "exp(coef)"]; lci <- cx$conf.int[1, "lower .95"]; uci <- cx$conf.int[1, "upper .95"]

  step_dat <- survfit_to_step_cns(fit, legend_labs)
  st_df <- step_dat$steps[time <= max_x]
  cn_df <- step_dat$cens[time <= max_x]
  col_map <- setNames(pal_colors, legend_labs)

  stat_md <- sprintf("Log-rank <i>P</i> = %s<br>HR = %.2f (%.2f\u2013%.2f)", format_p_fig5(km_p), hr, lci, uci)

  p_curve <- ggplot(st_df, aes(x = time, y = surv, color = strata, fill = strata)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.72) +
    geom_point(data = cn_df, aes(x = time, y = surv, color = strata),
               shape = 3, size = 0.78, stroke = 0.36, alpha = 0.72, show.legend = FALSE) +
    ggtext::geom_richtext(
      data = data.table(x = max_x * 0.03, y = 0.06, lbl = stat_md),
      aes(x = x, y = y, label = lbl),
      inherit.aes = FALSE, hjust = 0, vjust = 0,
      size = 2.35, family = cns_font, color = pal_cns$ink,
      fill = NA, label.color = NA, label.padding = unit(0, "pt"),
      lineheight = 1.25
    ) +
    scale_color_manual(values = col_map, name = legend_title_md) +
    scale_fill_manual(values = col_map, guide = "none") +
    scale_x_continuous(breaks = ticks) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    coord_cartesian(xlim = c(-max_x * 0.04, max_x * 1.04), ylim = c(0, 1.03), expand = FALSE) +
    labs(title = cns_panel_title(panel_letter, title_md, letter_pt = 11.2, title_pt = 8.5),
         y = ylab_txt, x = NULL) +
    theme_cns(base_size) +
    theme(
      plot.title           = ggtext::element_markdown(hjust = 0, margin = margin(t = 1, b = 3, l = 0)),
      legend.position      = "top",
      legend.justification = "left",
      legend.margin        = margin(0, 0, -4, 0),
      legend.key.width     = unit(4.5, "mm"),
      legend.key.height    = unit(2.8, "mm"),
      legend.spacing.x     = unit(2.0, "mm"),
      plot.margin          = margin(4, 9, 1, 5),
      axis.title.x         = element_blank()
    )

  sm <- summary(fit, times = ticks, extend = TRUE)
  tbl_df <- data.table(
    time   = sm$time,
    strata = factor(sub("^\\.g=", "", as.character(sm$strata)), levels = rev(legend_labs)),
    nrisk  = sm$n.risk
  )
  tbl_df[, y_pos := ifelse(strata == legend_labs[1], 2.0, 1.0)]

  y_labs_md <- vapply(rev(legend_labs), function(lb) {
    sprintf("<span style='color:%s;'><b>%s</b></span>", col_map[[lb]], lb)
  }, character(1))

  p_tbl <- ggplot(tbl_df, aes(x = time, y = y_pos, label = nrisk, color = strata)) +
    geom_text(size = 2.25, family = cns_font, fontface = "plain", show.legend = FALSE) +
    scale_color_manual(values = col_map) +
    scale_x_continuous(breaks = ticks) +
    scale_y_continuous(breaks = c(1.0, 2.0), labels = y_labs_md) +
    coord_cartesian(xlim = c(-max_x * 0.04, max_x * 1.04), ylim = c(0.35, 2.65), expand = FALSE, clip = "off") +
    labs(title = "Number at risk", x = xlab_txt, y = NULL) +
    theme_cns(base_size) +
    theme(
      plot.title   = ggtext::element_markdown(size = base_size - 1.3, face = "plain", color = pal_cns$slate, hjust = 0, margin = margin(t = 1, b = 1, l = 0)),
      axis.text.y  = ggtext::element_markdown(size = base_size - 1.2, margin = margin(r = 5)),
      axis.line    = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks   = element_blank(),
      axis.title.x = ggtext::element_markdown(size = base_size - 0.5, color = pal_cns$ink, margin = margin(t = 2, b = 1)),
      plot.margin  = margin(0, 9, 3, 5)
    )

  patchwork::wrap_elements(full = (p_curve / p_tbl + plot_layout(heights = c(3.65, 1.08))))
}

table_forest_cns <- function(d, panel_letter, title_md, x_lab, null_val,
                             left_header, hr_header, p_header, i2_header,
                             widths = c(1.62, 1.25, 0.92, 0.48, 0.33),
                             log_x = TRUE, x_range = NULL, x_breaks = NULL,
                             point_shape = 18, point_size = 2.8,
                             base_size = 8.0) {
  d <- copy(d)
  d[, row := rev(seq_len(.N))]
  zebra <- d[row %% 2L == 1L]
  n_rows <- nrow(d)
  y_sc <- scale_y_continuous(limits = c(0.40, n_rows + 0.62), expand = c(0, 0))
  bands <- geom_rect(data = zebra, aes(ymin = row - 0.48, ymax = row + 0.48),
                     xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = pal_cns$zebra)
  top_rule <- geom_hline(yintercept = n_rows + 0.58, color = pal_cns$ref_line, linewidth = 0.38)

  side_col_md <- function(col_name, hdr, x_pos = 0.01, hjust_val = 0) {
    ggplot(d, aes(y = row)) +
      bands + top_rule +
      ggtext::geom_richtext(
        aes(x = x_pos, label = .data[[col_name]]),
        hjust = hjust_val, size = 2.38, family = cns_font, color = pal_cns$ink,
        fill = NA, label.color = NA, label.padding = unit(0, "pt")
      ) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + y_sc +
      labs(title = hdr) +
      theme_void(base_family = cns_font) +
      theme(
        plot.title = ggtext::element_markdown(
          size = base_size - 0.6, color = pal_cns$subink,
          hjust = hjust_val, margin = margin(b = 3, l = 2, r = 2)
        ),
        plot.margin = margin(2, 1, 2, 1)
      )
  }

  p_mid <- ggplot(d, aes(x = effect, y = row)) +
    bands + top_rule +
    geom_vline(xintercept = null_val, linetype = "dashed", color = pal_cns$ref_line, linewidth = 0.40) +
    geom_errorbar(aes(xmin = low, xmax = high, color = col_group),
                  orientation = "y", width = 0.18, linewidth = 0.50) +
    geom_point(aes(color = col_group), shape = point_shape, size = point_size) +
    scale_color_identity() +
    y_sc + labs(title = " ", x = x_lab, y = NULL) +
    theme_cns(base_size) +
    theme(
      axis.line.y  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x  = ggtext::element_markdown(size = base_size - 0.9, color = pal_cns$subink),
      axis.title.x = ggtext::element_markdown(size = base_size - 0.4, margin = margin(t = 3)),
      plot.title   = ggtext::element_markdown(size = base_size - 0.6, margin = margin(b = 3, l = 0)),
      plot.margin  = margin(2, 5, 2, 5)
    )

  if (log_x) {
    p_mid <- p_mid + scale_x_log10(breaks = x_breaks, labels = function(x) format(x, trim = TRUE, drop0trailing = TRUE))
  } else if (!is.null(x_breaks)) {
    p_mid <- p_mid + scale_x_continuous(breaks = x_breaks)
  }
  if (!is.null(x_range)) p_mid <- p_mid + coord_cartesian(xlim = x_range)

  out <- (
    side_col_md("row_label", left_header, x_pos = 0.01, hjust_val = 0) |
    p_mid |
    side_col_md("stat_hr", hr_header, x_pos = 0.04, hjust_val = 0) |
    side_col_md("stat_p",  p_header,  x_pos = 0.96, hjust_val = 1) |
    side_col_md("stat_i2", i2_header, x_pos = 0.96, hjust_val = 1)
  ) +
    plot_layout(widths = widths) +
    plot_annotation(
      title = cns_panel_title(panel_letter, title_md, letter_pt = 11.2, title_pt = 8.5),
      theme = theme(
        plot.title  = ggtext::element_markdown(hjust = 0, margin = margin(b = 3, l = 0)),
        plot.margin = margin(3, 2, 2, 2)
      )
    )
  patchwork::wrap_elements(full = out)
}

# ------------------------------------------------------------------------------
# Load read-only cached datasets & stats
# ------------------------------------------------------------------------------
tc      <- readRDS(file.path(bulk_dir, "tcga_brca_clean_cache.rds"))
mb      <- readRDS(file.path(bulk_dir, "metabric_target_cache.rds"))
g20     <- readRDS(file.path(geo_dir,  "GSE20685_parsed.rds"))
g42     <- readRDS(file.path(geo_dir,  "GSE42568_parsed.rds"))
cl_tc   <- tc$clin[match(colnames(tc$tumor_log2), tc$clin$SAMPLE_ID)]
diff_dt <- fread(file.path(res_dir, "module4_differential_expression.csv"))
s13     <- fread(file.path(res_dir, "planA_pam50_subtypes_stats.csv"))

# ------------------------------------------------------------------------------
# Panel A: Tumor vs. Normal PPARG expression
# ------------------------------------------------------------------------------
df_tn <- rbind(
  data.table(Cohort = "GSE42568",
             G_simple = ifelse(grepl("normal", as.character(g42$pheno$title), ignore.case = TRUE), "Normal", "Tumor"),
             Expr = as.numeric(g42$expr_target["PPARG", ])),
  data.table(Cohort = "TCGA-BRCA", G_simple = "Normal", Expr = as.numeric(tc$normal_log2["PPARG", ])),
  data.table(Cohort = "TCGA-BRCA", G_simple = "Tumor",  Expr = as.numeric(tc$tumor_log2["PPARG", ]))
)
df_tn[, G_simple := factor(G_simple, levels = c("Normal", "Tumor"))]
vio_5a <- build_half_violin_dt(df_tn, "G_simple", "Expr", facet_col = "Cohort", width = 0.30, side = "right")

p5a_gse <- diff_dt[cohort == "GSE42568"  & gene_symbol == "PPARG" & comparison == "Tumor_vs_Normal"]$wilcox_pval[1]
p5a_tcg <- diff_dt[cohort == "TCGA_BRCA" & gene_symbol == "PPARG" & comparison == "Tumor_vs_Normal"]$wilcox_pval[1]
brk_5a <- df_tn[, .(ymax = max(Expr, na.rm = TRUE), ymin = min(Expr, na.rm = TRUE)), by = Cohort]
brk_5a[, `:=`(
  y_brk  = ymax + 0.08 * (ymax - ymin),
  y_tick = ymax + 0.03 * (ymax - ymin),
  y_txt  = ymax + 0.145 * (ymax - ymin),
  p_md   = c(paste0("<i>P</i> = ", format_p_fig5(p5a_gse)),
             paste0("<i>P</i> = ", format_p_fig5(p5a_tcg)))
)]

p5a <- ggplot(df_tn, aes(x = G_simple, y = Expr, fill = G_simple, color = G_simple)) +
  geom_polygon(data = vio_5a, aes(x = vx, y = vy, group = grp_id, fill = G_simple, color = G_simple),
               inherit.aes = FALSE, alpha = 0.28, linewidth = 0.26) +
  ggrastr::rasterise(
    geom_jitter(aes(x = as.numeric(G_simple) - 0.17), width = 0.055, size = 0.42, alpha = 0.24),
    dpi = 320
  ) +
  geom_boxplot(aes(x = as.numeric(G_simple) - 0.17), width = 0.14, outlier.shape = NA,
               alpha = 0.78, color = pal_cns$ink, linewidth = 0.34) +
  geom_segment(data = brk_5a, aes(x = 1, xend = 2, y = y_brk, yend = y_brk),
               inherit.aes = FALSE, color = pal_cns$subink, linewidth = 0.35) +
  geom_segment(data = brk_5a, aes(x = 1, xend = 1, y = y_tick, yend = y_brk),
               inherit.aes = FALSE, color = pal_cns$subink, linewidth = 0.35) +
  geom_segment(data = brk_5a, aes(x = 2, xend = 2, y = y_tick, yend = y_brk),
               inherit.aes = FALSE, color = pal_cns$subink, linewidth = 0.35) +
  ggtext::geom_richtext(
    data = brk_5a, aes(x = 1.5, y = y_txt, label = p_md),
    inherit.aes = FALSE, size = 2.35, family = cns_font, color = pal_cns$ink,
    fill = NA, label.color = NA, label.padding = unit(0, "pt")
  ) +
  facet_wrap(~ Cohort, scales = "free_y") +
  scale_x_discrete(limits = c("Normal", "Tumor")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  scale_fill_manual(values = c("Normal" = pal_cns$pparg_light, "Tumor" = pal_cns$cxcl8)) +
  scale_color_manual(values = c("Normal" = pal_cns$pparg, "Tumor" = pal_cns$cxcl8_dark)) +
  labs(
    title = cns_panel_title("A", "<i>PPARG</i> expression (Tumor vs. Normal)"),
    x = NULL,
    y = "Normalized log<sub>2</sub> expr."
  ) +
  theme_cns(8.0) +
  theme(
    legend.position = "none",
    plot.title      = ggtext::element_markdown(hjust = 0, margin = margin(t = 1, b = 4, l = 0)),
    strip.text      = ggtext::element_markdown(face = "bold", size = 8.0, color = pal_cns$subink, hjust = 0.5, margin = margin(b = 4)),
    panel.spacing.x = unit(4.5, "mm")
  )

# ------------------------------------------------------------------------------
# Panel B: PAM50 molecular subtypes (TCGA-BRCA)
# ------------------------------------------------------------------------------
pam_tc <- rbind(
  data.table(Gene = "CXCL8", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["CXCL8", ])),
  data.table(Gene = "PPARG", Subtype = gsub("^BRCA_", "", as.character(cl_tc$SUBTYPE)), Expr = as.numeric(tc$tumor_log2["PPARG", ]))
)[Subtype %in% c("Normal", "LumA", "LumB", "Her2", "Basal")]
pam_tc[, Subtype := factor(Subtype, levels = c("Normal", "LumA", "LumB", "Her2", "Basal"),
                           labels = c("Norm.", "LumA", "LumB", "HER2+", "Basal"))]
pam_tc[, Gene_md := factor(ifelse(Gene == "CXCL8", "<i>CXCL8</i>", "<i>PPARG</i>"),
                           levels = c("<i>CXCL8</i>", "<i>PPARG</i>"))]
vio_5b <- build_half_violin_dt(pam_tc, "Subtype", "Expr", facet_col = "Gene_md", width = 0.30, side = "right")

kw_cx <- s13[cohort == "TCGA_BRCA" & gene_symbol == "CXCL8"]$kw_pval[1]
w_cx  <- wilcox.test(pam_tc[Gene == "CXCL8" & Subtype == "Basal"]$Expr, pam_tc[Gene == "CXCL8" & Subtype == "LumA"]$Expr)$p.value
kw_pp <- s13[cohort == "TCGA_BRCA" & gene_symbol == "PPARG"]$kw_pval[1]
w_pp  <- s13[cohort == "TCGA_BRCA" & gene_symbol == "PPARG"]$wilcox_Basal_vs_LumA[1]

ann_5b <- data.table(
  Gene_md = factor(c("<i>CXCL8</i>", "<i>PPARG</i>"), levels = c("<i>CXCL8</i>", "<i>PPARG</i>")),
  ymax    = c(max(pam_tc[Gene == "CXCL8"]$Expr, na.rm = TRUE), max(pam_tc[Gene == "PPARG"]$Expr, na.rm = TRUE)),
  ymin    = c(min(pam_tc[Gene == "CXCL8"]$Expr, na.rm = TRUE), min(pam_tc[Gene == "PPARG"]$Expr, na.rm = TRUE)),
  kw_md   = c(paste0("KW <i>P</i> = ", format_p_fig5(kw_cx)),
              paste0("KW <i>P</i> = ", format_p_fig5(kw_pp))),
  w_md    = c(paste0("Basal vs. LumA: <i>P</i> ", ifelse(w_cx < 2.2e-16, "< 2.2e-16", paste0("= ", format_p_fig5(w_cx)))),
              paste0("Basal vs. LumA: <i>P</i> = ", format_p_fig5(w_pp)))
)
ann_5b[, `:=`(
  y_brk  = ymax + 0.08 * (ymax - ymin),
  y_tick = ymax + 0.03 * (ymax - ymin),
  y_w    = ymax + 0.18 * (ymax - ymin),
  y_kw   = ymax + 0.33 * (ymax - ymin)
)]

pam_pal <- c(
  "Norm." = "#91A3B0",
  "LumA"  = pal_cns$pparg,
  "LumB"  = pal_cns$pparg_light,
  "HER2+" = pal_cns$teal,
  "Basal" = pal_cns$cxcl8
)

p5b <- ggplot(pam_tc, aes(x = Subtype, y = Expr, fill = Subtype)) +
  geom_polygon(data = vio_5b, aes(x = vx, y = vy, group = grp_id, fill = Subtype, color = Subtype),
               inherit.aes = FALSE, alpha = 0.28, linewidth = 0.24) +
  geom_boxplot(aes(x = as.numeric(Subtype) - 0.14), width = 0.19,
               outlier.size = 0.32, outlier.alpha = 0.45, outlier.stroke = 0.2,
               alpha = 0.80, color = pal_cns$ink, linewidth = 0.34) +
  geom_segment(data = ann_5b, aes(x = 2, xend = 5, y = y_brk, yend = y_brk),
               inherit.aes = FALSE, color = pal_cns$subink, linewidth = 0.35) +
  geom_segment(data = ann_5b, aes(x = 2, xend = 2, y = y_tick, yend = y_brk),
               inherit.aes = FALSE, color = pal_cns$subink, linewidth = 0.35) +
  geom_segment(data = ann_5b, aes(x = 5, xend = 5, y = y_tick, yend = y_brk),
               inherit.aes = FALSE, color = pal_cns$subink, linewidth = 0.35) +
  ggtext::geom_richtext(
    data = ann_5b, aes(x = 3.35, y = y_w, label = w_md),
    inherit.aes = FALSE, size = 2.02, family = cns_font, color = pal_cns$ink,
    fill = NA, label.color = NA, label.padding = unit(0, "pt")
  ) +
  ggtext::geom_richtext(
    data = ann_5b, aes(x = 3.0, y = y_kw, label = kw_md),
    inherit.aes = FALSE, size = 2.12, family = cns_font, color = pal_cns$subink,
    fill = NA, label.color = NA, label.padding = unit(0, "pt")
  ) +
  facet_wrap(~ Gene_md, scales = "free_y") +
  scale_x_discrete(limits = c("Norm.", "LumA", "LumB", "HER2+", "Basal")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.42))) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = pam_pal) +
  scale_color_manual(values = pam_pal) +
  labs(
    title = cns_panel_title("B", "PAM50 subtypes (TCGA-BRCA)"),
    x = NULL,
    y = "log<sub>2</sub>(RSEM + 1)"
  ) +
  theme_cns(8.0) +
  theme(
    legend.position = "none",
    plot.title      = ggtext::element_markdown(hjust = 0, margin = margin(t = 1, b = 4, l = 0)),
    strip.text      = ggtext::element_markdown(face = "bold", size = 8.4, color = pal_cns$ink, hjust = 0.5, margin = margin(b = 4)),
    axis.text.x     = ggtext::element_markdown(size = 6.9, angle = 22, hjust = 1),
    panel.spacing.x = unit(4.5, "mm")
  )

# ------------------------------------------------------------------------------
# Panels C–F: Representative Kaplan-Meier survival curves (italic gene names)
# ------------------------------------------------------------------------------
idx_tc_erp <- which(cl_tc$ER_Xena == "Positive" & !is.na(cl_tc$OS_MONTHS) & cl_tc$OS_MONTHS > 0 & !is.na(cl_tc$OS_event))
df_tc_p <- data.table(time = cl_tc$OS_MONTHS[idx_tc_erp], status = cl_tc$OS_event[idx_tc_erp],
                      expr = as.numeric(tc$tumor_log2["PPARG", idx_tc_erp]))
df_tc_p[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_5c <- plot_km_cns(
  df_tc_p, "time", "status", "Group",
  panel_letter = "C", title_md = "TCGA ER+ OS (<i>PPARG</i>)",
  xlab_txt = "Months", ylab_txt = "Overall survival",
  legend_title_md = "<i>PPARG</i>", legend_labs = c("Low", "High"),
  pal_colors = c(pal_cns$cxcl8, pal_cns$pparg), break_by = 50, max_x = 250
)

idx_mb_erp <- which(mb$clin$ER_clean == "Positive" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_p <- data.table(time = mb$clin$OS_MONTHS[idx_mb_erp], status = mb$clin$OS_event[idx_mb_erp],
                      expr = as.numeric(mb$expr_target["PPARG", idx_mb_erp]))
df_mb_p[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_5d <- plot_km_cns(
  df_mb_p, "time", "status", "Group",
  panel_letter = "D", title_md = "METABRIC ER+ OS (<i>PPARG</i>)",
  xlab_txt = "Months", ylab_txt = "Overall survival",
  legend_title_md = "<i>PPARG</i>", legend_labs = c("Low", "High"),
  pal_colors = c(pal_cns$cxcl8, pal_cns$pparg), break_by = 100, max_x = 300
)

pd20 <- as.data.table(g20$pheno)
df_g20_c <- data.table(
  time = suppressWarnings(as.numeric(as.character(pd20[["follow_up_duration (years):ch1"]]))),
  status = suppressWarnings(as.integer(as.character(pd20[["event_death:ch1"]]))),
  expr = as.numeric(g20$expr_target["CXCL8", ])
)[!is.na(time) & time > 0 & !is.na(status)]
df_g20_c[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_5e <- plot_km_cns(
  df_g20_c, "time", "status", "Group",
  panel_letter = "E", title_md = "GSE20685 OS (<i>CXCL8</i>)",
  xlab_txt = "Years", ylab_txt = "Overall survival",
  legend_title_md = "<i>CXCL8</i>", legend_labs = c("Low", "High"),
  pal_colors = c(pal_cns$pparg, pal_cns$cxcl8), break_by = 3, max_x = 15
)

idx_mb_ern <- which(mb$clin$ER_clean == "Negative" & !is.na(mb$clin$OS_MONTHS) & mb$clin$OS_MONTHS > 0 & !is.na(mb$clin$OS_event))
df_mb_c <- data.table(time = mb$clin$OS_MONTHS[idx_mb_ern], status = mb$clin$OS_event[idx_mb_ern],
                      expr = as.numeric(mb$expr_target["CXCL8", idx_mb_ern]))
df_mb_c[, Group := ifelse(expr >= median(expr), "High", "Low")]
km_5f <- plot_km_cns(
  df_mb_c, "time", "status", "Group",
  panel_letter = "F", title_md = "METABRIC ER\u2212 OS (<i>CXCL8</i>)",
  xlab_txt = "Months", ylab_txt = "Overall survival",
  legend_title_md = "<i>CXCL8</i>", legend_labs = c("Low", "High"),
  pal_colors = c(pal_cns$pparg, pal_cns$cxcl8), break_by = 100, max_x = 300
)
rm(tc, mb, g20, g42, cl_tc); gc()

# ------------------------------------------------------------------------------
# Panel G: Standardized OS random-effects meta-analysis (CNS table-forest)
# ------------------------------------------------------------------------------
meta_hr <- fread(file.path(res_dir, "planA_standardized_os_meta.csv"))
hr_sub <- meta_hr[(gene_symbol == "CXCL8" & stratum %chin% c("Overall", "ER_neg")) |
                  (gene_symbol == "PPARG" & stratum %chin% c("Overall", "ER_pos"))]
hr_sub[, `:=`(
  row_label = sprintf(
    "<b><i>%s</i></b>\u2002%s <span style='color:%s;'>(%d cohorts, <i>n</i> = %s)</span>",
    gene_symbol, clean_subtype_cns(stratum), pal_cns$slate, k, format(total_n, big.mark = ",", trim = TRUE)
  ),
  effect    = RE_HR,
  low       = RE_LCI,
  high      = RE_UCI,
  col_group = ifelse(gene_symbol == "CXCL8", pal_cns$cxcl8, pal_cns$pparg),
  stat_hr   = sprintf("%.3f (%.3f\u2013%.3f)", RE_HR, RE_LCI, RE_UCI),
  stat_p    = sprintf("%.3g", RE_pval),
  stat_i2   = sprintf("%.0f%%", I2_pct)
)]

p5g <- table_forest_cns(
  hr_sub,
  panel_letter = "G",
  title_md     = "Standardized overall survival meta-analysis (random effects)",
  x_lab        = "Hazard ratio per SD (95% Wald CI)",
  null_val     = 1,
  left_header  = "<b>Gene &amp; stratum</b>",
  hr_header    = "<b>HR (95% CI)</b>",
  p_header     = "<b><i>P</i></b>",
  i2_header    = "<b><i>I</i><sup>2</sup></b>",
  widths       = c(1.60, 1.25, 0.90, 0.50, 0.35),
  log_x        = TRUE,
  x_range      = c(0.72, 1.50),
  x_breaks     = c(0.8, 1.0, 1.2, 1.4),
  point_shape  = 18,
  point_size   = 2.8
)

# ------------------------------------------------------------------------------
# Assemble full Fig 5
# ------------------------------------------------------------------------------
fig5_r1  <- wrap_elements(full = (p5a | p5b) + plot_layout(widths = c(0.95, 1.15)))
fig5_all <- wrap_plots(
  fig5_r1,
  wrap_elements(full = (km_5c | km_5d)),
  wrap_elements(full = (km_5e | km_5f)),
  p5g,
  ncol = 1,
  heights = c(57, 63, 63, 47)
)

save_cns_figure(fig5_all, "Fig5_Bulk_Expression_and_Survival_Meta", width_mm = 180, height_mm = 230, dpi_png = 350)
cat("CNS-Polished Fig 5 complete.\n")
