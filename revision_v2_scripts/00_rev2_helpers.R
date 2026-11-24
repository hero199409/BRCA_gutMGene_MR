# 00_rev2_helpers.R
# Shared visual system for Manuscript Revision v2 (Frontiers in Immunology target)
# Clean minimalist top-journal aesthetic: zero prose subtitles or metadata footers inside graphics.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggbeeswarm)
  library(ggrastr)
  library(ragg)
  library(patchwork)
  library(cowplot)
  library(survival)
  library(jpeg)
  library(png)
  library(magick)
  library(grid)
  library(scales)
})

resolve_project_root <- function() {
  env_root <- Sys.getenv("BRCA_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
    return(normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
root_dir  <- resolve_project_root()
cache_dir <- file.path(root_dir, "00_data/mr_cache")
bulk_dir  <- file.path(root_dir, "00_data/bulk_cohorts")
geo_dir   <- file.path(root_dir, "00_data/geo_cache")
sc_dir    <- file.path(root_dir, "00_data/scRNA_cache")
res_dir   <- file.path(root_dir, "02_results")
ihc_dir   <- file.path(root_dir, "03_figures/HPA_IHC_Images")

rev_dir     <- file.path(root_dir, "05_revision_v2")
rev_fig_dir <- file.path(rev_dir, "03_figures")
rev_prv_dir <- file.path(rev_fig_dir, "previews")
rev_ms_dir  <- file.path(rev_dir, "04_manuscript")
dir.create(rev_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rev_prv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rev_ms_dir,  recursive = TRUE, showWarnings = FALSE)

pal_rev <- list(
  cxcl8      = "#C83E2B",
  cxcl8_dark = "#8E1F13",
  cxcl8_soft = "#FEE2E2",
  pparg      = "#2B5B84",
  pparg_dark = "#1A365D",
  pparg_soft = "#E0F2FE",
  teal       = "#0F766E",
  teal_soft  = "#CCFBF1",
  amber      = "#B45309",
  slate      = "#475569",
  muted      = "#94A3B8",
  zebra      = "#F8FAFC",
  border     = "#CBD5E1",
  ink        = "#0F172A",
  subink     = "#334155"
)

theme_rev <- function(base_size = 8.5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line         = element_line(color = pal_rev$ink, linewidth = 0.45),
      axis.ticks        = element_line(color = pal_rev$subink, linewidth = 0.38),
      axis.ticks.length = unit(2.4, "pt"),
      axis.title        = element_text(face = "bold", size = base_size - 0.3, color = pal_rev$ink),
      axis.text         = element_text(size = base_size - 1.0, color = pal_rev$subink),
      panel.grid        = element_blank(),
      strip.background  = element_rect(fill = "#F1F5F9", color = pal_rev$border, linewidth = 0.38),
      strip.text        = element_text(face = "bold", size = base_size - 0.5, color = pal_rev$ink,
                                       margin = margin(t = 2.5, b = 2.5, l = 3, r = 3)),
      plot.title        = element_text(face = "bold", size = base_size + 1.0, color = pal_rev$ink,
                                       hjust = 0, margin = margin(b = 3)),
      plot.title.position = "plot",
      plot.subtitle     = element_text(size = base_size - 0.9, color = pal_rev$subink,
                                       hjust = 0, margin = margin(b = 3), lineheight = 1.12),
      plot.caption      = element_text(size = base_size - 1.2, color = pal_rev$slate,
                                       hjust = 0, margin = margin(t = 2), lineheight = 1.12),
      plot.caption.position = "plot",
      legend.title      = element_text(face = "bold", size = base_size - 0.8, color = pal_rev$ink),
      legend.text       = element_text(size = base_size - 1.0, color = pal_rev$subink),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      plot.margin       = margin(4, 6, 4, 6)
    )
}

theme_rev_grid <- function(base_size = 8.5, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.border      = element_rect(color = pal_rev$border, linewidth = 0.45, fill = NA),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "#F1F5F9", linewidth = 0.32),
      axis.ticks        = element_line(color = pal_rev$slate, linewidth = 0.35),
      axis.ticks.length = unit(2.2, "pt"),
      axis.title        = element_text(face = "bold", size = base_size - 0.3, color = pal_rev$ink),
      axis.text         = element_text(size = base_size - 1.0, color = pal_rev$subink),
      strip.background  = element_rect(fill = "#F1F5F9", color = pal_rev$border, linewidth = 0.38),
      strip.text        = element_text(face = "bold", size = base_size - 0.5, color = pal_rev$ink,
                                       margin = margin(t = 2.5, b = 2.5, l = 3, r = 3)),
      plot.title        = element_text(face = "bold", size = base_size + 1.0, color = pal_rev$ink,
                                       hjust = 0, margin = margin(b = 3)),
      plot.title.position = "plot",
      plot.subtitle     = element_text(size = base_size - 0.9, color = pal_rev$subink,
                                       hjust = 0, margin = margin(b = 3), lineheight = 1.12),
      plot.caption      = element_text(size = base_size - 1.2, color = pal_rev$slate,
                                       hjust = 0, margin = margin(t = 2), lineheight = 1.12),
      plot.caption.position = "plot",
      legend.title      = element_text(face = "bold", size = base_size - 0.8, color = pal_rev$ink),
      legend.text       = element_text(size = base_size - 1.0, color = pal_rev$subink),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      plot.margin       = margin(4, 6, 4, 6)
    )
}

save_rev_figure <- function(p, name, width_mm = 180, height_mm = 220, dpi_png = 350, out_dir = rev_fig_dir) {
  w_in <- width_mm / 25.4
  h_in <- height_mm / 25.4
  png_path <- file.path(out_dir, paste0(name, ".png"))
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  ggsave(png_path, p, width = w_in, height = h_in, dpi = dpi_png, device = ragg::agg_png, bg = "white")
  ggsave(pdf_path, p, width = w_in, height = h_in, device = cairo_pdf, bg = "white")
  cat(sprintf("Saved %-48s | %.0fx%.0f mm | PNG=%.2f MB | PDF=%.2f MB\n",
              name, width_mm, height_mm,
              file.info(png_path)$size / 1048576,
              file.info(pdf_path)$size / 1048576))
}

clean_subtype_rev <- function(x) {
  x <- gsub("BRCA_", "", x)
  x <- gsub("ER_pos", "ER+", x)
  x <- gsub("ER_neg", "ER\u2212", x)
  x <- gsub("HER2_pos", "HER2+", x)
  x
}

clean_lineage_rev <- function(x) {
  m <- c(
    "Mono/Macro"     = "Myeloid",
    "T_cells"        = "T Cells",
    "CD4Tconv"       = "CD4+ T",
    "CD8T"           = "CD8+ T",
    "CD8Tex"         = "Exh. CD8+ T",
    "Tprolif"        = "Prolif. T",
    "B"              = "B Cells",
    "Plasma"         = "Plasma",
    "Epithelial"     = "Norm. Epithelial",
    "Malignant"      = "Malignant",
    "Myofibroblasts" = "Myofibroblasts",
    "Fibroblasts"    = "Fibroblasts",
    "Endothelial"    = "Endothelial",
    "PVL"            = "PVL",
    "Mast"           = "Mast",
    "DC"             = "DC"
  )
  out <- unname(m[x])
  ifelse(is.na(out), x, out)
}

format_p_clean <- function(p) {
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

compute_directional_roc <- function(expr, event) {
  ok <- is.finite(expr) & event %in% c(0L, 1L)
  x <- expr[ok]; y <- event[ok]
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  ord <- order(x, decreasing = TRUE)
  xs <- x[ord]; ys <- y[ord]
  rle_x <- rle(xs)
  ends <- cumsum(rle_x$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  tp <- cumsum(vapply(seq_along(ends), function(i) sum(ys[starts[i]:ends[i]] == 1L), integer(1)))
  fp <- cumsum(vapply(seq_along(ends), function(i) sum(ys[starts[i]:ends[i]] == 0L), integer(1)))
  data.table(fpr = c(0, fp / n0), tpr = c(0, tp / n1))
}

survfit_to_step_rev <- function(fit, strata_levels) {
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

plot_km_rev <- function(df, time_col, status_col, group_col,
                        title_txt, xlab_txt, ylab_txt,
                        legend_title, legend_labs = c("Low", "High"),
                        pal_colors, break_by, max_x = NULL, base_size = 8.2) {
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

  step_dat <- survfit_to_step_rev(fit, legend_labs)
  st_df <- step_dat$steps[time <= max_x]
  cn_df <- step_dat$cens[time <= max_x]
  col_map <- setNames(pal_colors, legend_labs)

  stat_lbl <- sprintf("Log-rank P = %s\nHR = %.2f (%.2f\u2013%.2f)", format_p_clean(km_p), hr, lci, uci)

  p_curve <- ggplot(st_df, aes(x = time, y = surv, color = strata, fill = strata)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.14, color = NA) +
    geom_line(linewidth = 0.85) +
    geom_point(data = cn_df, aes(x = time, y = surv, color = strata),
               shape = 3, size = 1.1, stroke = 0.55, show.legend = FALSE) +
    annotate("text", x = max_x * 0.03, y = 0.07, label = stat_lbl,
             hjust = 0, vjust = 0, size = 2.35, lineheight = 1.18, color = pal_rev$ink) +
    scale_color_manual(values = col_map, name = legend_title) +
    scale_fill_manual(values = col_map, guide = "none") +
    scale_x_continuous(breaks = ticks) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    coord_cartesian(xlim = c(-max_x * 0.04, max_x * 1.04), ylim = c(0, 1.03), expand = FALSE) +
    labs(title = title_txt, y = ylab_txt, x = NULL) +
    theme_rev(base_size) +
    theme(
      legend.position = "top",
      legend.margin   = margin(0, 0, -4, 0),
      legend.key.width = unit(4.5, "mm"),
      plot.title      = element_text(size = base_size + 0.3, margin = margin(b = 2)),
      plot.margin     = margin(3, 8, 1, 5),
      axis.title.x    = element_blank()
    )

  sm <- summary(fit, times = ticks, extend = TRUE)
  tbl_df <- data.table(
    time   = sm$time,
    strata = factor(sub("^\\.g=", "", as.character(sm$strata)), levels = rev(legend_labs)),
    nrisk  = sm$n.risk
  )
  tbl_df[, y_pos := ifelse(strata == legend_labs[1], 2.0, 1.0)]

  p_tbl <- ggplot(tbl_df, aes(x = time, y = y_pos, label = nrisk, color = strata)) +
    geom_text(size = 2.25, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = col_map) +
    scale_x_continuous(breaks = ticks) +
    scale_y_continuous(breaks = c(1.0, 2.0), labels = rev(legend_labs)) +
    coord_cartesian(xlim = c(-max_x * 0.04, max_x * 1.04), ylim = c(0.35, 2.65), expand = FALSE, clip = "off") +
    labs(title = "At risk", x = xlab_txt, y = NULL) +
    theme_rev(base_size) +
    theme(
      plot.title   = element_text(size = base_size - 1.2, face = "bold", color = pal_rev$subink, hjust = 0, margin = margin(t = 1, b = 1)),
      axis.text.y  = element_text(size = base_size - 1.3, face = "bold", color = unname(col_map[rev(legend_labs)]), margin = margin(r = 6)),
      axis.line    = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks   = element_blank(),
      axis.title.x = element_text(size = base_size - 0.6, face = "bold", color = pal_rev$ink, margin = margin(t = 2, b = 1)),
      plot.margin  = margin(0, 8, 3, 5)
    )

  patchwork::wrap_elements(full = (p_curve / p_tbl + plot_layout(heights = c(3.6, 1.05))))
}

table_forest_rev <- function(d, title_txt, x_lab, null_val, stats_header, widths,
                             log_x = TRUE, x_range = NULL, x_breaks = NULL,
                             point_col = pal_rev$pparg, base_size = 8.2) {
  d <- copy(d)
  d[, row := rev(seq_len(.N))]
  zebra <- d[row %% 2L == 1L]
  y_sc <- scale_y_continuous(limits = c(0.42, nrow(d) + 0.58), expand = c(0, 0))
  bands <- geom_rect(data = zebra, aes(ymin = row - 0.48, ymax = row + 0.48),
                     xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = pal_rev$zebra)

  side_col <- function(col_name, hdr) {
    ggplot(d, aes(y = row)) + bands +
      geom_text(aes(x = 0.01, label = .data[[col_name]]), hjust = 0, size = 2.35, color = pal_rev$ink) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + y_sc +
      labs(title = hdr) + theme_void(base_family = "sans") +
      theme(plot.title = element_text(size = base_size - 0.8, face = "bold", color = pal_rev$subink, margin = margin(b = 2)),
            plot.margin = margin(1, 2, 1, 2))
  }

  has_col_var <- "col_group" %in% names(d)
  p_mid <- ggplot(d, aes(x = effect, y = row)) + bands +
    geom_vline(xintercept = null_val, linetype = "dashed", color = pal_rev$muted, linewidth = 0.45)
  if (has_col_var) {
    p_mid <- p_mid +
      geom_errorbar(aes(xmin = low, xmax = high, color = col_group), orientation = "y", width = 0.24, linewidth = 0.60) +
      geom_point(aes(color = col_group), size = 2.0) +
      scale_color_identity()
  } else {
    p_mid <- p_mid +
      geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0.24, linewidth = 0.60, color = point_col) +
      geom_point(size = 2.0, color = point_col)
  }
  p_mid <- p_mid + y_sc + labs(title = " ", x = x_lab, y = NULL) + theme_rev(base_size) +
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = base_size - 1.1),
          axis.title.x = element_text(size = base_size - 0.6),
          plot.title = element_text(size = base_size - 0.8, margin = margin(b = 2)),
          plot.margin = margin(1, 4, 1, 4))
  if (log_x) {
    if (!is.null(x_breaks)) {
      p_mid <- p_mid + scale_x_log10(breaks = x_breaks, labels = function(x) format(x, trim = TRUE, drop0trailing = TRUE))
    } else {
      p_mid <- p_mid + scale_x_log10(breaks = scales::breaks_log(n = 5))
    }
  } else if (!is.null(x_breaks)) {
    p_mid <- p_mid + scale_x_continuous(breaks = x_breaks)
  }
  if (!is.null(x_range)) p_mid <- p_mid + coord_cartesian(xlim = x_range)

  out <- (side_col("row_label", " ") | p_mid | side_col("stats", stats_header)) +
    plot_layout(widths = widths)
  if (!is.null(title_txt) && nzchar(title_txt)) {
    out <- out + plot_annotation(
      title = title_txt,
      theme = theme(plot.title = element_text(size = base_size + 0.8, face = "bold", color = pal_rev$ink, margin = margin(b = 2)),
                    plot.margin = margin(2, 2, 2, 2))
    )
  }
  patchwork::wrap_elements(full = out)
}
