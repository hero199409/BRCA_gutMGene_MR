# 25_0_cns_plot_helpers.R
# CNS / Nature Medicine 级顶刊绘图公共函数库
# 特点：100% 保留原有统计结果与子图序号，升级为：
#   - 斑马条纹森林图（含完整 Shape 图例）
#   - 带 95% CI 半透明阴影带 + 删失竖线 + 双指标(Log-rank P & HR) + 均匀间距风险表 KM 曲线
#   - 纯 ggplot2 原生 Raincloud 半小提琴+散点箱线图
#   - 圆角卡片、平滑三次 Hermite 曲线、HPA 高倍局部放大插图（Zoom-in Inset）

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  library(ggbeeswarm)
  library(ggrastr)
  library(ragg)
  library(patchwork)
  library(cowplot)
  library(survival)
  library(rms)
  library(circlize)
  library(gridGraphics)
  library(jpeg)
  library(png)
  library(magick)
  library(grid)
  library(scales)
})

cache_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/mr_cache")
bulk_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/bulk_cohorts")
geo_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/geo_cache")
sc_dir    <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data/scRNA_cache")
res_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")
pub_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "03_figures/Publication_PlanA_CNS")
ihc_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "03_figures/HPA_IHC_Images")
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)

# CNS 顶刊双轴语义色板（ER-/CXCL8 珊瑚红、ER+/PPARG 藏青蓝、协同绿、紫、琥珀）
pal_cns <- c(
  red       = "#E64B35",
  red_dark  = "#B91C1C",
  red_soft  = "#FFE4E6",
  blue      = "#3C5488",
  blue_dark = "#1E3A8A",
  blue_soft = "#DBEAFE",
  cyan      = "#4DBBD5",
  green     = "#00A087",
  green_soft= "#CCFBF1",
  coral     = "#F39B7F",
  slate     = "#64748B",
  amber     = "#D97706",
  purple    = "#7C3AED",
  null      = "#94A3B8"
)

# 兼容别名
pal_npg <- pal_cns

#' 开放式 L 型经典期刊坐标轴主题
theme_cns <- function(base_size = 9.5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line         = element_line(color = "#0F172A", linewidth = 0.48),
      axis.ticks        = element_line(color = "#334155", linewidth = 0.40),
      axis.ticks.length = unit(2.8, "pt"),
      axis.title        = element_text(face = "bold", size = base_size - 0.2, color = "#0F172A"),
      axis.text         = element_text(color = "#334155", size = base_size - 1.0),
      panel.grid        = element_blank(),
      strip.background  = element_rect(fill = "#F1F5F9", color = NA),
      strip.text        = element_text(face = "bold", size = base_size - 0.3, color = "#0F172A", margin = margin(t = 3.5, b = 3.5)),
      plot.title        = element_text(face = "bold", size = base_size + 1.4, color = "#0F172A", margin = margin(b = 6)),
      plot.subtitle     = element_blank(),
      legend.title      = element_text(face = "bold", size = base_size - 0.7, color = "#0F172A"),
      legend.text       = element_text(size = base_size - 1.0, color = "#334155"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      plot.margin       = margin(6, 10, 6, 8)
    )
}

#' 轻量级发丝网格主题（用于二维分面矩阵、GSEA、基因座图）
theme_cns_grid <- function(base_size = 9.2, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.border      = element_rect(color = "#CBD5E1", linewidth = 0.50, fill = NA),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "#F1F5F9", linewidth = 0.35),
      axis.ticks        = element_line(color = "#64748B", linewidth = 0.35),
      axis.title        = element_text(face = "bold", size = base_size - 0.3, color = "#0F172A"),
      axis.text         = element_text(color = "#334155", size = base_size - 1.0),
      strip.background  = element_rect(fill = "#F1F5F9", color = "#CBD5E1", linewidth = 0.45),
      strip.text        = element_text(face = "bold", size = base_size - 0.3, color = "#0F172A"),
      plot.title        = element_text(face = "bold", size = base_size + 1.2, color = "#0F172A", margin = margin(b = 6)),
      plot.subtitle     = element_blank(),
      legend.title      = element_text(face = "bold", size = base_size - 0.7, color = "#0F172A"),
      legend.text       = element_text(size = base_size - 1.0, color = "#334155"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      plot.margin       = margin(6, 8, 6, 6)
    )
}

save_cns <- function(p, name, w = 15, h = 8, dpi_png = 350) {
  png_path <- file.path(pub_dir, paste0(name, ".png"))
  pdf_path <- file.path(pub_dir, paste0(name, ".pdf"))
  ggsave(png_path, p, width = w, height = h, dpi = dpi_png, device = ragg::agg_png, bg = "white")
  ggsave(pdf_path, p, width = w, height = h, device = cairo_pdf, bg = "white")
  png_mb <- file.info(png_path)$size / 1048576
  pdf_mb <- file.info(pdf_path)$size / 1048576
  cat(sprintf("Exported %-46s | PNG=%.2f MB | PDF=%.2f MB\n", name, png_mb, pdf_mb))
  if (pdf_mb >= 10) warning(sprintf("PDF >= 10 MB (%s: %.2f MB)!", name, pdf_mb))
}

# 仅 Fig 1–3 使用的 Frontiers 双栏规格。保留其他图件调用的 save_cns 默认值。
save_frontiers_cns <- function(p, name, height_mm) {
  save_cns(p, name, w = 180 / 25.4, h = height_mm / 25.4, dpi_png = 600)
}

theme_frontiers_cns <- function() {
  theme_cns(9) + theme(
    plot.title = element_text(size = 10, margin = margin(b = 4)),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.subtitle = element_text(size = 8, color = "#334155", margin = margin(b = 4)),
    plot.caption = element_text(size = 8, color = "#334155", hjust = 0, margin = margin(t = 4)),
    panel.spacing = unit(3, "mm"),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(3.5, "mm"),
    plot.margin = margin(4, 6, 4, 6)
  )
}

# 详列数值时单独预留右侧列，避免随区间端点变化而遮盖点或截断。
forest_frontiers_cns <- function(d, title, subtitle, colors, shapes, legend_rows = 1) {
  stopifnot(nrow(d) > 0, !anyDuplicated(d$label), all(is.finite(d$OR_LCI95)),
            all(is.finite(d$OR_UCI95)), all(d$OR_LCI95 > 0))
  d <- copy(d)
  d[, label := factor(label, levels = rev(label))]
  lo <- min(log(d$OR_LCI95), 0)
  hi <- max(log(d$OR_UCI95), 0)
  span <- hi - lo
  ticks <- scales::breaks_extended(4)(exp(c(lo, hi)))
  ticks <- sort(c(1, ticks[ticks >= exp(lo) & ticks <= exp(hi) & abs(log(ticks)) >= 0.08 * span]))
  d[, annotation_x := exp(hi + 0.07 * span)]
  ggplot(d, aes(x = OR, y = label)) +
    geom_vline(xintercept = 1, color = "#94A3B8", linetype = "dashed", linewidth = 0.4) +
    geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95, color = tier),
                  orientation = "y", width = 0.24, linewidth = 0.55) +
    geom_point(aes(color = tier, shape = instrument), size = 2.2) +
    geom_text(aes(x = annotation_x, label = annotation), hjust = 0,
              size = 8 / ggplot2::.pt, color = "#0F172A", show.legend = FALSE) +
    scale_color_manual(values = colors, name = NULL, guide = guide_legend(order = 1, nrow = legend_rows)) +
    scale_shape_manual(values = shapes, name = NULL, guide = guide_legend(order = 2)) +
    scale_x_log10(breaks = ticks, labels = label_number(accuracy = NULL)) +
    coord_cartesian(xlim = exp(c(lo - 0.07 * span, hi + 0.90 * span))) +
    labs(title = title, subtitle = subtitle, x = "Odds ratio (log scale)", y = NULL) +
    theme_frontiers_cns() + theme(legend.box = "vertical")
}

clean_subtype_str <- function(x) {
  x <- gsub("BRCA_", "", x)
  x <- gsub("ER_pos", "ER+", x)
  x <- gsub("ER_neg", "ER\u2212", x)
  x <- gsub("HER2_pos", "HER2+", x)
  x
}

clean_lineage_str <- function(x) {
  m <- c(
    "Mono/Macro"     = "Mono/Macrophages",
    "T_cells"        = "T Cells",
    "CD4Tconv"       = "CD4+ T Cells",
    "CD8T"           = "CD8+ T Cells",
    "CD8Tex"         = "Exhausted CD8+ T",
    "Tprolif"        = "Prolif. T Cells",
    "B"              = "B Cells",
    "Plasma"         = "Plasma Cells",
    "Epithelial"     = "Normal Epithelial",
    "Malignant"      = "Malignant Cells",
    "Myofibroblasts" = "Myofibroblasts",
    "Fibroblasts"    = "Fibroblasts",
    "Endothelial"    = "Endothelial",
    "PVL"            = "PVL",
    "Mast"           = "Mast Cells",
    "DC"             = "Dendritic Cells",
    "SMC"            = "Smooth Muscle"
  )
  out <- unname(m[x])
  ifelse(is.na(out), x, out)
}

format_pval_expr <- function(p, prefix = "italic(P)") {
  if (is.na(p)) return(paste0(prefix, " == NA"))
  if (p < 1e-16 && p == 0) return(sprintf("%s < 2.2 %%*%% 10^{-16}", prefix))
  if (p < 0.001) {
    s <- formatC(p, format = "e", digits = 1)
    parts <- strsplit(s, "e")[[1]]
    mant <- as.numeric(parts[1])
    expo <- as.integer(parts[2])
    sprintf("%s == %.1f %%*%% 10^{%d}", prefix, mant, expo)
  } else {
    sprintf("%s == %.3f", prefix, p)
  }
}

#' 斑马条纹森林图（支持多形状图例与右对齐标注列）
plot_forest_cns <- function(dt, title_txt, xlim_vals, color_map,
                            text_x_mult = 1.035, base_size = 9.4,
                            show_shape_legend = FALSE) {
  d <- copy(dt)
  lvls <- levels(d$label)
  if (is.null(lvls)) lvls <- unique(as.character(d$label))
  d[, y_num := as.numeric(factor(label, levels = lvls))]
  zebra_dt <- data.table(y_num = seq_along(lvls))[y_num %% 2 == 1]

  shape_vals <- c("Primary" = 16, "Warning (2 SNPs)" = 17, "LD-Confounded (Excluded)" = 4)
  present_shapes <- intersect(names(shape_vals), unique(as.character(d$domain_flag)))
  shape_guide <- if (show_shape_legend && length(present_shapes) > 1) {
    guide_legend(order = 2, override.aes = list(size = 2.8, color = "#334155"))
  } else "none"

  ggplot(d, aes(x = OR, y = label, color = sig_tag)) +
    geom_rect(data = zebra_dt,
              aes(ymin = y_num - 0.48, ymax = y_num + 0.48, xmin = -Inf, xmax = Inf),
              inherit.aes = FALSE, fill = "#F8FAFC", color = NA) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
    geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.28, linewidth = 0.78) +
    geom_point(aes(shape = domain_flag), size = 2.95) +
    geom_text(aes(x = ifelse(OR_UCI95 < 1.01, pmax(OR_UCI95 * text_x_mult, 1.032), OR_UCI95 * text_x_mult),
                  label = anno_txt),
              hjust = 0, size = 2.65, color = "#0F172A", show.legend = FALSE) +
    scale_color_manual(values = color_map, guide = guide_legend(order = 1)) +
    scale_shape_manual(values = shape_vals, breaks = present_shapes, guide = shape_guide) +
    coord_cartesian(xlim = xlim_vals) +
    labs(title = title_txt, x = "Odds Ratio (95% CI)", y = NULL, color = NULL, shape = NULL) +
    theme_cns(base_size) +
    theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(t = -2))
}

#' 将 survfit 对象展开为严格右连续阶梯坐标，供 95% CI ribbon 与删失竖线使用
survfit_to_step_dt <- function(fit, strata_levels) {
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

    tt <- c(0, fit$time[idx])
    ss <- c(1, fit$surv[idx])
    lo <- c(1, ifelse(is.na(fit$lower[idx]), fit$surv[idx], fit$lower[idx]))
    hi <- c(1, ifelse(is.na(fit$upper[idx]), fit$surv[idx], fit$upper[idx]))

    # 展开成阶梯线段节点 (t_i, s_{i-1}) -> (t_i, s_i)
    m <- length(tt)
    sx <- numeric(2 * m - 1)
    sy <- numeric(2 * m - 1); slo <- numeric(2 * m - 1); shi <- numeric(2 * m - 1)
    sx[1] <- tt[1]; sy[1] <- ss[1]; slo[1] <- lo[1]; shi[1] <- hi[1]
    if (m > 1) {
      for (j in 2:m) {
        p1 <- 2 * j - 2; p2 <- 2 * j - 1
        sx[p1]  <- tt[j]; sx[p2]  <- tt[j]
        sy[p1]  <- ss[j - 1]; sy[p2]  <- ss[j]
        slo[p1] <- lo[j - 1]; slo[p2] <- lo[j]
        shi[p1] <- hi[j - 1]; shi[p2] <- hi[j]
      }
    }
    out_list[[k]] <- data.table(time = sx, surv = sy, lower = slo, upper = shi,
                                strata = factor(g_lbl, levels = strata_levels))
    c_mask <- fit$n.censor[idx] > 0
    if (any(c_mask)) {
      cens_list[[k]] <- data.table(time = fit$time[idx][c_mask], surv = fit$surv[idx][c_mask],
                                   strata = factor(g_lbl, levels = strata_levels))
    }
  }
  list(steps = rbindlist(out_list), cens = rbindlist(cens_list))
}

#' CNS 升级版 KM 生存曲线 + 95% CI 阴影带 + Log-rank P & HR + 均匀行距风险人数表
plot_km_cns <- function(df, time_col, status_col, group_col,
                        title_txt, xlab_txt, ylab_txt,
                        legend_title, legend_labs = c("Low", "High"),
                        pal_colors, pval_coord, break_by, max_x = NULL,
                        base_size = 9.0) {
  d <- copy(df)[!is.na(get(time_col)) & get(time_col) > 0 & !is.na(get(status_col))]
  d[, .t := as.numeric(get(time_col))]
  d[, .e := as.integer(get(status_col))]
  d[, .g := factor(as.character(get(group_col)), levels = legend_labs)]
  if (is.null(max_x)) max_x <- ceiling(max(d$.t, na.rm = TRUE) / break_by) * break_by
  ticks <- seq(0, max_x, by = break_by)

  fit <- survival::survfit(survival::Surv(.t, .e) ~ .g, data = as.data.frame(d))
  sd  <- survival::survdiff(survival::Surv(.t, .e) ~ .g, data = as.data.frame(d))
  km_p <- stats::pchisq(sd$chisq, df = 1, lower.tail = FALSE)
  cx   <- summary(survival::coxph(survival::Surv(.t, .e) ~ .g, data = as.data.frame(d)))
  hr   <- cx$conf.int[1, "exp(coef)"]
  lci  <- cx$conf.int[1, "lower .95"]
  uci  <- cx$conf.int[1, "upper .95"]

  step_dat <- survfit_to_step_dt(fit, legend_labs)
  st_df <- step_dat$steps[time <= max_x]
  cn_df <- step_dat$cens[time <= max_x]
  col_map <- setNames(pal_colors, legend_labs)

  p_expr <- format_pval_expr(km_p, prefix = '"Log-rank"~italic(P)')
  hr_txt <- sprintf("HR = %.2f (95%% CI: %.2f\u2013%.2f)", hr, lci, uci)

  p_curve <- ggplot(st_df, aes(x = time, y = surv, color = strata, fill = strata)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.14, color = NA) +
    geom_line(linewidth = 0.92) +
    geom_point(data = cn_df, aes(x = time, y = surv, color = strata), shape = 3, size = 1.35, stroke = 0.65, show.legend = FALSE) +
    annotate("text", x = pval_coord[1], y = pval_coord[2] + 0.085, label = p_expr,
             parse = TRUE, hjust = 0, size = 3.05, color = "#0F172A") +
    annotate("text", x = pval_coord[1], y = pval_coord[2], label = hr_txt,
             hjust = 0, size = 2.85, fontface = "bold", color = "#334155") +
    scale_color_manual(values = col_map, name = legend_title) +
    scale_fill_manual(values = col_map, guide = "none") +
    scale_x_continuous(breaks = ticks) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    coord_cartesian(xlim = c(-max_x * 0.06, max_x * 1.04), ylim = c(0, 1.03), expand = FALSE) +
    labs(title = title_txt, y = ylab_txt, x = NULL) +
    theme_cns(base_size) +
    theme(
      legend.position = "top",
      legend.margin   = margin(0, 0, -4, 0),
      plot.margin     = margin(6, 12, 2, 6),
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
    geom_text(size = 2.75, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = col_map) +
    scale_x_continuous(breaks = ticks) +
    scale_y_continuous(breaks = c(1.0, 2.0), labels = rev(legend_labs)) +
    coord_cartesian(xlim = c(-max_x * 0.06, max_x * 1.04), ylim = c(0.42, 2.58), expand = FALSE, clip = "off") +
    labs(title = "Number at risk", x = xlab_txt, y = NULL) +
    theme_cns(base_size) +
    theme(
      plot.title   = element_text(size = base_size - 0.8, face = "bold", color = "#334155", hjust = 0, margin = margin(t = 2, b = 3)),
      axis.text.y  = element_text(size = base_size - 1.0, face = "bold", color = unname(col_map[rev(legend_labs)]), margin = margin(r = 8)),
      axis.line    = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = base_size - 0.4, face = "bold", color = "#0F172A", margin = margin(t = 5, b = 2)),
      plot.margin  = margin(1, 12, 5, 6)
    )

  patchwork::wrap_elements(full = (p_curve / p_tbl + plot_layout(heights = c(4.1, 0.82))))
}

#' 原生 ggplot2 计算右半侧小提琴密度多边形（用于 Raincloud 图，零额外包依赖）
build_half_violin_dt <- function(dt, group_col, val_col, facet_col = NULL, width = 0.36, side = "right") {
  grp_cols <- if (is.null(facet_col)) group_col else c(facet_col, group_col)
  lvls <- levels(dt[[group_col]])
  if (is.null(lvls)) lvls <- unique(as.character(dt[[group_col]]))
  res <- dt[, {
    v <- get(val_col)
    v <- v[!is.na(v)]
    if (length(v) >= 5 && sd(v) > 0) {
      dns <- density(v, n = 128)
      sc  <- width / max(dns$y)
      x0  <- match(as.character(get(group_col)[1]), lvls)
      vx  <- if (side == "right") x0 + dns$y * sc else x0 - dns$y * sc
      # 闭合回基线
      data.table(vx = c(x0, vx, x0), vy = c(min(dns$x), dns$x, max(dns$x)))
    } else NULL
  }, by = grp_cols]
  res[, grp_id := do.call(paste, c(.SD, sep = "___")), .SDcols = grp_cols]
  res
}

#' 平滑三次 Hermite S 形连线
make_sigmoid_edges <- function(edges_df, direction = "x", n_pts = 60) {
  t <- seq(0, 1, length.out = n_pts)
  s <- 3 * t^2 - 2 * t^3
  res <- vector("list", nrow(edges_df))
  for (i in seq_len(nrow(edges_df))) {
    row <- edges_df[i]
    px <- if (direction == "x") row$x1 + t * (row$x2 - row$x1) else row$x1 + s * (row$x2 - row$x1)
    py <- if (direction == "x") row$y1 + s * (row$y2 - row$y1) else row$y1 + t * (row$y2 - row$y1)
    res[[i]] <- data.table(edge_id = i, x = px, y = py,
                           edge_type = if ("edge_type" %in% names(row)) row$edge_type else "flow",
                           is_core   = if ("is_core"   %in% names(row)) row$is_core   else FALSE)
  }
  rbindlist(res)
}

add_round_card <- function(p, xmin, xmax, ymin, ymax, fill, border_col = NA, lwd = 0.6, r_pt = 7) {
  g <- grid::roundrectGrob(
    r = unit(r_pt, "pt"),
    gp = grid::gpar(fill = fill, col = ifelse(is.na(border_col), NA, border_col), lwd = lwd)
  )
  p + annotation_custom(g, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

#' 构建带高倍局部放大镜插图 (Zoom-in Inset) 的 HPA IHC 病理图板
make_hpa_inset_panel <- function(img_path, crop_box, title_lbl, caption_txt,
                                 accent_col = "#3C5488", max_px = 1150) {
  mg_full <- magick::image_read(img_path)
  info <- magick::image_info(mg_full)
  W <- info$width; H <- info$height

  # 1) 全景圆形 TMA 缩略图
  mg_main <- magick::image_scale(mg_full, sprintf("%dx%d", max_px, max_px))
  tmp1 <- tempfile(fileext = ".jpg")
  magick::image_write(mg_main, tmp1, format = "jpeg", quality = 92)
  arr_main <- jpeg::readJPEG(tmp1)
  unlink(tmp1)

  # 2) 从 3000x3000 高清原图裁剪局部区域 (crop_box = c(xmin_frac, xmax_frac, ymin_frac, ymax_frac) in image pixel coords top-to-bottom)
  x0 <- round(crop_box[1] * W); cw <- round((crop_box[2] - crop_box[1]) * W)
  y0 <- round(crop_box[3] * H); ch <- round((crop_box[4] - crop_box[3]) * H)
  mg_zoom <- magick::image_crop(mg_full, sprintf("%dx%d+%d+%d", cw, ch, x0, y0))
  mg_zoom <- magick::image_scale(mg_zoom, "680x680")
  tmp2 <- tempfile(fileext = ".jpg")
  magick::image_write(mg_zoom, tmp2, format = "jpeg", quality = 94)
  arr_zoom <- jpeg::readJPEG(tmp2)
  unlink(tmp2)

  # 在 1:1 等比例坐标系 (xlim = [0, 1.54], ylim = [0, 1.02]) 中绘制，保证 TMA 圆片与局部放大框均严格 1:1 不变形
  g_main <- grid::rasterGrob(arr_main, interpolate = TRUE, width = unit(1, "npc"), height = unit(1, "npc"))
  g_zoom <- grid::rasterGrob(arr_zoom, interpolate = TRUE, width = unit(1, "npc"), height = unit(1, "npc"))

  bx1 <- crop_box[1]; bx2 <- crop_box[2]
  by1 <- 1 - crop_box[4]; by2 <- 1 - crop_box[3]

  gx1 <- 0.02; gx2 <- 0.96; gy1 <- 0.05; gy2 <- 0.99
  map_x <- function(x) gx1 + x * (gx2 - gx1)
  map_y <- function(y) gy1 + y * (gy2 - gy1)

  zx1 <- 1.02; zx2 <- 1.52; zy1 <- 0.25; zy2 <- 0.75

  p <- ggplot() +
    annotation_custom(g_main, xmin = gx1, xmax = gx2, ymin = gy1, ymax = gy2) +
    annotation_custom(g_zoom, xmin = zx1, xmax = zx2, ymin = zy1, ymax = zy2) +
    annotate("rect", xmin = map_x(bx1), xmax = map_x(bx2), ymin = map_y(by1), ymax = map_y(by2),
             fill = NA, color = accent_col, linewidth = 0.85) +
    annotate("segment", x = map_x(bx2), xend = zx1, y = map_y(by2), yend = zy2,
             color = accent_col, linewidth = 0.55, linetype = "dashed") +
    annotate("segment", x = map_x(bx2), xend = zx1, y = map_y(by1), yend = zy1,
             color = accent_col, linewidth = 0.55, linetype = "dashed") +
    annotate("rect", xmin = zx1, xmax = zx2, ymin = zy1, ymax = zy2,
             fill = NA, color = accent_col, linewidth = 1.05) +
    annotate("label", x = (zx1 + zx2) / 2, y = zy2 + 0.065, label = "High-Mag Inset",
             size = 2.35, fontface = "bold", fill = accent_col, color = "white",
             linewidth = 0, label.padding = unit(0.16, "lines")) +
    coord_fixed(ratio = 1, xlim = c(0, 1.55), ylim = c(0.02, 1.01), expand = FALSE) +
    labs(title = title_lbl, caption = caption_txt) +
    theme_void(base_family = "sans") +
    theme(
      plot.title   = element_text(face = "bold", size = 9.8, color = "#0F172A", hjust = 0, margin = margin(b = 2)),
      plot.caption = element_text(size = 8.3, face = "bold", color = "#334155", hjust = 0.5, margin = margin(t = 3)),
      plot.margin  = margin(4, 6, 4, 6)
    )
  p
}
