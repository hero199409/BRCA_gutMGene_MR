# 23_0_plot_helpers.R
# Nature / Cell / Frontiers 顶刊美学标准公共绘图库
# 开放式 L 型主题、零代码味标签映射、平滑贝塞尔网络、斜体科学计数法 KM + Risk Table、高分辨率栅格与 <10 MB 导出

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  library(ggrastr)
  library(ragg)
  library(patchwork)
  library(survival)
  library(survminer)
  library(rms)
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
pub_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "03_figures/Publication_PlanA")
ihc_dir   <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "03_figures/HPA_IHC_Images")
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)

# 统一双轴语义色板：ER-/CXCL8 珊瑚红、ER+/PPARG 藏青蓝、协同绿、石板蓝、中性灰
pal_npg <- c(
  red   = "#E64B35",
  blue  = "#3C5488",
  cyan  = "#4DBBD5",
  green = "#00A087",
  coral = "#F39B7F",
  slate = "#7B8FA1",
  amber = "#D97706",
  null  = "#94A3B8"
)

# 轻盈开放式 L 型经典坐标轴主题（去四边黑框与沉重网格）
theme_journal <- function(base_size = 9.5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line         = element_line(color = "#1E293B", linewidth = 0.45),
      axis.ticks        = element_line(color = "#334155", linewidth = 0.40),
      axis.ticks.length = unit(2.5, "pt"),
      axis.title        = element_text(face = "bold", size = base_size - 0.3, color = "#0F172A"),
      axis.text         = element_text(color = "#334155", size = base_size - 1.1),
      panel.grid        = element_blank(),
      strip.background  = element_rect(fill = "#F8FAFC", color = NA),
      strip.text        = element_text(face = "bold", size = base_size - 0.3, color = "#0F172A", margin = margin(t = 3, b = 3)),
      plot.title        = element_text(face = "bold", size = base_size + 1.2, color = "#0F172A", margin = margin(b = 6)),
      plot.subtitle     = element_blank(),
      legend.title      = element_text(face = "bold", size = base_size - 0.8, color = "#0F172A"),
      legend.text       = element_text(size = base_size - 1.1, color = "#334155"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      plot.margin       = margin(6, 10, 6, 8)
    )
}

# 用于 2D 分面热图/点图/曼哈顿峰图的轻盈极简发丝框主题
theme_panel_grid <- function(base_size = 9.2, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.border      = element_rect(color = "#CBD5E1", linewidth = 0.50, fill = NA),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "#F1F5F9", linewidth = 0.35),
      axis.ticks        = element_line(color = "#64748B", linewidth = 0.35),
      axis.title        = element_text(face = "bold", size = base_size - 0.4, color = "#0F172A"),
      axis.text         = element_text(color = "#334155", size = base_size - 1.1),
      strip.background  = element_rect(fill = "#F8FAFC", color = "#CBD5E1", linewidth = 0.45),
      strip.text        = element_text(face = "bold", size = base_size - 0.4, color = "#0F172A"),
      plot.title        = element_text(face = "bold", size = base_size + 1.1, color = "#0F172A", margin = margin(b = 6)),
      plot.subtitle     = element_blank(),
      legend.title      = element_text(face = "bold", size = base_size - 0.8, color = "#0F172A"),
      legend.text       = element_text(size = base_size - 1.1, color = "#334155"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      plot.margin       = margin(6, 8, 6, 6)
    )
}

save_pub <- function(p, name, w = 15, h = 8, dpi_png = 350) {
  png_path <- file.path(pub_dir, paste0(name, ".png"))
  pdf_path <- file.path(pub_dir, paste0(name, ".pdf"))
  ggsave(png_path, p, width = w, height = h, dpi = dpi_png, device = ragg::agg_png, bg = "white")
  ggsave(pdf_path, p, width = w, height = h, device = cairo_pdf, bg = "white")
  png_mb <- file.info(png_path)$size / 1048576
  pdf_mb <- file.info(pdf_path)$size / 1048576
  cat(sprintf("Exported %-46s | PNG=%.2f MB | PDF=%.2f MB\n", name, png_mb, pdf_mb))
  if (pdf_mb >= 10) warning(sprintf("PDF >= 10 MB (%s: %.2f MB)!", name, pdf_mb))
}

# 标签去代码味工具函数
clean_subtype_str <- function(x) {
  x <- gsub("BRCA_", "", x)
  x <- gsub("ER_pos", "ER+", x)
  x <- gsub("ER_neg", "ER\u2212", x)
  x <- gsub("HER2_pos", "HER2+", x)
  x
}

clean_lineage_str <- function(x) {
  m <- c(
    "Mono/Macro"    = "Mono/Macrophages",
    "T_cells"       = "T Cells",
    "CD4Tconv"      = "CD4+ T Cells",
    "CD8T"          = "CD8+ T Cells",
    "CD8Tex"        = "Exhausted CD8+ T",
    "Tprolif"       = "Prolif. T Cells",
    "B"             = "B Cells",
    "Plasma"        = "Plasma Cells",
    "Epithelial"    = "Normal Epithelial",
    "Malignant"     = "Malignant Cells",
    "Myofibroblasts"= "Myofibroblasts",
    "Fibroblasts"   = "Fibroblasts",
    "Endothelial"   = "Endothelial",
    "PVL"           = "PVL",
    "Mast"          = "Mast Cells",
    "DC"            = "Dendritic Cells",
    "SMC"           = "Smooth Muscle"
  )
  out <- unname(m[x])
  ifelse(is.na(out), x, out)
}

# 生成标准科学计数法 plotmath 斜体 P 值表达式
format_pval_expr <- function(p) {
  if (is.na(p)) return("italic(P) == NA")
  if (p < 0.001) {
    s <- formatC(p, format = "e", digits = 1)
    parts <- strsplit(s, "e")[[1]]
    mant <- as.numeric(parts[1])
    expo <- as.integer(parts[2])
    sprintf("italic(P) == %.1f %%*%% 10^{%d}", mant, expo)
  } else {
    sprintf("italic(P) == %.3f", p)
  }
}

plot_forest_panel <- function(dt, title_txt, xlim_vals, color_map, text_x_mult = 1.03, base_size = 9.5) {
  ggplot(dt, aes(x = OR, y = label, color = sig_tag)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.55) +
    geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.28, linewidth = 0.72) +
    geom_point(aes(shape = domain_flag), size = 2.75) +
    geom_text(aes(x = ifelse(OR_UCI95 < 1.01, pmax(OR_UCI95 * text_x_mult, 1.028), OR_UCI95 * text_x_mult),
                  label = anno_txt),
              hjust = 0, size = 2.65, color = "#1E293B", show.legend = FALSE) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = c("Primary" = 16, "Warning (2 SNPs)" = 17, "LD-Confounded (Excluded)" = 4), guide = "none") +
    coord_cartesian(xlim = xlim_vals) +
    labs(title = title_txt, x = "Odds Ratio (95% CI)", y = NULL, color = NULL) +
    theme_journal(base_size) +
    theme(legend.position = "bottom")
}

# 升级版 KM + 纯 ggplot 风险表（标准斜体 P = 8.4×10⁻⁴ + 宽松风险表底标距）
plot_km_with_risk <- function(df, time_col, status_col, group_col,
                              title_txt, xlab_txt, ylab_txt,
                              legend_title, legend_labs = c("Low", "High"),
                              pal_colors, pval_coord, break_by, max_x = NULL,
                              base_size = 9.2) {
  d <- copy(df)[!is.na(get(time_col)) & get(time_col) > 0 & !is.na(get(status_col))]
  d[, .t := as.numeric(get(time_col))]
  d[, .e := as.integer(get(status_col))]
  d[, .g := factor(as.character(get(group_col)), levels = legend_labs)]
  if (is.null(max_x)) max_x <- ceiling(max(d$.t, na.rm = TRUE) / break_by) * break_by
  ticks <- seq(0, max_x, by = break_by)

  fit <- survival::survfit(survival::Surv(.t, .e) ~ .g, data = as.data.frame(d))
  sd  <- survival::survdiff(survival::Surv(.t, .e) ~ .g, data = as.data.frame(d))
  km_p <- stats::pchisq(sd$chisq, df = 1, lower.tail = FALSE)
  p_expr <- format_pval_expr(km_p)

  gs <- survminer::ggsurvplot(
    fit, data = as.data.frame(d),
    pval = FALSE, risk.table = FALSE, break.time.by = break_by, xlim = c(0, max_x),
    palette = pal_colors, legend.title = legend_title, legend.labs = legend_labs,
    xlab = xlab_txt, ylab = ylab_txt, ggtheme = theme_journal(base_size)
  )

  p_curve <- gs$plot +
    annotate("text", x = pval_coord[1], y = pval_coord[2], label = p_expr,
             parse = TRUE, hjust = 0, size = 3.2, color = "#0F172A") +
    labs(title = title_txt, y = ylab_txt, x = NULL) +
    theme(
      legend.position = "top",
      legend.margin   = margin(0, 0, -4, 0),
      plot.margin     = margin(6, 12, 2, 6),
      axis.title.x    = element_blank()
    )
  p_curve$coordinates <- coord_cartesian(xlim = c(-max_x * 0.07, max_x * 1.05), ylim = c(0, 1.03), expand = FALSE)

  sm <- summary(fit, times = ticks, extend = TRUE)
  tbl_df <- data.table(
    time   = sm$time,
    strata = factor(sub("^\\.g=", "", as.character(sm$strata)), levels = rev(legend_labs)),
    nrisk  = sm$n.risk
  )
  col_map <- setNames(pal_colors, legend_labs)

  p_tbl <- ggplot(tbl_df, aes(x = time, y = strata, label = nrisk, color = strata)) +
    geom_text(size = 2.75, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = col_map) +
    scale_x_continuous(breaks = ticks) +
    scale_y_discrete(expand = expansion(add = c(0.60, 0.55))) +
    coord_cartesian(xlim = c(-max_x * 0.07, max_x * 1.05), expand = FALSE, clip = "off") +
    labs(title = "Number at risk", x = xlab_txt, y = NULL) +
    theme_journal(base_size) +
    theme(
      plot.title       = element_text(size = base_size - 0.8, face = "bold", color = "#334155", hjust = 0, margin = margin(t = 2, b = 4)),
      axis.text.y      = element_text(size = base_size - 1.0, face = "bold", color = unname(col_map[rev(legend_labs)]), margin = margin(r = 10)),
      axis.line        = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.title.x     = element_text(size = base_size - 0.5, face = "bold", color = "#0F172A", margin = margin(t = 10, b = 2)),
      plot.margin      = margin(1, 12, 6, 6)
    )

  patchwork::wrap_elements(full = (p_curve / p_tbl + plot_layout(heights = c(4.0, 0.72))))
}

# 将圆角软填充卡片作为原生 grid grob 注入 ggplot 区域（适合高级流程图与胶囊节点）
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

# 极简无边框图像嵌入（支持底部简洁图注与白边裁剪）
img_panel <- function(img_path, top_ttl = NULL, bottom_caption = NULL,
                      trim_bg = FALSE, max_px = 1100, show_border = FALSE) {
  if (grepl("\\.jpe?g$", img_path, ignore.case = TRUE)) {
    mg <- magick::image_read(img_path)
    if (!is.null(max_px)) mg <- magick::image_scale(mg, sprintf("%dx%d", max_px, max_px))
    tmp_jpg <- tempfile(fileext = ".jpg")
    magick::image_write(mg, tmp_jpg, format = "jpeg", quality = 92)
    arr <- jpeg::readJPEG(tmp_jpg)
    unlink(tmp_jpg)
  } else {
    if (trim_bg) {
      mg <- magick::image_read(img_path)
      mg <- magick::image_trim(mg, fuzz = 5)
      mg <- magick::image_border(mg, "white", "16x12")
      tmp_trim <- tempfile(fileext = ".png")
      magick::image_write(mg, tmp_trim, format = "png")
      arr <- png::readPNG(tmp_trim)
      unlink(tmp_trim)
    } else {
      arr <- png::readPNG(img_path)
    }
  }
  asp_r <- dim(arr)[1] / dim(arr)[2]
  g <- grid::rasterGrob(arr, interpolate = TRUE, width = unit(1, "npc"), height = unit(1, "npc"))
  ggplot() +
    annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    coord_fixed(ratio = asp_r) +
    labs(title = top_ttl, caption = bottom_caption) +
    theme_void(base_family = "sans") +
    theme(
      aspect.ratio = asp_r,
      plot.title   = element_text(face = "bold", size = 10.2, hjust = 0, color = "#0F172A", margin = margin(b = 4)),
      plot.caption = element_text(size = 8.6, face = "bold", color = "#1E293B", hjust = 0.5, margin = margin(t = 5)),
      panel.border = if (show_border) element_rect(color = "#E2E8F0", fill = NA, linewidth = 0.45) else element_blank(),
      plot.margin  = margin(4, 6, 4, 6)
    )
}
