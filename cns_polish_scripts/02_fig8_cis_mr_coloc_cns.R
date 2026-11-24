# 02_fig8_cis_mr_coloc_cns.R
# True Nature Genetics / Nature Medicine Grade Figure 8:
#   Row 1: A (Top-Left, 56%) Blood cis-eQTL MR forest plot (italic genes, focal highlights, tabular OR/CI/P)
#        + B (Top-Right, 44%) Plasma IL-8 pQTL MR 3-stratum dodged forest plot (Overall / ER+ / ER- unified axis)
#   Row 2: C (Mid, Full-Width) Vertical Mirrored Colocalization (PP.H4 up in teal, PP.H3 down in amber, italic genes)
#   Row 3: D1-D4 (Bottom, 2x2) Nature Genetics-style Regional Locus & Coloc Plots (shaded delta band, significance-tinted
#          points, diamond lead SNPs, borderless inline track tags, labeled exon-intron gene tracks)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggtext)
  library(ggrastr)
  library(ragg)
  library(cowplot)
  library(patchwork)
  library(grid)
  library(scales)
})

root_dir    <- Sys.getenv("BRCA_PROJECT_ROOT", unset = ".")
cache_dir   <- file.path(root_dir, "00_data/mr_cache")
res_dir     <- file.path(root_dir, "02_results")
cns_fig_dir <- file.path(root_dir, "06_figures_cns_polish/03_figures")
dir.create(cns_fig_dir, recursive = TRUE, showWarnings = FALSE)

pal <- list(
  cxcl8      = "#E64B35",
  cxcl8_soft = "#FEF0EE",
  pparg      = "#3C5488",
  pparg_soft = "#EEF3F9",
  teal       = "#00A087",
  teal_soft  = "#E6F5F2",
  amber      = "#E08E2B",
  amber_soft = "#FEF6E7",
  slate      = "#64748B",
  muted      = "#CBD5E1",
  zebra      = "#F8FAFC",
  ink        = "#0F172A",
  subink     = "#334155"
)

cns_title <- function(letter, subtitle_md) {
  sprintf("<span style='font-size:11pt; color:#0F172A'><b>%s</b></span>\u2002<span style='font-size:8.4pt; color:#1E293B'>%s</span>",
          letter, subtitle_md)
}

theme_ng <- function(base_size = 7.8, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line           = element_line(color = pal$ink, linewidth = 0.35),
      axis.ticks          = element_line(color = pal$subink, linewidth = 0.30),
      axis.ticks.length   = unit(1.8, "pt"),
      axis.title          = element_text(face = "plain", size = base_size - 0.1, color = pal$ink),
      axis.text           = element_text(size = base_size - 0.7, color = pal$subink),
      panel.grid          = element_blank(),
      panel.background    = element_rect(fill = "white", color = NA),
      plot.background     = element_rect(fill = "white", color = NA),
      plot.title          = element_markdown(size = base_size + 0.5, color = pal$ink, hjust = 0, margin = margin(b = 4)),
      plot.title.position = "plot",
      legend.background   = element_blank(),
      legend.key          = element_blank()
    )
}

clean_sub <- function(x) {
  x <- gsub("BRCA_", "", x)
  x <- gsub("ER_pos", "ER+", x)
  x <- gsub("ER_neg", "ER\u2212", x)
  x
}

read_res <- function(name) fread(file.path(res_dir, paste0(name, ".csv")))
subtype_order <- c("Overall", "ER+", "ER\u2212")

cat("=== Generating Nature Genetics-Grade Fig 8 ===\n")

# ---------------------------------------------------------------------------
# 1. Panel A: Blood cis-eQTL MR (13 nominal from 108 BCAC tests)
# ---------------------------------------------------------------------------
cis <- read_res("module3_eqtlgen_cis_only_mr_primary_full_denominator")
stopifnot(nrow(cis) == 108L, uniqueN(cis$id.exposure) == 36L)
cis <- cis[is.finite(pval) & pval < 0.05][order(gene_symbol, outcome_label)]
cis[, `:=`(
  label_md = sprintf("*%s* (%s, %d)", gene_symbol, clean_sub(outcome_label), nsnp),
  tier     = fifelse(FDR < 0.05, "FDR < 0.05", "Nominal P < 0.05"),
  inst     = fifelse(nsnp == 1L, "1 SNP", "\u22652 SNPs"),
  stats    = sprintf("%.2f (%.2f\u2013%.2f)  P=%.2g", OR, OR_LCI95, OR_UCI95, pval),
  is_cxcl8 = (gene_symbol == "CXCL8"),
  is_pparg = (gene_symbol == "PPARG")
)]
# Use 3 decimals for OR/CI to keep 100% precision with manuscript
cis[, stats := sprintf("%.3f [%.3f, %.3f]  P=%.3g", OR, OR_LCI95, OR_UCI95, pval)]
cis[, label_f := factor(label_md, levels = rev(label_md))]
cis[, y_num   := as.numeric(label_f)]

lo_a <- min(log(cis$OR_LCI95), log(0.50))
hi_a <- max(log(cis$OR_UCI95), log(1.95))
sp_a <- hi_a - lo_a
cis[, ann_x := exp(hi_a + 0.05 * sp_a)]

p8a <- ggplot(cis, aes(x = OR, y = y_num, color = tier)) +
  # Zebra & subtle focal highlights for CXCL8 and PPARG
  geom_rect(data = cis[y_num %% 2L == 1L & !is_cxcl8 & !is_pparg],
            aes(ymin = y_num - 0.45, ymax = y_num + 0.45), xmin = -Inf, xmax = Inf,
            inherit.aes = FALSE, fill = pal$zebra, color = NA) +
  geom_rect(data = cis[is_cxcl8 == TRUE],
            aes(ymin = y_num - 0.46, ymax = y_num + 0.46), xmin = -Inf, xmax = Inf,
            inherit.aes = FALSE, fill = pal$cxcl8_soft, color = NA) +
  geom_rect(data = cis[is_pparg == TRUE],
            aes(ymin = y_num - 0.46, ymax = y_num + 0.46), xmin = -Inf, xmax = Inf,
            inherit.aes = FALSE, fill = pal$pparg_soft, color = NA) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.38) +
  geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.24, linewidth = 0.55) +
  geom_point(aes(shape = inst), size = 2.05) +
  geom_text(aes(x = ann_x, label = stats), hjust = 0, size = 1.98, color = pal$ink, show.legend = FALSE) +
  scale_color_manual(values = c("FDR < 0.05" = pal$teal, "Nominal P < 0.05" = pal$pparg), name = NULL) +
  scale_shape_manual(values = c("1 SNP" = 16, "\u22652 SNPs" = 17), name = NULL) +
  scale_x_log10(breaks = c(0.6, 1.0, 1.5), labels = c("0.6", "1.0", "1.5")) +
  scale_y_continuous(breaks = cis$y_num, labels = cis$label_md, expand = expansion(add = 0.55)) +
  coord_cartesian(xlim = exp(c(lo_a - 0.04 * sp_a, hi_a + 1.18 * sp_a))) +
  labs(
    title = cns_title("A", "Blood <i>cis</i>-eQTL MR"),
    x     = "Odds ratio (log scale)",
    y     = NULL
  ) +
  theme_ng(7.6) +
  theme(
    axis.text.y      = element_markdown(size = 6.7, color = pal$ink, margin = margin(r = 3)),
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.spacing.x = unit(3.0, "mm"),
    legend.key.width = unit(3.5, "mm"),
    legend.margin    = margin(-2, 0, 0, -6),
    legend.text      = element_text(size = 6.3, margin = margin(r = 3)),
    plot.margin      = margin(4, 4, 2, 4)
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16, size = 2.2, linewidth = 0.55), order = 1),
    shape = guide_legend(override.aes = list(color = pal$slate, size = 2.1), order = 2)
  )

# ---------------------------------------------------------------------------
# 2. Panel B: Plasma IL-8 pQTL MR (Clean 3-Stratum Dodged Forest Plot)
# ---------------------------------------------------------------------------
pq       <- read_res("planA_pqtl_il8_mr")
pq_audit <- read_res("planA_pqtl_il8_iv_audit")
selection <- data.table(
  id.exposure = c("ebi-a-GCST90011994", "prot-a-749", "ebi-a-GCST90011994", "prot-a-749", "ebi-a-GCST004445", "ebi-a-GCST90000480", "prot-c-3447_64_2"),
  study       = c("SCALLOP", "INTERVAL", "SCALLOP", "INTERVAL", "YFS/FINRISK", "LBC", "KORA"),
  iv_scope    = c(rep("Cis_1Mb_CXCL8", 2), rep("All_IVs", 5)),
  thresh      = c(5e-8, rep(1e-5, 6))
)
pq <- pq[method %chin% c("Wald ratio", "Inverse variance weighted") & outcome_label %chin% c("BRCA_Overall", "BRCA_ER_pos", "BRCA_ER_neg")]
pq <- pq[selection, on = .(id.exposure, iv_scope, thresh), nomatch = 0]
audit_single <- unique(pq_audit[is_cis_GRCh37 == TRUE, .(id.exposure, thresh = thresh_used, SNP)])
audit_single <- audit_single[selection[iv_scope == "Cis_1Mb_CXCL8"], on = .(id.exposure, thresh), nomatch = 0]

pq[, row_key := paste(study, iv_scope)]
sel_labels <- selection[, paste(study, iv_scope)]
# Get representative nsnp per row (same across outcomes or range)
nsnp_map <- pq[, .(n_lbl = if (uniqueN(nsnp) == 1L) paste0("n=", nsnp[1]) else paste0("n=", min(nsnp), "\u2013", max(nsnp))), by = row_key]
pq <- merge(pq, nsnp_map, by = "row_key", sort = FALSE)

pretty_names <- c(
  paste0("SCALLOP *cis* (", audit_single[id.exposure == "ebi-a-GCST90011994", SNP], ", n=1)"),
  paste0("INTERVAL *cis* (", audit_single[id.exposure == "prot-a-749", SNP], ", n=1)"),
  paste0("SCALLOP all IVs (", nsnp_map[row_key == "SCALLOP All_IVs", n_lbl], ")"),
  paste0("INTERVAL all IVs (", nsnp_map[row_key == "INTERVAL All_IVs", n_lbl], ")"),
  paste0("YFS/FINRISK all IVs (", nsnp_map[row_key == "YFS/FINRISK All_IVs", n_lbl], ")"),
  paste0("LBC all IVs (", nsnp_map[row_key == "LBC All_IVs", n_lbl], ")"),
  paste0("KORA all IVs (", nsnp_map[row_key == "KORA All_IVs", n_lbl], ")")
)
pq[, row_f   := factor(row_key, levels = rev(sel_labels), labels = rev(pretty_names))]
pq[, y_base  := as.numeric(row_f)]
pq[, outcome := factor(clean_sub(outcome_label), levels = subtype_order)]
pq[, y_pos   := y_base + fifelse(outcome == "ER\u2212", 0.24, fifelse(outcome == "ER+", 0.0, -0.24))]
pq[, is_sig  := pval < 0.05]

out_pal <- c("Overall" = "#64748B", "ER+" = "#3C5488", "ER\u2212" = "#E64B35")
pq[, fill_hex := ifelse(is_sig, unname(out_pal[as.character(outcome)]), "white")]

# Highlight cis rows softly
cis_bg <- data.table(ymin = 5.52, ymax = 7.48)

p8b <- ggplot(pq, aes(x = OR, y = y_pos, color = outcome)) +
  geom_rect(data = cis_bg, aes(ymin = ymin, ymax = ymax), xmin = -Inf, xmax = Inf,
            inherit.aes = FALSE, fill = "#F8FAFC", color = NA) +
  annotate("segment", x = 0.72, xend = 2.05, y = 5.5, yend = 5.5, color = "#CBD5E1", linewidth = 0.40, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#94A3B8", linewidth = 0.38) +
  geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.14, linewidth = 0.50) +
  geom_point(aes(fill = fill_hex, shape = outcome), size = 1.95, stroke = 0.55) +
  scale_fill_identity() +
  scale_color_manual(values = out_pal, name = NULL) +
  scale_shape_manual(values = c("Overall" = 21, "ER+" = 22, "ER\u2212" = 23), name = NULL) +
  scale_x_log10(breaks = c(0.8, 1.0, 1.25, 1.5, 2.0), labels = c("0.8", "1.0", "1.25", "1.5", "2.0")) +
  scale_y_continuous(breaks = seq_along(pretty_names), labels = rev(pretty_names), expand = expansion(add = 0.55)) +
  coord_cartesian(xlim = c(0.75, 2.02)) +
  labs(
    title = cns_title("B", "Plasma IL-8 pQTL MR"),
    x     = "Odds ratio (log scale)",
    y     = NULL
  ) +
  theme_ng(7.6) +
  theme(
    axis.text.y     = element_markdown(size = 6.6, color = pal$ink, margin = margin(r = 3)),
    legend.position = "bottom",
    legend.margin   = margin(-2, 0, 0, 0),
    legend.text     = element_text(size = 6.6),
    plot.margin     = margin(4, 4, 2, 6)
  ) +
  guides(
    color = guide_legend(override.aes = list(fill = unname(out_pal), size = 2.3)),
    shape = guide_legend(override.aes = list(fill = unname(out_pal), size = 2.3))
  )

# ---------------------------------------------------------------------------
# 3. Panel C: Horizontal Mirrored Colocalization (PP.H4 UP, PP.H3 DOWN)
# ---------------------------------------------------------------------------
coloc <- read_res("module3_cis_positive_genes_coloc")[p12_prior == 1e-5]
coloc[, pair_md := sprintf("*%s*<br><span style='font-size:6pt;color:#64748B;'>(%s)</span>",
                           gene_symbol, clean_sub(outcome_label))]
coloc <- coloc[order(-(PP.H4 - PP.H3))]
coloc[, x_idx := seq_len(.N)]

mirror <- rbind(
  coloc[, .(x_idx, pair_md, prob =  PP.H4, hyp = "PP.H4 (shared variant)")],
  coloc[, .(x_idx, pair_md, prob = -PP.H3, hyp = "PP.H3 (distinct variants)")]
)
mirror[, lbl := fifelse(abs(prob) >= 0.05,
                        ifelse(abs(prob) > 0.99, ">0.99", sprintf("%.2f", abs(prob))),
                        "")]

# Highlight focal pairs (CXCL8 ER-, PPARG ER+)
focal_c <- coloc[gene_symbol %chin% c("CXCL8", "PPARG"), .(x_idx, fill_bg = ifelse(gene_symbol == "CXCL8", pal$cxcl8_soft, pal$pparg_soft))]

p8c <- ggplot(mirror, aes(x = x_idx, y = prob, fill = hyp, color = hyp)) +
  geom_rect(data = focal_c, aes(xmin = x_idx - 0.44, xmax = x_idx + 0.44), ymin = -Inf, ymax = Inf,
            inherit.aes = FALSE, fill = focal_c$fill_bg, color = NA) +
  geom_hline(yintercept = 0, color = pal$ink, linewidth = 0.42) +
  geom_hline(yintercept = c(0.5, -0.5), linetype = "dotted", color = "#CBD5E1", linewidth = 0.35) +
  geom_col(width = 0.44, alpha = 0.25, color = NA) +
  geom_segment(aes(x = x_idx, xend = x_idx, y = 0, yend = prob), linewidth = 0.85) +
  geom_point(shape = 21, color = "white", stroke = 0.35, size = 2.3) +
  geom_text(data = mirror[lbl != ""],
            aes(y = prob + sign(prob) * 0.09, vjust = ifelse(prob > 0, 0, 1), label = lbl),
            size = 2.05, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values  = c("PP.H4 (shared variant)" = pal$teal, "PP.H3 (distinct variants)" = pal$amber), name = NULL) +
  scale_color_manual(values = c("PP.H4 (shared variant)" = pal$teal, "PP.H3 (distinct variants)" = pal$amber), name = NULL) +
  scale_x_continuous(breaks = coloc$x_idx, labels = coloc$pair_md, expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = c(-1.0, -0.5, 0, 0.5, 1.0), labels = c("1.0", "0.5", "0.0", "0.5", "1.0"),
                     limits = c(-1.22, 0.78)) +
  labs(
    title = cns_title("C", "<i>cis</i>-eQTL and BCAC colocalization"),
    x     = NULL,
    y     = "Posterior prob."
  ) +
  theme_ng(7.6) +
  theme(
    axis.text.x      = element_markdown(size = 6.6, lineheight = 1.02, color = pal$ink, margin = margin(t = 2)),
    legend.position  = c(0.74, 0.89),
    legend.direction = "horizontal",
    legend.text      = element_text(size = 6.6, color = pal$ink),
    legend.key.size  = unit(3.2, "mm"),
    plot.margin      = margin(4, 6, 2, 4)
  )

# ---------------------------------------------------------------------------
# 4. Panel D1-D4: Nature Genetics-Style Regional Locus & Coloc Plots (GRCh37)
# ---------------------------------------------------------------------------
annotation <- readRDS(file.path(cache_dir, "planA_cns_loci_anno_cache.rds"))
region_jobs <- data.table(
  gene    = c("CXCL8", "IL1B", "PPARG", "CXCR2"),
  outcome = c("BRCA_ER_neg", "BRCA_ER_neg", "BRCA_ER_pos", "BRCA_Overall"),
  id      = c("ieu-a-1128", "ieu-a-1128", "ieu-a-1127", "ieu-a-1126"),
  title   = c("D1", "D2", "D3", "D4")
)

plot_region_ng <- function(job) {
  row    <- coloc[gene_symbol == job$gene & outcome_label == job$outcome]
  coords <- strsplit(row$region, "[:-]")[[1L]]
  xlim   <- as.numeric(coords[2:3]) / 1e6
  input  <- readRDS(file.path(cache_dir, sprintf("coloc_%s_%s_%s_%s.rds", row$id.exposure, job$id, coords[1L], coords[2L])))
  clean_f <- function(x) as.data.table(x)[!is.na(rsid) & !is.na(beta) & !is.na(se) & se > 0][!duplicated(rsid)]
  snps   <- merge(clean_f(input$e_reg), clean_f(input$o_reg), by = "rsid", suffixes = c(".e", ".o"))[
    is.finite(eaf.e) & eaf.e > 0 & eaf.e < 1]

  eq <- snps[, .(snp = rsid, x = position.e / 1e6, logp = -log10(pmax(p.e, .Machine$double.xmin)))]
  bc <- snps[, .(snp = rsid, x = position.o / 1e6, logp = -log10(pmax(p.o, .Machine$double.xmin)))]
  pk_eq <- eq[order(-logp, snp)][1]
  pk_bc <- bc[order(-logp, snp)][1]
  gap   <- abs(pk_eq$x - pk_bc$x) * 1000
  band_xmin <- min(pk_eq$x, pk_bc$x)
  band_xmax <- max(pk_eq$x, pk_bc$x)

  # Significance-tinted color scaling per track
  eq[, rel_p := logp / max(logp, 1)]
  bc[, rel_p := logp / max(logp, 1)]

  # Helper for clean borderless track
  make_track <- function(df_tr, pk_tr, hi_col, track_lbl, show_title = FALSE, pp_txt = NULL) {
    ymax <- max(df_tr$logp, na.rm = TRUE)
    p_tr <- ggplot(df_tr, aes(x = x, y = logp)) +
      # Subtle shaded gap band between lead SNPs
      annotate("rect", xmin = band_xmin, xmax = band_xmax, ymin = -Inf, ymax = Inf, fill = "#FEF3C7", alpha = 0.45) +
      geom_vline(xintercept = pk_eq$x, color = pal$pparg, linetype = "dashed", linewidth = 0.38, alpha = 0.85) +
      geom_vline(xintercept = pk_bc$x, color = pal$cxcl8, linetype = "dashed", linewidth = 0.38, alpha = 0.85) +
      ggrastr::rasterise(geom_point(aes(color = rel_p), size = 0.55, alpha = 0.78), dpi = 300) +
      scale_color_gradient(low = "#94A3B8", high = hi_col, guide = "none") +
      geom_point(data = pk_tr, shape = 23, size = 2.2, fill = hi_col, color = "white", stroke = 0.45) +
      annotate("text", x = xlim[1] + 0.015 * diff(xlim), y = ymax * 1.07, label = track_lbl,
               hjust = 0, vjust = 1, size = 2.15, fontface = "bold", color = hi_col) +
      scale_x_continuous(limits = xlim, breaks = scales::breaks_extended(3), expand = expansion(mult = 0.02)) +
      scale_y_continuous(breaks = function(lims) { b <- scales::breaks_pretty(n = 2)(c(0, lims[2])); b[b >= 0] },
                         expand = expansion(mult = c(0.05, 0.16))) +
      labs(
        title = if (show_title) sprintf("<b>%s</b>&ensp;<i>%s</i> (%s)",
                                        job$title, job$gene, clean_sub(job$outcome)) else NULL,
        x     = NULL,
        y     = expression(-log[10](italic(P)))
      ) +
      theme_ng(7.3) +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(size = 6.2),
        axis.title.y = element_text(size = 6.6),
        plot.title   = element_markdown(size = 7.6, hjust = 0, margin = margin(b = 2)),
        plot.margin  = margin(if (show_title) 3 else 0, 5, 1, 4)
      )
    p_tr
  }

  pp_str <- sprintf("PP.H4=%.2f", row$PP.H4)
  p_eq <- make_track(eq, pk_eq, pal$pparg, "eQTLGen", show_title = TRUE, pp_txt = pp_str)
  p_bc <- make_track(bc, pk_bc, pal$cxcl8, "BCAC",    show_title = FALSE)

  # Clean exon-intron gene track with italicized gene name
  key   <- paste(job$gene, sub("BRCA_", "", job$outcome), sep = "_")
  genes <- as.data.table(annotation$genes)[pair == key & is_target == TRUE]
  exons <- as.data.table(annotation$exons)[pair == key & gene_name == job$gene & start / 1e6 <= xlim[2] & end / 1e6 >= xlim[1]]
  g_mid <- if (nrow(genes) > 0) mean(c(pmax(genes$start[1] / 1e6, xlim[1]), pmin(genes$end[1] / 1e6, xlim[2]))) else mean(xlim)
  # Place gene label slightly offset so it never clips
  lbl_x <- if (g_mid > xlim[1] + 0.72 * diff(xlim)) g_mid - 0.04 * diff(xlim) else g_mid + 0.04 * diff(xlim)
  lbl_hj <- if (g_mid > xlim[1] + 0.72 * diff(xlim)) 1 else 0

  p_ex <- ggplot() +
    annotate("rect", xmin = band_xmin, xmax = band_xmax, ymin = 0, ymax = 1, fill = "#FEF3C7", alpha = 0.45) +
    geom_segment(data = genes, aes(x = pmax(start / 1e6, xlim[1]), xend = pmin(end / 1e6, xlim[2]), y = 0.60, yend = 0.60),
                 linewidth = 0.55, color = pal$ink) +
    geom_rect(data = exons, aes(xmin = pmax(start / 1e6, xlim[1]), xmax = pmin(end / 1e6, xlim[2]), ymin = 0.36, ymax = 0.84),
              fill = pal$ink) +
    annotate("text", x = lbl_x, y = 0.60, label = job$gene, fontface = "italic",
             hjust = lbl_hj, vjust = 0.5, size = 2.15, color = pal$ink) +
    scale_x_continuous(limits = xlim, breaks = scales::breaks_extended(3), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = sprintf("Chr%s (Mb, GRCh37)", coords[1]), y = NULL) +
    theme_ng(7.3) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.text.x  = element_text(size = 6.3),
      axis.title.x = element_text(size = 6.8, margin = margin(t = 1)),
      plot.margin  = margin(0, 5, 3, 4)
    )

  cowplot::plot_grid(p_eq, p_bc, p_ex, ncol = 1, align = "v", axis = "lr", rel_heights = c(1.00, 0.86, 0.44))
}

reg_plots <- lapply(seq_len(nrow(region_jobs)), function(i) plot_region_ng(region_jobs[i]))
p8d_grid  <- cowplot::plot_grid(plotlist = reg_plots, ncol = 2, nrow = 2)

# ---------------------------------------------------------------------------
# 5. Assemble Full Nature Genetics Figure 8
# ---------------------------------------------------------------------------
row1 <- cowplot::plot_grid(p8a, p8b, ncol = 2, rel_widths = c(0.56, 0.44), align = "h", axis = "tb")
fig8_cns <- cowplot::plot_grid(
  row1,
  p8c,
  p8d_grid,
  ncol = 1,
  rel_heights = c(0.32, 0.22, 0.46)
)

png_path <- file.path(cns_fig_dir, "Fig8_CisMR_pQTL_Coloc_and_Loci_CNS.png")
pdf_path <- file.path(cns_fig_dir, "Fig8_CisMR_pQTL_Coloc_and_Loci_CNS.pdf")
ggsave(png_path, fig8_cns, width = 185 / 25.4, height = 224 / 25.4, dpi = 350, device = ragg::agg_png, bg = "white")
ggsave(pdf_path, fig8_cns, width = 185 / 25.4, height = 224 / 25.4, device = cairo_pdf, bg = "white")
cat(sprintf("Saved Fig8_CisMR_pQTL_Coloc_and_Loci_CNS | 185x224 mm | PNG=%.2f MB | PDF=%.2f MB\n",
            file.info(png_path)$size / 1048576, file.info(pdf_path)$size / 1048576))
cat("Nature Genetics Fig 8 complete.\n")
