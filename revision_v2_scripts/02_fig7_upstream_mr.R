# 02_fig7_upstream_mr.R
# Figure 7 (merges Old Fig 1B-C & Fig 2C-D without repeating cis-MR):
# Upstream gut microbiota and plasma metabolite susceptibility MR (BCAC) + exact gutMGene PMID records.

local({
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)) else file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "05_revision_v2", "01_scripts")
  helper_path <- file.path(script_dir, "00_rev2_helpers.R")
  if (!file.exists(helper_path)) helper_path <- file.path("05_revision_v2", "01_scripts", "00_rev2_helpers.R")
  source(helper_path, local = FALSE)
})

cat("=== Generating Final Fig 7 (Upstream MR & gutMGene Records) ===\n")

read_res <- function(name) fread(file.path(res_dir, paste0(name, ".csv")))

nominal_mib <- read_res("mibiogen_bcac_mr_significant_p005")[gwas_name != "Marvinbryantia"]
nominal_met <- read_res("metabolite_bcac_mr_significant_p005")[
  id.exposure %chin% c("met-d-Acetate", "met-a-462", "met-a-385", "met-c-893", "met-a-446")]

prep_up <- function(d, col_lbl) {
  d <- copy(d)[order(outcome_label, pval)]
  d[, label := sprintf("%s (%s, %d SNPs)", get(col_lbl), clean_subtype_rev(outcome_label), nsnp)]
  d[, label := factor(label, levels = rev(label))]
  d[, tier := fcase(FDR < 0.05, "FDR < 0.05", FDR < 0.10, "0.05 \u2264 FDR < 0.10", default = "Nominal P < 0.05")]
  d[, inst := fifelse(nsnp == 2L, "2 SNPs", "\u22653 SNPs")]
  d[, stats := sprintf("%.2f [%.2f, %.2f] | P=%.2g", OR, OR_LCI95, OR_UCI95, pval)]
  d
}

nominal_met[, short_name := fcase(id.exposure == "met-a-385", "3-Phenylpropionate",
                                  id.exposure == "met-c-893", "Linoleic acid (18:2)",
                                  default = trait_name)]

up_cols <- c("FDR < 0.05" = pal_rev$teal, "0.05 \u2264 FDR < 0.10" = pal_rev$amber, "Nominal P < 0.05" = pal_rev$pparg)
up_shps <- c("\u22653 SNPs" = 16, "2 SNPs" = 17)

plot_up_forest <- function(d, title_txt, show_shape = TRUE) {
  lo <- min(log(d$OR_LCI95), 0); hi <- max(log(d$OR_UCI95), 0)
  sp <- hi - lo
  d[, ann_x := exp(hi + 0.06 * sp)]
  d[, y_num := as.numeric(label)]
  zb <- d[y_num %% 2L == 1L]

  ggplot(d, aes(x = OR, y = label, color = tier)) +
    geom_rect(data = zb, aes(ymin = y_num - 0.46, ymax = y_num + 0.46),
              xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = pal_rev$zebra) +
    geom_vline(xintercept = 1, linetype = "dashed", color = pal_rev$muted, linewidth = 0.45) +
    geom_errorbar(aes(xmin = OR_LCI95, xmax = OR_UCI95), orientation = "y", width = 0.24, linewidth = 0.60) +
    geom_point(aes(shape = inst), size = 2.15) +
    geom_text(aes(x = ann_x, label = stats), hjust = 0, size = 2.25, color = pal_rev$ink, show.legend = FALSE) +
    scale_color_manual(values = up_cols, name = NULL) +
    scale_shape_manual(values = up_shps, name = NULL, guide = if (show_shape) guide_legend() else "none") +
    scale_x_log10(breaks = c(0.6, 0.8, 1.0, 1.4, 1.8), labels = c("0.6", "0.8", "1.0", "1.4", "1.8")) +
    coord_cartesian(xlim = exp(c(lo - 0.06 * sp, hi + 0.92 * sp))) +
    labs(title = title_txt, x = "Odds ratio (log scale)", y = NULL) +
    theme_rev(8.2) +
    theme(legend.position = "bottom", legend.margin = margin(-3, 0, 0, 0),
          axis.text.y = element_text(size = 7.3))
}

p7a <- plot_up_forest(prep_up(nominal_mib, "gwas_name"), "A   Gut microbiota MR (MiBioGen vs. BCAC)", show_shape = TRUE)
p7b <- plot_up_forest(prep_up(nominal_met, "short_name"), "B   Plasma metabolite MR (OpenGWAS vs. BCAC)", show_shape = FALSE)

# 7C & 7D: Exact gutMGene v2.0 PMID records (no cross-table causal chaining)
gut_dir <- dirname(cache_dir)
microbe <- fread(file.path(gut_dir, "gutMGene_v2.0_Human_Gut_Microbe_Microbial_Metabolite_official_browse.csv"))
gene_tb <- fread(file.path(gut_dir, "gutMGene_v2.0_Human_Microbial_Metabolite_Host_Gene_official_browse.csv"))

make_pmid_mat <- function(d, row_key, col_key, row_names, col_names, title_txt) {
  ev <- d[get(row_key) %chin% row_names & get(col_key) %chin% col_names & !is.na(pmid),
          .(row = get(row_key), col = get(col_key), pmid)]
  ev <- ev[, .(evidence = paste(sort(unique(pmid)), collapse = "\n")), by = .(row, col)]
  ev <- ev[CJ(row = row_names, col = col_names, sorted = FALSE), on = .(row, col)]
  ev[, `:=`(row = factor(row, levels = rev(row_names)),
            col = factor(col, levels = col_names),
            found = !is.na(evidence))]
  ggplot(ev, aes(col, row)) +
    geom_tile(aes(fill = found), width = 0.94, height = 0.92, color = "white", linewidth = 0.8) +
    geom_text(aes(label = ifelse(found, evidence, "")), size = 2.35, color = "#123A3A", lineheight = 0.95) +
    scale_fill_manual(values = c(`TRUE` = "#CCFBF1", `FALSE` = "#F1F5F9"), guide = "none") +
    scale_x_discrete(position = "top", expand = expansion(add = 0.5)) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(title = title_txt, x = NULL, y = NULL) +
    theme_rev(8.2) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(face = "italic", size = 7.6),
          axis.text.x = element_text(face = "bold", size = 7.6))
}

met_names <- c("Acetate", "Butyrate", "Succinate")
p7c <- make_pmid_mat(microbe, "gut_microbiota", "metabolite", c("Dorea", "Blautia"), met_names, "C   gutMGene taxon\u2013metabolite PMIDs")
p7d <- make_pmid_mat(gene_tb, "gene", "metabolite", c("CXCL8", "IL1B", "FFAR2", "PPARG"), met_names, "D   gutMGene metabolite\u2013gene PMIDs")

fig7_all <- wrap_plots(
  wrap_elements(full = p7a),
  wrap_elements(full = p7b),
  wrap_elements(full = (p7c | p7d)),
  ncol = 1,
  heights = c(84, 68, 48)
)
save_rev_figure(fig7_all, "Fig7_Upstream_MR_and_gutMGene_Records", width_mm = 180, height_mm = 202, dpi_png = 350)
cat("Fig 7 complete.\n")
