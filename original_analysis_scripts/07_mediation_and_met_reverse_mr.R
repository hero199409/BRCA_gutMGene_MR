# 代谢物双向核查、FinnGen 独立队列检验及第一步中介 MR（MiBioGen -> 循环代谢物）。
project_root <- Sys.getenv("MR_PROJECT_ROOT", unset = Sys.getenv("BRCA_PROJECT_ROOT", unset = "."))
source(file.path(project_root, "01_scripts", "00_mr_extract_helpers.R"))
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

cache_dir <- Sys.getenv("MR_CACHE_DIR", unset = file.path(project_root, "00_data", "mr_cache"))
input_dir <- file.path(project_root, "02_results")
res_dir <- Sys.getenv("MR_OUTPUT_DIR", unset = input_dir)
bcac_ids <- names(mr_bcac_labels)
finn_ids <- c("finn-b-C3_BREAST_EXALLC", "finn-b-C3_BREAST")
methods <- c("mr_ivw", "mr_weighted_median", "mr_egger_regression", "mr_wald_ratio")
set.seed(20261122) # 复现加权中位数 bootstrap

read_primary <- function(stem) {
  name <- paste0(stem, "_bcac_mr_primary.csv")
  path <- file.path(res_dir, name)
  if (!file.exists(path)) path <- file.path(input_dir, name)
  fread(path)[id.outcome %in% bcac_ids & method %in% mr_primary_methods & is.finite(pval)]
}

# 取 BCAC 正向主检验名义候选。排除词串匹配误纳的 alpha-hydroxyisovalerate
# (met-a-616)，该 GWAS 性状与 gutMGene 的 valerate（戊酸）不同。
met_positive <- unique(read_primary("metabolite")[pval < 0.05, id.exposure])
met_ids <- setdiff(met_positive, "met-a-616")
taxa_positive <- unique(read_primary("mibiogen")[pval < 0.05, id.exposure])
if (!length(met_ids) || !length(taxa_positive)) stop("No candidates from forward BCAC MR.")

catalogue <- as.data.table(readRDS(file.path(project_root, "00_data", "opengwas_available_outcomes.rds")))
metadata <- catalogue[match(met_ids, id), .(id, trait, author, sample_size)]
if (anyNA(metadata$id) || anyNA(metadata$author)) stop("Missing metabolite mapping in the local OpenGWAS catalogue.")
print(metadata)

select_exposure <- function(exposure, ids, source) {
  if (!all(c("id.exposure", "SNP", "beta.exposure", "se.exposure", "pval.exposure") %in% names(exposure))) {
    stop("Invalid exposure cache: ", source, call. = FALSE)
  }
  exposure <- as.data.frame(exposure)
  exposure <- exposure[exposure$id.exposure %in% ids, , drop = FALSE]
  missing <- setdiff(ids, exposure$id.exposure)
  if (length(missing)) {
    stop("Missing instruments in ", source, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  exposure
}

collect_outcomes <- function(name, exposure, ids, proxies) {
  path <- file.path(cache_dir, paste0(name, ".rds"))
  snps <- exposure$SNP
  out <- mr_cache_rds(path, function() mr_fetch_outcomes(snps, ids, proxies))
  mr_select_outcomes(out, ids, snps, basename(path))
}
fit_mr <- function(harmonised, chosen_methods = methods) {
  as.data.table(TwoSampleMR::mr(harmonised, method_list = chosen_methods))
}
primary_results <- function(results) results[method %in% mr_primary_methods & is.finite(pval)]
reverse_ivw <- function(results, ids) {
  ivw <- results[method == "Inverse variance weighted" & is.finite(pval)]
  expected <- CJ(id.exposure = bcac_ids, id.outcome = ids)
  missing <- expected[!ivw, on = .(id.exposure, id.outcome)]
  if (nrow(missing)) {
    stop("Missing reverse IVW estimate for ", paste(paste(missing$id.exposure, missing$id.outcome, sep = "/"),
                                                       collapse = ", "), call. = FALSE)
  }
  ivw
}

# 反向阈值 P < 5e-8，F > 10，取得可直接复核的 BCAC -> taxa/代谢物估计。
reverse_cache <- file.path(cache_dir, "reverse_bcac_ivs.rds")
rev_exp <- mr_cache_rds(reverse_cache, function() {
  as.data.frame(rbindlist(lapply(bcac_ids, mr_fetch_iv, p_threshold = 5e-8), fill = TRUE))
})
rev_exp <- select_exposure(mr_filter_iv(rev_exp, 5e-8), bcac_ids, basename(reverse_cache))
rev_taxa_out <- collect_outcomes("reverse_mibiogen_outcome", rev_exp, taxa_positive, proxies = TRUE)
rev_taxa <- reverse_ivw(fit_mr(mr_harmonise(rev_exp, rev_taxa_out), "mr_ivw"), taxa_positive)
rev_taxa_p <- rev_taxa[, .(min_reverse_ivw_p = min(pval)), by = id.outcome]
print(rev_taxa_p)
taxa_ids <- rev_taxa_p[min_reverse_ivw_p >= 0.05, id.outcome]
if (!length(taxa_ids)) stop("No taxa pass the reverse IVW screen.")

# 反向 BCAC -> 代谢物。Shin HapMap2 按历史策略关闭代理，其余允许代理（r^2 >= 0.8）。
reverse_met_out <- collect_outcomes("reverse_met_outcome", rev_exp, met_ids,
                                    proxies = function(id) !startsWith(id, "met-a-"))
mr_rev_met <- fit_mr(mr_harmonise(rev_exp, reverse_met_out), methods[1:3])
rev_met <- reverse_ivw(mr_rev_met, met_ids)
print(rev_met[order(pval), .(id.exposure, id.outcome, nsnp, b, pval)])

# 第一阶段筛查。若将来候选中有反向名义显著的代谢物，保留其反向结果，
# 但在 taxa -> 代谢物中应用与菌群一致的反向筛选。
med_met_ids <- setdiff(met_ids, rev_met[pval < 0.05, unique(id.outcome)])
if (!length(med_met_ids)) stop("No metabolites pass the reverse IVW screen.")
taxa_cache <- file.path(cache_dir, "mibiogen_aligned_exp_ivs.rds")
taxa_exp <- mr_cache_rds(taxa_cache, function() {
  aligned <- fread(file.path(input_dir, "aligned_microbiome_gwas.csv"))
  ids <- unique(aligned[cohort == "MiBioGen_2021", id])
  as.data.frame(rbindlist(lapply(ids, mr_fetch_iv, p_threshold = 1e-5), fill = TRUE))
})
taxa_exp <- select_exposure(mr_filter_iv(taxa_exp, 1e-5), taxa_ids, basename(taxa_cache))
step1_out <- collect_outcomes("outcome_sig_met_for_taxa", taxa_exp, med_met_ids, proxies = TRUE)
step1_harm <- mr_harmonise(taxa_exp, step1_out)
step1 <- fit_mr(step1_harm)
step1_primary <- primary_results(step1)[, FDR := p.adjust(pval, method = "BH")]
expected <- CJ(id.exposure = taxa_ids, id.outcome = med_met_ids)
missing <- expected[!step1_primary, on = .(id.exposure, id.outcome)]
if (nrow(missing)) {
  message("Step-1 pairs without estimable primary MR: ",
          paste(paste(missing$id.exposure, missing$id.outcome, sep = "/"), collapse = ", "))
}
cat("Step-1 planned pairs:", nrow(expected), "| tested:", nrow(step1_primary),
    "| harmonised rows:", nrow(step1_harm), "| kept rows:", sum(step1_harm$mr_keep %in% TRUE), "\n")
print(step1_primary[pval < 0.05][order(pval), .(id.exposure, id.outcome, method, nsnp, b, pval, FDR)])

# 代谢物 -> FinnGen。聚合 RDS 优先，用原始逐性状 RDS/API 补足缺失的完整性状。
# Borges 先 P < 5e-8；LD clumping 后不足 3 枚 SNP 时回退至 P < 1e-5。
# Shin 和 Kettunen 使用 P < 1e-5。所有 SNP 均再次检查 F > 10。
met_cache <- file.path(cache_dir, "met_aligned_exp_ivs.rds")
exp_met <- if (file.exists(met_cache)) mr_cache_rds(met_cache) else data.frame()
missing_met <- setdiff(met_ids, exp_met$id.exposure)
if (length(missing_met)) {
  by_id <- lapply(missing_met, function(gwas_id) {
    path <- file.path(cache_dir, paste0("iv_met_", gwas_id, ".rds"))
    mr_cache_rds(path, function() {
      author <- metadata$author[match(gwas_id, metadata$id)]
      if (author == "Borges CM") {
        mr_fetch_iv(gwas_id, p_threshold = 5e-8, fallback_p = 1e-5)
      } else {
        mr_fetch_iv(gwas_id, p_threshold = 1e-5)
      }
    })
  })
  exp_met <- rbindlist(c(list(exp_met), by_id), fill = TRUE)
}
exp_met <- select_exposure(mr_filter_iv(exp_met, 1e-5), met_ids, basename(met_cache))
finn_met_out <- collect_outcomes("outcome_finngen_sig_met", exp_met, finn_ids, proxies = TRUE)
mr_finn_met <- fit_mr(mr_harmonise(exp_met, finn_met_out))
mr_finn_met[, OR := exp(b)]
print(primary_results(mr_finn_met)[order(pval), .(id.outcome, id.exposure, nsnp, OR, pval)])

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
# 保留原始 S6 全方法输出格式，另存每对唯一主方法与按整个矩阵校正的 BH FDR。
fwrite(step1, file.path(res_dir, "mediation_step1_taxa_to_metabolite.csv"))
fwrite(step1_primary[order(pval)], file.path(res_dir, "mediation_step1_taxa_to_metabolite_primary.csv"))
fwrite(mr_rev_met, file.path(res_dir, "metabolite_reverse_mr.csv"))
fwrite(mr_finn_met, file.path(res_dir, "metabolite_finngen_validation_mr.csv"))
