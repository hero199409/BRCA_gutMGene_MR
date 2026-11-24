# 02_align_gutmgene_opengwas.R
# 将本地 gutMGene v1/v2.0 人类菌群与代谢物清单，与缓存的 OpenGWAS 目录精确对齐（0 API 消耗）

source(file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "01_scripts/00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

data_dir <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "00_data")
res_dir  <- file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "02_results")

# 1. 读取 5 个 gutMGene 文件并构建标准化知识库
f_mg_v1   <- fread(file.path(data_dir, "Gut_Microbe_and_Gene-human.txt"), sep = "\t")
f_mm_v1   <- fread(file.path(data_dir, "Gut_Microbe_and_Metabolite-human.txt"), sep = "\t")
f_mm_v2   <- fread(file.path(data_dir, "gutMGene_v2.0_Human_Gut_Microbe_Microbial_Metabolite_official_browse.csv"))
f_mg2_v2  <- fread(file.path(data_dir, "gutMGene_v2.0_Human_Microbial_Metabolite_Host_Gene_official_browse.csv"))
f_trip_v2 <- fread(file.path(data_dir, "gutMGene_v2.0_Human_Microbe_Metabolite_Host_Gene_PubChem_CID_join.csv"))

clean_taxon <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\[|\\]", "", x)
  x
}

extract_genus <- function(x) {
  x <- clean_taxon(x)
  # 对二名法/菌株提取首词（注意 Candidatus 情况）
  ifelse(grepl("^Candidatus ", x, ignore.case = TRUE),
         word(x, 1, 2),
         word(x, 1))
}

# 汇总所有出现过的 gutMGene 人类微生物（及其来源表 & 是否存在下游宿主基因连接）
micro_v1_g <- data.table(
  gut_microbiota = clean_taxon(f_mg_v1[["Gut Microbiota"]]),
  rank = tolower(trimws(f_mg_v1[["Classification"]])),
  has_direct_gene = TRUE,
  has_metabolite = FALSE,
  has_triplet_gene = FALSE
)
micro_v1_m <- data.table(
  gut_microbiota = clean_taxon(f_mm_v1[["Gut Microbiota"]]),
  rank = tolower(trimws(f_mm_v1[["Classification"]])),
  has_direct_gene = FALSE,
  has_metabolite = TRUE,
  has_triplet_gene = FALSE
)
micro_v2_m <- data.table(
  gut_microbiota = clean_taxon(f_mm_v2[["gut_microbiota"]]),
  rank = tolower(trimws(f_mm_v2[["rank"]])),
  has_direct_gene = FALSE,
  has_metabolite = TRUE,
  has_triplet_gene = FALSE
)
micro_trip <- data.table(
  gut_microbiota = clean_taxon(f_trip_v2[["gut_microbiota"]]),
  rank = tolower(trimws(f_trip_v2[["gut_microbiota_rank"]])),
  has_direct_gene = FALSE,
  has_metabolite = TRUE,
  has_triplet_gene = TRUE
)

all_micro <- rbindlist(list(micro_v1_g, micro_v1_m, micro_v2_m, micro_trip), use.names = TRUE)
all_micro <- all_micro[!is.na(gut_microbiota) & nzchar(gut_microbiota)]
all_micro[, taxon_lower := tolower(gut_microbiota)]
all_micro[, genus_guess := tolower(extract_genus(gut_microbiota))]

micro_summary <- all_micro[, .(
  gut_microbiota_rep = gut_microbiota[1],
  ranks = paste(unique(rank[nzchar(rank)]), collapse = ";"),
  has_direct_gene = any(has_direct_gene),
  has_metabolite = any(has_metabolite),
  has_triplet_gene = any(has_triplet_gene)
), by = .(taxon_lower, genus_guess)]

cat("gutMGene 人类去重微生物条目数:", nrow(micro_summary),
    "| 去重属名候选数:", length(unique(micro_summary$genus_guess)), "\n")

# 汇总所有有基因连接的代谢物 + 所有菌群产出的代谢物
met_with_gene <- unique(c(
  tolower(trimws(f_mg2_v2$metabolite)),
  tolower(trimws(f_trip_v2$microbe_table_metabolite)),
  tolower(trimws(f_trip_v2$gene_table_metabolite))
))
met_with_gene <- met_with_gene[!is.na(met_with_gene) & nzchar(met_with_gene)]

met_all <- unique(c(
  met_with_gene,
  tolower(trimws(f_mm_v1$Metabolite)),
  tolower(trimws(f_mm_v2$metabolite))
))
met_all <- met_all[!is.na(met_all) & nzchar(met_all)]
cat("gutMGene 人类去重代谢物总数:", length(met_all), "| 其中已有下游宿主基因记录的代谢物:", length(met_with_gene), "\n")

# 2. 读取 OpenGWAS 目录并解析菌群和代谢物
ao <- readRDS(file.path(data_dir, "opengwas_available_outcomes.rds"))

# 2A. 解析 MiBioGen (Kurilshikov A)
mibiogen <- ao[author == "Kurilshikov A"]
mibiogen[, raw_inner := str_match(trait, "\\((.*)\\)")[, 2]]
mibiogen[, gwas_rank := word(raw_inner, 1)]
mibiogen[, gwas_name := trimws(str_replace(str_replace(raw_inner, "^(phylum|class|order|family|genus)\\s+", ""), "\\s+id\\.[0-9]+$", ""))]
mibiogen[, name_norm := tolower(gsub("group$|sensu stricto [0-9]+$|[0-9]+$", "", gwas_name))]
mibiogen[, name_norm := trimws(gsub("[^a-z0-9 ]+", " ", name_norm))]
mibiogen[, name_exact := tolower(trimws(gwas_name))]

# 2B. 解析 DMP (Lopera-Maya EA)
dmp <- ao[author == "Lopera-Maya EA"]
dmp[, raw_inner := str_match(trait, "\\((.*)\\)")[, 2]]
dmp_taxa <- dmp[grepl("^k_", raw_inner)]
# 取最末一级分类名 (s_... 或 g_... 等)
dmp_taxa[, last_tok := sapply(strsplit(raw_inner, "\\."), function(v) tail(v, 1))]
dmp_taxa[, gwas_rank := fifelse(grepl("^s_", last_tok), "species",
                        fifelse(grepl("^g_", last_tok), "genus",
                        fifelse(grepl("^f_", last_tok), "family",
                        fifelse(grepl("^o_", last_tok), "order",
                        fifelse(grepl("^c_", last_tok), "class", "phylum")))))]
dmp_taxa[, gwas_name := gsub("^[kpcofgs]_", "", last_tok)]
dmp_taxa[, gwas_name := gsub("_", " ", gwas_name)]
dmp_taxa[, name_exact := tolower(trimws(gwas_name))]
dmp_taxa[, name_norm := name_exact]

# 2C. 解析 FR02 (Qin Y)
fr02 <- ao[author == "Qin Y"]
fr02[, gwas_name := trimws(str_replace(trait, "\\s+abundance in stool$", ""))]
fr02[, gwas_rank := fifelse(grepl(" ", gwas_name), "species_or_sub", "higher_taxon")]
fr02[, name_exact := tolower(gwas_name)]
fr02[, name_norm := name_exact]

# 匹配函数
match_micro <- function(dt, cohort_name) {
  dt <- copy(dt)
  dt[, cohort := cohort_name]
  # 精确匹配 taxon_lower 或 genus_guess（针对 genus 级别）
  exact_taxa <- unique(micro_summary$taxon_lower)
  genus_taxa <- unique(micro_summary$genus_guess)
  
  dt[, matched_exact := name_exact %in% exact_taxa]
  dt[, matched_genus := (name_exact %in% genus_taxa) | (name_norm %in% genus_taxa)]
  
  # 对 matched 的条目补充 gutMGene 是否有宿主基因/代谢物支撑
  m_out <- dt[matched_exact == TRUE | matched_genus == TRUE]
  m_out[, has_direct_or_triplet_gene := sapply(seq_len(.N), function(i) {
    nm <- name_exact[i]; nm2 <- name_norm[i]
    any(micro_summary[(taxon_lower == nm | genus_guess == nm | genus_guess == nm2)]$has_direct_gene |
        micro_summary[(taxon_lower == nm | genus_guess == nm | genus_guess == nm2)]$has_triplet_gene)
  })]
  m_out[, has_metabolite := sapply(seq_len(.N), function(i) {
    nm <- name_exact[i]; nm2 <- name_norm[i]
    any(micro_summary[(taxon_lower == nm | genus_guess == nm | genus_guess == nm2)]$has_metabolite)
  })]
  m_out[, .(cohort, id, trait, gwas_rank, gwas_name, sample_size, matched_exact, matched_genus, has_direct_or_triplet_gene, has_metabolite)]
}

res_mib <- match_micro(mibiogen, "MiBioGen_2021")
res_dmp <- match_micro(dmp_taxa, "DMP_2022")
res_fr02 <- match_micro(fr02, "FR02_2022")

aligned_micro <- rbindlist(list(res_mib, res_dmp, res_fr02), use.names = TRUE)
fwrite(aligned_micro, file.path(res_dir, "aligned_microbiome_gwas.csv"))

cat("\n=== 菌群 GWAS 与 gutMGene 对齐结果 ===\n")
print(aligned_micro[, .(
  total_matched = .N,
  exact_match = sum(matched_exact),
  has_host_gene = sum(has_direct_or_triplet_gene),
  has_metabolite = sum(has_metabolite)
), by = cohort])

# 3. 代谢物 GWAS 对齐 (Shin / Borges / Kettunen)
met_gwas <- ao[author %in% c("Shin", "Borges CM", "Kettunen")]
met_gwas[, trait_clean := tolower(trimws(gsub("\\*$", "", trait)))]
# 匹配：精确匹配 或 常见别名（如 butyrate/3-hydroxybutyrate, acetate, indoleacetate/indole-3-acetate 等）
met_gwas[, matched_all_gutmgene := trait_clean %in% met_all]
met_gwas[, matched_gene_gutmgene := trait_clean %in% met_with_gene]

# 关键词模糊包含检索（针对短链脂肪酸、胆汁酸、吲哚类等常见漏网别名）
key_patterns <- c("butyrate", "propionate", "acetate", "valerate", "lactate", "succinate",
                  "deoxycholate", "lithocholate", "cholate", "ursodeoxycholate",
                  "indole", "kynurenine", "serotonin", "tryptophan", "trimethylamine", "betaine",
                  "equol", "urolithin", "hippurate", "benzoate", "phenylacetate", "cinnamate")
met_gwas[, keyword_hit := grepl(paste(key_patterns, collapse = "|"), trait_clean)]

aligned_met <- met_gwas[matched_all_gutmgene == TRUE | keyword_hit == TRUE,
                        .(id, trait, author, year, sample_size, matched_all_gutmgene, matched_gene_gutmgene, keyword_hit)]
fwrite(aligned_met, file.path(res_dir, "aligned_metabolite_gwas.csv"))

cat("\n=== 代谢物 GWAS 与 gutMGene / 核心菌群代谢物对齐结果 ===\n")
cat("精确命中 gutMGene 代谢物:", sum(aligned_met$matched_all_gutmgene),
    "| 其中自带宿主基因记录:", sum(aligned_met$matched_gene_gutmgene),
    "| 关键词命中(含重叠):", nrow(aligned_met), "\n")
print(aligned_met[order(-matched_gene_gutmgene, -matched_all_gutmgene)][1:min(30, .N)])
