# 01_test_opengwas_catalog.R
# 测试 OpenGWAS JWT 鉴权并检索乳腺癌 / 肠道菌群 / 代谢物 / eQTL 可用 GWAS 数据集

source(if (file.exists("scripts/00_config.R")) "scripts/00_config.R" else "00_config.R")
suppressPackageStartupMessages({
  library(ieugwasr)
  library(TwoSampleMR)
  library(data.table)
})

cat("=== 1. Testing User Auth & Single GWAS Info ===\n")
u <- tryCatch(ieugwasr::user(), error = function(e) e)
print(u)

g_test <- ieugwasr::gwasinfo("ieu-a-1126")
print(g_test[, c("id", "trait", "sample_size", "ncase", "ncontrol", "population", "year", "author")])

cat("\n=== 2. Pulling Available Outcomes Metadata (caching locally) ===\n")
ao <- TwoSampleMR::available_outcomes()
ao_dt <- as.data.table(ao)
saveRDS(ao_dt, "00_data/opengwas_available_outcomes.rds")
cat("Total available datasets:", nrow(ao_dt), "\n")

cat("\n--- Breast Cancer GWAS (BCAC / FinnGen / UKB) ---\n")
brca <- ao_dt[grepl("breast cancer|malignant neoplasm of breast", trait, ignore.case = TRUE)]
brca <- brca[order(-sample_size)]
print(brca[1:min(25, .N), .(id, trait, population, sample_size, ncase, ncontrol, author, year)])

cat("\n--- Gut Microbiome GWAS ---\n")
micro <- ao_dt[grepl("gut microbiota|microbiome|abundance|genus\\.|family\\.|order\\.|class\\.|phylum\\.", trait, ignore.case = TRUE) | grepl("ebi-a-GCST9001", id)]
cat("Matched microbiome-like entries:", nrow(micro), "\n")
print(micro[1:min(15, .N), .(id, trait, sample_size, author, year)])

cat("\n--- Metabolite GWAS counts by prefix ---\n")
met <- ao_dt[grepl("^met-", id) | grepl("metabolite|butyrate|propionate|bile acid|acetate|indole", trait, ignore.case = TRUE)]
cat("Matched metabolite entries:", nrow(met), "\n")
print(table(sub("-.*", "", met$id)))

cat("\n--- eQTLGen entries count ---\n")
eqtl <- ao_dt[grepl("^eqtl-a-", id)]
cat("eQTLGen entries:", nrow(eqtl), "\n")
