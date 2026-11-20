# Subtype-Divergent Gut Microbiota–Metabolite–Host Gene Axes in Breast Cancer (*CXCL8* in ER− and *PPARG* in ER+)

Reproducible analysis code accompanying the manuscript:

> **Multi-Omics Mendelian Randomization and Cross-Cohort Single-Cell Validation Reveal Subtype-Divergent Gut Microbiota–Metabolite–Host Gene Axes (*CXCL8* in ER-Negative and *PPARG* in ER-Positive Breast Cancer) Governing Tumor Microenvironment Remodeling and Immunochemotherapy Response**

---

## 1. Repository Structure

```text
.
├── .gitignore
├── LICENSE
├── README.md
└── scripts/
    ├── 00_config.example.R                           # Template config (copy to scripts/00_config.R)
    ├── 00_install_pkgs.R                             # R / Bioconductor package installation
    ├── 01_check_gutmgene.R                           # gutMGene v1.0/v2.0 parsing
    ├── 01_test_opengwas_catalog.R                    # MRC IEU OpenGWAS catalog snapshot
    ├── 02_align_gutmgene_opengwas.R                  # Taxa & metabolite alignment to OpenGWAS
    ├── 03_run_mibiogen_brca_mr.R                     # Stage 1: MiBioGen -> BCAC two-sample MR
    ├── 04_check_pos_taxa_gutmgene.R                  # Curated microbe-gene linkage audit
    ├── 05_mibiogen_finngen_and_reverse_mr.R          # Layer 2 FinnGen meta & bidirectional reverse MR
    ├── 06_run_metabolite_brca_mr.R                   # Stage 2: 47 plasma metabolites -> BCAC MR
    ├── 07_mediation_and_met_reverse_mr.R             # Step-1 mediation MR & metabolite reverse MR
    ├── 08_build_module3_candidate_genes.R            # 61 gutMGene candidate host genes extraction
    ├── 09_run_eqtlgen_all49_mr.R                     # Initial eQTLGen MR screening
    ├── 10_verify_pos_eqtl_cis_coloc_gtex.R           # Stage 3: +-1 Mb cis-audit, coloc.abf, GTEx v8
    ├── 10b_check_gtex_breast_single_eqtl.R           # GTEx v8 Breast_Mammary_Tissue single-variant check
    ├── 10c_gtex_breast_dyneqtl.R                     # GTEx v8 dynamic eQTL lookup
    ├── 11_module4_bulk_tcga_metabric_geo.R           # Baseline TCGA, METABRIC, GSE42568, GSE20685, GSE7390
    ├── 11b_fix_geo_and_cxcl8_audit.R                 # Entrez 3576 CXCL8 & GEO log2 auditing
    ├── 11c_update_module4_complete.R                 # 5-cohort survival meta-analysis
    ├── 12_module5_tme_gsea_local.R                   # immunedeconv (quanTIseq/MCP/ESTIMATE) & Hallmark GSEA
    ├── 13_module6_gdsc2_drugs_local.R                # DepMap 24Q4 x GDSC2 subtype-adjusted drug screen
    ├── 14_generate_draft_figures.R                   # Exploratory diagnostic graphics
    ├── 15_planA_deep_validation_and_hpa.R            # Multi-method sensitivity, HPA XML & IHC audit
    ├── 16_generate_publication_figures_planA.R       # Early Plan A publication figures
    ├── 17a_adipose_confounding_and_multivar_erpos.R  # ADIPOQ/PLIN1/StromalScore Cox deconfounding
    ├── 17b_tisch2_tumor_scrna_analysis.R             # Tumor scRNA-seq across GSE176078, EMTAB8107, GSE148673
    ├── 17c_extract_gse176078_targets.py              # Sparse MatrixMarket streaming for 100,064 cells
    ├── 18_module2_therapy_response_cohorts.R         # Endocrine, NACT & I-SPY2 GSE173839 therapy cohorts
    ├── 18b_module2_interactions_and_multivar.R       # Formal Biomarker x Arm interaction & DRFS models
    ├── 19_module3_dmp_fr02_and_pqtl.R                # DMP_2022/FR02_2022 microbiome & plasma IL-8 pQTL MR
    ├── 20_module4_os_meta_and_nomogram.R             # CXCL8 OS heterogeneity decomposition & ER+ nomogram
    ├── 21_generate_final_planA_figures.R             # Pre-revision composite Plan A figures
    ├── 22_export_supplementary_tables_upgrade.R      # Export 36 publication Supplementary Tables S1-S28b
    ├── 23_0_plot_helpers.R                           # Shared Frontiers plotting helpers & theme
    ├── 23_generate_frontiers_part1_figs1to5.R        # Final Figures 1-5 (300 dpi PNG & PDF)
    └── 24_generate_frontiers_part2_figs6toS3.R       # Final Figures 6-8 & Figures S1-S3 (300 dpi PNG & PDF)
```

---

## 2. Software Environment & Key Dependencies

Analyses were developed and verified under:
- **R `4.6.0`** with key packages: `TwoSampleMR` (`0.7.14`), `ieugwasr` (`1.2.1`), `MRPRESSO` (`1.0`), `coloc` (`5.2.3`), `meta` (`9.0.0`), `survival` (`3.8.6`), `rms` (`8.1.1`), `timeROC` (`0.4.1`), `immunedeconv` (`2.1.4`), `clusterProfiler` (`4.20.0`), `limma` (`3.67.3`), `GEOquery` (`2.80.1`), `rhdf5` (`2.56.0`), `org.Hs.eg.db` (`3.23.1`), `data.table` (`1.18.2.1`), `ggplot2` (`4.0.3`), `patchwork` (`1.3.2`), `ggrastr` (`1.0.2`), and `ragg` (`1.5.2`).
- **Python 3** (standard library for sparse MatrixMarket extraction in `scripts/17c_extract_gse176078_targets.py`).

---

## 3. Setup & Data Sources

```bash
git clone https://github.com/hero199409/BRCA_gutMGene_MR.git
cd BRCA_gutMGene_MR
```

1. **Configure Environment:**
   ```bash
   cp scripts/00_config.example.R scripts/00_config.R
   ```
   Obtain a personal JSON Web Token (JWT) from the MRC IEU OpenGWAS portal (`https://api.opengwas.io/`) and either set `export OPENGWAS_JWT="your_jwt_here"` or paste it into `scripts/00_config.R` (ignored by `.gitignore`).

2. **Prepare Public Input Data (`00_data/`):**
   - **gutMGene v1.0 / v2.0:** Place curated human microbe–metabolite–host gene files in `00_data/` (`http://bio-computing.hrbmu.edu.cn/gutmgene`).
   - **GWAS / eQTL / pQTL Summary Data:** Queried programmatically via `ieugwasr` / `TwoSampleMR` against `BCAC` (`ieu-a-1126/1127/1128`), `FinnGen` (`finn-b-C3_BREAST*`), `MiBioGen_2021`, `DMP_2022`, `FR02_2022`, plasma metabolites (`met-*`), `eQTLGen` (`eqtl-a-*`), and `SCALLOP`/`INTERVAL`/`YFS`/`KORA`/`LBC` IL-8 pQTLs.
   - **Tissue eQTL & Pathology:** `GTEx Portal v8` REST API (`Breast_Mammary_Tissue`) and `Human Protein Atlas` XMLs (`ENSG00000169429`, `ENSG00000132170`, `ENSG00000125538`, `ENSG00000126262`).
   - **Bulk, Therapy & Single-Cell Cohorts:** `TCGA-BRCA` (PanCancer Atlas 2018), `METABRIC`, `NCBI GEO` (`GSE42568`, `GSE20685`, `GSE7390`, `GSE9195`, `GSE17705`, `GSE25066`, `GSE32646`, `GSE173839`), `TISCH2` HDF5 (`BRCA_EMTAB8107`, `BRCA_GSE148673`, `BRCA_GSE176078`) plus raw `GSE176078` sparse counts, and `DepMap 24Q4` × `GDSC2`.

---

## 4. Execution Workflow

Run all commands from the repository root (`BRCA_gutMGene_MR/`) in sequential stage order:

```bash
# 0. Install required CRAN/Bioconductor/GitHub packages
Rscript scripts/00_install_pkgs.R

# 1. Stage 1-3: gutMGene curation, OpenGWAS alignment, Microbiome/Metabolite/cis-eQTL MR & Coloc
Rscript scripts/01_check_gutmgene.R
Rscript scripts/01_test_opengwas_catalog.R
Rscript scripts/02_align_gutmgene_opengwas.R
Rscript scripts/03_run_mibiogen_brca_mr.R
Rscript scripts/04_check_pos_taxa_gutmgene.R
Rscript scripts/05_mibiogen_finngen_and_reverse_mr.R
Rscript scripts/06_run_metabolite_brca_mr.R
Rscript scripts/07_mediation_and_met_reverse_mr.R
Rscript scripts/08_build_module3_candidate_genes.R
Rscript scripts/09_run_eqtlgen_all49_mr.R
Rscript scripts/10_verify_pos_eqtl_cis_coloc_gtex.R
Rscript scripts/10b_check_gtex_breast_single_eqtl.R
Rscript scripts/10c_gtex_breast_dyneqtl.R

# 2. Layer 3-5 External Validation (Bulk Survival, TME/GSEA, GDSC2, HPA, Adipose, scRNA-seq, Therapy)
Rscript scripts/11_module4_bulk_tcga_metabric_geo.R
Rscript scripts/11b_fix_geo_and_cxcl8_audit.R
Rscript scripts/11c_update_module4_complete.R
Rscript scripts/12_module5_tme_gsea_local.R
Rscript scripts/13_module6_gdsc2_drugs_local.R
Rscript scripts/15_planA_deep_validation_and_hpa.R
Rscript scripts/17a_adipose_confounding_and_multivar_erpos.R
python  scripts/17c_extract_gse176078_targets.py
Rscript scripts/17b_tisch2_tumor_scrna_analysis.R
Rscript scripts/18_module2_therapy_response_cohorts.R
Rscript scripts/18b_module2_interactions_and_multivar.R
Rscript scripts/19_module3_dmp_fr02_and_pqtl.R
Rscript scripts/20_module4_os_meta_and_nomogram.R

# 3. Export All 36 Publication Supplementary Tables (S1-S28b) & Final Figures (Fig 1-8, Fig S1-S3)
Rscript scripts/22_export_supplementary_tables_upgrade.R
Rscript scripts/23_generate_frontiers_part1_figs1to5.R
Rscript scripts/24_generate_frontiers_part2_figs6toS3.R
```

---

## 5. License

Distributed under the MIT License. See [`LICENSE`](LICENSE) for details.
