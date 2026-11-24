# CXCL8 and PPARG in breast cancer

Analysis code accompanying the **Frontiers in Immunology submission, Research Topic #88318, Final v3**, on CXCL8 in ER-negative and PPARG in ER-positive breast cancer.

## Release contents

| Directory or file | Purpose |
| --- | --- |
| `original_analysis_scripts/` | Upstream MR, colocalization, bulk and single-cell cohort analyses, therapy response, pharmacogenomics, plot caches and supplementary tables. |
| `revision_v2_scripts/` | Revised figures, panel traceability builder and I-SPY2 missing-pCR sensitivity analysis. |
| `cns_polish_scripts/` | Publication figure scripts and shared plotting theme. |
| `packaging/` | Python submission-package builder. |
| [Panel-to-data map](Panel_to_Data_Traceability.md) · [CSV](Panel_to_Data_Traceability.csv) | Panel-level source files, populations, models and generating scripts. |
| `docs/` | STROBE-MR checklist in Word and Markdown. |
| `scripts/` | Earlier release, including the Windows sparse-matrix extractor for GSE176078. |

## Configure a local workspace

Clone the repository and open a shell in its root. The commands below use the fresh clone as the analysis workspace. If your data are elsewhere, set `BRCA_PROJECT_ROOT` to the absolute directory containing `00_data/` and `02_results/`. Some scripts source helpers using the original analysis directory names. Populate these directories in a fresh workspace, then run analyses from the project root.

**PowerShell**

```powershell
$env:BRCA_PROJECT_ROOT = (Get-Location).Path
$folders = @{
    'original_analysis_scripts' = '01_scripts'
    'revision_v2_scripts' = '05_revision_v2/01_scripts'
    'cns_polish_scripts' = '06_figures_cns_polish/01_scripts'
}
foreach ($src in $folders.Keys) {
    $dest = Join-Path $env:BRCA_PROJECT_ROOT $folders[$src]
    New-Item -ItemType Directory -Force $dest | Out-Null
    Copy-Item "$src/*" $dest
}
```

**Bash**

```bash
export BRCA_PROJECT_ROOT="$PWD"
mkdir -p "$BRCA_PROJECT_ROOT"/{01_scripts,05_revision_v2/01_scripts,06_figures_cns_polish/01_scripts}
cp original_analysis_scripts/* "$BRCA_PROJECT_ROOT/01_scripts/"
cp revision_v2_scripts/* "$BRCA_PROJECT_ROOT/05_revision_v2/01_scripts/"
cp cns_polish_scripts/* "$BRCA_PROJECT_ROOT/06_figures_cns_polish/01_scripts/"
```

## Run the analysis

1. The earlier release used R 4.6.0. Start with `Rscript 01_scripts/00_install_pkgs.R` and install any additional libraries imported by the chosen module. Place public inputs under `00_data/` using the filenames expected by each script. Online MR extraction requires `MR_OFFLINE=FALSE` and a personal `OPENGWAS_JWT` in the environment. The default MR workflow reads local caches.
2. Run the numbered R scripts under `01_scripts/` in dependency order. Stages 01–10 cover gutMGene alignment, microbiome and metabolite MR, cis-eQTL MR and colocalization. Stages 11–22 cover bulk validation, TME/GSEA, single-cell analysis, treatment response, survival, pharmacogenomics and table export. When processing raw GSE176078 counts on Windows, run the earlier release's `scripts/17c_extract_gse176078_targets.py` from the project root before 17b. Scripts 23–26 build earlier figures and plot caches.
3. For the revised figures, use the panel-to-data map to identify inputs, then run the matching `02_fig*.R` in `05_revision_v2/01_scripts/` or `03_supp_figsS1toS4.R` for supplementary figures. For polished figures run the corresponding script in `06_figures_cns_polish/01_scripts/`.
4. For I-SPY2, stage 18 (`18_module2_therapy_response_cohorts.R`) creates `00_data/bulk_cohorts/planA_therapy_plot_cache.rds`. Run `Rscript 05_revision_v2/01_scripts/05_ispy2_missing_pcr_32scenarios.R` to assign all 32 binary pCR patterns to the five missing control-arm outcomes. The output is `05_revision_v2/ispy2_missing_pcr_32scenarios_sensitivity.csv`.
5. Run `python 05_revision_v2/01_scripts/04_build_traceability_table.py` to regenerate both panel-to-data files under `05_revision_v2/`.
6. For document packaging, run `python packaging/package_frontiers_submission_ready.py --help` from the clone. This step requires Pillow, PyMuPDF, XlsxWriter, Pandoc, manuscript chapters, figures, supplementary CSVs and the companion `packaging/frontiers_format/` resources from the manuscript workspace. Specify a new output directory with `--output` and use `--root` to override the project root.

## Public data access

These sources and accessions follow the submitted Data Availability Statement. Supplementary materials report screening results, denominators and access metadata.

| Source | Data or accessions |
| --- | --- |
| [gutMGene](http://bio-computing.hrbmu.edu.cn/gutmgene) | v1.0 and v2.0 curated human microbe–metabolite–host gene interactions. |
| [MRC IEU OpenGWAS](https://gwas.mrcieu.ac.uk/) | BCAC `ieu-a-1126`, `ieu-a-1127`, `ieu-a-1128`; FinnGen `finn-b-C3_BREAST`, `finn-b-C3_BREAST_EXALLC`; MiBioGen `ebi-a-GCST90016951` through `ebi-a-GCST90017039`; DMP_2022 and FR02_2022; metabolites `met-d-Acetate`, `met-c-893`, `met-a-385`, `met-a-446`, `met-a-462`; eQTLGen `eqtl-a-*`; SCALLOP `ebi-a-GCST90011994`; INTERVAL `prot-a-749`. |
| [GTEx](https://gtexportal.org/) | v8 Breast_Mammary_Tissue eQTLs. |
| [cBioPortal](https://www.cbioportal.org/) and [UCSC Xena](https://xenabrowser.net/) | Bulk cohorts, including TCGA-BRCA PanCancer Atlas 2018 and METABRIC. |
| [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) | GSE42568, GSE20685, GSE7390, GSE9195, GSE17705, GSE25066, GSE32646 and GSE173839. GSE176078 provides the Wu et al. tumor single-cell data. |
| [TISCH2](http://tisch.comp-genomics.org/) | BRCA_EMTAB8107, BRCA_GSE148673 and BRCA_GSE176078. |
| [Human Protein Atlas](https://www.proteinatlas.org/) | Pathology immunohistochemistry and single-cell expression. |
| [DepMap](https://depmap.org/) and [GDSC](https://www.cancerrxgene.org/) | DepMap 24Q4 cell-line data and GDSC2 drug sensitivity. |

Distributed under the [MIT License](LICENSE).
