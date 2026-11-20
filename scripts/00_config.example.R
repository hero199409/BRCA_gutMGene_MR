# Example configuration template for BRCA_gutMGene_MR
# Copy this file to `scripts/00_config.R` and configure your own OpenGWAS JWT.
# Never commit `scripts/00_config.R` containing real tokens to version control.

# 1. Optional HTTP/HTTPS proxy (leave commented out if direct connection works)
# Sys.setenv(http_proxy  = "http://127.0.0.1:YOUR_PROXY_PORT")
# Sys.setenv(https_proxy = "http://127.0.0.1:YOUR_PROXY_PORT")

# 2. MRC IEU OpenGWAS authentication JWT (obtain from https://api.opengwas.io/)
opengwas_token <- Sys.getenv("OPENGWAS_JWT", unset = "PASTE_YOUR_OPENGWAS_JWT_HERE")
Sys.setenv(OPENGWAS_JWT = opengwas_token)

# 3. Repository root and standard directory paths (relative to repo root by default)
proj_dir <- Sys.getenv("BRCA_PROJ_DIR", unset = if (dir.exists("00_data")) "." else "..")
data_dir <- file.path(proj_dir, "00_data")
res_dir  <- file.path(proj_dir, "02_results")
fig_dir  <- file.path(proj_dir, "03_figures")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(proj_dir, "tmp"), showWarnings = FALSE, recursive = TRUE)

options(timeout = 600)
message("[00_config] Environment initialized for project root: ", normalizePath(proj_dir, mustWork = FALSE))
