# Project settings. Inherit proxy settings and OpenGWAS credentials from the environment.
project_dir <- Sys.getenv("BRCA_PROJECT_ROOT", unset = Sys.getenv("BRCA_PROJECT_ROOT", unset = "."))
project_tmp <- file.path(project_dir, "tmp")
if (!dir.exists(project_tmp)) dir.create(project_tmp, recursive = TRUE)
Sys.setenv(TMPDIR = project_tmp, TEMP = project_tmp, TMP = project_tmp)
mr_offline <- tolower(Sys.getenv("MR_OFFLINE", unset = "1")) %in% c("1", "true", "yes")
# curl inherits the caller's proxy settings. Keep compatibility with existing callers.
curl_proxy_flag <- ""
