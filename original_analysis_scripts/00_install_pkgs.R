# 00_install_pkgs.R
# 安装本项目所需的缺失 R 包（仅使用预编译二进制源，避免缺乏 Rtools 编译失败）

Sys.setenv(TMPDIR = file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "tmp"))
Sys.setenv(TEMP = file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "tmp"))
Sys.setenv(TMP = file.path(Sys.getenv("BRCA_PROJECT_ROOT", unset = "."), "tmp"))

options(
  timeout = 600,
  pkgType = "win.binary",
  install.packages.check.source = "no",
  repos = c(
    mrcieu = "https://mrcieu.r-universe.dev",
    CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
  ),
  BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
)

cran_pkgs <- c("ieugwasr", "TwoSampleMR", "MendelianRandomization", "coloc", "oncoPredict")
bioc_pkgs <- c("GEOquery")

installed <- rownames(installed.packages())

for (pkg in cran_pkgs) {
  if (!(pkg %in% installed)) {
    message("=== Installing CRAN/r-universe binary: ", pkg, " ===")
    tryCatch({
      install.packages(pkg, dependencies = TRUE, type = "win.binary")
    }, error = function(e) {
      message("Failed to install ", pkg, ": ", conditionMessage(e))
    })
  } else {
    message("Already installed: ", pkg)
  }
}

for (pkg in bioc_pkgs) {
  if (!(pkg %in% installed)) {
    message("=== Installing Bioconductor binary: ", pkg, " ===")
    tryCatch({
      BiocManager::install(pkg, ask = FALSE, update = FALSE, type = "win.binary")
    }, error = function(e) {
      message("Failed to install ", pkg, ": ", conditionMessage(e))
    })
  } else {
    message("Already installed: ", pkg)
  }
}

# MRPRESSO 是纯 R 代码包，若 CRAN 无二进制可尝试 remotes::install_github("rondolab/MR-PRESSO")
if (!("MRPRESSO" %in% rownames(installed.packages()))) {
  message("=== Installing MRPRESSO from GitHub (pure R) ===")
  tryCatch({
    remotes::install_github("rondolab/MR-PRESSO", upgrade = "never", dependencies = TRUE)
  }, error = function(e) {
    message("Failed to install MRPRESSO: ", conditionMessage(e))
  })
}

message("\n=== Final Package Check ===")
check_pkgs <- c("TwoSampleMR", "ieugwasr", "MendelianRandomization", "MRPRESSO", "coloc", "GEOquery", "oncoPredict")
final_installed <- rownames(installed.packages())
print(setNames(check_pkgs %in% final_installed, check_pkgs))
