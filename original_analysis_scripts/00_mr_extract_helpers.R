# 本地缓存优先；仅 MR_OFFLINE=FALSE 时允许 TwoSampleMR 标准 API。
# 在线认证由 ieugwasr 从环境获取，无需读取项目配置文件。
mr_require_online <- function(operation) {
  offline <- tolower(trimws(Sys.getenv("MR_OFFLINE", unset = "true")))
  if (offline %in% c("true", "1", "yes", "on")) {
    stop("MR_OFFLINE enabled: missing local cache for ", operation, call. = FALSE)
  }
  if (!offline %in% c("false", "0", "no", "off")) {
    stop("Invalid MR_OFFLINE value; use TRUE or FALSE.", call. = FALSE)
  }
}

mr_cache_rds <- function(path, build) {
  if (file.exists(path)) {
    message("Reading MR cache: ", basename(path))
    return(readRDS(path))
  }
  # build 也可能完全通过其他本地缓存完成；API 在各提取入口拦截。
  result <- build()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  pending <- tempfile(pattern = ".mr-", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(pending), add = TRUE)
  saveRDS(result, pending)
  if (!file.exists(path) && !file.rename(pending, path) && !file.exists(path)) {
    stop("Cannot save MR cache: ", path, call. = FALSE)
  }
  result
}

mr_fetch_iv <- function(gwas_id, p_threshold, fallback_p = NULL) {
  mr_require_online(paste0("instruments of ", gwas_id))
  extract <- function(p) {
    message("Extracting instruments: ", gwas_id, " (P < ", p, ")")
    TwoSampleMR::extract_instruments(gwas_id, p1 = p, clump = TRUE, r2 = 0.001, kb = 10000)
  }
  result <- extract(p_threshold)
  if (!is.null(fallback_p) && (is.null(result) || nrow(result) < 3L)) {
    result <- extract(fallback_p)
    p_threshold <- fallback_p
  }
  if (is.null(result) || !nrow(result)) return(data.frame())
  mr_filter_iv(result, p_threshold)
}

mr_filter_iv <- function(dat, p_threshold) {
  if (!nrow(dat)) return(dat)
  dat <- as.data.frame(dat)
  f_stat <- (dat$beta.exposure / dat$se.exposure)^2
  keep <- is.finite(dat$pval.exposure) & dat$pval.exposure < p_threshold &
    is.finite(f_stat) & f_stat > 10 & is.finite(dat$se.exposure) & dat$se.exposure > 0
  dat$F_stat <- f_stat
  dat[keep, , drop = FALSE]
}

mr_fetch_outcomes <- function(snps, outcomes, proxies, chunk_size = 250L) {
  mr_require_online(paste0("outcome associations for ", paste(outcomes, collapse = ", ")))
  snps <- unique(snps[!is.na(snps) & nzchar(snps)])
  if (!length(snps) || !length(outcomes)) stop("No SNPs or outcomes requested.", call. = FALSE)
  chunks <- split(snps, ceiling(seq_along(snps) / chunk_size))
  parts <- lapply(outcomes, function(id) {
    use_proxy <- if (is.function(proxies)) proxies(id) else proxies
    message("Extracting outcome: ", id, " (proxies = ", use_proxy, ")")
    chunks_out <- lapply(chunks, function(chunk) {
      # 请求错误直接向上传播，避免把失败写成有效的部分缓存。API 正常返回的空块可以跳过。
      TwoSampleMR::extract_outcome_data(snps = chunk, outcomes = id, proxies = use_proxy,
                                        rsq = 0.8)
    })
    out <- data.table::rbindlist(chunks_out, fill = TRUE)
    if (!nrow(out)) stop("No associations returned for ", id, call. = FALSE)
    out
  })
  as.data.frame(data.table::rbindlist(parts, fill = TRUE))
}

mr_select_outcomes <- function(out, ids, snps, source) {
  if (!all(c("id.outcome", "SNP", "beta.outcome", "se.outcome") %in% names(out))) {
    stop("Invalid outcome cache: ", source, call. = FALSE)
  }
  out <- as.data.frame(out)
  out <- out[out$id.outcome %in% ids & out$SNP %in% snps, , drop = FALSE]
  missing <- setdiff(ids, out$id.outcome)
  if (length(missing)) {
    stop("No cached outcome rows matching current instruments in ", source, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  out
}

mr_harmonise <- function(exposure, outcome) {
  if (!nrow(exposure) || !nrow(outcome)) stop("No exposure or outcome rows to harmonise.")
  TwoSampleMR::harmonise_data(as.data.frame(exposure), as.data.frame(outcome), action = 2)
}

mr_bcac_labels <- c("ieu-a-1126" = "BRCA_Overall", "ieu-a-1127" = "BRCA_ER_pos",
                    "ieu-a-1128" = "BRCA_ER_neg")
mr_primary_methods <- c("Inverse variance weighted", "Wald ratio")
