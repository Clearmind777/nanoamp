# ---------------------------------------------------------------------------
# General utility functions
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

#' nanoamp version
#'
#' @return A version string.
#' @export
nanoamp_version <- function() "0.1.0"

log_msg <- function(level, ...) {
  msg <- paste0(...)
  cat(sprintf("[%s] %-5s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
  utils::flush.console()
}

log_info <- function(...) log_msg("INFO", ...)
log_warn <- function(...) log_msg("WARN", ...)
log_error <- function(...) log_msg("ERROR", ...)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

safe_div <- function(a, b) ifelse(b == 0, NA_real_, a / b)

reverse_complement <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

require_packages <- function(pkgs, strict = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    msg <- paste0(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Please install them, for example BiocManager::install(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), "))"
    )
    if (strict) stop(msg, call. = FALSE)
    log_warn(msg)
  }
  invisible(missing)
}

nanoamp_platform <- function() {
  sys <- tolower(Sys.info()[["sysname"]])
  os <- switch(sys, linux = "linux", windows = "windows", darwin = "macos", sys)
  arch <- tolower(R.version$arch)
  arch <- if (grepl("aarch64|arm64", arch)) {
    "arm64"
  } else if (grepl("x86_64|amd64", arch)) {
    "x86_64"
  } else {
    arch
  }
  paste(os, arch, sep = "-")
}

find_dependence_dir <- function(start = getwd()) {
  p <- normalizePath(start, mustWork = FALSE)
  repeat {
    for (name in c("03_dependence", "dependence")) {
      candidate <- file.path(p, name)
      if (dir.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  installed <- system.file("dependence", package = "nanoamp")
  if (nzchar(installed) && dir.exists(installed)) {
    return(normalizePath(installed, mustWork = TRUE))
  }
  NULL
}

nanoamp_dependence_dir <- function() {
  env <- Sys.getenv("NANOAMP_DEPENDENCE_DIR", unset = "")
  if (nzchar(env) && dir.exists(env)) return(normalizePath(env, mustWork = TRUE))
  find_dependence_dir()
}

nanoamp_tool_path <- function(tool, required = TRUE) {
  exe <- if (.Platform$OS.type == "windows" && !grepl("\\.exe$", tool)) {
    paste0(tool, ".exe")
  } else {
    tool
  }
  env_name <- paste0("NANOAMP_", toupper(gsub("[^A-Za-z0-9]", "_", tool)))
  env <- Sys.getenv(env_name, unset = "")
  candidates <- character(0)
  if (nzchar(env)) candidates <- c(candidates, env)
  dep <- nanoamp_dependence_dir()
  if (!is.null(dep)) {
    candidates <- c(candidates, file.path(dep, nanoamp_platform(), "bin", exe))
  }
  for (candidate in candidates) {
    if (!file.exists(candidate)) next
    prepared <- nanoamp_prepare_tool(candidate)
    if (!is.null(prepared)) return(prepared)
  }
  path <- Sys.which(tool)
  if (nzchar(path)) return(unname(path))
  if (required) {
    stop(sprintf(
      paste0(
        "External tool '%s' not found.\n",
        "Searched the %s variable, %s and PATH.\n",
        "Set %s, put the binary in <dependence>/%s/bin/, or add it to PATH."
      ),
      tool,
      env_name,
      if (is.null(dep)) "<no dependence directory>" else file.path(dep, nanoamp_platform(), "bin"),
      env_name,
      nanoamp_platform()
    ), call. = FALSE)
  }
  NULL
}

nanoamp_prepare_tool <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  if (file.access(path, 1) == 0) return(path)
  if (.Platform$OS.type != "windows") {
    try(Sys.chmod(path, mode = "0755"), silent = TRUE)
    if (file.access(path, 1) == 0) return(path)
  }
  log_warn("Bundled tool is not executable: ", path, "; falling back to PATH")
  NULL
}

nanoamp_tool_version <- function(tool, path = NULL) {
  path <- path %||% nanoamp_tool_path(tool, required = FALSE)
  if (is.null(path)) return(NA_character_)
  out <- tryCatch(
    system2(path, "--version", stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)
  )
  if (length(out) == 0) return(NA_character_)
  trimws(out[1])
}

check_external_tool <- function(tool) {
  path <- nanoamp_tool_path(tool, required = FALSE)
  if (is.null(path)) {
    nanoamp_tool_path(tool, required = TRUE)
  }
  path
}

write_tsv <- function(df, path) {
  data.table::fwrite(df, path, sep = "\t", quote = FALSE, na = "")
}

write_json <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

safe_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

format_op <- function(type, pos, ref, alt) {
  out <- character(length(pos))
  for (i in seq_along(pos)) {
    out[i] <- switch(
      type[i],
      snv = sprintf("%d%s>%s", pos[i], ref[i], alt[i]),
      ins = sprintf("%dins%s", pos[i], alt[i]),
      del = sprintf("%ddel%s", pos[i], ref[i]),
      delregion = sprintf("%ddel%s", pos[i], ref[i]),
      sprintf("%s%d%s>%s", type[i], pos[i], ref[i], alt[i])
    )
  }
  out
}

homopolymer_run <- function(ref_seq, pos) {
  L <- nchar(ref_seq)
  if (is.na(pos) || pos < 1 || pos > L) return(0L)
  base <- substr(ref_seq, pos, pos)
  if (!base %in% c("A", "C", "G", "T")) return(0L)
  i <- pos
  while (i > 1 && substr(ref_seq, i - 1, i - 1) == base) i <- i - 1
  j <- pos
  while (j < L && substr(ref_seq, j + 1, j + 1) == base) j <- j + 1
  as.integer(j - i + 1)
}

variant_homopolymer_run <- function(ref_seq, type, pos, ref) {
  L <- nchar(ref_seq)
  p <- as.integer(pos)
  if (type == "ins") {
    return(max(homopolymer_run(ref_seq, p), homopolymer_run(ref_seq, min(p + 1L, L))))
  }
  if (type %in% c("del", "delregion")) {
    span <- max(nchar(ref), 1L)
    positions <- seq.int(max(1L, p), min(L, p + span))
    return(max(vapply(positions, function(x) homopolymer_run(ref_seq, x), integer(1))))
  }
  homopolymer_run(ref_seq, p)
}

seq_context <- function(ref_seq, pos, ref, alt, width = 10) {
  L <- nchar(ref_seq)
  p <- max(1L, min(as.integer(pos), L))
  left <- substr(ref_seq, max(1L, p - width), p)
  right <- substr(ref_seq, min(L, p + 1L), min(L, p + width))
  ref_disp <- if (nzchar(ref)) ref else "-"
  alt_disp <- if (nzchar(alt)) alt else "-"
  paste0(left, "[", ref_disp, "/", alt_disp, "]", right)
}

default_params <- function() {
  list(
    top_n = 20L,
    min_reads = 3L,
    min_freq = 0.02,
    min_identity = 0.90,
    min_ref_coverage = 0.90,
    homopolymer = 4L,
    strand_bias = 0.90,
    identity_cutoff = 0.99,
    min_cluster_reads = 2L,
    max_msa_seqs = 100L,
    consensus_method = "decipher",
    aligner = "minimap2",
    threads = 4L,
    keep_intermediates = TRUE
  )
}
