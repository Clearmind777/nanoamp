# ---------------------------------------------------------------------------
# Command line interface
# ---------------------------------------------------------------------------

cli_usage <- function() {
  cat(
    "nanoamp - nanopore amplicon haplotype analysis\n\n",
    "Usage:\n",
    "  nanoamp call   --reads <fastq> --reference <fasta> --outdir <dir> [--mode A|B|C]\n",
    "  nanoamp batch  --sample-sheet <tsv> --outdir <dir> [--mode A|B|C]\n",
    "  nanoamp doctor\n",
    "  nanoamp help\n\n",
    "Examples:\n",
    "  nanoamp call --reads sample.fastq --reference target.fa --mode A --top-n 20 --outdir out\n",
    "  nanoamp doctor\n",
    sep = ""
  )
}

cli_call_options <- function() {
  list(
    optparse::make_option(c("--reads"), type = "character", help = "Input FASTQ"),
    optparse::make_option(c("--reference"), type = "character", help = "Target sequence FASTA"),
    optparse::make_option(c("--outdir"), type = "character", help = "Output directory"),
    optparse::make_option(c("--mode"), type = "character", default = "A",
                          help = "A reference-guided, B de novo, C exact [default A]"),
    optparse::make_option(c("--top-n"), type = "integer", default = 20,
                          help = "Number of top haplotypes to report [default 20]"),
    optparse::make_option(c("--min-reads"), type = "integer", default = 3,
                          help = "Minimum supporting reads per variant [default 3]"),
    optparse::make_option(c("--min-freq"), type = "double", default = 0.02,
                          help = "Minimum variant frequency [default 0.02]"),
    optparse::make_option(c("--min-identity"), type = "double", default = 0.90,
                          help = "Minimum read identity [default 0.90]"),
    optparse::make_option(c("--identity-cutoff"), type = "double", default = 0.99,
                          help = "Mode B clustering identity cutoff [default 0.99]"),
    optparse::make_option(c("--min-cluster-reads"), type = "integer", default = 2,
                          help = "Mode B minimum cluster size [default 2]"),
    optparse::make_option(c("--consensus-method"), type = "character", default = "decipher",
                          help = "Mode B consensus method: decipher or medoid"),
    optparse::make_option(c("--aligner"), type = "character", default = "minimap2",
                          help = "Alignment backend: minimap2 or r"),
    optparse::make_option(c("--threads"), type = "integer", default = 4,
                          help = "Number of threads [default 4]"),
    optparse::make_option(c("--ref-label"), type = "character", default = NULL,
                          help = "Reference label used in outputs"),
    optparse::make_option(c("--no-intermediates"), action = "store_true", default = FALSE,
                          help = "Do not keep BAM and other intermediate files")
  )
}

cli_batch_options <- function() {
  list(
    optparse::make_option(c("--sample-sheet"), type = "character"),
    optparse::make_option(c("--outdir"), type = "character"),
    optparse::make_option(c("--mode"), type = "character", default = "A"),
    optparse::make_option(c("--top-n"), type = "integer", default = 20),
    optparse::make_option(c("--threads"), type = "integer", default = 4),
    optparse::make_option(c("--min-reads"), type = "integer", default = 3),
    optparse::make_option(c("--min-freq"), type = "double", default = 0.02),
    optparse::make_option(c("--min-identity"), type = "double", default = 0.90),
    optparse::make_option(c("--identity-cutoff"), type = "double", default = 0.99),
    optparse::make_option(c("--min-cluster-reads"), type = "integer", default = 2),
    optparse::make_option(c("--consensus-method"), type = "character", default = "decipher"),
    optparse::make_option(c("--aligner"), type = "character", default = "minimap2"),
    optparse::make_option(c("--no-intermediates"), action = "store_true", default = FALSE)
  )
}

cli_cmd_call <- function(args) {
  opt <- optparse::parse_args(
    optparse::OptionParser(option_list = cli_call_options()), args = args
  )
  if (is.null(opt$reads) || is.null(opt$reference) || is.null(opt$outdir)) {
    cli_usage()
    stop("call requires --reads, --reference and --outdir", call. = FALSE)
  }
  run_haplotype_analysis(
    reads = opt$reads, reference = opt$reference, outdir = opt$outdir,
    mode = opt$mode, top_n = opt$`top-n`,
    min_reads = opt$`min-reads`, min_freq = opt$`min-freq`,
    min_identity = opt$`min-identity`,
    identity_cutoff = opt$`identity-cutoff`,
    min_cluster_reads = opt$`min-cluster-reads`,
    consensus_method = opt$`consensus-method`,
    aligner = opt$aligner,
    threads = opt$threads,
    keep_intermediates = !isTRUE(opt$`no-intermediates`),
    ref_label = opt$`ref-label`
  )
  invisible(TRUE)
}

cli_cmd_batch <- function(args) {
  opt <- optparse::parse_args(
    optparse::OptionParser(option_list = cli_batch_options()), args = args
  )
  if (is.null(opt$`sample-sheet`) || is.null(opt$outdir)) {
    cli_usage()
    stop("batch requires --sample-sheet and --outdir", call. = FALSE)
  }
  sheet <- data.table::fread(opt$`sample-sheet`, sep = "\t", header = TRUE)
  needed <- c("sample", "reads", "reference")
  if (!all(needed %in% names(sheet))) {
    stop(sprintf("sample-sheet must contain columns: %s", paste(needed, collapse = ", ")),
         call. = FALSE)
  }
  summary_rows <- list()
  for (i in seq_len(nrow(sheet))) {
    outdir <- file.path(opt$outdir, sheet$sample[i])
    status <- "ok"
    err <- ""
    tryCatch({
      run_haplotype_analysis(
        reads = sheet$reads[i], reference = sheet$reference[i], outdir = outdir,
        mode = opt$mode, top_n = opt$`top-n`, threads = opt$threads,
        min_reads = opt$`min-reads`, min_freq = opt$`min-freq`,
        min_identity = opt$`min-identity`, identity_cutoff = opt$`identity-cutoff`,
        min_cluster_reads = opt$`min-cluster-reads`,
        consensus_method = opt$`consensus-method`,
        aligner = opt$aligner,
        keep_intermediates = !isTRUE(opt$`no-intermediates`),
        ref_label = if ("ref_label" %in% names(sheet)) sheet$ref_label[i] else NULL
      )
      TRUE
    }, error = function(e) {
      status <<- "error"
      err <<- conditionMessage(e)
      FALSE
    })
    summary_rows[[i]] <- data.table::data.table(
      sample = sheet$sample[i], mode = opt$mode, outdir = outdir,
      status = status, error = err
    )
  }
  out <- data.table::rbindlist(summary_rows, use.names = TRUE)
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(out, file.path(opt$outdir, "batch_summary.tsv"),
                     sep = "\t", na = "")
  print(out)
  invisible(out)
}

cli_cmd_doctor <- function(args) {
  cat("nanoamp version:", nanoamp_version(), "\n")
  cat("R version:", R.version.string, "\n")
  cat("Rscript:", file.path(R.home("bin"), "Rscript"), "\n")
  cat("platform:", nanoamp_platform(), "\n")
  dep <- nanoamp_dependence_dir()
  cat("dependence directory:", if (is.null(dep)) "NOT FOUND" else dep, "\n")
  pkgs <- c("Biostrings", "IRanges", "Matrix", "Rsamtools", "ShortRead",
            "data.table", "optparse", "jsonlite", "readxl", "DECIPHER")
  for (p in pkgs) cat(sprintf("  %-12s %s\n", p, requireNamespace(p, quietly = TRUE)))
  for (tool in c("minimap2", "samtools")) {
    path <- nanoamp_tool_path(tool, required = FALSE)
    if (is.null(path)) {
      cat(sprintf("  %-12s NOT FOUND\n", tool))
    } else {
      cat(sprintf("  %-12s %s (%s)\n", tool, path, nanoamp_tool_version(tool, path)))
    }
  }
  invisible(TRUE)
}

#' nanoamp command line interface
#'
#' Run the nanoamp command line interface. This function is used by the
#' installed `nanoamp` executable and can also be called from R.
#'
#' @param args Character vector of command line arguments. Defaults to the
#'   arguments passed to `Rscript`.
#'
#' @return Invisibly returns the result of the selected command.
#' @export
nanoamp_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) >= 1 && args[1] == "--args") args <- args[-1]
  cmd <- if (length(args) >= 1 && !grepl("^--", args[1])) args[1] else "help"
  rest <- if (length(args) >= 1 && !grepl("^--", args[1])) args[-1] else args
  switch(
    cmd,
    call = cli_cmd_call(rest),
    batch = cli_cmd_batch(rest),
    doctor = cli_cmd_doctor(rest),
    help = cli_usage(),
    cli_usage()
  )
  invisible(NULL)
}
