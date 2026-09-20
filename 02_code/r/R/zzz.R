# Silence R CMD check notes for data.table non-standard evaluation.
utils::globalVariables(c(
  ".", ".I", ".N", ".SD",
  "alt", "ci_high", "ci_low", "cluster", "cluster_id", "count",
  "dp", "dp_minus", "dp_plus", "end", "Filter_Reason", "Filter_Status",
  "freq", "hp_run", "iend", "istart", "key", "len", "minus", "n_snv",
  "n_ins", "n_del", "ov", "passed_filter", "plus", "pos", "proportion",
  "rank", "ref", "ref_fwd", "ref_rev", "region_key", "sample_read_id",
  "Seq", "signature", "start", "support", "type",
  "cs", "DP4", "flag", "haplotype_id", "is_reference", "nm", "read_id",
  "ref_cov", "ref_end", "ref_span", "ref_start", "strand", "variants"
))

# ---------------------------------------------------------------------------
# Pairwise alignment provider.
#
# Bioconductor 3.19 moved pairwiseAlignment(), pattern(), subject(), aligned()
# and score() out of Biostrings into the pwalign package. Resolve the provider
# once at load time so the R-native alignment backend works on both old and
# new Bioconductor releases.
# ---------------------------------------------------------------------------
.pa_env <- new.env(parent = emptyenv())
.pa_optional_package <- "pwalign"

.pa_provider <- function(exported) {
  if ("pairwiseAlignment" %in% exported) {
    return(asNamespace("Biostrings"))
  }
  if (requireNamespace(.pa_optional_package, quietly = TRUE)) {
    return(asNamespace(.pa_optional_package))
  }
  stop(
    "Pairwise alignment is unavailable: this Biostrings build no longer ",
    "exports pairwiseAlignment() and the 'pwalign' package is not installed.\n",
    "Install it with: BiocManager::install(\"pwalign\")",
    call. = FALSE
  )
}

.onLoad <- function(libname, pkgname) {
  prov <- .pa_provider(getNamespaceExports("Biostrings"))
  for (fn in c("pairwiseAlignment", "pattern", "subject", "aligned", "score")) {
    assign(fn, get(fn, envir = prov), envir = .pa_env)
  }
  assign("provider", environmentName(prov), envir = .pa_env)
  invisible()
}

pa_pairwise_alignment <- function(...) .pa_env$pairwiseAlignment(...)
pa_pattern <- function(x) .pa_env$pattern(x)
pa_subject <- function(x) .pa_env$subject(x)
pa_aligned <- function(x) .pa_env$aligned(x)
pa_score <- function(x) .pa_env$score(x)
pa_provider_name <- function() .pa_env$provider
