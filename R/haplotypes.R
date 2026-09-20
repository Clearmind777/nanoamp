# ---------------------------------------------------------------------------
# Unified analysis dispatcher
# ---------------------------------------------------------------------------

#' Analyze nanopore amplicon haplotypes
#'
#' Align reads to a target sequence, correct sequencing errors, reconstruct
#' haplotypes, and report the top sequences with counts and proportions.
#'
#' @param reads Path to an input FASTQ file.
#' @param reference Path to a target sequence FASTA file.
#' @param outdir Output directory.
#' @param mode Analysis mode: `"A"` reference-guided correction (default),
#'   `"B"` de novo clustering, or `"C"` raw exact matching.
#' @param top_n Number of top haplotypes to report.
#' @param min_reads Minimum supporting reads for a candidate variant.
#' @param min_freq Minimum frequency for a candidate variant.
#' @param min_identity Minimum read identity to the reference.
#' @param min_ref_coverage Minimum fraction of the reference covered by a read.
#' @param homopolymer Homopolymer length threshold used for filtering.
#' @param strand_bias Strand bias filter threshold.
#' @param identity_cutoff Mode B clustering identity cutoff.
#' @param min_cluster_reads Mode B minimum cluster size.
#' @param max_msa_seqs Maximum number of sequences used for consensus alignment.
#' @param consensus_method Mode B consensus method, `"decipher"` or `"medoid"`.
#' @param aligner Alignment backend: `"minimap2"` (default) or `"r"` for the
#'   R-native pairwise alignment fallback.
#' @param use_samtools Use samtools for SAM/BAM conversion instead of
#'   Rsamtools. Rsamtools is recommended and is used by default.
#' @param threads Number of threads.
#' @param keep_intermediates Keep BAM and other intermediate files.
#' @param ref_label Optional reference label used in outputs.
#'
#' @return A list with `haplotypes`, `variants` and `qc` elements.
#' @export
run_haplotype_analysis <- function(reads, reference, outdir,
                                   mode = c("A", "B", "C"),
                                   top_n = 20L,
                                   min_reads = 3L, min_freq = 0.02,
                                   min_identity = 0.90, min_ref_coverage = 0.90,
                                   homopolymer = 4L, strand_bias = 0.90,
                                   identity_cutoff = 0.99,
                                   min_cluster_reads = 2L,
                                   max_msa_seqs = 100L,
                                   consensus_method = "decipher",
                                   aligner = c("minimap2", "r"),
                                   use_samtools = FALSE,
                                   threads = 4L,
                                   keep_intermediates = TRUE,
                                   ref_label = NULL) {
  mode <- toupper(match.arg(mode, c("A", "B", "C")))
  switch(
    mode,
    A = run_mode_a(
      reads, reference, outdir, top_n = top_n,
      min_reads = min_reads, min_freq = min_freq,
      min_identity = min_identity, min_ref_coverage = min_ref_coverage,
      homopolymer = homopolymer, strand_bias = strand_bias,
      aligner = aligner, use_samtools = use_samtools,
      threads = threads, keep_intermediates = keep_intermediates,
      ref_label = ref_label
    ),
    B = run_mode_b(
      reads, reference, outdir, top_n = top_n,
      identity_cutoff = identity_cutoff, min_cluster_reads = min_cluster_reads,
      min_identity = min_identity, min_ref_coverage = min_ref_coverage,
      max_msa_seqs = max_msa_seqs, consensus_method = consensus_method,
      aligner = aligner, use_samtools = use_samtools,
      threads = threads, keep_intermediates = keep_intermediates,
      ref_label = ref_label
    ),
    C = run_mode_c(
      reads, reference, outdir, top_n = top_n,
      keep_intermediates = keep_intermediates, ref_label = ref_label
    )
  )
}
