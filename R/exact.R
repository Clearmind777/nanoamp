# ---------------------------------------------------------------------------
# Mode C: raw exact matching (diagnostic mode)
# ---------------------------------------------------------------------------

run_mode_c <- function(reads_path, reference_path, outdir,
                       top_n = 20L, keep_intermediates = TRUE,
                       ref_label = NULL) {
  outdir <- ensure_dir(outdir)
  ref <- read_reference(reference_path)
  fq <- read_fastq(reads_path)
  log_info("Mode C: ", basename(reads_path), " -> ", ref$name, " (", nrow(fq), " reads)")
  if (nrow(fq) == 0) stop("Mode C: FASTQ is empty", call. = FALSE)

  ref_rc <- reverse_complement(ref$sequence)
  exact_fwd <- fq$sequence == ref$sequence
  exact_rev <- fq$sequence == ref_rc
  exact_any <- exact_fwd | exact_rev
  contains_ref <- grepl(ref$sequence, fq$sequence, fixed = TRUE) |
    grepl(ref_rc, fq$sequence, fixed = TRUE)

  counts <- fq[, .(count = .N), by = sequence][order(-count, sequence)]
  counts[, proportion := count / sum(count)]
  counts[, is_reference := sequence == ref$sequence | sequence == ref_rc]
  counts[, rank := seq_len(.N)]
  ci <- t(vapply(seq_len(nrow(counts)), function(i) {
    wilson_ci(counts$count[i], sum(counts$count))
  }, numeric(2)))
  counts[, `:=`(ci_low = ci[, 1], ci_high = ci[, 2], length = nchar(sequence))]

  write_tsv(counts[, .(rank, haplotype_id = sprintf("H%d", rank), count, proportion,
                       ci_low, ci_high, is_reference, length, sequence)],
            file.path(outdir, "haplotypes.tsv"))
  top <- utils::head(counts, top_n)
  write_fasta(stats::setNames(top$sequence, sprintf("H%d", top$rank)),
              file.path(outdir, "haplotypes.fasta"))
  write_tsv(data.table::data.table(
    Chr = character(0), Pos = integer(0), Ref = character(0), Alt = character(0),
    DP = integer(0), Ref_dp = integer(0), Alt_dp = integer(0), Freq = numeric(0),
    DP4 = character(0), Seq = character(0),
    Filter_Status = character(0), Filter_Reason = character(0)
  ), file.path(outdir, "variants.tsv"))

  qc <- list(
    mode = "C",
    reference_label = ref_label %||% ref$name,
    reference_length = ref$length,
    n_reads_total = nrow(fq),
    n_exact_forward = sum(exact_fwd),
    n_exact_reverse = sum(exact_rev),
    n_exact_either = sum(exact_any),
    exact_any_proportion = round(mean(exact_any), 6),
    n_contains_reference = sum(contains_ref),
    contains_reference_proportion = round(mean(contains_ref), 6),
    n_unique_raw_sequences = nrow(counts),
    top1_proportion = round(counts$proportion[1], 6),
    top1_is_reference = counts$is_reference[1]
  )
  write_tsv(build_qc_table(qc), file.path(outdir, "qc.tsv"))
  run_manifest(outdir, "C", list(top_n = top_n), ref, qc,
               extra = list(reads_md5 = safe_md5(reads_path)))
  invisible(list(haplotypes = counts, variants = NULL, qc = qc))
}
