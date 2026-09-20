# ---------------------------------------------------------------------------
# Mode B: de novo clustering and cluster consensus
# ---------------------------------------------------------------------------

greedy_cluster_from_distance <- function(dm, cutoff) {
  n <- nrow(dm)
  assigned <- rep(NA_integer_, n)
  k <- 0L
  for (i in seq_len(n)) {
    if (!is.na(assigned[i])) next
    k <- k + 1L
    assigned[i] <- k
    close <- which(is.na(assigned) & dm[i, ] <= cutoff)
    assigned[close] <- k
  }
  assigned
}

cluster_variant_patterns <- function(read_ids, read_vars, ref_len, identity_cutoff) {
  n <- length(read_ids)
  if (is.null(ref_len) || is.na(ref_len) || ref_len <= 0) ref_len <- 1L
  if (n == 1) {
    return(list(cluster = 1L, distance = matrix(0, 1, 1), method = "variant_matrix"))
  }
  rv <- data.table::as.data.table(read_vars)
  if (nrow(rv) == 0) {
    return(list(
      cluster = rep(1L, n),
      distance = matrix(0, n, n),
      method = "variant_matrix"
    ))
  }
  rv <- unique(rv[, .(read_id, key)])
  mi <- match(rv$read_id, read_ids)
  keep <- !is.na(mi)
  rv <- rv[keep]
  mi <- mi[keep]
  keys <- unique(rv$key)
  mj <- match(rv$key, keys)
  M <- Matrix::sparseMatrix(i = mi, j = mj, x = 1, dims = c(n, length(keys)))
  s <- Matrix::rowSums(M)
  shared <- as.matrix(M %*% Matrix::t(M))
  d <- outer(s, s, "+") - 2 * shared
  d <- d / ref_len
  clv <- greedy_cluster_from_distance(d, 1 - identity_cutoff)
  names(clv) <- read_ids
  list(cluster = clv, distance = d, method = "variant_matrix")
}

cluster_sequences <- function(seqs, identity_cutoff = 0.99, threads = 4L,
                              read_ids = NULL, read_vars = NULL, ref_len = NULL,
                              min_coverage = 0.5) {
  n <- length(seqs)
  if (n == 1) {
    return(list(cluster = 1L, distance = matrix(0, 1, 1), method = "single"))
  }
  cutoff <- 1 - identity_cutoff
  if (requireNamespace("DECIPHER", quietly = TRUE)) {
    ans <- tryCatch({
      x <- Biostrings::DNAStringSet(toupper(seqs))
      names(x) <- sprintf("r%06d", seq_len(n))
      d <- DECIPHER::DistanceMatrix(x, processors = threads, verbose = FALSE)
      dm <- as.matrix(d)
      exports <- getNamespaceExports("DECIPHER")
      if ("Clusterize" %in% exports) {
        cl <- DECIPHER::Clusterize(
          x, cutoff = cutoff, minCoverage = min_coverage,
          processors = threads, verbose = FALSE
        )
        clv <- as.integer(cl$cluster)
        names(clv) <- rownames(cl)
        clv <- clv[names(x)]
        method <- "DECIPHER::Clusterize"
      } else if ("IdClusters" %in% exports) {
        id_clusters <- getExportedValue("DECIPHER", "IdClusters")
        cl <- id_clusters(
          x, myDistMatrix = d, cutoff = cutoff,
          method = "complete", showPlot = FALSE, verbose = FALSE
        )
        clv <- if (is.null(dim(cl))) as.integer(cl) else as.integer(cl[, 1])
        if (!is.null(rownames(cl))) names(clv) <- rownames(cl)
        clv <- clv[names(x)]
        method <- "DECIPHER::IdClusters"
      } else {
        hc <- stats::hclust(stats::as.dist(dm), method = "complete")
        clv <- stats::cutree(hc, h = cutoff)
        names(clv) <- rownames(dm)
        clv <- clv[names(x)]
        method <- "DECIPHER::DistanceMatrix+hclust"
      }
      list(cluster = unname(clv), distance = dm, method = method)
    }, error = function(e) {
      log_warn("DECIPHER clustering failed; falling back to greedy clustering: ",
               conditionMessage(e))
      NULL
    })
    if (!is.null(ans)) return(ans)
  }
  if (!is.null(read_ids) && !is.null(read_vars)) {
    log_warn("DECIPHER unavailable; Mode B uses variant-pattern greedy clustering")
    return(cluster_variant_patterns(read_ids, read_vars, ref_len, identity_cutoff))
  }
  stop("DECIPHER unavailable and read_vars not supplied; cannot cluster for Mode B",
       call. = FALSE)
}

medoid_index <- function(dm, idx) {
  if (length(idx) == 1) return(idx)
  if (is.null(dm)) return(idx[1])
  sub <- dm[idx, idx, drop = FALSE]
  idx[which.min(rowSums(sub))]
}

majority_consensus <- function(aln) {
  m <- as.matrix(aln)
  n <- ncol(m)
  out <- character(n)
  for (j in seq_len(n)) {
    col <- m[, j]
    col <- col[col %in% c("A", "C", "G", "T", "-")]
    if (length(col) == 0) {
      out[j] <- ""
      next
    }
    tab <- table(col)
    top <- names(tab)[which.max(tab)]
    out[j] <- if (top == "-") "" else top
  }
  paste0(out, collapse = "")
}

build_cluster_consensus <- function(seqs, idx, dm,
                                    method = "decipher", max_msa_seqs = 100L, threads = 4L) {
  if (length(idx) == 1) return(seqs[idx])
  if (method %in% c("decipher", "auto") && requireNamespace("DECIPHER", quietly = TRUE)) {
    sub_idx <- idx
    if (length(idx) > max_msa_seqs) {
      med <- medoid_index(dm, idx)
      pos <- unique(round(seq(1, length(idx), length.out = max_msa_seqs)))
      sub_idx <- unique(c(med, idx[pos]))
    }
    ans <- tryCatch({
      x <- Biostrings::DNAStringSet(toupper(seqs[sub_idx]))
      names(x) <- sprintf("s%05d", seq_along(sub_idx))
      aln <- DECIPHER::AlignSeqs(
        x,
        processors = threads, verbose = FALSE
      )
      majority_consensus(aln)
    }, error = function(e) {
      log_warn("DECIPHER consensus failed; using medoid: ", conditionMessage(e))
      NULL
    })
    if (!is.null(ans) && nzchar(ans)) return(ans)
  }
  seqs[medoid_index(dm, idx)]
}

pairwise_diffs <- function(query, ref) {
  pa <- pa_pairwise_alignment(
    Biostrings::DNAString(query),
    Biostrings::DNAString(ref),
    type = "global", gapOpening = 5, gapExtension = 1
  )
  q <- as.character(pa_aligned(pa_pattern(pa)))
  s <- as.character(pa_aligned(pa_subject(pa)))
  alignment_to_ops(q, s)
}

run_mode_b <- function(reads_path, reference_path, outdir,
                       top_n = 20L,
                       identity_cutoff = 0.99, min_cluster_reads = 2L,
                       min_identity = 0.90, min_ref_coverage = 0.90,
                       max_msa_seqs = 100L, consensus_method = "decipher",
                       aligner = c("minimap2", "r"), use_samtools = FALSE,
                       threads = 4L, keep_intermediates = TRUE,
                       ref_label = NULL) {
  aligner <- match.arg(aligner, c("minimap2", "r"))
  outdir <- ensure_dir(outdir)
  ref <- read_reference(reference_path)
  n_total <- count_fastq_reads(reads_path)
  log_info("Mode B: ", basename(reads_path), " -> ", ref$name, " (", n_total, " reads)")

  prep <- prepare_alignment_data(
    reads_path, reference_path, outdir,
    aligner = aligner, threads = threads, use_samtools = use_samtools
  )
  aln <- prep$aln
  n_primary <- nrow(aln)
  if (n_primary == 0) stop("Mode B: no aligned reads", call. = FALSE)
  kept <- filter_alignment_reads(aln, min_identity, min_ref_coverage)
  if (nrow(kept) == 0) stop("Mode B: no reads left after filtering", call. = FALSE)

  read_vars <- prep$read_vars
  if (nrow(read_vars) > 0) read_vars <- read_vars[read_id %in% kept$read_id]
  seqs <- vapply(seq_len(nrow(kept)), function(i) {
    ops <- if (nrow(read_vars) > 0) read_vars[read_id == kept$read_id[i]] else empty_ops()
    apply_variants(ref$sequence, ops)
  }, character(1))

  cl <- cluster_sequences(
    seqs, identity_cutoff, threads,
    read_ids = kept$read_id, read_vars = read_vars, ref_len = ref$length
  )
  labels <- cl$cluster
  dm <- cl$distance
  cluster_ids <- sort(unique(labels))
  rows <- lapply(cluster_ids, function(cid) {
    idx <- which(labels == cid)
    consensus <- build_cluster_consensus(
      seqs, idx, dm, method = consensus_method,
      max_msa_seqs = max_msa_seqs, threads = threads
    )
    annotate <- length(idx) >= min_cluster_reads
    diffs <- if (annotate) pairwise_diffs(consensus, ref$sequence) else empty_ops()
    data.table::data.table(
      cluster_id = cid,
      count = length(idx),
      consensus = consensus,
      length = nchar(consensus),
      is_reference = if (annotate) identical(consensus, ref$sequence) else NA,
      n_snv = if (annotate) sum(diffs$type == "snv") else NA_integer_,
      n_ins = if (annotate) sum(diffs$type == "ins") else NA_integer_,
      n_del = if (annotate) sum(diffs$type == "del") else NA_integer_,
      variants = if (!annotate) "." else if (nrow(diffs) == 0) "." else paste(
        format_op(diffs$type, diffs$pos, diffs$ref, diffs$alt), collapse = ";"
      ),
      annotated = annotate,
      diff_ops = list(diffs)
    )
  })
  clusters <- data.table::rbindlist(rows, use.names = TRUE)
  clusters <- clusters[order(-count, cluster_id)]
  clusters[, passed_filter := count >= min_cluster_reads]
  clusters[, proportion := count / sum(count)]
  clusters[, rank := seq_len(.N)]
  ci <- t(vapply(seq_len(nrow(clusters)), function(i) {
    wilson_ci(clusters$count[i], sum(clusters$count))
  }, numeric(2)))
  clusters[, `:=`(ci_low = ci[, 1], ci_high = ci[, 2])]

  write_tsv(clusters[, .(rank, cluster_id, count, proportion, ci_low, ci_high,
                         passed_filter, length, is_reference,
                         n_snv, n_ins, n_del, variants)],
            file.path(outdir, "haplotypes.tsv"))
  top <- utils::head(clusters[passed_filter == TRUE], top_n)
  if (nrow(top) == 0) top <- utils::head(clusters, top_n)
  write_fasta(stats::setNames(top$consensus, sprintf("C%d_%s", top$cluster_id, top$variants)),
              file.path(outdir, "haplotypes.fasta"))

  var_rows <- data.table::rbindlist(lapply(seq_len(nrow(clusters)), function(i) {
    if (!isTRUE(clusters$annotated[i])) return(NULL)
    d <- clusters$diff_ops[[i]]
    if (nrow(d) == 0) return(NULL)
    d[, `:=`(cluster_id = clusters$cluster_id[i], count = clusters$count[i])]
    d
  }), use.names = TRUE)
  if (is.null(var_rows) || nrow(var_rows) == 0) {
    var_rows <- data.table::data.table(
      cluster_id = integer(0), count = integer(0), type = character(0),
      pos = integer(0), ref = character(0), alt = character(0)
    )
  } else {
    ref_seq <- ref$sequence
    var_rows[, Seq := mapply(
      function(p, r, a) seq_context(ref_seq, p, r, a),
      pos, ref, alt, SIMPLIFY = TRUE, USE.NAMES = FALSE
    )]
  }
  write_tsv(var_rows, file.path(outdir, "variants.tsv"))

  qc <- list(
    mode = "B",
    aligner = aligner,
    reference_label = ref_label %||% ref$name,
    reference_length = ref$length,
    n_reads_total = n_total,
    n_reads_primary = n_primary,
    n_reads_used = nrow(kept),
    mapping_rate = round(n_primary / max(n_total, 1), 6),
    mean_identity = round(mean(kept$identity), 6),
    mean_coverage = round(sum(kept$ref_span) / ref$length, 4),
    identity_cutoff = identity_cutoff,
    clustering_method = cl$method,
    consensus_method = consensus_method,
    decipher_version = if (requireNamespace("DECIPHER", quietly = TRUE)) {
      as.character(utils::packageVersion("DECIPHER"))
    } else {
      NA_character_
    },
    n_clusters = nrow(clusters),
    n_clusters_passed = sum(clusters$passed_filter),
    top1_proportion = round(clusters$proportion[1], 6),
    top1_is_reference = clusters$is_reference[1]
  )
  write_tsv(build_qc_table(qc), file.path(outdir, "qc.tsv"))
  run_manifest(outdir, "B", list(
    top_n = top_n, identity_cutoff = identity_cutoff,
    min_cluster_reads = min_cluster_reads, min_identity = min_identity,
    min_ref_coverage = min_ref_coverage, max_msa_seqs = max_msa_seqs,
    consensus_method = consensus_method, aligner = aligner, threads = threads
  ), ref, qc, extra = list(reads_md5 = safe_md5(reads_path)))

  if (!isTRUE(keep_intermediates) && !is.null(prep$bam)) {
    unlink(c(prep$bam, paste0(prep$bam, ".bai"), paste0(prep$bam, ".minimap2.log")))
  }
  invisible(list(haplotypes = clusters, variants = var_rows, qc = qc))
}
