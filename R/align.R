# ---------------------------------------------------------------------------
# minimap2 alignment and BAM parsing
# ---------------------------------------------------------------------------

align_reads <- function(reads_path, reference_path, out_bam,
                        threads = 4L, minimap2 = NULL, use_samtools = FALSE,
                        keep_sam = FALSE) {
  minimap2_bin <- minimap2 %||% nanoamp_tool_path("minimap2", required = TRUE)
  ensure_dir(dirname(out_bam))

  sam <- sub("\\.bam$", ".sam", out_bam)
  log_file <- paste0(out_bam, ".minimap2.log")
  status <- system2(
    minimap2_bin,
    c("-ax", "map-ont", "--cs", "-t", as.integer(threads),
      normalizePath(reference_path, mustWork = TRUE),
      normalizePath(reads_path, mustWork = TRUE)),
    stdout = sam, stderr = log_file
  )
  if (status != 0 || !file.exists(sam)) {
    stop(sprintf("minimap2 alignment failed with exit code %s", status), call. = FALSE)
  }

  if (isTRUE(use_samtools)) {
    samtools_bin <- nanoamp_tool_path("samtools", required = TRUE)
    status <- system2(
      samtools_bin,
      c("sort", "-@", as.integer(threads), "-o", out_bam, sam),
      stdout = FALSE, stderr = log_file
    )
    if (status != 0 || !file.exists(out_bam)) {
      stop(sprintf("samtools sort failed with exit code %s", status), call. = FALSE)
    }
    system2(samtools_bin, c("index", out_bam), stdout = FALSE, stderr = FALSE)
  } else {
    destination <- sub("\\.bam$", "", out_bam)
    bam <- Rsamtools::asBam(sam, destination, overwrite = TRUE, indexDestination = TRUE)
    if (!identical(
      normalizePath(bam, mustWork = FALSE),
      normalizePath(out_bam, mustWork = FALSE)
    )) {
      file.rename(bam, out_bam)
      if (file.exists(paste0(bam, ".bai"))) {
        file.rename(paste0(bam, ".bai"), paste0(out_bam, ".bai"))
      }
    }
  }
  if (!isTRUE(keep_sam)) unlink(sam)
  out_bam
}

parse_alignments <- function(bam) {
  flag <- Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = FALSE
  )
  param <- Rsamtools::ScanBamParam(
    what = c("qname", "flag", "rname", "pos", "cigar", "seq", "mapq"),
    tag = c("cs", "NM"),
    flag = flag
  )
  b <- Rsamtools::scanBam(bam, param = param)[[1]]
  if (length(b$qname) == 0) {
    return(data.table::data.table(
      read_id = character(0), flag = integer(0), ref_name = character(0),
      ref_start = integer(0), cigar = character(0), read_seq = character(0),
      mapq = integer(0), cs = character(0), nm = integer(0)
    ))
  }
  cs <- b$tag$cs
  if (is.null(cs)) cs <- rep(NA_character_, length(b$qname))
  nm <- b$tag$NM
  if (is.null(nm)) nm <- rep(NA_integer_, length(b$qname))
  data.table::data.table(
    read_id = as.character(b$qname),
    flag = as.integer(b$flag),
    ref_name = as.character(b$rname),
    ref_start = as.integer(b$pos),
    cigar = as.character(b$cigar),
    read_seq = toupper(as.character(b$seq)),
    mapq = as.integer(b$mapq),
    cs = as.character(cs),
    nm = as.integer(nm)
  )
}

count_bam_reads <- function(bam, mapped_only = FALSE) {
  flag <- if (mapped_only) Rsamtools::scanBamFlag(isUnmappedQuery = FALSE) else Rsamtools::scanBamFlag()
  as.integer(Rsamtools::countBam(bam, param = Rsamtools::ScanBamParam(flag = flag))$records)
}

alignment_read_strand <- function(flag) {
  ifelse(bitwAnd(as.integer(flag), 16L) > 0L, "-", "+")
}

cs_tokenize <- function(cs) {
  if (is.null(cs) || length(cs) == 0 || is.na(cs) || !nzchar(cs)) return(character(0))
  cs <- sub("^cs:Z:", "", cs)
  regmatches(cs, gregexpr("(:[0-9]+|\\*[A-Za-z]{2}|\\+[A-Za-z]+|-[A-Za-z]+)", cs))[[1]]
}

cs_ref_span <- function(cs) {
  tokens <- cs_tokenize(cs)
  if (!length(tokens)) return(0L)
  parts <- vapply(tokens, function(tok) {
    c0 <- substr(tok, 1, 1)
    if (c0 == ":") as.integer(sub("^:", "", tok))
    else if (c0 == "*") 1L
    else if (c0 == "-") nchar(tok) - 1L
    else 0L
  }, integer(1))
  as.integer(sum(parts))
}

prepare_alignment_stats <- function(aln, ref_len) {
  if (nrow(aln) == 0) {
    aln[, `:=`(ref_span = integer(0), identity = numeric(0),
               ref_end = integer(0), ref_cov = numeric(0),
               strand = character(0))]
    return(aln[])
  }
  aln <- data.table::copy(aln)
  aln[, strand := alignment_read_strand(flag)]
  aln[, ref_span := vapply(cs, cs_ref_span, integer(1))]
  aln[, nm := ifelse(is.na(nm), 0L, nm)]
  aln[, identity := 1 - nm / pmax(ref_span, 1L)]
  aln[, ref_end := ref_start + pmax(ref_span, 1L) - 1L]
  aln[, ref_cov := pmin(ref_span / ref_len, 1)]
  aln[]
}

filter_alignment_reads <- function(aln, min_identity = 0.90, min_ref_coverage = 0.90) {
  aln[!is.na(identity) & identity >= min_identity & ref_cov >= min_ref_coverage]
}

alignment_to_ops <- function(q_aligned, s_aligned) {
  q <- strsplit(q_aligned, "")[[1]]
  s <- strsplit(s_aligned, "")[[1]]
  n <- length(q)
  types <- character(0)
  poss <- integer(0)
  refs <- character(0)
  alts <- character(0)
  ref_pos <- 1L
  i <- 1L
  while (i <= n) {
    if (q[i] == s[i]) {
      if (q[i] != "-") ref_pos <- ref_pos + 1L
      i <- i + 1L
    } else if (q[i] != "-" && s[i] != "-") {
      types <- c(types, "snv")
      poss <- c(poss, ref_pos)
      refs <- c(refs, s[i])
      alts <- c(alts, q[i])
      ref_pos <- ref_pos + 1L
      i <- i + 1L
    } else if (q[i] == "-") {
      start <- ref_pos
      deleted <- character(0)
      while (i <= n && q[i] == "-" && s[i] != "-") {
        deleted <- c(deleted, s[i])
        ref_pos <- ref_pos + 1L
        i <- i + 1L
      }
      types <- c(types, "del")
      poss <- c(poss, start)
      refs <- c(refs, paste0(deleted, collapse = ""))
      alts <- c(alts, "")
    } else {
      anchor <- max(ref_pos - 1L, 0L)
      inserted <- character(0)
      while (i <= n && s[i] == "-" && q[i] != "-") {
        inserted <- c(inserted, q[i])
        i <- i + 1L
      }
      types <- c(types, "ins")
      poss <- c(poss, anchor)
      refs <- c(refs, "")
      alts <- c(alts, paste0(inserted, collapse = ""))
    }
  }
  if (length(types) == 0) return(empty_ops())
  data.table::data.table(type = types, pos = poss, ref = refs, alt = alts)
}

align_reads_r <- function(reads_path, reference_path, threads = 1L) {
  ref <- read_reference(reference_path)
  fq <- read_fastq(reads_path)
  n <- nrow(fq)
  if (n == 0) {
    return(list(
      aln = data.table::data.table(
        read_id = character(0), flag = integer(0), ref_name = character(0),
        ref_start = integer(0), cigar = character(0), read_seq = character(0),
        mapq = integer(0), cs = character(0), nm = integer(0),
        strand = character(0), ref_span = integer(0), identity = numeric(0),
        ref_end = integer(0), ref_cov = numeric(0)
      ),
      read_vars = data.table::data.table(
        type = character(0), pos = integer(0), ref = character(0),
        alt = character(0), read_id = character(0), strand = character(0),
        key = character(0)
      )
    ))
  }

  aln_rows <- vector("list", n)
  var_rows <- vector("list", n)
  for (i in seq_len(n)) {
    read <- fq$sequence[i]
    read_rc <- reverse_complement(read)
    a_fwd <- pa_pairwise_alignment(
      read, ref$sequence, type = "global", gapOpening = 5, gapExtension = 1
    )
    a_rev <- pa_pairwise_alignment(
      read_rc, ref$sequence, type = "global", gapOpening = 5, gapExtension = 1
    )
    if (pa_score(a_rev) > pa_score(a_fwd)) {
      aln <- a_rev
      strand <- "-"
    } else {
      aln <- a_fwd
      strand <- "+"
    }
    q <- as.character(pa_aligned(pa_pattern(aln)))
    s <- as.character(pa_aligned(pa_subject(aln)))
    ops <- alignment_to_ops(q, s)
    nm <- if (nrow(ops) == 0) {
      0L
    } else {
      sum(ops$type == "snv") +
        sum(nchar(ops$ref[ops$type == "del"])) +
        sum(nchar(ops$alt[ops$type == "ins"]))
    }
    aln_rows[[i]] <- data.table::data.table(
      read_id = fq$read_id[i],
      flag = ifelse(strand == "-", 16L, 0L),
      ref_name = ref$name,
      ref_start = 1L,
      cigar = NA_character_,
      read_seq = read,
      mapq = 60L,
      cs = NA_character_,
      nm = as.integer(nm),
      strand = strand,
      ref_span = ref$length,
      identity = 1 - nm / ref$length,
      ref_end = ref$length,
      ref_cov = 1
    )
    if (nrow(ops) > 0) {
      ops[, read_id := fq$read_id[i]]
      ops[, strand := strand]
      var_rows[[i]] <- ops
    }
  }
  aln <- data.table::rbindlist(aln_rows, use.names = TRUE)
  vars <- data.table::rbindlist(var_rows, use.names = TRUE)
  if (nrow(vars) == 0) {
    vars <- data.table::data.table(
      type = character(0), pos = integer(0), ref = character(0),
      alt = character(0), read_id = character(0), strand = character(0)
    )
  }
  vars[, key := variant_key(type, pos, ref, alt)]
  list(aln = aln, read_vars = vars)
}

prepare_alignment_data <- function(reads_path, reference_path, outdir,
                                   aligner = c("minimap2", "r"),
                                   threads = 4L, use_samtools = FALSE) {
  aligner <- match.arg(aligner, c("minimap2", "r"))
  if (identical(aligner, "r")) {
    res <- align_reads_r(reads_path, reference_path, threads = threads)
    return(list(aln = res$aln, read_vars = res$read_vars, bam = NULL, aligner = aligner))
  }
  ref <- read_reference(reference_path)
  bam <- file.path(outdir, "alignments.bam")
  align_reads(reads_path, reference_path, bam, threads = threads, use_samtools = use_samtools)
  aln <- prepare_alignment_stats(parse_alignments(bam), ref$length)
  list(aln = aln, read_vars = extract_read_variants(aln), bam = bam, aligner = aligner)
}
