# ---------------------------------------------------------------------------
# FASTA / FASTQ / table input and output
# ---------------------------------------------------------------------------

read_fastq <- function(path) {
  if (!file.exists(path)) stop(sprintf("FASTQ file not found: %s", path), call. = FALSE)
  fq <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, "rt")
    on.exit(close(con), add = TRUE)
    ShortRead::readFastq(con)
  } else {
    ShortRead::readFastq(path)
  }
  if (length(fq) == 0) {
    return(data.table::data.table(
      read_id = character(0), sequence = character(0), quality = character(0)
    ))
  }
  qmat <- methods::as(Biostrings::quality(fq), "matrix")
  qual_chr <- if (nrow(qmat) == 0) character(0) else {
    vapply(seq_len(nrow(qmat)), function(i) {
      q <- as.integer(qmat[i, ])
      q[is.na(q)] <- 0L
      q <- pmax(pmin(q, 93L), 0L)
      rawToChar(as.raw(q + 33L))
    }, character(1))
  }
  data.table::data.table(
    read_id = as.character(ShortRead::id(fq)),
    sequence = as.character(ShortRead::sread(fq)),
    quality = qual_chr
  )
}

count_fastq_reads <- function(path) {
  if (!file.exists(path)) stop(sprintf("FASTQ file not found: %s", path), call. = FALSE)
  n <- ShortRead::countFastq(path)
  as.integer(n$records[1])
}

read_reference <- function(path) {
  if (!file.exists(path)) stop(sprintf("Reference sequence not found: %s", path), call. = FALSE)
  x <- Biostrings::readDNAStringSet(path)
  if (length(x) == 0) stop(sprintf("Reference sequence is empty: %s", path), call. = FALSE)
  seq <- toupper(as.character(x[[1]]))
  list(
    name = names(x)[1] %||% "reference",
    sequence = seq,
    length = nchar(seq),
    path = normalizePath(path, mustWork = TRUE),
    md5 = safe_md5(path)
  )
}

write_fasta <- function(sequences, path) {
  if (length(sequences) == 0) {
    file.create(path)
    return(invisible(path))
  }
  x <- Biostrings::DNAStringSet(toupper(sequences))
  names(x) <- names(sequences)
  Biostrings::writeXStringSet(x, path)
  invisible(path)
}

read_company_variants <- function(path) {
  if (is.null(path) || is.na(path) || !file.exists(path)) return(NULL)
  df <- tryCatch(
    readxl::read_excel(path, sheet = 1),
    error = function(e) {
      log_warn("Could not read company variant table ", path, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  pick <- function(candidates) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) == 0) return(rep(NA, nrow(df)))
    df[[hit[1]]]
  }
  out <- data.table::data.table(
    company_pos = suppressWarnings(as.integer(
      pick(c("\u7a81\u53d8\u78b1\u57fa\u4f4d\u7f6e", "Pos", "position"))
    )),
    company_ref = as.character(pick(c("\u8f93\u51fa\u78b1\u57fa", "Ref", "ref"))),
    company_alt = as.character(pick(c("\u53d8\u5f02\u78b1\u57fa", "Alt", "alt"))),
    company_type = as.character(pick(c("\u53d8\u5f02\u7c7b\u578b", "type"))),
    company_freq = suppressWarnings(as.numeric(
      pick(c("\u53d8\u5f02\u6bd4\u4f8b(%)", "Freq", "freq"))
    ))
  )
  out[]
}
