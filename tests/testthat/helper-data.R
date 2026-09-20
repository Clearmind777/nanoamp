write_test_fastq <- function(seqs, path, prefix = "r") {
  ids <- sprintf("@%s%04d", prefix, seq_along(seqs))
  quals <- strrep("I", nchar(seqs))
  lines <- as.vector(rbind(ids, seqs, "+", quals))
  writeLines(lines, path)
}

write_test_ref <- function(seq, path) {
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(stats::setNames(seq, "ref")), path
  )
}

make_random_seq <- function(n, seed = 1L) {
  set.seed(seed)
  paste0(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
}
