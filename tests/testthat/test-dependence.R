test_that("tool resolver reports a usable platform", {
  expect_match(nanoamp_platform(), "^(linux|windows|macos)-")
})

test_that("tool resolver can find minimap2 when available", {
  path <- nanoamp_tool_path("minimap2", required = FALSE)
  if (!is.null(path)) {
    expect_true(file.exists(path))
  } else {
    succeed("minimap2 is not available on this platform")
  }
})

test_that("R-native aligner runs on a small synthetic dataset", {
  skip_if_not_installed("Biostrings")
  td <- tempfile("nanoamp_r_aligner_")
  dir.create(td)
  ref <- make_random_seq(300, seed = 42)
  mut <- ref
  substr(mut, 150, 150) <- ifelse(substr(mut, 150, 150) == "A", "T", "A")
  seqs <- c(rep(ref, 6), rep(mut, 4))
  fq <- file.path(td, "reads.fastq")
  ref_fa <- file.path(td, "ref.fa")
  write_test_fastq(seqs, fq)
  write_test_ref(ref, ref_fa)

  res <- run_mode_a(
    fq, ref_fa, file.path(td, "out"),
    top_n = 5, min_reads = 2, min_freq = 0.2,
    min_identity = 0.8, aligner = "r"
  )
  expect_true(nrow(res$haplotypes) >= 1)
  expect_equal(sum(res$haplotypes$count), 10)
})
