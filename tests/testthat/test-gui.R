test_that("GUI server runs a Mode C analysis", {
  skip_if_not_installed("shiny")

  td <- tempfile("nanoamp_gui_")
  dir.create(td, recursive = TRUE)
  ref <- make_random_seq(300, seed = 7)
  seqs <- c(ref, reverse_complement(ref), ref)
  fq <- file.path(td, "reads.fastq")
  ref_fa <- file.path(td, "ref.fa")
  write_test_fastq(seqs, fq)
  write_test_ref(ref, ref_fa)

  outdir <- file.path(td, "out")

  suppressWarnings(
    shiny::testServer(nanoamp:::nanoamp_gui_server, {
      session$setInputs(
        reads = list(datapath = fq, name = "reads.fastq"),
        reference = list(datapath = ref_fa, name = "ref.fa"),
        outdir = outdir, mode = "C", top_n = 5,
        min_reads = 3, min_freq = 0.02, min_identity = 0.9,
        identity_cutoff = 0.99, min_cluster_reads = 2,
        consensus_method = "decipher", aligner = "r",
        threads = 2, keep_intermediates = FALSE
      )
      session$setInputs(run = 1)
      expect_match(output$status, "Analysis complete")
    })
  )

  expect_true(file.exists(file.path(outdir, "haplotypes.tsv")))
  expect_true(file.exists(file.path(outdir, "variants.tsv")))
  expect_true(file.exists(file.path(outdir, "qc.tsv")))
})
