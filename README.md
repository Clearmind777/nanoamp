# nanoamp

`nanoamp` is an R package for analyzing Oxford Nanopore reads from PCR
amplicons. It aligns reads to a target sequence, corrects sequencing errors,
reconstructs haplotypes, and reports the most abundant sequences with counts
and proportions.

The package is designed for questions such as:

- How many reads match the intended PCR product exactly?
- What other sequences are present, and at what proportions?
- Which variants are real, and which are nanopore sequencing errors?
- Which haplotype carries which combination of variants?

## Installation

### 1. Install dependencies

```r
install.packages(c(
  "Biostrings", "Rsamtools", "ShortRead", "IRanges", "Matrix",
  "data.table", "optparse", "jsonlite", "readxl"
))

# Recommended for Mode B (de novo clustering)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DECIPHER")

# Required for aligner = "r" on Bioconductor >= 3.19, which moved
# pairwiseAlignment() out of Biostrings
BiocManager::install("pwalign")

# Optional GUI
install.packages(c("shiny", "DT"))
```

### 2. Install `nanoamp`

From GitHub:

```r
remotes::install_github("Clearmind777/nanoamp")
```

From a local checkout:

```bash
R CMD INSTALL .
```

During development:

```r
devtools::install()
```

### 3. Install external tools

`minimap2` is resolved from the `NANOAMP_MINIMAP2` environment variable, from
`$NANOAMP_DEPENDENCE_DIR/<os>-<arch>/bin/`, and finally from `PATH`. The
package ships no binaries, so a normal installation just needs `minimap2` on
`PATH`. `samtools` is optional: SAM to BAM conversion uses
`Rsamtools::asBam()` by default.

Detailed installation and `PATH` configuration instructions for Linux and
Windows are in [inst/docs/INSTALL_DEPENDENCIES.md](inst/docs/INSTALL_DEPENDENCIES.md).

```bash
minimap2 --version
# optional:
# samtools --version
```

Check everything from R:

```r
library(nanoamp)
nanoamp_cli("doctor")
```

## Quick start

```r
library(nanoamp)

res <- run_haplotype_analysis(
  reads     = "sample.fastq",
  reference = "target.fa",
  outdir    = "results/sampleA",
  mode      = "A",
  top_n     = 20
)

# Top haplotypes
res$haplotypes

# Candidate variants
res$variants

# QC metrics
res$qc
```

Basic input requirements:

- `reads`: FASTQ or FASTQ.GZ, single-end nanopore reads;
- `reference`: FASTA containing the intended amplicon sequence;
- `outdir`: output directory (created automatically).

## Analysis modes

### Mode A: reference-guided correction (recommended)

Mode A aligns reads to the target sequence, discovers candidate variants,
treats differences that do not pass the variant filters as sequencing errors,
and groups reads by their corrected sequence.

Use Mode A when:

- a reliable target sequence is available;
- you need quantitative haplotype proportions;
- you want to distinguish real variants from nanopore errors.

### Mode B: de novo clustering (exploratory)

Mode B clusters reads with `DECIPHER::Clusterize` and builds a polished
consensus for each cluster using `DECIPHER::AlignSeqs` followed by majority
voting.

Use Mode B when:

- no reliable reference is available;
- you want a data-driven overview of the main sequence groups;
- you accept that haplotypes differing by less than the sequencing error rate
  may not be resolved.

If `DECIPHER` is unavailable, Mode B falls back to variant-pattern greedy
clustering and records this in `qc.tsv`.

### Mode C: raw exact matching (diagnostic)

Mode C counts raw reads that match the reference exactly on either strand. It
is useful for demonstrating the effect of nanopore errors, but it is not
recommended for quantitative haplotype analysis.

## Parameters

Default parameters can be inspected with:

```r
nanoamp_defaults()
```

| Parameter | Default | Description |
|---|---:|---|
| `top_n` | 20 | Number of top haplotypes to report |
| `min_reads` | 3 | Minimum supporting reads for a candidate variant |
| `min_freq` | 0.02 | Minimum variant frequency |
| `min_identity` | 0.90 | Minimum read identity to the reference |
| `min_ref_coverage` | 0.90 | Minimum fraction of the reference covered by a read |
| `homopolymer` | 4 | Homopolymer length threshold for filtering |
| `strand_bias` | 0.90 | Strand bias threshold |
| `identity_cutoff` | 0.99 | Mode B clustering identity cutoff |
| `min_cluster_reads` | 2 | Mode B minimum cluster size |
| `max_msa_seqs` | 100 | Maximum sequences per consensus alignment |
| `consensus_method` | `"decipher"` | `"decipher"` or `"medoid"` |
| `aligner` | `"minimap2"` | `"minimap2"` or `"r"` (R-native fallback) |
| `use_samtools` | `FALSE` | Use samtools instead of Rsamtools for SAM to BAM |
| `threads` | 4 | Number of threads |
| `keep_intermediates` | `TRUE` | Keep BAM and other intermediate files |

## Output files

```text
outdir/
|-- haplotypes.tsv
|-- haplotypes.fasta
|-- variants.tsv
|-- qc.tsv
|-- run_manifest.json
`-- alignments.bam(.bai)     # Modes A and B, when keep_intermediates = TRUE
```

### haplotypes.tsv

| Column | Description |
|---|---|
| `rank` | Rank by supporting read count |
| `haplotype_id` / `cluster_id` | Haplotype or cluster identifier |
| `count` | Supporting reads |
| `proportion` | Fraction of assigned reads |
| `ci_low`, `ci_high` | 95% Wilson confidence interval |
| `is_reference` | Whether the sequence matches the reference |
| `n_snv`, `n_ins`, `n_del` | Number of variants |
| `length` | Haplotype length |
| `variants` | Variant description; `.` means no variant |

### variants.tsv

Mode A uses a company-compatible layout:

```text
Chr  Pos  Ref  Alt  DP  Ref_dp  Alt_dp  Freq  DP4  Seq  Filter_Status  Filter_Reason
```

- `Freq` is a fraction between 0 and 1;
- `-` in `Ref` or `Alt` represents an insertion or deletion;
- `Filter_Status` is `PASS` or `FILTERED`.

### qc.tsv

Two columns, `metric` and `value`, including read counts, mapping rate, mean
identity, coverage, clustering method, consensus method, and DECIPHER version.

## Command line interface

The package ships a CLI built from the same R code. It is available as the
`nanoamp` command after installation, and as `nanoamp_cli()` from R.

```bash
nanoamp doctor

nanoamp call \
  --reads sample.fastq \
  --reference target.fa \
  --mode A \
  --top-n 20 \
  --outdir results/sampleA

nanoamp batch \
  --sample-sheet samples.tsv \
  --mode A \
  --outdir results/batch
```

The batch sample sheet is a TSV with at least:

```text
sample	reads	reference
```

Optional columns: `ref_label`.

Use `--aligner r` to select the R-native alignment backend on platforms
without minimap2.

### Install the `nanoamp` command

```bash
sh "$(Rscript --vanilla -e 'cat(system.file("scripts", "install_cli.sh", package = "nanoamp"))')" ~/.local/bin
export PATH="$HOME/.local/bin:$PATH"
nanoamp doctor
```

Alternatively, call the CLI directly from R:

```r
library(nanoamp)
nanoamp_cli(c("call", "--reads", "sample.fastq", "--reference", "target.fa",
              "--outdir", "results/sampleA"))
```

## Graphical user interface

Launch the Shiny GUI:

```r
library(nanoamp)
nanoamp_gui()
```

The GUI provides file pickers, mode selection, advanced parameters, a run
button, a captured log, interactive haplotype/variant tables, QC output and
download buttons.

On Windows the shipped launcher can be used instead:

```bat
Rscript -e "library(nanoamp); nanoamp_gui()"
```

To obtain the Shiny app object without starting a server:

```r
app <- nanoamp_gui_app()
```

### External tools and the R-native backend

`aligner = "minimap2"` (the default) needs a `minimap2` executable. On a
machine without minimap2, use the R-native backend instead:

```r
run_haplotype_analysis(..., aligner = "r")
```

The R-native backend uses Biostrings pairwise alignment and requires no
external tool. It is slower and is intended for small and medium amplicons.

`samtools` is optional: SAM -> BAM conversion uses `Rsamtools::asBam()` by
default. Set `use_samtools = TRUE` only if you explicitly want the samtools
path.

## Tests

```r
# Unit tests
devtools::test()

# Full package check
devtools::check()
```

The package is verified with `R CMD check` and currently passes with
`Status: OK`.

## Troubleshooting

| Symptom | Solution |
|---|---|
| `minimap2` not found | Install minimap2 and add it to `PATH`, or set `NANOAMP_MINIMAP2` |
| `samtools` not found | Usually not needed: `Rsamtools` is the default. Install samtools only if `use_samtools = TRUE` |
| Mode B is slow | Reduce `max_msa_seqs`, increase `threads`, or use `mode = "A"` |
| Mode B cannot separate close haplotypes | This is expected below the sequencing error rate; use Mode A |
| `DECIPHER` not installed | Mode B falls back to greedy clustering; install DECIPHER for better results |
| `pairwiseAlignment` is not an exported object from Biostrings | Bioconductor >= 3.19 moved it to `pwalign`; install it with `BiocManager::install("pwalign")` |
| All proportions are low in Mode C | Nanopore reads contain errors; use Mode A |

## License

MIT.
