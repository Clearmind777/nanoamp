# Installing dependencies for nanoamp

This guide explains how to install and configure `minimap2`, `samtools` and
the R packages required by `nanoamp` on Linux and Windows.

## 1. What is needed

| Dependency | Type | Required for |
|---|---|---|
| `minimap2` | external command | Modes A and B (read alignment) |
| `samtools` | optional external command | Compatibility fallback; Rsamtools is used by default |
| R packages | R packages | Core analysis |
| `DECIPHER` | optional R package | Mode B de novo clustering |
| `shiny`, `bslib`, `DT` | optional R packages | GUI |

Mode C (raw exact matching) does not need `minimap2`. `samtools` is never
required because `Rsamtools::asBam()` handles SAM to BAM conversion.

## 1.1 Tool resolution and the R-native fallback

`nanoamp` looks for tools in this order:

1. `NANOAMP_MINIMAP2` / `NANOAMP_SAMTOOLS` environment variables;
2. `<dependence>/<os>-<arch>/bin/`, where `<dependence>` is
   `$NANOAMP_DEPENDENCE_DIR` or a `dependence/` (or `03_dependence/`) directory
   in the current working tree;
3. `PATH`.

The package itself ships no binaries, so on most installations `minimap2` is
simply found on `PATH`.

On platforms without a minimap2 binary (Windows, ARM), use the R-native
backend:

```r
run_haplotype_analysis(..., aligner = "r")
```

`samtools` is optional because `Rsamtools::asBam()` converts SAM to BAM by
default. Use `use_samtools = TRUE` only if you need the samtools path.

## 2. Linux

### Option A: conda / mamba (recommended)

```bash
conda create -n nanoamp -c conda-forge -c bioconda minimap2 samtools
conda activate nanoamp

which minimap2
minimap2 --version

which samtools
samtools --version
```

Then install R packages inside the same R environment (see section 4).

### Option B: system packages

Debian / Ubuntu:

```bash
sudo apt-get update
sudo apt-get install -y minimap2 samtools
```

CentOS / Rocky / AlmaLinux:

```bash
sudo dnf install -y minimap2 samtools
```

### Verify

```bash
which minimap2
which samtools

Rscript -e 'library(nanoamp); nanoamp_cli("doctor")'
```

## 3. Windows

There are no official Windows binaries for `minimap2` or `samtools`, and
conda-forge / bioconda do not provide win-64 builds for them. Choose one of the
following options.

### Option A: R-native backend (recommended)

Use the R-native alignment backend, which needs no external tool:

```r
run_haplotype_analysis(..., aligner = "r")
```

This is the simplest Windows setup and works for small and medium amplicons.

### Option B: WSL2 (recommended when minimap2 speed is needed)

1. Install WSL2 and Ubuntu.
2. Follow the Linux installation instructions in section 2 inside WSL.
3. Run nanoamp inside WSL, pointing it to the data files.

### Option C: third-party Windows binaries (optional)

If you have third-party builds, place them here:

```text
dependence\windows-x86_64\bin\minimap2.exe
dependence\windows-x86_64\bin\samtools.exe    # optional
```

nanoamp resolves these files automatically. `samtools.exe` is optional because
`Rsamtools` handles SAM to BAM conversion by default.

### Configure PATH (only for Option C)

1. Open **System Properties -> Environment Variables**.
2. Edit the `Path` variable.
3. Add the folder that contains `minimap2.exe` (and optionally `samtools.exe`).
4. Click OK and **restart RStudio / terminal**.
5. Verify in R:

```r
Sys.which("minimap2")
Sys.which("samtools")   # optional
library(nanoamp)
nanoamp_cli("doctor")
```

Alternatively, skip PATH entirely by placing the binaries under
`dependence\windows-x86_64\bin\` or by setting `NANOAMP_MINIMAP2`.

### Windows pitfalls

- **Do not rely on `conda install minimap2 samtools`**: there is no win-64
  build for these packages.
- **PATH not refreshed**: restart RStudio after changing `PATH`.
- **Spaces or non-ASCII characters in paths**: prefer `C:\tools\...`.
- **Windows SmartScreen**: allow the downloaded binaries if prompted.
- **Multiple R installations**: check `Rscript -e 'cat(R.home())'` and make
  sure the package is installed into the R you actually use.

## 4. R packages

Required:

```r
install.packages(c(
  "Biostrings", "Rsamtools", "ShortRead", "IRanges", "Matrix",
  "data.table", "optparse", "jsonlite", "readxl"
))
```

If the Bioconductor packages are not available from CRAN:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "Rsamtools", "ShortRead", "IRanges"))
```

Optional:

```r
BiocManager::install("DECIPHER")          # Mode B
install.packages(c("shiny", "DT"))  # GUI
```

## 5. Verification checklist

```r
library(nanoamp)
nanoamp_cli("doctor")
```

Expected output:

```text
nanoamp version: 0.1.0
R version: ...
Rscript: ...
  Biostrings   TRUE
  ...
  DECIPHER     TRUE
  minimap2     /path/to/minimap2
  samtools     /path/to/samtools
```

What to check:

- `minimap2` shows a path, not `NOT FOUND`;
- `samtools` is optional and may show `NOT FOUND` unless
  `use_samtools = TRUE`;
- R packages show `TRUE`;
- `DECIPHER` may be `FALSE`: Mode B still works with a fallback, but DECIPHER
  is recommended.

## 6. Dependency reduction status

The following improvements are already implemented:

1. `Rsamtools::asBam()` performs SAM to BAM conversion by default, so the
   `samtools` command is optional;
2. `aligner = "r"` provides an R-native pairwise alignment backend for small
   and medium datasets, and for Windows / ARM platforms without minimap2;
3. `minimap2` remains the recommended backend for large datasets.

Set `use_samtools = TRUE` only if you explicitly need the samtools path.
