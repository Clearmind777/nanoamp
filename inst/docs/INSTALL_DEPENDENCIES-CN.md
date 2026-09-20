# nanoamp 依赖安装与配置

本文说明如何在 Linux 和 Windows 上安装并配置 `minimap2`、`samtools`
以及 `nanoamp` 需要的 R 包。

## 1. 需要哪些依赖

| 依赖 | 类型 | 用途 |
|---|---|---|
| `minimap2` | 外部命令 | 方案 A/B 的 reads 比对 |
| `samtools` | 可选外部命令 | 兼容后备；默认使用 Rsamtools |
| R 包 | R 包 | 核心分析 |
| `DECIPHER` | 可选 R 包 | 方案 B 从头聚类 |
| `shiny`、`bslib`、`DT` | 可选 R 包 | GUI |

方案 C（原始精确匹配）不需要 `minimap2`。`samtools` 从来不是必需依赖，
因为 `Rsamtools::asBam()` 已经负责 SAM→BAM。

## 1.1 工具解析顺序与 R 内后备

`nanoamp` 查找工具的顺序：

1. 环境变量 `NANOAMP_MINIMAP2` / `NANOAMP_SAMTOOLS`；
2. `<dependence>/<os>-<arch>/bin/`，其中 `<dependence>` 是
   `$NANOAMP_DEPENDENCE_DIR`，或当前工作目录树里的 `dependence/`（或
   `03_dependence/`）目录；
3. `PATH`。

R 包本身不附带二进制，因此绝大多数情况下只要把 `minimap2` 放到 `PATH` 即可。

没有 minimap2 二进制的平台（Windows、ARM）可以使用 R 内后端：

```r
run_haplotype_analysis(..., aligner = "r")
```

`samtools` 不是必需依赖：默认用 `Rsamtools::asBam()` 完成 SAM→BAM。
只有显式设置 `use_samtools = TRUE` 时才走 samtools。

## 2. Linux

### 方式 A：conda / mamba（推荐）

```bash
conda create -n nanoamp -c conda-forge -c bioconda minimap2 samtools
conda activate nanoamp

which minimap2
minimap2 --version

which samtools
samtools --version
```

然后在这个 R 环境中安装 R 包（见第 4 节）。

### 方式 B：系统包管理器

Debian / Ubuntu：

```bash
sudo apt-get update
sudo apt-get install -y minimap2 samtools
```

CentOS / Rocky / AlmaLinux：

```bash
sudo dnf install -y minimap2 samtools
```

### 验证

```bash
which minimap2
which samtools

Rscript -e 'library(nanoamp); nanoamp_cli("doctor")'
```

## 3. Windows

`minimap2` 和 `samtools` 都没有官方 Windows 二进制，conda-forge / bioconda
也不提供 win-64 构建。请从以下方案中选择。

### 方式 A：R 内后端（推荐）

使用 R 内比对后端，不需要任何外部工具：

```r
run_haplotype_analysis(..., aligner = "r")
```

这是 Windows 上最简单的方案，适合中小扩增子。

### 方式 B：WSL2（需要 minimap2 速度时推荐）

1. 安装 WSL2 和 Ubuntu；
2. 在 WSL 内按第 2 节的 Linux 步骤安装 minimap2；
3. 在 WSL 中运行 nanoamp，并指向相应数据文件。

### 方式 C：第三方 Windows 二进制（可选）

如果你有第三方 Windows 构建，放到：

```text
dependence\windows-x86_64\bin\minimap2.exe
dependence\windows-x86_64\bin\samtools.exe    # 可选
```

nanoamp 会自动解析这些文件。`samtools.exe` 是可选的，因为默认用
`Rsamtools` 完成 SAM→BAM。

### 在 Windows 配置 PATH（仅方式 C 需要）

1. 打开 **系统属性 -> 环境变量**；
2. 编辑 `Path`；
3. 加入包含 `minimap2.exe`（以及可选的 `samtools.exe`）的目录；
4. 确定后**重启 RStudio / 终端**；
5. 在 R 中验证：

```r
Sys.which("minimap2")
Sys.which("samtools")   # 可选
library(nanoamp)
nanoamp_cli("doctor")
```

也可以跳过 PATH，直接把二进制放到
`dependence\windows-x86_64\bin\`，或者设置 `NANOAMP_MINIMAP2`。

### Windows 常见坑

- **不要依赖 `conda install minimap2 samtools`**：这两个包没有 win-64 构建；
- **PATH 未刷新**：修改 PATH 后必须重启 RStudio；
- **路径有空格或中文**：建议放到 `C:\tools\...`；
- **Windows SmartScreen**：如果提示拦截下载的 exe，需要手动允许；
- **多个 R 版本**：用 `Rscript -e 'cat(R.home())'` 确认当前 R，
  并确保 R 包装到了同一个 R 中。

## 4. R 包

必需：

```r
install.packages(c(
  "Biostrings", "Rsamtools", "ShortRead", "IRanges", "Matrix",
  "data.table", "optparse", "jsonlite", "readxl"
))
```

如果 Bioconductor 包无法从 CRAN 安装：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "Rsamtools", "ShortRead", "IRanges"))
```

可选：

```r
BiocManager::install("DECIPHER")             # 方案 B
install.packages(c("shiny", "DT"))  # GUI
```

## 5. 验证清单

```r
library(nanoamp)
nanoamp_cli("doctor")
```

期望输出：

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

需要确认：

- `minimap2` 显示路径而不是 `NOT FOUND`；
- `samtools` 是可选的，除非 `use_samtools = TRUE`，否则显示 `NOT FOUND` 也没关系；
- R 包显示 `TRUE`；
- `DECIPHER` 可以是 `FALSE`：方案 B 会自动降级，但建议安装。

## 6. 依赖缩减现状

以下改进已经实现：

1. 默认用 `Rsamtools::asBam()` 完成 SAM→BAM，`samtools` 命令变成可选；
2. `aligner = "r"` 提供 R 内成对比对后端，适合中小数据，以及没有
   minimap2 的 Windows / ARM 平台；
3. `minimap2` 仍然是大数据量下的推荐后端。

只有显式设置 `use_samtools = TRUE` 时才需要 samtools。
