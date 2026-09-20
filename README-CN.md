# nanoamp

`nanoamp` 是一个 R 包，用于分析纳米孔 PCR 产物的测序数据。它会把 reads 比对到
目的序列，校正测序错误，重建单倍型，并输出数量最多、比例最高的序列。

这个包面向以下问题：

- 有多少 reads 与预期 PCR 产物完全一致？
- 还存在哪些序列？各自比例是多少？
- 哪些变异是真实变异，哪些只是纳米孔测序错误？
- 哪个单倍型携带了哪一组变异？

## 安装

### 1. 安装依赖

```r
install.packages(c(
  "Biostrings", "Rsamtools", "ShortRead", "IRanges", "Matrix",
  "data.table", "optparse", "jsonlite", "readxl"
))

# 方案 B（从头聚类）推荐安装
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DECIPHER")

# 在 Bioconductor >= 3.19 上，pairwiseAlignment() 已移出 Biostrings，
# 使用 aligner = "r" 时必须安装 pwalign
BiocManager::install("pwalign")

# 图形界面（可选）
install.packages(c("shiny", "DT"))
```

### 2. 安装 `nanoamp`

从 GitHub 安装：

```r
remotes::install_github("Clearmind777/nanoamp")
```

从本地源码安装：

```bash
R CMD INSTALL .
```

开发时：

```r
devtools::install()
```

### 3. 安装外部工具

`minimap2` 的解析顺序是：环境变量 `NANOAMP_MINIMAP2` →
`$NANOAMP_DEPENDENCE_DIR/<os>-<arch>/bin/` → `PATH`。R 包本身不附带二进制，
因此常规安装只需把 `minimap2` 放到 `PATH` 上。`samtools` 是可选的：默认由
`Rsamtools::asBam()` 完成 SAM→BAM。

Linux 与 Windows 下的详细安装与 `PATH` 配置说明见
[inst/docs/INSTALL_DEPENDENCIES-CN.md](inst/docs/INSTALL_DEPENDENCIES-CN.md)。

```bash
minimap2 --version
# 可选：
# samtools --version
```

在 R 中检查环境：

```r
library(nanoamp)
nanoamp_cli("doctor")
```

## 快速开始

```r
library(nanoamp)

res <- run_haplotype_analysis(
  reads     = "sample.fastq",
  reference = "target.fa",
  outdir    = "results/sampleA",
  mode      = "A",
  top_n     = 20
)

# 主要单倍型
res$haplotypes

# 候选变异
res$variants

# QC 指标
res$qc
```

输入要求：

- `reads`：FASTQ 或 FASTQ.GZ，单端纳米孔 reads；
- `reference`：包含目的扩增子序列的 FASTA；
- `outdir`：输出目录（不存在会自动创建）。

## 分析方案

### 方案 A：参考序列引导校正（推荐）

方案 A 先把 reads 比对到目的序列，发现候选变异，把未通过变异过滤的差异当作
测序错误，最后按校正后的序列对 reads 分组。

以下情况建议用方案 A：

- 有可靠的参考序列；
- 需要定量的单倍型比例；
- 需要区分真实变异与纳米孔错误。

### 方案 B：从头聚类（探索性）

方案 B 用 `DECIPHER::Clusterize` 对 reads 聚类，再用
`DECIPHER::AlignSeqs` 加多数投票，为每个簇构建 polished 一致性序列。

以下情况建议用方案 B：

- 没有可靠参考序列；
- 想要数据驱动地了解主要序列分组；
- 接受“差异小于测序错误率的单倍型可能无法分开”。

如果 `DECIPHER` 不可用，方案 B 会退化为按变异模式的贪心聚类，并在
`qc.tsv` 中记录。

### 方案 C：原始精确匹配（诊断用）

方案 C 统计在正链或反链上与参考序列完全一致的原始 reads。它适合演示纳米孔错误
的影响，但不建议用于定量单倍型分析。

## 参数

查看默认参数：

```r
nanoamp_defaults()
```

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `top_n` | 20 | 输出的单倍型数量 |
| `min_reads` | 3 | 候选变异的最小支持 reads 数 |
| `min_freq` | 0.02 | 最小变异频率 |
| `min_identity` | 0.90 | reads 与参考序列的最小一致性 |
| `min_ref_coverage` | 0.90 | 单条 read 覆盖参考序列的最小比例 |
| `homopolymer` | 4 | 同聚物长度过滤阈值 |
| `strand_bias` | 0.90 | 链偏好阈值 |
| `identity_cutoff` | 0.99 | 方案 B 聚类一致性阈值 |
| `min_cluster_reads` | 2 | 方案 B 最小簇大小 |
| `max_msa_seqs` | 100 | 每次一致性比对使用的最大序列数 |
| `consensus_method` | `"decipher"` | `"decipher"` 或 `"medoid"` |
| `aligner` | `"minimap2"` | `"minimap2"` 或 `"r"`（R 内后端） |
| `use_samtools` | `FALSE` | 用 samtools 而不是 Rsamtools 完成 SAM→BAM |
| `threads` | 4 | 线程数 |
| `keep_intermediates` | `TRUE` | 保留 BAM 等中间文件 |

## 输出文件

```text
outdir/
|-- haplotypes.tsv
|-- haplotypes.fasta
|-- variants.tsv
|-- qc.tsv
|-- run_manifest.json
`-- alignments.bam(.bai)     # 方案 A/B，keep_intermediates = TRUE 时生成
```

### haplotypes.tsv

| 列名 | 说明 |
|---|---|
| `rank` | 按支持 reads 数排名 |
| `haplotype_id` / `cluster_id` | 单倍型或簇编号 |
| `count` | 支持 reads 数 |
| `proportion` | 占已分配 reads 的比例 |
| `ci_low`、`ci_high` | 95% Wilson 置信区间 |
| `is_reference` | 是否与参考序列一致 |
| `n_snv`、`n_ins`、`n_del` | 变异数量 |
| `length` | 单倍型长度 |
| `variants` | 变异描述；`.` 表示无变异 |

### variants.tsv

方案 A 使用与公司一致的列布局：

```text
Chr  Pos  Ref  Alt  DP  Ref_dp  Alt_dp  Freq  DP4  Seq  Filter_Status  Filter_Reason
```

- `Freq` 是 0~1 的比例；
- `Ref` 或 `Alt` 中的 `-` 表示插入或缺失；
- `Filter_Status` 为 `PASS` 或 `FILTERED`。

### qc.tsv

两列：`metric` 和 `value`，包含 reads 数、比对率、平均一致性、覆盖度、聚类方法、
一致性方法和 DECIPHER 版本等。

## 命令行版本

R 包自带基于同一份 R 代码的 CLI，安装后以 `nanoamp` 命令提供，也可以在 R 中调用
`nanoamp_cli()`。

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

批量样本表是至少包含以下列的 TSV：

```text
sample	reads	reference
```

可选列：`ref_label`。

没有 minimap2 的平台可以用 `--aligner r` 选择 R 内比对后端。

### 安装 `nanoamp` 命令

```bash
sh "$(Rscript --vanilla -e 'cat(system.file("scripts", "install_cli.sh", package = "nanoamp"))')" ~/.local/bin
export PATH="$HOME/.local/bin:$PATH"
nanoamp doctor
```

也可以直接在 R 中调用 CLI：

```r
library(nanoamp)
nanoamp_cli(c("call", "--reads", "sample.fastq", "--reference", "target.fa",
              "--outdir", "results/sampleA"))
```

## 图形界面

启动 Shiny 界面：

```r
library(nanoamp)
nanoamp_gui()
```

界面提供文件选择、模式选择、高级参数、运行按钮、日志窗口、单倍型/变异交互表格、
QC 结果和下载按钮。

Windows 下也可以用自带启动器：

```bat
Rscript -e "library(nanoamp); nanoamp_gui()"
```

如果只需要 Shiny app 对象而不启动服务器：

```r
app <- nanoamp_gui_app()
```

### 外部工具与 R 内后端

`aligner = "minimap2"`（默认）需要一个 `minimap2` 可执行文件。没有 minimap2 的
机器可以改用 R 内后端：

```r
run_haplotype_analysis(..., aligner = "r")
```

R 内后端使用 Biostrings 成对比对，不需要外部工具；速度较慢，适合中小扩增子。

`samtools` 不是必需依赖：默认用 `Rsamtools::asBam()` 完成 SAM→BAM。
只有显式设置 `use_samtools = TRUE` 才会走 samtools。

## 测试

```r
# 单元测试
devtools::test()

# 完整包检查
devtools::check()
```

本包已通过 `R CMD check`，当前状态为 `Status: OK`。

## 常见问题

| 现象 | 解决方式 |
|---|---|
| 找不到 `minimap2` | 安装 minimap2 并加入 `PATH`，或设置 `NANOAMP_MINIMAP2` |
| 找不到 `samtools` | 通常不需要：默认用 `Rsamtools`；只有 `use_samtools = TRUE` 才需要 |
| 方案 B 太慢 | 减小 `max_msa_seqs`、增加 `threads`，或改用 `mode = "A"` |
| 方案 B 分不开相近单倍型 | 差异低于测序错误率时是正常现象，建议改用方案 A |
| 没有安装 `DECIPHER` | 方案 B 会退化为贪心聚类；建议安装 DECIPHER |
| 报错 `pairwiseAlignment` is not an exported object from Biostrings | Bioconductor >= 3.19 已把它移到 `pwalign`，用 `BiocManager::install("pwalign")` 安装 |
| 方案 C 的比例都很低 | 纳米孔 reads 含错误，建议改用方案 A |

## 许可证

MIT。
