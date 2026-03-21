# CUT&RUN Pipeline — Complete Documentation

**Version:** 1.0  
**Reference:** Tulloch et al. eLife 2025 (doi:10.7554/eLife.107565)  
**Pipeline:** Replicates CUT&RUN analysis with alignment to mm39 (GRCm39)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Directory Structure & File Relationships](#2-directory-structure--file-relationships)
3. [Configuration Reference](#3-configuration-reference)
4. [Pipeline Modules (Step-by-Step)](#4-pipeline-modules-step-by-step)
5. [Snakefile Rules & Functions](#5-snakefile-rules--functions)
6. [Python Scripts Reference](#6-python-scripts-reference)
7. [Shell Scripts Reference](#7-shell-scripts-reference)
8. [Dependencies](#8-dependencies)
9. [Flexibility & Customization](#9-flexibility--customization)
10. [Stringency & Validation](#10-stringency--validation)
11. [Data Flow Diagram](#11-data-flow-diagram)
12. [Troubleshooting](#12-troubleshooting)

---

## 1. Overview

### Purpose

This pipeline processes CUT&RUN (and CUT&Tag) sequencing data from raw paired-end FASTQ files through quality control, alignment, deduplication, peak calling, and coverage visualization. It is modular: each stage can be run independently or as part of a full run.

### Workflow Summary

| Stage | Tool | Input | Output |
|-------|------|-------|--------|
| 1. QC | repair.sh, FastQC, Trim Galore | Raw FASTQ (R1, R2) | FastQC reports (with repair stats), trimmed FASTQ |
| 2. Align | BWA mem | Trimmed FASTQ | Sorted BAM |
| 3. Dedup | samtools markdup | BAM | Deduplicated BAM + index |
| 4. Peaks | MACS3 | Dedup BAM (treatment + control) | narrowPeak, summits BED |
| 5. BigWig | DeepTools bamCoverage | Dedup BAM | RPKM-normalized BigWig |

### Key Design Decisions

- **Results location:** Outputs go to `results/` in the **parent** of `fastq_dir`, so each dataset keeps its own results.
- **Sample discovery:** Samples are auto-discovered from subdirectories in `fastq_dir`; no hardcoded sample list.
- **Path handling:** Tools run on files copied to `/tmp` to avoid issues with paths containing spaces or `&`.
- **Validation:** Each module validates inputs before running; corrupt FASTQs are skipped.

---

## 2. Directory Structure & File Relationships

### Expected Input Structure

```
fastq_dir/
├── SampleA/
│   ├── SampleA_R1.fastq.gz
│   └── SampleA_R2.fastq.gz
├── SampleB/
│   ├── SampleB_R1.fastq.gz
│   └── SampleB_R2.fastq.gz
└── ...
```

**Naming convention:** `{sample_name}/{sample_name}_R1.fastq.gz` and `{sample_name}_R2.fastq.gz`

### Output Structure (results/)

```
results/                          # Parent of fastq_dir + /results
├── .setup_done                   # Setup completion marker
├── .check_qc_ok                  # QC inputs validated
├── .samples_discovered.txt       # Auto-discovered sample list
├── .check_align_ok               # Align inputs validated
├── .check_dedup_ok               # Dedup inputs validated
├── .check_peaks_ok               # Peaks inputs validated
├── .check_bigwig_ok              # BigWig inputs validated
├── validated/
│   └── samples.txt               # Valid samples (after gzip/format check)
├── fastq_inputs/
│   └── {sample}/
│       ├── {sample}_R1.fastq.gz  # Symlink to original
│       └── {sample}_R2.fastq.gz
├── fastq_ready/
│   └── {sample}/
│       ├── {sample}_R1.fastq.gz  # Repaired (paired) FASTQ
│       └── {sample}_R2.fastq.gz
├── fastqc/
│   ├── {sample}_R1_fastqc.html   # Includes Paired Read Repair section
│   ├── {sample}_R1_fastqc.zip
│   ├── {sample}_R2_fastqc.html  # Includes Paired Read Repair section
│   └── {sample}_R2_fastqc.zip
├── trimmed/
│   ├── {sample}_R1_val_1.fq.gz
│   ├── {sample}_R2_val_2.fq.gz
│   ├── {sample}_R1.fastq.gz_trimming_report.txt
│   └── {sample}_R2.fastq.gz_trimming_report.txt
├── bam/
│   ├── {sample}.bam
│   ├── {sample}.dedup.bam
│   └── {sample}.dedup.bam.bai
├── peaks/
│   ├── {sample}_peaks.narrowPeak
│   └── {sample}_summits.bed
└── bigwig/
    └── {sample}.rpkm.bw
```

### Pipeline Directory (cutrun_pipeline/)

```
cutrun_pipeline/
├── Snakefile              # Main workflow definition
├── config.yaml            # Configuration (paths, parameters)
├── environment.yaml       # Conda environment
├── run_pipeline.sh        # Full pipeline
├── run_qc.sh              # QC module only
├── run_align.sh           # Align module only
├── run_dedup.sh           # Dedup module only
├── run_peaks.sh           # Peaks module only
├── run_bigwig.sh          # BigWig module only
├── unlock.sh              # Clear Snakemake lock
├── index_mm39.sh          # Build BWA index
├── setup_wsl.sh           # WSL setup
├── scripts/
│   ├── check_inputs.py         # Input validation per module
│   ├── validate_fastq.py       # FASTQ integrity + format check
│   └── inject_repair_report.py # Inject repair stats into FastQC HTML
├── logs/
│   ├── repair/
│   ├── fastqc/
│   ├── trim_galore/
│   ├── bwa/
│   ├── dedup/
│   ├── macs3/
│   └── deeptools/
└── reference/             # Optional reference files
```

### File Dependency Graph

```
fastq_dir/{sample}/{sample}_R1.fastq.gz, _R2.fastq.gz
    │
    ├─► check_qc_inputs ─► .samples_discovered.txt
    │
    ├─► validate ─► validated/samples.txt
    │
    ├─► link_fastq ─► fastq_inputs/{sample}/*.fastq.gz (symlinks)
    │
    └─► repair_paired ─► fastq_ready/{sample}/*.fastq.gz (paired, orphaned removed)
            │
            ├─► fastqc ─► fastqc/*.html (with repair report), *.zip
            │
            └─► trim_galore ─► trimmed/*_val_1.fq.gz, *_val_2.fq.gz
            │
            └─► bwa_align ─► bam/{sample}.bam
                    │
                    └─► mark_duplicates ─► bam/{sample}.dedup.bam
                            │
                            ├─► macs3_callpeak ─► peaks/*_peaks.narrowPeak (TF samples only)
                            │
                            └─► bam_coverage ─► bigwig/*.rpkm.bw
```

---

## 3. Configuration Reference

### config.yaml — All Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|--------------|
| `fastq_dir` | string | (required at runtime) | Path to FASTQ folder. Relative to `cutrun_pipeline/` or absolute. Overridden by run scripts. |
| `genome_fa` | string | — | Path to mm39 FASTA (for samtools faidx). |
| `genome_bwa` | string | — | Path to mm39 FASTA for BWA (same as genome_fa). Must have BWA index (.bwt, .ann, etc.). |
| `blacklist_bed` | string | `""` | Optional BED file of regions to exclude from BigWig. Empty = no blacklist. |
| `control_sample` | string | `"IgG-2"` | Sample name used as control for MACS3 peak calling (e.g. IgG). |
| `bbmap_xmx` | string | `""` | bbmap Java heap. Empty = auto-detect (50% of available RAM, min 1g, max 32g). Override if needed: `"4g"`. |
| `samples` | dict | (optional) | Fallback sample list if auto-discovery finds nothing. Format: `SampleName: ["SampleName"]`. |
| `macs3_pvalue` | float | `0.05` | MACS3 peak p-value threshold. |
| `macs3_extsize` | int | `200` | MACS3 extension size (bp). |
| `macs3_nolambda` | bool | `true` | Disable MACS3 local lambda. |
| `threads` | int | `4` | Default threads for parallel steps. |

### Runtime Override

All run scripts pass `fastq_dir` to Snakemake:

```bash
snakemake all_qc --config fastq_dir="/path/to/FASTQ" -j 4
```

This overrides `config.yaml` for that run.

---

## 4. Pipeline Modules (Step-by-Step)

### Module 1: QC (all_qc)

**Target:** `all_qc`

**Steps (order):**
1. **setup** — Create output directories.
2. **check_qc_inputs** — Discover samples, validate R1/R2 exist, write `.samples_discovered.txt`.
3. **validate** (checkpoint) — Run gzip integrity + FASTQ format check; write `validated/samples.txt` with valid samples only.
4. **link_fastq** — Create symlinks in `fastq_inputs/` for each valid sample.
5. **repair_paired** — bbmap repair.sh: fix paired FASTQs with unequal read counts; output to `fastq_ready/`.
6. **fastqc** — Quality control on repaired R1/R2; injects repair stats (reads before/after, orphaned removed) into HTML reports.
7. **trim_galore** — Adapter and quality trimming.

**Tool arguments (repair.sh):**
- `in1` / `in2` — Input R1/R2 from `fastq_inputs/`.
- `out1` / `out2` — Output to `fastq_ready/`. Orphaned reads (no matching mate) are excluded.

**Tool arguments (Trim Galore):**
- `--paired` — Paired-end mode.
- `--output_dir` — Output directory (tmp, then moved).
- `--cores` — Threads from config.
- Adapter: auto-detected (Nextera, Illumina, etc.).
- Quality: Phred 20 cutoff (default).

**Tool arguments (FastQC):**
- `-o` — Output directory.
- `-t` — Threads (2 per sample).

**Repair report injection:** After FastQC, `inject_repair_report.py` adds a "Paired Read Repair" section to each HTML report, showing reads before/after repair, orphaned R1/R2 removed, and total reads removed. This documents the repair step for each sample.

### Module 2: Alignment (all_align)

**Target:** `all_align`

**Steps:**
1. **check_align_inputs** — Validate genome, BWA index, trimmed FASTQs.
2. **bwa_align** — BWA mem alignment, samtools sort.

**Tool arguments (BWA):**
- `bwa mem -t {threads}` — Align with specified threads.
- `samtools sort -@ {threads}` — Sort by coordinate.

### Module 3: Deduplication (all_dedup)

**Target:** `all_dedup`

**Steps:**
1. **check_dedup_inputs** — Validate BAM files exist.
2. **mark_duplicates** — samtools markdup with fixmate.

**Tool arguments (samtools):**
- `samtools sort -n` — Sort by read name (required for fixmate).
- `samtools fixmate -m` — Add mate score tags.
- `samtools sort` — Sort by coordinate.
- `samtools markdup -r -s` — Remove duplicates (`-r`), use secondary alignments (`-s`).
- `samtools index` — Create BAI index.

### Module 4: Peak Calling (all_peaks)

**Target:** `all_peaks`

**Steps:**
1. **check_peaks_inputs** — Validate control + treatment dedup BAMs.
2. **macs3_callpeak** — MACS3 peak calling (TF samples only; control excluded from treatment list).

**Tool arguments (MACS3):**
- `-t` — Treatment BAM.
- `-c` — Control BAM.
- `-n` — Output name prefix.
- `-p` — P-value (from config).
- `--extsize` — Extension size (from config).
- `--nolambda` — If `macs3_nolambda: true`.
- `--format BAM` — Input format.

**Outputs:** `*_peaks.narrowPeak`, `*_summits.bed`

### Module 5: BigWig (all_bigwig)

**Target:** `all_bigwig`

**Steps:**
1. **check_bigwig_inputs** — Validate dedup BAMs exist.
2. **bam_coverage** — DeepTools bamCoverage.

**Tool arguments (bamCoverage):**
- `--normalizeUsing RPKM` — RPKM normalization.
- `--binSize 10` — 10 bp bins.
- `--blackListFileName` — If blacklist BED is set.

---

## 5. Snakefile Rules & Functions

### Global Variables

| Variable | Source | Purpose |
|----------|--------|---------|
| `PIPELINE_DIR` | `config.yaml` location | Pipeline root. |
| `RESDIR` | `parent(fastq_dir)/results` | Results directory (relative). |
| `CONTROL` | `config["control_sample"]` | Control sample name. |
| `BLACKLIST` | `config["blacklist_bed"]` | Blacklist path. |

### Functions

| Function | Purpose |
|----------|---------|
| `get_valid_samples()` | Read `validated/samples.txt` (used when checkpoint not yet run). |
| `_get_valid_samples()` | Read from validate checkpoint output (used after checkpoint). |
| `_aggregate_qc_input(wildcards)` | Aggregate inputs for `all_qc`. |
| `_aggregate_align_input(wildcards)` | Aggregate inputs for `check_align_inputs`. |
| `_aggregate_align_target(wildcards)` | Aggregate targets for `all_align`. |
| `_aggregate_dedup_input(wildcards)` | Aggregate inputs for `check_dedup_inputs`. |
| `_aggregate_dedup_target(wildcards)` | Aggregate targets for `all_dedup`. |
| `_aggregate_peaks_input(wildcards)` | Aggregate inputs for `check_peaks_inputs`. |
| `_aggregate_peaks_target(wildcards)` | Aggregate targets for `all_peaks`. |
| `_aggregate_bigwig_input(wildcards)` | Aggregate inputs for `check_bigwig_inputs`. |
| `_aggregate_bigwig_target(wildcards)` | Aggregate targets for `all_bigwig`. |
| `_aggregate_all_input(wildcards)` | Aggregate inputs for full pipeline. |

### ruleorder

```python
ruleorder: mark_duplicates > bwa_align
```

Resolves ambiguity when a path like `IgG-2.dedup.bam` could match either `bwa_align` (sample=IgG-2.dedup) or `mark_duplicates` (sample=IgG-2). Prefer `mark_duplicates`.

---

## 6. Python Scripts Reference

### scripts/check_inputs.py

**Purpose:** Validate inputs for each module. Exits with code 1 on failure.

**Usage:**
```bash
check_inputs.py <module> [args...]
```

| Module | Args | Validates |
|--------|------|-----------|
| `qc` | `fastq_dir config.yaml [samples_out]` | FASTQ dir exists, samples have R1/R2; discovers samples; writes `samples_out`. |
| `align` | `genome_bwa samples.txt trimmed_dir` | Genome + BWA index exist; trimmed FASTQs exist for all samples. |
| `dedup` | `samples.txt bam_dir` | BAM files exist for all samples. |
| `peaks` | `samples.txt bam_dir control` | Control BAM + treatment dedup BAMs exist. |
| `bigwig` | `samples.txt bam_dir` | Dedup BAMs exist for all samples. |

**Functions:**
- `get_samples_from_file(samples_file)` — Read sample list from file.
- `get_samples_from_config(config_path)` — Read samples from config (fallback).
- `discover_samples(fastq_dir)` — Find subdirs with `{name}_R1.fastq.gz` and `{name}_R2.fastq.gz`.
- `check_qc(...)` — QC validation + discovery.
- `check_align(...)` — Align validation.
- `check_dedup(...)` — Dedup validation.
- `check_peaks(...)` — Peaks validation.
- `check_bigwig(...)` — BigWig validation.

### scripts/inject_repair_report.py

**Purpose:** Inject paired-read repair statistics into FastQC HTML reports. Adds a "Paired Read Repair" section at the top of each report showing reads before/after repair, orphaned R1/R2 removed, and total reads removed.

**Usage:**
```bash
inject_repair_report.py --sample <name> --input-r1 <path> --input-r2 <path> \
  --output-r1 <path> --output-r2 <path> --repair-log <path> \
  --fastqc-html <R1_html> [--fastqc-html-r2 <R2_html>]
```

**Behavior:** Counts reads in gzipped FASTQs (before and after repair), parses the repair log for orphaned counts, builds an HTML module, and injects it before the first FastQC module in both R1 and R2 reports. Called automatically by the `fastqc` rule after FastQC completes.

### scripts/validate_fastq.py

**Purpose:** Check FASTQ integrity (gzip CRC, format). Writes valid sample names only. Always exits 0.

**Usage:**
```bash
validate_fastq.py <fastq_dir> <sample1> <sample2> ... --output <out.txt>
```

**Functions:**
- `gzip_integrity_ok(path)` — Run `gzip -t`; returns False on CRC/truncation.
- `validate_fastq(path, max_records=1800000)` — Parse first ~1.8M records; check `@` and `+` lines.

**Behavior:** Invalid samples are skipped; only valid samples are written to output.

---

## 7. Shell Scripts Reference

| Script | Purpose | Usage |
|--------|---------|-------|
| `run_pipeline.sh` | Run all modules in sequence | `./run_pipeline.sh <fastq_dir> [jobs]` |
| `run_qc.sh` | QC only | `./run_qc.sh <fastq_dir> [jobs]` |
| `run_align.sh` | Align only | `./run_align.sh <fastq_dir> [jobs]` |
| `run_dedup.sh` | Dedup only | `./run_dedup.sh <fastq_dir> [jobs]` |
| `run_peaks.sh` | Peaks only | `./run_peaks.sh <fastq_dir> [jobs]` |
| `run_bigwig.sh` | BigWig only | `./run_bigwig.sh <fastq_dir> [jobs]` |
| `unlock.sh` | Clear Snakemake lock | `./unlock.sh` |
| `index_mm39.sh` | Build BWA index | `./index_mm39.sh` (reads genome from config) |
| `setup_wsl.sh` | WSL setup | `./setup_wsl.sh` |

**Common Snakemake options:** `-j N` (jobs), `-k` (keep going on failure), `--latency-wait 60`, `--rerun-incomplete`.

---

## 8. Dependencies

### Conda (environment.yaml)

| Package | Purpose |
|---------|---------|
| python >= 3.9 | Scripts |
| snakemake >= 7.0 | Workflow |
| bbmap | Paired-read repair (repair.sh) |
| fastqc | QC |
| trim-galore | Trimming |
| bwa | Alignment |
| samtools | BAM handling, dedup |
| macs3 | Peak calling |
| deeptools | BigWig |

**Channels:** bioconda, conda-forge, defaults

**Note:** Bioconda tools (bwa, samtools, macs3) are Linux-only; use WSL on Windows.

### External Requirements

- **Reference genome:** mm39 (GRCm39) FASTA + BWA index.
- **BWA index:** Run `index_mm39.sh` or `bwa index <genome>`.
- **Optional:** Blacklist BED for BigWig.

---

## 9. Flexibility & Customization

### What You Can Change

| Item | How |
|------|-----|
| `fastq_dir` | Pass at runtime or edit config. |
| Genome path | Edit `genome_fa` / `genome_bwa` in config. |
| Control sample | Edit `control_sample` in config. |
| MACS3 p-value | Edit `macs3_pvalue` in config. |
| MACS3 extsize | Edit `macs3_extsize` in config. |
| MACS3 nolambda | Edit `macs3_nolambda` in config. |
| Threads | Edit `threads` in config or pass `-j N` to Snakemake. |
| Blacklist | Set `blacklist_bed` in config. |
| Sample list | Rely on auto-discovery or set `samples` in config as fallback. |

### What Is Fixed

- FASTQ naming: `{sample}_R1.fastq.gz`, `{sample}_R2.fastq.gz`.
- Directory layout: `fastq_dir/{sample}/`.
- Results location: `parent(fastq_dir)/results/`.
- Trim Galore: paired, retain_unpaired, auto-adapters, Phred 20.
- BWA: mem algorithm.
- Dedup: samtools markdup with fixmate.
- BigWig: RPKM, binSize 10.

---

## 10. Stringency & Validation

### FASTQ Validation

- **gzip -t:** CRC and truncation check.
- **Format:** First ~1.8M records checked for `@` (header) and `+` (separator).
- **Outcome:** Invalid samples skipped; only valid samples in `validated/samples.txt`.

### Per-Module Checks

- **QC:** FASTQ dir exists; discovered samples have R1 and R2.
- **Align:** Genome and BWA index exist; trimmed FASTQs exist for all valid samples.
- **Dedup:** BAM files exist for all valid samples.
- **Peaks:** Control BAM exists; treatment dedup BAMs exist; at least one TF sample.
- **BigWig:** Dedup BAMs exist for all valid samples.

### Path Handling

- Files copied to `/tmp` before running tools to avoid spaces/`&` in paths.
- `cp -L` used to dereference symlinks and copy real files.

---

## 11. Data Flow Diagram

```
                    ┌─────────────────────────────────────────┐
                    │           fastq_dir (user input)         │
                    │  {sample}/{sample}_R1.fastq.gz, _R2      │
                    └────────────────────┬────────────────────┘
                                         │
                    ┌────────────────────▼────────────────────┐
                    │  check_qc_inputs → .samples_discovered   │
                    │  validate → validated/samples.txt       │
                    │  link_fastq → fastq_inputs/ (symlinks)  │
                    │  repair_paired → fastq_ready/           │
                    └────────────────────┬────────────────────┘
                                         │
              ┌──────────────────────────┼──────────────────────────┐
              │                          │                          │
              ▼                          ▼                          │
    ┌─────────────────┐      ┌──────────────────┐                  │
    │     fastqc      │      │   trim_galore    │                  │
    │  fastqc/*.html   │      │ trimmed/*.fq.gz  │                  │
    │  (+repair rpt)  │      └────────┬─────────┘                  │
    └─────────────────┘                                               │
                                       │                            │
                                       ▼                            │
                            ┌──────────────────┐                    │
                            │    bwa_align     │                    │
                            │   bam/*.bam      │                    │
                            └────────┬─────────┘                    │
                                     │                              │
                                     ▼                              │
                          ┌──────────────────────┐                   │
                          │  mark_duplicates     │                   │
                          │  bam/*.dedup.bam    │                   │
                          └──────────┬──────────┘                   │
                                      │                              │
                    ┌─────────────────┼─────────────────┐            │
                    │                 │                 │            │
                    ▼                 ▼                 │            │
          ┌──────────────┐  ┌─────────────────┐        │            │
          │ macs3_callpeak│  │  bam_coverage   │        │            │
          │ peaks/*.narrowPeak│ bigwig/*.rpkm.bw │      │            │
          └──────────────┘  └─────────────────┘        │            │
           (TF samples only)  (all samples)             │            │
                                                        │            │
                    Control sample excluded from TF list │            │
```

---

## 12. Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| `mv: cannot stat '..._fastqc.html'` | FastQC output naming | Fixed: outputs are `{sample}_R1_fastqc.html`. |
| `Input file names supplied with whitespace(s)` | Paths with spaces/`&` | Pipeline copies to `/tmp`; re-run. |
| `repair.sh` exit 127 | bbmap not in PATH | Run `conda install -c bioconda bbmap`. If that fails: `conda env update -f environment.yaml --prune`, or `mamba install -c bioconda bbmap`, or install [BBMap](https://sourceforge.net/projects/bbmap/) manually and add to PATH. |
| `repair.sh` "Unknown parameter CUT" | Path contains `&` (e.g. `CUT&Tag`) | Pipeline copies to `/tmp` before repair; re-run. |
| `Read 1 output is truncated...` | R1/R2 read count mismatch | bbmap repair.sh runs before Trim Galore; orphaned reads removed. Repair stats in FastQC HTML. |
| `pigz: abort: write error` | I/O (OneDrive, WSL, disk) | Use non-OneDrive path; ensure sync complete. |
| `cp: error reading... Input/output error` | File lock/corruption | Run `gzip -t file.fastq.gz`; re-download or exclude. |
| WSL exits unexpectedly | Memory, OneDrive, sleep | Add `.wslconfig` with `memory=8GB`; `wsl --shutdown`. |
| `Directory cannot be locked` | Snakemake lock | Run `./unlock.sh`. |
| `No samples found` | Wrong fastq_dir or layout | Check `fastq_dir/{sample}/{sample}_R1.fastq.gz` exists. |
| `Control BAM not found` | Wrong control_sample | Set `control_sample` in config to match sample name. |
| `repair.sh` OutOfMemoryError | Java heap auto-detect fails (e.g. WSL) | Set `bbmap_xmx: "4g"` in config.yaml. |

### BBMap installation (when conda fails)

Conda may fail to install bbmap due to solver conflicts or channel order. Try in order:

1. `conda env update -f environment.yaml --prune`
2. `mamba install -c bioconda bbmap` (install mamba first: `conda install -c conda-forge mamba`)
3. `conda install -c bioconda -c conda-forge -c defaults bbmap`
4. Manual: download BBMap, extract, add `bbmap/` to PATH.

---

## Appendix: Quick Reference

### Run Commands

```bash
# Full pipeline
./run_pipeline.sh "../../Mu CUT&Tag_OTX2_2023/FASTQ" 4

# Individual modules
./run_qc.sh "../../Mu CUT&Tag_OTX2_2023/FASTQ" 4
./run_align.sh "../../Mu CUT&Tag_OTX2_2023/FASTQ" 4
./run_dedup.sh "../../Mu CUT&Tag_OTX2_2023/FASTQ" 4
./run_peaks.sh "../../Mu CUT&Tag_OTX2_2023/FASTQ" 4
./run_bigwig.sh "../../Mu CUT&Tag_OTX2_2023/FASTQ" 4
```

### Key Output Files

- **QC:** `results/fastqc/*.html` (includes Paired Read Repair section), `results/trimmed/*_val_1.fq.gz`
- **Align:** `results/bam/*.bam`
- **Dedup:** `results/bam/*.dedup.bam`
- **Peaks:** `results/peaks/*_peaks.narrowPeak`
- **BigWig:** `results/bigwig/*.rpkm.bw`

---

*Document generated for CUT&RUN pipeline. For questions, refer to the Snakefile, config.yaml, and scripts.*
