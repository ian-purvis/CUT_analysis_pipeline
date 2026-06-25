# CUT&RUN Pipeline — Tulloch et al. eLife 2025

Pipeline reproducing the CUT&RUN analysis from **Tulloch et al. eLife 2025** ([doi:10.7554/eLife.107565](https://doi.org/10.7554/eLife.107565)), with alignment to **mm39 (GRCm39)** instead of mm10. This pipeline is compatible with Windows Subsystem for Linux (WSL).

## Data Source

- **Harvard Dataverse**: [doi:10.7910/DVN/TW0ZQL](https://doi.org/10.7910/DVN/TW0ZQL)
- Use `download_cutrun_from_metadata.ps1` in `E14 retinas/` to download FASTQ files

## Pipeline Steps

| Step | Tool | Version | Description |
|------|------|---------|-------------|
| QC | FASTQC | 0.11.9 | Quality control (includes repair report in HTML) |
| Repair | bbmap repair.sh | - | Fix paired FASTQs with unequal read counts |
| Trim | Trim Galore | 0.6.6 | Adapter trimming |
| Align | BWA mem | 0.7.17 | Alignment to **mm39** |
| Dedup | Samtools | 1.15.1 | Remove duplicates |
| Peaks | MACS3 | 3.x | Peak calling (p=0.05, ext=200, no lambda) |
| Coverage | DeepTools | 3.5 | RPKM BigWig tracks |

## Requirements

- **WSL (Ubuntu)** on Windows — BWA, samtools, MACS3 are Linux-only in bioconda
- Or native Linux/macOS

## Setup

### 0. Install Miniconda (if needed)

```bash
cd /path/to/cutrun_pipeline
bash Miniconda3-latest-Linux-x86_64.sh
# Follow prompts, then: source ~/.bashrc
```

### 1. Create conda environment

```bash
cd /path/to/cutrun_pipeline
conda env create -f environment.yaml
conda activate cutrun
```

Or run `bash setup_wsl.sh` (checks for conda first). If you already have the cutrun env and get `repair.sh: command not found`, run `conda install -c bioconda bbmap` to add bbmap.

### 2. Set up mm39 reference genome

Download or copy the mm39 (GRCm39) genome, then:

```bash
# Edit config.yaml to set genome_bwa path, then:
bash index_mm39.sh
```

Store the genome outside cloud-synced folders (e.g. `C:\genomes\mouse` or `~/genomes/mouse`) for reliable I/O.

### 3. Configure paths

Copy the example config and edit for your setup:

```bash
cp config.yaml.example config.yaml
```

Edit `config.yaml` for genome and samples:

- `genome_fa` / `genome_bwa`: absolute path to mm39 FASTA
- `blacklist_bed`: optional; use `""` to skip
- `control_sample`: IgG control for MACS3 peak calling
- `samples`: optional fallback; by default, samples are **auto-discovered** from `fastq_dir` (subdirs with `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`)
- `bbmap_xmx`: Java heap for repair.sh (e.g. `"4g"`). Empty = auto-detect. Set if you get OutOfMemoryError.
- `results_dir`: Optional. Override where outputs are written.

**`fastq_dir`** is specified at runtime when running the pipeline (see Run Pipeline below). Results are written to `results/` in the **parent** of `fastq_dir`, keeping each dataset's outputs separate (unless `results_dir` overrides).

Corrupted/incomplete FASTQs are auto-detected and skipped.

## Converting Windows paths to WSL

Windows `C:\` → `/mnt/c/` in WSL. Use `wslpath -u "C:\path"` to convert.

Paths with `&` (e.g. `CUT&RUN`) must be quoted in bash.

**Path handling:** The pipeline normalizes paths (strips trailing slashes, avoids double slashes) so that re-runs do not unnecessarily re-execute completed jobs. Config paths (`fastq_dir`, `genome_bwa`, etc.) are normalized automatically.

## Run Pipeline

**Required:** Pass `fastq_dir` (path to FASTQ folder) as the first argument. Results go to `results/` in that folder's parent.

### OneDrive / cloud-synced projects

When the pipeline and data live on OneDrive (or similar), use **`run_from_local.sh`** so Snakemake metadata (`.snakemake/`) is stored on the local WSL filesystem. This reduces I/O errors and random WSL exits. FASTQ and results stay on OneDrive.

```bash
# Recommended for OneDrive: metadata in ~/cutrun_work, data on OneDrive
./run_from_local.sh "../../Mu CUT&Tag_OTX2_2023/FASTQ" 4 all

# Single module: ./run_from_local.sh "<fastq_dir>" 4 all_qc
```

### Standard run (non-cloud or if run_from_local fails)

```bash
# Modular (recommended)
./run_qc.sh ../../Mu\ CUT&Tag_OTX2_2023/FASTQ 4
./run_align.sh ../../Mu\ CUT&Tag_OTX2_2023/FASTQ 4
./run_dedup.sh ../../Mu\ CUT&Tag_OTX2_2023/FASTQ 4
./run_peaks.sh ../../Mu\ CUT&Tag_OTX2_2023/FASTQ 4
./run_bigwig.sh ../../Mu\ CUT&Tag_OTX2_2023/FASTQ 4

# Or full pipeline
./run_pipeline.sh ../../Mu\ CUT&Tag_OTX2_2023/FASTQ 4
```

For paths with spaces or `&`, use quotes: `./run_qc.sh "../../Mu CUT&Tag_OTX2_2023/FASTQ" 4`

**Note:** Trailing slashes on `fastq_dir` are automatically removed to avoid path inconsistencies and unnecessary re-runs.

Each module validates inputs and runs upstream steps if needed. Use `./unlock.sh` if you get lock errors after interruption.

## Outputs

Outputs are written to `results/` in the parent of `fastq_dir` (keeps datasets separate).

| Output | Location |
|--------|----------|
| FASTQC | `results/fastqc/` (HTML reports include Paired Read Repair section: reads before/after, orphaned reads removed) |
| Trimmed FASTQ | `results/trimmed/` |
| Deduplicated BAM | `results/bam/*.dedup.bam` |
| MACS3 peaks | `results/peaks/*_peaks.narrowPeak` |
| RPKM BigWig | `results/bigwig/*.rpkm.bw` |

## Samples

- **TF samples**: FoxN4-1/2/3, Mybl1-1/2/3, Pax6-1/2/3, H3K4me3
- **Control**: IgG-2 (for MACS3 peak calling)

## Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| `repair.sh` exit 127 (command not found) | bbmap not installed or env created before bbmap was added | See [BBMap installation](#bbmap-installation) below. |
| `repair.sh` "Unknown parameter CUT" or path split at `&` | Path contains `&` (e.g. `Mu CUT&Tag_OTX2_2023`) | Pipeline now copies to `/tmp` before repair; re-run. |
| `inject_repair_report.py: No such file or directory` or `can't find '__main__' module in '...OneDrive'` | Path contains `&` (e.g. `CUT&RUN_Tulloch_2025`) | Fixed: script path is now quoted; re-run. |
| `repair.sh` OutOfMemoryError | bbmap memory autodetect fails (e.g. WSL) | Set `bbmap_xmx: "4g"` in config.yaml to override auto-detect; re-run. |
| `mv: cannot stat '..._fastqc.html': No such file or directory` | FastQC output names differ from expected | Fixed in pipeline: outputs are `{sample}_R1_fastqc.html` (not `{sample}_R1.fastq.gz_fastqc.html`). |
| `Input file names supplied with whitespace(s)` (Trim Galore) | Paths with spaces/`&` (e.g. `Mu CUT&Tag_OTX2_2023`) | Pipeline now copies to `/tmp` with `cp -L` before running. Re-run. |
| `Read 1 output is truncated... please check your paired-end input files` | R1 and R2 have different read counts | Pipeline runs bbmap repair.sh before Trim Galore; orphaned reads are removed. Repair stats appear in the FastQC HTML report. |
| `mv: error writing '...': Input/output error` | Cross-filesystem move (e.g. /tmp → OneDrive) failing | Pipeline now uses temp dirs *under* results (same filesystem as OneDrive), so final `mv` is a rename. Re-run. |
| `[Errno 28] No space left on device` | Disk full | Free space on C: (and OneDrive cache if applicable). |
| `pigz: abort: write error on <stdout>` | I/O error (OneDrive sync, WSL, disk) | Run from a non–OneDrive path, or ensure files are fully synced. |
| `cp: error reading '...': Input/output error` | File locked, corrupt, or OneDrive sync issue | Verify file: `gzip -t file.fastq.gz`. Re-download or exclude sample. |
| WSL exits unexpectedly | Memory pressure, OneDrive, or sleep | Use `./run_from_local.sh` so metadata is local; add `.wslconfig` with `memory=8GB`, run `wsl --shutdown`. |
| `Failed to set marker file... Errno 22 Invalid argument` | Path too long (OneDrive + nested dirs) | Use `./run_from_local.sh` so metadata lives in `~/cutrun_work`. Pipeline still runs. |
| Jobs re-run unnecessarily / "Set of input files has changed" | Path format changed (trailing slash, double slash) | Pipeline now normalizes paths. Use consistent `fastq_dir` (with or without trailing slash—both work). |

### BBMap installation

BBMap provides `repair.sh` for fixing paired-end FASTQs with unequal read counts. If `conda install -c bioconda bbmap` fails (solver conflicts, channel issues), try:

1. **Update environment from file:**
   ```bash
   conda activate cutrun
   conda env update -f environment.yaml --prune
   ```

2. **Use mamba** (often resolves conflicts better):
   ```bash
   conda install -c conda-forge mamba
   mamba install -c bioconda bbmap
   ```

3. **Explicit channel order:**
   ```bash
   conda install -c bioconda -c conda-forge -c defaults bbmap
   ```

4. **Manual install** (if conda still fails): Download [BBMap](https://sourceforge.net/projects/bbmap/), extract, and add the `bbmap/` directory to your PATH so `repair.sh` is found when the cutrun env is activated.

## Documentation

A full pipeline reference is in `PIPELINE_DOCUMENTATION.md` (Markdown). To generate a Word document:

```bash
pip install python-docx
python generate_pipeline_doc.py
```

This creates `PIPELINE_DOCUMENTATION.docx`. Open in Word and use **File → Save As → PDF** for a PDF.

## Reference

Tulloch AJ, Delgado RN, Catta-Preta R, Cepko CL (2025). *Massively parallel reporter assay for mapping gene-specific regulatory regions at single nucleotide resolution.* eLife. doi:10.7554/eLife.107565
