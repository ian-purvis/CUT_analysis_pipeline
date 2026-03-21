"""
CUT&RUN Processing Pipeline (Modular)
Replicates Tulloch et al. eLife 2025 (https://doi.org/10.7554/eLife.107565)

Temp dirs: All processing uses temp dirs under results (same filesystem as outputs
to avoid cross-fs move errors). Stale temp is cleaned at rule start; trap EXIT
removes on normal exit. Outputs go to target = parent(fastq_dir)/results.
"""

configfile: "config.yaml"
config["fastq_dir"] = (config.get("fastq_dir") or "").strip().rstrip("/") or "."

import os as _os
import subprocess
import sys

def _norm_path(p):
    if not p or not isinstance(p, str):
        return p or ""
    return _os.path.normpath(p.strip().rstrip("/")) if p.strip() else ""

config["genome_fa"] = _norm_path(config.get("genome_fa") or "")
config["genome_bwa"] = _norm_path(config.get("genome_bwa") or "")
config["blacklist_bed"] = _norm_path(config.get("blacklist_bed") or "")

PIPELINE_DIR = getattr(workflow, "basedir", None) or _os.path.dirname(_os.path.abspath("config.yaml"))
_fastq_abs = _os.path.normpath(_os.path.join(PIPELINE_DIR, config["fastq_dir"]))
_default_resdir = _os.path.normpath(_os.path.relpath(_os.path.join(_os.path.dirname(_fastq_abs), "results"), PIPELINE_DIR))
_override = _norm_path(config.get("results_dir") or "")
RESDIR = _os.path.normpath(_override) if _override and _os.path.isabs(_override) else (_os.path.normpath(_override) if _override else _default_resdir)
CONTROL = config.get("control_sample", "IgG-2")
BLACKLIST = config.get("blacklist_bed", "")
ruleorder: mark_duplicates > bwa_align

def get_bbmap_xmx():
    xmx = config.get("bbmap_xmx", "").strip() if config.get("bbmap_xmx") else ""
    if xmx:
        return "-Xmx" + xmx if not str(xmx).startswith("-Xmx") else xmx
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    kb = int(line.split()[1])
                    mb = max(1024, min(32768, kb * 50 // 100 // 1024))
                    # Round to nearest 2GB to avoid param churn and unnecessary re-runs
                    gb = max(1, round(mb / 1024))
                    return f"-Xmx{gb}g"
    except (OSError, ValueError):
        pass
    return "-Xmx4g"

def get_valid_samples():
    path = _os.path.join(getattr(workflow, "basedir", None) or PIPELINE_DIR, RESDIR, "validated", "samples.txt")
    if not _os.path.exists(path):
        return []
    with open(path) as f:
        return [s.strip() for s in f if s.strip()]

rule setup:
    output: touch(RESDIR + "/.setup_done")
    params: resdir = RESDIR,
    run:
        shell("mkdir -p '{params.resdir}/validated' '{params.resdir}/fastq_inputs' '{params.resdir}/fastq_ready' '{params.resdir}/fastqc' '{params.resdir}/trimmed' '{params.resdir}/bam' '{params.resdir}/peaks' '{params.resdir}/bigwig' logs/fastqc logs/trim_galore logs/repair logs/bwa logs/dedup logs/macs3 logs/deeptools reference; touch {output}")

rule check_qc_inputs:
    input: _setup = RESDIR + "/.setup_done",
    output: touch(RESDIR + "/.check_qc_ok"), samples = RESDIR + "/.samples_discovered.txt",
    run:
        cmd = [sys.executable, "scripts/check_inputs.py", "qc", config["fastq_dir"], "config.yaml", output.samples]
        subprocess.run(cmd, check=True, cwd=_os.path.abspath("."))

checkpoint validate:
    input: _check = RESDIR + "/.check_qc_ok", samples = RESDIR + "/.samples_discovered.txt",
    output: samples = RESDIR + "/validated/samples.txt",
    run:
        with open(input.samples) as f:
            slist = [s.strip() for s in f if s.strip()]
        if slist:
            fastq_dir = _os.path.normpath(_os.path.join(_os.path.abspath("."), config["fastq_dir"]))
            subprocess.run([sys.executable, "scripts/validate_fastq.py", fastq_dir] + slist + ["--output", output.samples], check=True, cwd=_os.path.abspath("."))
        else:
            shell("touch {output.samples}")

rule link_fastq:
    input:
        _setup = RESDIR + "/.setup_done",
        _valid = RESDIR + "/validated/samples.txt",
        r1 = lambda w: _os.path.join(config["fastq_dir"], w.sample, f"{w.sample}_R1.fastq.gz"),
        r2 = lambda w: _os.path.join(config["fastq_dir"], w.sample, f"{w.sample}_R2.fastq.gz"),
    output: r1 = RESDIR + "/fastq_inputs/{sample}/{sample}_R1.fastq.gz", r2 = RESDIR + "/fastq_inputs/{sample}/{sample}_R2.fastq.gz",
    run:
        workdir = _os.path.abspath(".")
        outdir = _os.path.join(workdir, RESDIR, "fastq_inputs", wildcards.sample)
        _os.makedirs(outdir, exist_ok=True)
        rel = _os.path.relpath(_os.path.join(workdir, config["fastq_dir"], wildcards.sample), outdir)
        for dst, fn in [(_os.path.join(outdir, f"{wildcards.sample}_R1.fastq.gz"), f"{wildcards.sample}_R1.fastq.gz"), (_os.path.join(outdir, f"{wildcards.sample}_R2.fastq.gz"), f"{wildcards.sample}_R2.fastq.gz")]:
            if _os.path.lexists(dst): _os.remove(dst)
            _os.symlink(_os.path.join(rel, fn), dst)

rule repair_paired:
    input: r1 = RESDIR + "/fastq_inputs/{sample}/{sample}_R1.fastq.gz", r2 = RESDIR + "/fastq_inputs/{sample}/{sample}_R2.fastq.gz",
    output: r1 = RESDIR + "/fastq_ready/{sample}/{sample}_R1.fastq.gz", r2 = RESDIR + "/fastq_ready/{sample}/{sample}_R2.fastq.gz",
    params: resdir = RESDIR, xmx = get_bbmap_xmx(), tmpdir = RESDIR + "/.tmp_repair_{sample}",
    log: "logs/repair/{sample}.log"
    shell:
        "TMP='{params.tmpdir}'; rm -rf \"$TMP\"; trap 'rm -rf \"$TMP\"' EXIT; "
        "mkdir -p \"$TMP\" '{params.resdir}/fastq_ready/{wildcards.sample}'; "
        "cp -L '{input.r1}' '{input.r2}' \"$TMP\"/; "
        "repair.sh {params.xmx} in1=\"$TMP/{wildcards.sample}_R1.fastq.gz\" in2=\"$TMP/{wildcards.sample}_R2.fastq.gz\" out1=\"$TMP/R1_repaired.fastq.gz\" out2=\"$TMP/R2_repaired.fastq.gz\" >> {log} 2>&1; "
        "mv \"$TMP/R1_repaired.fastq.gz\" '{output.r1}' && mv \"$TMP/R2_repaired.fastq.gz\" '{output.r2}' && rm -rf \"$TMP\""

rule fastqc:
    input:
        r1 = RESDIR + "/fastq_ready/{sample}/{sample}_R1.fastq.gz", r2 = RESDIR + "/fastq_ready/{sample}/{sample}_R2.fastq.gz",
        r1_input = RESDIR + "/fastq_inputs/{sample}/{sample}_R1.fastq.gz", r2_input = RESDIR + "/fastq_inputs/{sample}/{sample}_R2.fastq.gz",
        repair_log = "logs/repair/{sample}.log",
    output: html1 = RESDIR + "/fastqc/{sample}_R1_fastqc.html", html2 = RESDIR + "/fastqc/{sample}_R2_fastqc.html", zip1 = RESDIR + "/fastqc/{sample}_R1_fastqc.zip", zip2 = RESDIR + "/fastqc/{sample}_R2_fastqc.zip",
    params: resdir = RESDIR, tmpdir = RESDIR + "/.tmp_fastqc_{sample}",
    threads: 2
    log: "logs/fastqc/{sample}.log"
    run:
        workdir = _os.path.abspath(".")
        script = _os.path.join(workdir, "scripts", "inject_repair_report.py")
        shell(
            "TMP='{params.tmpdir}'; rm -rf \"$TMP\"; trap 'rm -rf \"$TMP\"' EXIT; "
            "mkdir -p \"$TMP\" '{params.resdir}/fastqc'; cp -L '{input.r1}' '{input.r2}' \"$TMP\"/; "
            "fastqc -o \"$TMP\" -t {threads} \"$TMP/{wildcards.sample}_R1.fastq.gz\" \"$TMP/{wildcards.sample}_R2.fastq.gz\" >> {log} 2>&1; "
            "mv \"$TMP/{wildcards.sample}_R1.fastq.gz_fastqc.html\" '{output.html1}' && mv \"$TMP/{wildcards.sample}_R1.fastq.gz_fastqc.zip\" '{output.zip1}' && "
            "mv \"$TMP/{wildcards.sample}_R2.fastq.gz_fastqc.html\" '{output.html2}' && mv \"$TMP/{wildcards.sample}_R2.fastq.gz_fastqc.zip\" '{output.zip2}' && rm -rf \"$TMP\""
        )
        shell(
            f"python '{script}' --sample {wildcards.sample} --input-r1 '{_os.path.join(workdir, input.r1_input)}' --input-r2 '{_os.path.join(workdir, input.r2_input)}' "
            f"--output-r1 '{_os.path.join(workdir, input.r1)}' --output-r2 '{_os.path.join(workdir, input.r2)}' --repair-log '{_os.path.join(workdir, input.repair_log)}' "
            f"--fastqc-html '{_os.path.join(workdir, output.html1)}' --fastqc-html-r2 '{_os.path.join(workdir, output.html2)}'"
        )

rule trim_galore:
    input: r1 = RESDIR + "/fastq_ready/{sample}/{sample}_R1.fastq.gz", r2 = RESDIR + "/fastq_ready/{sample}/{sample}_R2.fastq.gz",
    output: r1 = RESDIR + "/trimmed/{sample}_R1_val_1.fq.gz", r2 = RESDIR + "/trimmed/{sample}_R2_val_2.fq.gz", report1 = RESDIR + "/trimmed/{sample}_R1.fastq.gz_trimming_report.txt", report2 = RESDIR + "/trimmed/{sample}_R2.fastq.gz_trimming_report.txt",
    params: resdir = RESDIR, tmpdir = RESDIR + "/.tmp_trim_{sample}",
    threads: config["threads"]
    log: "logs/trim_galore/{sample}.log"
    shell:
        "TMP='{params.tmpdir}'; rm -rf \"$TMP\"; trap 'rm -rf \"$TMP\"' EXIT; "
        "mkdir -p \"$TMP\" '{params.resdir}/trimmed'; cp -L '{input.r1}' '{input.r2}' \"$TMP\"/; "
        "trim_galore --paired --output_dir \"$TMP\" --cores {threads} \"$TMP/{wildcards.sample}_R1.fastq.gz\" \"$TMP/{wildcards.sample}_R2.fastq.gz\" >> {log} 2>&1; "
        "mv \"$TMP/{wildcards.sample}_R1_val_1.fq.gz\" \"$TMP/{wildcards.sample}_R2_val_2.fq.gz\" \"$TMP/{wildcards.sample}_R1.fastq.gz_trimming_report.txt\" \"$TMP/{wildcards.sample}_R2.fastq.gz_trimming_report.txt\" '{params.resdir}/trimmed/' && rm -rf \"$TMP\""

def _aggregate_qc_input(wildcards):
    samples_file = checkpoints.validate.get(samples=RESDIR + "/validated/samples.txt").output[0]
    path = _os.path.join(PIPELINE_DIR, samples_file)
    try:
        with open(path) as f:
            slist = [s.strip() for s in f if s.strip()]
    except FileNotFoundError:
        return [samples_file]
    return [samples_file] + expand(RESDIR + "/trimmed/{s}_R1_val_1.fq.gz", s=slist)

rule all_qc:
    input: _aggregate_qc_input,

# ===========================================================================
# MODULE 2: Alignment
# ===========================================================================

def _aggregate_align_input(wildcards):
    samples = get_valid_samples()
    if not samples:
        return [RESDIR + "/validated/samples.txt"]
    return [RESDIR + "/validated/samples.txt"] + expand(RESDIR + "/trimmed/{s}_R1_val_1.fq.gz", s=samples)

rule check_align_inputs:
    input: _aggregate_align_input,
    output: touch(RESDIR + "/.check_align_ok"),
    run:
        script = _os.path.join(PIPELINE_DIR, "scripts", "check_inputs.py")
        samples_path = _os.path.join(PIPELINE_DIR, RESDIR, "validated", "samples.txt")
        trim_path = _os.path.join(PIPELINE_DIR, RESDIR, "trimmed")
        subprocess.run([sys.executable, script, "align", config["genome_bwa"], samples_path, trim_path], check=True, cwd=PIPELINE_DIR)

rule bwa_align:
    input:
        _check = RESDIR + "/.check_align_ok",
        r1 = RESDIR + "/trimmed/{sample}_R1_val_1.fq.gz",
        r2 = RESDIR + "/trimmed/{sample}_R2_val_2.fq.gz",
    output: bam = RESDIR + "/bam/{sample}.bam",
    params: genome = config["genome_bwa"],
    threads: config["threads"]
    log: "logs/bwa/{sample}.log"
    shell:
        "bwa mem -t {threads} '{params.genome}' '{input.r1}' '{input.r2}' | samtools sort -@ {threads} -o '{output.bam}' - >> {log} 2>&1"

def _aggregate_align_target(wildcards):
    return expand(RESDIR + "/bam/{s}.bam", s=get_valid_samples()) if get_valid_samples() else [RESDIR + "/validated/samples.txt"]

rule all_align:
    input: _aggregate_align_target,

# ===========================================================================
# MODULE 3: Deduplication
# ===========================================================================

def _aggregate_dedup_input(wildcards):
    samples = get_valid_samples()
    if not samples:
        return [RESDIR + "/validated/samples.txt"]
    return [RESDIR + "/validated/samples.txt"] + expand(RESDIR + "/bam/{s}.bam", s=samples)

rule check_dedup_inputs:
    input: _aggregate_dedup_input,
    output: touch(RESDIR + "/.check_dedup_ok"),
    run:
        subprocess.run([sys.executable, "scripts/check_inputs.py", "dedup", _os.path.join(PIPELINE_DIR, RESDIR, "validated", "samples.txt"), _os.path.join(PIPELINE_DIR, RESDIR, "bam")], check=True, cwd=PIPELINE_DIR)

rule mark_duplicates:
    input:
        _check = RESDIR + "/.check_dedup_ok",
        bam = RESDIR + "/bam/{sample}.bam",
    output: bam = RESDIR + "/bam/{sample}.dedup.bam",
    params: tmpdir = RESDIR + "/.tmp_dedup_{sample}",
    threads: config["threads"]
    log: "logs/dedup/{sample}.log"
    shell:
        "TMP='{params.tmpdir}'; rm -rf \"$TMP\"; trap 'rm -rf \"$TMP\"' EXIT; mkdir -p \"$TMP\"; "
        "samtools sort -n -@ {threads} -o \"$TMP/n.bam\" '{input.bam}'; "
        "samtools fixmate -m \"$TMP/n.bam\" \"$TMP/f.bam\"; samtools sort -@ {threads} -o \"$TMP/s.bam\" \"$TMP/f.bam\"; "
        "samtools markdup -r -s \"$TMP/s.bam\" '{output.bam}' >> {log} 2>&1 && samtools index '{output.bam}' && rm -rf \"$TMP\""

def _aggregate_dedup_target(wildcards):
    return expand(RESDIR + "/bam/{s}.dedup.bam", s=get_valid_samples()) if get_valid_samples() else [RESDIR + "/validated/samples.txt"]

rule all_dedup:
    input: _aggregate_dedup_target,

# ===========================================================================
# MODULE 4: Peak Calling
# ===========================================================================

def _aggregate_peaks_input(wildcards):
    samples = get_valid_samples()
    tf = [s for s in samples if s != CONTROL]
    if not tf:
        return [RESDIR + "/validated/samples.txt"]
    return [RESDIR + "/validated/samples.txt"] + expand(RESDIR + "/bam/{s}.dedup.bam", s=tf) + [RESDIR + f"/bam/{CONTROL}.dedup.bam"]

rule check_peaks_inputs:
    input: _aggregate_peaks_input,
    output: touch(RESDIR + "/.check_peaks_ok"),
    run:
        subprocess.run([sys.executable, "scripts/check_inputs.py", "peaks", _os.path.join(PIPELINE_DIR, RESDIR, "validated", "samples.txt"), _os.path.join(PIPELINE_DIR, RESDIR, "bam"), CONTROL], check=True, cwd=PIPELINE_DIR)

rule macs3_callpeak:
    input:
        _check = RESDIR + "/.check_peaks_ok",
        bam = RESDIR + "/bam/{sample}.dedup.bam",
    output:
        peaks = RESDIR + "/peaks/{sample}_peaks.narrowPeak",
        summits = RESDIR + "/peaks/{sample}_summits.bed",
    params:
        control = RESDIR + f"/bam/{CONTROL}.dedup.bam",
        outdir = RESDIR + "/peaks",
        pval = config.get("macs3_pvalue", 0.05),
        extsize = config.get("macs3_extsize", 200),
        nolambda = "--nolambda" if config.get("macs3_nolambda", True) else "",
    log: "logs/macs3/{sample}.log"
    shell:
        "if [ \"{wildcards.sample}\" = \"{CONTROL}\" ]; then touch '{output.peaks}' '{output.summits}'; "
        "else macs3 callpeak -t '{input.bam}' -c '{params.control}' -n '{wildcards.sample}' --outdir '{params.outdir}' "
        "-p {params.pval} --extsize {params.extsize} {params.nolambda} --format BAM >> {log} 2>&1 && "
        "rm -f '{params.outdir}/{wildcards.sample}_peaks.xls' '{params.outdir}/{wildcards.sample}_model.r'; fi"

def _aggregate_peaks_target(wildcards):
    samples = get_valid_samples()
    tf = [s for s in samples if s != CONTROL]
    return expand(RESDIR + "/peaks/{s}_peaks.narrowPeak", s=tf) if tf else [RESDIR + "/validated/samples.txt"]

rule all_peaks:
    input: _aggregate_peaks_target,

# ===========================================================================
# MODULE 5: BigWig
# ===========================================================================

def _aggregate_bigwig_input(wildcards):
    samples = get_valid_samples()
    if not samples:
        return [RESDIR + "/validated/samples.txt"]
    return [RESDIR + "/validated/samples.txt"] + expand(RESDIR + "/bam/{s}.dedup.bam", s=samples)

rule check_bigwig_inputs:
    input: _aggregate_bigwig_input,
    output: touch(RESDIR + "/.check_bigwig_ok"),
    run:
        subprocess.run([sys.executable, "scripts/check_inputs.py", "bigwig", _os.path.join(PIPELINE_DIR, RESDIR, "validated", "samples.txt"), _os.path.join(PIPELINE_DIR, RESDIR, "bam")], check=True, cwd=PIPELINE_DIR)

rule bam_coverage:
    input:
        _check = RESDIR + "/.check_bigwig_ok",
        bam = RESDIR + "/bam/{sample}.dedup.bam",
    output: bw = RESDIR + "/bigwig/{sample}.rpkm.bw",
    params:
        blacklist = f"--blackListFileName '{BLACKLIST}'" if BLACKLIST else "",
    threads: config["threads"]
    log: "logs/deeptools/{sample}.log"
    shell:
        "bamCoverage -b '{input.bam}' -o '{output.bw}' --normalizeUsing RPKM --binSize 10 {params.blacklist} -p {threads} >> {log} 2>&1"

def _aggregate_bigwig_target(wildcards):
    return expand(RESDIR + "/bigwig/{s}.rpkm.bw", s=get_valid_samples()) if get_valid_samples() else [RESDIR + "/validated/samples.txt"]

rule all_bigwig:
    input: _aggregate_bigwig_target,
