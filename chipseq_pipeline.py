#!/usr/bin/env python3
"""
ChIP-seq Peak Calling Pipeline (Ubuntu 24.04)
- QC: FastQC + MultiQC (optional; runs if found)
- Index: bowtie2-build (skips if index exists)
- Align: bowtie2 (paired-end), then samtools sort & index
- Peaks: MACS3 (tries model; falls back to --nomodel --extsize 147 if needed)
- Annotate: R/ChIPseeker using your EXACT GFF/GTF (TxDb from GFF)

Usage (BY4742 example):
  python chipseq_pipeline.py \
    --fasta ref/BY4742.fa \
    --gff   ref/BY4742.gff3 \
    --r1 fastq/TF_R1.fastq.gz --r2 fastq/TF_R2.fastq.gz \
    --input-r1 fastq/Input_R1.fastq.gz --input-r2 fastq/Input_R2.fastq.gz \
    --sample TF_BY4742 --genome-size 1.2e7 --threads 8 --outdir run_BY4742

Minimal (no control):
  python chipseq_pipeline.py \
    --fasta ref/YJM789.fa --gff ref/YJM789.gff3 \
    --r1 fastq/TF_R1.fastq.gz --r2 fastq/TF_R2.fastq.gz \
    --sample TF_YJM --outdir run_YJM

Notes:
- Set --genome-size for your organism (yeast ~1.2e7).
- FASTQ must be paired-end (R1/R2). Control (Input) is optional.
- Requires: bowtie2, samtools, macs3; optionally fastqc, multiqc; R + packages.
"""

import argparse, os, sys, shutil, subprocess, textwrap, time, pathlib, re

def sh(cmd, log_path, cwd=None, env=None, check=True):
    """Run a shell command, tee stdout/stderr into a log file, and also echo to console."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "a") as lf:
        lf.write(f"\n$ {cmd}\n")
        lf.flush()
        proc = subprocess.Popen(cmd, shell=True, cwd=cwd, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, bufsize=1)
        for line in proc.stdout:
            sys.stdout.write(line)
            lf.write(line)
        proc.wait()
        lf.write(f"[exit_code={proc.returncode}]\n")
        lf.flush()
        if check and proc.returncode != 0:
            raise RuntimeError(f"Command failed (exit {proc.returncode}). See log: {log_path}")
        return proc.returncode

def which(tool):
    return shutil.which(tool) is not None

def require_tools(tools):
    missing = [t for t in tools if not which(t)]
    if missing:
        raise SystemExit(f"ERROR: Missing required tools: {', '.join(missing)}\n"
                         f"Install via conda:\n  conda install -c bioconda -c conda-forge {' '.join(missing)}")

def derive_index_basename(fasta_path, out_ref_dir):
    stem = pathlib.Path(fasta_path).stem
    return str(pathlib.Path(out_ref_dir) / stem)

def bowtie2_index_exists(index_base):
    # bowtie2-build creates files like <base>.1.bt2 etc.
    expected = [f"{index_base}.{i}.bt2" for i in [1,2,3,4]] + [f"{index_base}.rev.1.bt2", f"{index_base}.rev.2.bt2"]
    return all(os.path.exists(p) for p in expected)

def make_dirs(base):
    subdirs = ["logs", "qc", "ref", "bam", "peaks", "annotation", "cmds"]
    for sd in subdirs:
        os.makedirs(os.path.join(base, sd), exist_ok=True)
    return {sd: os.path.join(base, sd) for sd in subdirs}

def write_cmd(snapshot_path, cmd):
    with open(snapshot_path, "a") as f:
        f.write(cmd.strip() + "\n")

def write_r_script(r_path):
    r_code = r"""#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(rtracklayer)
  library(GenomeInfoDb)
  library(AnnotationDbi)
  library(org.Sc.sgd.db)  # OK if not perfect; we still get structural annotation
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: annotate_chipseq.R <narrowPeak> <gff_path> <outdir> <tss_up> <tss_down>")
}
peak_file <- args[1]
gff_path  <- args[2]
outdir    <- args[3]
tss_up    <- as.integer(args[4])
tss_down  <- as.integer(args[5])

message("Reading peaks: ", peak_file)
peaks <- readPeakFile(peak_file)

message("Building TxDb from: ", gff_path)
txdb <- txdbmaker::makeTxDbFromGFF(gff_path)  # robust for fresh GenomicFeatures versions

message("Annotating peaks...")
peakAnno <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(tss_up, tss_down), annoDb = "org.Sc.sgd.db")

# Save table
anno_df <- as.data.frame(peakAnno)
out_csv <- file.path(outdir, "chipseq_annotation.csv")
readr::write_csv(anno_df, out_csv)
message("Wrote: ", out_csv)

# Plots
pdf(file.path(outdir, "annotation_bar.pdf"), width=7, height=5)
print(plotAnnoBar(peakAnno))
dev.off()

pdf(file.path(outdir, "dist_to_TSS.pdf"), width=7, height=5)
print(plotDistToTSS(peakAnno))
dev.off()

pdf(file.path(outdir, "peak_heatmap.pdf"), width=7, height=5)
print(plotAnnoPie(peakAnno))
dev.off()

message("Annotation done.")
"""
    with open(r_path, "w") as f:
        f.write(r_code)

def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
        ChIP-seq Peak Calling Pipeline (paired-end)

        REQUIRED:
          --fasta      Path to reference FASTA (matching your strain, e.g., BY4742.fa)
          --gff        Path to matching GFF3/GTF annotation (same assembly as FASTA)
          --r1 --r2    Treatment FASTQ R1/R2 (ChIP)

        OPTIONAL:
          --input-r1/--input-r2  Control (Input DNA) FASTQs
          --genome-size          Effective genome size (default: 1.2e7 for yeast)
          --sample               Sample name (used in output file names)
          --threads              Bowtie2/samtools threads
          --outdir               Output directory (will be created)
          --nomodel              Force MACS3 --nomodel --extsize 147
          --skip-qc              Skip FastQC/MultiQC
          --tss-up/down          TSS window for annotation (default: -1000 / 1000)
        """)
    )
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--r1", required=True, help="ChIP R1 fastq.gz")
    ap.add_argument("--r2", required=True, help="ChIP R2 fastq.gz")
    ap.add_argument("--input-r1", default=None, help="Input DNA R1 fastq.gz (optional)")
    ap.add_argument("--input-r2", default=None, help="Input DNA R2 fastq.gz (optional)")
    ap.add_argument("--sample", default="chipseq_sample")
    ap.add_argument("--genome-size", default="1.2e7")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--outdir", default="chipseq_run")
    ap.add_argument("--nomodel", action="store_true", help="Use MACS3 --nomodel --extsize 147")
    ap.add_argument("--skip-qc", action="store_true")
    ap.add_argument("--tss-up", type=int, default=-1000)
    ap.add_argument("--tss-down", type=int, default=1000)
    args = ap.parse_args()

    # Check tools
    require_tools(["bowtie2", "samtools", "macs3"])
    # QC tools are optional
    qc_tools = [t for t in ["fastqc", "multiqc"] if which(t)]
    if not args.skip_qc and len(qc_tools) == 0:
        print("[info] FastQC/MultiQC not found; QC will be skipped.")
        args.skip_qc = True

    # Layout
    out = make_dirs(args.outdir)
    logs = out["logs"]
    refd = out["ref"]
    bamd = out["bam"]
    peaksd = out["peaks"]
    annd = out["annotation"]
    cmds = out["cmds"]

    # Copy/reference files (optional): we just record absolute paths
    fasta = os.path.abspath(args.fasta)
    gff   = os.path.abspath(args.gff)
    r1    = os.path.abspath(args.r1)
    r2    = os.path.abspath(args.r2)
    input_r1 = os.path.abspath(args.input_r1) if args.input_r1 else None
    input_r2 = os.path.abspath(args.input_r2) if args.input_r2 else None

    # 0) QC
    if not args.skip_qc:
        qc_cmd = f"fastqc {r1} {r2} -o {out['qc']}"
        write_cmd(os.path.join(cmds, "00_fastqc.cmd"), qc_cmd)
        sh(qc_cmd, os.path.join(logs, "00_fastqc.log"))
        if input_r1 and input_r2:
            qc_cmd2 = f"fastqc {input_r1} {input_r2} -o {out['qc']}"
            write_cmd(os.path.join(cmds, "00_fastqc.cmd"), qc_cmd2)
            sh(qc_cmd2, os.path.join(logs, "00_fastqc_input.log"))
        if which("multiqc"):
            mqc_cmd = f"multiqc {out['qc']} -o {out['qc']}"
            write_cmd(os.path.join(cmds, "01_multiqc.cmd"), mqc_cmd)
            sh(mqc_cmd, os.path.join(logs, "01_multiqc.log"))

    # 1) bowtie2-build
    index_base = derive_index_basename(fasta, refd)
    if not bowtie2_index_exists(index_base):
        build_cmd = f"bowtie2-build {fasta} {index_base}"
        write_cmd(os.path.join(cmds, "02_bowtie2_build.cmd"), build_cmd)
        sh(build_cmd, os.path.join(logs, "02_bowtie2_build.log"))
    else:
        print("[info] bowtie2 index already exists; skipping build.")

    # 2) Alignment (treatment)
    sam_treat = os.path.join(bamd, f"{args.sample}.sam")
    bam_treat = os.path.join(bamd, f"{args.sample}_sorted.bam")
    align_cmd = (
        f"bowtie2 --very-sensitive -x {index_base} -1 {r1} -2 {r2} "
        f"-p {args.threads} -S {sam_treat}"
    )
    write_cmd(os.path.join(cmds, "03_align_treatment.cmd"), align_cmd)
    sh(align_cmd, os.path.join(logs, "03_align_treatment.log"))

    # Convert, sort, index (treatment)
    sort_cmd = f"samtools view -@ {args.threads} -bS {sam_treat} | samtools sort -@ {args.threads} -o {bam_treat}"
    write_cmd(os.path.join(cmds, "04_sort_treatment.cmd"), sort_cmd)
    sh(sort_cmd, os.path.join(logs, "04_sort_treatment.log"))
    os.remove(sam_treat) if os.path.exists(sam_treat) else None
    idx_cmd = f"samtools index {bam_treat}"
    write_cmd(os.path.join(cmds, "05_index_treatment.cmd"), idx_cmd)
    sh(idx_cmd, os.path.join(logs, "05_index_treatment.log"))

    # 3) Alignment (control) if provided
    bam_control = None
    if input_r1 and input_r2:
        sam_ctrl = os.path.join(bamd, f"{args.sample}_input.sam")
        bam_ctrl = os.path.join(bamd, f"{args.sample}_input_sorted.bam")
        ctrl_align_cmd = (
            f"bowtie2 --very-sensitive -x {index_base} -1 {input_r1} -2 {input_r2} "
            f"-p {args.threads} -S {sam_ctrl}"
        )
        write_cmd(os.path.join(cmds, "06_align_control.cmd"), ctrl_align_cmd)
        sh(ctrl_align_cmd, os.path.join(logs, "06_align_control.log"))

        ctrl_sort_cmd = f"samtools view -@ {args.threads} -bS {sam_ctrl} | samtools sort -@ {args.threads} -o {bam_ctrl}"
        write_cmd(os.path.join(cmds, "07_sort_control.cmd"), ctrl_sort_cmd)
        sh(ctrl_sort_cmd, os.path.join(logs, "07_sort_control.log"))
        os.remove(sam_ctrl) if os.path.exists(sam_ctrl) else None
        ctrl_idx_cmd = f"samtools index {bam_ctrl}"
        write_cmd(os.path.join(cmds, "08_index_control.cmd"), ctrl_idx_cmd)
        sh(ctrl_idx_cmd, os.path.join(logs, "08_index_control.log"))
        bam_control = bam_ctrl

    # 4) MACS3 peak calling
    name = args.sample
    peaks_prefix = os.path.join(peaksd, name)
    base_macs_cmd = f"macs3 callpeak -t {bam_treat} " + (f"-c {bam_control} " if bam_control else "") + \
                    f"-f BAM -g {args.genome_size} -n {name} --outdir {peaksd} -q 0.05"

    if args.nomodel:
        base_macs_cmd += " --nomodel --extsize 147"

    write_cmd(os.path.join(cmds, "09_macs3.cmd"), base_macs_cmd)
    macs_log = os.path.join(logs, "09_macs3.log")
    sh(base_macs_cmd, macs_log, check=False)  # don't hard-fail—check for model warning

    # Fallback: re-run with --nomodel if model-building failure found and user didn't force --nomodel
    fallback_needed = False
    if not args.nomodel:
        if os.path.exists(macs_log):
            with open(macs_log) as f:
                text = f.read()
            if re.search(r"MACS3 needs at least 100 paired peaks", text, re.IGNORECASE):
                fallback_needed = True

    if fallback_needed:
        fb_cmd = base_macs_cmd + " --nomodel --extsize 147"
        write_cmd(os.path.join(cmds, "10_macs3_fallback.cmd"), fb_cmd)
        sh(fb_cmd, os.path.join(logs, "10_macs3_fallback.log"))

    # 5) R annotation (ChIPseeker)
    # pick the latest narrowPeak
    narrowpeak = os.path.join(peaksd, f"{name}_peaks.narrowPeak")
    if not os.path.exists(narrowpeak):
        # sometimes MACS uses slightly different naming; look for *peaks.narrowPeak
        candidates = sorted([p for p in pathlib.Path(peaksd).glob(f"{name}*peaks.narrowPeak")], key=os.path.getmtime)
        if candidates:
            narrowpeak = str(candidates[-1])

    if not os.path.exists(narrowpeak):
        print(f"[warn] No narrowPeak file found for {name}. Skipping annotation.")
    else:
        r_script = os.path.join(args.outdir, "annotation", "annotate_chipseq.R")
        write_r_script(r_script)
        os.chmod(r_script, 0o755)
        r_cmd = f"Rscript {r_script} {narrowpeak} {gff} {annd} {args.tss_up} {args.tss_down}"
        write_cmd(os.path.join(cmds, "11_annotate.cmd"), r_cmd)
        # We don't enforce check=True for R so the rest of pipeline still leaves artifacts if R lacks pkgs.
        sh(r_cmd, os.path.join(logs, "11_annotate.log"), check=False)

    # 6) Write a small run summary
    with open(os.path.join(args.outdir, "RUN_SUMMARY.txt"), "w") as rs:
        rs.write(textwrap.dedent(f"""
        ChIP-seq pipeline run summary
        =============================
        Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}
        Sample:    {args.sample}

        Inputs:
          FASTA: {fasta}
          GFF:   {gff}
          ChIP:  {r1}
                 {r2}
          Input: {input_r1 or '(none)'}
                 {input_r2 or '(none)'}

        Threads: {args.threads}
        Genome size: {args.genome_size}
        MACS3 nomodel: {'YES' if args.nomodel else 'NO (fallback on need)'}

        Outputs:
          BAM/TBI:   {bamd}
          Peaks:     {peaksd}
          Logs:      {logs}
          Annotation:{annd}

        Commands recorded under: {cmds}
        """).strip()+"\n")

    print("\n✅ Done. Check logs/ for detailed outputs and annotation/ for CSV + plots.")

if __name__ == "__main__":
    main()

