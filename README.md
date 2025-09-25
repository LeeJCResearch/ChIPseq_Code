# ChIP-seq Peak Calling & Annotation Pipeline (Python + R)

This repository contains a single Python script, `chipseq_pipeline.py`, that performs:

* **QC**: FastQC → MultiQC (optional)
* **Indexing**: `bowtie2-build` for your **exact** strain’s genome
* **Alignment** (paired-end): Bowtie2 **`--very-sensitive`**
* **BAM processing**: `samtools view/sort/index`
* **Peak calling**: MACS3 (auto-fallback to `--nomodel --extsize 147` if needed)
* **Annotation**: R/ChIPseeker using a **TxDb built from your exact GFF**
  (strain-correct, not S288C-only)

It works well for *S. cerevisiae* strains like **BY4742**, **YJM789**, **RM11**, etc., as long as the **FASTA and GFF** match the assembly you aligned to.

---

## Quick Start

### 1) Clone & layout your files

Recommended repo layout (you can copy this as-is):

```
your-repo/
├── chipseq_pipeline.py
├── README.md
├── ref/                 # put reference files here
│   ├── BY4742.fa        # reference FASTA (strain-specific)
│   └── BY4742.gff3      # matching GFF3 or GTF annotation (same assembly!)
├── fastq/               # put FASTQs here (gzipped)
│   ├── TF_R1.fastq.gz   # treatment R1
│   ├── TF_R2.fastq.gz   # treatment R2
│   ├── Input_R1.fastq.gz  # (optional) input/control R1
│   └── Input_R2.fastq.gz  # (optional) input/control R2
└── runs/                # pipeline creates run folders here (optional)
```

**Naming tips (makes life easier):**

* FASTQs: `SAMPLE_R1.fastq.gz`, `SAMPLE_R2.fastq.gz` (paired-end).
* One **FASTA** and one **GFF3/GTF** per strain/assembly, with **matching** contig names.
* Keep extensions standard: `.fa`/`.fasta`, `.gff3`/`.gtf`, `.fastq.gz`.

> If your contigs look like `gi|...|gb|JRIR01000001.1|` in files, that’s fine—annotation handles it as long as **FASTA and GFF** come from the same source (same accession strings).

---

### 2) Create the software environment

**Option A: Conda (recommended)**

```bash
conda create -n chipseq_env -c bioconda -c conda-forge \
  bowtie2 samtools macs3 fastqc multiqc r-base r-essentials \
  r-biocmanager r-stringr r-readr r-ggplot2
conda activate chipseq_env

# R packages for annotation (run inside R)
R -q <<'RSCRIPT'
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("ChIPseeker","GenomicFeatures","GenomeInfoDb","rtracklayer","org.Sc.sgd.db","txdbmaker"))
install.packages(c("readr","ggplot2"))
RSCRIPT
```

**Option B: System packages**
You can install bowtie2, samtools, R, etc. via apt, but MACS3 and Bioconductor are much smoother in conda.

---

### 3) Run the pipeline

Basic example (**with Input/control**):

```bash
python chipseq_pipeline.py \
  --fasta ref/BY4742.fa \
  --gff   ref/BY4742.gff3 \
  --r1 fastq/TF_R1.fastq.gz --r2 fastq/TF_R2.fastq.gz \
  --input-r1 fastq/Input_R1.fastq.gz --input-r2 fastq/Input_R2.fastq.gz \
  --sample TF_BY4742 \
  --genome-size 1.2e7 \
  --threads 8 \
  --outdir runs/run_BY4742
```

Minimal example (**no Input/control**):

```bash
python chipseq_pipeline.py \
  --fasta ref/YJM789.fa \
  --gff   ref/YJM789.gff3 \
  --r1 fastq/TF_R1.fastq.gz --r2 fastq/TF_R2.fastq.gz \
  --sample TF_YJM789 \
  --genome-size 1.2e7 \
  --threads 8 \
  --outdir runs/run_YJM789
```

> **Bowtie2 sensitivity:** The script uses `bowtie2 --very-sensitive` by default.

---

## What the script does (stage by stage)

1. **QC (optional)**

   * Runs **FastQC** on all FASTQs (treatment + control if provided)
   * Aggregates with **MultiQC**
   * Outputs: `outdir/qc/*` and logs in `outdir/logs/`

2. **Indexing**

   * Builds a `bowtie2-build` index from your `--fasta` into `outdir/ref/`
   * Skips if index files already exist for that FASTA

3. **Alignment** (paired)

   * **Treatment**: `bowtie2 --very-sensitive -x index -1 R1 -2 R2`
   * Converts to BAM, sorts, and indexes with `samtools`
   * **Control** (optional): same as treatment
   * Outputs: sorted BAMs in `outdir/bam/`, logs in `outdir/logs/`

4. **Peak calling (MACS3)**

   * `macs3 callpeak -t TREAT.bam [-c CONTROL.bam] -f BAM -g <genome_size> -q 0.05`
   * If MACS3 can’t build a model (common with TFs/low peaks), it **auto-retries** with:

     ```
     --nomodel --extsize 147
     ```
   * Outputs: `outdir/peaks/<sample>_peaks.narrowPeak`, summits, model files

5. **Annotation (R/ChIPseeker)**

   * Builds a **TxDb directly from your `--gff`** (via `txdbmaker::makeTxDbFromGFF`)
   * Annotates peaks with **TSS window** (default: -1000..+1000)
   * Saves a CSV table and several PDFs in `outdir/annotation/`

6. **Reproducibility**

   * Every command is copied into `outdir/cmds/*.cmd`
   * Full console outputs are in `outdir/logs/*.log`
   * A text summary is saved at `outdir/RUN_SUMMARY.txt`

---

## Command-line options

```
--fasta <path>        Reference FASTA (must match your strain/assembly)
--gff <path>          Matching GFF3/GTF (same assembly & contig names as FASTA)
--r1/--r2 <paths>     Treatment paired FASTQs (gz ok)
--input-r1/--input-r2 Optional control (Input DNA) paired FASTQs
--sample <name>       Sample label used in output file names (default: chipseq_sample)
--genome-size <str>   MACS3 genome size (yeast ≈ 1.2e7). Accepts notation like 1.2e7.
--threads <int>       Threads for Bowtie2/Samtools (default: 8)
--outdir <path>       Output directory (created if missing)
--nomodel             Force MACS3 --nomodel --extsize 147 (skip model building)
--skip-qc             Don’t run FastQC/MultiQC
--tss-up <int>        ChIPseeker upstream TSS window (default: -1000)
--tss-down <int>      ChIPseeker downstream TSS window (default: 1000)
```

---

## Output structure

```
outdir/
├── annotation/
│   ├── annotate_chipseq.R          # generated helper script
│   ├── chipseq_annotation.csv      # annotated peaks
│   ├── annotation_bar.pdf
│   ├── dist_to_TSS.pdf
│   └── peak_heatmap.pdf
├── bam/
│   ├── <sample>_sorted.bam
│   └── <sample>_sorted.bam.bai
├── cmds/
│   └── *.cmd                       # exact commands for every step
├── logs/
│   └── *.log                       # full stdout/stderr of each command
├── peaks/
│   ├── <sample>_peaks.narrowPeak
│   ├── <sample>_summits.bed
│   └── <sample>_model.r            # if model built
├── qc/                              # FastQC & MultiQC outputs
└── ref/                             # bowtie2 index files for your FASTA
RUN_SUMMARY.txt
```

---

## Best practices & naming conventions

* **FASTA & GFF must match.** Use reference + annotation from the **same source** and version (same contig names/accessions).
  Example: `BY4742.fa` with `BY4742.gff3` from the same release.

* **FASTQ naming:** `SAMPLE_R1.fastq.gz`, `SAMPLE_R2.fastq.gz`. Avoid spaces; keep gzipped.

* **Sample names:** Pass a clear `--sample` (e.g., `TF_BY4742`). It becomes the prefix of peak files and BAMs.

* **Genome size:** Use `--genome-size 1.2e7` for budding yeast. (Adjust for other organisms.)

* **Alignment preset:** We use `bowtie2 --very-sensitive` (good default for TF ChIP-seq).

---

## Troubleshooting

* **MACS3 “needs at least 100 paired peaks” / model fail**

  * The pipeline auto-retries with `--nomodel --extsize 147`.
    You can force this from the start with `--nomodel`.

* **Zero peaks or poor alignment**

  * Confirm your FASTA and GFF are for the **same strain/assembly** you expect.
  * Check `logs/03_align_treatment.log` for alignment rate and warnings.
  * Inspect QC in `qc/` (if enabled).

* **Annotation says “Unknown ID type”**

  * That’s OK; ChIPseeker still adds **structural** annotations (promoter/exon/intron/TSS distances).
  * Symbol mapping to **S288C** genes may not apply for non-S288C strains; we prioritize strain-correct coordinates using your GFF.

* **R package errors**

  * Re-open the conda env and reinstall the Bioconductor packages shown above.
  * Make sure `txdbmaker` is installed (the script calls `txdbmaker::makeTxDbFromGFF`).

* **Permissions on cluster/shared systems**

  * Write outputs into a user-writable location (e.g., `runs/`).
  * If `fastqc` / `multiqc` are missing, use `--skip-qc` or install them.

---

## Reproducibility & provenance

* Every command is saved under `outdir/cmds/*.cmd`.
* Every console output is captured under `outdir/logs/*.log`.
* A human-readable `RUN_SUMMARY.txt` records inputs, tool settings, and outputs.

---

## Examples

**BY4742 with Input DNA:**

```bash
python chipseq_pipeline.py \
  --fasta ref/BY4742.fa \
  --gff   ref/BY4742.gff3 \
  --r1 fastq/TF_R1.fastq.gz --r2 fastq/TF_R2.fastq.gz \
  --input-r1 fastq/Input_R1.fastq.gz --input-r2 fastq/Input_R2.fastq.gz \
  --sample TF_BY4742 \
  --genome-size 1.2e7 \
  --threads 8 \
  --outdir runs/run_BY4742
```

**YJM789 without Input DNA:**

```bash
python chipseq_pipeline.py \
  --fasta ref/YJM789.fa \
  --gff   ref/YJM789.gff3 \
  --r1 fastq/TF_R1.fastq.gz --r2 fastq/TF_R2.fastq.gz \
  --sample TF_YJM789 \
  --genome-size 1.2e7 \
  --threads 8 \
  --outdir runs/run_YJM789
```

---

## Notes & extensions

* For multiple samples/strains, run the script once per sample with different `--sample`, `--fasta`, `--gff`, and FASTQs.
* If you want replicates + IDR or differential binding (e.g., DiffBind), open an issue—we can extend this to a Snakemake workflow.

---

## Citation

If you use this pipeline, please cite the underlying tools:

* **Bowtie2**: Langmead & Salzberg (2012)
* **SAMtools**: Li et al. (2009)
* **MACS3**: Zhang et al. (2008) / updated MACS3
* **ChIPseeker**: Yu et al. (2015)
* **Bioconductor/TxDb/GenomeInfoDb/rtracklayer**: respective packages

---

### Contact

Questions or bugs? Please open an issue with:

* Your command line
* The `RUN_SUMMARY.txt`
* Relevant `logs/*.log`

Happy peak-calling! 🧬✨
