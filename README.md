# ChIP-seq Peak Calling Pipeline

This repository provides a reproducible pipeline to process **paired-end ChIP-seq data** for *Saccharomyces cerevisiae* (BY4742, YJM789, RM11, etc.) from raw FASTQ files to **annotated peaks**.

The pipeline performs:

1. **QC** with FastQC and MultiQC
2. **Alignment** with Bowtie2 (`--very-sensitive`)
3. **Sorting + Indexing** with Samtools
4. **Peak Calling** with MACS3 (paired-end mode, `-f BAMPE`)
5. **Annotation** with R/Bioconductor (`ChIPseeker`, `TxDbMaker`, `org.Sc.sgd.db`)

All steps log their outputs to `logs/`, and results are saved under `runs/<strain>/`.

---

## File Organization

```text
chipseq_project/
├── data/                  # Raw FASTQ files
│   └── STRAIN_NAME/
│       ├── TF_R1.fastq.gz
│       ├── TF_R2.fastq.gz
│       ├── CTRL_R1.fastq.gz
│       └── CTRL_R2.fastq.gz
├── refs/                  # Reference files
│   └── STRAIN_NAME/
│       ├── genome.fasta   # reference genome (FASTA)
│       └── genome.gff3    # matching annotation (GFF3/GTF)
├── scripts/               # Pipeline + R annotator
│   ├── chipseq_pipeline.py
│   └── annotate_peaks.R
├── runs/                  # Pipeline outputs (bam, peaks, etc.)
├── logs/                  # Per-step log files
└── README.md
```

Replace `STRAIN_NAME` with e.g. **BY4742**, **YJM789**, or **RM11**.

---

## Installation

### Conda/Mamba (recommended)

```bash
mamba create -n chipseq_env -y -c conda-forge -c bioconda \
  python=3.11 macs3 bowtie2 samtools bedtools fastqc multiqc \
  r-base=4.4 r-tidyverse r-data.table r-stringr \
  bioconductor-genomicranges bioconductor-genomicfeatures bioconductor-rtracklayer \
  bioconductor-annotationdbi bioconductor-genomeinfodb \
  bioconductor-org.sc.sgd.db bioconductor-chipseeker bioconductor-txdbmaker
```

---

## Usage

```bash
python scripts/chipseq_pipeline.py \
  --strain STRAIN_NAME \
  --ref-fasta refs/STRAIN_NAME/genome.fasta \
  --gff refs/STRAIN_NAME/genome.gff3 \
  --r1 data/STRAIN_NAME/TF_R1.fastq.gz \
  --r2 data/STRAIN_NAME/TF_R2.fastq.gz \
  --control-r1 data/STRAIN_NAME/CTRL_R1.fastq.gz \
  --control-r2 data/STRAIN_NAME/CTRL_R2.fastq.gz \
  --name TF_STRAIN_NAME \
  --threads 8 \
  --macs-genome 1.2e7 \
  --bowtie2-very-sensitive \
  --paired-end
```

Outputs:

* QC → `runs/run<STRAIN_NAME>/qc/`
* BAM → `runs/run<STRAIN_NAME>/bam/`
* Peaks → `runs/run<STRAIN_NAME>/peaks/`
* Annotation → `annotation/TF_STRAIN_NAME_annot.csv`
* Logs → `logs/`

---

## Common Errors & Fixes

### 1. Conda solver conflict (`libmamba Could not solve for environment specs`)

* Cause: mismatch between R version and Bioconductor packages.
* Fix: pin `r-base=4.4` and include `-c conda-forge -c bioconda`.

### 2. MACS3 model building error

```
WARNING: MACS3 needs at least 100 paired peaks...
```

* Fix: use fallback mode:

  ```bash
  macs3 callpeak -t TF.bam -c CTRL.bam -f BAMPE -g 1.2e7 -n TF --nomodel --extsize 147
  ```

### 3. Annotation error: `Unknown ID type`

* Cause: strain-specific GFF (e.g., BY4742 contigs).
* Fix: annotation will still assign structural features (promoter, exon, intron). Symbol mapping is only guaranteed for S288C (org.Sc.sgd.db).

### 4. Empty BAM or 0% alignment

* Cause: wrong reference genome (FASTA) or corrupted FASTQ.
* Fix: recheck that `bowtie2-build` used the same FASTA as alignment.

---

## Attribution & Citations

This pipeline was co-written with assistance from **ChatGPT (OpenAI, 2025)**.
If you use this code, please include proper software citations:

* **FastQC**: Andrews S. FastQC: A Quality Control Tool for High Throughput Sequence Data. (2010).
* **MultiQC**: Ewels P, et al. *MultiQC: summarize analysis results for multiple tools and samples in a single report.* Bioinformatics. 2016.
* **Bowtie2**: Langmead B, Salzberg SL. *Fast gapped-read alignment with Bowtie 2.* Nat Methods. 2012.
* **Samtools**: Li H, et al. *The Sequence Alignment/Map format and SAMtools.* Bioinformatics. 2009.
* **MACS3**: Zhang Y, et al. *Model-based Analysis of ChIP-Seq (MACS).* Genome Biol. 2008. (MACS2/3 are maintained versions).
* **ChIPseeker**: Yu G, et al. *ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization.* Bioinformatics. 2015.
* **TxDbMaker**: Bioconductor project (GenomicFeatures / txdbmaker).
* **org.Sc.sgd.db**: Bioconductor annotation package (S. cerevisiae S288C genome).

---

## Reproducibility Checklist

Before running the pipeline, ensure:

1. **Environment**

   * Export Conda env:

     ```bash
     conda env export --name chipseq_env > environment.yml
     ```
   * Record software versions:

     ```bash
     fastqc --version
     bowtie2 --version
     samtools --version
     macs3 --version
     Rscript -e "sessionInfo()"
     ```

2. **References**

   * Record source and version of reference FASTA and GFF3 (e.g., *SGD release 2022*).
   * Ensure FASTA and GFF3 come from the same release.

3. **Data**

   * Store raw FASTQ files in `data/STRAIN_NAME/`.
   * Do not rename reads beyond `R1` / `R2` convention.

4. **Pipeline Outputs**

   * Keep all logs in `logs/`.
   * Keep QC reports (`runs/run<STRAIN_NAME>/qc/`) for auditing.

5. **Publication**

   * Cite all tools listed above.
   * State strain background (BY4742, YJM789, RM11, etc.) and reference genome used.

