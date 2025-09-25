#Code written in conjunction with openAI ChatGPT. Uses Python and R code.
#I reccomend running this code in a conda environment -- Code listed below to run in terminal.
#Also, make sure to install the correct packages in your R environment -- Code listed below

#####-- _Conda Code_
conda create -n chipseq_env -c bioconda -c conda-forge \
  bowtie2 samtools macs3 fastqc multiqc r-base r-essentials \
  r-biocmanager r-stringr r-readr r-ggplot2
conda activate chipseq_env
#####--

#####-- _R CODE_
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("ChIPseeker","GenomicFeatures","GenomeInfoDb","rtracklayer","org.Sc.sgd.db","txdbmaker"))
install.packages(c("readr","ggplot2"))
#####--
