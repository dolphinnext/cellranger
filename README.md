Cell Ranger (v7.0.0) is a set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis. 


Steps:
  1. Cellranger count takes FASTQ files performs alignment, filtering, barcode counting, and UMI counting. It uses the Chromium cellular barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis.  cellranger count also processes Feature Barcoding data alongside Gene Expression reads.
  2. Cellranger aggr aggregates outputs from multiple runs of cellranger count, normalizing those runs to the same sequencing depth and then recomputing the feature-barcode matrices and analysis on the combined data. The aggr pipeline can be used to combine data from multiple samples into an experiment-wide feature-barcode matrix and analysis.

Please check following web site for detailed information:  https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count

Inputs: 
  - Reads
  - Mate (single or pair)

Citation:

If you use DolphinNext in your research, please cite: 
Yukselen, O., Turkyilmaz, O., Ozturk, A.R. et al. DolphinNext: a distributed data processing platform for high throughput genomics. BMC Genomics 21, 310 (2020). https://doi.org/10.1186/s12864-020-6714-x


Program Versions:
  - Cellranger-v7.0.0
  - fastqc=0.11.8
  - multiqc=1.7
  - Cellranger-atac-2.0.0.
  - bcl2fastq2=2.20.0.422
  - samtools=1.3
  - igvtools=2.5.3
  - bedtools=2.30.0

Containers:
  - Docker: dolphinnext/cellranger:6.0
  - Singularity: https://galaxyweb.umassmed.edu/pub/dolphinnext_singularity/UMMS-Biocore-singularity-cellranger-v6.simg


