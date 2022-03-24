# RNA_sequencing_code
## Input

Reads files: fastq

Genome: fasta

Annotation: gtf

## Pipelines

### 0 QC fastqc

### 1 STAR index & mapping

### 2 featurecounts quantification

### 3 DEG_pre-process

### 4 DEG_DESeq2

### 5 Visualization (heatmap, PCA, correlation)

### Useful manual links



(RNA sequencing: the teenage years)[https://www.nature.com/articles/s41576-019-0150-2]

STAR : https://github.com/alexdobin/STAR

Manul of STAR: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

Ensembl reference genome download:  http://www.ensembl.org/Gallus_gallus/Info/Index

Featurecounts manual: https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf


DESeq2 tutorial: https://lashlock.github.io/compbio/R_presentation.html
