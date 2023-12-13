
###Basic commands to run cellranger for single-cell or single-nuclei RNA-sequencing Data
#References: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
#Reference Genome:https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest? 

datadirectory="/.../fastq/directory/"
reference="/diskmnt/Datasets/Reference/Cellranger-2020-A/refdata-gex-GRCh38-2020-A"
cellranger="/diskmnt/Software/cellranger-6.0.1/bin/cellranger"

cellranger count --id=${sample} --fastqs=${datadirectory} --sample=${sample} --include-introns --localmem=200 --localcores=30 --transcriptome=${reference}