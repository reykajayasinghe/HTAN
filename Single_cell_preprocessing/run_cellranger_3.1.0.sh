
###Basic commands to run cellranger for single-cell or single-nuclei RNA-sequencing Data
#References: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
#Reference Genome:https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest? 

datadirectory="/.../fastq/directory/"
reference="refdata-cellranger-GRCh38-3.0.0"
cellranger="/../Software/cellranger-3.1.0/cellranger-cs/3.1.0/bin/count"

${cellranger} --id ${sample} -lib1 --fastqs ${datadirectory} --transcriptome ${reference} --localcores=50 --localmem=400