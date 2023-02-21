
###Basic commands to run cellranger-arc
##applies to versions 1.0.1 and 2.0.0
#References: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc
#Reference Genome: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

datadirectory="/.../Cellranger-arc/preprocessing/"
reference="refdata-cellranger-arc-GRCh38-2020-A"

cellranger-arc count --id ${SAMPLE} --reference ${reference} --libraries ${datadirectory}/${SAMPLE}/library.csv --disable-ui --localcores 16 --localmem 64