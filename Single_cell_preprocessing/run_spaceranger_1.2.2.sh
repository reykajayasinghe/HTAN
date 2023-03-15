
###Basic commands to run spaceranger for spatial transcriptomics Data
#References: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger

sample="SAMPLEID"
sample2="SAMPLEID_2"
TIF_image="TIF_location"
SLIDE_SERIAL_ID="serial_ID"
AREA="slide_area" #https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/slide-info
datadirectory="/.../fastq/directory/"
reference="refdata-gex-GRCh38-2020-A"

spaceranger count --id=${sample} --transcriptome=${reference} --fastqs=${datadirectory} --sample=${sample2} --image=${TIF_IMAGE} --slide=${SLIDE_SERIAL_ID} --area=${AREA} --reorient-images --localcores=32 --localmem=150
