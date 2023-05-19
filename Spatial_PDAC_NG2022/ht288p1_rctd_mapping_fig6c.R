###MAP ANNOTATED SCRNA from HT288P1 to ST HT288P1 sample
#Reference Figure 6C https://www.nature.com/articles/s41588-022-01157-1

##Base code from Liang-Bo Wang
##Modified by Reyka Jayasinghe (reyka@wustl.edu)
conda activate spotlight
working directory: /diskmnt/Datasets/mmy_scratch/PDAC/ST_analysis

library(devtools)
devtools::install_github("https://github.com/MarcElosua/SPOTlight")
devtools::install_git("https://github.com/MarcElosua/SPOTlight")

library(here)
library(arrow)
library(Seurat)
library(SPOTlight)
library(scatterpie)
library(patchwork)
library(ggrastr)
library(tidyverse)

scrna=readRDS("HT288P1_v2.rds")
newmetadata<-read.table("ht288p1_metadata.txt", row.names=1,sep = '\t', header = TRUE)
scrna <- AddMetaData(
  object = scrna,
  metadata = newmetadata)    


##Note certain annotations are changed to PDAC in celltypespecific_v3 metadata
##EMT/Other etc.

DefaultAssay(scrna) = "SCT"

Seurat::Idents(scrna) = "celltypespecific_v3"

scrna_cluster_markers_selected = Seurat::FindAllMarkers(
    object = scrna, 
    assay = "SCT",
    slot = "data",
    verbose = TRUE, 
    only.pos = TRUE,  # Only positive markers
    logfc.threshold = 1,
    min.pct = 0.5  # Most cells should express the markers
)

set.seed(202011)

sample="HT288P1"
st_merged = readRDS(paste0("HT288P1_H4_U1_processed.rds"))

spotlight_ls = spotlight_deconvolution(
    se_sc = scrna,
    counts_spatial = st_merged@assays$Spatial@counts,
    clust_vr = "celltypespecific_v3",
    cluster_markers = scrna_cluster_markers_selected,
    cl_n = 100,
    hvg = 3000,
    ntop = NULL,
    transf = "uv",
    method = "nsNMF",
    min_cont = 0.09
)

decon_mat = spotlight_ls[[2]]
cell_types_all = colnames(decon_mat)[which(colnames(decon_mat) != "res_ss")]

st_merged@meta.data = cbind(st_merged@meta.data[, 1:8], decon_mat)

#Cell type proportion
cell_types_all %>%
    walk(function(cell_type) {
        SpatialFeaturePlot(
            st_merged,
            features = cell_type, 
            pt.size.factor = 2,
            alpha = c(0, 1)
        ) &
            # ggplot2::scale_fill_viridis_c(limits = c(0, 1) , option = "plasma")
            ggplot2::scale_fill_gradientn(
                colours = heat.colors(10, rev = TRUE),
                limits = c(0, 1)
            )
        
        ggsave(glue::glue('cell_type_proportion_on_slides.{cell_type}.pdf'), width = 8, height = 5)
    })

#Deconvolution Pie Chart
pdf('HT288P1_H4_U1.pdf', width = 8, height = 6, useDingbats = FALSE)
SPOTlight::spatial_scatterpie(
    se_obj = st_merged,
    cell_types_all = cell_types_all,
    slice = 'HT288P1_H4_U1',
    img_path = '/outs/spatial/tissue_lowres_image.png',
    pie_scale = 0.4
)
dev.off()
