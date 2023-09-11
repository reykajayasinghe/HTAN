##Reanalysis of CAF Subsets January-September 2023

#Load in relevant packages
library(Seurat)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)

working_dir="PDAC_dir/"

#Load in seurat object from Cui Zhou et al study 
all_ourstudy = readRDS(paste0(working_dir,"Everything_merge_v6.rds"))

#Add Metadata
celltypes<-read.table(file=paste0(working_dir,"Spatial_PDAC_Metadata_Celltypes.txt"), sep="\t", header = TRUE)

rownames(celltypes)<-celltypes$Barcodes
celltypes$Barcodes<-NULL

all_ourstudy <- AddMetaData(
  object = all_ourstudy,
  metadata = celltypes)


###Marker Genes from Publications of interest

#CAF Markers
Elyada_gen<-c("COL1A1","FAP","PDPN","DCN") #Mouse/Human
Wang_gen<-c("COL1A1","LUM")
Werba_gen<-c("DCN")

CAF_unique<-c(Elyada_gen,Wang_gen,Werba_gen)%>% unique()

#iCAF Markers
Elyada_iCAF<-c("PDGFRA","CXCL12","CFD")
Wang_iCAF<-c("APOD","PTGDS","C7","C3","MGP","FBLN1")
Galbo_iCAF<-c("CFD", "C3", "CXCL14", "CXCL12")
Li_iCAF<-c("GSN","APOD","CFD","PTGDS")
Werba_iCAF<-c("C3", "C7", "CFD", "PTGDS")
iCAF_unique<-c(Elyada_iCAF,Wang_iCAF,Werba_iCAF,Galbo_iCAF,Li_iCAF)%>% unique()

#myCAF Markers
Elyada_myCAF<-c("ACTA2","TAGLN","MMP11","MYL9","HOPX","POSTN","TPM1","TPM2")
Wang_myCAF<-c("CRIP1","MYH11","DSTN","TAGLN")
Galbo_myCAF<-c("ACTA2","MYH11", "MCAM", "TAGLN","MYLK")
Li_myCAF<-c("TAGLN","MYH11")
Werba_myCAF<-c("ACTA2", "MMP11")
myCAF_unique<-c(Elyada_myCAF,Wang_myCAF,Werba_myCAF,Galbo_myCAF,Li_myCAF)%>% unique()

Pericyte<-c( "RGS5", "CSPG4", "NOTCH3")

allmarkers<-c(CAF_unique,iCAF_unique,myCAF_unique)

Elyadal<-c(Elyada_gen,Elyada_iCAF,Elyada_myCAF) %>% unique()
Wangl<-c(Wang_gen,Wang_iCAF,Wang_myCAF)%>% unique()
Werbal<-c(Werba_gen,Werba_iCAF,Werba_myCAF)%>% unique()
Galbol<-c(Galbo_iCAF,Galbo_myCAF)%>% unique()
Lil<-c(Li_iCAF,Li_myCAF)%>% unique()

lt=list(Elyada=Elyadal,
        Wang=Wangl,
        Werba=Werbal,
        Galbo=Galbol,
        Li=Lil)

data=list_to_matrix(lt)

options(repr.plot.width=6, repr.plot.height=2)

dotplot=data %>% melt %>% mutate(`Var1` = factor(`Var1`, levels=allmarkers)) %>% 
    filter(!value == 0) %>% ggplot(aes(x=Var1, y=Var2, size = value)) +
    geom_point(alpha=0.8) + theme_minimal()+
    theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5))

print(dotplot)

pdf("CAF_Markers_dotplot.pdf", w=6, h=2)
print(dotplot)
dev.off()


##Recluster CAF object with pericytes

Idents(all_ourstudy)<-all_ourstudy@meta.data$cell_type_specific
all<-subset(all_ourstudy,idents=c("Pericyte","iCAF","CD133+_iCAF","myCAF","apCAF","CXCR4+_iCAF"))
Idents(all)<-all@meta.data$tissue
all_tumor_only<-subset(all,idents=c("Tumor"))
all_tumor_only <- SCTransform(all_tumor_only, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
all_tumor_only <- RunPCA(all_tumor_only, npcs = 30)
all_tumor_only <- RunUMAP(all_tumor_only, reduction = "pca", dims = 1:30)
all_tumor_only <- FindNeighbors(all_tumor_only, reduction = "pca", dims = 1:30,force.recalc=TRUE)
all_tumor_only <- FindClusters(all_tumor_only, resolution = 0.5)
all_tumor_only <- RunUMAP(all_tumor_only, dims = 1:30)

saveRDS(all_tumor_only,"CAFs_pericytes.rds")

#Reannotate CAF's by cluster 
obj<-readRDS(paste0(working_dir,"CAFs_pericytes.rds"))

#Load CAF object for tumor only object 
CAF=readRDS(paste0(working_dir,"CAF_v3_tumor_only.rds"))

Idents(CAF)<-CAF@meta.data$anno
CAF<-RenameIdents(CAF,
"apCAF_10" = "apCAF",
"CD133+_iCAF_16" = "CD133+_iCAF",
"CXCR4+_iCAF_9" = "CXCR4+_iCAF",
"iCAF_0" = "iCAF",
"iCAF_1" = "myCAF_1",
"iCAF_12"= "iCAF",
"iCAF_13"= "iCAF",
"iCAF_14"= "iCAF",
"iCAF_15"= "iCAF",
"iCAF_17"= "iCAF",
"iCAF_18"= "iCAF",
"iCAF_19"= "iCAF",
"iCAF_4"= "iCAF",
"iCAF_5"= "iCAF",
"iCAF_6"= "iCAF",
"iCAF_7"= "iCAF", 
"myCAF_2"= "myCAF_2", 
"myCAF_8"= "myCAF_2",
"myCAF_11"= "myCAF_split")
CAF@meta.data$anno2<-Idents(CAF)

#Transfer annotation to pericyte + CAF object
metadata<-CAF@meta.data %>% select(anno2)
obj_metadata<-obj@meta.data %>% select(cell_type_specific) %>% filter(cell_type_specific %in% c("Pericyte"))
colnames(obj_metadata)<-c("anno2")
allmeta<-rbind(obj_metadata,metadata)
obj <- AddMetaData(
  object = obj,
  metadata = allmeta)

#Figure related to CAF_pericyte object with specific cluster annotation added

my_cols<-c("CXCR4+_iCAF"="#448aff",
           "iCAF"="#2a9d8f",
           "CD133+_iCAF"="#355070",
           "myCAF_1"="#ff9800",
           "myCAF_2"="#f44336",
           "myCAF_split"="#ad1457",
           "Pericyte"="#fbc3bc",
           "apCAF"="#b18bda")

pdf("CAF_umap.pdf", w=6, h=6)
DimPlot(caf_pericyte,group.by=c("anno2"),cols=my_cols,raster=TRUE,pt.size=2)&coord_equal()
dev.off()


#Create figure for marker annotation and subset CAF annotation
allmarkers<-c(CAF_unique,iCAF_unique,myCAF_unique,Pericyte)

celltypeorder<-c('iCAF','CXCR4+_iCAF','CD133+_iCAF','myCAF','apCAF','Pericyte',
    'cDC1','Macrophage','Exhausted_T','Treg','PDAC','cDC2','Plasma','NK',
    'CD8_T','CD8_T_cytotoxic','Neutrophil','Endothelial','Tuft',
    'Proliferating_T','CD4_T','Monocyte','B','Erythrocyte','PanIN',
    'Islet','Duct_like_2','Mast','Duct_like_1','Acinar_REG+','Acinar','ADM')

all_ourstudy@meta.data= all_ourstudy@meta.data %>% mutate(
    `cell_type_specific` = factor(`cell_type_specific`, levels=celltypeorder))

b<-DotPlot(all_ourstudy,features=allmarkers,scale=TRUE,group.by="cell_type_specific",cols=c("azure","deeppink4"))&
  theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&
  ggtitle("CAF Markers")

annorder<-c('myCAF_1','myCAF_2','myCAF_split','iCAF','CXCR4+_iCAF','CD133+_iCAF','apCAF','Pericyte')
allmarkers<-c(iCAF_unique,myCAF_unique,Pericyte)
caf_pericyte@meta.data= caf_pericyte@meta.data %>% mutate(
    `anno2` = factor(`anno2`, levels=annorder))
a<-DotPlot(caf_pericyte,group.by=c("anno2"),features=allmarkers,scale=TRUE,cols=c("azure","deeppink4"))&
  theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5))&
  ggtitle("CAF Markers")

pdf("CAF_Markers.pdf", w=10, h=8)
print(b)
dev.off()

pdf("CAF_breakdown_Markers.pdf", w=9, h=3.5)
print(a)
dev.off()


####ST Related Analysis

#Remake spatial plots 
st=readRDS(paste0(working_dir,"/PDAC/HT259P1/H1/HT259P1_H1A2_U1_processed_label_transferred.rds"))

b<-SpatialFeaturePlot(st, features = "POSTN",image.alpha=0,slot="data", pt.size.factor = 1)+
scale_fill_gradientn(colors=c("white","red4"))

d<-SpatialFeaturePlot(st, features = "C7",image.alpha=0,slot="data", pt.size.factor = 1)+
scale_fill_gradientn(colors=c("white","red4"))

c<-SpatialFeaturePlot(st, features = "C7",alpha=0, image.alpha=1,slot="data", pt.size.factor = 1)+
scale_fill_gradientn(colors=c("white","red4"))

pdf(paste("CAF_ST.pdf", sep=""), width=7, height=5)
ggarrange(c,b,d, ncol = 3, nrow = 1)
dev.off()

