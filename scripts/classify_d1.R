#!/usr/bin/env Rscript
#.libPaths("/nfs/users/nfs_f/fr7/anaconda2/envs/r_env/lib/R/library")
m <- modules::use("/Users/fr7/git_repos/single_cell/scripts/SC.R")
library(argparse)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(gridExtra)
library(stats)
library(randomForest)
library(loomR)

cols <- c(
  "UNDIFFERENTIATED.1" = "springgreen4",
  "UNDIFFERENTIATED.2" = "steelblue" ,  
  "UNDIFFERENTIATED.3" =  "goldenrod",
  "ENTEROCYTE.1" = "orange",
  "ENTEROCYTE.2" = "darksalmon",
  "ENTEROCYTE.3" = "darkolivegreen4",
  "ENTEROCYTE.4" = "magenta1",
  "ENTEROCYTE.ISG15" = "red",
  "ENTEROCYTE.ACE2" = "darkorchid4",
  "GOBLET" = "cornflowerblue",
  "TUFT" = "darkmagenta",
  "EE" = "chocolate"
)

dir <- '/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/day1_classify'
control <- readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/seurat_object.rds')
samples_directory<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/QC/combined'
meta_data<-read.table('/Users/fr7/git_repos/single_cell/metadata/samples.txt',header=T)
condition<-'infected.24'

#find the IDs of the samples that we want
samples<-meta_data[which(meta_data$condition == condition & meta_data$experiment == '4'),'sample_id']

#load seurat objects
objects<-list()
for (sample in samples){
  object<-readRDS(file.path(samples_directory,paste0(sample,".rds")))
  objects[[sample]]<-object
}

#merge all replicates into one object
day.1 <- merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects), project = "day1" )

#normalize, scale, run PCA
day.1<- m$normalize_and_run_pca(day.1,40,dir)
day.1<-m$do_clustering(day.1,35,resolution=0.5,min.dist=0.05)

#plots
markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Isg15","Ace2","Ptprc")
p1<-m$QC_violins(day.1, dir)
p2<-m$QC_scatters(day.1,dir)
p3<-m$plot_UMAPS(day.1,dir)
p4<-m$marker_violins(day.1,dir,markers_to_plot)
p5<-m$markers_feature_plot(day.1,dir,markers_to_plot)

#markers
DefaultAssay(day.1)<-"RNA"
markers<-m$find_all_markers(day.1,dir)
marker_dots<-m$marker_dotplot(day.1,markers,dir)

#transfer labels
DefaultAssay(control)<-"integrated"
control$celltype<-control@active.ident
anchors <- FindTransferAnchors(reference =  control, query = day.1, dims = 1:30 )
predictions <- TransferData(anchorset = anchors, refdata = control$celltype, dims = 1:30)
day.1<-AddMetaData(day.1, metadata = predictions)
Idents(day.1)<-"predicted.id"

#


pdf(file.path(dir,"UMAP_clusters.pdf"))
DimPlot(day.1,reduction="umap",label=TRUE,label.size=2,repel=TRUE, cols=cols)+NoLegend()
dev.off()

pdf(file.path(dir,"UMAP_batches.pdf"))
DimPlot(day.1,reduction="umap",group.by = "orig.ident")
dev.off()

saveRDS(day.1,file.path(dir,'seurat_object.rds'))

day.1<-readRDS(file.path(dir,'seurat_object.rds'))

#select the interesting cluster
p<-DimPlot(day.1,reduction="umap",label=TRUE,label.size=2,repel=TRUE, cols=cols)+NoLegend()
select.cells <- CellSelector(plot = p)
Idents(day.1, cells = select.cells) <- "NewCells"

#refind all markers
markers<-m$find_all_markers(day.1,dir)
marker_dots<-m$marker_dotplot(day.1,markers,dir)


newcells.markers <- FindMarkers(day.1, ident.1 = "NewCells", min.pct = 0.3,logfc.threshold = 0.5,only.pos=TRUE)

day.1<-FindVariableFeatures(day.1)
day.1@graphs<-list()
day.1.loom <- as.loom(day.1, filename = file.path(dir,"day.1.loom"), verbose = TRUE)
