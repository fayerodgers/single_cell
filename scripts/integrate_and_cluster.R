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
library(loomR)

parser <- ArgumentParser()
parser$add_argument("--results_directory",help="directory for files to be written to")
parser$add_argument("--samples_directory",help="directory to find saved, QC'd seurat objects")
parser$add_argument("--metadata", help="File with sample metadata")
parser$add_argument("--condition", help="Integrate and cluster samples from this condition")

args <- parser$parse_args()

#for now
args$results_directory<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control'
args$samples_directory<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/QC/combined'
args$metadata<-'/Users/fr7/git_repos/single_cell/metadata/samples.txt'
args$condition<-'control'


dir <- args$results_directory
samples_dir <- args$samples_directory
meta_data<-read.table(args$metadata,header=T)

#find the IDs of the samples that we want
samples<-meta_data[which(meta_data$treatment == 'control' & meta_data$experiment == '4'),'sample_id']

#load seurat objects
objects<-list()
for (sample in samples){
  object<-readRDS(file.path(samples_dir,paste0(sample,".rds")))
  objects[[sample]]<-object
}

#normalize and find variable features
for (i in names(objects)){
  objects[[i]] <- NormalizeData(objects[[i]])
  objects[[i]] <- FindVariableFeatures(objects[[i]], selection.method = "vst", nfeatures = 2000)
}

#integrate and cluster
combined<-m$integrate_datasets(objects,40)
p<-m$elbow_plot(combined,40,"pca", dir)
combined<-m$do_clustering(combined,35,resolution=0.5,min.dist=0.05)
saveRDS(combined,file.path(dir,'seurat_object.rds'))

#plots
markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Isg15","Ace2","Ptprc")
p1<-m$QC_violins(combined, dir)
p2<-m$QC_scatters(combined,dir)
p3<-m$plot_UMAPS(combined,dir)
p4<-m$marker_violins(combined,dir,markers_to_plot)
p5<-m$markers_feature_plot(combined,dir,markers_to_plot)

#select the Chga expressing cells as their own cluster
#p<-FeaturePlot(combined,features="Chga")
#combined <- CellSelector(plot = p, object = combined, ident = "ENTEROENDOCRINE")

#rename/merge clusters

new.clusters.ids <-list()
new.clusters.ids[["5"]] <- "UNDIFFERENTIATED.1" 
new.clusters.ids[["9"]] <-  "UNDIFFERENTIATED.1"  
new.clusters.ids[["6"]] <-  "UNDIFFERENTIATED.2"
new.clusters.ids[["4"]] <- "UNDIFFERENTIATED.3"   
new.clusters.ids[["1"]] <- "ENTEROCYTE.1"
new.clusters.ids[["2"]] <- "ENTEROCYTE.1"
new.clusters.ids[["0"]] <- "ENTEROCYTE.2"
new.clusters.ids[["7"]] <-  "ENTEROCYTE.2"
new.clusters.ids[["8"]] <-  "ENTEROCYTE.3"
new.clusters.ids[["11"]] <- "ENTEROCYTE.4"
new.clusters.ids[["10"]] <- "ENTEROCYTE.ISG15"
new.clusters.ids[["12"]] <- "ENTEROCYTE.ACE2"
new.clusters.ids[["3"]] <-  "GOBLET" 
new.clusters.ids[["13"]] <- "TUFT" 
new.clusters.ids[["14"]] <- "EE"
  

combined <- RenameIdents(combined,new.clusters.ids)

DefaultAssay(combined)<-"RNA"
markers<-m$find_all_markers(combined,dir)
marker_dots<-m$marker_dotplot(combined,markers,dir)

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

#UMAP for figure
pdf(file.path(dir,"UMAP_clusters.pdf"))
DimPlot(combined,reduction="umap",label=TRUE,label.size=2,repel=TRUE, cols=cols)+NoLegend()
dev.off()

pdf(file.path(dir,"UMAP_batches.pdf"))
DimPlot(combined,reduction="umap",group.by = "orig.ident")
dev.off()

DefaultAssay(combined)<-"RNA"
markers<-m$find_all_markers(combined,dir)
marker_dots<-m$marker_dotplot(combined,markers,dir)
saveRDS(combined,file.path(dir,'seurat_object.rds'))

#also save the object in loom format for scanpy
#need to rerun findvariable features to export as loom
combined<-FindVariableFeatures(combined)
combined.loom <- as.loom(combined, filename = file.path(dir,"control.loom"), verbose = TRUE)

