#!/usr/bin/env Rscript
#.libPaths("/nfs/users/nfs_f/fr7/anaconda2/envs/r_env/lib/R/library")
m <- modules::use("./SC.R")
library(argparse)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(gridExtra)

parser <- ArgumentParser()
parser$add_argument("--results_directory",help="directory for files to be written to")
parser$add_argument("--samples", help="File with a list of sample IDs to analyse, one per line")
parser$add_argument("--metadata", help="File with sample metadata")
parser$add_argument("--cellranger_version", help="cellranger131, cellranger211 or cellranger302")
parser$add_argument("--annotation_version", help="mouse reference annotation used in Cellranger",default='mm10-3_0_0')
parser$add_argument("--mito_cutoff", type="integer", help="Filter out cells with % reads aligning to mitochondrial genes greater than this (default: no cutoff)",default=100)
parser$add_argument("--nfeatures_cutoff", type="integer",help="Filter out cells with fewer features than this (default: no filter)", default=0)
parser$add_argument("--ncells_cutoff",type="integer",help="Only include features (genes) that are expressed in more than this number of cells (default: no filter)",default=0)
parser$add_argument("--variable_features",type="integer",help="Take top x variable genes (default: 2000)",default=2000)
parser$add_argument("--dimensions_to_assess", type="integer",help="Number of prinicpal components to calculate (default: 30)",default=30)
parser$add_argument("--dimensions_to_analyse",type="integer",help="Number of prinicpal components to include in UMAP and neighbour finding (default: 25)", default=25)
parser$add_argument("--resolution",type="integer",help="Resolution for identifying clusters (default: 1)",default=1)

args <- parser$parse_args()

#enable paralellisation
#plan("multiprocess", workers = 4)
#options(future.globals.maxSize = 9437184000)

#for now
args$cellranger_version<-'cellranger302'
args$metadata<-'/Users/fr7/git_repos/single_cell/metadata/samples.txt'
args$samples<-'/Users/fr7/git_repos/single_cell/experiment_4/samples_to_analyse.txt'
args$results_directory<-'/Users/fr7/git_repos/single_cell/experiment_4'
#args$mito_cutoff<-0
args$ncells_cutoff<-10

evaluate<-function(text){
  x<-eval(parse(text=text))
  return(x)
}


samples<-scan(args$samples,what=character())
meta_data<-read.table(args$metadata,header=T)
for (sample in samples){
  dir.create(file.path(args$results_directory,sample))
}


#TODO write a readme into the data dir with parameters used for the analysis

#QC: for each sample, identify clusters of dying cells and remove them.
seurat_objects<-lapply(samples,m$normalize_data,
                       args$cellranger_version,
                       meta_data,
                       args$variable_features,
                       args$ncells_cutoff,
                       args$mito_cutoff,
                       args$nfeatures_cutoff,
                       args$annotation_version)
names(seurat_objects)<-samples

#plots of variable features - check we capture the highly variable genes with 2000 features.
for (i in names(seurat_objects)){
  out_dir<-file.path(args$results_directory,i)
  p<-m$variable_feature_plot(seurat_objects[[i]],out_dir)
}


#scale data and run pca
seurat_objects<-lapply(seurat_objects,m$run_pca,args$dimensions_to_assess)

#elbow plots to select no. of PCs to use in clustering
elbow_plots<-c()
for (i in names(seurat_objects)){
  out_dir<-file.path(args$results_directory,i)
  p<-m$elbow_plot(seurat_objects[[i]],args$dimensions_to_assess,"pca", out_dir)
  elbow_plots<-c(elbow_plots,p)
}

#clustering
seurat_objects<-lapply(seurat_objects,m$do_clustering,args$dimensions_to_analyse,5)

markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4")

for (i in names(seurat_objects)){
  out_dir<-file.path(args$results_directory,i)
  p1<-m$QC_violins(seurat_objects[[i]], out_dir)
  p2<-m$QC_scatters(seurat_objects[[i]],out_dir)
  p3<-m$plot_UMAPS(seurat_objects[[i]],out_dir)
  p4<-m$marker_violins(seurat_objects[[i]],out_dir,markers_to_plot)
  markers<-m$find_all_markers(seurat_objects[[i]],out_dir)
  saveRDS(seurat_objects[[i]],file.path(out_dir,"seurat_object_unfiltered.rds"))
}

#inspect each sample and select clusters to remove and doublet thresholds
clusters_to_remove<-read.table(file.path(args$results_directory,"clusters_to_remove.txt"),header=TRUE)
clusters_to_remove$cluster<-(clusters_to_remove$cluster)

#doublet filtering: can not pass variables to subset, so have to do this bit manually.
seurat_objects$`4672STDY8112878`<-subset(seurat_objects$`4672STDY8112878`,  subset = nCount_RNA < 90000)
seurat_objects$`4672STDY8112879`<-subset(seurat_objects$`4672STDY8112879`,  subset = nCount_RNA < 80000)
seurat_objects$`4672STDY8112880`<-subset(seurat_objects$`4672STDY8112880`,  subset = nCount_RNA < 100000)
seurat_objects$`4672STDY8112881`<-subset(seurat_objects$`4672STDY8112881`,  subset = nCount_RNA < 125000)
seurat_objects$`4672STDY8112882`<-subset(seurat_objects$`4672STDY8112882`,  subset = nCount_RNA < 80000)
seurat_objects$`4672STDY8112883`<-subset(seurat_objects$`4672STDY8112883`,  subset = nCount_RNA < 100000)
seurat_objects$`4672STDY8112884`<-subset(seurat_objects$`4672STDY8112884`,  subset = nCount_RNA < 80000)
seurat_objects$`4672STDY8112885`<-subset(seurat_objects$`4672STDY8112885`,  subset = nCount_RNA < 100000)
seurat_objects$`4672STDY8112974`<-subset(seurat_objects$`4672STDY8112974`,  subset = nCount_RNA < 100000)
seurat_objects$`4672STDY8112975`<-subset(seurat_objects$`4672STDY8112975`,  subset = nCount_RNA < 100000)
seurat_objects$`4672STDY8112976`<-subset(seurat_objects$`4672STDY8112976`,  subset = nCount_RNA < 100000)
seurat_objects$`4672STDY8112977`<-subset(seurat_objects$`4672STDY8112977`,  subset = nCount_RNA < 100000)
seurat_objects$`4672STDY8113070`<-subset(seurat_objects$`4672STDY8113070`,  subset = nCount_RNA < 100000)
seurat_objects$`4672STDY8112979`<-subset(seurat_objects$`4672STDY8112979`,  subset = nCount_RNA < 100000)
seurat_objects$`4672STDY8112981`<-subset(seurat_objects$`4672STDY8112981`,  subset = nCount_RNA < 80000)
seurat_objects$`4672STDY8112980`<-subset(seurat_objects$`4672STDY8112980`,  subset = nCount_RNA < 125000)

#rescale and find variable features
for (i in names(seurat_objects)){
  out_dir<-file.path(args$results_directory,i)
  seurat_objects[[i]]<-NormalizeData(seurat_objects[[i]])
  seurat_objects[[i]] <- FindVariableFeatures(seurat_objects[[i]], selection.method = "vst", nfeatures = 2000)
  saveRDS(seurat_objects[[i]],file.path(out_dir,"seurat_object_postQC.rds"))
}

#integrate datasets

filtered<-c()
for (i in names(samples)){
  temp<-readRDS(file.path(args$results_directory,i,"seurat_object_postQC.rds"))
  filtered<-c(filtered,temp)
}
names(filtered)<-names(samples)

combined<-m$integrate_datasets(seurat_objects,30)







#replot and save objects post filtering
for (i in names(seurat_objects)){
  dir.create(file.path(args$results_directory,i,"postQC"))
  out_dir<-file.path(args$results_directory,i,"postQC")
  p1<-m$QC_violins(seurat_objects[[i]], out_dir)
  p2<-m$QC_scatters(seurat_objects[[i]],out_dir)
  p3<-m$plot_UMAPS(seurat_objects[[i]],out_dir)
  p4<-m$marker_violins(seurat_objects[[i]],out_dir,markers_to_plot)
#  markers<-m$find_all_markers(seurat_objects[[i]],out_dir)
  saveRDS(seurat_objects[[i]],file.path(out_dir,"seurat_object_filtered.rds"))
}




#Integrate datasets into one object
combined<-m$integrate_datasets(seurat_objects,args$dimensions_to_assess)
combined<-m$do_clustering(combined,args$dimensions_to_analyse,args$resolution)


DefaultAssay(combined) <- "RNA"
saveRDS(combined,paste0(args$results_directory,"/seurat_object.rds"))

#plots to select which clusters to remove
mito.vln<-VlnPlot(combined, features = c("percent.mito"), cols=rep("red",length(levels(combined$seurat_clusters))), pt.size = 0,do.return=T)
mito.vln<-mito.vln + theme(axis.text.x=element_text(size=8)) + labs(x="",y="% mitochondrial genes")
nfeatures.vln<-VlnPlot(combined, features = c("nFeature_RNA"), cols=rep("forestgreen",length(levels(combined$seurat_clusters))), pt.size = 0,do.return=T)
nfeatures.vln<-nfeatures.vln + theme(axis.text.x=element_text(size=8)) + labs(x="",y="Number of genes")
counts.vln<-VlnPlot(combined, features = c("nCount_RNA"), cols=rep("tan1",length(levels(combined$seurat_clusters))), pt.size = 0,do.return=T)
counts.vln<-counts.vln + theme(axis.text.x=element_text(size=8)) + labs(x="",y="Number of molecules")
p<-DimPlot(combined,reduction="umap",split.by = "time.status")
dim.plots<-m$plot_UMAPS(combined)
markers<-FeaturePlot(combined, features = c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Il33","Ly6a"), min.cutoff="q9",pt.size = 0.1)

combined <- subset(combined, subset = nCount_RNA < 150000)
combined <- subset(combined, idents = c('1','6','11','29'), invert = TRUE)
#re-do clustering
elbow<-ElbowPlot(combined,ndims=args$dimensions_to_assess)
args$dimensions_to_analyse <- 20
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst",nfeatures = 2000)
combined <- m$run_pca(combined,20)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 1)

#

print_plot<-function(dir,file_name,plot){
  pdf(paste0(dir,file_name))
  print(plot)
  dev.off()
}



#Printing plots
pdf(paste0(args$results_directory,"/elbowplot.pdf"))
ElbowPlot(combined,ndims=args$dimensions_to_assess)
dev.off()

mito.vln<-VlnPlot(combined, features = c("percent.mito"), cols=rep("red",length(samples)), pt.size = 0,group.by='orig.ident',do.return=T)
mito.vln<-mito.vln + theme(axis.text.x=element_text(size=8)) + labs(x="",y="% mitochondrial genes", title = args$cellranger_version)

nfeatures.vln<-VlnPlot(combined, features = c("nFeature_RNA"), cols=rep("forestgreen",length(samples)), pt.size = 0,group.by='orig.ident',do.return=T)
nfeatures.vln<-nfeatures.vln + theme(axis.text.x=element_text(size=8)) + labs(x="",y="Number of genes", title = args$cellranger_version)

pdf(paste0(args$results_directory,"/percentmito.pdf"),12,6)
print(mito.vln)
dev.off()

pdf(paste0(args$results_directory,"/nfeatures.pdf"),12,6)
print(nfeatures.vln)
dev.off()

pdf(paste0(args$results_directory,"/umaps.pdf"),20,20)
print(dim.plots)
dev.off()

pdf(paste0(args$results_directory,"/infection_umap.pdf"),20,20)
print(p)
dev.off()

pdf(paste0(args$results_directory,"/markers.pdf"),20,20)
print(markers)
dev.off()

#find cluster markers
clusters<-levels(combined$seurat_clusters)
markers<-lapply(clusters,m$find_markers,combined,args$results_directory)
names(markers)<-clusters

#dotplot of top 3 marker genes for each cluster
DefaultAssay(combined) <- "RNA"
markers_to_plot<-c()
for (cluster in clusters){
  x<-evaluate(paste0("markers$'",cluster,"'"))
  markers_to_plot<-c(markers_to_plot, rownames(x)[0:3] )
}
markers_to_plot<-unique(markers_to_plot)
p<-DotPlot(combined, features = rev(markers_to_plot), dot.scale = 4)
p<-p+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
pdf(paste0(args$results_directory,'/markers_dots.pdf'),15,10)
print(p)
dev.off()

#plot cell numbers per cluster
proportions<-as.data.frame(prop.table(table(Idents(combined), combined$orig.ident)))
names(proportions)<-c("cluster","sample","proportion")
proportions<-merge(proportions,meta_data,by.x = "sample",by.y="sample_id")
proportions$time.status<-paste0(proportions$time,"-",proportions$treatment)
proportion_plots<-lapply(clusters, m$plot_proportions, proportions,"time.status")
all_plots<-plot_grid(plotlist=proportion_plots,ncol=4)
pdf(paste0(args$results_directory,'/proportion_plots.pdf'),20,20)
print(all_plots)
dev.off()

#find regulated genes- TODO- make it so that we're comparing within a time point.
combined_24h<-subset(combined,idents="24hpi")
Idents(combined_24h)<-'seurat_clusters'
combined_72h<-subset(combined,idents="72hpi")
Idents(combined_72h)<-'seurat_clusters'
regulated_24<-lapply(clusters,m$find_regulated_genes,combined_24h,args$results_directory,'24hpi')
names(regulated_24)<-clusters
regulated_72<-lapply(clusters,m$find_regulated_genes,combined_72h,args$results_directory,'72hpi')
names(regulated_72)<-clusters



