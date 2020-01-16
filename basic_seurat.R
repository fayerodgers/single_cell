#!/usr/bin/env Rscript
.libPaths("/nfs/users/nfs_f/fr7/anaconda2/envs/r_env/lib/R/library")
m <- modules::use("./SC.R")
library(argparse)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)


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
#args$cellranger_version<-'cellranger302'
#args$metadata<-'/Users/fr7/git_repos/single_cell/metadata/samples.txt'
#args$samples<-'/Users/fr7/git_repos/single_cell/experiment_4/samples_to_analyse.txt'
#args$results_directory<-'/Users/fr7/git_repos/single_cell/experiment_4'
#args$mito_cutoff<-0
#args$ncells_cutoff<-3

evaluate<-function(text){
  x<-eval(parse(text=text))
  return(x)
}


samples<-scan(args$samples,what=character())
meta_data<-read.table(args$metadata,header=T)
for (sample in samples){
  dir.create(sample)
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
#don't do this for now
seurat_objects<-lapply(seurat_objects,m$run_pca,args$dimensions_to_assess)
#seurat_objects<-lapply(seurat_objects,m$do_clustering,args$dimensions_to_analyse,args$resolution)
QC_plots<-lapply(seurat_objects,m$QC_violins)
#markers<-lapply(seurat_objects,m$find_all_markers)
names(seurat_objects)<-samples

#inspect each library and identify clusters to be removed
#clusters_to_remove<-read.table('/Users/fr7/git_repos/single_cell/experiment_4_mt_filter/clusters_to_remove.txt',header=T)

#filtered_seurat_objects<-c()

#seurat_objects<-lapply(seurat_objects,m$remove_mito_clusters,clusters_to_remove)


#Integrate datasets into one object
combined<-m$integrate_datasets(seurat_objects,args$dimensions_to_assess)
combined<-m$do_clustering(combined,args$dimensions_to_analyse,args$resolution)
saveRDS(combined,paste0(args$results_directory,"/seurat_object.rds"))
#make new slot for time~infection status
combined[["time.status"]]<-paste0(combined$time,".",combined$infection_status)
saveRDS(combined,paste0(args$results_directory,"/seurat_object.rds"))

p<-DimPlot(combined,reduction="umap",split.by = "time.status")
dim.plots<-m$plot_UMAPS(combined)
markers<-FeaturePlot(combined, features = c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Il33","Ly6a"), min.cutoff="q9",pt.size = 0.1)

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



