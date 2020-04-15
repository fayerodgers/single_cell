#!/usr/bin/env Rscript
#.libPaths("/nfs/users/nfs_f/fr7/anaconda2/envs/r_env/lib/R/library")
m <- modules::use("/Users/fr7/git_repos/single_cell/SC.R")
library(argparse)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(gridExtra)
library(stats)
library(DoubletDecon)

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
args$results_directory<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL'
#args$mito_cutoff<-0
#args$ncells_cutoff<-10

#find samples and meta data
samples<-scan(args$samples,what=character())
meta_data<-read.table(args$metadata,header=T)

#make output directories
dir<-file.path(args$results_directory,'QC')
dir.create(dir)
for (sample in samples){
  dir.create(file.path(dir,sample))
}

#QC:get samples into seurat_objects
seurat_objects<-lapply(samples,m$normalize_data,
                       args$cellranger_version,
                       meta_data,
                       args$variable_features,
                       args$ncells_cutoff,
                       args$mito_cutoff,
                       args$nfeatures_cutoff,
                       args$annotation_version)
names(seurat_objects)<-samples

#plot UMIs for each sample
umi_frequency_plots<-lapply(names(seurat_objects),m$umi_frequency_plot,seurat_objects,dir,28000) #try various to get the right threshold for each sample

#also plot %mito and ngenes for each sample 
mito_frequency_plots<-lapply(names(seurat_objects),m$mito_frequency_plot,seurat_objects,dir)
all<-m$print_combined(mito_frequency_plots,dir,'mito_frequency_plots')

ngenes_plots<-lapply(names(seurat_objects),m$genes_frequency_plot,seurat_objects,dir)
all<-m$print_combined(ngenes_plots,dir,'ngenes_plots')


#Inspect UMI plots and select a cutoff for each sample
#Apply UMI filters
umi_filters<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/UMI_filters.txt',header=T)
#plot all selected filters
all_p<-list()
for (sample in names(seurat_objects)){
  threshold<-umi_filters[which(umi_filters$sample == sample),'umi_threshold' ]
  p<-m$umi_frequency_plot(sample,seurat_objects,dir, threshold)
  all_p[[sample]]<-p
}
all<-m$print_combined(all_p,dir,'umi_with_thresholds')

#apply the filters
umi_filtered<-list()
umi_filtered$`4672STDY8112878`<-subset(seurat_objects$`4672STDY8112878`,  subset = nCount_RNA > 19000)
umi_filtered$`4672STDY8112879`<-subset(seurat_objects$`4672STDY8112879`,  subset = nCount_RNA > 25000)
umi_filtered$`4672STDY8112880`<-subset(seurat_objects$`4672STDY8112880`,  subset = nCount_RNA > 22500)
umi_filtered$`4672STDY8112881`<-subset(seurat_objects$`4672STDY8112881`,  subset = nCount_RNA > 26000)
umi_filtered$`4672STDY8112882`<-subset(seurat_objects$`4672STDY8112882`,  subset = nCount_RNA > 25000)
umi_filtered$`4672STDY8112883`<-subset(seurat_objects$`4672STDY8112883`,  subset = nCount_RNA > 23000)
umi_filtered$`4672STDY8112884`<-subset(seurat_objects$`4672STDY8112884`,  subset = nCount_RNA > 20000)
umi_filtered$`4672STDY8112885`<-subset(seurat_objects$`4672STDY8112885`,  subset = nCount_RNA > 23000)
umi_filtered$`4672STDY8112974`<-subset(seurat_objects$`4672STDY8112974`,  subset = nCount_RNA > 21500)
umi_filtered$`4672STDY8112975`<-subset(seurat_objects$`4672STDY8112975`,  subset = nCount_RNA > 22000)
umi_filtered$`4672STDY8112976`<-subset(seurat_objects$`4672STDY8112976`,  subset = nCount_RNA > 24000)
umi_filtered$`4672STDY8112977`<-subset(seurat_objects$`4672STDY8112977`,  subset = nCount_RNA > 26500)
umi_filtered$`4672STDY8113070`<-subset(seurat_objects$`4672STDY8113070`,  subset = nCount_RNA > 28000)
umi_filtered$`4672STDY8112979`<-subset(seurat_objects$`4672STDY8112979`,  subset = nCount_RNA > 24000)
umi_filtered$`4672STDY8112981`<-subset(seurat_objects$`4672STDY8112981`,  subset = nCount_RNA > 25000)
umi_filtered$`4672STDY8112980`<-subset(seurat_objects$`4672STDY8112980`,  subset = nCount_RNA > 23000)

#new directory for the UMI filtered samples
dir<-file.path(dir,'umi_filtered')
dir.create(dir)
for (sample in samples){
  dir.create(file.path(dir,sample))
}

#redo all plots post umi filtering
umi_frequency_plots<-lapply(names(umi_filtered),m$umi_frequency_plot,umi_filtered,dir,0)
all<-m$print_combined(umi_frequency_plots,dir,'UMI_frequency')

mito_frequency_plots<-lapply(names(umi_filtered),m$mito_frequency_plot,umi_filtered,dir)
all<-m$print_combined(mito_frequency_plots,dir,'mito_frequency')

genes_frequency_plots<-lapply(names(umi_filtered),m$genes_frequency_plot,umi_filtered,dir)
all<-m$print_combined(genes_frequency_plots,dir,'genes_frequency')


#apply a % mito cutoff too
sub<-function(x) subset(x,subset=percent.mito<=30)
mito_filtered<-lapply(umi_filtered,sub)

#new directory for the mito filtered samples
dir<-file.path(dir,'mito_filtered')
dir.create(dir)
for (sample in samples){
  dir.create(file.path(dir,sample))
}

#redo all plots post mito filtering
umi_frequency_plots<-lapply(names(umi_filtered),m$umi_frequency_plot,umi_filtered,dir,0)
all<-m$print_combined(umi_frequency_plots,dir,'UMI_frequency')

mito_frequency_plots<-lapply(names(umi_filtered),m$mito_frequency_plot,umi_filtered,dir)
all<-m$print_combined(mito_frequency_plots,dir,'mito_frequency')

genes_frequency_plots<-lapply(names(umi_filtered),m$genes_frequency_plot,umi_filtered,dir)
all<-m$print_combined(genes_frequency_plots,dir,'genes_frequency')


#run PCA and cluster each sample
#scale data and run pca
mito_filtered<-lapply(mito_filtered,m$run_pca,args$dimensions_to_assess)

#elbow plots to select no. of PCs to use in clustering
elbow_plots<-c()
for (i in names(mito_filtered)){
  out_dir<-file.path(dir,i)
  p<-m$elbow_plot(mito_filtered[[i]],args$dimensions_to_assess,"pca", out_dir)
  elbow_plots<-c(elbow_plots,p)
}

#clustering
mito_filtered<-lapply(mito_filtered,m$do_clustering,25,1)
markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Isg15")
for (i in names(mito_filtered)){
  out_dir<-file.path(dir,i)
  p1<-m$QC_violins(mito_filtered[[i]], out_dir)
  p2<-m$QC_scatters(mito_filtered[[i]],out_dir)
  p3<-m$plot_UMAPS(mito_filtered[[i]],out_dir)
  p4<-m$marker_violins(mito_filtered[[i]],out_dir,markers_to_plot)
# markers<-m$find_all_markers(mito_filtered[[i]],out_dir)
  saveRDS(mito_filtered[[i]],file.path(out_dir,"seurat_object.rds"))
}

#new directory for doublet filtering
dir<-file.path(dir,'doublet_filtering')
dir.create(dir)
for (sample in samples){
  dir.create(file.path(dir,sample))
}

#Doublet decon
doublet_filtered<-list()
ndoublets<-c()
for (i in names(mito_filtered)){
  out_dir<-file.path(dir,i)
  files<-m$prepare_for_doublet_decon(mito_filtered[[i]],out_dir,i)
}

#starting here today- reload seurat objects
mito_filtered<-list()
dir<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/QC/umi_filtered/mito_filtered'
for (sample in samples){
  i<-readRDS(file.path(dir,sample,'seurat_object.rds'))
  mito_filtered[[sample]]<-i
}
dir<-file.path(dir,'doublet_filtering')
###

#run this for a few different values of rhop and decide on the best value for each sample

for (i in names(mito_filtered)){
  out_dir<-file.path(dir,i)
  expression_file<-read.table(file.path(out_dir,paste0(i,"_expression")),row.names=1)
  groups_file<-read.table(file.path(out_dir,paste0(i,"_groups")),header=FALSE,row.names = 1)
  results<-m$run_doublet_decon(expression_file,groups_file,i,out_dir,0.75,TRUE)
}

doublet_filters<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/doublet_filters.txt',header=T)
results<-list()
for (i in names(mito_filtered)){
  out_dir<-file.path(dir,i)
  expression_file<-read.table(file.path(out_dir,paste0(i,"_expression")),row.names=1)
  groups_file<-read.table(file.path(out_dir,paste0(i,"_groups")),header=FALSE,row.names = 1)
  filter<-doublet_filters[which(doublet_filters$sample == i),'rhop']
  results[[i]]<-m$run_doublet_decon(expression_file,groups_file,i,out_dir,filter,TRUE)
}

#extract doublets
after_doublets<-list()

for (i in names(results)){
  doublets<-row.names(results[[i]]$DRS_doublet_table[which(results[[i]]$DRS_doublet_table$isADoublet == 'TRUE'),])
  print(c(i,length(doublets)))
  after_doublets[[i]]<-subset(mito_filtered[[i]],cells=doublets,invert=TRUE)
}

after_doublets<-lapply(after_doublets,m$run_pca,args$dimensions_to_assess)
after_doublets<-lapply(after_doublets,m$do_clustering,args$dimensions_to_analyse,1)
for (i in names(after_doublets)){
  out_dir<-file.path(dir,i)
  p1<-m$QC_violins(after_doublets[[i]], out_dir)
  p2<-m$QC_scatters(after_doublets[[i]],out_dir)
  p3<-m$plot_UMAPS(after_doublets[[i]],out_dir)
  p4<-m$marker_violins(after_doublets[[i]],out_dir,markers_to_plot)
  #  markers<-m$find_all_markers(postQC[[i]],out_dir)
  saveRDS(after_doublets[[i]],file.path(out_dir,"after_doublets.rds"))
}
