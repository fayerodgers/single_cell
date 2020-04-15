#!/usr/bin/env Rscript
#.libPaths("/nfs/users/nfs_f/fr7/anaconda2/envs/r_env/lib/R/library")
m <- modules::use("./SC.R")
library(argparse)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(gridExtra)
library(stats)

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

#evaluate<-function(text){
#  x<-eval(parse(text=text))
#  return(x)
#}


samples<-scan(args$samples,what=character())
meta_data<-read.table(args$metadata,header=T)
for (sample in samples){
  dir.create(file.path(args$results_directory,sample))
}


#TODO write a readme into the data dir with parameters used for the analysis

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
umi_frequency_plots<-lapply(names(seurat_objects),m$umi_frequency_plot,seurat_objects,args$results_directory,28000)
all_umi_frequency_plots<-plot_grid(plotlist=umi_frequency_plots,ncol=4,nrow=4)
pdf(file.path(args$results_directory,'all_umi_frequency_plots_preQC.pdf'),30,30)
print(all_umi_frequency_plots)
dev.off()

#Plot % mito for each sample
mito_frequency_plots<-lapply(names(seurat_objects),m$mito_frequency_plot,seurat_objects,args$results_directory)
all_mito_frequency_plots<-plot_grid(plotlist=mito_frequency_plots,ncol=4,nrow=4)
pdf(file.path(args$results_directory,'all_mito_frequency_plots_preQC.pdf'),30,30)
print(all_mito_frequency_plots)
dev.off()

#Plot n genes for each sample
genes_frequency_plots<-lapply(names(seurat_objects),m$genes_frequency_plot,seurat_objects,args$results_directory)
all_genes_frequency_plots<-plot_grid(plotlist=genes_frequency_plots,ncol=4,nrow=4)
pdf(file.path(args$results_directory,'all_genes_frequency_plots_preQC.pdf'),30,30)
print(all_genes_frequency_plots)
dev.off()

#Inspect UMI plots and select a cutoff for each sample
#Apply UMI filters
umi_filters<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/UMI_filters.txt',header=T)
#plot all selected filters
all_p<-list()
for (sample in names(seurat_objects)){
  threshold<-umi_filters[which(umi_filters$sample == sample),'umi_threshold' ]
  p<-m$umi_frequency_plot(sample,seurat_objects,args$results_directory, threshold)
  all_p[[sample]]<-p
}
all_threshold_plots<-plot_grid(plotlist=all_p,ncol=4,nrow=4)
pdf(file.path(args$results_directory,'all_umi_frequency_plots_thresholds.pdf'),30,30)
print(all_threshold_plots)
dev.off()

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
  
#replot
dir.create(file.path(args$results_directory,'umi_filter'))
for (sample in samples){
  dir.create(file.path(args$results_directory,'umi_filter',sample))
}
dir<-file.path(args$results_directory,'umi_filter')
#plot UMIs for each sample
umi_frequency_plots<-lapply(names(umi_filtered),m$umi_frequency_plot,umi_filtered,dir,0)
all_umi_frequency_plots<-plot_grid(plotlist=umi_frequency_plots,ncol=4,nrow=4)
pdf(file.path(args$results_directory,'all_umi_frequency_plots_post_UMIfilter.pdf'),30,30)
print(all_umi_frequency_plots)
dev.off()

#Plot % mito for each sample
mito_frequency_plots<-lapply(names(umi_filtered),m$mito_frequency_plot,umi_filtered,dir)
all_mito_frequency_plots<-plot_grid(plotlist=mito_frequency_plots,ncol=4,nrow=4)
pdf(file.path(args$results_directory,'all_mito_frequency_plots_post_UMIfilter.pdf'),30,30)
print(all_mito_frequency_plots)
dev.off()

#Plot n genes for each sample
genes_frequency_plots<-lapply(names(umi_filtered),m$genes_frequency_plot,umi_filtered,dir)
all_genes_frequency_plots<-plot_grid(plotlist=genes_frequency_plots,ncol=4,nrow=4)
pdf(file.path(args$results_directory,'all_genes_frequency_plots_post_UMIfilter.pdf'),30,30)
print(all_genes_frequency_plots)
dev.off()


#apply a % mito cutoff too
sub<-function(x) subset(x,subset=percent.mito<=30)
postQC<-lapply(umi_filtered,sub)



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

#maybe
filtered<-c()
for (i in names(samples)){
  temp<-readRDS(file.path(args$results_directory,i,"seurat_object_postQC.rds"))
  filtered<-c(filtered,temp)
}
names(filtered)<-names(samples)

#integrate datasets
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

combined<-m$do_clustering(combined,args$dimensions_to_analyse,args$resolution)
saveRDS(combined,file.path(out_dir,"seurat_object_filtered.rds"))

#maybe
combined<-readRDS('/Users/fr7/git_repos/single_cell/experiment_4/QC_bespoke/1/local/combined_clustered.rds')

print_plot<-function(dir,file_name,plot){
  pdf(paste0(dir,file_name))
  print(plot)
  dev.off()
}

#find cluster markers
markers<-m$find_all_markers(combined,args$results_directory)
#maybe
markers<-read.table(args$results_directory,header=T,stringsAsFactors=FALSE)

#dotplot of top 3 marker genes for each cluster
marker_dots<-m$marker_dotplot(combined,markers,args$results_directory)

#rename clusters based on marker genes
cluster_ids<-read.table(file.path(args$results_directory,"cluster_ids.txt"),header=TRUE,stringsAsFactors=FALSE)
new.cluster.ids<-cluster_ids$id
names(new.cluster.ids)<-cluster_ids$cluster
combined <- RenameIdents(combined,new.cluster.ids)

#plot cell numbers per cluster
m$plot_cluster_distributions(combined,args$results_directory,meta_data)

#retrieve a table of 
numbers<-as.data.frame(table(Idents(combined), combined$orig.ident))
colnames(numbers)<-c("cluster","sample","cells")
totals<-as.data.frame(table(combined$orig.ident))
colnames(totals)<-c("sample","total_cells")
clusters<-levels(numbers$cluster)
df<-merge(numbers,meta_data,by.x="sample",by.y="sample_id")
df<-merge(df,totals,by="sample")
df$cluster<-as.character(df$cluster)

for (cluster in clusters){
  temp<-df[which(df$cluster==cluster),]
  temp24<-temp[which(temp$time=="24hpi"),]
  model24<-glm(cells ~ treatment + offset(log(total_cells)),family=poisson(link=log),data=temp24)
  temp72<-temp[which(temp$time=="72hpi"),]
  model72<-glm(cells ~ treatment + offset(log(total_cells)),family=poisson(link=log),data=temp72)
  sink(file.path(args$results_directory,paste0(cluster,".24.glm.txt")))
  print(temp24)
  print(summary(model24))
  sink()
  sink(file.path(args$results_directory,paste0(cluster,".72.glm.txt")))
  print(temp72)
  print(summary(model72))
  sink()  
}




test2_24<-test2[which(test2$time=="24hpi"),]
test2_72<-test2[which(test2$time=="72hpi"),]
model24<-glm(cells ~ treatment + offset(log(total_cells)),family=poisson(link=log),data=test2_24)
model72<-glm(cells ~ treatment + offset(log(total_cells)), data=test2_72,family=poisson(link=log))
model72_x<-glm(cells ~ 1 + offset(log(total_cells)), data=test2_72)

#find regulated genes
combined_24h<-subset(combined,idents="24hpi")
Idents(combined_24h)<-'seurat_clusters'
combined_72h<-subset(combined,idents="72hpi")
Idents(combined_72h)<-'seurat_clusters'
regulated_24<-lapply(clusters,m$find_regulated_genes,combined_24h,args$results_directory,'24hpi')
names(regulated_24)<-clusters
regulated_72<-lapply(clusters,m$find_regulated_genes,combined_72h,args$results_directory,'72hpi')
names(regulated_72)<-clusters



