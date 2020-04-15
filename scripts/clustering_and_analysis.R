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
parser$add_argument("--resolution",help="Resolution for identifying clusters (default: 1)",default=1)
parser$add_argument("--metadata", help="File with sample metadata")
args <- parser$parse_args()

args$results_directory<-'~/git_repos/single_cell/experiment_4/QC_bespoke'
args$resolution<-1
args$metadata<-'~/git_repos/single_cell/metadata/samples.txt'

resolution<-as.numeric(args$resolution)
meta_data<-read.table(args$metadata,header=T)
dir.create(file.path(args$results_directory,args$resolution,'local'))
out_dir<-file.path(args$results_directory,args$resolution,'local')
combined<-readRDS(file.path(args$results_directory,"combined_seurat_object.rds"))
p<-m$elbow_plot(combined,30,"pca", out_dir)
combined<-m$do_clustering(combined,30,resolution)
saveRDS(combined,file.path(out_dir,"combined_clustered.rds"))


cluster_ids<-read.table('~/git_repos/single_cell/experiment_4/QC_bespoke/1/local/cluster_ids.txt',header=T)
new.cluster.ids <-as.character(cluster_ids$id)
names(new.cluster.ids)<-levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)

markers_to_plot<-c("Aqp8","Krt20","Muc2","Reg4","Chga","Lgr5","Hopx","Wfdc2","Dclk1","Cdk4")
p1<-m$QC_violins(combined, out_dir)
p2<-m$QC_scatters(combined,out_dir)
p3<-m$plot_UMAPS(combined,out_dir)
p4<-m$marker_violins(combined,out_dir,markers_to_plot)
p5<-m$split_UMAP(combined,out_dir)
#markers<-m$find_all_markers(combined,out_dir)
markers<-read.table('~/git_repos/single_cell/experiment_4/QC_bespoke/markers.txt',header=T)
markers$gene<-as.character(markers$gene)
p6<-m$marker_dotplot(combined,markers,out_dir)
p7<-m$plot_cluster_distributions(combined,out_dir,meta_data)
p8<-m$markers_feature_plot(combined,out_dir,markers_to_plot)


#find regulated genes
clusters<-levels(combined)
out_dir<-file.path(args$results_directory,args$resolution,'local','24h')
dir.create(out_dir)
combined_24h<-subset(combined, subset = time ==  '24hpi')
regulated_24<-lapply(clusters,m$find_regulated_genes,combined_24h,out_dir,'24hpi')
names(regulated_24)<-clusters
p9<-m$plot_all_regulated(combined_24h,out_dir,regulated_24)

combined_72h<-subset(combined, subset = time ==  '72hpi')
regulated_72<-lapply(clusters,m$find_regulated_genes,combined_72h,out_dir,'72hpi')
names(regulated_72)<-clusters
out_dir<-file.path(args$results_directory,args$resolution,'local','72h')
dir.create(out_dir)
p10<-m$plot_all_regulated(combined_72h,out_dir,regulated_72)

richard<-c("Per1","Per2","Clock","Arntl","Cry1","Cry2","Timeless","Muc2","Reg4")
DefaultAssay(combined)<"RNA"
for (gene in richard){
  p<-FeaturePlot(combined,gene,min.cutoff="q9",pt.size = 0.1)
  pdf(paste0(gene,".pdf"))
  print(p)
  dev.off()
}
