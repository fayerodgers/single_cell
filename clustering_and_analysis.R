#!/usr/bin/env Rscript
m <- modules::use("./SC.R")
library(argparse)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(gridExtra)

parser <- ArgumentParser()
parser$add_argument("--results_directory",help="directory for files to be written to")
parser$add_argument("--resolution",type="integer",help="Resolution for identifying clusters (default: 1)",default=1)
parser$add_argument("--metadata", help="File with sample metadata")
args <- parser$parse_args()


meta_data<-read.table(args$metadata,header=T)
dir.create(file.path(args$results_directory,args$resolution))
out_dir<-file.path(args$results_directory,args$resolution)
combined<-readRDS(file.path(args$results_directory,"combined_seurat_object.rds"))
p<-m$elbow_plot(combined,30,"pca", out_dir)
combined<-m$do_clustering(combined,30,args$resolution)

markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4")
p1<-m$QC_violins(combined, out_dir)
p2<-m$QC_scatters(combined,out_dir)
p3<-m$plot_UMAPS(combined,out_dir)
p4<-m$marker_violins(combined,out_dir,markers_to_plot)
p5<-m$split_UMAP(combined,out_dir)
markers<-m$find_all_markers(combined,out_dir)
p6<-m$marker_dotplot(combined,markers,out_dir)
p7<-m$plot_cluster_distributions(combined,out_dir,meta_data)
p8<-m$markers_feature_plot(combined,out_dir,markers_to_plot)

#find regulated genes
combined_24h<-subset(combined, subset = time ==  '24hpi')
regulated_24<-lapply(clusters,m$find_regulated_genes,combined_24h,args$results_directory,'24hpi')
names(regulated_24)<-clusters
combined_72h<-subset(combined, subset = time ==  '72hpi')
regulated_72<-lapply(clusters,m$find_regulated_genes,combined_72h,args$results_directory,'72hpi')
names(regulated_72)<-clusters
                    

