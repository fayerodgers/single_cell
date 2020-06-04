#!/usr/bin/env Rscript
#.libPaths("/nfs/users/nfs_f/fr7/anaconda2/envs/r_env/lib/R/library")
m <- modules::use("/Users/fr7/git_repos/single_cell/scripts/SC.R")
library(argparse)
library(Seurat)
library(ggplot2)

control<-readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/seurat_object.rds')
d1<-readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/day1_classify/seurat_object.rds')
d3<-readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/day3_classify/seurat_object.rds')

dir<-'~/temp'
  
control.markers<-m$find_all_markers(control,dir)
d1.markers<-m$find_all_markers(d1,dir)
d3.markers<-m$find_all_markers(d3,dir)

all.markers<-rbind(control.markers,d1.markers,d3.markers)

my.markers<-unique(all.markers[which((all.markers$cluster == "ENTEROCYTE.ISG15" |
                              all.markers$cluster == "ENTEROCYTE.ACE2"  |
                              all.markers$cluster == "NewCells"  ) &
                              all.markers$p_val_adj < 0.01
                              ), "gene"])

DefaultAssay(control)<-"integrated"

pdf(file.path(dir,'control.heatmap.pdf'),30,10)
DoHeatmap(control,features=my.markers)
dev.off()

pdf(file.path(dir,'d1.heatmap.pdf'),30,10)
DoHeatmap(d1,features=my.markers)
dev.off()

pdf(file.path(dir,'d3.heatmap.pdf'),30,10)
DoHeatmap(d3,features=my.markers)
dev.off()

