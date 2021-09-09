#!/usr/bin/env Rscript
#.libPaths("/nfs/users/nfs_f/fr7/anaconda2/envs/r_env/lib/R/library")
m <- modules::use("/Users/fr7/git_repos/single_cell/scripts/SC.R")
library(argparse)
library(Seurat)
library(ggplot2)
library(DescTools)
library(plyr)


cols<- c(
  "Und.1" = "#a6cee3",
  "Und.2" = "#fb9a99" ,  
  "Und.3" =  "#b2df8a",
  "Entero.1" = "#1f78b4",
  "Entero.2" = "#ff7f00",
  "Entero.3" = "#cab2d6",
  "Entero.Isg15" = "#e31a1c",
  "Goblet" = "#33a02c",
  "Tuft" = "#fdbf6f",
  "Ee" = "#b15928"
)

control<-readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/seurat_object.rds')
d1<-readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/day1_classify/seurat_object.rds')
d3<-readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/day3_classify/seurat_object.rds')

dir<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated'
metadata<-read.table('/Users/fr7/git_repos/single_cell/metadata/samples.txt', header=T)
figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'
clusters<-levels(control@active.ident)

numbers.control<-as.data.frame(table(Idents(control), control$orig.ident))
numbers.d1<-as.data.frame(table(Idents(d1), d1$orig.ident))
numbers.d3<-as.data.frame(table(Idents(d3), d3$orig.ident))

proportions.control<-as.data.frame(prop.table(table(Idents(control), control$orig.ident),margin=2))
proportions.d1<-as.data.frame(prop.table(table(Idents(d1), d1$orig.ident),margin=2))
proportions.d3<-as.data.frame(prop.table(table(Idents(d3), d3$orig.ident),margin=2))

proportions<-rbind(proportions.control,proportions.d1,proportions.d3)
colnames(proportions) <- c("celltype","sample","proportion")

proportions<-merge(proportions,metadata,by.x = "sample",by.y="sample_id")

#change to nice label names
proportions$condition.i <- revalue(proportions$condition.i,c("control" = "C", "d1" = "D1", "d3" = "D3"))

g<-ggplot(proportions,aes(condition.i,proportion,fill=celltype)) +
#  geom_dotplot(binaxis = "y", stackdir = "center")+ 
  geom_boxplot()+
  scale_fill_manual(values=cols)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_wrap(~celltype, scales="free_y", ncol=2)+
  NoLegend()+
  ylab("Proportion of cells")+
  xlab("")

pdf(file.path(dir,"cluster.sizes.pdf"),5,10)
print(g)
dev.off()

entero<-proportions[which(proportions$celltype %like% "%Entero%"),]
entero.cols<-cols[grep("Entero",names(cols))]

g<-ggplot(entero,aes(condition.i,proportion,fill=celltype)) +
  geom_boxplot()+
  scale_fill_manual(values=entero.cols)+
  facet_wrap(~celltype, scales="free_y", ncol=2)+
  theme(axis.text.y = element_text(size=5), strip.text.x = element_text(size=20)) + 
  NoLegend()+
  ylab("Proportion of cells")+
  xlab("")

pdf(file.path(dir,"enterocytes.pdf"),6,6)
print(g)
dev.off()

pdf(file.path(figures_dir,"fig.4c.pdf"),2.5,3)
print(g)
dev.off()

nonentero<-proportions[which(! proportions$celltype %like% "%Entero%"),]
nonentero.cols<-cols[grep("Entero",names(cols),invert=TRUE)]



g<-ggplot(nonentero,aes(condition.i,proportion,fill=celltype)) +
  geom_boxplot()+
  scale_fill_manual(values=nonentero.cols)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) + 
  facet_wrap(~celltype, scales="free_y", ncol=2)+
  NoLegend()+
  ylab("Proportion of cells")+
  xlab("")

pdf(file.path(dir,"non-enterocytes.pdf"),6,6)
print(g)
dev.off()