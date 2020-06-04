#!/usr/bin/env Rscript
#.libPaths("/nfs/users/nfs_f/fr7/anaconda2/envs/r_env/lib/R/library")
m <- modules::use("/Users/fr7/git_repos/single_cell/scripts/SC.R")
library(argparse)
library(Seurat)
library(ggplot2)
library(DescTools)

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

control<-readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/seurat_object.rds')
d1<-readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/day1_classify/seurat_object.rds')
d3<-readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/day3_classify/seurat_object.rds')

dir<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated'
metadata<-read.table('/Users/fr7/git_repos/single_cell/metadata/samples.txt', header=T)

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

entero<-proportions[which(proportions$celltype %like% "%ENTEROCYTE%"),]
entero.cols<-cols[grep("ENTEROCYTE",names(cols))]

g<-ggplot(entero,aes(condition.i,proportion,fill=celltype)) +
  geom_boxplot()+
  scale_fill_manual(values=entero.cols)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) + 
  facet_wrap(~celltype, scales="free_y", ncol=2)+
  NoLegend()+
  ylab("Proportion of cells")+
  xlab("")

pdf(file.path(dir,"enterocytes.pdf"),6,6)
print(g)
dev.off()


nonentero<-proportions[which(! proportions$celltype %like% "%ENTEROCYTE%"),]
nonentero.cols<-cols[grep("ENTEROCYTE",names(cols),invert=TRUE)]



g<-ggplot(nonentero,aes(condition.i,proportion,fill=celltype)) +
  geom_boxplot()+
  scale_fill_manual(values=nonentero.cols)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) + 
  facet_wrap(~celltype, scales="free_y", ncol=2)+
  NoLegend()+
  ylab("Proportion of cells")+
  xlab("")

pdf(file.path(dir,"enterocytes.pdf"),6,6)
print(g)
dev.off()