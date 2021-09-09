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
library(randomForest)
library(loomR)
library(ggpubr)


cols.2 <- c(
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


dir <- '/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/day3_classify'
control <- readRDS('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/seurat_object.rds')
samples_directory<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/QC/combined'
meta_data<-read.table('/Users/fr7/git_repos/single_cell/metadata/samples.txt',header=T)
condition<-'infected.72'
figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'

#find the IDs of the samples that we want
samples<-meta_data[which(meta_data$condition == condition & meta_data$experiment == '4'),'sample_id']

#load seurat objects
objects<-list()
for (sample in samples){
  object<-readRDS(file.path(samples_dir,paste0(sample,".rds")))
  objects[[sample]]<-object
}

#merge all replicates into one object
day.3 <- merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects), project = "day3" )

#normalize, scale, run PCA
day.3<- m$normalize_and_run_pca(day.3,40,dir)
day.3<-m$do_clustering(day.3,35,resolution=0.5,min.dist=0.05)

#plots
markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Isg15","Ace2","Ptprc")
p1<-m$QC_violins(day.3, dir)
p2<-m$QC_scatters(day.3,dir)
p3<-m$plot_UMAPS(day.3,dir)
p4<-m$marker_violins(day.3,dir,markers_to_plot)
p5<-m$markers_feature_plot(day.3,dir,markers_to_plot)

#markers
DefaultAssay(day.1)<-"RNA"
markers<-m$find_all_markers(day.3,dir)
marker_dots<-m$marker_dotplot(day.3,markers,dir)

#transfer labels
DefaultAssay(control)<-"integrated"
control$celltype<-control@active.ident
anchors <- FindTransferAnchors(reference =  control, query = day.3, dims = 1:30 )
predictions <- TransferData(anchorset = anchors, refdata = control$celltype, dims = 1:30)
day.3<-AddMetaData(day.3, metadata = predictions)
Idents(day.3)<-"predicted.id"

#
#Idents(day.3)<-"seurat_clusters"
DimPlot(day.3,reduction="umap", cols=cols)

pdf(file.path(dir,"UMAP_clusters.pdf"))
DimPlot(day.3,reduction="umap",label=TRUE,label.size=2,repel=TRUE, cols=cols)+NoLegend()
dev.off()

pdf(file.path(dir,"UMAP_batches.pdf"))
DimPlot(day.3,reduction="umap",group.by = "orig.ident")
dev.off()

saveRDS(day.3,file.path(dir,'seurat_object.rds'))

day.3<-readRDS(file.path(dir,'seurat_object.rds'))

p<-DimPlot(day.3,reduction="umap",label=TRUE,label.size=2,repel=TRUE, cols=cols)+NoLegend()
select.cells <- CellSelector(plot = p)
Idents(day.3, cells = select.cells) <- "NewCells"

#refind all markers
markers<-m$find_all_markers(day.3,dir)
marker_dots<-m$marker_dotplot(day.3,markers,dir)

day.3<-FindVariableFeatures(day.3)
day.3@graphs<-list()
day.3.loom <- as.loom(day.3, filename = file.path(dir,"day.3.loom"), verbose = TRUE)

fig.4a.d3 <-DimPlot(d3,reduction="umap",label=TRUE,label.size=1,repel=TRUE, cols=cols,pt.size=0.1)+
  NoLegend() +
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5))

pdf(file.path(figures_dir,"fig.4a.d3.pdf"),2.5,2.5)
fig.4a.d3
dev.off()

F20 <-DimPlot(day.3,reduction="umap", cols=cols.2,pt.size=0.3)+
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5))

leg<-get_legend(F20)

pdf(file.path(figures_dir,"umap.legend.pdf"),2.5,2.8)
as_ggplot(leg)
dev.off()

F20 <-DimPlot(day.3,reduction="umap", cols=cols.2,pt.size=0.1)+
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5))+
  NoLegend()

pdf(file.path(figures_dir,"F20.pdf"),2.5,2.5)
F20
dev.off()

###VlnPlots for NewCells###
new.genes<-c("Muc2", "Reg4","Ido1","Fcgbp","Clca1","Fer1l6","Zg16","Sct","Isg15","Sytl2")
levels(d3)<-c("NewCells",levels(control))
figs<-list()
for (gene in new.genes){
  fig <- VlnPlot(d3,features=gene,pt.size=0,cols=cols)+
    theme(text = element_text(size=2),
          axis.text.y = element_text(size=4,margin=margin(-100,0,0,0),hjust=-100),
          axis.text.x=element_blank(),
          legend.position = "none",
          axis.title.x =element_blank(),
          axis.title.y =element_blank(),
          axis.line = element_line(size = 0.1),
          axis.ticks=element_blank(),
          plot.title=element_text(size=5, face="plain")
    )
  
  fig$layers[[1]]$aes_params$size = 0.1
  figs[[gene]]<-fig
}

fig.4d.d3<-plot_grid(plotlist=figs, ncol=10)

pdf(file.path(figures_dir,"fig.4d.d3.pdf"),7,1)
fig.4d.d3
dev.off()

#get legend


#for the legend
p<- VlnPlot(d1,features=gene,pt.size=0,cols=cols)+
  theme(text = element_text(size=5),
        axis.text.y = element_text(size=4,margin=margin(-100,0,0,0),hjust=-100),
        axis.text.x=element_blank(),
        axis.title.x =element_blank(),
        axis.title.y =element_blank(),
        axis.line = element_line(size = 0.1),
        axis.ticks=element_blank(),
        plot.title=element_text(size=5, face="plain"))



leg <- get_legend(p)

pdf(file.path(figures_dir,"fig.4d.legend.pdf"),1,5)
as_ggplot(leg)
dev.off()
