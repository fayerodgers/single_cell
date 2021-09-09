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
library(loomR)
library(RColorBrewer)
library("scales")
library(ggpubr)

parser <- ArgumentParser()
parser$add_argument("--results_directory",help="directory for files to be written to")
parser$add_argument("--samples_directory",help="directory to find saved, QC'd seurat objects")
parser$add_argument("--metadata", help="File with sample metadata")
parser$add_argument("--condition", help="Integrate and cluster samples from this condition")

args <- parser$parse_args()

#for now
args$results_directory<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control'
args$samples_directory<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/QC/combined'
args$metadata<-'/Users/fr7/git_repos/single_cell/metadata/samples.txt'
args$condition<-'control'
figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'

dir <- args$results_directory
samples_dir <- args$samples_directory
meta_data<-read.table(args$metadata,header=T)

#find the IDs of the samples that we want
samples<-meta_data[which(meta_data$treatment == 'control' & meta_data$experiment == '4'),'sample_id']

#load seurat objects
objects<-list()
for (sample in samples){
  object<-readRDS(file.path(samples_dir,paste0(sample,".rds")))
  objects[[sample]]<-object
}

#normalize and find variable features
for (i in names(objects)){
  objects[[i]] <- NormalizeData(objects[[i]])
  objects[[i]] <- FindVariableFeatures(objects[[i]], selection.method = "vst", nfeatures = 2000)
}

#integrate and cluster
combined<-m$integrate_datasets(objects,40)
p<-m$elbow_plot(combined,40,"pca", dir)
combined<-m$do_clustering(combined,35,resolution=0.7,min.dist=0.05)
saveRDS(combined,file.path(dir,'seurat_object.rds'))

#plots
markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Isg15","Reg3g")
p1<-m$QC_violins(combined, dir)
p2<-m$QC_scatters(combined,dir)
p3<-m$plot_UMAPS(combined,dir)
p4<-m$marker_violins(combined,dir,markers_to_plot)
p5<-m$markers_feature_plot(combined,dir,markers_to_plot)

#select the Chga expressing cells as their own cluster
#p<-FeaturePlot(combined,features="Chga")
#combined <- CellSelector(plot = p, object = combined, ident = "ENTEROENDOCRINE")

#rename/merge clusters

new.clusters.ids <-list()
new.clusters.ids[["3"]] <-  "Und.1"
new.clusters.ids[["7"]] <-  "Und.2"
new.clusters.ids[["4"]] <-  "Und.3"
new.clusters.ids[["8"]] <-  "Entero.1"
new.clusters.ids[["0"]] <-  "Entero.1"
new.clusters.ids[["5"]] <-  "Entero.1"
new.clusters.ids[["10"]] <- "Entero.1"
new.clusters.ids[["1"]] <-  "Entero.2"
new.clusters.ids[["2"]] <-  "Entero.2"
new.clusters.ids[["9"]] <-  "Entero.2"
new.clusters.ids[["11"]] <- "Entero.2"
new.clusters.ids[["14"]] <- "Entero.2"
new.clusters.ids[["15"]] <- "Entero.3"
new.clusters.ids[["13"]] <- "Entero.Isg15"
new.clusters.ids[["6"]] <-  "Goblet"
new.clusters.ids[["12"]] <- "Goblet"
new.clusters.ids[["16"]] <- "Tuft"
new.clusters.ids[["17"]] <- "Ee"
combined <- RenameIdents(combined,new.clusters.ids)

saveRDS(combined,file.path(dir,'seurat_object.rds'))

DefaultAssay(combined)<-"RNA"
markers<-m$find_all_markers(combined,dir)
#maybe
#markers<-read.table(file.path(dir,"markers.txt"),stringsAsFactors = FALSE)
markers<-markers[which(markers$p_val_adj < 0.05 ),]
markers<-markers[order(-abs(markers$avg_logFC)),]

marker_dots<-m$marker_dotplot(combined,markers,figures_dir)

F10<-marker_dots +
        theme(text = element_text(size=5),
         axis.text=element_text(size=5),
         legend.title = element_text(size = 5), 
         legend.text = element_text(size = 5),
         legend.margin=margin(0,0,0,0), 
         legend.box.margin=margin(-10,-10,-10,-10) )+
        scale_color_gradient(name = "Average expression\n(scaled)",   low = "#56B1F7",
                             high = "#132B43")+
        guides(color = guide_colourbar(barwidth = 0.5, barheight = 2), 
         size = guide_legend(keywidth = 0.5, keyheight = 0.5))+
         labs(x=NULL,y=NULL)+
        scale_size_continuous(name = "Percentage of cells\nexpressing marker",range = c(0.1,1.5))
  

pdf(file.path(figures_dir,"F10.pdf"),5,2.5)
F10
dev.off()



cols.2 <- c(
  "Und.1" = "#a6cee3",
  "Und.2" = "#fb9a99" ,  
  "Und.3" =  "#b2df8a",
  "Entero.1" = "#1f78b4",
  "Entero.2" = "#ff7f00",
  "Entero.3" = "#cab2d6",
#  "Entero.4" = "#ffff99",
#  "Entero.5" = "#6a3d9a",
  "Entero.Isg15" = "#e31a1c",
  "Goblet" = "#33a02c",
  "Tuft" = "#fdbf6f",
  "Ee" = "#b15928"
)

#UMAP for figure
pdf(file.path(dir,"UMAP_clusters.pdf"))
DimPlot(combined,reduction="umap",label=TRUE,label.size=2,repel=TRUE, cols=cols.2)+NoLegend()
dev.off()

F6 <-DimPlot(combined,reduction="umap",label=TRUE,label.size=1,repel=TRUE, cols=cols.2,pt.size=0.1)+
         NoLegend() +
         theme(axis.text=element_text(size=5), axis.title=element_text(size=5))

pdf(file.path(figures_dir,"F6.pdf"),2.5,2.5)
F6
dev.off()

umap.for.grant <-DimPlot(combined,reduction="umap", cols=cols.2,pt.size=0.3)+
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5))

leg<-get_legend(umap.for.grant)

pdf(file.path(figures_dir,"F6.legend.pdf"),2.5,2.7)
as_ggplot(leg)
dev.off()

umap.for.grant <-DimPlot(combined,reduction="umap", cols=cols.2,pt.size=0.1)+
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5))+
  NoLegend()

pdf(file.path(figures_dir,"F6.pdf"),2.5,2.5)
umap.for.grant
dev.off()

pdf(file.path(figures_dir,"UMAP_batches.pdf"))
DimPlot(combined,reduction="umap",group.by = "orig.ident")
dev.off()

DefaultAssay(combined)<-"RNA"
markers<-m$find_all_markers(combined,figures_dir)
marker_dots<-m$marker_dotplot(combined,markers,figures_dir)

saveRDS(combined,file.path(dir,'seurat_object.rds'))

#also save the object in loom format for scanpy
#need to rerun findvariable features to export as loom
combined<-FindVariableFeatures(combined)
combined.loom <- as.loom(combined, filename = file.path(dir,"control.loom"), verbose = TRUE)

#extract cycling cells
cycling<-subset(combined,idents=c("UNDIFFERENTIATED.1","UNDIFFERENTIATED.2","UNDIFFERENTIATED.3"))
cycling<-FindVariableFeatures(cycling)
cycling.loom<-as.loom(cycling, filename = file.path(dir,"cycling.loom"), verbose = TRUE)

DefaultAssay(cycling)<-"RNA"
cycling_dir<-file.path(dir,"cycling")
dir.create(cycling_dir)

#subcluster cycling cells
cycling <- m$normalize_and_run_pca(cycling,50,cycling_dir)
cycling <- m$do_clustering(cycling,35,resolution=0.2, min.dist = 0.5)

p1<-m$QC_violins(cycling, cycling_dir)
p2<-m$QC_scatters(cycling,cycling_dir)
p3<-m$plot_UMAPS(cycling,cycling_dir)
p4<-m$marker_violins(cycling,cycling_dir,markers_to_plot)
p5<-m$markers_feature_plot(cycling,cycling_dir,markers_to_plot)

cycling_markers<-m$find_all_markers(cycling,cycling_dir)
marker_dots<-m$marker_dotplot(cycling,cycling_markers,cycling_dir)
cycling@graphs<-list()
cycling.loom<-as.loom(cycling, filename = file.path(cycling_dir,"cycling.loom"), verbose = TRUE)
saveRDS(cycling,file.path(cycling_dir,"seurat_object.rds"))

#extract enterocytes
enterocytes<-subset(combined,idents=c("ENTEROCYTE.1","ENTEROCYTE.2","ENTEROCYTE.3","ENTEROCYTE.4","ENTEROCYTE.ISG15","ENTEROCYTE.ACE2"))
enterocytes<-FindVariableFeatures(enterocytes)
enterocytes.loom<-as.loom(enterocytes, filename = file.path(dir,"enterocytes.loom"), verbose = TRUE)

#GO term plots
go<-read.csv(file.path(dir,"curated_go_terms.csv"))

#enterocyte.isg15 go term plot
go.isg15<-go[which(go$cluster == "ENTEROCYTE.ISG15"),]

fig.2c<-ggplot(go.isg15,x=cluster,y=PathwayName)+
  geom_point(aes(x=cluster,y=PathwayName, size=ngenes, col=PathwayPadj))+
  scale_color_gradient(limits = c(0.0001,0.01),
                       oob=squish, 
                       name = "Pathway p-value\n(adjusted)", 
                       breaks = c(0.0001,0.01),
                       labels = c("< 0.0001","> 0.01"),
                       guide = guide_colourbar())+
  scale_size_continuous(name="Number of\ngenes", range = c(0.1,2))+
  theme(text = element_text(size=5),
        axis.text.x=element_text(size=5, angle=45,hjust=1),
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-10,-10) )+
  labs(x=NULL,y=NULL)+
  guides(color = guide_colourbar(barwidth = 0.5, barheight = 2), 
         size = guide_legend(keywidth = 0.5, keyheight = 0.5))
             
pdf(file.path(figures_dir,'fig.2c.pdf'),2.5,1.8)
print(fig.2c)
dev.off()

#isg.15.genes<-scan(file = file.path(dir,'isg15.genes.txt'), sep = "\n", what = character())
#isg.15.genes.selected<-isg.15.genes[isg.15.genes != "B2m" & isg.15.genes != "H2-D1" & isg.15.genes != "Ifih1" & isg.15.genes != "H2-T23"]

isg.15.genes.selected<-c("Ddx60","Ifit1","Isg15","Oasl2","Stat1","Zbp1")
  
figs<-list()
for (gene in isg.15.genes.selected){
  fig <- VlnPlot(combined,features=gene,pt.size=0,cols=cols.2)+
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

F11<-plot_grid(plotlist=figs, ncol=3)
                 
pdf(file.path(figures_dir,"F11.pdf"),4,2.5)
F11
dev.off()

#for the legend
p<- VlnPlot(combined,features=gene,pt.size=0,cols=cols)+
  theme(text = element_text(size=5),
        axis.text.y = element_text(size=4,margin=margin(-100,0,0,0),hjust=-100),
        axis.text.x=element_blank(),
        axis.title.x =element_blank(),
        axis.title.y =element_blank(),
        axis.line = element_line(size = 0.1),
        axis.ticks=element_blank(),
        plot.title=element_text(size=5, face="plain"))
        

  
leg <- get_legend(p)

pdf(file.path(figures_dir,"fig.2d.legend.pdf"),1,5)
as_ggplot(leg)
dev.off()


#enterocyte.ace2 go term plot
go.ace2<-go[which(go$cluster == "ENTEROCYTE.ACE2"),]

fig.2e<-ggplot(go.ace2,x=cluster,y=PathwayName)+
  geom_point(aes(x=cluster,y=PathwayName, size=ngenes, col=PathwayPadj))+
  scale_color_gradient(limits = c(0.0001,0.01),
                       oob=squish, 
                       name = "Pathway p-value\n(adjusted)", 
                       breaks = c(0.0001,0.01),
                       labels = c("< 0.0001","> 0.01"),
                       guide = guide_colourbar())+
  scale_size_continuous(name="Number of\ngenes", range = c(0.1,2))+
  theme(text = element_text(size=5),
        axis.text.x=element_text(size=5, angle=45,hjust=1),
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-10,-10) )+
  labs(x=NULL,y=NULL)+
  guides(color = guide_colourbar(barwidth = 0.5, barheight = 2), 
         size = guide_legend(keywidth = 0.5, keyheight = 0.5))

pdf(file.path(figures_dir,'fig.2e.pdf'),1.8,1.8)
print(fig.2e)
dev.off()


ace2.genes<-c("Reg3g","Reg3b","Ace2","Fabp6")
figs<-list()
for (gene in ace2.genes){
  fig <- VlnPlot(combined,features=gene,pt.size=0,cols=cols.2)+
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

F12<-plot_grid(plotlist=figs, ncol=2)

pdf(file.path(figures_dir,"F12.pdf"),4,2.5)
F12
dev.off()

###Reg3g/Reg3b/Ace2 feature plot###



###VlnPlots for NewCells###
new.genes<-c("Muc2", "Reg4","Ido1","Fcgbp","Clca1","Fer1l6","Zg16","Sct","Isg15","Sytl2")
figs<-list()
for (gene in new.genes){
  fig <- VlnPlot(combined,features=gene,pt.size=0,cols=cols)+
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

fig.4d.control<-plot_grid(plotlist=figs, ncol=10)

pdf(file.path(figures_dir,"fig.4d.control.pdf"),7,1)
fig.4d.control
dev.off()

# finding GO terms for all clusters
markers<-read.table(file.path(dir,"markers.txt"),stringsAsFactors = FALSE)


