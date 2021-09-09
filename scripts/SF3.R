#!/usr/bin/env Rscript
#.libPaths("/nfs/users/nfs_f/fr7/anaconda2/envs/r_env/lib/R/library")
m    <- modules::use("/Users/fr7/git_repos/single_cell/scripts/SC.R")
mono <- modules::use("/Users/fr7/git_repos/single_cell/scripts/monocle.R") 
library(argparse)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(gridExtra)
library(stats)
library(loomR)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(plyr)
library(patchwork)
library(clusterProfiler)
library(biomaRt)
library(org.Mm.eg.db)
library(inauguration)
library(monocle3)

figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'
files_dir<-'/Users/fr7/git_repos/single_cell/SupplementaryFiles'
merge_dir <- '/Users/fr7/git_repos/single_cell/experiment_4/FINAL/merge'
dir <- '/Users/fr7/git_repos/single_cell/experiment_4/FINAL/merge/undiff'
all.samples <- readRDS(file.path(merge_dir,'seurat_object.rds'))

### isolate the undiff cells and cluster them alone ###
Idents(all.samples)<-all.samples$cell.type
undiff <- subset(all.samples, idents = c("Und.1","Und.2","Und.3"))
undiff <- m$normalize_and_run_pca(undiff,40,dir)
undiff <- m$do_clustering(undiff,30,resolution=0.4,min.dist = 0.05)

#undiff.control <- subset(undiff,subset = condition == "Control")
#undiff.markers<-m$find_all_markers(undiff.control,undiff_dir)

new.undiff.ids <- list()
new.undiff.ids[["0"]] <- "S-Stem/TA cells"
new.undiff.ids[["4"]] <- "S-Stem/TA cells"
new.undiff.ids[["1"]] <- "G2M-Stem/TA cells"
new.undiff.ids[["3"]] <- "G2M-Stem/TA cells"
new.undiff.ids[["2"]] <- "Progenitor Entero.1"
new.undiff.ids[["7"]] <- "Progenitor Entero.1"
new.undiff.ids[["5"]] <- "Progenitor Entero.2"
new.undiff.ids[["6"]] <- "DSC"
undiff <- RenameIdents(undiff,new.undiff.ids)

# split UMAP of these cells

cols.undiff<-c(
  "G2M-Stem/TA cells"               =  "#a6cee3",
  "S-Stem/TA cells"                 =  "#fb9a99",
  "Progenitor Entero.1"             =  "#C7E9C0",
  "Progenitor Entero.2"             =  "#006D2C",
  "DSC"                             =  "#5c1a33"
)

cols <- c(
  "Und.1"                 =  "#a6cee3",
  "Und.2"                 =  "#b2df8a",
  "Und.3"                 =  "#fb9a99",
  "Early entero.1"        =  "#1f78b4",
  "Early entero.2"        =  "#ff7f00",
  "Enterocyte"            =  "#cab2d6",
  "Entero.Isg15"          =  "#e31a1c",
  "Entero.AMP"            =  "#6a3d9a",
  "Goblet"                =  "#33a02c",
  "Tuft"                  =  "#fdbf6f",
  "Ee"                    =  "#b15928"
)

cols.all <- c(cols.undiff,cols)

SF3A.split <- DimPlot(undiff,reduction="umap", cols=cols.undiff,pt.size=0.3, split.by = "condition")+
      theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 7, hjust = -0.01))

SFleg<-get_legend(SF3A.split)

pdf(file.path(figures_dir,"SF3A.legend.pdf"),2.5,3)
as_ggplot(leg)
dev.off()

SF3A.split <- SF3A.split + NoLegend()

pdf(file.path(figures_dir,"SF3A.split.pdf"),7.5,2.5)
print(SF3A.split)
dev.off()

phase.cols <- c(
  "S"   = "yellow",
  "G2M" = "magenta",
  "G1"  = "#3e236e" # dark purple
)
  
SF3B <- DimPlot(undiff,reduction="umap", group.by = "Phase",pt.size=0.3, split.by = "condition", cols = phase.cols)+
  theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 7, hjust = -0.01))

leg<-get_legend(SF3B)

pdf(file.path(figures_dir,"SF3B.legend.pdf"),2.5,3)
as_ggplot(leg)
dev.off()

SF3B <- SF3B + NoLegend()

pdf(file.path(figures_dir,"SF3B.pdf"),7.5,2.5)
print(SF3B)
dev.off()

# subset controls and calculate markers
undiff.control <- subset(undiff, subset = infection_status == "control")

# plot stem/TA markers
genes_to_plot <- c("Lgr5", "Ascl2", "Smoc2", "Lrig1", "Sox9", "Myc", "Hopx", "Cd44", "Birc5", "Pcna", "Mki67","Top2a","Ube2c","Hes1","Spdef", "Atoh1")
genes_to_plot_2<-c("F3","Hmgb2","Stmn1","Bm1","Text","Hope","Lrig1","Smoc2","Afap1l1","Agr3","Cnn3","Dach1","Slco3a1","Sorbs2","Tns3","Vdr","Aqp4","Cdca7","Cdk6","Clca4a", "Clca4b","Kcnq1","Nav1","Soat1")
SF3A.genes <- FeaturePlot(undiff.control,features = genes_to_plot, pt.size = 0.05, cols = c('grey',"#A50026"), combine = FALSE)
SF3A.genes.2 <- FeaturePlot(undiff.control,features = genes_to_plot_2, pt.size = 0.05, cols = c('grey',"#A50026"), combine = FALSE)
add_formatting<-function(gg){
  gg<-gg+theme(axis.text=element_text(size=12), axis.title=element_text(size=12), plot.title = element_text(size=12,face='bold.italic'))
  return(gg)
}
SF3A.genes.fig <- lapply(SF3A.genes,add_formatting)
SF3A.fig <- plot_grid(plotlist=SF3A.genes.fig, ncol = 4)

pdf(file.path(figures_dir,"SF3A.genes.2.pdf"),15,12)
print(SF3A.fig)
dev.off()



markers<-m$find_all_markers(undiff.control,dir)
markers<-markers[which(markers$p_val_adj < 0.05 ),]

# supplementary file
text<-"Cluster-defining marker genes for subclustered undifferentiated populations (control cells, Wilcox test)"
markers.file <- markers[,c("cluster","gene","pct.1","pct.2","avg_logFC","p_val","p_val_adj")]
colnames(markers.file) <- c("Cluster","Gene","Proportion of cells in cluster expressing marker", "Proportion of all other cells expressing marker", "Average log FoldChange", "p value", "Adjusted p value")
markers.file.text<-mrna$text_matrix(markers.file,text)
markers.file.text<-mrna$write_results_text(markers.file.text,files_dir,"SFile3.3.tsv")


markers<-markers[order(-abs(markers$avg_logFC)),]

# select top 5 markers from each cluster to plot 
# make sure we also plot known markers, even if not in top 5
#known_markers <- read.table('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/known_markers.txt',header=TRUE,stringsAsFactors = FALSE)
clusters<-levels(undiff.control@active.ident)
markers_to_plot<-c()
for (cluster in clusters){
  x<-markers[which(markers$cluster == cluster),'gene']
#  y<-known_markers[which(known_markers$cluster == cluster),'marker']
  markers_to_plot<-c(markers_to_plot, x[0:5])
}
markers_to_plot<-unique(markers_to_plot)
cols.2<-rev(brewer.pal(n=11,name="RdYlBu"))

SF3C <- DotPlot(undiff.control, features = rev(markers_to_plot), dot.scale = 4) +
  theme(text = element_text(size=7),
        axis.text=element_text(size=7),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(0,0,-10,-10) )+
  scale_color_gradientn(name = "Average expression\n(scaled)",   colours = cols.2)+
  guides(color = guide_colourbar(barwidth = 0.5, barheight = 2), 
         size = guide_legend(keywidth = 0.5, keyheight = 0.5))+
  labs(x=NULL,y=NULL)+
  scale_y_discrete(limits = rev(levels(undiff)))+    # reverse the order of the y axis
  scale_size_continuous(name = "Percentage of cells\nexpressing marker",range = c(0.1,1.5))


pdf(file.path(figures_dir,"SF3C.pdf"),5,2.2)
print(SF3C)
dev.off()


# add new undiff idents to the original Seurat object
other.cells <- subset(all.samples, idents = c("Und.1","Und.2","Und.3"), invert = TRUE)
all.samples <- AddMetaData(all.samples,
                           metadata = unlist(list(Idents(undiff), Idents(other.cells))), # no idea why this works, but this seems to be how to concatenate factor vectors.
                           col.name = 'cell.type.2')

Idents(all.samples)<-all.samples$cell.type.2
saveRDS(all.samples,file.path(dir,'seurat_object.rds'))

# move to monocle (control cells only)
all.samples.control <- subset(all.samples, subset = condition == 'Control')
cds <- mono$seurat_to_monocle(all.samples.control)

# normalize and pre-process
cds <- preprocess_cds(cds, num_dim = 100)

# remove batch effects 
cds <- align_cds(cds, alignment_group = "orig.ident")

# dimensional reduction
cds <- reduce_dimension(cds)

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster",label_cell_groups=FALSE)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "partition",label_cell_groups=FALSE)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_cluster",label_cell_groups=FALSE)

SF3D <- plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_cluster", reduction_method = "UMAP",label_cell_groups=FALSE, show_trajectory_graph = FALSE)
SF3D <- SF3D + scale_color_manual(values = cols.all)+
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5), legend.title = element_blank())

leg <- get_legend(SF3D)

pdf(file.path(figures_dir,"SF3D.legend.pdf"),2.5,4)
as_ggplot(leg)
dev.off()

SF3D <- SF3D + NoLegend()
pdf(file.path(figures_dir,"SF3D.pdf"),2.5,2.5)
print(SF3D)
dev.off()

SF3E <- plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "phase", reduction_method = "UMAP",label_cell_groups=FALSE, show_trajectory_graph = FALSE)
SF3E <- SF3E +
  scale_color_manual(values = phase.cols)+
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5))+
  NoLegend()

pdf(file.path(figures_dir,'SF3E.pdf'),2.5,2.5)
print(SF3E)
dev.off()

cds <- order_cells(cds)
SF3F <- plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "pseudotime", reduction_method = "UMAP",label_cell_groups=FALSE, label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE)
SF3F <- SF3F +
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5))

  
pdf(file.path(figures_dir,'SF3F.pdf'),3,2.5)
print(SF3F)
dev.off()

# export for PAGA
all.samples.control <- FindVariableFeatures(all.samples.control)
control.loom <- as.loom(all.samples.control, filename = file.path(merge_dir,"control.loom"), verbose = TRUE)


# plots for rebuttal

# B1 probes
pdf(file.path(figures_dir,"B1_probes.pdf"),10,5)
VlnPlot(all.samples,features=c("Lgr5","Ascl2","Isg15"),pt.size=0)+theme(axis.text=element_text(size=5), axis.title.x=element_blank(), axis.title.y = element_text(size=5))
dev.off()

# B2 probes
pdf(file.path(figures_dir,"B3_probes.pdf"),10,5)
VlnPlot(all.samples,features=c("Reg4","Muc2","Krt20"),pt.size=0)+theme(axis.text=element_text(size=5), axis.title.x=element_blank(), axis.title.y = element_text(size=5))
dev.off()

# B4 probes
pdf(file.path(figures_dir,"B4_probes.pdf"),13,5)
VlnPlot(all.samples,features=c("Stmn1","Mki67","Ube2c","Top2a"),pt.size=0, ncol = 4)+theme(axis.text=element_text(size=5), axis.title.x=element_blank(), axis.title.y = element_text(size=5))
dev.off()

# B5 probes
pdf(file.path(figures_dir,"B5_probes.pdf"),10,5)
VlnPlot(all.samples,features=c("Hmgb2","Birc5","Car1"),pt.size=0)+theme(axis.text=element_text(size=5), axis.title.x=element_blank(), axis.title.y = element_text(size=5))
dev.off()

# Try to recluster the S-Stem/TA cells
test <- subset(undiff, idents = c("S-Stem/TA cells"))
test <- m$normalize_and_run_pca(test,40,dir)
test <- m$do_clustering(test,30,resolution=0.4,min.dist = 0.05)

DimPlot(test,reduction="umap",pt.size=0.3, split.by = "condition")+
theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 7, hjust = -0.01))

markers.test<-m$find_all_markers(test,dir)
markers.test<-markers.test[which(markers.test$p_val_adj < 0.05 ),]
