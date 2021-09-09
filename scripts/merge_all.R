#!/usr/bin/env Rscript
#.libPaths("/nfs/users/nfs_f/fr7/anaconda2/envs/r_env/lib/R/library")
m <- modules::use("/Users/fr7/git_repos/single_cell/scripts/SC.R")
mrna <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")
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

parser <- ArgumentParser()
parser$add_argument("--results_directory",help="directory for files to be written to")
parser$add_argument("--samples_directory",help="directory to find saved, QC'd seurat objects")
parser$add_argument("--metadata", help="File with sample metadata")
parser$add_argument("--condition", help="Integrate and cluster samples from this condition")

args <- parser$parse_args()

#for now
args$results_directory<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/merge'
args$samples_directory<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/QC/combined'
args$metadata<-'/Users/fr7/git_repos/single_cell/metadata/samples.txt'
#args$condition<-'control'
figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'
files_dir<-'/Users/fr7/git_repos/single_cell/SupplementaryFiles'

dir <- args$results_directory
samples_dir <- args$samples_directory
meta_data<-read.table(args$metadata,header=T)

#find the IDs of the samples that we want
samples<-meta_data[which(meta_data$experiment == '4'),'sample_id']

#load seurat objects
objects<-list()
for (sample in samples){
  object<-readRDS(file.path(samples_dir,paste0(sample,".rds")))
  objects[[sample]]<-object
}

all.samples <- merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects), project = "all" )

all.samples <- m$normalize_and_run_pca(all.samples,40,dir)
all.samples <- m$do_clustering(all.samples,30,resolution=0.6,min.dist = 0.05)

markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Ube2c","Dclk1","Cdk4","Isg15","Reg3g")
p1<-m$QC_violins(all.samples, dir)
p2<-m$QC_scatters(all.samples,dir)
p3<-m$plot_UMAPS(all.samples,dir)
p4<-m$marker_violins(all.samples,dir,markers_to_plot)
p5<-m$markers_feature_plot(all.samples,dir,markers_to_plot)

# DimPlot(all.samples,cells.highlight = WhichCells(all.samples,idents ="0"))
# Initial cell types
new.clusters.ids <-list()
new.clusters.ids[["3"]] <-  "Und.1"
new.clusters.ids[["11"]] <- "Und.2"
new.clusters.ids[["4"]] <-  "Und.3"
new.clusters.ids[["8"]] <-  "Early entero.1"
new.clusters.ids[["1"]] <-  "Early entero.2"
new.clusters.ids[["5"]] <-  "Early entero.2"
new.clusters.ids[["0"]] <-  "Enterocyte"
new.clusters.ids[["2"]] <-  "Enterocyte"
new.clusters.ids[["9"]] <-  "Enterocyte"
new.clusters.ids[["12"]] <- "Enterocyte"
new.clusters.ids[["7"]] <-  "Entero.Isg15"
new.clusters.ids[["15"]] <- "Entero.AMP"
new.clusters.ids[["6"]] <-  "Goblet"
new.clusters.ids[["10"]] <- "Goblet"
new.clusters.ids[["13"]] <- "Ee"
new.clusters.ids[["14"]] <- "Tuft"
all.samples <- RenameIdents(all.samples,new.clusters.ids)
# stash idents as orig.cell.type
all.samples$cell.type <- all.samples@active.ident

saveRDS(all.samples,file.path(dir,'seurat_object.rds'))
# all.samples <- readRDS(file.path(dir,'seurat_object.rds'))

# extract only control samples to find cluster markers
control.samples<-subset(all.samples, subset = infection_status == "control")
markers<-m$find_all_markers(control.samples,dir)


### Fig1B - UMAP ###

cols <- c(
  "Und.1"                 =  "#a6cee3",
  "Und.2"                 =  "#b2df8a",
  "Und.3"                 =  "#fb9a99",
  "Early entero.1"        =  "#1f78b4",
  "Early entero.2"        =  "#ff7f00",
  "Enterocyte"            =  "#cab2d6",
  "Entero.Isg15"          =  "#e31a1c",
  "Entero.AMP"             =  "#6a3d9a",
  "Goblet"                =  "#33a02c",
  "Tuft"                  =  "#fdbf6f",
  "Ee"                    =  "#b15928"
)

F1B <-DimPlot(all.samples,reduction="umap", cols=cols,pt.size=0.3)+
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5))

leg<-get_legend(F1B)

pdf(file.path(figures_dir,"F1B.legend.pdf"),2.5,3)
as_ggplot(leg)
dev.off()

F1B <-DimPlot(all.samples,reduction="umap", cols=cols,pt.size=0.1)+
  theme(axis.text=element_text(size=5), axis.title=element_text(size=5))+
  NoLegend()

pdf(file.path(figures_dir,"F1B.pdf"),2.5,2.5)
print(F1B)
dev.off()

# add new metadata
Idents(all.samples)<- all.samples$time.status
new.condition.ids <-list()
new.condition.ids[["24hpi.control"]]  <- "Control"
new.condition.ids[["72hpi.control"]]  <- "Control"
new.condition.ids[["24hpi.infected"]] <- "D1 infected"
new.condition.ids[["72hpi.infected"]] <- "D3 infected"
all.samples <- RenameIdents(all.samples,new.condition.ids)
all.samples$condition <- all.samples@active.ident
Idents(all.samples) <- all.samples$cell.type

F1B.split <-DimPlot(all.samples,reduction="umap", cols=cols,pt.size=0.1,split.by="condition")+
  theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 11, hjust = -0.01))+
  NoLegend()

pdf(file.path(figures_dir,"F1B.split.pdf"),7.5,2.5)
print(F1B.split)
dev.off()

### SF2E - split UMAP ###
Idents(all.samples) <- all.samples$infection_status
infection_status.caps<-list()
infection_status.caps[["control"]] <- "Control"
infection_status.caps[["infected"]] <- "Infected"
all.samples <- RenameIdents(all.samples,infection_status.caps)
all.samples$infection_status_caps<- all.samples@active.ident
Idents(all.samples) <- all.samples$cell.type

d1.samples <- subset(all.samples,subset = time == "24hpi")
SF2E.d1 <-  DimPlot(d1.samples,reduction="umap", cols=cols,pt.size=0.1,split.by="infection_status_caps")+
  theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 11, hjust = -0.01))+
  NoLegend()

pdf(file.path(figures_dir,"SF2E.d1.pdf"),5,2.5)
print(SF2E.d1)
dev.off()

d3.samples <- subset(all.samples,subset = time == "72hpi")
SF2E.d3 <-  DimPlot(d3.samples,reduction="umap", cols=cols,pt.size=0.1,split.by="infection_status_caps")+
  theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 11, hjust = -0.01))+
  NoLegend()

pdf(file.path(figures_dir,"SF2E.d3.pdf"),5,2.5)
print(SF2E.d3)
dev.off()


### Fig SF2B/C ###

Idents(all.samples) <- all.samples$orig.ident
new.mouse.ids<-list()
new.mouse.ids[["4672STDY8112878"]] <- "C1 (D1)"
new.mouse.ids[["4672STDY8112879"]] <- "C2 (D1)"
new.mouse.ids[["4672STDY8112880"]] <- "C3 (D1)"
new.mouse.ids[["4672STDY8112881"]] <- "C4 (D1)"
new.mouse.ids[["4672STDY8112882"]] <- "I1 (D1)"
new.mouse.ids[["4672STDY8112883"]] <- "I2 (D1)"
new.mouse.ids[["4672STDY8112884"]] <- "I3 (D1)"
new.mouse.ids[["4672STDY8112885"]] <- "I4 (D1)"
new.mouse.ids[["4672STDY8112974"]] <- "C5 (D3)"
new.mouse.ids[["4672STDY8112975"]] <- "C6 (D3)"
new.mouse.ids[["4672STDY8112976"]] <- "C7 (D3)"
new.mouse.ids[["4672STDY8112977"]] <- "C8 (D3)"
new.mouse.ids[["4672STDY8113070"]] <- "I5 (D3)"
new.mouse.ids[["4672STDY8112979"]] <- "I6 (D3)"
new.mouse.ids[["4672STDY8112980"]] <- "I7 (D3)"
new.mouse.ids[["4672STDY8112981"]] <- "I8 (D3)"
all.samples <- RenameIdents(all.samples,new.mouse.ids)
all.samples$mouse.id <- all.samples@active.ident
Idents(all.samples) <- all.samples$cell.type

batch.cols.2 <- c(
  "C1 (D1)" = "tomato1",
  "C2 (D1)" = "tomato2",
  "C3 (D1)" = "tomato3",
  "C4 (D1)" = "tomato4",
  "I1 (D1)" = "palegreen1",
  "I2 (D1)" = "palegreen2",
  "I3 (D1)" = "palegreen3",
  "I4 (D1)" = "palegreen4",
  "C5 (D3)" = "steelblue1",
  "C6 (D3)" = "steelblue2",
  "C7 (D3)" = "steelblue3",
  "C8 (D3)" = "steelblue4",
  "I5 (D3)" = "tan1",
  "I6 (D3)" = "tan2",
  "I7 (D3)" = "tan3",
  "I8 (D3)" = "tan4"
)

SF2B <-VlnPlot(all.samples,features = "percent.mito", group.by = "mouse.id", pt.size = 0, cols = batch.cols.2)+
      ylab("Percentage mitochondrial genes")+
      theme(axis.text=element_text(size=5), axis.title.x=element_blank(), axis.title.y = element_text(size=5))+
      ggtitle(NULL)+
      NoLegend()

pdf(file.path(figures_dir,"SF2B.pdf"),5,2)
print(SF2B)
dev.off()


SF2D <-DimPlot(all.samples, group.by = "mouse.id", reduction="umap",pt.size=0.05,cols = batch.cols.2,  ncol = 2)+
  theme(axis.text=element_text(size=5), axis.title=element_text(size=2))

leg<-get_legend(SF2D)

pdf(file.path(figures_dir,"SF2D.legend.pdf"),2.5,3.5)
as_ggplot(leg)
dev.off()

SF2D <-DimPlot(all.samples, group.by = "mouse.id", reduction="umap",pt.size=0.05,cols = batch.cols.2,  ncol = 2)+
  theme(axis.text=element_text(size=24), axis.title=element_text(size=24))+
  NoLegend()

pdf(file.path(figures_dir,"SF2D.pdf"),7,7)
print(SF2D)
dev.off()


### Fig 1C - marker dot plot ### 

# maybe read in markers
# markers<-read.table(file.path(dir,"markers.txt"),stringsAsFactors = FALSE)

# select sig markers for supplementary file 
markers.file<-markers[which(markers$p_val_adj < 0.05 ),]
text<-"Cluster-defining marker genes for caecum IEC populations (control cells, Wilcox test)"
markers.file <- markers.file[,c("cluster","gene","pct.1","pct.2","avg_logFC","p_val","p_val_adj")]
colnames(markers.file) <- c("Cluster","Gene","Proportion of cells in cluster expressing marker", "Proportion of all other cells expressing marker", "Average log FoldChange", "p value", "Adjusted p value")
markers.file.text<-mrna$text_matrix(markers.file,text)
markers.file.text<-mrna$write_results_text(markers.file.text,files_dir,"SFile3.1.tsv")

# Order by logFC
markers<-markers[which(markers$p_val_adj < 0.05 ),]
markers<-markers[order(-abs(markers$avg_logFC)),]

# select top 5 markers from each cluster to plot 
# make sure we also plot known markers, even if not in top 5
known_markers <- read.table('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/known_markers.txt',header=TRUE,stringsAsFactors = FALSE)
clusters<-levels(all.samples@active.ident)
markers_to_plot<-c()
for (cluster in clusters){
  x<-markers[which(markers$cluster == cluster),'gene']
  y<-known_markers[which(known_markers$cluster == cluster),'marker']
  markers_to_plot<-c(markers_to_plot, x[0:5],y)
}
markers_to_plot<-unique(markers_to_plot)
cols.2<-rev(brewer.pal(n=11,name="RdYlBu"))

F1C <- DotPlot(all.samples, features = rev(markers_to_plot), dot.scale = 4) +
       theme(text = element_text(size=9),
          axis.text=element_text(size=9),
          axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
          legend.title = element_text(size = 9), 
          legend.text = element_text(size = 7),
          legend.margin=margin(0,0,0,0), 
          legend.box.margin=margin(0,0,-10,-10) )+
        scale_color_gradientn(name = "Average expression\n(scaled)",   colours = cols.2)+
        guides(color = guide_colourbar(barwidth = 0.5, barheight = 2), 
          size = guide_legend(keywidth = 0.5, keyheight = 0.5))+
        labs(x=NULL,y=NULL)+
        scale_y_discrete(limits = rev(levels(all.samples)))+    # reverse the order of the y axis
        scale_size_continuous(name = "Percentage of cells\nexpressing marker",range = c(0.1,1.5))


pdf(file.path(figures_dir,"F1C.pdf"),10,2.5)
print(F1C)
dev.off()

### F1D - quantify cluster sizes - C, D1, D3 ###
proportions.condition<-as.data.frame(prop.table(table(Idents(all.samples), all.samples$orig.ident),margin=2), stringsAsFactors = FALSE)
colnames(proportions.condition) <- c("celltype","sample","proportion")
proportions.condition<-merge(proportions.condition,meta_data,by.x = "sample",by.y="sample_id")

# change to nice label names
proportions.condition$condition <- revalue(proportions.condition$condition,c("control.24" = "C", "infected.24" = "D1", "control.72" = "C", "infected.72" = "D3"))
proportions.condition$condition<-as.character(proportions.condition$condition)

# Not used but for our info
F1D.box<-ggplot(proportions.condition,aes(condition,proportion,fill=celltype)) +
  #geom_dotplot(binaxis = "y", stackdir = "center")+ 
  geom_boxplot()+
  scale_fill_manual(values=cols)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_wrap(~celltype, scale = "free_y", ncol=11)+
  NoLegend()+
  ylab("Proportion of cells")+
  xlab("")


pdf(file.path(figures_dir,"F1D.box.pdf"),10,2.5)
print(F1D.box)
dev.off()

# calculate mean and standard deviation
calculate_mean <- function(celltype, condition, proportions.table){
  values <- proportions.table[which(proportions.table$celltype == celltype & proportions.table$condition == condition),"proportion"]
  return(mean(values))
}

calculate_sd <- function(celltype, condition, proportions.table){
  values <-  proportions.table[which(proportions.table$celltype == celltype & proportions.table$condition == condition),"proportion"]
  return(sd(values))
}

all.means <- unique(proportions.condition[,c("celltype", "condition")])
all.means$celltype <- factor(all.means$celltype, levels = levels(all.samples))
for (row in 1:nrow(all.means)){
  all.means$mean[row] <- calculate_mean(all.means$celltype[row], all.means$condition[row], proportions.condition)
  all.means$sd[row]   <- calculate_sd(all.means$celltype[row], all.means$condition[row], proportions.condition)
}

F1D<-ggplot(all.means,aes(condition,mean,fill=celltype)) +
  #geom_dotplot(binaxis = "y", stackdir = "center")+ 
  geom_col()+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2)+
  scale_fill_manual(values=cols)+
  theme(axis.text.x = element_text()) + 
  facet_wrap(~celltype, ncol=11)+
  NoLegend()+
  ylab("Proportion of cells")+
  xlab("")+
  ylim(0,0.6)


pdf(file.path(figures_dir,"F1D.pdf"),10,2)
print(F1D)
dev.off()

### t tests of proportions ###
run_t_test <- function(celltype, condition, proportions.table){
  infected <- proportions.table[which(proportions.table$celltype == celltype  & proportions.table$condition == condition),"proportion"]
  control <- proportions.table[which(proportions.table$celltype == celltype  & proportions.table$condition == "C"),"proportion"]
  test <- t.test(x = control, y = infected, alternative = "two.sided")
  return(test$p.value)
}


clusters <- levels(all.samples)
# consider all controls together
t.tests.d1 <- lapply(clusters, run_t_test, "D1", proportions.condition)
names(t.tests.d1) <- clusters
t.tests.d3 <- lapply(clusters, run_t_test, "D3", proportions.condition)
names(t.tests.d3) <- clusters

# time matched controls only
proportions.condition.24 <- proportions.condition[which(proportions.condition$time == "24hpi"),]
tm.t.tests.d1 <- lapply(clusters, run_t_test, "D1", proportions.condition.24)
names(tm.t.tests.d1) <- clusters
proportions.condition.72 <- proportions.condition[which(proportions.condition$time == "72hpi"),]
tm.t.tests.d3 <- lapply(clusters, run_t_test, "D3", proportions.condition.72)
names(tm.t.tests.d3) <- clusters

### Fig SF2F Violin plots of interesting genes ###

violin_genes<-c(
  "Isg15",
  "Ifi27l2b",
  "Ddx60",
  "Ifit1",
  "Ifit1bl1",
  "Irf7",
  "Lypd8",
  "Ifit3",
  "Rsad2",
  "Zbp1",
  "Oasl1",
  "Reg3g",
  "Reg3b"
)



### stacked violin plot code from https://rpubs.com/DarrenVan/628853 ###

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-2, 0, -5, 0), "pt"), 
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... , do.return = TRUE)
  p$layers[[1]]$aes_params$size = 0
  p <- p +
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 7, angle = 90), 
          axis.text.y = element_text(size = 5), 
          plot.margin = plot.margin,
          plot.title = element_blank(),
          axis.line = element_line(size = 0.3)
       )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, size = 7), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

### end of code from https://rpubs.com/DarrenVan/628853 ###

# use 72h samples (and all controls) for violin plots
samples.c.72i <- subset(all.samples,  subset = time.status %in% c('24hpi.control','72hpi.control','72hpi.infected'))
F1D <- StackedVlnPlot(obj = samples.c.72i, features = violin_genes, split.by = "infection_status")

pdf(file.path(figures_dir,"F1D.pdf"),2,6)
print(F1D)
dev.off()

### SF4A - expression of alarmins ###

cols.heat <- colorRampPalette(colors = c("#313695","#FFFFBF","#A50026"))(200)
SF4A <- FeaturePlot(all.samples,features = c("Il33","Il25","Tslp","Ifna1","Ifnb1","Ifng","Ifnl2","Ifnl3"),split.by = "condition", pt.size = 0.1, cols = c('grey',"#A50026"))+
        theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 7, hjust = -0.01))

pdf(file.path(figures_dir,"SF4A.pdf"),12,20)
print(SF4A)
dev.off()

# can't get a legend when using split.by, so write a function for these feature plots

FeaturePlotSplitScale <- function(seurat_object, feature, split_by, ...){
  all.cells <- colnames(seurat_object)
  groups <- levels(seurat_object@meta.data[,split_by])
  
  # extract min and max values for the scale
  
  minimum <- min(seurat_object[['RNA']]@data[feature, ])
  maximum <- max(seurat_object[['RNA']]@data[feature, ])
  if (maximum == 0){
    maximum = 0.1 # arbitrary value to make all 0s grey
  }
  ggs<-list()
  for (group in groups){
    subset_indx<- seurat_object@meta.data[,split_by] == group
    subset_cells<- all.cells[subset_indx]
    new_obj <- subset(all.samples, cells = subset_cells )
    gg <- FeaturePlot(new_obj, features = feature, ...)+
          scale_color_gradient(limits = c(minimum, maximum), low = 'grey', high = "#A50026", oob = squish, name = feature)+
          ggtitle(group)+
          theme(axis.text=element_text(size=9), axis.title=element_text(size=9),legend.text = element_text(size = 9),legend.title = element_text(size = 11, face="bold.italic"), plot.title = element_text(face = "italic", size = 9, hjust = -0.01))+
          guides(color = guide_colourbar(barwidth = 1, barheight = 5))
    ggs[[group]] <- gg
  }
  return(ggs)
}

immune_genes <- c("Il33","Il25","Tslp","Ifna1","Ifnb1","Ifng","Ifnl2","Ifnl3")
for (immune_gene in immune_genes){
  ggs<-FeaturePlotSplitScale(all.samples,immune_gene,"condition", pt.size = 0.1)
  fig <- wrap_plots(ggs, guides='collect')
  pdf(file.path(figures_dir,paste0("SF11.",immune_gene,".pdf")),7.5,2.5)
  print(fig)
  dev.off()
}

### Feature plots of (selected) regulated genes - SF4B ###
de.genes <- scan('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/merge/degenes.txt',sep = "\n", what = character())

SF4B.1 <- FeaturePlot(all.samples,features = de.genes[1:6],split.by = "condition", pt.size = 0.1, cols = cols.heat)+
  theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 7, hjust = -0.01))
pdf(file.path(figures_dir,"SF4B.1.pdf"),16,20)
print(SF4B.1)
dev.off()

SF4B.2 <- FeaturePlot(all.samples,features = de.genes[7:12],split.by = "condition", pt.size = 0.1, cols = cols.heat)+
  theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 7, hjust = -0.01))
pdf(file.path(figures_dir,"SF4B.2.pdf"),16,20)
print(SF4B.2)
dev.off()

SF4B.3 <- FeaturePlot(all.samples,features = de.genes[13:18],split.by = "condition", pt.size = 0.1, cols = cols.heat)+
  theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 7, hjust = -0.01))
pdf(file.path(figures_dir,"SF4B.3.pdf"),16,20)
print(SF4B.3)
dev.off()

SF4B.4 <- FeaturePlot(all.samples,features = de.genes[19:24],split.by = "condition", pt.size = 0.1, cols = cols.heat)+
  theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 7, hjust = -0.01))
pdf(file.path(figures_dir,"SF4B.4.pdf"),16,20)
print(SF4B.4)
dev.off()

SF4B.5 <- FeaturePlot(all.samples,features = de.genes[25:26],split.by = "condition", pt.size = 0.1, cols = cols.heat)+
  theme(axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text = element_text(face = "italic", size = 7, hjust = -0.01))
pdf(file.path(figures_dir,"SF4B.5.pdf"),16,20)
print(SF4B.5)
dev.off()

### GO enrichment of cluster marker genes SF2G ###

# clusterprofiler - first convert gene symbols to Entrez IDs

genesymbol_to_entrez <-function(gene_list){
  ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl",host = "asia.ensembl.org")
  filters <- ('external_gene_name')
  attributes <- c('ensembl_gene_id','entrezgene_id','entrezgene_accession','external_gene_name')
  query <- getBM(attributes = attributes, filters = filters, values = gene_list, mart = ensembl)
  return(query)
}

retrieve_markers<-function(cluster,markers){
  cluster_markers <- markers[which(markers$cluster == cluster & markers$p_val_adj < 0.05), "gene"]
  entrez_markers <- genesymbol_to_entrez(cluster_markers)
  return(entrez_markers)
}

enrich_gos<-function(markers_table){
  go_enrichment <- enrichGO(markers_table$entrezgene_id,  OrgDb=org.Mm.eg.db, qvalueCutoff=1, pvalueCutoff=1, ont = "BP")
  return(go_enrichment)
}

marker_ids<-lapply(clusters, retrieve_markers, markers)
names(marker_ids)<-clusters
go_enrichments<-lapply(marker_ids, enrich_gos)
names(go_enrichments)<-clusters

# function to plot terms
plot_gos <- function(cluster, go_enrichments){
  n<-length(rownames(data.frame(go_enrichments[[cluster]])[which(data.frame(go_enrichments[[cluster]])$p.adjust< 0.05),]))
  pdf(file.path(figures_dir,paste0(cluster,".godots.pdf")),20,40)
  print(dotplot(go_enrichments[[cluster]],showCategory=n))
  dev.off()
  return(cluster)
}

# plot terms for each cluster 
clusters <- sapply(names(go_enrichments), plot_gos, go_enrichments)

# supplementary file of all sig GO terms for each cluster
all.terms <- data.frame()
for (cluster in clusters){
  terms <- data.frame(go_enrichments[[cluster]])
  sig.terms <- terms[which(terms$p.adjust< 0.05 & terms$Count >=5),]
  sig.terms$Cluster <- cluster
  all.terms <- rbind(all.terms,sig.terms)
}
all.terms <- all.terms[,c("Cluster","ID", "Description", "GeneRatio","pvalue","p.adjust")]
colnames(all.terms) <- c("Cluster","GO term ID", "Description", "Gene ratio", "p value", "Adjusted p value")
text <- "GO term (BP) enrichment in IEC cell type cluster marker genes"
cluster.go.text<-mrna$text_matrix(all.terms,text)
cluster.go.text<-mrna$write_results_text(cluster.go.text,files_dir,"SFile3.2.tsv")


pdf("Isg15.pdf", 3, 3)
FeaturePlot(d3.samples,features = c("Isg15"), pt.size = 0.1, cols = c('grey',"#A50026"))+
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9), plot.title = element_text(face = "italic", size = 13, hjust = -0.01))
dev.off()

pdf("Ddx60.pdf", 3, 3)
FeaturePlot(d3.samples,features = c("Ddx60"), pt.size = 0.1, cols = c('grey',"#A50026"))+
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9), plot.title = element_text(face = "italic", size = 13, hjust = -0.01))
dev.off()

pdf("Irf7.pdf", 3, 3)
FeaturePlot(d3.samples,features = c("Irf7"), pt.size = 0.1, cols = c('grey',"#A50026"))+
  theme(axis.text=element_text(size=9), axis.title=element_text(size=9), plot.title = element_text(face = "italic", size = 13, hjust = -0.01))
dev.off()



