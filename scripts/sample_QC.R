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
library(DoubletDecon)

parser <- ArgumentParser()
parser$add_argument("--results_directory",help="directory for files to be written to")
parser$add_argument("--samples", help="File with a list of sample IDs to analyse, one per line")
parser$add_argument("--metadata", help="File with sample metadata")
parser$add_argument("--cellranger_version", help="cellranger131, cellranger211 or cellranger302")
parser$add_argument("--annotation_version", help="mouse reference annotation used in Cellranger",default='mm10-3_0_0')
parser$add_argument("--mito_cutoff", type="integer", help="Filter out cells with % reads aligning to mitochondrial genes greater than this (default: no cutoff)",default=100)
parser$add_argument("--nfeatures_cutoff", type="integer",help="Filter out cells with fewer features than this (default: no filter)", default=0)
parser$add_argument("--ncells_cutoff",type="integer",help="Only include features (genes) that are expressed in more than this number of cells (default: no filter)",default=0)

args <- parser$parse_args()

#enable paralellisation
#plan("multiprocess", workers = 4)
#options(future.globals.maxSize = 9437184000)

#for now
args$cellranger_version<-'cellranger302'
args$metadata<-'/Users/fr7/git_repos/single_cell/metadata/samples.txt'
args$samples<-'/Users/fr7/git_repos/single_cell/experiment_4/samples_to_analyse.txt'
args$results_directory<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL'
figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'

#find samples and meta data
samples<-scan(args$samples,what=character())
meta_data<-read.table(args$metadata,header=T)

#make output directories
dir<-file.path(args$results_directory,'QC')
dir.create(dir)
for (sample in samples){
  dir.create(file.path(dir,sample))
}

#QC:get samples into seurat objects
seurat_objects<-lapply(samples,m$normalize_data,
                       args$cellranger_version,
                       meta_data,
                       2000,
                       args$ncells_cutoff,
                       args$mito_cutoff,
                       args$nfeatures_cutoff,
                       args$annotation_version)
names(seurat_objects)<-samples

#plot UMIs for each sample
umi_frequency_plots<-lapply(names(seurat_objects),m$umi_frequency_plot,seurat_objects,dir,19000) #try various to get the right threshold for each sample

SF2A.UMI <- m$umi_frequency_plot("4672STDY8112878", seurat_objects, dir, 19000)
SF2A.UMI <- SF2A.UMI + ylab("Density") + ggtitle("") + theme(axis.text=element_text(size=7), axis.title=element_text(size=9))
pdf(file.path(figures_dir,"SF2A.UMI.pdf"), 2.5,2.5)
print(SF2A.UMI)
dev.off()

#also plot %mito and ngenes for each sample 
mito_frequency_plots<-lapply(names(seurat_objects),m$mito_frequency_plot,seurat_objects,dir)
all<-m$print_combined(mito_frequency_plots,dir,'mito_frequency_plots')

SF2A.MITO.PRE <- m$mito_frequency_plot("4672STDY8112878", seurat_objects, dir)
SF2A.MITO.PRE <- SF2A.MITO.PRE + ylab("Density") + xlab("Percentage mitochondrial features") + ggtitle("") + theme(axis.text=element_text(size=7), axis.title=element_text(size=9))
pdf(file.path(figures_dir,"SF2A.MITO.PRE.pdf"), 2.5,2.5)
print(SF2A.MITO.PRE)
dev.off()

ngenes_plots<-lapply(names(seurat_objects),m$genes_frequency_plot,seurat_objects,dir)
all<-m$print_combined(ngenes_plots,dir,'ngenes_plots')


#Inspect UMI plots and select a cutoff for each sample
#Apply UMI filters
umi_filters<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/UMI_filters.txt',header=T)
#plot all selected filters
all_p<-list()
for (sample in names(seurat_objects)){
  threshold<-umi_filters[which(umi_filters$sample == sample),'umi_threshold' ]
  p<-m$umi_frequency_plot(sample,seurat_objects,dir, threshold)
  all_p[[sample]]<-p
}
all<-m$print_combined(all_p,dir,'umi_with_thresholds')

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

#new directory for the UMI filtered samples
dir<-file.path(dir,'umi_filtered')
dir.create(dir)
for (sample in samples){
  dir.create(file.path(dir,sample))
}

#redo all plots post umi filtering
umi_frequency_plots<-lapply(names(umi_filtered),m$umi_frequency_plot,umi_filtered,dir,0)
all<-m$print_combined(umi_frequency_plots,dir,'UMI_frequency')

mito_frequency_plots<-lapply(names(umi_filtered),m$mito_frequency_plot,umi_filtered,dir)
all<-m$print_combined(mito_frequency_plots,dir,'mito_frequency')

SF2A.MITO.POST <- m$mito_frequency_plot("4672STDY8112878", umi_filtered, dir)
SF2A.MITO.POST <- SF2A.MITO.POST + ylab("Density") + xlab("Percentage mitochondrial features") + ggtitle("") + theme(axis.text=element_text(size=7), axis.title=element_text(size=9))
pdf(file.path(figures_dir,"SF2A.MITO.POST.pdf"), 2.5,2.5)
print(SF2A.MITO.POST)
dev.off()


genes_frequency_plots<-lapply(names(umi_filtered),m$genes_frequency_plot,umi_filtered,dir)
all<-m$print_combined(genes_frequency_plots,dir,'genes_frequency')


#apply a % mito cutoff too
sub<-function(x) subset(x,subset=percent.mito<=30)
mito_filtered<-lapply(umi_filtered,sub)

#new directory for the mito filtered samples
dir<-file.path(dir,'mito_filtered')
dir.create(dir)
for (sample in samples){
  dir.create(file.path(dir,sample))
  out_dir<-file.path(dir,sample)
  saveRDS(mito_filtered[[sample]],file.path(out_dir,"seurat_object.rds"))
}

#redo all plots post mito filtering
umi_frequency_plots<-lapply(names(mito_filtered),m$umi_frequency_plot,mito_filtered,dir,0)
all<-m$print_combined(umi_frequency_plots,dir,'UMI_frequency')

mito_frequency_plots<-lapply(names(mito_filtered),m$mito_frequency_plot,mito_filtered,dir)
all<-m$print_combined(mito_frequency_plots,dir,'mito_frequency')

genes_frequency_plots<-lapply(names(mito_filtered),m$genes_frequency_plot,mito_filtered,dir)
all<-m$print_combined(genes_frequency_plots,dir,'genes_frequency')

#maybe
mito_filtered<-list()
dir<-file.path(args$results_directory,'QC','umi_filtered','mito_filtered')
for (sample in samples){
  mito_filtered[[sample]]<-readRDS(file.path(dir,sample,'seurat_object.rds'))
}

#merge all samples
all<- merge(mito_filtered[[1]], y = mito_filtered[-1], add.cell.ids = names(mito_filtered), project = "all" )
dir<-file.path(args$results_directory,'QC','combined')

#standard clustering workflow
all <- CellCycleScoring(all, s.features = s.genes, g2m.features = g2m.genes)
all$CC.Difference <- all$S.Score - all$G2M.Score
all<-NormalizeData(all)
all<-FindVariableFeatures(all, selection.method = "vst",nfeatures = 2000)
all<-ScaleData(all, verbose = FALSE, vars.to.regress = c("orig.ident","CC.Difference","nCount_RNA"))
all<- RunPCA(all, verbose = FALSE, npcs=40)
all <- m$do_clustering(all,pcdimensions=35,resolution=0.2,min.dist=0.01)

p<-m$plot_UMAPS(all,dir)
markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Isg15","Ptprc","Reg3g")
p<-m$marker_violins(all,dir,markers_to_plot)
p<-m$markers_feature_plot(all,dir,markers_to_plot)
markers<-m$find_all_markers(all,dir)
marker_dots<-m$marker_dotplot(all,markers,dir)

#identify and merge clusters
new.cluster.ids<-c("0","0","2","3","4","5","6","7","8")
names(new.cluster.ids)<-levels(all$seurat_clusters)
all <- RenameIdents(all,new.cluster.ids)

#recalculate markers
markers<-m$find_all_markers(all,dir)
marker_dots<-m$marker_dotplot(all,markers,dir)

#doublet decon
files<-m$prepare_for_doublet_decon(all,dir,'all')
#run DD 10 times and take the intersection.
results_list<-list()
for (i in 1:10){
  results<-m$run_doublet_decon(files$newExpressionFile,files$newGroupsFile,'all',dir,1,FALSE)
  results_list[[i]]<-results
}

all_doublets <- row.names(results_list[[1]]$Final_doublets_groups)
for (results in results_list){
  doublets<-row.names(results$Final_doublets_groups)
  all_doublets <- intersect(all_doublets,doublets) 
}


all.1<-subset(all,cells=all_doublets,invert=TRUE)
all.1 <- m$do_clustering(all.1,pcdimensions=35,resolution=0.3,min.dist=0.01)

p<-m$plot_UMAPS(all.1,dir)
p<-m$marker_violins(all.1,dir,markers_to_plot)
p<-DimPlot(all.1,cells.highlight=WhichCells(all.1,ident="9"))
markers<-m$find_all_markers(all.1,dir)
marker_dots<-m$marker_dotplot(all.1,markers,dir)

#Also extract the immune cell cluster
immune_cell_markers<-c("Ptprc","Thy1","Cd3e","Cd3d","Insr","Fcgr1","Fcgr3","Fcgr2b","Epcam","Cd80","Cd86")
fibroblast_markers<-c("Vim","Col6a2","Col1a2","Dpt")

pdf(file.path(dir,"immune_cluster_violins.pdf"),20,10)
VlnPlot(all.1,features=immune_cell_markers,pt.size=0,ncol=5)
dev.off()

pdf(file.path(dir,"fibroblast_cluster_violins.pdf"),20,10)
VlnPlot(all.1,features=fibroblast_markers,pt.size=0,ncol=4)
dev.off()


clusters<-levels(all.1)
clusters<-clusters[clusters != "8"]
all.2<-subset(all.1,idents=clusters)

#separate back into original datasets 
Idents(all.2)<-"orig.ident"
for (sample in samples){
  object<-subset(all.2,idents=sample)
  saveRDS(object,file.path(dir,paste0(sample,".rds")))
}

##############


