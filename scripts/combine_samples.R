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

parser <- ArgumentParser()
parser$add_argument("--results_directory",help="directory for files to be written to")
parser$add_argument("--samples", help="File with a list of sample IDs to analyse, one per line")
parser$add_argument("--metadata", help="File with sample metadata")
args <- parser$parse_args()
args$samples<-'/Users/fr7/git_repos/single_cell/experiment_4/samples_to_analyse.txt'
args$metadata<-'/Users/fr7/git_repos/single_cell/metadata/samples.txt'
dir<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/combined'
samples<-scan(args$samples,what=character())
meta_data<-read.table(args$metadata,header=T)
sample_dir<-'/Users/fr7/git_repos/single_cell/experiment_4/FINAL/QC/umi_filtered/mito_filtered'

#maybe
#combined<-readRDS(file.path(dir,"seurat_object.rds"))

seurat_objects<-list()
for (sample in samples){
  object<-readRDS(file.path(sample_dir,sample,'seurat_object.rds'))
  seurat_objects[[sample]]<-object
}

#scale and find variable features
for (i in names(seurat_objects)){
  seurat_objects[[i]]<-NormalizeData(seurat_objects[[i]])
  seurat_objects[[i]] <- FindVariableFeatures(seurat_objects[[i]], selection.method = "vst", nfeatures = 2000)
}

dir<-file.path(dir,'QC')
dir.create(dir)

#first integrate all objects, to check for the immune cell cluster and doublets
combined<-m$integrate_datasets(seurat_objects,30)
p<-m$elbow_plot(combined,30,"pca", dir)
combined<-m$do_clustering(combined,25,0.4)
markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Isg15","Cd3e","Ptprc")
p1<-m$QC_violins(combined, dir)
p2<-m$QC_scatters(combined,dir)
p3<-m$plot_UMAPS(combined,dir)
p4<-m$marker_violins(combined,dir,markers_to_plot)
saveRDS(combined,file.path(dir,'seurat_object.rds'))

#doublet decon
DefaultAssay(combined) <- "RNA"
files<-m$prepare_for_doublet_decon(combined,dir,'combined')
results<-m$run_doublet_decon(files$newExpressionFile,files$newGroupsFile,'combined',dir,1.04,TRUE)
doublets<-row.names(results$Final_doublets_groups)
test<-subset(combined,cells=doublets,invert=TRUE)

dir<-file.path(dir,'post_doublet_decon')
dir.create(dir)

test <- RunUMAP(test, reduction = "pca", dims = 1:30,min.dist=0.01)
test <- FindNeighbors(test, reduction = "pca", dims = 1:30)
test <- FindClusters(test, resolution = 3)

#remove immune cell cluster
test<-subset(test,idents="45",inverse=TRUE)

p3<-m$plot_UMAPS(test,dir)
p4<-m$marker_violins(test,dir,markers_to_plot)
DimPlot(test,reduction="umap",cells.highlight = WhichCells(test, idents = ""))
FeaturePlot(test,features=c("Dclk1"))
tufts<-subset(test,idents=c("67"))
tufts <- RunUMAP(tufts, reduction = "pca", dims = 1:25)
tufts <- FindNeighbors(tufts, reduction = "pca", dims = 1:25)
tufts <- FindClusters(tufts, resolution = 1)
p3<-m$plot_UMAPS(tufts,dir)
tuft.markers<-m$find_all_markers(tufts,dir)
FeaturePlot(tufts,features=c("Eif3a"))
FeaturePlot(test,features=c("Eif3a","Dclk1","Cdk4"))
p1<-m$QC_violins(tufts, dir)


###
markers<-m$find_all_markers(combined,dir)
marker_dots<-m$marker_dotplot(combined,markers,dir)



combined_r1<-m$do_clustering(combined,25,1)
saveRDS(combined_r1,file.path(dir,"seurat_object.rds"))
markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Isg15")
p1<-m$QC_violins(combined_r1, dir)
p2<-m$QC_scatters(combined_r1,dir)
p3<-m$plot_UMAPS(combined_r1,dir)
p4<-m$marker_violins(combined_r1,dir,markers_to_plot)
markers<-m$find_all_markers(combined_r1,dir)
marker_dots<-m$marker_dotplot(combined_r1,markers,dir)

immune_cell_markers<-c("Ptprc","Thy1","Cd3e","Cd3d","Insr","Fcgr1","Fcgr3","Fcgr2b","Epcam","Cd80","Cd86")
fibroblast_markers<-c("Vim","Col6a2","Col1a2","Dpt")

pdf(file.path(dir,"immune_cluster_violins.pdf"),20,10)
VlnPlot(combined_r1,features=immune_cell_markers,pt.size=0,ncol=5)
dev.off()

pdf(file.path(dir,"fibroblast_cluster_violins.pdf"),20,10)
VlnPlot(combined_r1,features=fibroblast_markers,pt.size=0,ncol=4)
dev.off()

#cluster 22 identified as contaminant immune cluster
immune<-WhichCells(combined_r1, idents = "22")
pdf(file.path(dir,"immune_highlighted.pdf"),6,5)
DimPlot(combined_r1,reduction="umap",cells.highlight = immune)
dev.off()

#Remove cluster 22
clusters<-levels(combined_r1)
clusters<-clusters[clusters != "22"]
combined_r1_immune_removed<-subset(combined_r1,idents=clusters)

#separate back into original datasets
seurat_objects_immune_removed<-list()
Idents(combined_r1_immune_removed)<-"orig.ident"
for (sample in samples){
  object<-subset(combined_r1_immune_removed,idents=sample)
  seurat_objects_immune_removed[[sample]]<-object
}

#count how many cells have been lost in each sample
for (sample in samples){
  original<-length(seurat_objects[[sample]]$orig.ident)
  current<-length(seurat_objects_immune_removed[[sample]]$orig.ident)
  print( c(sample,original,current) )
}

#rescale and integrate the control samples
control_objects<-list()
for (i in names(seurat_objects_immune_removed)){
  if (unique(seurat_objects_immune_removed[[i]]$infection_status == 'control')){
    control_objects[[i]]<-seurat_objects_immune_removed[[i]]
  }
}

for (i in names(control_objects)){
  DefaultAssay(control_objects[[i]]) <- "RNA"
  control_objects[[i]]<-NormalizeData(control_objects[[i]])
  control_objects[[i]] <- FindVariableFeatures(control_objects[[i]], selection.method = "vst", nfeatures = 2000)
}

combined_control<-m$integrate_datasets(control_objects,30)
ElbowPlot(combined_control,ndims = 30,reduction = "pca")
dir<-file.path(dir,'control')
dir.create(dir)

resolution<-1.5
#testing
dir<-file.path(dir,resolution)
dir.create(dir)
combined_control <- RunUMAP(combined_control, reduction = "pca", dims = 1:25,min.dist=0.001)
combined_control <- FindNeighbors(combined_control, reduction = "pca", dims = 1:25)
combined_control <- FindClusters(combined_control, resolution = 3)
p3<-m$plot_UMAPS(combined_control,dir)
p4<-m$marker_violins(combined_control,dir,markers_to_plot)

FeaturePlot(combined_control,features=c("Chga","Dclk1","Chgb"))
DimPlot(combined_control,reduction="umap",cells.highlight = WhichCells(combined_control,idents="36"))

combined_control<-m$do_clustering(combined_control,25,resolution)
dir<-file.path(dir,resolution)
dir.create(dir)

p1<-m$QC_violins(combined_control, dir)
p2<-m$QC_scatters(combined_control,dir)
p3<-m$plot_UMAPS(combined_control,dir)
p4<-m$marker_violins(combined_control,dir,markers_to_plot)
markers<-m$find_all_markers(combined_control,dir)
marker_dots<-m$marker_dotplot(combined_control,markers,dir)

#######

combined.2<-m$do_clustering(combined,25,1)
combined<-m$do_clustering(combined,25,1)

p3<-m$plot_UMAPS(combined,dir)
p4<-m$marker_violins(combined,dir,markers_to_plot)

combined <- RunUMAP(combined, reduction = "pca", dims = 1:10)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:10)
combined <- FindClusters(combined, resolution = 0.15)
p3<-m$plot_UMAPS(combined,dir)
p4<-m$marker_violins(combined,dir,markers_to_plot)



####### 


combined_r2<-m$do_clustering(combined,25,2)
combined_r0.2<-m$do_clustering(combined,25,0.2)
combined_r0.1<-m$do_clustering(combined,25,0.1)
saveRDS(combined,file.path(dir,"seurat_object.rds"))
markers_to_plot<-c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Isg15")
p1<-m$QC_violins(combined, dir)
p2<-m$QC_scatters(combined,dir)
p3<-m$plot_UMAPS(combined,dir)
p4<-m$marker_violins(combined,dir,markers_to_plot)

#find cluster markers
markers<-m$find_all_markers(combined,dir)
marker_dots<-m$marker_dotplot(combined,markers,dir)
#plot cell numbers per cluster
m$plot_cluster_distributions(combined,dir,meta_data)

#find regulated genes
Idents(combined)<-'time'
combined_24h<-subset(combined,idents="24hpi")
Idents(combined_24h)<-'seurat_clusters'
combined_72h<-subset(combined,idents="72hpi")
Idents(combined_72h)<-'seurat_clusters'
clusters<-levels(combined_24h)
regulated_24<-lapply(clusters,m$find_regulated_genes,combined_24h,dir,'24hpi')
names(regulated_24)<-clusters
regulated_72<-lapply(clusters,m$find_regulated_genes,combined_72h,dir,'72hpi')
names(regulated_72)<-clusters

p1<-m$plot_all_regulated(combined_24h, dir, regulated_24)
p2<-m$plot_all_regulated(combined_72h, dir, regulated_72)

#highlight immune cell cluster
immune<-WhichCells(combined, idents = "5")
DimPlot(combined,reduction="umap",cells.highlight = immune)
pdf(file.path(dir,"immune_highlighted.pdf"),6,5)
DimPlot(combined,reduction="umap",cells.highlight = immune)
dev.off()

pdf(file.path(dir,"immune_cluster_violins.pdf"),10,10)
VlnPlot(combined,features=c("Ptprc","Cd3e","Cd3d","Insr","Epcam","Krt20","Rgs1","Emb","Ccl5"),pt.size=0,ncol=3)
dev.off()

#remove those cells
combined<-subset(combined,idents=c("0","1","2","3","4","6"))

#separate back into original datasets and reintegrate
seurat_objects<-list()
Idents(combined)<-"orig.ident"
for (sample in samples){
  object<-subset(combined,idents=sample)
  seurat_objects[[sample]]<-object
  
}

for (i in names(seurat_objects)){
  DefaultAssay(seurat_objects[[i]]) <- "RNA"
  seurat_objects[[i]]<-NormalizeData(seurat_objects[[i]])
  seurat_objects[[i]] <- FindVariableFeatures(seurat_objects[[i]], selection.method = "vst", nfeatures = 2000)
}
combined.2<-m$integrate_datasets(seurat_objects,30)
combined.2<-m$do_clustering(combined,25,1)
combined<-m$do_clustering(combined,25,1)

p3<-m$plot_UMAPS(combined,dir)
p4<-m$marker_violins(combined,dir,markers_to_plot)
ElbowPlot(combined,ndims = 30,reduction = "pca")
combined <- RunUMAP(combined, reduction = "pca", dims = 1:25)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:25)
combined <- FindClusters(combined, resolution = 0.15)
p3<-m$plot_UMAPS(combined,dir)
p4<-m$marker_violins(combined,dir,markers_to_plot)

#go with 0.15 initially
pdf(file.path(dir,"batch_effects.pdf"),6,5)
DimPlot(combined,reduction="umap",group.by="orig.ident")
dev.off()

#name the clusters
cluster_ids<-read.table(file.path(dir,"cluster_ids.txt"),header=TRUE,stringsAsFactors=FALSE)
new.cluster.ids<-cluster_ids$id
names(new.cluster.ids)<-cluster_ids$cluster
combined <- RenameIdents(combined,new.cluster.ids)
p3<-m$plot_UMAPS(combined,dir)
p4<-m$marker_violins(combined,dir,markers_to_plot)

#highlight cell types
enterocytes<-WhichCells(combined, idents = c("ENTEROCYTE.1","ENTEROCYTE.2","ENTEROCYTE.3","ENTEROCYTE.4"))
DimPlot(combined,reduction="umap",cells.highlight = enterocytes)
pdf(file.path(dir,"enterocyte_highlighted.pdf"),6,5)
DimPlot(combined,reduction="umap",cells.highlight = enterocytes)
dev.off()

goblets<-WhichCells(combined, idents = c("GOBLET"))
pdf(file.path(dir,"goblets_highlighted.pdf"),6,5)
DimPlot(combined,reduction="umap",cells.highlight = goblets)
dev.off()

TA<-WhichCells(combined, idents = c("TA"))
pdf(file.path(dir,"TA_highlighted.pdf"),6,5)
DimPlot(combined,reduction="umap",cells.highlight = TA)
dev.off()

pdf(file.path(dir,"labelled_umap.pdf"),6,5)
DimPlot(combined,reduction="umap")
dev.off()

markers_to_plot<-c("Muc2","Chga","Dclk1","Aqp8","Krt20","Cdk4", "Lgr5" )
pdf(file.path(dir,"known_violins.pdf"),20,20)
VlnPlot(combined,features=markers_to_plot,pt.size=0,ncol=3)
dev.off()

#find cluster markers
markers<-m$find_all_markers(combined,dir)
marker_dots<-m$marker_dotplot(combined,markers,dir)
#plot cell numbers per cluster
m$plot_cluster_distributions(combined,dir,meta_data)

#virus markers
virus_markers<-c("Isg15","Ifit1","Ifit1bl1","Ddx60","Oasl1","Ifit3")
p<-m$markers_feature_plot(combined,dir,virus_markers)

pdf(file.path(dir,"virus_violins.pdf"),10,10)
VlnPlot(combined, features = virus_markers, pt.size=0, ncol=3)
dev.off()

enterocyte.4<-WhichCells(combined, idents = c("ENTEROCYTE.4"))
pdf(file.path(dir,"enterocyte.4_highlighted.pdf"),5,5)
DimPlot(combined,reduction="umap",cells.highlight = enterocyte.4,split.by = "time.status",ncol=2,sizes.highlight = 0.2)
dev.off()


m$plot_cluster_distributions(combined,dir,meta_data)


#retrieve a table of 
numbers<-as.data.frame(table(Idents(combined), combined$orig.ident))
colnames(numbers)<-c("cluster","sample","cells")
totals<-as.data.frame(table(combined$orig.ident))
colnames(totals)<-c("sample","total_cells")
clusters<-levels(numbers$cluster)
df<-merge(numbers,meta_data,by.x="sample",by.y="sample_id")
df<-merge(df,totals,by="sample")
df$cluster<-as.character(df$cluster)

for (cluster in clusters){
  temp<-df[which(df$cluster==cluster),]
  temp24<-temp[which(temp$time=="24hpi"),]
  model24<-glm(cells ~ treatment + offset(log(total_cells)),family=poisson(link=log),data=temp24)
  temp72<-temp[which(temp$time=="72hpi"),]
  model72<-glm(cells ~ treatment + offset(log(total_cells)),family=poisson(link=log),data=temp72)
  sink(file.path(dir,paste0(cluster,".24.glm.txt")))
  print(temp24)
  print(summary(model24))
  sink()
  sink(file.path(dir,paste0(cluster,".72.glm.txt")))
  print(temp72)
  print(summary(model72))
  sink()  
}


#find regulated genes
Idents(combined)<-'time'
combined_24h<-subset(combined,idents="24hpi")
Idents(combined_24h)<-'seurat_clusters'
combined_24h <- RenameIdents(combined_24h,new.cluster.ids)
combined_72h<-subset(combined,idents="72hpi")
Idents(combined_72h)<-'seurat_clusters'
combined_72h <- RenameIdents(combined_72h,new.cluster.ids)
clusters<-levels(combined_24h)
regulated_24<-lapply(clusters,m$find_regulated_genes,combined_24h,dir,'24hpi')
names(regulated_24)<-clusters
regulated_72<-lapply(clusters,m$find_regulated_genes,combined_72h,dir,'72hpi')
names(regulated_72)<-clusters

p1<-m$plot_all_regulated(combined_24h, dir, regulated_24)
p2<-m$plot_all_regulated(combined_72h, dir, regulated_72)

sig_24h<-c("2010109I03Rik","Adm","Aldh1b1","Anxa1","Atf3","Bok","Car2","Ccdc71l","Cd38","Cd9","Ces1g","Clec2h","Cox7a1","Cyp2c55","Dstn","Edn1","F3","Fhl2","Fos","Gadd45a","Gdf15","Higd1a","Hmox1","Hs3st1","Hsp90aa1","Hspa1a","Hspa1b","Hspe1","Il18","Il1rn","Itga2","Jund","Lmo7","Mal","Mast4","Mt1","Mt2","Pim3","Plaur","Plet1","Reg3b","S100a11","S100a14","S100a16","S100a6","Saa1","Sdcbp2","Selenop","Sgk1","Slc17a4","Slc37a2","Sprr1a","Sprr2h","Tagln2","Tfrc","Tjp2","Tnip3","Txnip") 
DotPlot(combined_24h, features = sig_24h, cols = c("blue", "red"), dot.scale = 4, split.by = "infection_status")+RotatedAxis()
sig_72h<-c("2010109I03Rik","Acox1","Akr1c19","Apol9a","Apol9b","Aqp8","Atp5g1","B3galt5","Bst2","Ces2a","Cox7a1","Ddx60","Fgfbp1","Gm26917","Gprc5a","Gsdmc2","Gsdmc4","Hist1h1c","Hmgcs2","Hmox1","Ifi27","Ifi27l2a","Ifi27l2b","Ifi35","Ifi47","Ifit1","Ifit1bl1","Ifit2","Ifit3","Ifit3b","Irf7","Irgm1","Isg15","Krt19","Lgals3bp","Lgals9","Ly6a","Ly6e","Lypd8","Mal","Mrpl30","Oas1a","Oasl2","Phf11d","Plac8","Psmb8","Psmb9","Rbp2","Rnf213","Rtp4","S100a11","S100g","Saa1","Selenop","Spink4","Stat1","Trim15","Trim30a","Trim30d","Trpv6","Usp18","Xaf1","Zbp1","mt-Atp6","mt-Cytb","mt-Nd4l")
DotPlot(combined_72h, features = sig_72h, cols = c("blue", "red"), dot.scale = 4, split.by = "infection_status")+RotatedAxis()



