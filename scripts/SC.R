import("Seurat")
import("ggplot2")
import("utils")
import("DESeq2")
import("grDevices")
import("cowplot")
import("DoubletDecon")

evaluate<-function(text){
  x<-eval(parse(text=text))
  return(x)
}

normalize_data<-function(sample,cellranger_version,meta_data,variable_features,ncells,nmito,nfeatures,annotation){
  data <- Read10X(data.dir = paste0("/Users/fr7/git_repos/single_cell/count_matrices/",annotation,"/",cellranger_version,"/",sample,"/filtered_feature_bc_matrix"),unique.features = TRUE)
  seurat_object <- CreateSeuratObject(counts=data, project = paste0(sample),min.cells=ncells,min.features=nfeatures)
  seurat_object[["infection_status"]] <- as.character(meta_data[which(meta_data$sample_id == sample),'treatment']) 
  seurat_object[["experiment"]] <- as.character(meta_data[which(meta_data$sample_id == sample),'experiment'])
  seurat_object[["cellranger_version"]]<-cellranger_version
  seurat_object[["percent.mito"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object[["time"]]<-as.character(meta_data[which(meta_data$sample_id == sample),'time']) 
  seurat_object[["time.status"]]<-paste0(seurat_object$time,".",seurat_object$infection_status)
  seurat_object <- subset(seurat_object, subset.name = "percent.mito", high.threshold = nmito)
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst",nfeatures = variable_features)
  return(seurat_object)
}

variable_feature_plot<-function(seurat_object,out_dir){
  top10 <- head(VariableFeatures(seurat_object), 10)
  p1 <- VariableFeaturePlot(seurat_object)
  p1 <- LabelPoints(p1, points = top10, repel = TRUE)
  pdf(file.path(out_dir,"variable_features_plot.pdf"))
  print(p1)
  dev.off()
  return(p1)
}

run_pca<-function(seurat_object,ndims){
  seurat_object <- ScaleData(seurat_object, verbose = FALSE)
  seurat_object <- RunPCA(seurat_object, npcs = ndims, verbose = FALSE)
  return(seurat_object)
}

normalize_and_run_pca<-function(seurat_object,npcs,out_dir){
	seurat_object <- NormalizeData(seurat_object)
	seurat_object<- FindVariableFeatures(seurat_object, selection.method = "vst",nfeatures = 2000)
	seurat_object <- ScaleData(seurat_object, verbose = FALSE)
	seurat_object <- RunPCA(seurat_object, npcs = npcs, verbose = FALSE)
	p <- elbow_plot(seurat_object,npcs,"pca",out_dir)
	return(seurat_object)
}


elbow_plot<-function(seurat_object,ndims,reduction,out_dir){
  p1<-ElbowPlot(seurat_object,ndims = ndims,reduction = reduction)
  pdf(file.path(out_dir,"elbowplot.pdf"))
  print(p1)
  dev.off()
  return(p1)
}

integrate_datasets<-function(list_of_seurat_objects,ndims){        
  anchors<-FindIntegrationAnchors(list_of_seurat_objects, dims = 1:35)
  combined<- IntegrateData(anchorset = anchors, dims = 1:35)
  DefaultAssay(combined) <- "integrated"
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = ndims, verbose = FALSE)
  return(combined)
}  

do_clustering<-function(combined,pcdimensions,resolution,min.dist){
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:pcdimensions, min.dist=min.dist)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:pcdimensions)
  combined <- FindClusters(combined, resolution = resolution)
  return(combined)
}

plot_UMAPS<-function(seurat_object,out_dir){
  p1<-DimPlot(seurat_object,reduction="umap")
  p2<-DimPlot(seurat_object,reduction="umap",group.by = "orig.ident")
  p3<-FeaturePlot(seurat_object,features="percent.mito")
  p4<-FeaturePlot(seurat_object,features = "nFeature_RNA")
  all<-CombinePlots(plots=list(p1,p2,p3,p4),ncol=2)
  pdf(file.path(out_dir,"UMAPs.pdf"),20,20)
  print(all)
  dev.off()
  return(all)
}

split_UMAP<-function(seurat_object,out_dir){
  p1<-DimPlot(seurat_object,reduction="umap",split.by="time.status")
  pdf(file.path(out_dir,"split_UMAP.pdf"),50,20)
  print(p1)
  dev.off()
  return(p1)
}

find_markers<-function(cluster,seurat_object,data_dir){
  DefaultAssay(seurat_object) <- "RNA"
  markers<-FindMarkers(seurat_object, ident.1 = cluster, print.bar = FALSE)
  file<-paste0(data_dir,"/cluster.",cluster,".markers.tsv")
  write.table(markers, file = file, sep ='\t', eol = '\n') 
  return(markers)
}

find_regulated_genes<-function(cluster,seurat_object,data_dir,timepoint){
  x<-subset(seurat_object, idents =  cluster) 
  Idents(x)<-"infection_status"
  DefaultAssay(x) <- "RNA"
  infection.response<-try(FindMarkers(x,ident.1="control",ident.2="infected",print.bar=FALSE),silent=TRUE)
  if (is.data.frame(infection.response)){
   file<-paste0(data_dir,"/",timepoint,".cluster.",cluster,".regulated.tsv")
   write.table(infection.response, file = file, sep ='\t', eol = '\n') 
  }
  else{
    print(infection.response)
  }
  return(infection.response)
}

#plotting functions
get_dotplot<-function(data,x,y,fill,metric){
  g<-ggplot(data,aes(x,y, colour=fill)) +
    geom_jitter(width = 0.1, size = 3) +
    labs(title= metric, x= deparse(substitute(x)) , y = metric)
  return(g)
}

boxplot<-function(data,x,y,fill,metric){
  g<-ggplot(data,aes(x,y, colour=fill)) +
    geom_boxplot() +
    labs(title= metric, x= deparse(substitute(x)) , y = metric)
  return(g)
}

boxplot<-function(data,x,y,fill,metric){
  g<-ggplot(data,aes(x,y, colour=fill)) +
    geom_boxplot() +
    labs(title= metric, x= deparse(substitute(x)) , y = metric) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(g)
}

plot_proportions<-function(cluster,proportions_table,x){
    my_data<-proportions_table[which(proportions_table$cluster == cluster),]
    x<-my_data[,x]
 #   plot<-get_dotplot(my_data,x,my_data$n,NULL,paste0("cluster.",cluster))
    plot<-boxplot(my_data,x,my_data$n,NULL,paste0(cluster))
    return(plot)
}

plot_cluster_distributions<-function(seurat_object,out_dir,meta_data){
  clusters<-levels(seurat_object@active.ident)
  numbers<-as.data.frame(table(Idents(seurat_object), seurat_object$orig.ident))
  proportions_by_sample<-as.data.frame(prop.table(table(Idents(seurat_object), seurat_object$orig.ident),margin=2))
  proportions_by_cluster<-as.data.frame(prop.table(table(Idents(seurat_object), seurat_object$orig.ident),margin=1))
  y<-list(numbers,proportions_by_sample,proportions_by_cluster)
  names(y)<-c("number","proportions_by_sample","proportions_by_cluster")
  for (x in names(y)){
    names(y[[x]])<-c("cluster","sample","n")
    y[[x]]<-merge(y[[x]],meta_data,by.x = "sample",by.y="sample_id")
    y[[x]]$time.status<-paste0(y[[x]]$time,"-",y[[x]]$treatment)
    plots<-lapply(clusters, plot_proportions, y[[x]],"time.status")
    all_plots<-plot_grid(plotlist=plots,ncol=8)
    pdf(file.path(out_dir,paste0(x,'.pdf')),30,5)
    print(all_plots)
    dev.off()
  }
}
  

set_colour_by_experiment<-function(seurat_object){
  if (unique(seurat_object$experiment) == '1'){col="red"}
  if (unique(seurat_object$experiment) == '2'){col="forestgreen"}
  if (unique(seurat_object$experiment == '3.1') | unique(seurat_object$experiment) == '3.2'){col="cornflowerblue"}
  if (unique(seurat_object$experiment) == '4'){col="orange"}
  return(col)
}

set_colour_by_time_condition<-function(seurat_object){
  if (unique(seurat_object$time.status) == '24hpi.control'){col="red"}
  if (unique(seurat_object$time.status) == '24hpi.infected'){col="forestgreen"}
  if (unique(seurat_object$time.status) == '72hpi.control'){col="cornflowerblue"}
  if (unique(seurat_object$time.status) == '72hpi.infected'){col="orange"}
  return(col)
}

violin<-function(seurat_object,feature,group,ymax){
  col<-set_colour_by_experiment(seurat_object)
  vln<-VlnPlot(seurat_object, features = c(feature), cols=c(col), pt.size = 0,group.by = group, do.return=T,y.max=ymax)
  vln<-vln + theme(text = element_text(size=30), axis.text.x=element_text(angle=0,hjust=0.5,size=25),axis.text.y=element_text(size=25),legend.position = "none") + labs(x="",y="",title=paste0(seurat_object$infection_status)) 
  return(vln)
}

scatter<-function(seurat_object,feature1,feature2,xlim=NULL,ylim=NULL){
  scatter_plot<-FeatureScatter(seurat_object,feature1=feature1,feature2 = feature2)
  scatter_plot<-scatter_plot+ theme(text = element_text(size=30), axis.text.x=element_text(angle=0,hjust=0.5,size=25),axis.text.y=element_text(size=25),legend.position = "none") + labs(title=paste0(seurat_object$orig.ident)) 
#  if (!is.null(xlim)) scatter_plot <- scatter_plot + xlim(xlim)
#  if (!is.null(ylim)) scatter_plot <- scatter_plot + ylim(ylim)
  return(scatter_plot)
}

find_all_markers<-function(seurat_object,out_dir){
  markers<-FindAllMarkers(seurat_object,min.pct = 0.3,logfc.threshold = 0.5,only.pos=TRUE)
  write.table(markers, file = file.path(out_dir,'markers.txt'), sep ='\t', eol = '\n') 
  return(markers)
}

QC_violins<-function(seurat_object,out_dir){
  p1<-VlnPlot(seurat_object,features=c("percent.mito","nFeature_RNA","nCount_RNA"),pt.size=0,ncol=1)
  pdf(file.path(out_dir,"QC_violins.pdf"),20,10)
  print(p1)
  dev.off()
  p2<-VlnPlot(seurat_object,features=c("percent.mito","nFeature_RNA","nCount_RNA"),pt.size=0,ncol=1,log=TRUE)
  pdf(file.path(out_dir,"QC_violins_logy.pdf"),20,10)
  print(p2)
  dev.off()
  p3<-VlnPlot(seurat_object,features=c("percent.mito","nFeature_RNA","nCount_RNA"),pt.size=0,ncol=1,sort=TRUE)
  pdf(file.path(out_dir,"QC_violins_sorted.pdf"),20,10)
  print(p3)
  dev.off()
}

QC_scatters<-function(seurat_object,out_dir){
  p1<-scatter(seurat_object,"nCount_RNA","nFeature_RNA")
  p2<-scatter(seurat_object,"nCount_RNA","percent.mito")
  p3<-scatter(seurat_object,"nFeature_RNA","percent.mito")
  p4<-CombinePlots(plots = list(p1,p2,p3),ncol=3)
  pdf(file.path(out_dir,"QC_scatters.pdf"),35,10)
  print(p4)
  dev.off()
  return(p4)
}

marker_violins<-function(seurat_object,out_dir,markers){
  p1<-VlnPlot(seurat_object,features=markers,pt.size=0,ncol=1)
  pdf(file.path(out_dir,"violin_markers.pdf"),20,50)
  print(p1)
  dev.off()
}
  
remove_mito_clusters<-function(seurat_object,to_remove_table){
  sample<-levels(seurat_object$orig.ident)
  to_remove<-to_remove_table[which(to_remove_table$sample == sample),'cluster']
  try(seurat_object<-subset(seurat_object, idents = to_remove, invert = TRUE))
  return(seurat_object)
}

counts_filter<-function(seurat_object,thresholds){
  sample<-levels(seurat_object$orig.ident)
  print(sample)
  threshold<-thresholds[which(thresholds$sample == sample),'threshold']
  print(threshold)
  seurat_object
  try(seurat_object<-subset(seurat_object, subset = nCount_RNA < threshold))
  return(seurat_object)
}

marker_dotplot<-function(seurat_object,markers,out_dir){
  DefaultAssay(seurat_object) <- "RNA"
  clusters<-levels(seurat_object@active.ident)
  markers_to_plot<-c()
  for (cluster in clusters){
    x<-markers[which(markers$cluster == cluster),'gene']
    markers_to_plot<-c(markers_to_plot, x[0:5] )
  }
  markers_to_plot<-unique(markers_to_plot)
  p<-DotPlot(seurat_object, features = rev(markers_to_plot), dot.scale = 4)
  p<-p+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
  pdf(file.path(out_dir,"markers_dots.pdf"),15,5)
  print(p)
  dev.off()
  return(p)
}
  
markers_feature_plot<-function(seurat_object,out_dir,markers){
  markers<-FeaturePlot(seurat_object, features = markers, min.cutoff="q9",pt.size = 0.1)
  pdf(file.path(out_dir,"markers_feature_plot.pdf"),20,20)
  print(markers)
  dev.off()
}

plot_regulated<-function(cluster,seurat_object,regulated){
  temp<-subset(seurat_object,idents = cluster)
  Idents(temp) <- "infection_status"
  avg <- log1p(AverageExpression(temp, verbose = FALSE)$RNA)
  avg$gene <- rownames(avg)
  genes.to.label<-rownames(regulated[[cluster]][which(regulated[[cluster]]$p_val_adj<0.05 & abs(regulated[[cluster]]$avg_logFC)>0.5),])
  p1 <- ggplot(avg, aes(control, infected)) + geom_point() + ggtitle(paste0(cluster))
  if (length(genes.to.label) > 0){
    p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  }
  return(p1)
}

plot_all_regulated<-function(seurat_object,out_dir,regulated){
  clusters<-levels(seurat_object)
  plots<-lapply(clusters,plot_regulated,seurat_object,regulated)
  all_plots<-plot_grid(plotlist=plots,ncol=8)
  pdf(file.path(out_dir,'regulated_plots.pdf'),25,5)
  print(all_plots)
  dev.off() 
  return(all_plots)
}

umi_frequency_plot<-function(sample_id,seurat_objects,out_dir, vline){
  UMIs<-seurat_objects[[sample_id]]@meta.data$nCount_RNA
  UMI.df<-data.frame(UMIs=UMIs)
  myfill<-set_colour_by_time_condition(seurat_objects[[sample_id]])
  p<-ggplot(UMI.df,aes(UMIs)) + geom_density(fill=myfill) + ggtitle(sample_id) 
  p<- p + geom_vline(xintercept = vline, linetype="dotted", size=1)
  pdf(file.path(out_dir,sample_id,paste0('UMI_frequency_plot_',vline,'.pdf')))
  print(p)
  dev.off()
  return(p)
}

mito_frequency_plot<-function(sample_id,seurat_objects,out_dir){
  mito<-seurat_objects[[sample_id]]@meta.data$percent.mito
  mito.df<-data.frame(Percent_mito=mito)
  myfill<-set_colour_by_time_condition(seurat_objects[[sample_id]])
  p<-ggplot(mito.df,aes(Percent_mito)) + geom_density(fill=myfill) + ggtitle(sample_id) 
  pdf(file.path(out_dir,sample_id,'Mito_frequency_plot.pdf'))
  print(p)
  dev.off()
  return(p)
}

genes_frequency_plot<-function(sample_id,seurat_objects,out_dir){
  genes<-seurat_objects[[sample_id]]@meta.data$nFeature_RNA
  genes.df<-data.frame(nGenes=genes)
  myfill<-set_colour_by_time_condition(seurat_objects[[sample_id]])
  p<-ggplot(genes.df,aes(nGenes)) + geom_density(fill=myfill) + ggtitle(sample_id) 
  pdf(file.path(out_dir,sample_id,'nGenes_frequency_plot.pdf'))
  print(p)
  dev.off()
  return(p)
}

print_combined<-function(plotlist,out_dir,title){
  all<-plot_grid(plotlist=plotlist,ncol=4,nrow=4)
  pdf(file.path(out_dir,paste0(title,'.pdf')),30,30)
  print(all)
  dev.off()
  return(all)
}

prepare_for_doublet_decon<-function(seurat_object,out_dir,sample){
  
  newFiles=Improved_Seurat_Pre_Process(seurat_object, num_genes=50, write_files=FALSE)
  write.table(newFiles$newExpressionFile, file.path(out_dir, paste0(sample, "_expression")), sep="\t")
  write.table(newFiles$newFullExpressionFile, file.path(out_dir, paste0(sample, "_fullExpression")), sep="\t")
  write.table(newFiles$newGroupsFile, file.path(out_dir, paste0(sample, "_groups")), sep="\t", col.names = F)
 
  return(newFiles) 
}

run_doublet_decon<-function(expression_file,groups_file,sample,out_dir,rhop,heatmap){
  pdf(file.path(out_dir, paste0(sample,".rhop.",rhop,'.pdf')))
  results<-try(Main_Doublet_Decon(rawDataFile=expression_file, 
                             groupsFile=groups_file, 
                             filename=sample, 
                             location=out_dir,
                             fullDataFile=NULL, 
                             removeCC=FALSE, 
                             species="mmu", 
                             rhop=rhop, 
                             write=TRUE, 
                             PMF=TRUE, 
                             useFull=FALSE, 
                             heatmap=heatmap,
                             centroids=TRUE,
                             num_doubs=100, 
                             only50=FALSE,
                             min_uniq=4,
                             nCores=-1))
  dev.off()
  return(results)
}
  


