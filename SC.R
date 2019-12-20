import("Seurat")
import("ggplot2")
import("utils")
import("DESeq2")
import("grDevices")

evaluate<-function(text){
  x<-eval(parse(text=text))
  return(x)
}

normalize_data<-function(sample,cellranger_version,meta_data,variable_features,ncells,nmito,nfeatures,annotation){
  data <- Read10X(data.dir = paste0("./count_matrices/",annotation,"/",cellranger_version,"/",sample,"/filtered_feature_bc_matrix"),unique.features = TRUE)
  seurat_object <- CreateSeuratObject(counts=data, project = paste0(sample),min.cells=ncells,min.features=nfeatures)
  seurat_object[["infection_status"]] <- as.character(meta_data[which(meta_data$sample_id == sample),'treatment']) 
  seurat_object[["experiment"]] <- as.character(meta_data[which(meta_data$sample_id == sample),'experiment'])
  seurat_object[["cellranger_version"]]<-cellranger_version
  seurat_object[["percent.mito"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object[["time"]]<-as.character(meta_data[which(meta_data$sample_id == sample),'time']) 
  seurat_object <- subset(seurat_object, subset.name = "percent.mito", high.threshold = nmito)
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst",nfeatures = variable_features)
  return(seurat_object)
}

run_pca<-function(seurat_object,ndims){
  seurat_object <- ScaleData(seurat_object, verbose = FALSE)
  seurat_object <- RunPCA(seurat_object, npcs = ndims, verbose = FALSE)
  return(seurat_object)
}

integrate_datasets<-function(list_of_seurat_objects,ndims){        
  anchors<-FindIntegrationAnchors(list_of_seurat_objects,dims=1:ndims)
  combined<- IntegrateData(anchorset = anchors, dims = 1:ndims)
  DefaultAssay(combined) <- "integrated"
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = ndims, verbose = FALSE)
  return(combined)
}  

do_clustering<-function(combined,pcdimensions,resolution){
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:pcdimensions)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:pcdimensions)
  combined <- FindClusters(combined, resolution = resolution)
  return(combined)
}

plot_UMAPS<-function(seurat_object){
  p1<-DimPlot(seurat_object,reduction="umap")
  p2<-DimPlot(seurat_object,reduction="umap",group.by = "orig.ident")
  p3<-FeaturePlot(seurat_object,features="percent.mito")
  p4<-FeaturePlot(seurat_object,features = "nFeature_RNA")
  all<-CombinePlots(plots=list(p1,p2,p3,p4),ncol=3)
  return(all)
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
  DefaultAssay(seurat_object) <- "RNA"
  infection.response<-try(FindMarkers(x,ident.1="control",ident.2="infected",print.bar=FALSE,test.use="DESeq2"),silent=TRUE)
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

plot_proportions<-function(cluster,proportions_table,x){
    my_data<-proportions_table[which(proportions_table$cluster == cluster),]
    x<-my_data[,x]
    print(x)
    plot<-get_dotplot(my_data,x,my_data$proportion,NULL,paste0("cluster.",cluster))
    return(plot)
  }

set_colour_by_experiment<-function(seurat_object){
  if (seurat_object$experiment == '1'){col="red"}
  if (seurat_object$experiment == '2'){col="forestgreen"}
  if (seurat_object$experiment == '3.1' | seurat_object$experiment == '3.2'){col="cornflowerblue"}
  if (seurat_object$experiment == '4'){col="orange"}
  return(col)
}

violin<-function(seurat_object,feature,group,ymax){
  col<-set_colour_by_experiment(seurat_object)
  vln<-VlnPlot(seurat_object, features = c(feature), cols=c(col), pt.size = 0,group.by = group, do.return=T,y.max=ymax)
  vln<-vln + theme(text = element_text(size=30), axis.text.x=element_text(angle=0,hjust=0.5,size=25),axis.text.y=element_text(size=25),legend.position = "none") + labs(x="",y="",title=paste0(seurat_object$infection_status)) 
  return(vln)
}

scatter<-function(seurat_object,feature1,feature2,xmin,xmax,ymin,ymax){
  col<-set_colour_by_experiment(seurat_object)
  scatter_plot<-FeatureScatter(seurat_object,feature1=feature1,feature2 = feature2,cols=c(col))
  scatter_plot<-scatter_plot+ xlim(xmin,xmax) + ylim(ymin,ymax)+ theme(text = element_text(size=30), axis.text.x=element_text(angle=0,hjust=0.5,size=25),axis.text.y=element_text(size=25),legend.position = "none") + labs(title=paste0(seurat_object$orig.ident)) 
  return(scatter_plot)
}

find_all_markers<-function(seurat_object){
  markers<-FindAllMarkers(seurat_object,min.pct = 0.3,logfc.threshold = 0.75)
  write.table(markers, file = paste0(levels(seurat_object$orig.ident),"/markers.txt"), sep ='\t', eol = '\n') 
  return(markers)
}

QC_violins<-function(seurat_object){
  p<-VlnPlot(seurat_object,features=c("percent.mito","nFeature_RNA","nCount_RNA"),pt.size=0)
  pdf(paste0(levels(seurat_object$orig.ident),"/plots.pdf"),12,6)
  print(p)
  dev.off()
  return(p)
}
  
remove_mito_clusters<-function(seurat_object,to_remove_table){
  sample<-levels(seurat_object$orig.ident)
  to_remove<-to_remove_table[which(to_remove_table$sample == sample),'remove']
  seurat_object<-subset(seurat_object, idents = to_remove, invert = TRUE)
  return(seurat_object)
}
