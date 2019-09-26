#use_python("/Users/fr7/miniconda2/bin/python",required=TRUE)

import("Seurat")
#import("reticulate")
#import("Matrix")
#import("dplyr")
#import("metap")
#import("ggplot2")
#import("cowplot")


normalize_data<-function(sample,cellranger_version,meta_data,variable_features,ncells,nmito,nfeatures){
  data.dir = paste0("./count_matrices/",cellranger_version,"/",sample)
  print(data.dir)
  data <- Read10X(data.dir = paste0("./count_matrices/",cellranger_version,"/",sample),unique.features = TRUE)
  seurat_object <- CreateSeuratObject(counts=data, project = paste0(sample),min.cells=ncells,min.features=nfeatures)
  seurat_object[["infection_status"]] <- as.character(meta_data[which(meta_data$sample_id == sample),'treatment']) 
  seurat_object[["experiment"]] <- as.character(meta_data[which(meta_data$sample_id == sample),'experiment'])
  seurat_object[["cellranger_version"]]<-cellranger_version
  seurat_object[["percent.mito"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object <- SubsetData(seurat_object, subset.name = "percent.mito", high.threshold = nmito)
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst",nfeatures = variable_features)
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
  all<-CombinePlots(plots=list(p1,p2,p3,p4),ncol=2)
  return(all)
}

find_markers<-function(seurat_object,cluster){
  markers<-FindMarkers(seurat_object, ident.1 = cluster, print.bar = FALSE)
  file<-paste0("cluster.",cluster,".markers.tsv")
  write.table(markers, file = file, sep ='\t', eol = '\n') 
  return(markers)
}

find_regulated_genes<-function(seurat_object,cluster){
  x<-subset(seurat_object, idents =  cluster) 
  Idents(x)<-"infection_status"
  infection.response<-FindMarkers(x,ident.1="control",ident.2="infected",print.bar=FALSE)
  file<-paste0("cluster.",cluster,".regulated.tsv")
  write.table(infection.response, file = file, sep ='\t', eol = '\n') 
  return(infection.response)
}