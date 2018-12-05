#For a look at the raw data (before we do any filtering).
#This script contains the parameters used for experiment 2.
library(Seurat)
library(Matrix)
library(dplyr)
library(metap)

#Define the output directory where figures will be put
out<-"/Users/fr7/Documents/single_cell/Experiment_2/pngs/"

#Define a file name for the R data structure
data_file<-"/Users/fr7/Documents/single_cell/Experiment_2/experiment_2_raw.Rdata"

#Read the 10X data directories and set up Seurat objects
control_files <- c("24134_4672STDY7140017","24134_4672STDY7140018")
infected_files <- c("24134_4672STDY7140019","24134_4672STDY7140020") 

samples <- c()
i=1
for (file in control_files){
  data <- Read10X(data.dir = paste0("/Users/fr7/Documents/single_cell/Experiment_2/",file,"/"))
  colnames(data)=paste0("m",i,"_",colnames(data)) #need to add a unique tag to each run-cell names might not be unique between runs
  mouse <- CreateSeuratObject(raw.data=data, project = paste0("m",i) )
  mouse@meta.data$infection_status <- "CONTROL"
  samples <- c(samples, mouse)
  i = i+1
}

for (file in infected_files){
  data <- Read10X(data.dir = paste0("/Users/fr7/Documents/single_cell/Experiment_2/",file,"/"))
  colnames(data)=paste0("m",i,"_",colnames(data))
  mouse <- CreateSeuratObject(raw.data=data, project = paste0("m",i) )
  mouse@meta.data$infection_status <- "INFECTED"
  samples <- c(samples, mouse)
  i = i+1
}

#A function to calculate percentage of mitochondrial reads; normalise and scale data; find variable genes
normalise_and_scale<-function(sample){
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = sample@data), value = TRUE)
  percent.mito <- Matrix::colSums(sample@raw.data[mito.genes, ])/Matrix::colSums(sample@raw.data)
  sample <- AddMetaData(sample, metadata = percent.mito, col.name = "percent.mito")
  sample <- NormalizeData(sample)
  sample <- ScaleData(sample)
  sample <- FindVariableGenes(sample, do.plot=F)
  return(sample)
}

samples<-lapply(samples, FUN = normalise_and_scale)

#take variable genes (dispersion=variance/mean) for each data set.
genes.use <- c()
for (sample in samples){
  genes.use <- unique(c(genes.use,sample@var.genes))
  genes.use <- intersect(genes.use, rownames(sample@scale.data))
}

#Run CCA
experiment<- RunMultiCCA(samples, genes.use = genes.use, num.ccs = 40)
pdf(paste0(out,"cc_saturation_raw.pdf"))
p1 <- MetageneBicorPlot(experiment, grouping.var = "infection_status", dims.eval = 1:40, display.progress = FALSE) #use this plot to decide how many CCs to include downstream. Choose 30.
dev.off()

#STOP HERE! Decide how many CCs to include in the analysis.
################

experiment <- AlignSubspace(experiment, reduction.type = "cca", grouping.var = "infection_status", dims.align = 1:30)
experiment <- FindClusters(experiment, reduction.type = "cca.aligned", resolution = 1, dims.use = 1:30) #Might want to play with the resolution.
experiment <- RunTSNE(experiment, reduction.use = "cca.aligned", dims.use = 1:30, do.fast = T)
pdf(paste0(out,'tsne_infection_status_raw.pdf'))
p2 <- TSNEPlot(experiment, do.return = T, pt.size = 1.5, group.by = "infection_status")
p2
dev.off()
pdf(paste0(out,'tsne_clusters_raw.pdf'))
p3 <- TSNEPlot(experiment, do.label = T, do.return = T, pt.size = 1.5)
p3
dev.off()
pdf(paste0(out,'markers_raw.pdf'))
p4 <- FeaturePlot(experiment, features.plot = c("Alpi","Muc2","Chga","Lgr5","Dclk1","Cdk4","Il33","Prdx6","Rps26"), min.cutoff="q9", cols.use = c("lightgrey", "blue"), pt.size = 1.5) 
p4
dev.off()

####################
#Stop here! Examine the clusters to identify cell types
#rename cells according to identified clusters

new.ident<-c("ENTEROCYTE_A","ENTEROCYTE_B","ENTEROCYTE_C","TA_A","ENTEROCYTE_D","GOBLET","TA_B","ENTEROCYTE_E","ENTEROCYTE_F","TUFT")#This will change depending on the clusters identified.
for (i in 0:9){ experiment<-RenameIdent(experiment,old.ident.name = i,new.ident.name = new.ident[i+1])}
#Stash the new identities as 'cell_type'
experiment<-StashIdent(experiment,save.name="cell_type")
pdf(paste0(out,'named_tsne_clusters_raw.pdf'))
p5<- TSNEPlot(experiment, do.label = T, do.return = T, pt.size = 1.5)
p5
dev.off()

#Plots of nGene and percent.mito
pdf(paste0(out,'ngene_mito_by_mouse_raw.pdf'))
p6<-VlnPlot(experiment,features.plot=c('nGene','percent.mito'),group.by='orig.ident',point.size.use = 0.5, cols.use = BlackAndWhite())
p6
dev.off()

pdf(paste0(out,'ngene_mito_by_cluster_raw.pdf'))
p6<-VlnPlot(experiment,features.plot=c('nGene','percent.mito'),point.size.use = 0.5,do.sort=T,x.lab.rot = T, cols.use = BlackAndWhite())
p6
dev.off()

#Get cell numbers for each cluster
cell_types<-c("ENTEROCYTE_A","ENTEROCYTE_B","ENTEROCYTE_C","ENTEROCYTE_D","ENTEROCYTE_E","ENTEROCYTE_F","TA_A","TA_B","GOBLET","TUFT") 
control<-c()
infected<-c()
control_percent<-c()
infected_percent<-c()
control_sum = nrow(experiment@meta.data[which(experiment@meta.data$infection_status=='CONTROL'),])
infected_sum = nrow(experiment@meta.data[which(experiment@meta.data$infection_status=='INFECTED'),])
for (cell_type in cell_types){
  x = nrow(experiment@meta.data[which(experiment@meta.data$cell_type==cell_type & experiment@meta.data$infection_status=='CONTROL'),])
  x_percent = round((x*100/control_sum),1)
  y = nrow(experiment@meta.data[which(experiment@meta.data$cell_type==cell_type & experiment@meta.data$infection_status=='INFECTED'),])
  y_percent = round((y*100/infected_sum),1)
  control<-c(control,x)
  infected<-c(infected,y)
  control_percent<-c(control_percent,x_percent)
  infected_percent<-c(infected_percent,y_percent)
}
control<-c(control,control_sum)
infected<-c(infected,infected_sum)
control_percent<-c(control_percent,round((control_sum*100/control_sum),1))
infected_percent<-c(infected_percent,round((infected_sum*100/infected_sum),1))

df<-data.frame(control,infected,control_percent,infected_percent)
colnames(df)=c("CONTROL(n)","INFECTED(n)","CONTROL(%)","INFECTED(%)")
rownames(df)=c(cell_types,"SUM")

#save the data structure
saveRDS(experiment,data_file)




