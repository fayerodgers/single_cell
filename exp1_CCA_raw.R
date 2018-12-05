#For a look at the raw data (before we do any filtering)
library(Seurat)
library(Matrix)
library(dplyr)
library(metap)

#Define the output directory where figures will be put
out<-"/Users/fr7/Documents/single_cell/Experiment_1/pngs/"

#Define a file name for the R data structure
data_file<-"/Users/fr7/Documents/single_cell/Experiment_1/experiment_1_raw.Rdata"

#Read the 10X data directories and set up Seurat objects
#Control
control.data <- Read10X(data.dir = "/Users/fr7/Documents/single_cell/Experiment_1/22414_7/")
control <- CreateSeuratObject(raw.data=control.data, project = "CONTROL")
control@meta.data$infection_status <- "CONTROL"

#Infected
infected.data <- Read10X(data.dir = "/Users/fr7/Documents/single_cell/Experiment_1/22414_8/")
infected <- CreateSeuratObject(raw.data=infected.data, project = "INFECTED")
infected@meta.data$infection_status <- "INFECTED"

#A function to calculate % mitochondrial reads, normalise data, scale data and find variable genes.
normalise_and_scale<-function(sample){
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = sample@data), value = TRUE)
  percent.mito <- Matrix::colSums(sample@raw.data[mito.genes, ])/Matrix::colSums(sample@raw.data)
  sample <- AddMetaData(sample, metadata = percent.mito, col.name = "percent.mito")
  sample <- NormalizeData(sample)
  sample <- ScaleData(sample)
  sample <- FindVariableGenes(sample, do.plot=F)
  return(sample)
}

samples <- list(control,infected)
samples<-lapply(samples, FUN = normalise_and_scale)

#take variable genes (dispersion=variance/mean) for each data set. 
for (sample in samples){
  genes.use <- unique(c(genes.use,sample@var.genes))
  genes.use <- intersect(genes.use, rownames(sample@scale.data))
}

#Run CCA
experiment<- RunCCA(samples[[1]],samples[[2]], genes.use = genes.use, num.cc = 30, add.cell.id1 ='control', add.cell.id2 = 'infected')
pdf(paste0(out,"cc_saturation_raw.pdf"))
p1 <- MetageneBicorPlot(experiment, grouping.var = "infection_status", dims.eval = 1:30, display.progress = FALSE) #use this plot to decide how many CCs to include downstream. Choose 15.
p1
dev.off()

#STOP HERE! Decide how many CCs to include in the analysis.
################

experiment <- AlignSubspace(experiment, reduction.type = "cca", grouping.var = "infection_status", dims.align = 1:15)
experiment <- FindClusters(experiment, reduction.type = "cca.aligned", resolution = 1, dims.use = 1:15)
experiment <- RunTSNE(experiment, reduction.use = "cca.aligned", dims.use = 1:15, do.fast = T)
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

new.ident<-c("ENTEROCYTE_A","TA_A","ENTEROCYTE_B","GOBLET","ENTEROCYTE_C","ENTEROCYTE_D","TA_B")
for (i in 0:6){ experiment<-RenameIdent(experiment,old.ident.name = i,new.ident.name = new.ident[i+1])}
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
cell_types<-c("ENTEROCYTE_A","ENTEROCYTE_B","ENTEROCYTE_C","ENTEROCYTE_D","TA_A","TA_B","GOBLET")
control<-c()
infected<-c()
control_percent<-c()
infected_percent<-c()
control_sum = nrow(experiment@meta.data[which(experiment@meta.data$infection_status=='CONTROL'),])
infected_sum = nrow(experiment@meta.data[which(experiment@meta.data$infection_status=='INFECTED'),])
for (cell_type in cell_types){
  x = nrow(experimentw@meta.data[which(experiment@meta.data$cell_type==cell_type & experiment@meta.data$infection_status=='CONTROL'),])
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
