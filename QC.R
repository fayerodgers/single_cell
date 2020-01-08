library(cowplot) #for plot_grid
m <- modules::use("SC.R")

evaluate<-function(text){
  x<-eval(parse(text=text))
  return(x)
}

setwd('/Users/fr7/git_repos/single_cell')
CR <- read.table('metadata/cellranger_metrics.txt', header = T)
meta_data <- read.table('metadata/samples.txt',header=T)
samples<-meta_data$sample_id  #analyse all samples

#join cellranger table with meta data
CR <- merge(CR,meta_data,no.dups=TRUE)
#make a new column for mouse_id-sample_id (for plotting)
CR$mouse_id.sample_id <- paste0(CR$mouse_id,'-',CR$sample_id)
#extract all CellRangerv3, mm10-3_0_0 samples
CR.3<-CR[which(CR$cellranger_versiom == '3.0.2' & CR$transcriptome == 'mm10-3_0_0'),]

#dot plots for interesting metrics
metrics<-c("Number_of_Reads","Sequencing_Saturation" ,"Reads_Mapped_to_Genome", "Estimated_Number_of_Cells" , "Fraction_Reads_in_Cells" , "Mean_Reads_per_Cell" , "Median_Genes_per_Cell" , "Total_Genes_Detected" , "Median_UMI_Counts_per_Cell" )
for (metric in metrics){
  y<-evaluate(paste0("CR.3$",metric))
  plot <- m$get_dotplot(CR.3,CR.3$experiment,y,CR.3$QC,metric)
  assign(paste0(metric), plot)
  pdf(paste0("QC_plots/all_data/",metric,".pdf") )
  print(plot)
  dev.off()
}

#make one figure with all plots
plots<-lapply(metrics,evaluate)
all_plots<-plot_grid(plotlist=plots,ncol=3)
pdf('QC_plots/all_data/all_QC_plots.pdf',12,8)
print(all_plots)
dev.off()

#generate seurat objects
seurat_objects<-lapply(samples,m$normalize_data,'cellranger302', meta_data, 2000,0,100,0,'mm10-3_0_0')

#violin plots of %mito
mito_violins<-lapply(seurat_objects,m$violin,"percent.mito","orig.ident",100)
all_mito<-plot_grid(plotlist=mito_violins,ncol=9,nrow=4)
pdf('QC_plots/all_data/all_mito.pdf',60,30)
print(all_mito)
dev.off()

#violin plots of nfeatures
nfeatures_violins<-lapply(seurat_objects,m$violin,"nFeature_RNA","orig.ident",8500)
all_nfeatures<-plot_grid(plotlist=nfeatures_violins,ncol=9,nrow=4)
pdf('QC_plots/all_data/all_nfeatures.pdf',60,30)
print(all_nfeatures)
dev.off()

#scatter plots
scatter_plots<-lapply(seurat_objects,m$scatter,"nFeature_RNA","percent.mito",0,8500,0,100)
all_scatters<-plot_grid(plotlist=scatter_plots,ncol=9,nrow=4)
pdf('QC_plots/all_data/all_scatters.pdf',60,30)
print(all_scatters)
dev.off()
