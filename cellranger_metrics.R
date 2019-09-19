library(ggplot2)

setwd('/Users/fr7/git_repos/single_cell')
CR <- read.table('cellranger_metrics.tsv', header = T)
meta_data <- read.table('samples.txt',header=T)

#join cellranger table with meta data
CR <- merge(CR,meta_data,no.dups=TRUE)
#make a new column for mouse_id-sample_id (for plotting)
CR$mouse_id.sample_id <- paste0(CR$mouse_id,'-',CR$sample_id)

write.table(CR,'cellranger_metadata.tsv',sep="\t")

#select only samples that have been run with versions 2 AND 3.
CR.2.3 <- CR[0,]
for (i in unique(CR$sample_id)){
  versions <- CR[which(CR$sample_id == i), ]
  if ( '2.1.1' %in% versions$cellranger_versiom & '3.0.2' %in% versions$cellranger_versiom  ){
    CR.2.3 <- rbind(CR.2.3, CR[which(CR$sample_id == i & CR$cellranger_versiom == '3.0.2'), ] ,  CR[which(CR$sample_id == i & CR$cellranger_versiom == '2.1.1'), ] )
  }
}

#plotting function
get_dotplot<-function(data,x,y,fill,metric){
  g<-ggplot(data,aes(x,y, colour=fill)) +
  geom_jitter(width = 0.1, size = 3) +
  labs(title= metric, x= deparse(substitute(x)) , y = metric)
  return(g)
}

metrics<-c("Number_of_Reads","Sequencing_Saturation" ,"Reads_Mapped_to_Genome", "Estimated_Number_of_Cells" , "Fraction_Reads_in_Cells" , "Mean_Reads_per_Cell" , "Median_Genes_per_Cell" , "Total_Genes_Detected" , "Median_UMI_Counts_per_Cell" )

for (metric in metrics){
  y<-eval(parse(text= paste0("CR.2.3$",metric)))
  plot <- get_dotplot(CR.2.3,CR.2.3$cellranger_versiom,y,CR.2.3$mouse_id.sample_id,metric)
  assign(paste0(metric), plot)
# print(plot)
  pdf(paste0("QC_plots/",metric,".pdf") )
  print(plot)
  dev.off()
}




