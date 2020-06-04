library('tximport')
library('DESeq2')
library("pheatmap")
library("vsn")
library("RColorBrewer")
library("ggplot2")
library("reshape2")
library("biomaRt")
library("gProfileR")
library("plyr")
library("gridExtra")
library("EnhancedVolcano")
library("dplyr")
m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")

###### things to edit ######
#paths and metadata
maindir<-'/Users/fr7/git_repos/single_cell/experiment_4/bulk_libraries'
metadata<-read.table(file.path(maindir,'metadata.txt'),header=TRUE,stringsAsFactors = TRUE)
tx2gene <- read.table('~/git_repos/Trichuris_transwells/transcripts2genes.E98.tsv', header = T)

#remove samples with missing data
metadata<-metadata[which( !(metadata$sample %in% c("4672STDY8348322","4672STDY8348329","4672STDY8348338")) ),]

#add a new column for time-status
metadata$time.status<-paste0(metadata$time,'_',metadata$infection_status)
#metadata$time.status<-factor(metadata$time.status, levels=c("24h_control", "24h_infected", "72h_control", "72h_infected"))
metadata<-metadata[order(metadata$time.status),]

###### QC: import the table of read numbers and counts and plot ######

read_counts <- read.table(file.path(maindir,'kallisto_counts','counts_table.tsv'),header=FALSE)
names(read_counts) <- c("sample","total_reads","kallisto_counts")
read_counts<-melt(read_counts,id.vars = "sample")
p<-ggplot(data=read_counts, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") + 
  theme_minimal() +
  xlab('Sample') + 
  ylab('Count') +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=20), axis.text.x=element_text(angle=90))

pdf(file.path(maindir,'read_counts.pdf'),30,10)
print(p)
dev.off()

###### QC: PCA plot of all samples ######

#find the correct samples
sorting<-c('bulk','sort')
time.status<-c('24h_control','24h_infected','72h_control','72h_infected')
design_formula<-"~sorting+time.status"

subdir<-paste0('time.status.',paste(time.status,collapse ="."),'.',paste(sorting,collapse ="."))
dir.create(file.path(maindir, subdir))
dir<-file.path(maindir,subdir)
dir.all<-dir
metadata.all<-metadata[which(metadata$sorting %in% sorting & metadata$time.status %in% time.status),]

#write a table so we know which samples were used in the analysis
write.table(metadata.all,file.path(dir,'metadata.txt'),sep="\t")

#setup dds object
dds.all<-m$setup_dds_from_kallisto(metadata.all,maindir,as.formula(design_formula),tx2gene)

#plot PCAs
for (i in design){
  pca<-m$plot_pca(dds.all,i,dir)
}

###### compare 24h samples ######

#find the correct samples
sorting<-c('bulk','sort')
time.status<-c('24h_control','24h_infected')
design_formula<-"~sorting+time.status"

subdir<-paste0('time.status.',paste(time.status,collapse ="."),'.',paste(sorting,collapse ="."))
dir.create(file.path(maindir, subdir))
dir<-file.path(maindir,subdir)
dir.24<-dir
metadata.24<-metadata[which(metadata$sorting %in% sorting & metadata$time.status %in% time.status),]

#write a table so we know which samples were used in the analysis
write.table(metadata.24,file.path(dir,'metadata.txt'),sep="\t")

#setup dds object
dds.24<-m$setup_dds_from_kallisto(metadata.24,maindir,as.formula(design_formula),tx2gene)

#relevel factors
dds.24$time.status <- relevel(dds.24$time.status, ref = "24h_control")

#plot PCAs
for (i in design){
  pca<-m$plot_pca(dds.24,i,dir)
}

#get results - Wald test
dds.24<-DESeq(dds.24)
res.24 <- results(dds.24)
res.24<-res.24[order(res.24$padj),]

#biomart
res.24<-m$mart(res.24,0.05)
res.24<-res.24[order(res.24$padj),]

#write and plot results
res.24<-m$write_results(res.24,dir,"Wald")
volcano.24<-m$volcano_plot(res.24,dir,"volcano")

###### compare 72h samples ######

#find the correct samples
sorting<-c('bulk','sort')
time.status<-c('72h_control','72h_infected')
design_formula<-"~sorting+time.status"

subdir<-paste0('time.status.',paste(time.status,collapse ="."),'.',paste(sorting,collapse ="."))
dir.create(file.path(maindir, subdir))
dir<-file.path(maindir,subdir)
dir.72<-dir
metadata.72<-metadata[which(metadata$sorting %in% sorting & metadata$time.status %in% time.status),]

#write a table so we know which samples were used in the analysis
write.table(metadata.72,file.path(dir,'metadata.txt'),sep="\t")

#setup dds object
dds.72<-m$setup_dds_from_kallisto(metadata.72,maindir,as.formula(design_formula),tx2gene)

#relevel factors
dds.72$time.status <- relevel(dds.72$time.status, ref = "72h_control")

#plot PCAs
for (i in design){
  pca<-m$plot_pca(dds.72,i,dir)
}

#get results - Wald test
dds.72<-DESeq(dds.72)
res.72 <- results(dds.72)
res.72<-res.72[order(res.72$padj),]

#biomart
res.72<-m$mart(res.72,0.05)
res.72<-res.72[order(res.72$padj),]

#write and plot results
res.72<-m$write_results(res.72,dir,"Wald")
volcano.72<-m$volcano_plot(res.72,dir,"volcano")

###### heatmaps across all timepoints ######

#retrieve files of interesting genes - extract genes with abs(log2FC) > 1
genes.24<-scan(file = file.path(dir.24,'immune.txt'), sep = "\n", what = character())
genes.24<-res.24[res.24$external_gene_name %in% genes.24, c("ensembl_gene_id","external_gene_name","log2FoldChange" )]
genes.24<-genes.24[which(abs(genes.24$log2FoldChange) > 1), ]

genes.72<-scan(file = file.path(dir.72,'immune.txt'), sep = "\n", what = character())
genes.72<-res.72[res.72$external_gene_name %in% genes.72, c("ensembl_gene_id","external_gene_name","log2FoldChange" )]
genes.72<-genes.72[which(abs(genes.72$log2FoldChange) > 1), ]

genes<-rbind(genes.24,genes.72)

#extract relevant fpkms
fpkms<-fpkm(dds.all)[genes$ensembl_gene_id,]

rownames(fpkms)<-genes$external_gene_name
#colnames(fpkms)<-metadata.all$time.status
annotation<-data.frame(metadata.all$time.status)
row.names(annotation)<-metadata.all$sample

dir<-dir.all

#heatmap of replicates
pdf(file.path(dir,"heatmap_replicates.pdf"))
pheatmap(log((fpkms + 1), base=2),cluster_cols=F,treeheight_row = 0, treeheight_col = 0 , annotation_col = annotation , show_colnames = F)
dev.off()

#heatmap of means
control.24h <- as.character(metadata.all[which(metadata.all$time.status == "24h_control"), "sample" ])
infected.24h <- as.character(metadata.all[which(metadata.all$time.status == "24h_infected"), "sample" ])
control.72h <- as.character(metadata.all[which(metadata.all$time.status == "72h_control"), "sample" ])
infected.72h <- as.character(metadata.all[which(metadata.all$time.status == "72h_infected"), "sample" ])

#new dataframe with means
fpkm.means<-data.frame("control.24h" = rowMeans(fpkms[,control.24h]),
                       "infected.24h" = rowMeans(fpkms[,infected.24h]),
                       "control.72h" = rowMeans(fpkms[,control.72h]),
                       "infected.72h" = rowMeans(fpkms[,infected.72h])
                       )

pdf(file.path(dir,"heatmap_means.pdf"))
pheatmap(log((fpkm.means + 1), base=2),cluster_cols=F,treeheight_row = 0, treeheight_col = 0 ,cellwidth = 20, cellheight = 10)
dev.off()

#just plot the 24h genes
dir<-dir.24

fpkms.24<-fpkms[genes.24$external_gene_name,c(control.24h, infected.24h) ]
fpkms.24.means<-fpkm.means[genes.24$external_gene_name,c("control.24h","infected.24h")]

#24h replicates
pdf(file.path(dir,"heatmap_replicates.pdf"))
pheatmap(log((fpkms.24 + 1), base=2),cluster_cols=F,treeheight_row = 0, treeheight_col = 0 , annotation_col = annotation , show_colnames = F, cellheight= 10, annotation_names_col = T, annotation_legend = F)
dev.off()

#24h means
pdf(file.path(dir,"heatmap_means.pdf"))
pheatmap(log((fpkms.24.means + 1), base=2),cluster_cols=F,treeheight_row = 0, treeheight_col = 0 , cellwidth = 20, cellheight = 10)
dev.off()

#just plot the 72h genes
dir<-dir.72

fpkms.72<-fpkms[genes.72$external_gene_name,c(control.72h, infected.72h) ]
fpkms.72.means<-fpkm.means[genes.72$external_gene_name,c("control.72h","infected.72h")]

#72h replicates
pdf(file.path(dir,"heatmap_replicates.pdf"))
pheatmap(log((fpkms.72 + 1), base=2),cluster_cols=F,treeheight_row = 0, treeheight_col = 0 , annotation_col = annotation , show_colnames = F, cellheight= 10, annotation_names_col = T, annotation_legend = F)
dev.off()

#72h means
pdf(file.path(dir,"heatmap_means.pdf"))
pheatmap(log((fpkms.72.means + 1), base=2),cluster_cols=F,treeheight_row = 0, treeheight_col = 0 , cellwidth = 20, cellheight = 10)
dev.off()



