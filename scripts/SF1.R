#!/usr/bin/env Rscript
m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(scales)
library(stringr)
library(goseq)
library("clusterProfiler")
library("pathview")
library("fgsea")
library("GO.db")
library("msigdbr")
library("reshape2")

### paths and metadata ###
maindir_whole<-'/Users/fr7/timecourse/kallisto_analysis'
metadata_whole<-read.table(file.path(maindir_whole,'metadata.txt'),header=T,stringsAsFactors = TRUE)
metadata_whole$time<-as.factor(metadata_whole$time)

maindir_epi<-'/Users/fr7/git_repos/single_cell/experiment_4/bulk_libraries'
metadata_epi<-read.table(file.path(maindir_epi,'metadata.txt'),header=TRUE,stringsAsFactors = TRUE)

figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'
files_dir<-'/Users/fr7/git_repos/single_cell/SupplementaryFiles'

tx2gene <- read.table('~/git_repos/Trichuris_transwells/transcripts2genes.E98.tsv', header = T)

### analysis - whole caecum samples D0,4,7 ###

# d0, d4 and d7 only, high dose infection and uninfected only. Remove the d4 and d7 controls to analyse as time series
metadata<-metadata_whole[which(((metadata_whole$time == '0' & metadata_whole$infection == 'N') |
                            (metadata_whole$time == '4' & metadata_whole$infection == 'HD') |
                            (metadata_whole$time == '7' & metadata_whole$infection == 'HD')) &
                           (metadata_whole$tissue == 'LIA' | metadata_whole$tissue == 'LIE' | metadata_whole$tissue == 'LIC') &
                           (metadata_whole$batch == '409850') 
),]

dir.time<-file.path(maindir_whole,'as_time_course')
dir.create(dir.time)
design<-c('tissue','time')
design_formula<-"~tissue+time"

# write a table so we know which samples were used in the analysis
write.table(metadata,file.path(dir.time,'metadata.txt'),sep="\t")

# relevel factors
metadata$time<-relevel(metadata$time, ref="0")

# setup dds.whole object
dds.whole<-m$setup_dds_from_kallisto(metadata,maindir_whole,as.formula(design_formula),tx2gene)

# plot PCAs
for (i in design){
  pca<-m$plot_pca(dds.whole,i,dir.time)
}

reduced_formula<-"~tissue"
dds.whole<-m$deseq2_lrt(dds.whole,as.formula(reduced_formula))

results_d4<-results(dds.whole,contrast=c("time","4","0"))
results_d4<-m$mart(results_d4,0.05)
results_d4<-results_d4[order(results_d4$padj),]
results_d4<-m$write_results(results_d4,dir.time,"LRT_d4")
volcano<-m$volcano_plot(results_d4,dir.time,'volcano_d4')

results_d7<-results(dds.whole,contrast=c("time","7","0"))
results_d7<-m$mart(results_d7,0.05)
results_d7<-results_d7[order(results_d7$padj),]
results_d7<-m$write_results(results_d7,dir.time,"LRT_d7")
volcano<-m$volcano_plot(results_d7,dir.time,"volcano_d7")

# write nice supplementary files
res.d4.sig<-results_d4[which(results_d4$padj < 0.05),]
text<-"Differential expression of genes in whole caecum across a 7 day infection time course (likelihood ratio test). Fold changes show changes at 4 days post infection (relative to day 0 uninfected controls)."
res.d4.text<-m$text_matrix(res.d4.sig,text)
#res.d4.text<-m$write_results_text(res.d4.text,files_dir,"SF.1.tsv")

res.d7.sig<-results_d7[which(results_d7$padj < 0.05),]
text<-"Differential expression of genes in whole caecum across a 7 day infection time course (likelihood ratio test). Fold changes show changes at 7 days post infection (relative to day 0 uninfected controls)."
res.d7.text<-m$text_matrix(res.d7.sig,text)
#res.d7.text<-m$write_results_text(res.d7.text,files_dir,"SF.2.tsv")

# merge the d4 and d7 results into one table
SFile1.1 <- merge(x = res.d4.sig, y = res.d7.sig, by = colnames(res.d4.sig)[! colnames(res.d4.sig) %in% c("log2FoldChange","lfcSE")])
SFile1.1 <- SFile1.1[order(SFile1.1$padj),]
colnames(SFile1.1) <- c("Ensembl gene ID", "Gene name", "Description", "Mean of normalized counts", "LRT statistic", "LRT p value", "Adjusted p value", "Log2Fold Change (D4 v D0)", "Standard error (D4 v D0)", "Log2Fold Change (D7 v D0)", "Standard error (D7 v D0)" )
text<-"Differential expression of genes in whole caecum across a 7 day infection time course (likelihood ratio test)."
SFile1.1.text<-m$text_matrix(SFile1.1,text)
SFile1.1.text<-m$write_results_text(SFile1.1.text,files_dir,"SFile1.1.tsv")

# GO enrichment

tx.lengths<-rowMedians(dds.whole@assays@data$avgTxLength)
names(tx.lengths)<-names(dds.whole)

# remove genes that are not expressed
results_d7.expressed<-results_d7[!(results_d7$baseMean == 0),]

sig.d7<-as.integer(results_d7.expressed$padj<0.05 & !is.na(results_d7.expressed$padj))
names(sig.d7)<-results_d7.expressed$ensembl_gene_id
tx.lengths<-tx.lengths[names(tx.lengths) %in% names(sig.d7) ]

pwf <- nullp(sig.d7, "mm10", "ensGene",bias.data=tx.lengths)

go <- goseq(pwf, "mm10","ensGene", test.cats=c("GO:BP"))
go$p.adjust <- p.adjust(go$over_represented_pvalue, method = "BH")
go <- go[which(go$p.adjust < 0.05),]
go <- go[order(go$p.adjust),]
go <- m$write_results(go,dir.time,"GO.BP")

# write nice supplementary file
colnames(go) <- c("GO ID", "Over represented p value", "Under represented p value", "Num DE genes in cat", "Num genes in cat", "GO term name", "Ontology", "BH-adjusted p value")
text <- "GO term (BP) enrichment in significantly regulated genes in whole caecum infection time course"
go.text<-m$text_matrix(go,text)
go.text<-m$write_results_text(go.text,files_dir,"SFile1.2.tsv")

### analysis - epithelium samples D1,3 ###

# remove samples with missing data
metadata_epi<-metadata_epi[which( !(metadata_epi$sample %in% c("4672STDY8348322","4672STDY8348329","4672STDY8348338")) ),]

# add a new column for time-status
metadata_epi$time.status<-as.factor(paste0(metadata_epi$time,'_',metadata_epi$infection_status))
metadata_epi<-metadata_epi[order(metadata_epi$time.status),]

###### QC: import the table of read numbers and counts and plot ######

read_counts <- read.table(file.path(maindir_epi,'kallisto_counts','counts_table.tsv'),header=FALSE)
names(read_counts) <- c("sample","total_reads","kallisto_counts")
read_counts<-melt(read_counts,id.vars = "sample")
p<-ggplot(data=read_counts, aes(x=sample,y=value,fill=variable)) +
  geom_bar(stat="identity",position = "dodge") + 
  theme_minimal() +
  xlab('Sample') + 
  ylab('Count') +
  theme(axis.title = element_text(size=20), axis.text  = element_text(size=20), axis.text.x=element_text(angle=90))

pdf(file.path(maindir_epi,'read_counts.pdf'),30,10)
print(p)
dev.off()

# QC: PCA plot of all samples

# find the correct samples
design<-c('sorting','time.status')
sorting<-c('bulk','sort')
time.status<-c('24h_control','24h_infected','72h_control','72h_infected')
design_formula<-"~sorting+time.status"

subdir<-paste0('time.status.',paste(time.status,collapse ="."),'.',paste(sorting,collapse ="."))
dir.create(file.path(maindir_epi, subdir))
dir.epi<-file.path(maindir_epi,subdir)
metadata_epi<-metadata_epi[which(metadata_epi$sorting %in% sorting & metadata_epi$time.status %in% time.status),]

# write a table so we know which samples were used in the analysis
write.table(metadata_epi,file.path(dir.epi,'metadata.txt'),sep="\t")

# setup dds object
dds.all.epi<-m$setup_dds_from_kallisto(metadata_epi,maindir_epi,as.formula(design_formula),tx2gene)

# plot PCAs
for (i in design){
  pca<-m$plot_pca(dds.all.epi,i,dir.epi)
}

# compare 24h samples 

# find the correct samples
sorting<-c('bulk','sort')
time.status<-c('24h_control','24h_infected')
design_formula<-"~sorting+time.status"

subdir<-paste0('time.status.',paste(time.status,collapse ="."),'.',paste(sorting,collapse ="."))
dir.create(file.path(maindir_epi, subdir))
dir<-file.path(maindir_epi,subdir)
dir.24<-dir
metadata.24<-metadata_epi[which(metadata_epi$sorting %in% sorting & metadata_epi$time.status %in% time.status),]

# write a table so we know which samples were used in the analysis
write.table(metadata.24,file.path(dir.24,'metadata.txt'),sep="\t")

# setup dds object
dds.24<-m$setup_dds_from_kallisto(metadata.24,maindir_epi,as.formula(design_formula),tx2gene)

# relevel factors
dds.24$time.status <- relevel(dds.24$time.status, ref = "24h_control")

# plot PCAs
for (i in design){
  pca<-m$plot_pca(dds.24,i,dir)
}

# get results - Wald test
dds.24 <- DESeq(dds.24)
res.24 <- results(dds.24)
res.24 <- res.24[order(res.24$padj),]

# biomart
res.24<-m$mart(res.24,0.05)
res.24<-res.24[order(res.24$padj),]

# write and plot results
res.24<-m$write_results(res.24,dir.24,"Wald")
volcano.24<-m$volcano_plot(res.24,dir.24,"volcano")

# write supplementary file
res.24.sig<-res.24[which(res.24$padj < 0.05),]
colnames(res.24.sig) <- c("Ensembl gene ID", "Gene name", "Description", "Mean of normalized counts", "Log2Fold Change", "Log2FoldChange Standard error", "Wald statistic", "Wald p value", "Adjusted p value")
text<-"Differential expression of genes in caecal IECs 1 day post Trichuris infection (Wald test; control (D1) v infected (D1))"
res.24.text<-m$text_matrix(res.24.sig,text)
res.24.text<-m$write_results_text(res.24.text,files_dir,"SFile2.1.tsv")

# GO - 24h 
# get transcript length data from dds object

tx.lengths<-rowMedians(dds.24@assays@data$avgTxLength)
names(tx.lengths)<-names(dds.24)

# remove genes that are not expressed
res.24.expressed<-res.24[!(res.24$baseMean == 0),]

sig.24<-as.integer(res.24.expressed$padj<0.05 & !is.na(res.24.expressed$padj))
names(sig.24)<-res.24.expressed$ensembl_gene_id
tx.lengths<-tx.lengths[names(tx.lengths) %in% names(sig.24) ]

pwf <- nullp(sig.24, "mm10", "ensGene",bias.data=tx.lengths)

go.24 <- goseq(pwf, "mm10","ensGene", test.cats=c("GO:BP"))
go.24$p.adjust <- p.adjust(go.24$over_represented_pvalue, method = "BH")
go.24<-m$write_results(go.24,dir.24,"GO.BP") # nothing significant after adjustment

# write nice supplementary file
text<-"GO term (BP) enrichment in significantly regulated genes in caecum epithelial cells 1 day post infection"
go.24.text<-m$text_matrix(go.24,text)
go.24.text<-m$write_results_text(go.24.text,files_dir,"SF.7.tsv")


# compare 72h samples

# find the correct samples
sorting<-c('bulk','sort')
time.status<-c('72h_control','72h_infected')
design_formula<-"~sorting+time.status"

subdir<-paste0('time.status.',paste(time.status,collapse ="."),'.',paste(sorting,collapse ="."))
dir.create(file.path(maindir_epi, subdir))
dir<-file.path(maindir_epi,subdir)
dir.72<-dir
metadata.72<-metadata_epi[which(metadata_epi$sorting %in% sorting & metadata_epi$time.status %in% time.status),]

# write a table so we know which samples were used in the analysis
write.table(metadata.72,file.path(dir.72,'metadata.txt'),sep="\t")

# setup dds object
dds.72<-m$setup_dds_from_kallisto(metadata.72,maindir_epi,as.formula(design_formula),tx2gene)

# relevel factors
dds.72$time.status <- relevel(dds.72$time.status, ref = "72h_control")

# plot PCAs
for (i in design){
  pca<-m$plot_pca(dds.72,i,dir)
}

# get results - Wald test
dds.72<-DESeq(dds.72)
res.72 <- results(dds.72)
res.72<-res.72[order(res.72$padj),]

# biomart
res.72<-m$mart(res.72,0.05)
res.72<-res.72[order(res.72$padj),]

# write and plot results
res.72<-m$write_results(res.72,dir,"Wald")
volcano.72<-m$volcano_plot(res.72,dir,"volcano")

# write supplementary file
res.72.sig<-res.72[which(res.72$padj < 0.05),]
colnames(res.72.sig) <- c("Ensembl gene ID", "Gene name", "Description", "Mean of normalized counts", "Log2Fold Change", "Log2FoldChange Standard error", "Wald statistic", "Wald p value", "Adjusted p value")
text<-"Differential expression of genes in caecal IECs 3 days post Trichuris infection (Wald test; control (D3) v infected (D3))"
res.72.text<-m$text_matrix(res.72.sig,text)
res.72.text<-m$write_results_text(res.72.text,files_dir,"SFile2.2.tsv")

# GO - 72h

tx.lengths<-rowMedians(dds.72@assays@data$avgTxLength)
names(tx.lengths)<-names(dds.72)

# remove genes that are not expressed
res.72.expressed<-res.72[!(res.72$baseMean == 0),]

sig.72<-as.integer(res.72.expressed$padj<0.05 & !is.na(res.72.expressed$padj))
names(sig.72)<-res.72.expressed$ensembl_gene_id
tx.lengths<-tx.lengths[names(tx.lengths) %in% names(sig.72) ]

pwf <- nullp(sig.72, "mm10", "ensGene",bias.data=tx.lengths)

go.72 <- goseq(pwf, "mm10","ensGene", test.cats=c("GO:BP"))
go.72$p.adjust <- p.adjust(go.72$over_represented_pvalue, method = "BH")
go.72<-m$write_results(go.72,dir.72,"GO.BP")

# write nice supplementary file
colnames(go.72) <- c("GO ID", "Over represented p value", "Under represented p value", "Num DE genes in cat", "Num genes in cat", "GO term name", "Ontology", "BH-adjusted p value")
text<-"GO term (BP) enrichment in significantly regulated genes in caecum IECs 3 days post T. muris infection"
go.72.text<-m$text_matrix(go.72,text)
go.72.text<-m$write_results_text(go.72.text,files_dir,"SFile2.3.tsv")


### GO plot - SF1A ###

# paths
d1.go<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/bulk_libraries/time.status.24h_control.24h_infected.bulk.sort/GO.BP_results.tsv',header=T, stringsAsFactors = FALSE)
d3.go<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/bulk_libraries/time.status.72h_control.72h_infected.bulk.sort/GO.BP_results.tsv',header=T, stringsAsFactors = FALSE)
whole.go<-read.table('/Users/fr7/timecourse/kallisto_analysis/as_time_course/GO.BP_results.tsv',header=T, stringsAsFactors = FALSE)

d1.go<-go.24
d3.go<-go.72
whole.go<-go

# add new column for adjusted p val
d1.go$p.adjust <- p.adjust(d1.go$over_represented_pvalue, method = "BH")
d3.go$p.adjust <- p.adjust(d3.go$over_represented_pvalue, method = "BH")
whole.go$p.adjust <- p.adjust(whole.go$over_represented_pvalue, method = "BH")

# define GO terms to plot

go_to_plot<-c(
  "GO:0002221",
# "GO:0042742", # defense response to bacteria
  "GO:0009615",
  "GO:0051707",
  "GO:0045087",
  "GO:0035456",
  "GO:0044419",
  "GO:0034340",
  "GO:0035455"
)

lookup_ngenes<-function(go,go.table){
  ngenes<-go.table[which(go.table$category == go),'numDEInCat']
  return(ngenes)
}

lookup_pval<-function(go,go.table){
  pval<-go.table[which(go.table$category == go),'p.adjust']
  return(pval)
}

lookup_description<-function(go,go.table){
  desc<-go.table[which(go.table$category == go),'term']
  return(desc)
}

# d1
go.df.d1 <- data.frame("pathway_id" = go_to_plot, stringsAsFactors = FALSE )
go.df.d1$time <- 'Epithelium, D1'
go.df.d1$ngenes <- sapply(go_to_plot,lookup_ngenes,d1.go)
go.df.d1$pval <- sapply(go_to_plot,lookup_pval,d1.go)
go.df.d1$description<-sapply(go_to_plot,lookup_description,d1.go)

# d3
go.df.d3 <- data.frame("pathway_id" = go_to_plot, stringsAsFactors = FALSE )
go.df.d3$time <- 'Epithelium, D3'
go.df.d3$ngenes <- sapply(go_to_plot,lookup_ngenes,d3.go)
go.df.d3$pval <- sapply(go_to_plot,lookup_pval,d3.go)
go.df.d3$description<-sapply(go_to_plot,lookup_description,d3.go)

# whole
go.df.whole <- data.frame("pathway_id" = go_to_plot, stringsAsFactors = FALSE )
go.df.whole$time <- 'Complete caecum'
go.df.whole$ngenes <- sapply(go_to_plot,lookup_ngenes,whole.go)
go.df.whole$pval <- sapply(go_to_plot,lookup_pval,whole.go)
go.df.whole$description<-sapply(go_to_plot,lookup_description,whole.go)

go.df <- rbind(go.df.d1,go.df.d3,go.df.whole)

# plot

cols<-brewer.pal(n=11,name="RdYlBu")

# remove rows with 0 genes
go.df <- go.df[! go.df$ngenes == 0, ]

SF1A<-ggplot(go.df,x=time,y=description)+
  geom_point(aes(x=time,y=description, size=ngenes, col=pval))+
  scale_color_gradientn(limits = c(0.0001,0.05),
                       colours=cols,
                       oob=squish, 
                       name = "Pathway p-value\n(adjusted)", 
                       breaks = c(0.0001,0.05),
                       labels = c("< 0.0001","> 0.05"),
                       guide = guide_colourbar())+
  scale_size_continuous(name="Number of\ngenes", range = c(0.1,2), breaks = c(5,100,200))+
  theme(text = element_text(size=8),
        axis.text.x=element_text(size=6, angle=45, hjust = 1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-5,-10,-10) )+
  labs(x=NULL,y=NULL)+
  guides(color = guide_colourbar(barwidth = 0.5, barheight = 2), 
         size = guide_legend(keywidth = 0.5, keyheight = 0.5))

pdf(file.path(figures_dir,'SF1A.pdf'),3,2)
print(SF1A)
dev.off()

### Heatmap whole caecum - SF1B ###

heatmap_genes<-('/Users/fr7/git_repos/single_cell/experiment_4/bulk_libraries/genes_to_plot.tsv')
genes<-scan(file = heatmap_genes, sep = "\n", what = character())
genes<-genes[! genes %in% c("Fabp6","Apoa1","Clca4b")]

genes<-res.d7.sig[which(res.d7.sig$external_gene_name %in% genes & abs(res.d7.sig$log2FoldChange)>1.5),]
fpkms<-fpkm(dds.whole)[genes$ensembl_gene_id,]
rownames(fpkms)<-genes$external_gene_name

annotation <- data.frame("Time" = metadata$time, "Sample" = metadata$sample,  stringsAsFactors = FALSE)
d7.samples <- as.character(annotation[which(annotation$Time == "7"),"Sample"]) # required later for ordering
annotation <- annotation[order(annotation$Time),]
rownames(annotation) <- annotation$Sample
annotation <- droplevels(annotation)
col_order <- annotation$Sample
fpkms <- fpkms[,col_order]
annotation <- annotation["Time"]
annotation$Time <- as.character(annotation$Time)

annotation$Time[annotation$Time == "0"] <- "D0"
annotation$Time[annotation$Time == "4"] <- "D4"
annotation$Time[annotation$Time == "7"] <- "D7"

cols.heat <- colorRampPalette(colors = c("#313695","#FFFFBF","#A50026"))(200)

# heatmap of replicates
# order by average expression at D7
d7.fpkms <- fpkms[,d7.samples]
d7.means <- sort(apply(d7.fpkms, 1, mean), decreasing = TRUE )
fpkms <- fpkms[names(d7.means),]
#fpkms <- fpkms[,col_order]

pdf(file.path(figures_dir,"SF1B.pdf"),15,7)
pheatmap(log((fpkms + 1), base=2),
         col = cols.heat, 
         cluster_cols=F,
         treeheight_row = 0, 
         treeheight_col = 0, 
         annotation_col = annotation, 
         show_colnames = F, 
         cluster_rows = T,
         cellwidth = 20, 
         cellheight = 10, 
      #  fontsize = 4, 
      #  breaks = 0:12,
      #  annotation_legend = FALSE, 
      #  legend = FALSE),
         legend_breaks = 0:12
      #  cex=1.0
)
dev.off()

### Heatmap whole caecum - SF4C ###

gene_names <- c("Il33","Il25","Tslp","Ifna1","Ifnb1","Ifng","Ifnl2","Ifnl3")
genes <- m$genename_to_ensembl(gene_names)
fpkms<-fpkm(dds.whole)[genes$ensembl_gene_id,]
rownames(fpkms)<-genes$external_gene_name

annotation <- data.frame("Time" = metadata$time, "Sample" = metadata$sample,  stringsAsFactors = FALSE)
d7.samples <- as.character(annotation[which(annotation$Time == "7"),"Sample"]) # required later for ordering
annotation <- annotation[order(annotation$Time),]
rownames(annotation) <- annotation$Sample
annotation <- droplevels(annotation)
col_order <- annotation$Sample
fpkms <- fpkms[,col_order]
annotation <- annotation["Time"]
annotation$Time <- as.character(annotation$Time)

annotation$Time[annotation$Time == "0"] <- "D0"
annotation$Time[annotation$Time == "4"] <- "D4"
annotation$Time[annotation$Time == "7"] <- "D7"

cols.heat <- colorRampPalette(colors = c("#313695","#FFFFBF","#A50026"))(200)

# heatmap of replicates
# order by average expression at D7
d7.fpkms <- fpkms[,d7.samples]
d7.means <- sort(apply(d7.fpkms, 1, mean), decreasing = TRUE )
fpkms <- fpkms[names(d7.means),]
#fpkms <- fpkms[,col_order]

pdf(file.path(figures_dir,"SF4C.pdf"),15,7)
pheatmap(log((fpkms + 1), base=2),
         col = cols.heat, 
         cluster_cols=F,
         treeheight_row = 0, 
         treeheight_col = 0, 
         annotation_col = annotation, 
         show_colnames = F, 
         cluster_rows = T,
         cellwidth = 20, 
         cellheight = 10, 
         #  fontsize = 4, 
         #  breaks = 0:12,
         #  annotation_legend = FALSE, 
         #  legend = FALSE),
         legend_breaks = 0:12
         #  cex=1.0
)
dev.off()

## alternative version of SF4C for the Wellcome grant response ##
gene_names <- c("Il33","Tslp","Il25","Ifna1","Ifnb1")
genes <- m$genename_to_ensembl(gene_names)
fpkms<-fpkm(dds.whole)[genes$ensembl_gene_id,]
rownames(fpkms)<-genes$external_gene_name

annotation <- data.frame("Time" = metadata$time, "Sample" = metadata$sample,  stringsAsFactors = FALSE)
d0.d7.samples <- as.character(annotation[which(annotation$Time == "7" | annotation$Time == "0" ),"Sample"]) 
annotation <- annotation[order(annotation$Time),]
rownames(annotation) <- annotation$Sample


# select D0 and D7 only
fpkms<-fpkms[,d0.d7.samples]
annotation<-annotation[which(annotation$Sample %in% d0.d7.samples),]
annotation <- droplevels(annotation)

col_order <- annotation$Sample
fpkms <- fpkms[,col_order]
annotation <- annotation["Time"]
annotation$Time <- as.character(annotation$Time)

annotation$Time[annotation$Time == "0"] <- "D0"
annotation$Time[annotation$Time == "7"] <- "D7"

cols.heat <- colorRampPalette(colors = c("#313695","#FFFFBF","#A50026"))(200)


fpkms <- fpkms[gene_names,]

pdf(file.path(figures_dir,"ALT.SF4C.pdf"),10,7)
pheatmap(log((fpkms + 1), base=2),
         col = cols.heat, 
         cluster_cols=F,
         treeheight_row = 0, 
         treeheight_col = 0, 
         annotation_col = annotation, 
         show_colnames = F, 
         cluster_rows = F,
         cellwidth = 20, 
         cellheight = 10, 
         #  fontsize = 4, 
         #  breaks = 0:12,
         #  annotation_legend = FALSE, 
         #  legend = FALSE),
         legend_breaks = 0:12
         #  cex=1.0
)
dev.off()

### Heatmap epithelium - SF1C ###

genes<-scan(file = heatmap_genes, sep = "\n", what = character())
genes.24<-res.24.sig[which(res.24.sig$external_gene_name %in% genes & abs(res.24.sig$log2FoldChange)>1.5),c("external_gene_name","ensembl_gene_id")]
genes.72<-res.72.sig[which(res.72.sig$external_gene_name %in% genes & abs(res.72.sig$log2FoldChange)>1.5),c("external_gene_name","ensembl_gene_id")]
genes.epi<-unique(rbind(genes.24,genes.72))
genes.epi<-genes.epi[! genes.epi$external_gene_name %in% c("Reg3g","Reg3b","Apoa1","Anpep","Fabp6","Olfm4","Cxcl13","Clca4b"),]

fpkms.epi<-fpkm(dds.all.epi)[genes.epi$ensembl_gene_id,]
rownames(fpkms.epi)<-genes.epi$external_gene_name

annotation<-data.frame("Condition"= metadata_epi$time.status, stringsAsFactors = FALSE)
row.names(annotation)<-metadata_epi$sample
annotation$Condition <- as.character(annotation$Condition)
annotation$Condition[annotation$Condition == "24h_control"] <- "D1 Control"
annotation$Condition[annotation$Condition == "24h_infected"] <- "D1 Infected"
annotation$Condition[annotation$Condition == "72h_control"] <- "D3 Control"
annotation$Condition[annotation$Condition == "72h_infected"] <- "D3 Infected"

d3.infected.samples <- as.character(metadata_epi[which(metadata_epi$time.status == "72h_infected"),"sample"])
d3.infected.fpkms <- fpkms.epi[,d3.infected.samples]
d3.infected.means <- sort(apply(d3.infected.fpkms, 1, mean), decreasing = TRUE)
fpkms.epi<-fpkms.epi[names(d3.infected.means), ]

# heatmap of replicates
pdf(file.path(figures_dir,"SF1C.orderbyd3i.pdf"),15,7)
pheatmap(log((fpkms.epi + 1), base=2),
         col = cols.heat, 
         cluster_cols=F,
         treeheight_row = 0, 
         treeheight_col = 0, 
         annotation_col = annotation, 
         show_colnames = F, 
         cluster_rows = F,
         cellwidth = 20, 
         cellheight = 10, 
         #  fontsize = 4, 
         #  breaks = 0:12,
         #  annotation_legend = FALSE, 
         #  legend = FALSE),
          legend_breaks = 0:9,
         #  cex=1.0
)
dev.off()

## Alternative version of SF1C for Wellcome grant response ##

genes<-c("Isg15", "Ifit1", "Ifit1bl1", "Ang4", "Ddx60", "Ifit3", "Ifi27l2a", "Oasl2", "Usp18", "Irf7")
genes.72 <- m$genename_to_ensembl(genes)

fpkms.72<-fpkm(dds.72)[genes.72$ensembl_gene_id,]
rownames(fpkms.72)<-genes.72$external_gene_name

annotation<-data.frame("Condition"= metadata.72$time.status, stringsAsFactors = FALSE)
row.names(annotation)<-metadata.72$sample
annotation$Condition <- as.character(annotation$Condition)
annotation$Condition[annotation$Condition == "72h_control"] <- "D3 Control"
annotation$Condition[annotation$Condition == "72h_infected"] <- "D3 Infected"

d3.infected.samples <- as.character(metadata.72[which(metadata.72$time.status == "72h_infected"),"sample"])
d3.infected.fpkms <- fpkms.72[,d3.infected.samples]
d3.infected.means <- sort(apply(d3.infected.fpkms, 1, mean), decreasing = TRUE)
fpkms.72<-fpkms.72[names(d3.infected.means), ]

# heatmap of replicates
pdf(file.path(figures_dir,"ALT.SF1C.pdf"),8,7)
pheatmap(log((fpkms.72 + 1), base=2),
         col = cols.heat, 
         cluster_cols=F,
         treeheight_row = 0, 
         treeheight_col = 0, 
         annotation_col = annotation, 
         show_colnames = F, 
         cluster_rows = F,
         cellwidth = 20, 
         cellheight = 10, 
         #  fontsize = 4, 
         #  breaks = 0:12,
         #  annotation_legend = FALSE, 
         #  legend = FALSE),
         legend_breaks = 0:9,
         #  cex=1.0
)
dev.off()


### Heatmap epithelium - SF4D ###

gene_names <- c("Il33","Il25","Tslp","Ifna1","Ifnb1","Ifng","Ifnl2","Ifnl3")
genes.epi <- m$genename_to_ensembl(gene_names)
fpkms<-fpkm(dds.whole)[genes.epi$ensembl_gene_id,]
rownames(fpkms)<-genes.epi$external_gene_name

fpkms.epi<-fpkm(dds.all.epi)[genes.epi$ensembl_gene_id,]
rownames(fpkms.epi)<-genes.epi$external_gene_name

annotation<-data.frame("Condition"= metadata_epi$time.status, stringsAsFactors = FALSE)
row.names(annotation)<-metadata_epi$sample
annotation$Condition <- as.character(annotation$Condition)
annotation$Condition[annotation$Condition == "24h_control"] <- "D1 Control"
annotation$Condition[annotation$Condition == "24h_infected"] <- "D1 Infected"
annotation$Condition[annotation$Condition == "72h_control"] <- "D3 Control"
annotation$Condition[annotation$Condition == "72h_infected"] <- "D3 Infected"

d3.infected.samples <- as.character(metadata_epi[which(metadata_epi$time.status == "72h_infected"),"sample"])
d3.infected.fpkms <- fpkms.epi[,d3.infected.samples]
d3.infected.means <- sort(apply(d3.infected.fpkms, 1, mean), decreasing = TRUE)
fpkms.epi<-fpkms.epi[names(d3.infected.means), ]

# heatmap of replicates
pdf(file.path(figures_dir,"SF4D.pdf"),15,7)
pheatmap(log((fpkms.epi + 1), base=2),
         col = cols.heat, 
         cluster_cols=F,
         treeheight_row = 0, 
         treeheight_col = 0, 
         annotation_col = annotation, 
         show_colnames = F, 
         cluster_rows = F,
         cellwidth = 20, 
         cellheight = 10, 
         #  fontsize = 4, 
         #  breaks = 0:12,
         #  annotation_legend = FALSE, 
         #  legend = FALSE),
         legend_breaks = 0:9,
         #  cex=1.0
)
dev.off()


### pathway enrichment - F1A ###

# IECs D1 
# rank all genes by logFC
ensembl_entrez_all<-m$ensembl_to_entrez(res.24$ensembl_gene_id)
entrez.res.df.all<-merge(ensembl_entrez_all,res.24, by='ensembl_gene_id',all=TRUE)
# limit to genes that are expressed
entrez.res.df.all<-entrez.res.df.all[!(entrez.res.df.all$baseMean == 0),]
# remove NAs
entrez.res.df.all <- entrez.res.df.all[!is.na(entrez.res.df.all$entrezgene_id),]
entrez.res.df.all <- entrez.res.df.all[!is.na(entrez.res.df.all$pvalue),]
# new column for -log10(pvalue)
entrez.res.df.all$logp <- ifelse(entrez.res.df.all$log2FoldChange < 0, log10(entrez.res.df.all$pvalue), -log10(entrez.res.df.all$pvalue))

# rank all genes by logFC
logp_allgenes<-entrez.res.df.all$logp
names(logp_allgenes)<-entrez.res.df.all$entrezgene_id
h_gene_sets<-msigdbr(species = "Mus musculus", category = "H")
msigdbr_list = split(x = h_gene_sets$entrez_gene, f = h_gene_sets$gs_name)
fgseaRes <- fgsea(msigdbr_list, 
                  logp_allgenes, 
                  minSize=15, 
                  maxSize = 500, 
                  nperm=1000) 
fgsea.sig<-as.data.frame(fgseaRes[which(fgseaRes$padj<0.05),])

# write fgsea sig to file
text<-''
temp<-m$text_matrix(fgsea.sig,text)
temp<-m$write_results_text(temp,dir.24,"msigdbr_pathways.tsv")

pathways<-fgsea.sig$pathway
enrichment_plots<-lapply(pathways,m$plot_enrichment_score,msigdbr_list,logp_allgenes,dir.24)
F1A.IEC.D1<-m$plot_enrichment_score("HALLMARK_INTERFERON_ALPHA_RESPONSE",msigdbr_list,logp_allgenes,dir.24)
F1A.IEC.D1<-F1A.IEC.D1 + ggtitle("IECs D1 - Interferon alpha response") + theme(plot.title = element_text(size = 15))
pdf(file.path(figures_dir,"F1A.IEC.D1.pdf"),5,2)
print(F1A.IEC.D1)
dev.off()

# IECs D3
# rank all genes by logFC
ensembl_entrez_all<-m$ensembl_to_entrez(res.72$ensembl_gene_id)
entrez.res.df.all<-merge(ensembl_entrez_all,res.72, by='ensembl_gene_id',all=TRUE)
# limit to genes that are expressed
entrez.res.df.all<-entrez.res.df.all[!(entrez.res.df.all$baseMean == 0),]
# remove NAs
entrez.res.df.all <- entrez.res.df.all[!is.na(entrez.res.df.all$entrezgene_id),]
entrez.res.df.all <- entrez.res.df.all[!is.na(entrez.res.df.all$pvalue),]
# new column for -log10(pvalue)
entrez.res.df.all$logp <- ifelse(entrez.res.df.all$log2FoldChange < 0, log10(entrez.res.df.all$pvalue), -log10(entrez.res.df.all$pvalue))

# rank all genes by logFC
logp_allgenes<-entrez.res.df.all$logp
names(logp_allgenes)<-entrez.res.df.all$entrezgene_id
h_gene_sets<-msigdbr(species = "Mus musculus", category = "H")
msigdbr_list = split(x = h_gene_sets$entrez_gene, f = h_gene_sets$gs_name)
fgseaRes <- fgsea(msigdbr_list, 
                  logp_allgenes, 
                  minSize=15, 
                  maxSize = 500, 
                  nperm=1000) 
fgsea.sig<-as.data.frame(fgseaRes[which(fgseaRes$padj<0.05),])

# write fgsea sig to file
text<-''
temp<-m$text_matrix(fgsea.sig,text)
temp<-m$write_results_text(temp,dir.72,"msigdbr_pathways.tsv")

pathways<-fgsea.sig$pathway
enrichment_plots<-lapply(pathways,m$plot_enrichment_score,msigdbr_list,logp_allgenes,dir.72)
F1A.IEC.D3<-m$plot_enrichment_score("HALLMARK_INTERFERON_ALPHA_RESPONSE",msigdbr_list,logp_allgenes,dir.72)
F1A.IEC.D3<-F1A.IEC.D3 + ggtitle("IECs D3 - Interferon alpha response") + theme(plot.title = element_text(size = 15))
pdf(file.path(figures_dir,"F1A.IEC.D3.pdf"),5,2)
print(F1A.IEC.D3)
dev.off()

# Caecum
# rank all genes by logFC
ensembl_entrez_all<-m$ensembl_to_entrez(results_d7$ensembl_gene_id)
entrez.res.df.all<-merge(ensembl_entrez_all,results_d7, by='ensembl_gene_id',all=TRUE)
# limit to genes that are expressed
entrez.res.df.all<-entrez.res.df.all[!(entrez.res.df.all$baseMean == 0),]
# remove NAs
entrez.res.df.all <- entrez.res.df.all[!is.na(entrez.res.df.all$entrezgene_id),]
entrez.res.df.all <- entrez.res.df.all[!is.na(entrez.res.df.all$pvalue),]
# new column for -log10(pvalue)
entrez.res.df.all$logp <- ifelse(entrez.res.df.all$log2FoldChange < 0, log10(entrez.res.df.all$pvalue), -log10(entrez.res.df.all$pvalue))

# rank all genes by logFC
logp_allgenes<-entrez.res.df.all$logp
names(logp_allgenes)<-entrez.res.df.all$entrezgene_id
h_gene_sets<-msigdbr(species = "Mus musculus", category = "H")
msigdbr_list = split(x = h_gene_sets$entrez_gene, f = h_gene_sets$gs_name)
fgseaRes <- fgsea(msigdbr_list, 
                  logp_allgenes, 
                  minSize=15, 
                  maxSize = 500, 
                  nperm=1000) 
fgsea.sig<-as.data.frame(fgseaRes[which(fgseaRes$padj<0.05),])

# write fgsea sig to file
text<-''
temp<-m$text_matrix(fgsea.sig,text)
temp<-m$write_results_text(temp,dir.time,"msigdbr_pathways.tsv")

pathways<-fgsea.sig$pathway
enrichment_plots<-lapply(pathways,m$plot_enrichment_score,msigdbr_list,logp_allgenes,dir.time)
F1A.CAECUM<-m$plot_enrichment_score("HALLMARK_INTERFERON_ALPHA_RESPONSE",msigdbr_list,logp_allgenes,dir.time)
F1A.CAECUM<-F1A.CAECUM + ggtitle("Caecum - Interferon alpha response") + theme(plot.title = element_text(size = 15))
pdf(file.path(figures_dir,"F1A.CAECUM.pdf"),5,2)
print(F1A.CAECUM)
dev.off()