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
library("scales")
library("cowplot")
library("stringr")
library("clusterProfiler")
library("goseq")
library("fgsea")
library("GO.db")
library("msigdbr")
library("pathview")


m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")

###### things to edit ######
#paths and metadata
maindir<-'/Users/fr7/git_repos/single_cell/experiment_4/bulk_libraries'
figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'
files_dir<-'/Users/fr7/git_repos/single_cell/SupplementaryFiles'
metadata<-read.table(file.path(maindir,'metadata.txt'),header=TRUE,stringsAsFactors = TRUE)
tx2gene <- read.table('~/git_repos/Trichuris_transwells/transcripts2genes.E98.tsv', header = T)
heatmap_genes<-'/Users/fr7/git_repos/single_cell/experiment_4/bulk_libraries/genes_to_plot.tsv'

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
design<-c('sorting','time.status')
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

#write supplementary file
res.24.sig<-res.24[which(res.24$padj < 0.05),]
text<-"Differential expression of genes in caecal epithelial cells 1 day post Trichuris infection (Wald test; control (d1) v infected (d1))"
res.24.text<-m$text_matrix(res.24.sig,text)
res.24.text<-m$write_results_text(res.24.text,files_dir,"SF.5.tsv")

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

#write supplementary file
res.72.sig<-res.72[which(res.72$padj < 0.05),]
text<-"Differential expression of genes in caecal epithelial cells 3 days post Trichuris infection (Wald test; control (d3) v infected (d3))"
res.72.text<-m$text_matrix(res.72.sig,text)
res.72.text<-m$write_results_text(res.72.text,files_dir,"SF.6.tsv")

#### bar graph of number of regulated genes ####
up.24 <- length(rownames(res.24[which(res.24$padj < 0.05 & res.24$log2FoldChange > 1),]))
down.24 <- length(rownames(res.24[which(res.24$padj < 0.05 & res.24$log2FoldChange < -1),]))
up.72 <- length(rownames(res.72[which(res.72$padj < 0.05 & res.72$log2FoldChange > 1),]))
down.72 <- length(rownames(res.72[which(res.72$padj < 0.05 & res.72$log2FoldChange < -1),]))

regulated_genes <- data.frame(
  "time"= c("D1","D1","D3","D3"),
  "regulation"=c("Up","Down","Up","Down"),
  "ngenes"=c(up.24,down.24,up.72,down.72)
)

fig.1d<-ggplot(regulated_genes,x=time,y=ngenes)+
  geom_col(aes(time,ngenes,fill=regulation))+
  # scale_fill_manual(values=cbPalette)+
  #  facet_grid(~eval(substitute(facet), df),scales="free", space="free_x", labeller=labeller(labels))+
  scale_fill_manual(values=c('black','grey'))+
  theme(text = element_text(size=8),legend.title = element_text(size = 5), legend.text = element_text(size = 5),
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,-10,-10,-10) )+
  ylab('Regulated genes')+
  xlab('')+
  labs(fill='Regulation\n(Log2FC > 1)')+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))


pdf(file.path(dir,'ngenes.pdf'))
print(fig.1d)
dev.off()

pdf(file.path(figures_dir,'fig.1d.pdf'),1,2)
print(fig.1d)
dev.off()


#### GSEA #####

# get transcript length data from dds object

tx.lengths<-rowMedians(dds.24@assays@data$avgTxLength)
names(tx.lengths)<-names(dds.24)

# GO - 24h

# remove genes that are not expressed
res.24.expressed<-res.24[!(res.24$baseMean == 0),]

sig.24<-as.integer(res.24.expressed$padj<0.05 & !is.na(res.24.expressed$padj))
names(sig.24)<-res.24.expressed$ensembl_gene_id
tx.lengths<-tx.lengths[names(tx.lengths) %in% names(sig.24) ]

pwf <- nullp(sig.24, "mm10", "ensGene",bias.data=tx.lengths)

go.24 <- goseq(pwf, "mm10","ensGene", test.cats=c("GO:BP"))

go.24<-m$write_results(go.24,dir.24,"GO.BP")
# write nice supplementary file
text<-"GO term (BP) enrichment in significantly regulated genes in caecum epithelial cells 1 day post infection"
go.24.text<-m$text_matrix(go.24,text)
go.24.text<-m$write_results_text(go.24.text,files_dir,"SF.7.tsv")

# kegg - 24h
# search_kegg_organism('mmu', by='kegg_code')
# use biomart to convert Ensembl gene IDs to Entrez

entrez.df<-m$ensembl_to_entrez(res.24[which(res.24$padj<0.05),'ensembl_gene_id'])
entrez_genes <- entrez.df$entrezgene_id[!is.na(entrez.df$entrezgene_id)]
kk <- enrichKEGG(gene = entrez_genes, organism = 'mmu') # no significant KEGG pathways at 24h

# pathways - 24h
# rank all genes by logFC
ensembl_entrez_all<-m$ensembl_to_entrez(res.24$ensembl_gene_id)
entrez.res.df.all<-merge(ensembl_entrez_all,res.24, by='ensembl_gene_id',all=TRUE)
#limit to genes that are expressed
entrez.res.df.all<-entrez.res.df.all[!(entrez.res.df.all$baseMean == 0),]
#remove NAs
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
F25<-m$plot_enrichment_score("HALLMARK_INTERFERON_ALPHA_RESPONSE",msigdbr_list,logp_allgenes,dir.24)
F25<-F25 + ggtitle("Interferon alpha response (epithelium, D1)")
pdf(file.path(figures_dir,"F25.pdf"),5,2)
print(F25)
dev.off()


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

go.72<-m$write_results(go.72,dir.72,"GO.BP")

# write nice supplementary file
text<-"GO term (BP) enrichment in significantly regulated genes in caecum epithelial cells 3 days post infection"
go.72.text<-m$text_matrix(go.72,text)
go.72.text<-m$write_results_text(go.72.text,files_dir,"SF.8.tsv")

# kegg - 72h
# search_kegg_organism('mmu', by='kegg_code')
# use biomart to convert Ensembl gene IDs to Entrez

entrez.df<-m$ensembl_to_entrez(res.72[which(res.72$padj<0.05),'ensembl_gene_id'])
entrez.df <- entrez.df[!is.na(entrez.df$entrezgene_id),]
kk <- enrichKEGG(gene = entrez.df$entrezgene_id, organism = 'mmu')
kk<-m$write_results(kk,dir.72,"KEGG")

# plot KEGG pathways
entrez.res.df<-merge(entrez.df,res.72)
logFC<-entrez.res.df$log2FoldChange
names(logFC)<-entrez.res.df$entrezgene_id

keggs<-kk$ID
draw<-lapply(keggs,m$draw_kegg_pathway,logFC,dir.72)

# pathways - 72h
# rank all genes by logFC
ensembl_entrez_all<-m$ensembl_to_entrez(res.72$ensembl_gene_id)
entrez.res.df.all<-merge(ensembl_entrez_all,res.72, by='ensembl_gene_id',all=TRUE)
#limit to genes that are expressed
entrez.res.df.all<-entrez.res.df.all[!(entrez.res.df.all$baseMean == 0),]
#remove NAs
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
temp<-m$write_results_text(temp,dir.72,"msigdbr_pathways")

pathways<-fgsea.sig$pathway
enrichment_plots<-lapply(pathways,m$plot_enrichment_score,msigdbr_list,logp_allgenes,dir.72)
F26<-m$plot_enrichment_score("HALLMARK_INTERFERON_ALPHA_RESPONSE",msigdbr_list,logp_allgenes,dir.72)
F26<-F26 + ggtitle("Interferon alpha response (epithelium, D3)")
pdf(file.path(figures_dir,"F26.pdf"),5,2)
print(F26)
dev.off()

# heatmaps
genes<-scan(file = heatmap_genes, sep = "\n", what = character())
coul <- colorRampPalette(colors = c("blue","white","red"))(15)

# 24h - choose
genes.24<-res.24.sig[which(res.24.sig$external_gene_name %in% genes & abs(res.24.sig$log2FoldChange)>4.5),]
fpkms.24<-fpkm(dds.24)[genes.24$ensembl_gene_id,]

rownames(fpkms.24)<-genes.24$external_gene_name
#colnames(fpkms)<-metadata.all$time.status
annotation<-data.frame(metadata.24$time.status)
row.names(annotation)<-metadata.24$sample

#heatmap of replicates
pdf(file.path(dir.24,"heatmap_replicates.pdf"))
pheatmap(log((fpkms.24 + 1), base=2),color = coul, cluster_cols=F,treeheight_row = 0, treeheight_col = 0 , annotation_col = annotation , show_colnames = F, cluster_rows = T, cellwidth = 5, cellheight = 5, fontsize = 5, annotation_legend = FALSE)
dev.off()

# 72h - choose
genes.72<-res.72.sig[which(res.72.sig$external_gene_name %in% genes & abs(res.72.sig$log2FoldChange)>1.5),]
fpkms.72<-fpkm(dds.72)[genes.72$ensembl_gene_id,]

rownames(fpkms.72)<-genes.72$external_gene_name
#colnames(fpkms)<-metadata.all$time.status
annotation<-data.frame(metadata.72$time.status)
row.names(annotation)<-metadata.72$sample

#heatmap of replicates
pdf(file.path(dir.72,"heatmap_replicates.pdf"),3,4)
pheatmap(log((fpkms.72 + 1), base=2),col = coul, cluster_cols=T,treeheight_row = 0, treeheight_col = 0 , annotation_col = annotation , show_colnames = F, cluster_rows = T, cellwidth = 5, cellheight = 5, fontsize = 5, annotation_legend = FALSE, legend = FALSE)
dev.off()


#### GO term dot plot - deprecated ####

#BP GO terms with an adjusted pathway p value < 0.01 at at least one time point
#were selected and manually curated for redundancy

go<-read.csv(file.path(dir.all,'curated_go_terms.csv'),header=TRUE, sep =",")

#restrict to terms that are significant in at least one time point 
sig_terms<-go[which(go$PathwayPadj<0.01),"PathwayId"]
go.sig<-go[which(go$PathwayId %in% sig_terms),]


fig.1e<-ggplot(go.sig,x=time,y=PathwayName)+
  geom_point(aes(x=time,y=PathwayName, size=ngenes, col=PathwayPadj))+
  scale_color_gradient(limits = c(0.0001,0.01),
                       oob=squish, 
                       name = "Pathway p-value\n(adjusted)", 
                       breaks = c(0.0001,0.01),
                       labels = c("< 0.0001","> 0.01"),
                       guide = guide_colourbar())+
  scale_size_continuous(name="Number of\ngenes", range = c(0.1,2))+
  theme(text = element_text(size=5),
        axis.text.x=element_text(size=7),
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-10,-10) )+
  labs(x=NULL,y=NULL)+
  guides(color = guide_colourbar(barwidth = 0.5, barheight = 2), 
         size = guide_legend(keywidth = 0.5, keyheight = 0.5))

pdf(file.path(figures_dir,'fig.1e.pdf'),2.5,1.8)
print(fig.1e)
dev.off()



###### heatmaps across all timepoints - deprecated ######

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

#with 24h + 72h control samples combined - choose this one for figure
fpkm.means.1<-data.frame("C" = rowMeans(fpkms[,c(control.24h,control.72h)]),
                       "D1" = rowMeans(fpkms[,infected.24h]),
                       "D3" = rowMeans(fpkms[,infected.72h])
)


#GO term- gene associations for annotations
go.gene<-read.table(file.path(dir,"goterm-gene.txt"))
colnames(go.gene)<-c("term","gene")

gos.to.plot<-go$PathwayId
genes.to.plot<-rownames(fpkm.means.1)
go.annotation<-data.frame(matrix(ncol=length(gos.to.plot),nrow=length(genes.to.plot)))
colnames(go.annotation)<-gos.to.plot
rownames(go.annotation)<-genes.to.plot

for (goterm in gos.to.plot){
  temp.list<-c()
  for (gene in genes.to.plot){
    if (length(rownames(go.gene[which(go.gene$term == goterm & go.gene$gene == gene),]) > 0) ){
      temp.list<-c(temp.list,1)
    }
    else{
      temp.list<-c(temp.list,NA)
    }
  }
  go.annotation[[goterm]]<-temp.list
}

#remove colons from col names
names(go.annotation)<-str_replace_all(names(go.annotation), c(":" = "."))
#make a temp merged table for ordering for heat map
go.joined<-merge(fpkm.means.1,go.annotation,by=0)
go.joined<-go.joined[with(go.joined, order(GO.0051607,-D3)), ]

#separate tables again
fpkm.means.1<-go.joined[,c("C","D1","D3")]
rownames(fpkm.means.1)<-go.joined$Row.names

#replace colnames with full names in go.annotation
full.names<-go[go$PathwayId %in% gos.to.plot,"PathwayName"]
colnames(go.annotation)<-full.names

#remove columns that are all NA
go.annotation <- go.annotation[,colSums(is.na(go.annotation))<nrow(go.annotation)]

pheatmap(log((fpkm.means.1 + 1), base=2),
                 cluster_cols=F,
                 cluster_rows=F,
                 treeheight_row = 0,
                 treeheight_col = 0,
                 cellwidth = 5,
                 cellheight = 5,
                 angle_col=90, 
                 annotation_row = go.annotation, 
                 annotation_legend = FALSE, 
                 annotation_names_row = TRUE,
                 fontsize = 5,
                 ann_colors="black",
                filename=file.path(figures_dir,"fig.1f.pdf"))

#consider extra ordering

pdf(file.path(figures_dir,"fig.1f.pdf"))
fig.1f
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




