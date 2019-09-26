#!/usr/bin/env Rscript
m <- modules::use("SC.R")
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

parser <- ArgumentParser()
parser$add_argument("--results_directory",help="directory for files to be written to")
parser$add_argument("--samples", help="File with a list of sample IDs to analyse, one per line")
parser$add_argument("--metadata", help="File with sample metadata")
parser$add_argument("--cellranger_version", help="cellranger131, cellranger211 or cellranger302")
parser$add_argument("--mito_cutoff", type="integer", help="Filter out cells with % reads aligning to mitochondrial genes greater than this (default: no cutoff)",default=100)
parser$add_argument("--nfeatures_cutoff", type="integer",help="Filter out cells with fewer features than this (default: no filter)", default=0)
parser$add_argument("--ncells_cutoff",type="integer",help="Only include features (genes) that are expressed in more than this number of cells (default: no filter)",default=0)
parser$add_argument("--variable_features",type="integer",help="Take top x variable genes (default: 2000)",default=2000)
parser$add_argument("--dimensions_to_assess", type="integer",help="Number of prinicpal components to calculate (default: 30)",default=30)
parser$add_argument("--dimensions_to_analyse",type="integer",help="Number of prinicpal components to include in UMAP and neighbour finding (default: 25)", default=25)
parser$add_argument("--resolution",type="integer",help="Resolution for identifying clusters (default: 1)",default=1)

args <- parser$parse_args()

samples<-scan(args$samples,what=character())
meta_data<-read.table(args$metadata,header=T)

seurat_objects<-lapply(samples,m$normalize_data,
                       args$cellranger_version,
                       meta_data,
                       args$variable_features,
                       args$ncells_cutoff,
                       args$mito_cutoff,
                       args$nfeatures_cutoff)
combined<-m$integrate_datasets(seurat_objects,args$dimensions_to_assess)
combined<-m$do_clustering(combined,args$dimensions_to_analyse,args$resolution)
dim.plots<-m$plot_UMAPS(combined)
markers<-FeaturePlot(combined, features = c("Aqp8","Krt20","Muc2","Chga","Lgr5","Dclk1","Cdk4","Il33","Ly6a"), min.cutoff="q9",pt.size = 0.1)

#save data structure
saveRDS(combined,paste0(args$results_directory,"/seurat_object.rds"))

#Printing plots
pdf(paste0(args$results_directory,"/elbowplot.pdf"))
ElbowPlot(combined,ndims=args$dimensions_to_assess)
dev.off()

mito.vln<-VlnPlot(combined, features = c("percent.mito"), cols=rep("red",length(samples)), pt.size = 0,group.by='orig.ident',do.return=T)
mito.vln<-mito.vln + theme(axis.text.x=element_text(size=8)) + labs(x="",y="% mitochondrial genes", title = args$cellranger_version)

nfeatures.vln<-VlnPlot(combined, features = c("nFeature_RNA"), cols=rep("forestgreen",length(samples)), pt.size = 0,group.by='orig.ident',do.return=T)
nfeatures.vln<-nfeatures.vln + theme(axis.text.x=element_text(size=8)) + labs(x="",y="Number of genes", title = args$cellranger_version)

pdf(paste0(args$results_directory,"/percentmito.pdf"),12,6)
print(mito.vln)
dev.off()

pdf(paste0(args$results_directory,"/nfeatures.pdf"),12,6)
print(nfeatures.vln)
dev.off()

pdf(paste0(args$results_directory,"/umaps.pdf"),20,20)
print(dim.plots)
dev.off()

pdf(paste0(args$results_directory,"/markers.pdf"),20,20)
print(markers)
dev.off()

