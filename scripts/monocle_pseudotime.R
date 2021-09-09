library(Seurat)
library(monocle3)

### params to edit ###

# directory with Seurat object
seurat_dir<-file.path('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/merge/undiff.entero1')

# output directory
dir<-file.path('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/merge/undiff.entero1')

figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'
###

# load the seurat object
seurat_data<-readRDS(file.path(seurat_dir,'seurat_object.rds'))

# maybe - only run analysis on control cells
seurat_data <- subset(seurat_data, subset = condition == "Control")

### Building the necessary parts for a basic cds

# gene annotations

gene_annotation <- as.data.frame(rownames(seurat_data@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat_data@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# cell information

cell_metadata <- as.data.frame(seurat_data@assays[["RNA"]]@counts@Dimnames[[2]], row.names = seurat_data@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
cell_metadata$orig.ident<-seurat_data@meta.data$orig.ident
cell_metadata$cell.type<-seurat_data@meta.data$cell.type
cell_metadata$phase<-seurat_data@meta.data$Phase
cell_metadata$time<-seurat_data@meta.data$time
cell_metadata$seurat_cluster<-seurat_data@active.ident

# counts sparse matrix

New_matrix <- seurat_data@assays[["RNA"]]@counts
#New_matrix <- New_matrix[rownames(seurat_data@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

# make cds object

cds <- new_cell_data_set(expression_matrix,
                        cell_metadata = cell_metadata,
                        gene_metadata = gene_annotation)


# normalize and pre-process

cds <- preprocess_cds(cds, num_dim = 100)

# remove batch effects 

cds <- align_cds(cds, alignment_group = "orig.ident")

# dimensional reduction

cds <- reduce_dimension(cds)

pdf(file.path(dir,'by_seurat_cluster.pdf'))
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_cluster", reduction_method = "UMAP",label_cell_groups=FALSE)
dev.off()

pdf(file.path(dir,'by_cell_type.pdf'))
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type", reduction_method = "UMAP",label_cell_groups=FALSE)
dev.off()

#cluster

cds <- cluster_cells(cds, reduction_method = "UMAP")

pdf(file.path(dir,'by_monocle_cluster.pdf'))
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cluster",label_cell_groups=FALSE)
dev.off()

pdf(file.path(dir,'by_partition.pdf'))
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "partition",label_cell_groups=FALSE)
dev.off()

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "seurat_cluster",
           label_groups_by_cluster=TRUE,
           label_leaves=TRUE,
           label_branch_points=TRUE)

pdf(file.path(dir,'by_phase.pdf'))
plot_cells(cds,
           color_cells_by = "phase",
           label_cell_groups = FALSE
)
dev.off()

cds <- order_cells(cds)

pdf(file.path(dir,'by_psuedotime.pdf'))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()

stem_features <- c(
  "Lgr5", "Lrig1", "Sox9", "Myc", "Hopx", "Ascl2", "Smoc2", "Rfp43", "Cd44", "Surviving", "Pcna", "Ki67", "Spdef", "Atoh1", "Hes1"
)

pdf(file.path(dir,'undiff_genes.pdf'))
plot_cells(cds,
           genes = stem_features,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()


cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(cds_pr_test_res[cds_pr_test_res$morans_I>0.6 & cds_pr_test_res$q_value < 0.05,] )

plot_cells(cds, genes=c("Top2a", "Ube2c"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

gene_module_df <- as.data.frame(find_gene_modules(cds[pr_deg_ids,], resolution=c(0.1,10^seq(-6,-1))))
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$seurat_cluster)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

plot_cells(cds,
           genes=gene_module_df[gene_module_df$module=="2","id"],
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

stem_cds <- cds[rowData(cds)$gene_short_name %in% c("Ube2c","Top2a","Lgr5"),
                       ]

plot_genes_in_pseudotime(stem_cds,
                         color_cells_by="seurat_cluster",
                         min_expr=0.5)

# goblets

goblets<-choose_cells(cds)
goblet_test_res <- graph_test(goblets, neighbor_graph="principal_graph", cores=4)
goblet_ids <- row.names(goblet_test_res[goblet_test_res$morans_I>0.6 & goblet_test_res$q_value < 0.05 ,])
gene_modules_goblets <- find_gene_modules(goblets[goblet_ids,], resolution=0.001)
agg_mat <- aggregate_gene_expression(goblets, gene_modules_goblets)
module_dendro <- hclust(dist(agg_mat))
gene_modules_goblets$module <- factor(gene_modules_goblets$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(goblets,
           genes=gene_modules_goblets,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# enterocytes

enterocytes<-choose_cells(cds)
enterocyte_test_res <- graph_test(enterocytes, neighbor_graph="principal_graph", cores=4)
enterocyte_ids <- row.names(enterocyte_test_res[enterocyte_test_res$morans_I>0.6 & enterocyte_test_res$q_value < 0.05 ,])
gene_modules_enterocytes <- find_gene_modules(enterocytes[enterocyte_ids,])
cell_group_df <- data.frame(cell=row.names(colData(enterocytes)), 
                            cell_group=clusters(enterocytes)[colnames(enterocytes)])
agg_mat <- aggregate_gene_expression(enterocytes, gene_modules_enterocytes, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Cluster ", colnames(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")
module_dendro <- hclust(dist(agg_mat))
gene_modules_enterocytes$module <- factor(gene_modules_enterocytes$module, 
                                      levels = row.names(agg_mat)[module_dendro$order])

plot_cells(enterocytes,
           genes=gene_modules_enterocytes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# cycling

cyclings<-choose_cells(cds)
cycling_test_res <- graph_test(cyclings, neighbor_graph="principal_graph", cores=4)
cycling_ids <- row.names(cycling_test_res[cycling_test_res$morans_I>0.4 & cycling_test_res$q_value < 0.05 ,])
gene_modules_cyclings <- find_gene_modules(cyclings[cycling_ids,],resolution=0.0001)
agg_mat <- aggregate_gene_expression(cyclings, gene_modules_cyclings)
module_dendro <- hclust(dist(agg_mat))
gene_modules_cyclings$module <- factor(gene_modules_cyclings$module, 
                                          levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cyclings,
           genes=gene_modules_cyclings,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

cyclings <- cluster_cells(cyclings, reduction_method = "UMAP",resolution=0.0005)
plot_cells(cyclings, label_groups_by_cluster=FALSE,  color_cells_by = "cluster")
cyclings<-learn_graph(cyclings)
plot_cells(cyclings, label_groups_by_cluster=FALSE,  color_cells_by = "cluster")

cycling_test_res <- graph_test(cyclings, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(cycling_test_res[cycling_test_res$morans_I>0.6 & cycling_test_res$q_value < 0.05 ,])
gene_module_df <- as.data.frame(find_gene_modules(cyclings[pr_deg_ids,], resolution=1e-2))
cell_group_df <- data.frame(cell=row.names(colData(cyclings)), 
                                cell_group=clusters(cyclings)[colnames(cyclings)])
agg_mat <- aggregate_gene_expression(cyclings, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Cluster ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
