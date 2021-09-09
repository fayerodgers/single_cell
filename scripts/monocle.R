import("monocle3")
import("Seurat")


seurat_to_monocle <- function(seurat_data){

	# gene annotations

   gene_annotation <- as.data.frame(rownames(seurat_data@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat_data@reductions[["pca"]]@feature.loadings))
#  gene_annotation <- as.data.frame(rownames(seurat_data@assays[["RNA"]]@counts), row.names = rownames(seurat_data@assays[["RNA"]]@counts))
  colnames(gene_annotation) <- "gene_short_name"

	# cell information

	cell_metadata <- as.data.frame(seurat_data@assays[["RNA"]]@counts@Dimnames[[2]], row.names = seurat_data@assays[["RNA"]]@counts@Dimnames[[2]])
	colnames(cell_metadata) <- "barcode"
	cell_metadata$orig.ident<-seurat_data@meta.data$orig.ident
	cell_metadata$cell.type<-seurat_data@meta.data$cell.type
	cell_metadata$phase<-seurat_data@meta.data$Phase
	cell_metadata$time<-seurat_data@meta.data$time
	cell_metadata$seurat_cluster<-seurat_data@active.ident
	#cell_metadata$condition<-seurat_data$condition

	# counts sparse matrix

	New_matrix <- seurat_data@assays[["RNA"]]@counts
	New_matrix <- New_matrix[rownames(seurat_data@reductions[["pca"]]@feature.loadings), ]
	expression_matrix <- New_matrix

	# make cds object

	cds <- new_cell_data_set(expression_matrix,
                        cell_metadata = cell_metadata,
                        gene_metadata = gene_annotation)
	return(cds)
}
