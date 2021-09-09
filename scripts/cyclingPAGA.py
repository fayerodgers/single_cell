import numpy as np 
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, facecolor='white')


#pre-processing function
def pre_process(adata):
	#normalise to 10 000 reads per cell
	sc.pp.normalize_total(adata, target_sum=1e4)	
	#log
	sc.pp.log1p(adata)
	#find variable genes
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	#filter to take only variable genes 
	#adata = adata[:, adata.var.highly_variable]	
	#regress out batch effects
	sc.pp.regress_out(adata, ['orig_ident'])
	#scale data
	sc.pp.scale(adata, max_value=10)
	#calculate PCA
	sc.tl.pca(adata, svd_solver='arpack', n_comps=50)	
	return adata



#define marker genes for plotting
genes=['Muc2','Aqp8','Krt20','Isg15','Reg3g','Chga','Top2a','Lgr5','Smoc2','Ascl2','Cdk4','Rps15a']

cycling = sc.read_loom("/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/cycling.loom", sparse=True)

sc.pp.filter_genes(cycling, min_counts=1)
cycling = pre_process(cycling)

sc.pl.pca_variance_ratio(cycling, log=True, n_pcs = 50, save = ".cycling.png")

sc.pp.neighbors(cycling, n_pcs=35, n_neighbors=50)

sc.tl.umap(cycling)
sc.pl.umap(cycling, color=genes, save=".cycling.png")

sc.tl.leiden(cycling, resolution = 1)
sc.pl.umap(cycling, color = ['leiden', 'ClusterName'], save = ".clusters.cycling.png")

sc.tl.rank_genes_groups(cycling, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(cycling, n_genes=20, sharey=False,save=".cycling.png")

genes.append('leiden')
genes.append('ClusterName')

sc.tl.paga(cycling, groups='leiden')
sc.pl.paga(cycling, color=genes, threshold = 0.45, save = ".cycling.png")

sc.tl.draw_graph(cycling, init_pos='paga')
sc.pl.draw_graph(cycling, color=genes, save = ".cycling.png")

genes=['Npm1','Ncl','Top2a','Cdk1','Papss2','Mpp6','Cdkn3','Cdc20','Hepacam2','Fcgbp']

stems=['Top2a','Cdk1','Tubb5','Smc2','Smc4','Pbk','Ube2c','Tuba1b','Mki67','Cdca8','Aurkb','Prc1','Hist1h2ap','Birc5','Kif11','Hmgb2','Nusap1','Cks1b','H2afx','Cenpe']
