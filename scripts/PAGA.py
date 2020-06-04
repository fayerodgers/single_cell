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
	sc.tl.pca(adata, svd_solver='arpack', n_comps=100)	
	return adata



#define marker genes for plotting
genes=['Muc2','Aqp8','Krt20','Isg15','Reg3g','Chga','Top2a','Lgr5','Smoc2','Ascl2','Fer1l6']

control = sc.read_loom("/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/control.loom", sparse=True)
day1 = sc.read_loom("/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/day1_classify/day.1.loom", sparse=True)
day3 = sc.read_loom("/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/day3_classify/day.3.loom", sparse=True)

#pre-process
control = pre_process(control)
day1 = pre_process(day1)
day3 = pre_process(day3)

#plots to select number of components for clustering
sc.pl.pca_variance_ratio(control, log=True, n_pcs = 100, save = ".control.png")
sc.pl.pca_variance_ratio(day1, log=True, n_pcs = 100, save = ".day1.png")
sc.pl.pca_variance_ratio(day3, log=True, n_pcs = 100, save = ".day3.png")

#compute the neighbourhood graphs
sc.pp.neighbors(control, n_pcs=35, n_neighbors=4) #might want to play with n_neighbours
sc.pp.neighbors(day1, n_pcs=35, n_neighbors=4)
sc.pp.neighbors(day3, n_pcs=35, n_neighbors=4)

#plot UMAPs
sc.tl.umap(control)
sc.pl.umap(control, color=genes, save=".control.png")

sc.tl.umap(day1)
sc.pl.umap(day1, color=genes, save=".day1.png")

sc.tl.umap(day3)
sc.pl.umap(day3, color=genes, save=".day3.png")


#clustering
sc.tl.leiden(control, resolution = 1)
sc.pl.umap(control, color = ['leiden', 'ClusterName'], save = ".clusters.control.png")

sc.tl.leiden(day1, resolution = 1)
sc.pl.umap(day1, color = ['leiden', 'ClusterName'], save = ".clusters.day1.png")

sc.tl.leiden(day3, resolution = 1)
sc.pl.umap(day3, color = ['leiden', 'ClusterName'], save = ".clusters.day3.png")

#denoise the graph
#sc.tl.diffmap(control)
#sc.pp.neighbors(control, n_neighbors=4, use_rep='X_diffmap')
#sc.tl.draw_graph(control)
#sc.pl.draw_graph(control, color = 'ClusterName')

#find markers
sc.tl.rank_genes_groups(control, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(control, n_genes=20, sharey=False,save=".control.png")

sc.tl.rank_genes_groups(day1, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(day1, n_genes=20, sharey=False,save=".day1.png")

sc.tl.rank_genes_groups(day3, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(day3, n_genes=20, sharey=False,save=".day3.png")

#get tables of marker genes

genes.append('leiden')
genes.append('ClusterName')

#run PAGA

sc.tl.paga(control, groups='leiden')
sc.pl.paga(control, color=genes, threshold = 0.1, save = ".control.png")

sc.tl.paga(day1, groups='leiden')
sc.pl.paga(day1, color=genes, threshold = 0.1, save = ".day1.png")

sc.tl.paga(day3, groups='leiden')
sc.pl.paga(day3, color=genes, threshold = 0.1, save = ".day3.png")

#draw graphs

sc.tl.draw_graph(control, init_pos='paga')
sc.pl.draw_graph(control, color=genes, save = ".control.png")

sc.tl.draw_graph(day1, init_pos='paga')
sc.pl.draw_graph(day1, color=genes, save = ".day1.png")

sc.tl.draw_graph(day3, init_pos='paga')
sc.pl.draw_graph(day3, color=genes, save = ".day3.png")

late_entero_markers = ['Aqp8','Rhoc','Car4',"Cpn1",'Selenop','Clca4a','Irf7','Rbp2','S100a14','Capg','Dgat2']
sc.pl.draw_graph(control, color=late_entero_markers, save = "late_entero_markers.control.png")

early_entero_markers = ['Gsdmc2','Slc20a1','Mpp6','Gsdmc4','Papss2','Slc37a2','Pycard','Gstm1','Hmox1','Mt1','Mt2','Sult1a1']
sc.pl.draw_graph(control, color=early_entero_markers, save = "early_entero_markers.control.png")


goblet_markers = ['Fxyd3','Klk1','Qsox1','Smim14','Agr2','Ramp1','Muc2','Fcgbp','Tff3',"Sytl2",'Hepacam2','Tspan13','Reg4']
sc.pl.draw_graph(control, color=goblet_markers, save = ".goblet_markers.control.png")

#choose root cells

control.uns['iroot'] = np.flatnonzero(control.obs['leiden']  == '1')[0]
day1.uns['iroot'] = np.flatnonzero(day1.obs['leiden']  == '18')[0]
day3.uns['iroot'] = np.flatnonzero(day3.obs['leiden']  == '22')[0]

#pseudotime
sc.tl.dpt(control)
sc.pl.draw_graph(control, color=['leiden', 'dpt_pseudotime'], legend_loc='on data', save = "dpt_pseudotime.control.png")

sc.tl.dpt(day1)
sc.pl.draw_graph(day1, color=['leiden', 'dpt_pseudotime'], legend_loc='on data', save = "dpt_pseudotime.day1.png")

sc.tl.dpt(day3)
sc.pl.draw_graph(day3, color=['leiden', 'dpt_pseudotime'], legend_loc='on data', save = "dpt_pseudotime.day3.png")


####

