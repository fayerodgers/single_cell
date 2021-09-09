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
	adata = adata[:, adata.var.highly_variable]
    #remove genes with 0 counts	
	#adata = adata[:, adata.X.sum(axis=0) > 0]
	#regress out batch effects
	sc.pp.regress_out(adata, ['orig_ident'])
	#scale data
	sc.pp.scale(adata, max_value=10)
	#calculate PCA
	sc.tl.pca(adata, svd_solver='arpack', n_comps=100)	
	return adata



#define marker genes for plotting
genes=['Muc2','Aqp8','Krt20','Isg15','Reg3g','Chga','Top2a','Lgr5','Smoc2','Ascl2','Fer1l6']

control = sc.read_loom("/Users/fr7/git_repos/single_cell/experiment_4/FINAL/merge/control.loom", sparse=True)

#pre-process
control = pre_process(control)

#plots to select number of components for clustering
sc.pl.pca_variance_ratio(control, log=True, n_pcs = 100)


#compute the neighbourhood graphs
sc.pp.neighbors(control, n_pcs=35, n_neighbors=4) #might want to play with n_neighbours

#plot UMAPs
sc.tl.umap(control)
sc.pl.umap(control, color=genes, save=".control.png")
sc.pl.umap(control, color='ClusterName', save=".control.png")

#clustering
sc.tl.leiden(control, resolution = 0.3)
sc.pl.umap(control, color = ['leiden', 'ClusterName'], save = ".clusters.control.png")


#denoise the graph
sc.tl.diffmap(control)
sc.pp.neighbors(control, n_neighbors=4, use_rep='X_diffmap')
sc.tl.draw_graph(control, init_pos='paga')
sc.pl.draw_graph(control, color = 'ClusterName')

#find markers
sc.tl.rank_genes_groups(control, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(control, n_genes=10, sharey=False,save=".control.png")



#get tables of marker genes

genes.append('leiden')
genes.append('ClusterName')

genes = ['Ube2c','Tubb5','Top2a','Gsdmc4','Aldoa','Gsdmc2','Krt18','Slc12a2','Gpx2','Cdkn3','Cdc20','Ccnb2','Cpn1','Tmigd1','Rbp2','Sytl2','Hepacam2','Muc2','ClusterName', 'Aqp8', 'Krt20']

#run PAGA
new_colors = np.array(control.uns['ClusterName_colors'])
new_colors[[0]] = "#5c1a33" #DSC
new_colors[[1]] = "#1f78b4" #Early entero.1
new_colors[[2]] = "#006D2C" #Early entero.2
new_colors[[3]] = "#b15928" #Ee
new_colors[[4]] = "#6a3d9a" #Entero.AMP
new_colors[[5]] = "#e31a1c" #Entero.Isg15
new_colors[[6]] = "#cab2d6" #Enterocyte
new_colors[[7]] = "#a6cee3" #G2M-Stem/TA cells
new_colors[[8]] = "#33a02c" #Goblet
new_colors[[9]] = "#C7E9C0" #Progenitor Entero.1
new_colors[[10]] = "#ff7f00" #Progenitor Entero.2
new_colors[[11]] = "#fb9a99" #S-Stem/TA cells
new_colors[[12]] = "#fdbf6f" #Tuft


control.uns['ClusterName_colors'] = new_colors

sc.tl.paga(control, groups='ClusterName')
sc.pl.paga(control, threshold = 0.3, fontoutline=1, edge_width_scale=0.3, fontsize=4.5, save = ".control.png", frameon=False)

sc.tl.umap(control, init_pos='paga', min_dist = 0.001, spread = 2 )
sc.pl.umap(control, color=genes, save = ".control.png")


#draw graphs

sc.tl.draw_graph(control, init_pos='paga')
sc.pl.draw_graph(control, color=genes, save = ".control.png")



late_entero_markers = ['Aqp8','Rhoc','Car4',"Cpn1",'Selenop','Clca4a','Irf7','Rbp2','S100a14','Capg','Dgat2']
sc.pl.draw_graph(control, color=late_entero_markers, save = "late_entero_markers.control.png")

early_entero_markers = ['Gsdmc2','Slc20a1','Mpp6','Gsdmc4','Papss2','Slc37a2','Pycard','Gstm1','Hmox1','Mt1','Mt2','Sult1a1']
sc.pl.draw_graph(control, color=early_entero_markers, save = "early_entero_markers.control.png")


goblet_markers = ['Fxyd3','Klk1','Qsox1','Smim14','Agr2','Ramp1','Muc2','Fcgbp','Tff3',"Sytl2",'Hepacam2','Tspan13','Reg4']
sc.pl.draw_graph(control, color=goblet_markers, save = ".goblet_markers.control.png")

#choose root cells

control.uns['iroot'] = np.flatnonzero(control.obs['leiden']  == '3')[0]


#pseudotime
sc.tl.diffmap(control,n_comps=35)
sc.tl.dpt(control)
sc.pl.umap(control, color=['dpt_pseudotime'], legend_loc='on data', save = "dpt_pseudotime.control.png")


####

