###Import dependent python packages
import igraph
import louvain
import numpy as np
import pandas as pd
import scrublet as scr
import scanpy as sc
import matplotlib.pyplot as pl
from matplotlib import rcParams
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi_save=300)
## RNA velocity
import scvelo as scv
import loompy as lp
import re
scv.logging.print_version()
pl.rcParams['pdf.fonttype'] = 42

#######doublets detection##########
###### use CS10_YS as an example######
adata = sc.read_10x_mtx(
        './CS10_YS/outs/filtered_feature_bc_matrix/',
        var_names='gene_symbols',
        cache=True)

adata.var_names_make_unique()
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata

#####calculate highly variable_genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ['n_counts'])
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata)

######remove doublet cells detected by scrublet######

adata_dr = adata.copy()
scrub = scr.Scrublet(adata_dr.raw.X)
adata_dr.obs['doublet_scores'], adata_dr.obs['predicted_doublets'] = scrub.scrub_doublets()
    #take doublet cells based on expected doublets based on number of cells (1k cells = 1% increase)
expected_doublets = np.int(len(adata_dr)/100000*len(adata_dr))
doublet_cells = adata_dr.obs['doublet_scores'].sort_values(ascending=False).head(n=expected_doublets).index
adata_dr.obs['man_doublets'] = False
adata_dr.obs.loc[doublet_cells,'man_doublets'] = True

###rewrite 'doublets' into meta.data
adata_dr.obs['man_doublets'].value_counts()
adata.obs["doublets"] = adata_dr.obs['man_doublets']
adata.obs["doublets"] = adata.obs["doublets"].astype("category")

pd.DataFrame(adata_dr.obs).to_csv("CS10_YS_doublet_annot.csv")

###########################################################
####construct object according to the seurat results#######
results_file = './Ery.integrated.h5ad' 
adata = sc.read_csv('./Ery.integrated.expression.csv')####genes as row and cell_is as row
meta_info =  pd.read_csv('./Ery.integrated.meta.csv',index_col = 0)
meta_info
cell_embedding = pd.read_csv('./Ery.integrated.cell.embeddings.csv',index_col=0)
cell_embedding = np.array(cell_embedding)
adata.obs["cluster"] = meta_info["defined"]
adata.obs["n_counts"] = meta_info["nCount_RNA"]
adata.obs["samples"] = meta_info["samples"]
adata.obs["phase"] = meta_info["Phase"]
adata.var_names_make_unique()
adata
sc.pl.highest_expr_genes(adata, n_top=20)
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ['n_counts'])
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata)

adata.obsm["X_umap"] = cell_embedding

sc.pl.umap(adata, color=["cluster"])
adata.write(results_file)
###################################
#######Figure1-Sfigure1######
load("./Ery.integrated.h5ad")
adata = Ery.integrated.RData
cluster_col = ["#1F77B4","#FF7F0E","#8C564B","#D62728"]
stage_col = [("#461652","#FAE40E","#53C7DC"]
phase_col = [("#0E1674","#FF4A46","#4C9A25"]

pl.rcParams['pdf.fonttype'] = 42
sc.settings.set_figure_params(dpi=300)
sc.pl.umap(adata, color=["cluster"],palette= cluster_col,
         save="UMAP_clusters.pdf")
sc.pl.umap(adata, color=["stage"],palette= cluster_col,
        save="UMAP_stage.pdf")
from matplotlib import colors
cmap = colors.LinearSegmentedColormap.from_list('custom red', 
                                             [(0,    '#d8dcd6'),
                                              (0.2, '#d8dcd6'),
                                              (0.3, '#ffe36e'),
                                              (1,    '#e50000')], N=256)
                                              sc.settings.set_figure_params(dpi=300)
sc.pl.umap(adata, color=["GYPA","HBE1","HBG2","HBB"],color_map=cmap,vmax=4,
           ncols = 2,save="UMAP_marker_genes.pdf")
           
######################
#######velocity analysis#######
#######load the loom file of each sample##########
######use CS10_YS as example##########
CS10_loom = sc.read_loom("./loom_files/CS10_YS.loom")
CS10_loom.var_names_make_unique()
tmp = [re.sub(pattern= "CS10_YS:", repl= "YS_", string= x) for x in CS10_loom.obs_names]
tmp = [re.sub(pattern= "x", repl= "_1", string= x) for x in tmp]
CS10_loom.obs.index = tmp
#############combined loom files#######
loom_combine = CS10_loom.concatenate(CS11_loom, batch_key="stage")
tmp = [x[:-2] for x in loom_combine.obs_names]
loom_combine.obs.index = tmp
#################
set(adata.var_names).issubset(loom_combine.var_names)
vars_com = np.unique(adata.var_names.intersection(loom_combine.var_names))
loom_combine = loom_combine[adata.obs_names,:]
loom_combine = loom_combine[:, loom_combine.var_names.isin(vars_com)]
adata_vl = adata[:,adata.var_names.isin(vars_com)]
common_obs = adata_vl.obs_names.intersection(loom_combine.obs_names)
common_vars = adata_vl.var_names.intersection(loom_combine.var_names)
print(len(common_obs), len(common_vars))
loom_combine = loom_combine[:, adata_vl.var_names]
adata_vl = scv.utils.merge(adata = adata_vl, ldata = loom_combine)
scv.pp.filter_and_normalize(adata_vl, min_shared_counts=30, n_top_genes=1000)
scv.pp.moments(adata_vl, n_pcs=50, n_neighbors=50)
scv.tl.recover_dynamics(adata_vl)
scv.tl.velocity(adata_vl, mode='dynamical')
scv.tl.velocity_graph(adata_vl)
scv.settings.set_figure_params('scvelo', dpi_save=300)
scv.pl.velocity_embedding_stream(adata_vl, basis='umap',color="cluster",save = "Velocity_cluster.pdf",
                                smooth = 0.8, min_mass=1,legend_loc='right margin',figsize=[3.6,3.6],
                                palette = ["#1F77B4","#FF7F0E","#8C564B","#D62728"],
                                 outline_color = "white")
##########################
######PAGA analysis#######
sc.tl.paga(adata, groups='cluster',model="v1.0")
sc.settings.set_figure_params(dpi=300)
sc.pl.paga(adata,fontsize = 5, frameon = True,  edge_width_scale = 5, threshold = 0.02,  node_size_power=0.5,
           color = "cluster",save="cluster_PAGA.pdf")

#############################
#######Sfigure5##############
load("./pre_UCB_ery.h5ad")					
load("./Ery.integrated.h5ad")					
adata = pre_UCB_ery	
adata_ref = Ery.integrated.h5ad   

var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]

sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.ingest(adata, adata_ref, obs='cluster')

pd.DataFrame(adata.obs).to_csv("./pre_UCB_ery_predict_meta.csv")
pd.DataFrame(adata.obsm["X_umap"]).to_csv("./pre_UCB_ery_predict_X_UMAP.csv")

adata.write('/pre_UCB_ery_predicted.h5ad')

sc.pl.umap(adata, color=["cluster"],
           palette= cluster_col,save="pre_UCB_ery_cluster.pdf"
          )  