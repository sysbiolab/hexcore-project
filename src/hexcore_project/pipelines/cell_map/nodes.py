import cell2location
import numpy as np
from hexcore_project.pipelines import utils
import scanpy as sc

import matplotlib  
matplotlib.use('Agg')  # Backend não interativo

def setup_visium_data(adata_ref_with_signatures, sample_name):
    sample_list = {
        'normal': 'https://datasets.cellxgene.cziscience.com/3574c8a3-2bf8-4aa6-97b8-3c4a7751757a.h5ad',
        'tumor': 'https://datasets.cellxgene.cziscience.com/6cffa6f2-fd1b-49ab-841e-5fcbeab669df.h5ad'
    }
    
    inf_aver = adata_ref_with_signatures.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' \
                                           for i in adata_ref_with_signatures.uns['mod']['factor_names']]]
    inf_aver.columns = adata_ref_with_signatures.uns['mod']['factor_names']

    adata_viz = sc.read(
        f'./data/raw/adata_viz/{sample_name}.h5ad',
        backup_url=sample_list[sample_name]
    )

    adata_viz.obs['sample'] = list(adata_viz.uns['spatial'].keys())[0]
    adata_viz.var['SYMBOL'] = adata_viz.var['feature_name']

    # find mitochondria-encoded (MT) genes
    adata_viz.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_viz.var['SYMBOL']]

    # remove MT genes for spatial mapping (keeping their counts in the object)
    adata_viz.obsm['MT'] = adata_viz[:, adata_viz.var['MT_gene'].values].X.toarray()
    adata_viz = adata_viz[:, ~adata_viz.var['MT_gene'].values]
    
    intersect = np.intersect1d(adata_viz.var_names, inf_aver.index)
    
    adata_viz = adata_viz[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    return adata_viz, inf_aver