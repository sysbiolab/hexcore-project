import cell2location
import numpy as np
from hexcore_project.pipelines import utils
import scanpy as sc

import matplotlib  
matplotlib.use('Agg')  # Backend não interativo

def setup_visium_data(adata_ref_with_signatures, sample):    
    inf_aver = adata_ref_with_signatures.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' \
                                           for i in adata_ref_with_signatures.uns['mod']['factor_names']]]
    inf_aver.columns = adata_ref_with_signatures.uns['mod']['factor_names']

    adata_viz = sc.read(
        f"./data/raw/adata_viz/{sample['name']}.h5ad",
        backup_url=sample['url']
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

def run_cell2location_map_model(adata_viz, inf_aver, params):
    cell2location.models.Cell2location.setup_anndata(adata=adata_viz, batch_key=params["batch_key"])
    # create and train the model
    model_cell2location = cell2location.models.Cell2location(
        adata_viz, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=params['N_cells_per_location'],
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=params['detection_alpha']
    )

    model_cell2location.train(max_epochs=params['max_epochs'],
        # train using full data (batch_size=None)
        batch_size=params['batch_size'],
        # use all data points in training because
        # we need to estimate cell abundance at all locations
        train_size=params['train_size']
    )

    adata_viz = model_cell2location.export_posterior(
        adata_viz, sample_kwargs={'num_samples': params['num_samples'], 'batch_size': model_cell2location.adata.n_obs}
    )

    return model_cell2location, adata_viz, utils.plot_history(model=model_cell2location, iter_start=1000)