from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
import cell2location
import numpy as np
from hexcore_project.pipelines import utils

import matplotlib  
matplotlib.use('Agg')  # Backend não interativo

def setup_adata_ref_signature(multi_tissue_tumor_microenvironment_atlas, n_samples_subset):
    
    # rename genes to ENSEMBL ID for correct matching between single cell and spatial data
    multi_tissue_tumor_microenvironment_atlas.var['SYMBOL'] = multi_tissue_tumor_microenvironment_atlas.var['feature_name']

    adata_ref = multi_tissue_tumor_microenvironment_atlas.copy()

    # Filtra todos os subtipos presentes no tecido de mama e depois disso seleciona de todo o dataset esses subtipos
    adata_ref = adata_ref[adata_ref.obs['author_cell_type'].isin(adata_ref[adata_ref.obs['tissue'] == 'breast'].obs['author_cell_type'])]

    # Filtra as células pelo numero de amostras definido
    if n_samples_subset > 0:
        indices = np.random.choice(adata_ref.n_obs, n_samples_subset, replace=False)  # Seleciona índices aleatórios
        adata_ref = adata_ref[indices, :]  

    # Convertendo adata_ref.raw.X de float32 para int e movendo para X, mantendo a versão original em 'autor_expression'.
    # Criar uma cópia de raw.X em uma nova camada antes da conversão
    adata_ref.layers["normalized_expression"] = adata_ref.X.copy() 

    # Converter para int32
    adata_ref.layers["counts"] = adata_ref.raw.X.astype(np.int32)

    return [adata_ref, utils.plot_genes_expression_distribution(adata_ref, "normalized_expression", "counts")]

def filter_genes_for_signature_model(adata_ref, adata_ref_params):
    # Fonte: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/spatial/cell2location_lymph_node_spatial_tutorial.html
    # Filtra os genes com base nos parâmetros

    selected = filter_genes(adata_ref, 
                            cell_count_cutoff=adata_ref_params['cell_count_cutoff'], 
                            cell_percentage_cutoff2=adata_ref_params['cell_percentage_cutoff2'], 
                            nonz_mean_cutoff=adata_ref_params['nonz_mean_cutoff'])

    # Filtra o AnnData com base nos genes selecionados
    adata_ref = adata_ref[:, selected].copy()

    return [adata_ref, utils.plot_filter_genes(adata_ref, 
                                         adata_ref_params['cell_count_cutoff'], 
                                         adata_ref_params['cell_percentage_cutoff2'], 
                                         adata_ref_params['nonz_mean_cutoff'])]

def create_signatures_model(adata_ref, params):
    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                            layer='counts', # Use layer recovered from raw data 
                            # 10X reaction / sample / batch
                            batch_key=params['batch_key'],
                            # cell type, covariate used for constructing signatures
                            labels_key=params['labels_key'],
                            # multiplicative technical effects (platform, 3' vs 5', donor effect)
                            categorical_covariate_keys=[params['categorical_covariate_keys']]
                           )
    model_cell_type_signatures = RegressionModel(adata_ref)
    model_cell_type_signatures.train(max_epochs=params['epochs'])

    adata_ref = model_cell_type_signatures.export_posterior(
        adata_ref, sample_kwargs={'num_samples': params['num_samples'], 'batch_size': params['batch_size']}
    )

    return adata_ref, model_cell_type_signatures, utils.plot_history(model=model_cell_type_signatures, iter_start=20)