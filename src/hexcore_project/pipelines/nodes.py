from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
import cell2location
import numpy as np

n_samples = 1000  # Número de células que você quer selecionar


import matplotlib  
matplotlib.use('Agg')  # Backend não interativo

def setup_adata_ref_signature(multi_tissue_tumor_microenvironment_atlas, n_samples_subset):
    
    # rename genes to ENSEMBL ID for correct matching between single cell and spatial data
    multi_tissue_tumor_microenvironment_atlas.var['SYMBOL'] = multi_tissue_tumor_microenvironment_atlas.var['feature_name']

    del multi_tissue_tumor_microenvironment_atlas.raw

    adata_ref = multi_tissue_tumor_microenvironment_atlas.copy()

    # Filtra todos os subtipos presentes no tecido de mama e depois disso seleciona de todo o dataset esses subtipos
    adata_ref = adata_ref[adata_ref.obs['author_cell_type'].isin(adata_ref[adata_ref.obs['tissue'] == 'breast'].obs['author_cell_type'])]

    # Filtra as células pelo numero de amostras
    if n_samples_subset > 0:
        indices = np.random.choice(adata_ref.n_obs, n_samples_subset, replace=False)  # Seleciona índices aleatórios
        adata_ref = adata_ref[indices, :]  

    # Desnormaliza a matriz de expressão gênica
    adata_ref_desnorm = desnormalize_log1p(adata_ref).copy()

    return [adata_ref_desnorm, plot_genes_expression_distribution(adata_ref, adata_ref_desnorm)]

def filter_genes_for_signature_model(adata_ref, adata_ref_params):
    # fonte: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/spatial/cell2location_lymph_node_spatial_tutorial.html
    # Filtra os genes com base nos parâmetros

    selected = filter_genes(adata_ref, 
                            cell_count_cutoff=adata_ref_params['cell_count_cutoff'], 
                            cell_percentage_cutoff2=adata_ref_params['cell_percentage_cutoff2'], 
                            nonz_mean_cutoff=adata_ref_params['nonz_mean_cutoff'])
    
    # Aplica os filtros de genes caso a função filter_genes não esteja funcionando
    # Filtro hardcoded por conta de a função filter_genes não estar rolando

    gene_filter = (
        ((adata_ref.X > 0).sum(axis=0) / adata_ref.n_obs > 0.05) &  # Expressão > 0 em pelo menos 5% das células
        (np.array(adata_ref.X.mean(axis=0)).flatten() > 1.1) &       # Média de expressão > 1.1
        ((adata_ref.X > 0).sum(axis=0) / adata_ref.n_obs > 0.0005)   # Expressão > 0 em pelo menos 0.05% das células
    )

    # adata_ref = adata_ref[:, gene_filter].copy()

    # Filtra o AnnData com base nos genes selecionados
    adata_ref = adata_ref[:, selected].copy()

    return [adata_ref, plot_filter_genes(adata_ref, 
                                         adata_ref_params['cell_count_cutoff'], 
                                         adata_ref_params['cell_percentage_cutoff2'], 
                                         adata_ref_params['nonz_mean_cutoff'])]

def create_signatures_model(adata_ref, params):
    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
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

    return adata_ref, model_cell_type_signatures, model_cell_type_signatures.plot_history(20)