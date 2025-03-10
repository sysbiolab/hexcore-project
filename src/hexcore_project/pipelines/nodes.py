from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
import cell2location
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse

n_samples = 1000  # Número de células que você quer selecionar


import matplotlib  
matplotlib.use('Agg')  # Backend não interativo

def plot_filter_genes(adata, cell_count_cutoff, cell_percentage_cutoff2, nonz_mean_cutoff):
    
    adata.var["n_cells"] = np.array((adata.X > 0).sum(0)).flatten()
    adata.var["nonz_mean"] = np.array(adata.X.sum(0)).flatten() / adata.var["n_cells"]

    cell_count_cutoff = np.log10(cell_count_cutoff)
    cell_count_cutoff2 = np.log10(adata.shape[0] * cell_percentage_cutoff2)
    nonz_mean_cutoff = np.log10(nonz_mean_cutoff)

    gene_selection = (np.array(np.log10(adata.var["n_cells"]) > cell_count_cutoff2)) | (
        np.array(np.log10(adata.var["n_cells"]) > cell_count_cutoff)
        & np.array(np.log10(adata.var["nonz_mean"]) > nonz_mean_cutoff)
    )
    gene_selection = adata.var_names[gene_selection]
    adata_shape = adata[:, gene_selection].shape

    fig, ax = plt.subplots()
    ax.hist2d(
        np.log10(adata.var["nonz_mean"]),
        np.log10(adata.var["n_cells"]),
        bins=100,
        norm=matplotlib.colors.LogNorm(),
        range=[[0, 0.5], [1, 4.5]],
    )
    ax.axvspan(0, nonz_mean_cutoff, ymin=0.0, ymax=(cell_count_cutoff2 - 1) / 3.5, color="darkorange", alpha=0.3)
    ax.axvspan(
        nonz_mean_cutoff,
        np.max(np.log10(adata.var["nonz_mean"])),
        ymin=0.0,
        ymax=(cell_count_cutoff - 1) / 3.5,
        color="darkorange",
        alpha=0.3,
    )
    plt.vlines(nonz_mean_cutoff, cell_count_cutoff, cell_count_cutoff2, color="darkorange")
    plt.hlines(cell_count_cutoff, nonz_mean_cutoff, 1, color="darkorange")
    plt.hlines(cell_count_cutoff2, 0, nonz_mean_cutoff, color="darkorange")
    plt.xlabel("Mean non-zero expression level of gene (log)")
    plt.ylabel("Number of cells expressing gene (log)")
    plt.title(f"Gene filter: {adata_shape[0]} cells x {adata_shape[1]} genes")
    
    return plt

def desnormalize_log1p(adata: AnnData) -> AnnData:
    """
    Aplica a transformação inversa de log1p (exp(x) - 1) para recuperar valores aproximados das contagens originais.

    Parâmetros:
    - adata: AnnData contendo a matriz de expressão gênica normalizada com log1p.

    Retorna:
    - Um novo objeto AnnData com a matriz desnormalizada.
    """
    adata = adata.copy()  # Evita modificar o objeto original
    
    if sp.issparse(adata.X):
        adata.X = sp.csr_matrix(np.expm1(adata.X.toarray()))  # Converte para denso antes da transformação
    else:
        adata.X = np.expm1(adata.X)
    
    return adata

def plot_genes_expression_distribution(adata_ref: AnnData, andata_ref_disnorm: AnnData) -> plt:
    """
    Plota a distribuição de expressão gênica para um conjunto de genes para dois objetos AnnData.

    Parâmetros:
    - adata_ref: AnnData contendo a matriz de expressão gênica original.
    - andata_ref_disnorm: AnnData contendo a matriz de expressão gênica desnormalizada.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Converte para array denso, se for uma matriz esparsa
    if scipy.sparse.issparse(adata_ref.X):
        gene_expression_values_ref = adata_ref.X.toarray().flatten()
    else:
        gene_expression_values_ref = adata_ref.X.flatten()

    if scipy.sparse.issparse(andata_ref_disnorm.X):
        gene_expression_values_ref_disnorm = andata_ref_disnorm.X.toarray().flatten()
    else:
        gene_expression_values_ref_disnorm = andata_ref_disnorm.X.flatten()

    # Cria o histograma para adata_ref
    sns.histplot(gene_expression_values_ref, bins=100, kde=True, ax=axes[0])
    axes[0].set_xlabel("Nível de expressão gênica")
    axes[0].set_ylabel("Frequência")
    axes[0].set_title("Distribuição de todas as expressões gênicas (Dataset Original)")
    axes[0].set_yscale("log")  # Escala logarítmica no eixo Y para melhor visualização

    # Cria o histograma para andata_ref_disnorm
    sns.histplot(gene_expression_values_ref_disnorm, bins=100, kde=True, ax=axes[1])
    axes[1].set_xlabel("Nível de expressão gênica")
    axes[1].set_ylabel("Frequência")
    axes[1].set_title("Distribuição de todas as expressões gênicas (Desnormalizada)")
    axes[1].set_yscale("log")  # Escala logarítmica no eixo Y para melhor visualização

    plt.tight_layout()
    return plt

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