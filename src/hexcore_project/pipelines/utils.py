import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib


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
        adata.X.data = np.expm1(adata.X.data)  # Aplica a transformação diretamente nos dados esparsos
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

    # Obtém os valores de expressão gênica sem converter para array denso
    gene_expression_values_ref = adata_ref.X.data if sp.issparse(adata_ref.X) else adata_ref.X.flatten()
    gene_expression_values_ref_disnorm = andata_ref_disnorm.X.data if sp.issparse(andata_ref_disnorm.X) else andata_ref_disnorm.X.flatten()

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
