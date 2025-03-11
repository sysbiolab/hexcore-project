from scipy.sparse import diags
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib

def plot_history(model, iter_start=0, iter_end=-1, ax=None):
        r"""Plot training history
        Parameters
        ----------
        iter_start
            omit initial iterations from the plot
        iter_end
            omit last iterations from the plot
        ax
            matplotlib axis
        """
        if ax is None:
            ax = plt.gca()
        if iter_end == -1:
            iter_end = len(model.history_["elbo_train"])

        ax.plot(
            np.array(model.history_["elbo_train"].index[iter_start:iter_end]),
            np.array(model.history_["elbo_train"].values.flatten())[iter_start:iter_end],
            label="train",
        )
        ax.legend()
        ax.set_xlim(0, len(model.history_["elbo_train"]))
        ax.set_xlabel("Training epochs")
        ax.set_ylabel("-ELBO loss")
        return plt

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

# ------------------------- deprecated
def desnormalize_log1p(adata: AnnData, layer_src: str, layer_trg: str) -> AnnData:
    """
    Normalização aplicada pelo autor: https://github.com/bioinformatics-inca/MyeloidDerivedCells_in_TME/blob/822912bfd6487b86cca875324c0ab231405d7423/scRNA-seq/02_Integration/Integration.py#L33C1-L44C2
    Aplica a transformação inversa de log1p (exp(x) - 1) para recuperar valores aproximados das contagens originais.

    Parâmetros:
    - adata: AnnData contendo a matriz de expressão gênica normalizada.

    Retorna:
    - Um novo objeto AnnData com a matriz desnormalizada.
    """
    adata = adata.copy()
    
    # Armazenar a matriz original normalizada em uma camada
    adata.layers["author_expression"] = adata.X.copy()

    X_log = adata.X
    # Reverter a transformação log1p mantendo esparsidade
    X_norm = X_log.expm1()

    # Estimar os fatores de escala por célula (cada linha representa uma célula)
    scaling_factors = adata.obs["nCount_RNA"].values / 1e4

    # Criar uma matriz diagonal esparsa com os fatores de escala
    scaling_matrix = diags(scaling_factors)
    X_counts = scaling_matrix @ X_norm
    adata.X = X_counts
    
    return adata

def plot_genes_expression_distribution(adata_ref: AnnData, layer_src: str, layer_trg: str) -> plt:
    """
    Plota a distribuição de expressão gênica para diferentes layers do anndata.

    Parâmetros:
    - adata_ref: AnnData contendo a matriz de expressão gênica original.
    - layer_src: Nome da camada contendo a matriz de expressão gênica original.
    - layer_trg: Nome da camada contendo a matriz de expressão gênica desnormalizada.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Obtém os valores de expressão gênica sem converter para array denso
    gene_expression_values_src = adata_ref.layers[layer_src].data
    gene_expression_values_target = adata_ref.layers[layer_trg].data

    # Cria o histograma para adata_ref
    sns.histplot(gene_expression_values_src, bins=100, kde=True, ax=axes[0])
    axes[0].set_xlabel("Nível de expressão gênica")
    axes[0].set_ylabel("Frequência")
    axes[0].set_title("Distribuição de todas as expressões gênicas (Normalizada pelo autor)")
    axes[0].set_yscale("log")  # Escala logarítmica no eixo Y para melhor visualização

    # Cria o histograma para andata_ref_disnorm
    sns.histplot(gene_expression_values_target, bins=100, kde=True, ax=axes[1])
    axes[1].set_xlabel("Nível de expressão gênica")
    axes[1].set_ylabel("Frequência")
    axes[1].set_title("Distribution of all gene expressions (Contagens originais)")
    axes[1].set_yscale("log")  # Escala logarítmica no eixo Y para melhor visualização

    plt.tight_layout()
    return plt
