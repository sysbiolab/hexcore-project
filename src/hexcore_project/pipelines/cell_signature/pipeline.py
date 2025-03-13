from kedro.pipeline import Pipeline, node
from hexcore_project.pipelines.cell_signature.nodes import setup_adata_ref_signature, create_signatures_model, filter_genes_for_signature_model


def create_pipeline(**kwargs):
    return Pipeline([
        node(setup_adata_ref_signature, 
             ["multi_tissue_tumor_microenvironment_atlas", "params:adata_ref_filter.n_samples_subset"], 
             ["adata_ref", "adata_ref_gene_expression_distributions"], 
             name="setup_reference_signature_dataset"),
        node(filter_genes_for_signature_model,
             ["adata_ref", "params:adata_ref_filter"],
             ["adata_ref_filtered", "adata_ref_filtered_genes"],
             name="filter_genes_for_signature_model"),
        node(
            create_signatures_model,
            ["adata_ref_filtered", "params:model_cell_type_signatures"],
            ["adata_ref_with_signatures", "model_cell_type_signatures", "signature_model_elbo_loss"],
            name="create_signature_model"
        )
    ] )