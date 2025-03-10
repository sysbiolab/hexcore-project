from kedro.pipeline import Pipeline, node
from hexcore_project.pipelines.nodes import setup_adata_ref_signature, create_signatures_model


def create_pipeline(**kwargs):
    return Pipeline([
        node(setup_adata_ref_signature, 
             ["multi_tissue_tumor_microenvironment_atlas", "params:adata_ref_filter"], 
             "adata_ref", 
             name="setup_reference_signature_dataset"),
        node(
            create_signatures_model,
            ["adata_ref", "params:model_cell_type_signatures"],
            ["adata_ref_with_signatures", "model_cell_type_signatures"],
            name="create_signatures_model"
        )
    ] )