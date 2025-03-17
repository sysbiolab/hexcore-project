from kedro.pipeline import Pipeline, node
from hexcore_project.pipelines.cell_map.nodes import setup_visium_data, run_cell2location_map_model

def create_pipeline(sample: str, **kwargs):
    return Pipeline([
        node(setup_visium_data, 
             ["adata_ref_with_signatures", f"params:adata_viz_{sample}"], 
             [f"adata_viz_{sample}_tissue_filtered", f"inference_average_{sample}"], 
             name=f"setup_visium_data_{sample}"),
        node(run_cell2location_map_model,
             [f"adata_viz_{sample}_tissue_filtered", f"inference_average_{sample}", "params:model_cell_map"],
             [f"map_model_{sample}", f"adata_viz_{sample}_cell_abundance", f"map_model_elbo_loss_{sample}"],
             name=f"run_cell2location_{sample}_model")
    ] )
