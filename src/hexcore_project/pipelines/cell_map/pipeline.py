from kedro.pipeline import Pipeline, node
from hexcore_project.pipelines.cell_map.nodes import setup_visium_data


sample = 'normal'

def create_pipeline(**kwargs):
    return Pipeline([
        node(setup_visium_data, 
             ["adata_ref_with_signatures", f"{sample}"], 
             ["adata_viz_normal_tissue_filtered", "inference_average"], 
             name=f"setup_visium_data_{sample}"),
    ] )