"""Project pipelines."""

from kedro.pipeline import Pipeline

from hexcore_project.pipelines.cell_signature import pipeline as cell_signature_pipeline
from hexcore_project.pipelines.cell_map import pipeline as cell_map_pipeline

def register_pipelines() -> dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from pipeline names to ``Pipeline`` objects.
    """
    cell_signature = cell_signature_pipeline.create_pipeline()
    cell_map_normal = cell_map_pipeline.create_pipeline(sample='normal')
    cell_map_turmor = cell_map_pipeline.create_pipeline(sample='tumor')
    
    return {
        "__default__": Pipeline([cell_signature, cell_map_normal, cell_map_turmor]),
        "cell_signature": Pipeline([cell_signature]),
        "cell_map": Pipeline([cell_map_normal, cell_map_turmor]),
    }
