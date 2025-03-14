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
    cell_map = cell_map_pipeline.create_pipeline()

    return {"__default__": Pipeline([cell_signature, cell_map])}
