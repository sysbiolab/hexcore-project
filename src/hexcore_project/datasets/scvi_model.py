
from kedro.io import AbstractDataset
from pathlib import PurePosixPath
import fsspec
from typing import Any, Dict
import scanpy as sc
import logging
import os
from scvi.model.base import BaseModelClass

from kedro.io.core import (
    Version,
    get_filepath_str,
    get_protocol_and_path,
)

logger = logging.getLogger(__name__)

class ScviModel(AbstractDataset):
    def __init__(self,         
        *,
        filepath: str | None = None,
        tags: list[str] = [],
        metadata: dict[str, Any] | None = None,) -> None:
        """Creates a new instance of ScviModel to load / save model for given filepath.

        Args:
            filepath: The location of the h5 file to load / save data.
        """
        # parse the path and protocol (e.g. file, http, s3, etc.)
        protocol, path = get_protocol_and_path(filepath)
        self._protocol = protocol
        self.path = path
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol)
        self.metadata = metadata
        self.tags = tags
    
    def _load(self) -> BaseModelClass:
        return 
    
    def _save(self, model: BaseModelClass):
        """Saves image data to the specified filepath"""
        # TODO: versioning and add metadata to the model
        os.makedirs(os.path.dirname(self.path), exist_ok=True)

        return model.save(f"{self.path}", overwrite=True)

    def _describe(self) -> Dict[str, Any]:
        """Returns a dict that describes the attributes of the dataset."""
        return dict(filepath=self._filepath, protocol=self._protocol)