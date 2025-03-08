
from kedro.io import AbstractDataset
from pathlib import PurePosixPath
import fsspec
from typing import Any, Dict
import scanpy as sc
import logging

from kedro.io.core import (
    Version,
    get_filepath_str,
    get_protocol_and_path,
)

logger = logging.getLogger(__name__)

class AnnDataset(AbstractDataset):
    def __init__(self,         
        *,
        filepath: str,
        url: str,
        tags: list[str] = [],
        metadata: dict[str, Any] | None = None,) -> None:
        """Creates a new instance of AnnDataset to load / save data for given filepath or url.

        Args:
            filepath: The location of the h5 file to load / save data.
            url: The url to download the h5 file from.
        """
        # parse the path and protocol (e.g. file, http, s3, etc.)
        protocol, path = get_protocol_and_path(filepath)
        self._protocol = protocol
        self.url = url
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol)
        self.metadata = metadata
        self.tags = tags
    
    def _load(self) -> Any:
        print(self._filepath)
        return sc.read(
            self._filepath,
            backup_url=self.url
        )
    
    def _save(self, data) -> Any:
        """Saves image data to the specified filepath"""
        adata = data.copy()
        return adata.write(self._filepath)

    def _describe(self) -> Dict[str, Any]:
        """Returns a dict that describes the attributes of the dataset."""
        return dict(filepath=self._filepath, protocol=self._protocol)