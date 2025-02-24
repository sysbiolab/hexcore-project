
from kedro.io import AbstractDataset
from pathlib import PurePosixPath
import fsspec
from typing import Any, Dict
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
        tags: list[str] = [],
        metadata: dict[str, Any] | None = None,) -> None:
        """Creates a new instance of ImageDataset to load / save image data for given filepath.

        Args:
            filepath: The location of the image file to load / save data.
        """
        # parse the path and protocol (e.g. file, http, s3, etc.)
        protocol, path = get_protocol_and_path(filepath)
        self._protocol = protocol
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol)
        self.metadata = metadata
        self.tags = tags
    
    def _load(self) -> Any:
        load_path = get_filepath_str(self._filepath, self._protocol)
        with self._fs.open(load_path) as file:
            data = file.read()
            for tag in self.tags:
                data = data.replace(f'#{tag}'.encode(), b'')
            return edn_format.loads_all(data, debug=True)
    
    def _save(self, data) -> Any:
        """Saves image data to the specified filepath"""
        return None

    def _describe(self) -> Dict[str, Any]:
        """Returns a dict that describes the attributes of the dataset."""
        return dict(filepath=self._filepath, protocol=self._protocol)