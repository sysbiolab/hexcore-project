"""Project settings. There is no need to edit this file unless you want to change values
from the Kedro defaults. For further information, including these default values, see
https://docs.kedro.org/en/stable/kedro_project_setup/settings.html."""

# Instantiated project hooks.
# For example, after creating a hooks.py and defining a ProjectHooks class there, do
# from pandas_viz.hooks import ProjectHooks

# Hooks are executed in a Last-In-First-Out (LIFO) order.
# HOOKS = (ProjectHooks(),)

# Installed plugins for which to disable hook auto-registration.
# DISABLE_HOOKS_FOR_PLUGINS = ("kedro-viz",)

# Class that manages storing KedroSession data.
# from kedro.framework.session.store import BaseSessionStore
# SESSION_STORE_CLASS = BaseSessionStore
# Keyword arguments to pass to the `SESSION_STORE_CLASS` constructor.
# SESSION_STORE_ARGS = {
#     "path": "./sessions"
# }

# Directory that holds configuration.
# CONF_SOURCE = "conf"


# Class that manages how configuration is loaded.
from kedro.config import OmegaConfigLoader  # noqa: E402
from kedro_viz.integrations.kedro.sqlite_store import SQLiteStore
from pathlib import Path
import os
import socket

import scvi
import torch


# scvi.settings.num_threads = 31
# scvi.settings.batch_size = 512
# scvi.settings.progress_bar_style = "rich"
# scvi.settings.dl_num_workers = 16

# torch.set_float32_matmul_precision('high')

# torch.cuda.empty_cache()  # Limpa cache de memória antiga
# torch.cuda.memory_allocated()  # Verifica memória usada
# torch.cuda.memory_reserved()  # Verifica memória reservada
# torch.cuda.set_per_process_memory_fraction(0.95, device=0)


hostname = socket.gethostname()

os.environ["KEDRO_SQLITE_STORE_USERNAME"] = hostname

CONFIG_LOADER_CLASS = OmegaConfigLoader
# Keyword arguments to pass to the `CONFIG_LOADER_CLASS` constructor.
CONFIG_LOADER_ARGS = {
      "base_env": "base",
      "default_run_env": "local",
#       "config_patterns": {
#           "spark" : ["spark*/"],
#           "parameters": ["parameters*", "parameters*/**", "**/parameters*"],
#       }
}


SESSION_STORE_CLASS = SQLiteStore
SESSION_STORE_ARGS = {
    "path": str(Path(__file__).parents[2] / "data/kedro-viz"),
}

# Class that manages Kedro's library components.
# from kedro.framework.context import KedroContext
# CONTEXT_CLASS = KedroContext

# Class that manages the Data Catalog.
# from kedro.io import DataCatalog
# DATA_CATALOG_CLASS = DataCatalog
