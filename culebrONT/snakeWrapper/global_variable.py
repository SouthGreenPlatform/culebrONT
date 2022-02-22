import sys
from pathlib import Path
from .usefull_function import package_name
from ..global_variables import *

# Hack for build docs with unspecified path install
args = str(sys.argv)
if "sphinx" in args:
    INSTALL_PATH = Path(f"/Path/to/{package_name()}_install/")
else:
    INSTALL_PATH = Path(__file__).resolve().parent.parent
SNAKEFILE = INSTALL_PATH.joinpath("snakefiles", "Snakefile")
INSTALL_MODE = INSTALL_PATH.joinpath(".mode.txt")
SNAKEMAKE_SCRIPTS = INSTALL_PATH.joinpath("snakemake_scripts")
DEFAULT_PROFILE = INSTALL_PATH.joinpath("default_profile")
GIT_CONFIG_PATH = INSTALL_PATH.joinpath("install_files", "config.yaml")

GIT_TOOLS_PATH = INSTALL_PATH.joinpath("install_files", "tools_path.yaml")
USER_TOOLS_PATH = Path(f"~/.config/{package_name()}/tools_path.yaml").expanduser()
ARGS_TOOLS_PATH = Path(f"~/.config/{package_name()}/tools_path_args.yaml").expanduser()

GIT_CLUSTER_CONFIG = DEFAULT_PROFILE.joinpath("cluster_config.yaml")
USER_CLUSTER_CONFIG = Path(f"~/.config/{package_name()}/cluster_config.yaml").expanduser()
ARGS_CLUSTER_CONFIG = Path(f"~/.config/{package_name()}/cluster_config_args.yaml").expanduser()
