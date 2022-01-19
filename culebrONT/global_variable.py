import sys
from pathlib import Path

DOCS = "https://culebront-pipeline.readthedocs.io/en/latest/"
GIT_URL = "https://github.com/SouthGreenPlatform/CulebrONT_pipeline"

# Hack for build docs with unspecified path install
args = str(sys.argv)
if "sphinx" in args:
    CULEBRONT_PATH = Path("/Path/to/culebrONT_install/")
else:
    CULEBRONT_PATH = Path(__file__).resolve().parent
CULEBRONT_SNAKEFILE = CULEBRONT_PATH.joinpath("snakefiles", "Snakefile")
CULEBRONT_MODE = CULEBRONT_PATH.joinpath(".mode.txt")
CULEBRONT_SCRIPTS = CULEBRONT_PATH.joinpath("snakemake_scripts")
CULEBRONT_PROFILE = CULEBRONT_PATH.joinpath("default_profile")
CULEBRONT_CONFIG_PATH = CULEBRONT_PATH.joinpath("install_files", "config.yaml")

CULEBRONT_TOOLS_PATH = CULEBRONT_PATH.joinpath("install_files", "tools_path.yaml")
CULEBRONT_USER_TOOLS_PATH = Path("~/.config/CulebrONT/tools_path.yaml").expanduser()
CULEBRONT_ARGS_TOOLS_PATH = Path("~/.config/CulebrONT/tools_path_args.yaml").expanduser()

CULEBRONT_CLUSTER_CONFIG = CULEBRONT_PROFILE.joinpath("cluster_config.yaml")
CULEBRONT_USER_CLUSTER_CONFIG = Path("~/.config/CulebrONT/cluster_config.yaml").expanduser()
CULEBRONT_ARGS_CLUSTER_CONFIG = Path("~/.config/CulebrONT/cluster_config_args.yaml").expanduser()


AVAIL_ASSEMBLY = ("CANU", "FLYE", "MINIASM", "RAVEN", "SMARTDENOVO", "SHASTA")
AVAIL_CORRECTION = ("NANOPOLISH", "MEDAKA", "PILON")
AVAIL_POLISHING = ("RACON")
AVAIL_QUALITY = ("BUSCO", "QUAST", "BLOBTOOLS", "ASSEMBLYTICS", "KAT", "FLAGSTATS")

ALLOW_FASTQ_EXT = (".fastq", ".fq", ".fq.gz", ".fastq.gz")

SINGULARITY_URL_FILES = [('https://itrop.ird.fr/culebront_utilities/singularity_build/Singularity.culebront_tools.sif',
              f'{CULEBRONT_PATH}/containers/Singularity.culebront_tools.sif'),
             ('https://itrop.ird.fr/culebront_utilities/singularity_build/Singularity.report.sif',
              f'{CULEBRONT_PATH}/containers/Singularity.report.sif')]
