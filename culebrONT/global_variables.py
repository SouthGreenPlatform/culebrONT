from pathlib import Path

DOCS = "https://culebront-pipeline.readthedocs.io/en/latest/"
GIT_URL = "https://github.com/SouthGreenPlatform/culebrONT"

INSTALL_PATH = Path(__file__).resolve().parent
SINGULARITY_URL_FILES = [('https://itrop.ird.fr/culebront_utilities/singularity_build/Singularity.culebront_tools.sif',
                          f'{INSTALL_PATH}/containers/Singularity.culebront_tools.sif'),
                         ('https://itrop.ird.fr/culebront_utilities/singularity_build/Singularity.report.sif',
                          f'{INSTALL_PATH}/containers/Singularity.report.sif')
                         ]

DATATEST_URL_FILES = ("Data-Xoo-sub.zip", "https://itrop.ird.fr/culebront_utilities/Data-Xoo-sub.zip")


AVAIL_ASSEMBLY = ("CANU", "FLYE", "MINIASM", "RAVEN", "SMARTDENOVO", "SHASTA")
AVAIL_CORRECTION = ("NANOPOLISH", "MEDAKA", "PILON")
AVAIL_POLISHING = ("RACON")
AVAIL_QUALITY = ("BUSCO", "QUAST", "BLOBTOOLS", "ASSEMBLYTICS", "KAT", "FLAGSTATS", "MERQURY")

ALLOW_FASTQ_EXT = (".fastq", ".fq", ".fq.gz", ".fastq.gz")

