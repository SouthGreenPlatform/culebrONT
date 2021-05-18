#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Sebastien Ravel

##################################################
# Modules
##################################################
import snakemake
import sys
import os
from pathlib import Path
import argparse
from datetime import datetime
from pprint import pprint as pp

from AdditionalScripts.toolbox import existant_file, welcome_args, get_last_version, get_version
###################################

try:
    CULEBRONT_PATH = Path(os.environ["CULEBRONT"])
except KeyError:
    CULEBRONT_PATH = Path(__file__).resolve().parent
SNAKEFILE = CULEBRONT_PATH.joinpath("Snakefile")
CLUSTER_CONFIG = CULEBRONT_PATH.joinpath("cluster_config.yaml")
SLURM_STATUS = CULEBRONT_PATH.joinpath("slurm_status.py")
SLURM_WRAPPER = CULEBRONT_PATH.joinpath("slurm_wrapper.py")
version_CulebrONT = get_version(CULEBRONT_PATH)

##################################################
# build_parser function use with sphinxcontrib.autoprogram
##################################################
def build_parser():
    # Change the description HERE
    description_tools = f"""
    Welcome to CulebrONT  !
    Created on November 2019
    version: {version_CulebrONT}
    @author: Julie Orjuela (IRD), Aurore Comte (IRD), Sebastien Ravel (CIRAD), Florian Charriat (CIRAD), Bao Tram Vi (IRD), François Sabot (IRD) and Sebastien Cunnac (IRD)
    @email: julie.orjuela@ird.fr, aurore@comte.ird.fr

    #     .-+.
    #   `odNNh
    #   +NNd:
    #  .Nh.   ---:`
    #  -Nm`  ````./
    #   oNm/ ```-o-
    #    .yNmo/+/.     .oooooo.               oooo             .o8                  .oooooo.   ooooo      ooo ooooooooooooo
    #    `-+ydNmo.    d8P'  `Y8b              `888            "888                 d8P'  `Y8b  `888b.     `8' 8'   888   `8
    #  `/s+../oNNN:  888          oooo  oooo   888   .ooooo.   888oooo.  oooo d8b 888      888  8 `88b.    8       888
    #  +s/     `hNm  888          `888  `888   888  d88' `88b  d88' `88b `888""8P 888      888  8   `88b.  8       888
    #  os:  ////hNN` 888           888   888   888  888ooo888  888   888  888     888      888  8     `88b.8       888
    #  -so- .//sNNs  `88b    ooo   888   888   888  888    .o  888   888  888     `88b    d88'  8       `888       888
    #   `/oo:.yNd/    `Y8bood8P'   `V88V"V8P' o888o `Y8bod8P'  `Y8bod8P' d888b     `Y8bood8P'  o8o        `8      o888o
    #     -yooo+.
    #   `yNs.`-/oo:
    #   dNo` ....+s+
    #   :shmdddhhy+:

    Please cite our github https://github.com/SouthGreenPlatform/CulebrONT_pipeline
    Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)
    and GPLv3 Intellectual property belongs to IRD and authors.
    """
    parser_mandatory = argparse.ArgumentParser(add_help=False)
    mandatory = parser_mandatory.add_argument_group('Input mandatory infos for running')

    # Write mandatory arguments HERE
    mandatory.add_argument('-c', '--config', metavar="path/to/file/config_file", type=existant_file, required=True,
                           dest='configfile', help='path to file params file')
    # END Write mandatory arguments HERE

    parser_other = argparse.ArgumentParser(
            parents=[parser_mandatory],
            add_help=False,
            prog=Path(__file__).name,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=description_tools,
            epilog=get_last_version(version_CulebrONT)
    )

    optional = parser_other.add_argument_group('Input infos not mandatory')
    optional.add_argument('-v', '--version', action='version', version=version_CulebrONT,
                          help=f'Use if you want to know which version of {Path(__file__).name} you are using')
    optional.add_argument('-h', '--help', action='help', help=f'show this help message and exit')
    optional.add_argument('-d', '--debug', action='store_true', help='enter verbose/debug mode')

    # Write optional arguments HERE
    optional.add_argument('-n', '--dry_run', action='store_true', help='run with dry run mode')
    # optional.add_argument('-p', '--profile', action='store_true', help='run profile instead of cluster')
    # optional.add_argument('-p', '--profile', metavar="profile", type=str, required=False, default = None,
                           # dest='profile', help='path to file params file')


    # END Write optional arguments HERE

    return parser_other


@welcome_args(version_CulebrONT, build_parser())
def main():
    prog_args = build_parser().parse_args()
    # Main Code HERE

    # run CulebrONT!!
    # if not prog_args.profile:
    status = snakemake.snakemake(snakefile=SNAKEFILE,
                                 lock=False,
                                 use_conda=True,
                                 use_singularity=True,
                                 cores=1,
                                 verbose=False,
                                 latency_wait=60000000,
                                 keepgoing=True,
                                 restart_times=1,
                                 force_incomplete=True,
                                 configfiles=[prog_args.configfile.as_posix()],
                                 cluster=f"python3 {CULEBRONT_PATH}/slurm_wrapper.py {prog_args.configfile.as_posix()} {CLUSTER_CONFIG.as_posix()}",
                                 cluster_config=CLUSTER_CONFIG.as_posix(),
                                 cluster_status=f"python3 {CULEBRONT_PATH}/slurm_status.py",
                                 dryrun=prog_args.dry_run
                                     )
    # elif prog_args.profile:
        # print("NOT AVAIL NOW !!!!!!!")
        # exit()
        # from snakemake import get_argument_parser, get_profile_file
        # argv=None
        # # reparse args while inferring config file from profile
        # parser2 = get_argument_parser(prog_args.profile)
        # args = parser2.parse_args(argv)

        # pp(parser2)
        # pp(str(args).split(","))
        # def adjust_path(f):
            # if os.path.exists(f) or os.path.isabs(f):
                # return f
            # else:
                # return get_profile_file(prog_args.profile, f, return_default=True)

        # # update file paths to be relative to the profile
        # # (if they do not exist relative to CWD)

        # if args.jobscript:
            # args.jobscript = adjust_path(args.jobscript)
        # if args.cluster:
            # args.cluster = adjust_path(args.cluster)
            # # args.cluster = Path("/shared/home/sravel/.config/snakemake/IFB_CulebrONT/config.yaml").as_posix()
        # if args.cluster_sync:
            # args.cluster_sync = adjust_path(args.cluster_sync)
        # if args.cluster_status:
            # args.cluster_status = adjust_path(args.cluster_status)
        # if args.report_stylesheet:
            # args.report_stylesheet = adjust_path(args.report_stylesheet)

        # print(Path(args.cluster).resolve())
        # # print(Path(args.cluster).resolve().open(mode='r').read())
        # # exit()
        # print(f"""\n\n
                  # args.jobscript         {args.jobscript}
                  # args.cluster           {args.cluster},
                  # args.cluster_sync      {args.cluster_sync}
                  # args.cluster_status    {args.cluster_status}
                  # args.report_stylesheet {args.report_stylesheet}\n\n""")
        # status = None
        # status = snakemake.snakemake(snakefile=SNAKEFILE,
                             # lock=False,
                             # use_conda=True,
                             # use_singularity=True,
                             # cores=1,
                             # verbose=False,
                             # configfiles=[prog_args.configfile.as_posix()],
                             # # jobscript=args.jobscript,
                             # cluster=args.jobscript,
                             # cluster_config=CLUSTER_CONFIG.as_posix(),
                             # cluster_status=args.cluster_status,
                             # dryrun=prog_args.dry_run
                             # )

    # # SLURM JOBS WITHOUT PROFILES
    # snakemake --nolock --use-conda --use-singularity --cores -p --verbose -s $CULEBRONT_PATH/Snakefile \
    # --latency-wait 60 --keep-going --restart-times 1 --rerun-incomplete  \
    # --configfile $config \
    # --cluster "python3 $CULEBRONT_PATH/slurm_wrapper.py $config $CULEBRONT_PATH/cluster_config.yaml" \
    # --cluster-config $CULEBRONT_PATH/cluster_config.yaml \
    # --cluster-status "python3 $CULEBRONT_PATH/slurm_status.py"
    # USING PROFILES
    # snakemake --nolock --use-singularity --use-conda --cores -p --verbose -s Snakefile --configfile config.yaml \
    # --latency-wait 60 --keep-going --restart-times 1 --rerun-incomplete --cluster-config cluster_config.yaml --profile slurm-culebrONT

    if status: # translate "success" into shell exit code of 0
       return 0
    return 1


if __name__ == '__main__':
    main()
