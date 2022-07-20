#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
from snakemake.logging import logger
from snakemake.utils import validate
import re
from .global_variables import *
from .snakeWrapper import *


class CulebrONT(SnakeWrapper):
    """
    to read file config
    """

    def __init__(self, workflow, config, ):
        super().__init__(workflow, config)
        # workflow is available only in __init__
        # print("\n".join(list(workflow.__dict__.keys())))
        # print(workflow.__dict__)

        # Initialisation of CulebrONT attributes
        self.assembly_tools_activated = []
        self.polishing_tools_activated = []
        self.correction_tools_activated = []
        self.quality_tools_activated = []
        self.quality_step = []
        self.last_steps_list = []
        self.pipeline_stop = None

        self.fastq_files_list = []
        self.fastq_files_ext = []
        self.fastq_gzip = None

        self.illumina_files_list = []
        self.illumina_files_ext = []
        self.illumina_gzip = None

        self.R1 = ""
        self.R2 = ""

        self.add_circular_name = None
        self.TMP = None
        self.TCM = None
        self.TAG = None

        self.draft_to_correction = None
        self.draft_to_correction_index_fai = None
        self.draft_to_correction_index_mmi = None

        self.nb_racon_rounds = None
        self.nb_pilon_rounds = None

        self.__check_config_dic()
        self.__cleaning_for_rerun()

        # write corrected config to output:
        self.write_config(f"{self.get_config_value('DATA', 'OUTPUT')}/config_corrected.yaml")

        # using schemas to check mandatory value from yaml format
        try:
            validate(self.config, INSTALL_PATH.joinpath("schemas/config.schema.yaml").resolve().as_posix())
        except Exception as e:
            raise ValueError(
                f"{e}\n\nCONFIG FILE CHECKING STRUCTURE FAIL : you need to verify {self.path_config} KEYS:VALUES: {str(e)[30:76]}\n")

    def __split_illumina(self):
        R1 = []
        R2 = []
        for fastq in self.illumina_files_list:
            if '_R1' in fastq:
                R1.append(fastq)
            elif '_R2' in fastq:
                R2.append(fastq)
        return R1, R2

    def __cleaning_for_rerun(self):
        """Cleaning Report and dags if we rerun the snakemake a second time
        TODO: improve check if final report run maybe check html file"""
        if Path(self.get_config_value('DATA', 'OUTPUT')).joinpath("FINAL_REPORT/CulebrONT_report.html").exists():
            logger.warning(
                f"\nWARNING: If you want to rerun culebront a second time please delete the {self.get_config_value('DATA', 'OUTPUT')}FINAL_REPORT and REPORT repertories for each sample\n")

    def __check_tools_config(self, tool, mandatory=()):
        """Check if path is a file if not empty
        :return absolute path file"""
        tool_OK = False
        if tool in ["FLAGSTATS"]:
            tool = "SAMTOOLS"

        def check_singularity(sif):
            level1 = "SINGULARITY"
            path_file = self.tools_config[level1][sif]
            if re.findall("shub://SouthGreenPlatform/CulebrONT_pipeline", path_file, flags=re.IGNORECASE):
                raise ValueError(
                    f'CONFIG FILE CHECKING FAIL : shub download is obsolete')
            else:
                path_file = eval(f"f'{path_file}'")
                path = Path(path_file).resolve().as_posix()
                if path and path_file:
                    if not Path(path).exists() or not Path(path).is_file():
                        raise FileNotFoundError(
                            f'CONFIG FILE CHECKING FAIL : please check tools_config.yaml in the {level1} section, {tool} file "{path}" {"does not exist" if not Path(path).exists() else "is not a valid file"}')
                    else:
                        self.tools_config[level1][tool] = path

        # If only singularity tools
        if not self.use_env_modules and self.use_singularity:
            check_singularity("REPORT")     # check singularity image for Report
            check_singularity("TOOLS")      # check singularity image for conda (contains a lot of conda envs)
            tool_OK = True

        # If only envmodule
        if self.use_env_modules and not self.use_singularity:
            envmodule_key = self.tools_config["ENVMODULE"][tool]
            if not envmodule_key:
                raise ValueError(
                    f'CONFIG FILE CHECKING FAIL : please check tools_config.yaml in the "ENVMODULE" section, {tool} is empty')
            tool_OK = True

        # If envmodule and singularity
        if self.use_env_modules and self.use_singularity:
            raise ValueError(
                f"CONFIG FILE CHECKING FAIL : Use env-module or singularity but don't mix them")

        if len(mandatory) > 0 and not tool_OK:
            raise FileNotFoundError(
                f'CONFIG FILE CHECKING FAIL : please check tools_config.yaml in the  {tool} params, please append Singularity or module load, is mandatory for tool: {" ".join(mandatory)}')

    def __build_tools_activated(self, level1, allow, mandatory=False):
        tools_activate = []
        for tool, activated in self.config[level1].items():
            if tool in allow:
                boolean_activated = var_2_bool(level1, tool, activated)
                if boolean_activated:
                    tools_activate.append(tool)
                    self.set_config_value(level1=level1, level2=tool, value=boolean_activated)
                    self.__check_tools_config(tool, [tool])
            else:
                raise ValueError(f'CONFIG FILE CHECKING FAIL : {level1} {tool} not allow on CulebrONT"')
        if len(tools_activate) == 0 and mandatory:
            raise ValueError(f"CONFIG FILE CHECKING FAIL : you need to set True for at least one {level1} from {allow}")
        return tools_activate

    def __build_quality_step_list(self, only_last=False):
        last_steps_list = []
        suffix = ""
        if bool(self.config["FIXSTART"]):
            suffix = "sFIX"
        if self.correction_tools_activated:
            step = "CORRECTION"
            for corrector in self.correction_tools_activated:
                last_steps_list.append(f"{corrector}{suffix if step in self.pipeline_stop else ''}")
            if only_last:
                return last_steps_list

        if self.polishing_tools_activated:
            step = "POLISHING"
            for polisher in self.polishing_tools_activated:
                last_steps_list.append(f"{polisher}{suffix if step in self.pipeline_stop else ''}")
            if only_last:
                return last_steps_list
        if self.assembly_tools_activated:
            step = "ASSEMBLY"
            last_steps_list.append(f"ASSEMBLY{suffix if step in self.pipeline_stop else ''}")
            if only_last:
                return last_steps_list
        return last_steps_list

    def __get_last_step(self):
        if self.correction_tools_activated:
            return "CORRECTION"
        if self.polishing_tools_activated:
            return "POLISHING"
        if self.assembly_tools_activated:
            return "ASSEMBLY"

    def __check_config_dic(self):
        """Configuration file checking"""
        # check tools activation
        self.assembly_tools_activated = self.__build_tools_activated("ASSEMBLY", AVAIL_ASSEMBLY, True)
        self.polishing_tools_activated = self.__build_tools_activated("POLISHING", AVAIL_POLISHING)
        self.correction_tools_activated = self.__build_tools_activated("CORRECTION", AVAIL_CORRECTION)
        self.quality_tools_activated = self.__build_tools_activated("QUALITY", AVAIL_QUALITY)
        self.pipeline_stop = self.__get_last_step()

        # check mandatory directory
        self._check_dir_or_string(level1="DATA", level2="OUTPUT")
        self._check_dir_or_string(level1="DATA", level2="FASTQ", mandatory=self.assembly_tools_activated)

        # check if fastq file for assembly
        self.fastq_files_list, fastq_files_list_ext = get_files_ext(self.get_config_value('DATA', 'FASTQ'),
                                                                    ALLOW_FASTQ_EXT)
        if not self.fastq_files_list:
            raise ValueError(
                f"CONFIG FILE CHECKING FAIL : you need to append at least on fastq with extension on {ALLOW_FASTQ_EXT}")
        # check if all fastq have the same extension
        if len(fastq_files_list_ext) > 1:
            raise ValueError(
                f"CONFIG FILE CHECKING FAIL : Please use only the same format for assembly FASTQ data, not: {fastq_files_list_ext}")
        else:
            self.fastq_files_ext = fastq_files_list_ext[0]
        # check if fastq are gzip
        if "gz" in self.fastq_files_ext:
            self.fastq_gzip = True

        # CHECK ILLUMINA
        need_illumina = [var_2_bool("QUALITY", "KAT", self.get_config_value("QUALITY", "KAT")),
                         var_2_bool("QUALITY", "MERQURY", self.get_config_value("QUALITY", "MERQURY")),
                         var_2_bool("CORRECTION", "PILON", self.get_config_value("CORRECTION", "PILON")),
                         var_2_bool("QUALITY", "FLAGSTATS", self.get_config_value("QUALITY", "FLAGSTATS"))]
        if True in need_illumina:
            self._check_dir_or_string(level1="DATA", level2="ILLUMINA", mandatory=["KAT", "PILON", "FLAGSTATS"])
            self.illumina_files_list, illumina_files_list_ext = get_files_ext(self.get_config_value('DATA', 'ILLUMINA'),
                                                                              ALLOW_FASTQ_EXT)

            if not self.illumina_files_list:
                raise ValueError(
                    f"CONFIG FILE CHECKING FAIL : you need to append at least on fastq illumina with extension on {ALLOW_FASTQ_EXT}")
            # check if all fastq have the same extension
            if len(illumina_files_list_ext) > 1:
                raise ValueError(
                    f"CONFIG FILE CHECKING FAIL : Please use only the same format for ILLUMINA FASTQ data, not: {illumina_files_list_ext}")
            else:
                self.illumina_files_ext = illumina_files_list_ext[0]
            # check if fastq are gzip
            if "gz" in self.illumina_files_ext:
                self.illumina_gzip = True

            self.R1, self.R2 = self.__split_illumina()

        # build quality_step and last_steps_list; function test if FIXSTART to append only to last step
        self.last_steps_list = self.__build_quality_step_list(only_last=True)
        self.quality_step = self.__build_quality_step_list(only_last=False)

        # check files if QUAST
        if bool(self.config["QUALITY"]["QUAST"]):
            self._check_file_or_string(level1="DATA", level2="REF")
        genome_pb = convert_genome_size(self.get_config_value('DATA', 'GENOME_SIZE'))
        self.set_config_value(level1="params", level2="QUAST", level3="GENOME_SIZE_PB", value=genome_pb)

        # check files if MERQURY
        if bool(self.config["QUALITY"]["MERQURY"]):
            self._check_file_or_string(level1="DATA", level2="REF")
        genome_pb = convert_genome_size(self.get_config_value('DATA', 'GENOME_SIZE'))
        self.set_config_value(level1="params", level2="MERQURY", value={"GENOME_SIZE_PB":genome_pb})

        # check files if MAUVE
        if bool(self.config["MSA"]["MAUVE"]):
            self._check_file_or_string(level1="DATA", level2="REF", mandatory=["MAUVE"])
            # Make sure running Mauve makes sense
            if len(self.assembly_tools_activated) < 2 and len(self.quality_step) < 2:
                raise ValueError(
                    "CONFIG FILE CHECKING ERROR : MAUVE is irrelevant if you have a single final assembly as your config file implies (need more than one assembler and/or a correction step). Mauve will not be run !! ")

        # check BUSCO database if activate
        if bool(self.config["QUALITY"]["BUSCO"]):
            self._check_dir_or_string(level1='params', level2='BUSCO', level3='DATABASE', mandatory=["BUSCO"], check_string=True)

        # check DIAMOND database if activate
        if bool(self.config["QUALITY"]["BLOBTOOLS"]):
            self._check_file_or_string(level1='params', level2='DIAMOND', level3='DATABASE', mandatory=["BLOBTOOLS", 'DIAMOND'])

        # check if NANOPOLISH activate, if true compare fastq and fast5 files
        if "NANOPOLISH" in self.correction_tools_activated:
            self._check_dir_or_string(level1="DATA", level2="FAST5", mandatory=["NANOPOLISH"])
            fast5_dir_list = get_dir(self.config['DATA']['FAST5'])
            fastq_files_list, _ = get_files_ext(self.config['DATA']['FASTQ'], ALLOW_FASTQ_EXT, add_ext=False)

            if set(fastq_files_list) - set(fast5_dir_list):
                raise ValueError(
                    f"CONFIG FILE CHECKING ERROR : You don't have a fast5 repository for each of your fastq file (they should have the same name). This can raise a problem if you choose to use Nanopolish. Please check your data.:\n\t- fast5_dir_list:{fast5_dir_list}\n\t- fastq_files_list: {fastq_files_list}\n\n")

        # check Medaka config if activate
        if bool(self.config['CORRECTION']['MEDAKA']):
            if not bool(self.config['params']['MEDAKA']['MEDAKA_TRAIN_WITH_REF']):
                self._check_file_or_string(level1='params', level2='MEDAKA', level3='MEDAKA_MODEL_PATH', mandatory=["MEDAKA"], check_string=True)
            else:
                self._check_file_or_string(level1='DATA', level2='REF', mandatory=["MEDAKA"])

        # Check racon round
        if bool(self.config['POLISHING']['RACON']) and not (0 < int(self.config['params']['RACON']['RACON_ROUNDS']) < 10):
            raise ValueError(f"CONFIG FILE CHECKING ERROR : You have activated RACON, but RACON_ROUNDS is invalid, 0 < RACON_ROUNDS={self.config['params']['RACON']['RACON_ROUNDS']} < 10 . \n")

        ##############################
        # check workflow compatibility
        if bool(self.config['QUALITY']['ASSEMBLYTICS']):
            self._check_file_or_string(level1='DATA', level2='REF', mandatory=["ASSEMBLYTICS"])

        if not bool(self.config['POLISHING']['RACON']) and bool(self.config['ASSEMBLY']['MINIASM']):
            logger.warning(
                f"\nWARNING: RACON is automatically launched (2 rounds by default) for minipolish if MINIASM is activated !! . \n")

        # check size of genome
        if int(convert_genome_size(self.get_config_value('DATA', 'GENOME_SIZE'))) >= 50000000 and self.get_config_value('MSA', 'MAUVE'):
            logger.warning(
                f"WARNING: CONFIG FILE CHECKING WARNING : MAUVE if fixed to FALSE because genome size  >= 50 mb !! \n")
            self.set_config_value('MSA', 'MAUVE', False)

        # Make sure running fixstart makes sense
        if not bool(self.config['CIRCULAR']) and bool(self.config['FIXSTART']):
            raise ValueError(
                f"CONFIG FILE CHECKING ERROR : FIXSTART is irrelevant if you have not activated CIRCULAR. FIXSTART will not be run !! \n")

        # if you want to run mauve fixstart has to be activated
        if (bool(self.config['CIRCULAR']) and bool(self.config['MSA']['MAUVE'])) and not bool(self.config['FIXSTART']):
            raise ValueError(
                f"CONFIG FILE CHECKING ERROR : FIXSTART must be activated for MAUVE if CIRCULAR is TRUE on config file.  !! \n")

        # check for BLOBTOOLS
        if bool(self.config["QUALITY"]["BLOBTOOLS"]):
            self._check_file_or_string(level1='params', level2='BLOBTOOLS', level3='NAMES', mandatory=["BLOBTOOLS"])
            self._check_file_or_string(level1='params', level2='BLOBTOOLS', level3='NODES', mandatory=["BLOBTOOLS"])

        #############################################
        # Build variables name for files
        if bool(self.config['CIRCULAR']):
            self.add_circular_name = "CIRCULARISED"
            self.TMP = "TMP"
            self.TCM = "FALSE"
        else:
            self.add_circular_name = ""
            self.TMP = ""
            self.TCM = ""

        if not bool(self.config['POLISHING']['RACON']) and bool(self.config['CIRCULAR']):
            self.TAG = "TMPTAG"
        else:
            self.TAG = ""

        # ############################### DEF ####################################
        self.nb_racon_rounds = '2' if not bool(self.config['POLISHING']['RACON']) else str(
            self.config['params']['RACON']['RACON_ROUNDS'])
        self.nb_pilon_rounds = '1' if not bool(self.config['CORRECTION']['PILON']) else str(
            self.config['params']['PILON']['PILON_ROUNDS'])
