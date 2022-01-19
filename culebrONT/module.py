#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from snakemake.logging import logger
from snakemake.utils import validate
from snakemake.io import load_configfile
import yaml
from collections import OrderedDict
import pprint
import re
import culebrONT
import click
from culebrONT.global_variable import *
from culebrONT.usefull_function import get_install_mode


def get_dir(path):
    """List of directory included on folder"""
    return [elm.name for elm in Path(path).glob("*") if elm.is_dir()]


def get_files_ext(path, extensions, add_ext=True):
    """List of files with specify extension included on folder

    Arguments:
        path (str): a path to folder
        extensions (list or tuple): a list or tuple of extension like (".py")
        add_ext (bool): if True (default), file have extension

    Returns:
        :class:`list`: List of files name with or without extension , with specify extension include on folder
        :class:`list`: List of  all extension found

    Examples:
        >>> version = get_version("/path/to/install/culebront")
        >>> print(version)
            1.3.0
     """
    if not (extensions, (list, tuple)) or not extensions:
        raise ValueError(f'ERROR CulebrONT: "extensions" must be a list or tuple not "{type(extensions)}"')
    tmp_all_files = []
    all_files = []
    files_ext = []
    for ext in extensions:
        tmp_all_files.extend(Path(path).glob(f"**/*{ext}"))

    for elm in tmp_all_files:
        ext = "".join(elm.suffixes)
        if ext not in files_ext:
            files_ext.append(ext)
        if add_ext:
            all_files.append(elm.as_posix())
        else:
            if len(elm.suffixes) > 1:

                all_files.append(Path(elm.stem).stem)
            else:
                all_files.append(elm.stem)
    return all_files, files_ext


def convert_genome_size(size):
    mult = dict(K=10 ** 3, M=10 ** 6, G=10 ** 9, T=10 ** 12, N=1)
    search = re.search(r'^(\d+\.?\d*)\s*(.*)$', size)
    if not search or len(search.groups()) != 2:
        raise ValueError(
            f"CONFIG FILE CHECKING FAIL : not able to convert genome size please only use int value with{' '.join(mult.keys())} upper or lower, N or empty is bp size")
    else:
        value, unit = search.groups()
        if not unit:
            return int(value)
        elif unit and unit.upper() not in mult.keys():
            raise ValueError(
                f"CONFIG FILE CHECKING FAIL : '{unit}' unit value not allow or not able to convert genome size please only use int value with{' '.join(mult.keys())} upper or lower, N or empty is bp size")
        else:
            return int(float(value) * mult[unit.upper()])


class CulebrONT(object):
    """
    to read file config
    """

    def __init__(self, workflow, config):

        # workflow is available only in __init__
        self.snakefile = CULEBRONT_SNAKEFILE
        self.tools_config = None

        if not workflow.overwrite_configfiles:
            raise ValueError("ERROR CulebrONT: You need to use --configfile option to snakemake command line")
        else:
            self.path_config = workflow.overwrite_configfiles[0]
       
        if not workflow.overwrite_clusterconfig and get_install_mode() == "cluster":
            self.cluster_config = load_configfile(CULEBRONT_PROFILE.joinpath("cluster_config.yaml"))
        elif not workflow.overwrite_clusterconfig and get_install_mode() == "local":
            self.cluster_config = None
        else:
            self.cluster_config = workflow.overwrite_clusterconfig

        self.load_tool_configfile()
        # print("\n".join(list(workflow.__dict__.keys())))
        # print(workflow.__dict__)

        # --- Verification Configuration Files --- #
        self.config = config
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

        self.use_env_modules = workflow.use_env_modules
        self.use_conda = workflow.use_conda
        self.use_singularity = workflow.use_singularity

        self.__check_config_dic()
        self.__cleaning_for_rerun()

        # With good config to output:
        self.write_config(f"{self.config['DATA']['OUTPUT']}/config_corrected.yaml")

        # using schemas to check mandatory value from yaml format
        try:
            validate(self.config, culebrONT.CULEBRONT_PATH.joinpath("schemas/config.schema.yaml").resolve().as_posix())
        except Exception as e:
            raise ValueError(
                f"{e}\n\nCONFIG FILE CHECKING STRUCTURE FAIL : you need to verify {self.path_config} KEYS:VALUES: {str(e)[30:76]}\n")

    def load_tool_configfile(self):
        if CULEBRONT_USER_TOOLS_PATH.exists() and not CULEBRONT_ARGS_TOOLS_PATH.exists():
            self.tools_config = load_configfile(CULEBRONT_USER_TOOLS_PATH)
        elif CULEBRONT_ARGS_TOOLS_PATH.exists():
            self.tools_config = load_configfile(CULEBRONT_ARGS_TOOLS_PATH)
            CULEBRONT_ARGS_TOOLS_PATH.unlink()
        else:
            self.tools_config = load_configfile(CULEBRONT_TOOLS_PATH)

    def __split_illumina(self):
        R1 = []
        R2 = []
        for fastq in self.illumina_files_list:
            if '_R1' in fastq:
                R1.append(fastq)
            elif '_R2' in fastq:
                R2.append(fastq)
        return R1, R2

    def get_config_value(self, section, key, subsection=None, type_value=None):
        # TODO add type_value to load a good type
        if subsection:
            return self.config[section][subsection][key]
        else:
            return self.config[section][key]

    def set_config_value(self, section, key, value, subsection=None):
        if subsection:
            self.config[section][subsection][key] = value
        else:
            self.config[section][key] = value

    def write_config(self, path):
        p = Path(path).parent
        p.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as config_open:
            config_open.write(self.export_use_yaml)

    @property
    def export_use_yaml(self):
        """Use to print a dump config.yaml with corrected parameters"""

        def represent_dictionary_order(self, dict_data):
            return self.represent_mapping('tag:yaml.org,2002:map', dict_data.items())

        def setup_yaml():
            yaml.add_representer(OrderedDict, represent_dictionary_order)

        setup_yaml()
        return yaml.dump(self.config, default_flow_style=False, sort_keys=False, indent=4)

    def __cleaning_for_rerun(self):
        """Cleaning Report and dags if we rerun the snakemake a second time
        TODO: improve check if final report run maybe check html file"""
        if Path(self.config["DATA"]["OUTPUT"]).joinpath("FINAL_REPORT/CulebrONT_report.html").exists():
            logger.warning(
                f"\nWARNING: If you want to rerun culebront a second time please delete the {self.config['DATA']['OUTPUT']}FINAL_REPORT and REPORT repertories for each sample\n")

    def __check_dir(self, section, key, mandatory=[], subsection=None):
        """Check if path is a directory if not empty
            resolve path on config

        Arguments:
            section (str): the first level on config.yaml
            key (str): the final level on config.yaml
            mandatory (list tuple): a list or tuple with tools want mandatory info
            subsection (str): the second level on config.yaml (ie 3 level)
        Returns:
            :class:`list`: List of files name with or without extension , with specify extension include on folder
            :class:`list`: List of all extension found
        Raises:
            NotADirectoryError: If config.yaml data `path` does not exist.
        """
        path_value = self.get_config_value(section=section, key=key, subsection=subsection)
        if path_value:
            path = Path(path_value).resolve().as_posix() + "/"
            if (not Path(path).exists() or not Path(path).is_dir()) and key not in ["OUTPUT"]:
                raise NotADirectoryError(
                    f'CONFIG FILE CHECKING FAIL : in the "{section}" section, {f"subsection {subsection}" if subsection else ""}, {key} directory "{path}" {"does not exist" if not Path(path).exists() else "is not a valid directory"}')
            else:
                self.set_config_value(section, key, path, subsection)
        elif len(mandatory) > 0:
            raise NotADirectoryError(
                f'CONFIG FILE CHECKING FAIL : in the "{section}" section, {f"subsection {subsection}" if subsection else ""}, {key} directory "{path_value}" {"does not exist" if not Path(path_value).exists() else "is not a valid directory"} but is mandatory for tool: {" ".join(mandatory)}')

    def __check_file(self, section, key, mandatory=[], subsection=None):
        """Check if path is a file if not empty
        :return absolute path file"""
        path_value = self.get_config_value(section=section, key=key, subsection=subsection)
        path = Path(path_value).resolve().as_posix()
        if path_value != "":
            if not Path(path).exists() or not Path(path).is_file():
                raise FileNotFoundError(
                    f'CONFIG FILE CHECKING FAIL : in the {section} section, {f"subsection {subsection}" if subsection else ""},{key} file "{path}" {"does not exist" if not Path(path).exists() else "is not a valid file"}')
            else:
                self.set_config_value(section, key, path, subsection)
        elif len(mandatory) > 0:
            raise FileNotFoundError(
                f'CONFIG FILE CHECKING FAIL : in the "{section}" section, {f"subsection {subsection}" if subsection else ""},{key} file "{path_value}" {"does not exist" if not Path(path_value).exists() else "is not a valid file"} but is mandatory for tool: {" ".join(mandatory)}')

    def __check_dir_or_string (self, section, key, mandatory=[], subsection=None, check_string=False):
        """ similar to check_dir"""
        path_value = self.get_config_value(section=section, key=key, subsection=subsection)
        if path_value:
            path = Path(path_value).resolve().as_posix() + "/"
            # it is a path
            if path_value != "" and "/" in path_value:
                if (not Path(path).exists() or not Path(path).is_dir()) and key not in ["OUTPUT"]:
                    raise NotADirectoryError(
                        f'CONFIG FILE CHECKING FAIL : in the "{section}" section, {f"subsection {subsection}" if subsection else ""}, {key} directory "{path}" {"does not exist" if not Path(path).exists() else "is not a valid directory"}')
                else:
                    self.set_config_value(section, key, path, subsection)
            # it is not a path
            elif path_value != "" and not "/" in path_value and check_string:
                self.set_config_value(section, key, path_value, subsection)
            # it is empty
            elif path_value == "" and not "/" in path_value and check_string:
                raise ValueError(
                    f'CONFIG FILE CHECKING FAIL : in the {section} section, {f"subsection {subsection}" if subsection else ""},{key} Value "{path_value}" is empty')
        elif len(mandatory) > 0:
            raise NotADirectoryError(
                f'CONFIG FILE CHECKING FAIL : in the "{section}" section, {f"subsection {subsection}" if subsection else ""}, {key} directory "{path_value}" {"does not exist" if not Path(path_value).exists() else "is not a valid directory"} but is mandatory for tool: {" ".join(mandatory)}')

    def __check_file_or_string(self, section, key, mandatory=[], subsection=None, check_string=False):
        """Check if path is a file if not empty
        :return absolute path file"""
        path_value = self.get_config_value(section=section, key=key, subsection=subsection)
        path = Path(path_value).resolve().as_posix()
        # it is a path
        if path_value != "" and "/" in path_value:
            if not Path(path).exists() or not Path(path).is_file():
                raise FileNotFoundError(
                    f'CONFIG FILE CHECKING FAIL : in the {section} section, {f"subsection {subsection}" if subsection else ""},{key} file "{path}" {"does not exist" if not Path(path).exists() else "is not a valid file"}')
            else:
                self.set_config_value(section, key, path, subsection)
        # it is not a path
        elif path_value != "" and not "/" in path_value and check_string:
            self.set_config_value(section, key, path_value, subsection)
        # it is empty
        elif path_value == "" and not "/" in path_value and check_string:
            raise ValueError(
                f'CONFIG FILE CHECKING FAIL : in the {section} section, {f"subsection {subsection}" if subsection else ""},{key} Value "{path_value}" is empty')
        elif len(mandatory) > 0:
            raise FileNotFoundError(
                f'CONFIG FILE CHECKING FAIL : in the "{section}" section, {f"subsection {subsection}" if subsection else ""}:{key} file "{path_value}" {"does not exist" if not Path(path_value).exists() else "is not a valid file"} but is mandatory for tool: {" ".join(mandatory)}')


    def __check_tools_config(self, tool, mandatory=[]):
        """Check if path is a file if not empty
        :return absolute path file"""
        tool_OK = False
        if tool in ["FLAGSTATS"]:
            tool = "SAMTOOLS"
        def check_singularity(tool):
            section = "SINGULARITY"
            path_file = self.tools_config[section][tool]
            if re.findall("shub://SouthGreenPlatform/CulebrONT_pipeline", path_file, flags=re.IGNORECASE):
                raise ValueError(
                    f'CONFIG FILE CHECKING FAIL : shub download is obsolete')
            else :
                path_file = eval(f"f'{path_file}'")
                path = Path(path_file).resolve().as_posix()
                if path and path_file:
                    if not Path(path).exists() or not Path(path).is_file():
                        raise FileNotFoundError(
                            f'CONFIG FILE CHECKING FAIL : please check tools_config.yaml in the {section} section, {tool} file "{path}" {"does not exist" if not Path(path).exists() else "is not a valid file"}')
                    else:
                        self.tools_config[section][tool] = path

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

    @property
    def string_to_dag(self):
        """ return command line for rule graph """
        return f"""snakemake -s {self.snakefile} {'--use-singularity' if self.use_singularity else ''} {'--use-envmodules' if self.use_env_modules else ''}  --rulegraph"""        

    def __var_2_bool(self, key, tool, to_convert):
        """convert to boolean"""
        if isinstance(type(to_convert), bool):
            return to_convert
        elif f"{to_convert}".lower() in ("yes", "true", "t"):
            return True
        elif f"{to_convert}".lower() in ("no", "false", "f"):
            return False
        else:
            raise TypeError(
                f'CONFIG FILE CHECKING FAIL : in the "{key}" section, "{tool}" key: "{to_convert}" is not a valide boolean')

    def __build_tools_activated(self, key, allow, mandatory=False):
        tools_activate = []
        for tool, activated in self.config[key].items():
            if tool in allow:
                boolean_activated = self.__var_2_bool(key, tool, activated)
                if boolean_activated:
                    tools_activate.append(tool)
                    self.config[key][tool] = boolean_activated
                    self.__check_tools_config(tool, [tool])
            else:
                raise ValueError(f'CONFIG FILE CHECKING FAIL : {key} {tool} not allow on CulebrONT"')
        if len(tools_activate) == 0 and mandatory:
            raise ValueError(f"CONFIG FILE CHECKING FAIL : you need to set True for at least one {key} from {allow}")
        return tools_activate

    def __build_quality_step_list(self, only_last=False):
        last_steps_list = []
        suffix = ""
        if bool(self.config["FIXSTART"]):
            suffix = "_STARTFIXED"
        if self.correction_tools_activated:
            step = "CORRECTION"
            for corrector in self.correction_tools_activated:
                last_steps_list.append(f"STEP_{step}_{corrector}{suffix if step in self.pipeline_stop else ''}")
            if only_last:
                return last_steps_list

        if self.polishing_tools_activated:
            step = "POLISHING"
            for polisher in self.polishing_tools_activated:
                last_steps_list.append(f"STEP_{step}_{polisher}{suffix if step in self.pipeline_stop else ''}")
            if only_last:
                return last_steps_list
        if self.assembly_tools_activated:
            step = "ASSEMBLY"
            last_steps_list.append(f"STEP_{step}{suffix if step in self.pipeline_stop else ''}")
            # for assembler in self.assembly_tools_activated:
            # last_steps_list.append(f"STEP_{step}_{assembler}{suffix}" )
            if only_last: return last_steps_list
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
        self.__check_dir(section="DATA", key="OUTPUT")
        self.__check_dir(section="DATA", key="FASTQ", mandatory=self.assembly_tools_activated)

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

        ##### CHECK KAT
        if bool(self.config["QUALITY"]["KAT"]) or bool(self.config["CORRECTION"]["PILON"]) or bool(
                self.config["QUALITY"]['FLAGSTATS']):
            self.__check_dir(section="DATA", key="ILLUMINA", mandatory=["KAT", "PILON", "FLAGSTATS"])
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
            self.__check_file(section="DATA", key="REF")
        genome_pb = convert_genome_size(self.get_config_value('DATA', 'GENOME_SIZE'))
        self.set_config_value(section="params", subsection="QUAST", key="GENOME_SIZE_PB", value=genome_pb)

        # check files if MAUVE
        if bool(self.config["MSA"]["MAUVE"]):
            self.__check_file(section="DATA", key="REF", mandatory=["MAUVE"])
            # Make sure running Mauve makes sense
            if len(self.assembly_tools_activated) < 2 and len(self.quality_step) < 2:
                raise ValueError(
                    "CONFIG FILE CHECKING ERROR : MAUVE is irrelevant if you have a single final assembly as your config file implies (need more than one assembler and/or a correction step). Mauve will not be run !! ")

        # check BUSCO database if activate
        if bool(self.config["QUALITY"]["BUSCO"]):
            #self.__check_dir(section='params', subsection='BUSCO', key='DATABASE', mandatory=["BUSCO"])
            self.__check_dir_or_string(section='params', subsection='BUSCO', key='DATABASE', mandatory=["BUSCO"], check_string=True)

        # check DIAMOND database if activate
        if bool(self.config["QUALITY"]["BLOBTOOLS"]):
            self.__check_file(section='params', subsection='DIAMOND', key='DATABASE',
                              mandatory=["BLOBTOOLS", 'DIAMOND'])

        # check if NANOPOLISH activate, if true compare fastq and fast5 files
        if "NANOPOLISH" in self.correction_tools_activated:
            self.__check_dir(section="DATA", key="FAST5", mandatory=["NANOPOLISH"])
            fast5_files_list = get_dir(self.config['DATA']['FAST5'])
            fastq_files_list, _ = get_files_ext(self.config['DATA']['FASTQ'], ALLOW_FASTQ_EXT, add_ext=False)

            if set(fastq_files_list) - set(fast5_files_list):
                raise ValueError(
                    f"CONFIG FILE CHECKING ERROR : You don't have a fast5 repository for each of your fastq file (they should have the same name). This can raise a problem if you choose to use Nanopolish. Please check your data.:\n\t- fast5_files_list:{fast5_files_list}\n\t- fastq_files_list: {fastq_files_list}\n\n")

        # check Medaka config if activate
        if bool(self.config['CORRECTION']['MEDAKA']):
            # check singularity image
            if not bool(self.config['params']['MEDAKA']['MEDAKA_TRAIN_WITH_REF']):
                #self.__check_file(section='params', subsection='MEDAKA', key='MEDAKA_MODEL_PATH', mandatory=["MEDAKA"])
                self.__check_file_or_string(section='params', subsection='MEDAKA', key='MEDAKA_MODEL_PATH',
                                            mandatory=["MEDAKA"], check_string=True)
            else:
                self.__check_file(section='DATA', key='REF', mandatory=["MEDAKA"])

        # Check racon round
        if bool(self.config['POLISHING']['RACON']) and not (
                0 < int(self.config['params']['RACON']['RACON_ROUNDS']) < 10):
            raise ValueError(
                f"CONFIG FILE CHECKING ERROR : You have activated RACON, but RACON_ROUNDS is invalid, 0 < RACON_ROUNDS={self.config['params']['RACON']['RACON_ROUNDS']} < 10 . \n")

        ##############################
        # check workflow compatibility
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

        # if you want run mauve fixstart has to be activated
        if bool(self.config['CIRCULAR']) and not bool(self.config['FIXSTART']):
            raise ValueError(
                f"CONFIG FILE CHECKING ERROR : FIXSTART must be activated if CIRCULAR is TRUE on config file.  !! \n")

        # check for BLOBTOOLS
        if bool(self.config["QUALITY"]["BLOBTOOLS"]):
            self.__check_file(section='params', subsection='BLOBTOOLS', key='NAMES', mandatory=["BLOBTOOLS"])
            self.__check_file(section='params', subsection='BLOBTOOLS', key='NODES', mandatory=["BLOBTOOLS"])

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

    def __repr__(self):
        return f"{self.__class__}({pprint.pprint(self.__dict__)})"
