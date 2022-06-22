from pathlib import Path
from .global_variable import *
from .usefull_function import get_install_mode
from snakemake.io import load_configfile
import yaml
from collections import OrderedDict
import pprint


class SnakeWrapper(object):
    """ test generic wrapper class"""
    def __init__(self, workflow, config):
        # workflow is available only in __init__
        # print("\n".join(list(workflow.__dict__.keys())))
        # print(workflow.__dict__)
        self.snakefile = SNAKEFILE
        self.tools_config = None
        self.path_config = workflow.overwrite_configfiles[0]
        self.config = config

        self.use_env_modules = workflow.use_env_modules
        self.use_conda = workflow.use_conda
        self.use_singularity = workflow.use_singularity
        # test if cluster_config.yaml pass to snakemake
        if not workflow.overwrite_clusterconfig and get_install_mode() == "cluster":
            self.cluster_config = load_configfile(DEFAULT_PROFILE.joinpath("cluster_config.yaml"))
        elif not workflow.overwrite_clusterconfig and get_install_mode() == "local":
            self.cluster_config = None
        else:
            self.cluster_config = workflow.overwrite_clusterconfig

        self.load_tool_configfile()

    def load_tool_configfile(self):
        """Test path of tools_path.yaml on default install path, home or argument"""
        if USER_TOOLS_PATH.exists() and not ARGS_TOOLS_PATH.exists():
            self.tools_config = load_configfile(USER_TOOLS_PATH)
        elif ARGS_TOOLS_PATH.exists():
            self.tools_config = load_configfile(ARGS_TOOLS_PATH)
            ARGS_TOOLS_PATH.unlink()
        else:
            self.tools_config = load_configfile(GIT_TOOLS_PATH)

    def get_config_value(self, level1, level2=None, level3=None):
        """get value on config_file"""
        # TODO add type_value to load a good type
        if level3:
            return self.config[level1][level2][level3]
        elif level2:
            return self.config[level1][level2]
        else:
            return self.config[level1]

    def set_config_value(self, level1, level2, value, level3=None):
        if level3:
            self.config[level1][level2][level3] = value
        else:
            self.config[level1][level2] = value

    def write_config(self, path):
        p = Path(path).parent
        p.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as config_open:
            config_open.write(self.export_use_yaml)

    def _check_dir_or_string(self, level1, level2, mandatory=(), level3=None, check_string=False):
        """ similar to check_dir"""
        path_value = self.get_config_value(level1=level1, level2=level2, level3=level3)
        if path_value:
            path = Path(path_value).resolve().as_posix() + "/"
            # it is a path
            if path_value != "" and "/" in path_value:
                if (not Path(path).exists() or not Path(path).is_dir()) and level2 not in ["OUTPUT"]:
                    raise NotADirectoryError(
                        f'CONFIG FILE CHECKING FAIL : in section:{level1}, {f"subsection:{level2} directory:{level3}" if level3 else f"directory:{level2}"}, "{path}" {"does not exist" if not Path(path).exists() else "is not a valid directory"}')
                else:
                    self.set_config_value(level1=level1, level2=level2, level3=level3, value=path)
            # it is not a path
            elif path_value and "/" not in path_value and check_string:
                self.set_config_value(level1=level1, level2=level2, level3=level3, value=path_value)
            # it is empty
            elif path_value and "/" not in path_value and check_string:
                raise ValueError(
                    f'CONFIG FILE CHECKING FAIL : in section:{level1}, {f"subsection:{level2} value:{level3}" if level3 else f"value:{level2}"}, "{path_value}" is empty')
        elif len(mandatory) > 0:
            raise NotADirectoryError(
                f'CONFIG FILE CHECKING FAIL : in section:{level1}, {f"subsection:{level2} directory:{level3}" if level3 else f"directory:{level2}"}, "{path_value}" {"does not exist" if not Path(path_value).exists() else "is not a valid directory"} but is mandatory for tool: "{",".join(mandatory)}"')

    def _check_file_or_string(self, level1, level2, mandatory=(), level3=None, check_string=False):
        """Check if path is a file if not empty
        :return absolute path file"""
        path_value = self.get_config_value(level1=level1, level2=level2, level3=level3)
        path = Path(path_value).resolve().as_posix()
        # it is a path
        if path_value != "" and "/" in path_value:
            if not Path(path).exists() or not Path(path).is_file():
                raise FileNotFoundError(
                    f'CONFIG FILE CHECKING FAIL : in section:{level1}, {f"subsection:{level2}, file {level3}" if level3 else f"file {level2}"}, "{path}" {"does not exist" if not Path(path).exists() else "is not a valid file"}')
            else:
                self.set_config_value(level1=level1, level2=level2, level3=level3, value=path)
        # it is not a path
        elif path_value != "" and "/" not in path_value and check_string:
            self.set_config_value(level1=level1, level2=level2, level3=level3, value=path_value)
        # it is empty
        elif path_value == "" and "/" not in path_value and check_string:
            raise ValueError(
                f'CONFIG FILE CHECKING FAIL : in section:{level1}, {f"subsection:{level2}, value {level3}" if level3 else f"value:{level2}"}, "{path_value}" is empty')
        elif len(mandatory) > 0:
            raise FileNotFoundError(
                f'CONFIG FILE CHECKING FAIL : in  section:{level1} , {f"subsection:{level2}, value {level3}" if level3 else f"value:{level2}"}, "{path_value}" {"does not exist" if not Path(path_value).exists() else "is not a valid file"} but is mandatory for tool: "{",".join(mandatory)}"')

    @property
    def export_use_yaml(self):
        """Use to print a dump config.yaml with corrected parameters"""
        def represent_dictionary_order(yamldef, dict_data):
            return yamldef.represent_mapping('tag:yaml.org,2002:map', dict_data.items())

        def setup_yaml():
            yaml.add_representer(OrderedDict, represent_dictionary_order)
        setup_yaml()
        return yaml.dump(self.config, default_flow_style=False, sort_keys=False, indent=4)

    @property
    def string_to_dag(self):
        """ return command line for rule graph """
        return f"""snakemake -s {self.snakefile} {'--use-singularity' if self.use_singularity else ''} {'--use-envmodules' if self.use_env_modules else ''}  --rulegraph"""

    def __repr__(self):
        return f"{self.__class__}({pprint.pprint(self.__dict__)})"
