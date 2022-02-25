import os
import re
import subprocess
import urllib.request
from shutil import rmtree
import click
from tqdm import tqdm
from .global_variable import *


def check_privileges():
    if not os.environ.get("SUDO_UID") and os.geteuid() != 0:
        click.secho(f"\n    ERROR : You need to run -r, --restore with sudo privileges or as root\n", fg="red")
    else:
        return True


class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(tuple_value):
    url, output_path = tuple_value
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def multiprocessing_download(urls_list, threads = 2):
    from multiprocessing.pool import ThreadPool
    return ThreadPool(threads).imap_unordered(download_url, urls_list)


INSTALL_PATH = Path(__file__).resolve().parent.parent


def package_name():
    return Path(INSTALL_PATH).stem


def get_install_mode():
    """detect install mode"""
    if INSTALL_PATH.joinpath(".mode.txt").exists():
        return INSTALL_PATH.joinpath(".mode.txt").open("r").readline().strip()
    else:
        return "notInstall"


def get_version():
    """Read VERSION file to know current version
    Returns:
        version: actual version read on the VERSION file
    Examples:
        version = get_version()
        print(version)
            1.3.0
    """
    with open(INSTALL_PATH.joinpath("VERSION"), 'r') as version_file:
        return version_file.readline().strip()


def get_last_version(url, current_version):
    """Function for know the last version of Git repo in website"""
    try:
        from urllib.request import urlopen
        from re import search
        import click
        module_mane = url.split('/')[-1]
        HTML = urlopen(f"{url}/tags").read().decode('utf-8')
        str_search = f"{url.replace('https://github.com', '')}/releases/tag/.*"
        lastRelease = search(str_search, HTML).group(0).split("/")[-1].split('"')[0]
        epilogTools = "\n"
        if str(current_version) != lastRelease:
            if lastRelease < str(current_version):
                epilogTools = click.style(f"\n    ** NOTE: This {module_mane} version ({current_version}) is higher than the production version ({lastRelease}), you are using a dev version\n\n", fg="yellow", bold=True)
            elif lastRelease > str(current_version):
                epilogTools = click.style(f"\n    ** NOTE: The Latest version of {module_mane} {lastRelease} is available at {url}/releases\n\n",fg="yellow", underline=True)
        return epilogTools
    except Exception as e:
        epilogTools = click.style(f"\n    ** ENABLE TO GET LAST VERSION, check internet connection\n{e}\n\n", fg="red")
        return epilogTools


def command_required_option_from_option(require_name, require_map):
    import click

    class CommandOptionRequiredClass(click.Command):
        def invoke(self, ctx):
            require = ctx.params[require_name]
            if require not in require_map:
                raise click.ClickException(click.style(f"Unexpected value for --'{require_name}': {require}\n", fg="red"))
            if ctx.params[require_map[require].lower()] is None:
                raise click.ClickException(click.style(f"With {require_name}={require} must specify option --{require_map[require]} path/to/modules\n", fg="red"))
            super(CommandOptionRequiredClass, self).invoke(ctx)
    return CommandOptionRequiredClass


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

     """
    if not (extensions, (list, tuple)) or not extensions:
        raise ValueError(f'ERROR: "extensions" must be a list or tuple not "{type(extensions)}"')
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


def var_2_bool(key, tool, to_convert):
    """convert to boolean"""
    if isinstance(type(to_convert), bool):
        return to_convert
    elif f"{to_convert}".lower() in ("yes", "true", "t"):
        return True
    elif f"{to_convert}".lower() in ("no", "false", "f"):
        return False
    else:
        raise TypeError(
            f'CONFIG FILE CHECKING FAIL : in the "{key}" section, "{tool}" key: "{to_convert}" is not a valid boolean')


def create_bash_completion():
    bashrc_file = Path("~/.bashrc").expanduser().as_posix()
    output = subprocess.run(
        ["bash", "-c", "echo ${BASH_VERSION}"], stdout=subprocess.PIPE
    )
    match = re.search(r"^(\d+)\.(\d+)\.\d+", output.stdout.decode())
    if match is not None:
        major, minor = match.groups()
    if major < "4" or (major == "4" and minor < "4"):
        click.secho(f"\n    WARNING Shell completion is not supported for Bash versions older than 4.4.", fg="red", nl=False)
        # raise RuntimeError("Shell completion is not supported for Bash versions older than 4.4.")
    else:
        if not Path(f"{INSTALL_PATH}/{package_name()}-complete.sh").exists():
            build_completion = f"_{package_name().upper()}_COMPLETE=bash_source {package_name()} > {INSTALL_PATH}/{package_name()}-complete.sh"
            os.system(build_completion)
            path_bash = None
            with open(bashrc_file, "r") as bash_file_read:
                for line in bash_file_read:
                    if f"{package_name().upper()}" in line:
                        path_bash = bash_file_read.readline().strip()
            if path_bash:
                load = f"{path_bash[2:]}"
                if f"{INSTALL_PATH}/{package_name()}-complete.sh" != load:
                    click.secho(f"\n    WARNING autocompletion for {package_name().upper()} already found on {bashrc_file}, with other path please fix the good:", fg="red", nl=False)
                    click.secho(f"\n    Load on bashrc: {load}\n    New install:    {INSTALL_PATH}/{package_name()}-complete.sh", fg='bright_red')
            else:
                with open(bashrc_file, "a") as bash_file_open:
                    append_bashrc = f"\n#Add autocompletion for {package_name().upper()}\n. {INSTALL_PATH}/{package_name()}-complete.sh"
                    bash_file_open.write(append_bashrc)
                    click.secho(f"\n    INSTALL autocompletion for {package_name().upper()} on {bashrc_file} with command {append_bashrc}\nUpdate with command:\nsource ~/.bashrc", fg="yellow")


def create_envmodules(modules_dir):
    from .. import MODULE_FILE, __version__
    modules_dir = Path(modules_dir)
    modules_dir.mkdir(parents=True, exist_ok=True)
    modules_dir.joinpath(f"{__version__}").open("w").write(MODULE_FILE)
    click.edit(require_save=False, extension='', filename=modules_dir.joinpath(f"{__version__}").as_posix())
    click.secho(f"\n    Success install module file for version {__version__} on path '{modules_dir}'", fg="yellow")


def clean_home():
    if Path(f"~/.config/{package_name()}/").expanduser().exists():
        rmtree(Path(f"~/.config/{package_name()}/").expanduser().as_posix())


def check_and_download_singularity():
    # check if already download if true pop from list
    from . import SINGULARITY_URL_FILES, INSTALL_PATH
    SINGULARITY_URL_FILES_DOWNLOAD = []
    for imgs_list in SINGULARITY_URL_FILES:
        url, path_install = imgs_list
        if not Path(path_install).exists():
            SINGULARITY_URL_FILES_DOWNLOAD.append(imgs_list)
            click.secho(f"    File: {path_install} is download.", fg="yellow", nl=True)
        else:
            click.secho(f"    File: {path_install} already download, done.", fg="yellow", nl=True)
    results = multiprocessing_download(SINGULARITY_URL_FILES_DOWNLOAD)
    for r in results:
        click.secho(r, fg="blue")
    click.secho(f"\n    WARNING please check if binding is active on your singularity configuration, see https://sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html", fg = "bright_red")


def convert_genome_size(size):
    import re
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

