import click
from click import Abort
from pathlib import Path
from shutil import rmtree, copyfile, unpack_archive
import culebrONT
import os
import re
from cookiecutter.main import cookiecutter
from culebrONT.global_variable import *
from culebrONT.usefull_function import command_required_option_from_option, multiprocessing_download, get_install_mode

required_options = {
    True: 'modules_dir',
    False: 'scheduler'
}

def create_bash_completion():
    bashrc_file = Path("~/.bashrc").expanduser().as_posix()
    import subprocess
    output = subprocess.run(
        ["bash", "-c", "echo ${BASH_VERSION}"], stdout=subprocess.PIPE
    )
    match = re.search(r"^(\d+)\.(\d+)\.\d+", output.stdout.decode())
    if match is not None:
        major, minor = match.groups()
    if major < "4" or (major == "4" and minor < "4"):
        click.secho(f"\n    WARNNING Shell completion is not supported for Bash versions older than 4.4.", fg="red", nl=False)
        # raise RuntimeError("Shell completion is not supported for Bash versions older than 4.4.")
    else:
        if not Path(f"{CULEBRONT_PATH}/culebrONT-complete.sh").exists():
            build_completion = f"_CULEBRONT_COMPLETE=bash_source culebrONT > {CULEBRONT_PATH}/culebrONT-complete.sh"
            os.system(build_completion)
        with open(bashrc_file, "r") as bash_file_read:
            if not [True for line in bash_file_read if "CULEBRONT" in line]:
                with open(bashrc_file, "a") as bash_file_open:
                    append_bashrc = f"\n#Add autocompletion for CULEBRONT\n. {CULEBRONT_PATH}/culebrONT-complete.sh"
                    bash_file_open.write(append_bashrc)
                    click.secho(f"\n    INSTALL autocompletion for CULEBRONT on {bashrc_file} with command {append_bashrc}\nUpdate with commande:\nsource ~/.bashrc",
                                fg="yellow")
            else:
                path_culebront = ""
                with open(bashrc_file, "r") as bash_file_read:
                    for line in bash_file_read:
                        if "CULEBRONT" in line:
                            path_bash = bash_file_read.readline().strip()
                            path_culebront = f"{line}{path_bash}"
                load = f"{path_bash[2:]}"
                if f"{CULEBRONT_PATH}/culebrONT-complete.sh" != load:
                    click.secho(
                        f"\n    WARNNING autocompletion for CULEBRONT  already found on {bashrc_file}, with other path please fix the good:",
                        fg="red", nl=False)
                    click.secho( f"\n    Load on bashrc: {load}\n    New install:    {CULEBRONT_PATH}/culebrONT-complete.sh", fg='bright_red')


def create_envmodules(modules_dir):
    from culebrONT import MODULE_FILE
    modules_dir = Path(modules_dir)
    modules_dir.mkdir(parents=True, exist_ok=True)
    modules_dir.joinpath(f"{culebrONT.__version__}").open("w").write(MODULE_FILE)
    click.edit(require_save=False, extension='', filename=modules_dir.joinpath(f"{culebrONT.__version__}").as_posix())
    click.secho(f"\n    Success install module file for version {culebrONT.__version__} on path '{modules_dir}'", fg="yellow")

def clean_home():
    if Path("~/.config/CulebrONT/").expanduser().exists():
        rmtree(Path("~/.config/CulebrONT/").expanduser().as_posix())

def check_and_download_singularity():
    # check if already download if true pop from list
    SINGULARITY_URL_FILES_DOWNLOAD = []
    for imgs_list in SINGULARITY_URL_FILES:
        url, path_install = imgs_list
        if not Path(path_install).exists():
            SINGULARITY_URL_FILES_DOWNLOAD.append(imgs_list)
        else:
            click.secho(f"    File: {path_install} already download, done.", fg="yellow", nl=True)
    results = multiprocessing_download(SINGULARITY_URL_FILES_DOWNLOAD)
    for r in results:
        click.secho(r, fg="blue")
    click.secho(f"\n    WARNNING please check if binding is active on your singularity configuration, see https://sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html", fg = "bright_red")


@click.command("install_cluster", short_help='Install CulebrONT on HPC cluster', context_settings=dict(max_content_width=800),
               cls=command_required_option_from_option('create_envmodule', required_options), no_args_is_help=True)
@click.option('--scheduler', '-s', default="slurm", type=click.Choice(['slurm', 'sge', 'lsf'], case_sensitive=False),
              prompt='Choose your HPC scheduler', show_default=True, help='Type the HPC scheduler')
@click.option('--env', '-e', default="modules", type=click.Choice(['modules', 'singularity'], case_sensitive=False),
              prompt='Choose mode for tools dependencies', show_default=True, help='Mode for tools dependencies ')
@click.option('--bash_completion/--no-bash_completion', is_flag=True, required=False, default=True, show_default=True,
              help='Allow bash completion of culebrONT commands on the bashrc file')
@click.option('--create_envmodule/--no-create_enmodules', is_flag=True, required=False, default=False, show_default=True,
              help='Create a env module file allowing load culebrONT in a cluster')
@click.option('--modules_dir', '-m', default=None,
              type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True, resolve_path=True),
              required=False, show_default=True, help='Directory used to save the module created by --create_envmodule parameter', is_eager=True)
def install_cluster(scheduler, env, bash_completion, create_envmodule, modules_dir):
    """Run installation for HPC cluster"""
    # add file for installation mode:
    mode_file = CULEBRONT_MODE
    # rm previous install (ie @julie cluster then local)
    clean_home()
    # build default profile path
    default_profile = CULEBRONT_PROFILE
    def fail():
        rmtree(default_profile, ignore_errors=True)
        mode_file.unlink(missing_ok=True)
        click.secho(f"\n    INSTALL FAIL remove  already install files: {default_profile} ", fg="red")
        exit()

    # test if install already run
    try:
        if default_profile.exists() and click.confirm(click.style(f'    Profile "{default_profile}" exist do you want to remove and continue?\n\n', fg="red"), default=False, abort=True):
            rmtree(default_profile, ignore_errors=True)
        default_profile.mkdir(exist_ok=True)

        default_cluster = CULEBRONT_PATH.joinpath("install_files",f"cluster_config_{scheduler.upper()}.yaml")
        copyfile(default_cluster, default_profile.joinpath("cluster_config.yaml"))

        # Download cookiecutter of scheduler
        click.secho(f"    You choose '{scheduler}' as scheduler. Download cookiecutter:\n    https://github.com/Snakemake-Profiles/{scheduler.lower()}.git", fg="yellow")
        cookiecutter(f'https://github.com/Snakemake-Profiles/{scheduler.lower()}.git',
                     checkout=None,
                     no_input=True,
                     extra_context={"profile_name": f'',
                                    "sbatch_defaults": "--export=ALL",
                                    "cluster_config": f"{default_profile}/cluster_config.yaml",
                                    "advanced_argument_conversion": 1,
                                    "cluster_name": ""
                                    },
                     replay=False, overwrite_if_exists=True,
                     output_dir=f'{default_profile}', config_file=None,
                     default_config=False, password=None, directory=None, skip_if_file_exists=True)

        try:
            # default slurm cookiecutter not contain all snakemake variables
            if scheduler == "slurm":
                extra_slurm = f"use-envmodules: {'true' if env == 'modules' else 'false' }\nuse-singularity: {'true' if env == 'singularity' else 'false' }\nrerun-incomplete: true\nprintshellcmds: true"
                with open(f"{default_profile}/config.yaml", "a") as config_file:
                    config_file.write(extra_slurm)

            # Edit cluster_config.yaml, config.yaml and tools_path.yaml
            if click.confirm(click.style(f'\n    Now Edit file cluster_config.yaml according to your HPC:\n    More info on on {DOCS} \n    (enter to continue)', fg="blue"), default=True, abort=True, show_default=True):
                click.edit(require_save=False, extension='.yaml', filename=f"{default_profile}/cluster_config.yaml")

            if click.confirm(click.style(f'\n    Now Edit file config.yaml according to your HPC:\n    More info on {DOCS}\n    (enter to continue)', fg="blue"), default=True, abort=True, show_default=True):
                click.edit(require_save=False, extension='.yaml', filename=f"{default_profile}/config.yaml")

            if env != "singularity":
                if click.confirm(click.style(f'\n    Now Edit file tools_path.yaml according to your HPC:\n    More info on {DOCS}\n    (enter to continue)', fg="blue"), default=True, abort=True, show_default=True):
                    click.edit(require_save=False, extension='.yaml', filename=f"{CULEBRONT_TOOLS_PATH}")
            click.secho(f"\n    Profile is success install on {default_profile}", fg="yellow")
        except Abort:
            fail()
    except Abort:
        click.secho(f"\n    Profile is already created, skipping {default_profile}", fg="yellow")
    except Exception as e:
        print(e)
        fail()
    try:
        # if singularity activation
        if env == 'singularity':
            # check if already download
            check_and_download_singularity()
        # if envmodule activation
        if create_envmodule:
            create_envmodules(modules_dir)
        # export to add bash completion
        if bash_completion:
            create_bash_completion()
        click.secho(f"\n    Congratulations, you have successfully installed culebrONT !!!\n\n", fg="green", bold=True)
        mode_file.open("w").write("cluster")
    except Exception as e:
        click.secho(f"\n    ERROR : an error was detected, please check {e}",  fg="red")
        fail()


@click.command("install_local", short_help='Install CulebrONT on local computer', context_settings=dict(max_content_width=800), no_args_is_help=False)
@click.option('--bash_completion/--no-bash_completion', is_flag=True, required=False, default=True, show_default=True,
              help='Allow bash completion of culebrONT commands on the bashrc file')
def install_local(bash_completion):
    """
    \b
    Run installation for local computer with download the singularity images.
    """
    # add file for installation mode:
    mode_file = CULEBRONT_MODE
    # rm previous install (ie @julie cluster then local)
    clean_home()
    # add path to download
    CULEBRONT_PATH.joinpath("containers").mkdir(exist_ok=True, parents=True)
    try:

        check_and_download_singularity()
        # export to add bash completion
        if bash_completion:
           create_bash_completion()
        click.secho(f"\n    Congratulations, you have successfully installed culebrONT !!!\n\n", fg="green", bold=True)
        mode_file.open("w").write("local")
    except Exception as e:
        mode_file.unlink(missing_ok=True)
        click.secho(f"\n    ERROR : an error was detected, please check {e}", fg="red")
        exit()


@click.command("test_install", short_help='Test culebrONT with data_test after install (local or cluster)', context_settings=dict(max_content_width=800), no_args_is_help=False)
@click.option('--data_dir', '-d', default=None,
              type=click.Path(exists=False, file_okay=False, dir_okay=True, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to download datatest and create config.yaml to run test')
def test_install(data_dir):
    # create dir test and configure config.yaml
    data_dir = Path(data_dir)
    data_config_path = data_dir.joinpath("data_test_config.yaml")
    data_dir.mkdir(parents=True, exist_ok=True)
    txt =  CULEBRONT_CONFIG_PATH.open("r").read().replace("./Data-Xoo-sub",f"{data_dir}/Data-Xoo-sub").replace("./CulebrONT_OUTPUT/",f"{data_dir}/CulebrONT_OUTPUT/")
    data_config_path.open("w").write(txt)

    # download data
    download_zip = data_dir.joinpath("Data-Xoo-sub.zip")
    results = multiprocessing_download([("https://itrop.ird.fr/culebront_utilities/Data-Xoo-sub.zip", download_zip.as_posix())], threads=1)
    for r in results:
        click.secho(r, fg="blue")
    unpack_archive(download_zip, data_dir)
    download_zip.unlink()

    # build commande line
    mode = get_install_mode()
    cmd = f"\n    culebrONT {'run_cluster' if mode == 'cluster' else 'run_local --threads 1'} --config {data_config_path}\n\n"
    click.secho(cmd, fg='bright_blue')