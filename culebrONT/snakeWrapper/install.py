import click
from click import Abort
from shutil import rmtree, unpack_archive
import subprocess
from cookiecutter.main import cookiecutter
from .usefull_function import *
from .snakeWrapper import *

required_options = {
    True: 'modules_dir',
    False: 'scheduler'
}


@click.command("install_cluster", short_help=f'Install {package_name()} on HPC cluster',
               context_settings=dict(max_content_width=800),
               cls=command_required_option_from_option('create_envmodule', required_options), no_args_is_help=True)
@click.option('--scheduler', '-s', default="slurm", type=click.Choice(['slurm'], case_sensitive=False),
              prompt='Choose your HPC scheduler', show_default=True, help='Type the HPC scheduler (for the moment, only slurm is available ! )')
@click.option('--env', '-e', default="modules", type=click.Choice(['modules', 'singularity'], case_sensitive=False),
              prompt='Choose mode for tools dependencies', show_default=True, help='Mode for tools dependencies ')
@click.option('--bash_completion/--no-bash_completion', is_flag=True, required=False, default=True, show_default=True,
              help=f'Allow bash completion of {package_name()} commands on the bashrc file')
@click.option('--create_envmodule/--no-create_enmodules', is_flag=True, required=False, default=False,
              show_default=True,
              help=f'Create a env module file allowing load {package_name()} in a cluster')
@click.option('--modules_dir', '-m', default=None,
              type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True, resolve_path=True),
              required=False, show_default=True,
              help='Directory used to save the module created by --create_envmodule parameter', is_eager=True)
def install_cluster(scheduler, env, bash_completion, create_envmodule, modules_dir):
    """Run installation for HPC cluster"""
    # rm previous install (ie @julie cluster then local)
    clean_home()
    # build default profile path
    default_profile = DEFAULT_PROFILE

    def fail():
        rmtree(default_profile, ignore_errors=True)
        INSTALL_MODE.unlink(missing_ok=True)
        click.secho(f"\n    INSTALL FAIL remove already install files: {default_profile} ", fg="red", err=True)
        raise SystemExit

    # test if install already run
    try:
        if default_profile.exists() and click.confirm(
                click.style(f'    Profile "{default_profile}" exist do you want to remove and continue?\n\n', fg="red"),
                default=False, abort=True):
            rmtree(default_profile, ignore_errors=True)
        default_profile.mkdir(exist_ok=True)

        default_cluster = INSTALL_PATH.joinpath("install_files", f"cluster_config_{scheduler.upper()}.yaml").open(
            "r").read()
        command = r"""sinfo -s | grep "*" | cut -d"*" -f1 """
        default_partition = subprocess.check_output(command, shell=True).decode("utf8").strip()
        if not default_partition:
            click.secho("    Error: Slurm was not found on your system !!", fg="red", err=True)
            fail()
            raise SystemExit
        default_profile.joinpath("cluster_config.yaml").open("w").write(
            default_cluster.replace("PARTITION", default_partition))

        # Download cookiecutter of scheduler
        click.secho(
            f"    You choose '{scheduler}' as scheduler. Download cookiecutter:\n    "
            f"https://github.com/Snakemake-Profiles/{scheduler.lower()}.git",
            fg="yellow")
        cookiecutter(f'https://github.com/Snakemake-Profiles/{scheduler.lower()}.git',
                     checkout=None,
                     no_input=True,
                     extra_context={"profile_name": f'',
                                    "sbatch_defaults": "--export=ALL",
                                    "cluster_config": "$HOME/.config/culebrONT/cluster_config.yaml",
                                    "advanced_argument_conversion": 1,
                                    "cluster_name": ""
                                    },
                     replay=False, overwrite_if_exists=True,
                     output_dir=f'{default_profile}', config_file=None,
                     default_config=False, password=None, directory=None, skip_if_file_exists=True)

        try:
            # default slurm cookiecutter not contain all snakemake variables
            if scheduler == "slurm":
                default_slurm = 'restart-times: 0\njobscript: "slurm-jobscript.sh"\ncluster: ' \
                                '"slurm-submit.py"\ncluster-status: "slurm-status.py"\nmax-jobs-per-second: ' \
                                '1\nmax-status-checks-per-second: 10\nlocal-cores: 1\njobs: 200\nlatency-wait: ' \
                                '1296000\n'
                extra_slurm = f"use-envmodules: {'true' if env == 'modules' else 'false'}\nuse-singularity: " \
                              f"{'true' if env == 'singularity' else 'false'}\nrerun-incomplete: " \
                              f"true\nprintshellcmds: true"
                with open(f"{default_profile}/config.yaml", "w") as config_file:
                    config_file.write(f"{default_slurm}{extra_slurm}")

                # adding a line in RESOURCES_MAPPING dico in slurm profile
                with open(f"{default_profile}/slurm-submit.py", "r") as slurm_submit_script:
                    search_text = '"nodes": ("nodes", "nnodes"),'
                    replace_text = '"nodes": ("nodes", "nnodes"),\n    "nodelist" : ("w", "nodelist"),'
                    data = slurm_submit_script.read()
                    data = data.replace(search_text, replace_text)
                with open(f"{default_profile}/slurm-submit.py", "w") as slurm_submit_script:
                    slurm_submit_script.write(data)

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
            git_tools_file = GIT_TOOLS_PATH.open("r").read()
            GIT_TOOLS_PATH.open("w").write(git_tools_file.replace("{INSTALL_PATH}", f"{INSTALL_PATH}"))
        # if envmodule activation
        if create_envmodule:
            create_envmodules(modules_dir)
        # export to add bash completion
        if bash_completion:
            create_bash_completion()
        click.secho(f"\n    Congratulations, you have successfully installed {package_name()} !!!\n\n", fg="green", bold=True)
        INSTALL_MODE.open("w").write("cluster")
    except Exception as e:
        click.secho(f"\n    ERROR : an error was detected, please check {e}", fg="red")
        fail()


@click.command("install_local", short_help=f'Install {package_name()} on local computer',
               context_settings=dict(max_content_width=800), no_args_is_help=False)
@click.option('--bash_completion/--no-bash_completion', is_flag=True, required=False, default=True, show_default=True,
              help=f'Allow bash completion of {package_name()} commands on the bashrc file')
def install_local(bash_completion):
    """
    \b
    Run installation for local computer downloading singularity images automatically.
    """
    # rm previous install (ie @julie cluster then local)
    clean_home()
    # add path to download
    INSTALL_PATH.joinpath("containers").mkdir(exist_ok=True, parents=True)
    try:
        check_and_download_singularity()
        # export to add bash completion
        if bash_completion:
            create_bash_completion()
        # Edit tools to write good path for snakefile
        git_tools_file = GIT_TOOLS_PATH.open("r").read()
        USER_TOOLS_PATH.parent.mkdir(exist_ok=True, parents=True)
        USER_TOOLS_PATH.open("w").write(git_tools_file.replace("{INSTALL_PATH}", f"{INSTALL_PATH}"))
        # good install
        click.secho(f"\n    Congratulations, you have successfully installed {package_name()} !!!\n\n", fg="green", bold=True)
        INSTALL_MODE.open("w").write("local")
    except Exception as e:
        INSTALL_MODE.unlink(missing_ok=True)
        click.secho(f"\n    ERROR : an error was detected, please check {e}", fg="red")
        raise SystemExit


@click.command("test_install", short_help=f'Test {package_name()} with data_test after install (local or cluster)',
               context_settings=dict(max_content_width=800), no_args_is_help=False)
@click.option('--data_dir', '-d', default=None,
              type=click.Path(exists=False, file_okay=False, dir_okay=True, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to download data test and create config.yaml to run test')
def test_install(data_dir):
    """Test_install function downloads a scaled data test, writes a configuration file adapted to it and proposes a command line already to run !!!"""
    # create dir test and configure config.yaml
    data_dir = Path(data_dir).resolve()
    click.secho(f"\n    Created data test dir {data_dir}\n", fg="yellow")
    data_dir.mkdir(parents=True, exist_ok=True)

    data_config_path = data_dir.joinpath("data_test_config.yaml")
    click.secho(f"    Created config file to run data test: {data_config_path}\n", fg="yellow")
    txt = GIT_CONFIG_PATH.open("r").read().replace("DATA_DIR", f"{data_dir}").replace(
        f"./{package_name()}_OUTPUT/", f"{data_dir}/{package_name()}_OUTPUT/")
    data_config_path.open("w").write(txt)

    # download data
    download_zip = data_dir.joinpath(DATATEST_URL_FILES[0])
    if not Path(download_zip.as_posix()[:-4]).exists():
        if not download_zip.exists():
            click.secho(f"    Download data test\n", fg="yellow")
            results = multiprocessing_download(
                [(DATATEST_URL_FILES[1], download_zip.as_posix())], threads=1)
            for r in results:
                click.secho(r, fg="blue")
        click.secho(f"    Extract archive {download_zip} to {data_dir.as_posix()}\n", fg="yellow")
        unpack_archive(download_zip.as_posix(), data_dir.as_posix())
        download_zip.unlink()

    # build command line
    click.secho(f"    Write command line to run workflow on data test:\n", fg="yellow")
    mode = get_install_mode()
    cmd = f"\n    {package_name()} {'run_cluster' if mode == 'cluster' else 'run_local --threads 1'} --config {data_config_path}\n\n"
    click.secho(cmd, fg='bright_blue')
