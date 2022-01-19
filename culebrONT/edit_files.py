import click
from shutil import copyfile
from culebrONT.global_variable import *
from culebrONT.usefull_function import get_install_mode


@click.command("edit_tools", short_help='Edit own tools version', no_args_is_help=False)
@click.option('--restore', '-r', is_flag=True, required=False, default=False, show_default=True, help='Restore default tools_config.yaml (from install)')
def edit_tools(restore):
    if restore:
        if CULEBRONT_USER_TOOLS_PATH.exists():
            CULEBRONT_USER_TOOLS_PATH.unlink()
            click.secho(f"\n    Success remove your own tools_path.yaml on path '{CULEBRONT_USER_TOOLS_PATH}'", fg="yellow")
            click.secho(f"    culebrONT used '{CULEBRONT_TOOLS_PATH}' at default now !!\n\n", fg="bright_red")
        else:
            click.secho(f"    culebrONT already used default tools_path.yaml '{CULEBRONT_TOOLS_PATH}'!!\n\n", fg="red")
    else:
        if not CULEBRONT_USER_TOOLS_PATH.exists():
            CULEBRONT_USER_TOOLS_PATH.parent.mkdir(parents=True, exist_ok=True)
            copyfile(CULEBRONT_TOOLS_PATH, CULEBRONT_USER_TOOLS_PATH)
        click.edit(require_save=True, extension='.yaml', filename=CULEBRONT_USER_TOOLS_PATH)
        click.secho(f"\n    Success install your own tools_path.yaml on path '{CULEBRONT_USER_TOOLS_PATH}'", fg="yellow")
        click.secho(f"    culebrONT used '{CULEBRONT_USER_TOOLS_PATH}' at default now !!\n\n", fg="bright_red")


@click.command("create_config", short_help='Create config.yaml for run', no_args_is_help=True)
@click.option('--configyaml', '-c', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create config.yaml')
def create_config(configyaml):
    configyaml = Path(configyaml)
    configyaml.parent.mkdir(parents=True, exist_ok=True)
    copyfile(CULEBRONT_CONFIG_PATH.as_posix(), configyaml.as_posix())
    click.edit(require_save=True, extension='.yaml', filename=configyaml)
    click.secho(f"\n    Success create config file on path '{configyaml}'\n    add to command:", fg="yellow")
    mode = get_install_mode()
    click.secho(f"    culebrONT {'run_cluster' if mode == 'cluster' else 'run_local'} --config {configyaml}\n\n", fg="bright_blue")


@click.command("create_cluster_config", short_help='Create cluster_config.yaml', no_args_is_help=True)
@click.option('--clusterconfig', '-c', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create cluster_config.yaml')
def create_cluster_config(clusterconfig):
    clusterconfig = Path(clusterconfig)
    clusterconfig.parent.mkdir(parents=True, exist_ok=True)
    copyfile(CULEBRONT_CLUSTER_CONFIG.as_posix(), clusterconfig.as_posix())
    click.edit(require_save=True, extension='.yaml', filename=clusterconfig)
    click.secho(f"\n    Success create cluster_config file on path '{clusterconfig}', add to command:", fg="yellow")
    click.secho(f"    culebrONT run_cluster --clusterconfig {clusterconfig}\n\n", fg="bright_blue")