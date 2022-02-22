import click
from shutil import copyfile
from .global_variable import *
from .usefull_function import get_install_mode, package_name


@click.command("edit_tools", short_help='Edit own tools version', no_args_is_help=False)
@click.option('--restore', '-r', is_flag=True, required=False, default=False, show_default=True, help='Restore default tools_config.yaml (from install)')
def edit_tools(restore):
    if restore:
        if USER_TOOLS_PATH.exists():
            USER_TOOLS_PATH.unlink()
            click.secho(f"\n    Success remove your own tools_path.yaml on path '{USER_TOOLS_PATH}'", fg="yellow")
            click.secho(f"    {package_name()} used '{GIT_TOOLS_PATH}' at default now !!\n\n", fg="bright_red")
        else:
            click.secho(f"    {package_name()} already used default tools_path.yaml '{GIT_TOOLS_PATH}'!!\n\n", fg="red")
    else:
        if not USER_TOOLS_PATH.exists():
            USER_TOOLS_PATH.parent.mkdir(parents=True, exist_ok=True)
            copyfile(GIT_TOOLS_PATH, USER_TOOLS_PATH)
        click.edit(require_save=True, extension='.yaml', filename=USER_TOOLS_PATH)
        click.secho(f"\n    Success install your own tools_path.yaml on path '{USER_TOOLS_PATH}'", fg="yellow")
        click.secho(f"    {package_name()} used '{USER_TOOLS_PATH}' at default now !!\n\n", fg="bright_red")


@click.command("create_config", short_help='Create config.yaml for run', no_args_is_help=True)
@click.option('--configyaml', '-c', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create config.yaml')
def create_config(configyaml):
    configyaml = Path(configyaml)
    configyaml.parent.mkdir(parents=True, exist_ok=True)
    copyfile(GIT_CONFIG_PATH.as_posix(), configyaml.as_posix())
    click.edit(require_save=True, extension='.yaml', filename=configyaml)
    click.secho(f"\n    Success create config file on path '{configyaml}'\n    add to command:", fg="yellow")
    mode = get_install_mode()
    click.secho(f"    {package_name()} {'run_cluster' if mode == 'cluster' else 'run_local'} --config {configyaml}\n\n", fg="bright_blue")


@click.command("edit_cluster_config", short_help='Edit cluster_config.yaml use by profile', no_args_is_help=False)
def edit_cluster_config():
    f"""Edit file on {USER_CLUSTER_CONFIG}"""
    USER_CLUSTER_CONFIG.parent.mkdir(parents=True, exist_ok=True)
    if not USER_CLUSTER_CONFIG.exists():
        copyfile(GIT_CLUSTER_CONFIG.as_posix(), USER_CLUSTER_CONFIG.as_posix())
    click.edit(require_save=True, extension='.yaml', filename=USER_CLUSTER_CONFIG.as_posix())
    click.secho(f"\n    See '{DOCS}' to adapt on your cluster", fg="yellow")
    click.secho(f"\n    Success edit cluster_config file on path '{USER_CLUSTER_CONFIG}'", fg="yellow")
