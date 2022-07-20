import click
from shutil import copyfile
from .global_variable import *
from .usefull_function import package_name
import os
import re
import subprocess

def rewrite_if_bind(snakemake_other):
    """
    Function to rewrite --bind params
    It modifies click.UNPROCESSED
    """
    bind_args = list(filter(re.compile(".*--bind.*").match, snakemake_other))  # Try to select --bind
    if bind_args:
        bind_args_rewrite = f'"--bind {bind_args[0].split(" ")[1]}"'
        snakemake_other_list = list(filter(lambda x: x not in bind_args[0], snakemake_other))  # remove value to rewrite
        snakemake_other_list.append(bind_args_rewrite)
        return snakemake_other_list
    else:
        return snakemake_other


def build_pdf(cmd_snakemake_base):
    dag_cmd_snakemake = f"{cmd_snakemake_base} --dag | dot -Tpdf > schema_pipeline_dag.pdf"
    click.secho(f"    {dag_cmd_snakemake}\n", fg='bright_blue')
    process = subprocess.run(dag_cmd_snakemake, shell=True, check=False, stdout=sys.stdout, stderr=sys.stderr)
    if int(process.returncode) >= 1:
        sys.exit(1)
    rulegraph_cmd_snakemake = f"{cmd_snakemake_base} --rulegraph | dot -Tpdf > schema_pipeline_global.pdf"
    click.secho(f"    {rulegraph_cmd_snakemake}\n", fg='bright_blue')
    process = subprocess.run(rulegraph_cmd_snakemake, shell=True, check=False, stdout=sys.stdout, stderr=sys.stderr)
    if int(process.returncode) >= 1:
        sys.exit(1)
    filegraph_cmd_snakemake = f"{cmd_snakemake_base} --filegraph | dot -Tpdf > schema_pipeline_files.pdf"
    click.secho(f"    {filegraph_cmd_snakemake}\n", fg='bright_blue')
    process = subprocess.run(filegraph_cmd_snakemake, shell=True, check=False, stdout=sys.stdout, stderr=sys.stderr)
    if int(process.returncode) >= 1:
        sys.exit(1)


@click.command("run_cluster", short_help='Run workflow on HPC', context_settings=dict(ignore_unknown_options=True, max_content_width=800), no_args_is_help=True)
@click.option('--config', '-c', type=click.Path(exists=True, file_okay=True, readable=True, resolve_path=True), required=True, show_default=True, help=f'Configuration file for run {package_name()}')
@click.option('--pdf', '-pdf', is_flag=True, required=False, default=False, show_default=True, help='Run snakemake with --dag, --rulegraph and --filegraph')
@click.argument('snakemake_other', nargs=-1, type=click.UNPROCESSED)
def run_cluster(config, pdf, snakemake_other):
    """
    \b
    Run snakemake command line with mandatory parameters.
    SNAKEMAKE_OTHER: You can also pass additional Snakemake parameters
    using snakemake syntax.
    These parameters will take precedence over Snakemake ones, which were
    defined in the profile.
    See: https://snakemake.readthedocs.io/en/stable/executing/cli.html

    Example:
        culebrONT run_cluster -c config.yaml --dry-run --jobs 200
    """
    profile = DEFAULT_PROFILE
    tools = GIT_TOOLS_PATH
    # get user arguments
    click.secho(f'    Profile file: {profile}', fg='yellow')
    click.secho(f'    Config file: {config}', fg='yellow')

    if USER_CLUSTER_CONFIG.exists():
        clusterconfig = USER_CLUSTER_CONFIG
    else:
        click.secho(f"    Please run command line '{package_name()} edit_cluster_config' before the first run of {package_name()} see {DOCS}", fg="red", err=True)
        exit()
    cmd_clusterconfig = f"--cluster-config {clusterconfig}"
    click.secho(f'    Cluster config file load: {clusterconfig}', fg='yellow')

    if tools != GIT_TOOLS_PATH.as_posix():
        ARGS_TOOLS_PATH.parent.mkdir(parents=True, exist_ok=True)
        copyfile(tools, ARGS_TOOLS_PATH)
    elif USER_TOOLS_PATH.exists():
        tools = USER_TOOLS_PATH
    click.secho(f'    Tools Path file: {tools}', fg='yellow')

    cmd_snakemake_base = f"snakemake --show-failed-logs -p -s {SNAKEFILE} --configfile {config} --profile {profile} {cmd_clusterconfig} {' '.join(rewrite_if_bind(snakemake_other))}"
    click.secho(f"\n    {cmd_snakemake_base}\n", fg='bright_blue')
    process = subprocess.run(cmd_snakemake_base, shell=True, check=False, stdout=sys.stdout, stderr=sys.stderr)
    if int(process.returncode) >= 1:
        sys.exit(1)

    if pdf:
        build_pdf(cmd_snakemake_base)


@click.command("run_local", short_help='Run a workflow on local computer (use singularity mandatory)', context_settings=dict(ignore_unknown_options=True, max_content_width=800),
               no_args_is_help=True)
@click.option('--config', '-c', type=click.Path(exists=True, file_okay=True, readable=True, resolve_path=True), required=True, help=f'Configuration file for run {package_name()}')
@click.option('--threads', '-t', type=int, required=True, help='Number of threads')
@click.option('--pdf', '-p', is_flag=True, required=False, default=False, help='Run snakemake with --dag, --rulegraph and --filegraph')
@click.argument('snakemake_other', nargs=-1, type=click.UNPROCESSED)
def run_local(config, threads, pdf, snakemake_other):
    """
    \b
    Run snakemake command line with mandatory parameters.
    SNAKEMAKE_OTHER: You can also pass additional Snakemake parameters
    using snakemake syntax.
    These parameters will take precedence over Snakemake ones, which were
    defined in the profile.
    See: https://snakemake.readthedocs.io/en/stable/executing/cli.html

    Example:

        culebrONT run_local -c config.yaml --threads 8 --dry-run

        culebrONT run_local -c config.yaml --threads 8 --singularity-args '--bind /mnt:/mnt'

        # in LOCAL using 6 threads for Canu assembly from the total 8 threads\n
        culebrONT run_local -c config.yaml --threads 8 --set-threads run_canu=6

    """
    click.secho(f'    Config file: {config}', fg='yellow')
    click.secho(f'    Tools config file: {USER_TOOLS_PATH}', fg='yellow')

    cmd_snakemake_base = f"snakemake --latency-wait 1296000 --cores {threads} --use-singularity --show-failed-logs --printshellcmds -s {SNAKEFILE} --configfile {config}  {' '.join(rewrite_if_bind(snakemake_other))}"
    click.secho(f"\n    {cmd_snakemake_base}\n", fg='bright_blue')
    process = subprocess.run(cmd_snakemake_base, shell=True, check=False, stdout=sys.stdout, stderr=sys.stderr)
    if int(process.returncode) >= 1:
        sys.exit(1)

    if pdf:
        build_pdf(cmd_snakemake_base)
