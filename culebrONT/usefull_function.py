from pathlib import Path
import urllib.request
from tqdm import tqdm
import os
import click


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


INSTALL_PATH = Path(__file__).resolve().parent


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
        epilogTools = ""
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

