"""
This module contains utility functions for the Research Workbench.

These functions should live in the 'util_functions.ipynb.' notebook. It is included here
for version control purposes and cannot be imported via the gnomad_qc module in
Researcher Workbench unless the gnomad_qc repo is cloned prior to use.
"""

import subprocess
import sys
from typing import Optional

from IPython.display import Javascript, display


def restart_kernel_with_gnomad_packages(
    qc_branch: Optional[str] = None,
    methods_branch: Optional[str] = None,
    directory="/home/jupyter/packages/",
) -> None:
    """
    Clone the gnomad_qcrepository, update the gnomAD package, and restart the kernel.

    :param qc_branch: Branch of the gnomad_qc repository to clone, defaults to None.
    :param methods_branch: Branch of the gnomad_methods repository to clone, defaults to None.
    :param directory: Package directory, defaults to "/home/jupyter/packages/".
    """
    subprocess.run(["mkdir", directory])

    def _clone_and_move_repo_command(
        repo: str,
        branch: Optional[str] = None,
        secondary_name=None,
        directory=directory,
        install_reqs=False,
    ) -> None:
        """
        Clone a git repository.

        :param repo: Broad git repository name
        :param branch: Feature branch, defaults to None and thus main
        :param secondary_name: Secondary name of the repository, defaults to None
        :param directory: Package directory, defaults to "/home/jupyter/packages/"
        :param install_reqs: Install requirements, defaults to False
        """
        clone_command = [
            "git",
            "clone",
            f"https://github.com/broadinstitute/{repo}.git",
        ]

        if branch:
            clone_command.extend(["--branch", branch, "--single-branch"])

        subprocess.run(clone_command, check=True)

        if install_reqs:
            req_install_command = ["pip", "install", "-r", f"{repo}/requirements.txt"]
            subprocess.run(req_install_command, check=True)

        # The way we import the gnomad packages is by adding the package to the python
        # path, e.g. from gnomad import sample_qc, so we need to move the secondary name
        # argument into the package directory or well need import with the syntax
        # from gnomad_methods.gnomad import sample_qc which is not what we want given
        # the way the gnomad package is structured. Using rsync allows us to move the
        # package to the package directory multiple times on the same cluster.
        rsync_command = [
            "rsync",
            "-r",
            "--update",
            "--remove-source-files",
            f"{repo}/{secondary_name}" if secondary_name else repo,
            directory,
        ]
        subprocess.run(rsync_command, check=True)

        # Remove the repository after moving the package to the package directory so
        # we can clone the same repository more than once, i.e. you push a commit to a
        # branch and do not want to spin up a new cluster, you can rerun this function.
        # This is necessary as we cannot overwrite when cloning.
        rm_command = ["rm", "-rf", repo]
        subprocess.run(rm_command, check=True)

    _clone_and_move_repo_command("gnomad_qc", qc_branch, "gnomad_qc")
    _clone_and_move_repo_command(
        "gnomad_methods", methods_branch, "gnomad", install_reqs=True
    )

    # Put installed packages first in the python path. This is necessary as we cannot
    # overwrite or delete the preinstalled gnomad package
    sys.path.insert(0, "/home/jupyter/packages")

    # Restart the kernel.
    display(Javascript("Jupyter.notebook.kernel.restart()"))
