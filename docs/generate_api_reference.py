#!/usr/bin/env python3

"""Generate API reference documentation for a package."""

import importlib
import os
import pathlib
import pkgutil
import re
import sys

import setuptools

REPOSITORY_ROOT_PATH = str(pathlib.Path(os.path.abspath(__file__)).parent.parent)

sys.path.insert(0, REPOSITORY_ROOT_PATH)

ROOT_PACKAGE_PATH = os.path.join(REPOSITORY_ROOT_PATH, "gnomad_qc")

DOCS_DIRECTORY = os.path.join(REPOSITORY_ROOT_PATH, "docs")

FLOWCHART_DIRECTORY = os.path.join(DOCS_DIRECTORY, "flowcharts")

EXCLUDE_PACKAGES = []

EXCLUDE_TOP_LEVEL_PACKAGES = [
    "gnomad_qc.v2", "gnomad_qc.v3", "gnomad_qc.example_notebooks"
]

PACKAGE_DOC_TEMPLATE = """{title}

{package_doc}

{flowchart}

.. toctree::
    :maxdepth: 2

    {module_links}
"""

MODULE_DOC_TEMPLATE = """{title}

{module_doc}

{flowchart}

{argparse_doc}

Module Functions
****************

.. gnomad_automodulesummary:: {module_name}

.. automodule:: {module_name}
    :exclude-members: get_script_argument_parser
"""

ARGPARSE_TEMPLATE = """
.. argparse::
   :ref: {module_name}.get_script_argument_parser
   :prog: {module_name}.py

"""

MERMAID_TEMPLATE = """
.. mermaid:: {mmd_path}
"""


def module_doc_path(module):
    """Get the path for a module's documentation file."""
    return os.path.join(
        DOCS_DIRECTORY,
        "api_reference",
        re.sub(r"\.py$", ".rst", os.path.relpath(module.__file__, ROOT_PACKAGE_PATH)),
    )


def module_flowchart_path(module, local=False):
    """
    Get the path for a module's mermaid flowchart file.

    If local is True, returns the absolute path to the file.
    If local is False, returns the path relative to the module's documentation file.

    :param module: Module to get flowchart path for.
    :param local: Whether to return the absolute path to the file.
        Default is False.
    :return: Path to the module's mermaid flowchart file.
    """
    mmd_path = os.path.join(
        FLOWCHART_DIRECTORY,
        re.sub(r"\.py$", ".mmd", os.path.relpath(module.__file__, ROOT_PACKAGE_PATH))
    )
    if local:
        return mmd_path
    else:
        return os.path.relpath(mmd_path, os.path.split(module_doc_path(module))[0])


def package_doc_path(package):
    """Get the path for a package's documentation file."""
    return os.path.join(
        DOCS_DIRECTORY,
        "api_reference",
        os.path.relpath(package.__path__[0], ROOT_PACKAGE_PATH),
        "index.rst",
    )


def package_flowchart_path(package, local=False):
    """
    Get the path for a package's mermaid flowchart file.

    If local is True, returns the absolute path to the file.
    If local is False, returns the path relative to the package's documentation file.

    :param package: Package to get flowchart path for.
    :param local: Whether to return the absolute path to the file.
        Default is False.
    :return: Path to the package's mermaid flowchart file.
    """
    mmd_path = os.path.join(
        FLOWCHART_DIRECTORY,
        os.path.relpath(package.__path__[0], ROOT_PACKAGE_PATH),
        "index.mmd"
    )
    if local:
        return mmd_path
    else:
        return os.path.relpath(mmd_path, os.path.split(package_doc_path(package))[0])


def write_file(path, contents):
    """Write a file after ensuring that the target directory exists."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as out:
        out.write(contents)


def format_title(title):
    """
    Format title for reST.

    reST requires header text to have an underline at least as long as the text.
    """
    underline = "=" * len(title)
    return f"{title}\n{underline}"


def write_module_doc(module_name):
    """Write API reference documentation file for a module."""
    module = importlib.import_module(module_name)

    if hasattr(module, "get_script_argument_parser"):
        argparse_doc = ARGPARSE_TEMPLATE.format(module_name=module_name)
    else:
        argparse_doc = ""

    if os.path.exists(module_flowchart_path(module, local=True)):
        flowchart = MERMAID_TEMPLATE.format(mmd_path=module_flowchart_path(module))
    else:
        flowchart = ""

    doc = MODULE_DOC_TEMPLATE.format(
        module_name=module_name,
        title=format_title(module_name),
        module_doc=module.__doc__ or "",
        flowchart=flowchart,
        argparse_doc=argparse_doc,
    )

    doc_path = module_doc_path(module)
    write_file(doc_path, doc)


def write_package_doc(package_name):
    """Write API reference documentation file for a package."""
    package = importlib.import_module(package_name)

    module_links = []

    for module in pkgutil.iter_modules(package.__path__):
        if module.name in EXCLUDE_PACKAGES:
            continue

        full_module_name = f"{package_name}.{module.name}"
        if module.ispkg:
            write_package_doc(full_module_name)
            module_links.append(f"{module.name} <{module.name}/index>")
        else:
            write_module_doc(full_module_name)
            module_links.append(f"{module.name} <{module.name}>")

    if os.path.exists(package_flowchart_path(package, local=True)):
        flowchart = MERMAID_TEMPLATE.format(mmd_path=package_flowchart_path(package))
    else:
        flowchart = ""

    doc = PACKAGE_DOC_TEMPLATE.format(
        title=format_title(package_name),
        package_doc=package.__doc__ or "",
        flowchart=flowchart,
        module_links="\n    ".join(module_links),
    )

    doc_path = package_doc_path(package)
    write_file(doc_path, doc)


if __name__ == "__main__":
    packages = setuptools.find_namespace_packages(
        where=REPOSITORY_ROOT_PATH, include=["gnomad_qc.*"]
    )
    top_level_packages = [
        pkg for pkg in packages if pkg.count(".") == 1 and pkg not in EXCLUDE_TOP_LEVEL_PACKAGES
    ]

    for pkg in top_level_packages:
        write_package_doc(pkg)

    root_doc = PACKAGE_DOC_TEMPLATE.format(
        title=format_title("gnomad_qc"),
        package_doc="",
        flowchart="",
        module_links="\n    ".join(
            f"{pkg.split('.')[1]} <{pkg.split('.')[1]}/index>"
            for pkg in top_level_packages
        ),
    )

    write_file(os.path.join(DOCS_DIRECTORY, "api_reference", "index.rst"), root_doc)
