import networkx as nx
import pkgutil
import inspect
import itertools
import os
import re
import sys
import pathlib
import importlib
import setuptools

import hail as hl

from gnomad.resources.resource_utils import (
    BaseResource,
    DataException,
    BaseVersionedResource,
)
from gnomad.utils.file_utils import check_file_exists_raise_error

# Need to add the repository root to the path before importing any gnomad_qc modules.
REPOSITORY_ROOT_PATH = str(pathlib.Path(os.path.abspath(__file__)).parent.parent)
sys.path.insert(0, REPOSITORY_ROOT_PATH)

from gnomad_qc.v4.resources.variant_qc import TRUTH_SAMPLES

ROOT_PACKAGE_PATH = os.path.join(REPOSITORY_ROOT_PATH, "gnomad_qc")
"""Path to the root package in the repository."""

DOCS_DIRECTORY = os.path.join(REPOSITORY_ROOT_PATH, "docs")
"""Path to the docs directory in the repository."""

FLOWCHART_DIRECTORY = os.path.join(DOCS_DIRECTORY, "flowcharts")
"""Path to the flowcharts directory in the docs."""

API_REFERENCE_PATH = "https://broadinstitute.github.io/gnomad_qc/api_reference"
"""Path to the API reference documentation on GitHub Pages."""

INCLUDE_TOP_LEVEL_PACKAGES = {"gnomad_qc.v4"}
"""Top-level packages to include in the flowchart construction."""

EXCLUDE_RESOURCES = {
    "get_checkpoint_path",
    "get_logging_path",
    "add_meta",
    "get_gnomad_v4_vds"
}
"""Resources functions to exclude from the flowchart resource nodes."""

EXCLUDE_RESOURCE_PARAMS = {"test", "model_id"}
"""Resource function parameters to exclude from the resource node label."""

RESOURCE_PACKAGES = ["gnomad_qc.v4.resources"]
"""Resource packages to include in the flowchart resource nodes."""

KWARG_OPTIONS_BY_RESOURCE = {
    "get_freq_comparison": {"method": ["contingency_table_test", "cmh_test"]},
    "relatedness": {"method": [None, "cuking", "pc_relate"]},
    "hgdp_tgp_updated_callstats": {
        "subset": [
            'added',
            'subtracted',
            'pop_diff',
            'join',
            'v3_release_an',
            'v3_pop_diff_an',
            'pre_validity_check'
        ]
    },
    "calling_intervals": {
        "interval_name": ["broad", "ukb", "union", "intersection"],
        "calling_interval_padding": [0, 50]
    },
    "get_binned_concordance": {
        "model_id": ["vqsr_quasi_adj_ts_ss_calling", "rf_ba64690f"],
        "truth_sample": list(TRUTH_SAMPLES.keys()),
    },
    "get_callset_truth_data": {
        "truth_sample": list(TRUTH_SAMPLES.keys()),
    },
    "get_rf_training": {
        "model_id": ["rf_ba64690f"],
    },
    "get_score_bins": {
        "model_id": ["vqsr_quasi_adj_ts_ss_calling", "rf_ba64690f"],
    },
    "get_variant_qc_result": {
        "model_id": ["vqsr_quasi_adj_ts_ss_calling", "rf_ba64690f"],
    },
    "get_rf_model_path": {
        "model_id": ["rf_ba64690f"],
    },
    "get_freq": {
        "intermediate_subset": [None, "pre_hom_alt_adjustment"]
    },
    "get_sample_qc": {
        "strat": ["all", "bi_allelic", "under_three_alt_alleles"]
    }
}
"""Resource function parameter options to include in the resource nodes."""

NODE_COLORS = {
    "step": "#98BFFF",  # Light Blue
    "func": "#E0D5E6",  # Light Purple
    "gnomad_methods": "#F9F2CE",  # Light Yellow
    "hail": "#FAEACE",  # Light Orange
    "resource": "#D6D6D6",  # Light Gray
    "validity_check": "#2C5D4A",  # Dark Green
}
"""
Types of nodes:

    - Pipeline step
    - Function call of function within script
    - Function call of function in gnomAD methods
    - Function call of function in hail
    - Resource
    - Validity checks
"""


class QCFlowchart:
    """
    Class to represent a directed graph (flowchart) of a quality control pipeline.
    """

    def __init__(self):
        """
        Initialize a QCFlowchart object.
        """
        self.graph = nx.DiGraph()

    def add_node(self, node, **kwargs):
        """
        Add a node to the graph.
        """
        self.graph.add_node(node, **kwargs)

    def add_edge(self, edge, **kwargs):
        """
        Add an edge to the graph.
        """
        self.graph.add_edge(edge, **kwargs)

    def get_graph(self):
        """
        Get the graph.
        """
        return self.graph


class ResourceNode:

    def __init__(self, res_name, res_obj, module):
        """
        Initialize a ResourceNode object.

        :param res_name: Function name.
        :param res_obj: Resource function object.
        """
        self.res_name = res_name
        self.res_obj = res_obj
        self.module = module

        self.is_resource = self.get_is_resource()
        self.is_module_func = self.get_is_module_func()

    def get_is_resource(self):
        """
        Check if the resource object is a resource.
        """
        return (
            issubclass(type(self.res_obj), BaseResource)
            or issubclass(type(self.res_obj), BaseVersionedResource)
        )

    def get_is_module_func(self):
        """
        Check if the resource object is a function defined in the current module.
        """
        return (
            inspect.isfunction(self.res_obj)
            and self.res_obj.__module__ == self.module.name
        )

    def get_signature(self):
        """
        Get the signature of the resource function.
        """
        return inspect.signature(self.res_obj)

    def get_resource_kwargs(self):
        """
        Get Dictionary of kwargs for the ResourceNode function.

        :return: Dictionary of kwargs for the ResourceNode function.
        """
        kwargs_options = KWARG_OPTIONS_BY_RESOURCE.get(self.res_name, {})
        for k, v in self.get_signature().parameters.items():
            if (k == "test") or (k == "overwrite"):
                kwargs_options[k] = [False]
            elif v.annotation.__name__ == "bool":
                kwargs_options[k] = [True, False]
            elif k == "data_type":
                kwargs_options[k] = ["exomes", "genomes", "joint"]

        return kwargs_options

    def get_resource_kwargs_combos(self):
        """
        Get all combinations of kwargs for the ResourceNode function.

        :return: List of all combinations of kwargs for the ResourceNode function.
        """
        all_options = []
        for k, v in self.get_resource_kwargs().items():
            all_options.append([(k, val) for val in v])

        return list(itertools.product(*all_options))

    @staticmethod
    def remove_excluded_parameters(bind_kwargs):
        """
        Remove excluded parameters from the resource function kwargs.

        :return: Dictionary of resource function kwargs with excluded parameters removed.
        """
        return {
            k: v for k, v in bind_kwargs if k not in EXCLUDE_RESOURCE_PARAMS
        }

    def get_resource_display_name(self, bind_kwargs=None):
        """
        Get the display name for the ResourceNode.

        :param bind_kwargs: Optional resource function kwargs to run the function with.
        :return: Display name for the ResourceNode.
        """
        if bind_kwargs is None:
            res_display = self.res_name
        else:
            bind_kwargs = self.remove_excluded_parameters(bind_kwargs)
            res_display = ",\n".join("{}={}".format(k, v) for k, v in bind_kwargs)

            if len(bind_kwargs) > 1:
                res_display = "\n" + res_display

        return f"{self.res_name}({res_display})"

    def get_resource_id(self, bind_kwargs=None):
        """
        Get the ID for the ResourceNode.

        :param bind_kwargs: Optional resource function kwargs to run the function with.
        :return: ID for the ResourceNode.
        """
        if bind_kwargs is None:
            return self.res_name
        else:
            res_arg_str = "_".join(
                "{}_{}".format(k, v)
                for k, v in self.remove_excluded_parameters(bind_kwargs)
            )
            return f"{self.res_name}_{res_arg_str}"

    def add_resource_to_flowchart(self, flowchart, resource, bind_kwargs=None):
        """
        Run a resource function and add its resource nodes to the flowchart.

        :param flowchart: QCFlowchart object.
        """
        if hasattr(resource, "versions"):
            paths = [r.path for r in resource.versions.values()]
        elif hasattr(resource, "path"):
            paths = [resource.path]
        elif isinstance(resource, str):
            paths = [resource]
        else:
            err_msg = "" if bind_kwargs is None else f"with {bind_kwargs}"
            raise ValueError(
                f"Resource {self.res_obj} {err_msg}did not return a valid resource."
            )

        for p in paths:
            if check_file_exists_raise_error(p, error_if_exists=False):
                flowchart.add_node(
                    p,
                    type="resource",
                    color=NODE_COLORS["resource"],
                    label=self.get_resource_display_name(bind_kwargs),
                    id=self.get_resource_id(bind_kwargs),
                    href=get_function_doc_link(self.module, self.res_name),
                    res_node=self,
                )

    def add_to_flowchart(self, flowchart):
        """
        Run a resource function and add its resource nodes to the flowchart.

        :param flowchart: QCFlowchart object.
        """
        if self.is_module_func:
            res_sig = self.get_signature()
            for bind_kwargs in self.get_resource_kwargs_combos():
                res_bind = res_sig.bind(**dict(bind_kwargs))

                try:
                    self.add_resource_to_flowchart(
                        flowchart,
                        self.res_obj(*res_bind.args, **res_bind.kwargs),
                        bind_kwargs=bind_kwargs,
                    )
                except (ValueError, DataException, KeyError):
                    continue
        else:
            self.add_resource_to_flowchart(flowchart, self.res_obj, bind_kwargs=None)


def parse_resource_package(flowchart, resource_package, full_package_name):
    """
    Parse a resource package and add its resource nodes to the QCFlowchart.

    :param flowchart: QCFlowchart object.
    :param resource_package: Package object for resources.
    :param full_package_name: Full package name.
    :return: None.
    """
    resource_modules = pkgutil.iter_modules(
        resource_package.__path__, prefix=full_package_name + '.'
    )
    for module in resource_modules:
        module_members = inspect.getmembers(importlib.import_module(module.name))
        for res_name, res_obj in module_members:
            # Skip private functions and excluded resources.
            if res_name.startswith("_") or res_name in EXCLUDE_RESOURCES:
                continue

            res_node = ResourceNode(res_name, res_obj, module)

            # Skip members that are not resources or functions defined in the current
            # module.
            if not res_node.is_resource and not res_node.is_module_func:
                continue

            res_node.add_to_flowchart(flowchart)

        if module.ispkg:
            parse_resource_package(flowchart, module)


def package_doc_path(package):
    """
    Get the relative path for a package's documentation HTML.
    """
    return os.path.join(
        os.path.relpath(package.__path__[0], ROOT_PACKAGE_PATH),
        "index.html",
    )


def module_doc_path(module):
    """
    Get the relative path for a module's documentation HTML.
    """
    return re.sub(
        r"\.py$", ".html", os.path.relpath(module.__file__, ROOT_PACKAGE_PATH)
    )


def get_module_doc_link(module):
    """
    Get the full link to a module's documentation.
    """
    module_f = importlib.import_module(module.name)
    if module.ispkg:
        rel_path = package_doc_path(module_f)
    else:
        rel_path = module_doc_path(module_f)

    return os.path.join(API_REFERENCE_PATH, rel_path)


def get_function_doc_link(module, func_name):
    """
    Get the full link to a function's documentation.
    """
    return f"{get_module_doc_link(module)}#{module.name}.{func_name}"


def main():
    hl.init()
    sys.path.insert(0, REPOSITORY_ROOT_PATH)

    packages = setuptools.find_namespace_packages(
        where=REPOSITORY_ROOT_PATH, include=["gnomad_qc.*"]
    )
    top_level_packages = [pkg for pkg in packages if pkg in INCLUDE_TOP_LEVEL_PACKAGES]

    fc = QCFlowchart()

    for package_name in top_level_packages:
        package = importlib.import_module(package_name)
        for module in pkgutil.iter_modules(package.__path__):
            full_module_name = f"{package_name}.{module.name}"
            module = importlib.import_module(full_module_name)
            if full_module_name in RESOURCE_PACKAGES:
                parse_resource_package(fc, module, full_module_name)

    # TODO: Add edges to the graph.
    # NOTE: The below print statements are just to see that the resource nodes are being
    # added to the graph correctly.
    print(fc.get_graph())
    print(fc.graph.nodes.data(True))


if __name__ == "__main__":
    main()
