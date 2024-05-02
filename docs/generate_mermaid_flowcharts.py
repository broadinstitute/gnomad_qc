import ast
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
import argparse
import pickle
from collections import defaultdict, Counter

import hail as hl

from gnomad.resources.resource_utils import (
    BaseResource,
    BaseVersionedResource,
)
from gnomad.utils.file_utils import check_file_exists_raise_error

from generate_api_reference import module_flowchart_path, package_flowchart_path

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

GNOMAD_API_REFERENCE_PATH = "https://broadinstitute.github.io/gnomad_methods/api_reference"
"""Path to the gnomAD methods API reference documentation on GitHub Pages."""

INCLUDE_TOP_LEVEL_PACKAGES = {"gnomad_qc.v4"}
"""Top-level packages to include in the flowchart construction."""

MODULES_TO_SKIP = {
    "cuKING",
    "resources",
    "subset",
    "create_release_utils",
    "insilico_predictors",
    "vep_context_ht",
    "vrs_annotation_batch",
    "make_var_annot_hists",
    "create_sample_qc_metadata_ht_genomes",
    "import_variant_qc_vcf"
}
"""Modules to skip when constructing the flowchart."""

RAW_RESOURCE_FUNC = "get_gnomad_v4_vds"
EXTERNAL_RESOURCES = {"get_gnomad_v4_vds", "calling_intervals"}
EXCLUDE_RESOURCES = {
    "get_checkpoint_path",
    "get_logging_path",
    "add_meta",
    "gnomad_v4_genotypes",
}
"""Resources functions to exclude from the flowchart resource nodes."""

EXCLUDE_FUNCTIONS = {
    "get_pipeline_resources",
    "get_logging_path",
    "check_resource_existence",
    "bi_allelic_expr",
    "filter_to_autosomes",
    "annotate_adj",
    "platform_table_to_dict",
}

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
        "strat": ["under_three_alt_alleles"]
    },
    "get_downsampling": {
        "subset": [None, "non_ukb"]
    },
    "info_vcf_path": {
        "info_method": ["AS", "quasi", "set_long_AS_missing"]
    }
}
"""Resource function parameter options to include in the resource nodes."""

KWARG_OPTIONS_BY_MODULE = {
    "gnomad_qc.v4.sample_qc.hard_filters": {
        "calling_interval_name": ["intersection"],  #["broad", "ukb", "intersection"],
        "calling_interval_padding": [50],  # [0, 50]
        "n_alt_alleles_strata_name": ["three"],
    },
    "gnomad_qc.v4.sample_qc.interval_qc": {
        "calling_interval_name": ["intersection"],  #["broad", "ukb", "intersection"],
        "calling_interval_padding": [50],  # [0, 50]
    },
    "gnomad_qc.v4.sample_qc.relatedness": {
        "relatedness_method": ["cuking", "pc_relate"]
    },
    "gnomad_qc.v4.variant_qc.evaluation": {
        "model_id": ["vqsr_quasi_adj_ts_ss_calling", "rf_ba64690f"]
    },
    "gnomad_qc.v4.variant_qc.final_filter": {
        "model_id": ["vqsr_quasi_adj_ts_ss_calling", "rf_ba64690f"]
    },
    "gnomad_qc.v4.variant_qc.final_filter_genomes": {
        "model_id": ["vqsr_alleleSpecificTrans"]
    },
    "gnomad_qc.v4.variant_qc.random_forest": {
        "model_id": ["rf_ba64690f"]
    },
    "gnomad_qc.v4.sample_qc.outlier_filtering": {
        "exclude_unreleasable_samples_all_steps": [False],
        "use_nearest_neighbors_approximation": [False],
        "apply_stratified_filtering": [False],
    },
    "gnomad_qc.v4.annotations.compute_coverage": {
        "calling_interval_name": ["intersection"],  # ["broad", "ukb", "intersection"],
        "calling_interval_padding": [50],  # [0, 50]
    },
}
"""Parameter options to pass to 'get_pipeline_resources' function by module."""

NODE_COLORS = {
    "step": "#98BFFF",  # Light Blue
    "func": "#E0D5E6",  # Light Purple
    "gnomad_methods": "#F9F2CE",  # Light Yellow
    "hail": "#FAEACE",  # Light Orange
    "resource": "#D6D6D6",  # Light Gray
    "main_resource": "#d9a56c",  # Light Brown
    "raw_resource": "#FFCCCB",  # Light Red
}
"""
Types of nodes:

    - Pipeline step
    - Function call of function within script
    - Function call of function in gnomAD methods
    - Function call of function in hail
    - Resource
    - Main resource (resource used in other pipelines)
"""


########################################################################################
# Functions to get the relative path for a package's or module's documentation HTML.
########################################################################################
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


def get_module_doc_link(module_name, relative=False, gnomad_methods=False):
    """
    Get the full link to a module's documentation.
    """
    module_f = importlib.import_module(module_name)
    if gnomad_methods:
        return os.path.join(
            GNOMAD_API_REFERENCE_PATH,
            "/".join(module_name.split(".")[1:]) + ".html",
        )

    if inspect.ismodule(module_f):
        rel_path = module_doc_path(module_f)
    else:
        rel_path = package_doc_path(module_f)

    if relative:
        return rel_path

    return os.path.join(API_REFERENCE_PATH, rel_path)


def get_function_doc_link(module_name, func_name, relative=False, gnomad_methods=False):
    """
    Get the full link to a function's documentation.
    """
    return f"{get_module_doc_link(module_name, relative=relative, gnomad_methods=gnomad_methods)}#{module_name}.{func_name}"


########################################################################################
# QC Flowchart class and Base QC Parser class to help construct the flowchart.
########################################################################################
class QCFlowchart:
    """
    Represents a directed graph (flowchart) of a quality control pipeline.

    This class facilitates the construction and manipulation of a graph that
    models the workflow of a quality control pipeline, allowing for the addition,
    removal, and query of nodes and edges with detailed attributes.
    """

    def __init__(self, graph=None):
        """
        Initialize a QCFlowchart object with an optional existing graph.

        :param graph: An existing directed graph if one is to be used, otherwise None.
        """
        self.graph = nx.DiGraph(graph) if graph is not None else nx.DiGraph()

    def _format_node(self, node, node_type=None, script_name=None):
        """
        Format the node name based on its type and script name for consistency in the graph.

        :param node: The base name of the node.
        :param node_type: The type of the node (optional).
        :param script_name: The script associated with the node (optional).
        :return: A string that represents the fully formatted node name.
        """
        if script_name:
            node = f"{script_name}::{node}"
        if node_type:
            node = f"{node_type}::{node}"
        return node

    def add_node(
        self,
        node,
        node_type,
        node_id,
        label,
        module_name,
        href=None,
        script_name=None,
        **kwargs
    ):
        """
        Add a node to the graph with specified attributes.

        :param node: The base name of the node.
        :param node_type: The type of the node.
        :param node_id: A unique identifier for the node.
        :param label: A human-readable label for the node.
        :param module_name: The name of the module where the node is defined.
        :param href: An optional hyperlink for the node.
        :param script_name: An optional script name associated with the node.
        :param kwargs: Additional keyword arguments for node attributes.
        """
        node_name = self._format_node(node, node_type, script_name)
        node_attr = {
            "name": node,
            "type": node_type,
            "id": node_id,
            "label": label,
            "module": module_name,
            "href": href or "",
            "script_name": script_name,
            **kwargs
        }
        self.graph.add_node(node_name, **node_attr)

    def add_edge(
        self,
        from_node,
        to_node,
        from_node_type=None,
        to_node_type=None,
        from_script_name=None,
        to_script_name=None,
        **kwargs
    ):
        """
        Add an edge between two nodes, optionally adding the nodes if they do not exist.

        :param from_node: The originating node's base name.
        :param to_node: The destination node's base name.
        :param from_node_type: The type of the originating node (optional).
        :param to_node_type: The type of the destination node (optional).
        :param from_script_name: The script associated with the originating node
            (optional).
        :param to_script_name: The script associated with the destination node
            (optional).
        :param kwargs: Additional keyword arguments for edge attributes.
        """
        from_node_name = self._format_node(from_node, from_node_type, from_script_name)
        to_node_name = self._format_node(to_node, to_node_type, to_script_name)

        if not self.graph.has_node(from_node_name):
            self.add_node(from_node, from_node_type, f"{from_node_type}_{from_node}", from_node, "", from_script_name)

        if not self.graph.has_node(to_node_name):
            self.add_node(to_node, to_node_type, f"{to_node_type}_{to_node}", to_node, "", to_script_name)

        self.graph.add_edge(from_node_name, to_node_name, **kwargs)

    def remove_node(self, node, node_type=None, script_name=None):
        """
        Remove a node from the graph and automatically reconnect its predecessors to its successors.

        :param node: The base name of the node.
        :param node_type: The type of the node (optional).
        :param script_name: The script associated with the node (optional).
        """
        node_name = self._format_node(node, node_type, script_name)
        predecessors = list(self.graph.predecessors(node_name))
        successors = list(self.graph.successors(node_name))
        self.graph.remove_node(node_name)

        for pred in predecessors:
            for succ in successors:
                self.graph.add_edge(pred, succ)

    def neighbors(self, node, node_type=None, script_name=None):
        """
        Get the neighbors of a node within the graph.

        :param node: The base name of the node.
        :param node_type: The type of the node (optional).
        :param script_name: The script associated with the node (optional).
        :return: A list of neighboring node names.
        """
        node_name = self._format_node(node, node_type, script_name)
        return list(self.graph.neighbors(node_name))

    def get_graph(self):
        """
        Get the internal graph object.

        :return: The internal NetworkX DiGraph object.
        """
        return self.graph

    def set_node_attributes(self, node, attrs):
        """
        Set attributes for a specific node in the graph.

        :param node: The node for which to set attributes.
        :param attrs: A dictionary of attributes to set.
        """
        nx.set_node_attributes(self.graph, {node: attrs})


class BaseQCParser:
    """
    Base class for parsing quality control pipeline resources and functions.

    This parser helps to manage and define how resources and functions are
    integrated and visualized in a QC Flowchart.
    """

    def __init__(self, name, obj, module):
        """
        Initialize the base QC parser.

        :param name: Name of the resource or function.
        :param obj: The actual resource or function object.
        :param module: Module where the resource or function is defined.
        """
        self.name = name
        self.obj = obj
        self.module = module
        self.is_resource = self.get_is_resource()
        self.is_module_func = self.get_is_module_func()
        self.type = None

    def get_is_resource(self):
        """
        Determine if the object is a subclass of BaseResource or BaseVersionedResource.

        :return: True if the object is a resource, otherwise False.
        """
        return issubclass(type(self.obj), (BaseResource, BaseVersionedResource))

    def get_is_module_func(self):
        """
        Check if the object is a function defined within the same module.

        :return: True if the object is a function defined in the current module,
            otherwise False.
        """
        return inspect.isfunction(self.obj) and self.obj.__module__ == self.module.name

    def get_signature(self):
        """
        Retrieve the signature of the function.

        :return: The signature object of the function.
        """
        return inspect.signature(self.obj)

    def get_resource_kwargs(self, kwargs_options=None):
        """
        Construct a dictionary of keyword arguments for PipelineStepNode function from the signature.

        :param kwargs_options: Dictionary to store and pass keyword arguments.
        :return: Dictionary with possible keyword arguments and their values.
        """
        if kwargs_options is None:
            kwargs_options = {}

        signature = self.get_signature()
        for k, v in signature.parameters.items():
            if k in ["test", "overwrite"]:
                kwargs_options[k] = [False]
            elif k not in kwargs_options:
                if v.annotation == bool:
                    kwargs_options[k] = [True, False]
                elif k == "data_type":
                    kwargs_options[k] = ["exomes", "genomes", "joint"]

        return kwargs_options

    def get_resource_kwargs_combos(self, kwargs_options=None):
        """
        Generate all combinations of keyword arguments for the resource function.

        :param kwargs_options: Predefined options for keyword arguments.
        :return: List of all combinations of keyword arguments.
        """
        if kwargs_options is None:
            kwargs_options = {}

        kwargs_options = self.get_resource_kwargs(kwargs_options)
        all_options = []
        for k, v in kwargs_options.items():
            all_options.append([(k, val) for val in v])

        return list(itertools.product(*all_options))

    def remove_excluded_parameters(self, bind_kwargs):
        """
        Filter out excluded parameters from the function's keyword arguments.

        :param bind_kwargs: Dictionary of bound keyword arguments.
        :return: Dictionary with excluded parameters removed.
        """
        if bind_kwargs is None:
            bind_kwargs = {}
        return {k: v for k, v in bind_kwargs.items() if k not in EXCLUDE_PARAMS}

    def get_resource_display_name(self, bind_kwargs=None):
        """
        Generate a display name for the resource node using the provided keyword arguments.

        :param bind_kwargs: Dictionary of resource function keyword arguments.
        :return: Formatted display name for the resource node.
        """
        bind_kwargs = self.remove_excluded_parameters(
            bind_kwargs) if bind_kwargs else {}
        parts = [f"{k}={v}" for k, v in bind_kwargs.items()]
        res_display = "\n".join(parts) if parts else ""
        return f"{self.name}({res_display})"

    def get_resource_id(self, bind_kwargs=None):
        """
        Generate an identifier for the resource node.

        :param bind_kwargs: Dictionary of resource function keyword arguments.
        :return: Unique identifier for the resource node.
        """
        bind_kwargs = self.remove_excluded_parameters(
            bind_kwargs) if bind_kwargs else {}
        parts = [f"{k}_{v}" for k, v in bind_kwargs.items()]
        name_reformat = self.name.replace('--', '').replace('-', '_')
        res_arg_str = "_".join(parts) if parts else ""
        return f"{self.type}_{name_reformat}{('_' + res_arg_str if res_arg_str else '')}"

    def get_doc_link(self, relative=True, add_function=True):
        """
        Construct a link to the documentation of the resource or module.

        :param relative: Flag to indicate if the link should be relative.
        :param add_function: Flag to include function details in the link.
        :return: URL string to the documentation.
        """
        if add_function:
            return get_function_doc_link(self.module.name, self.name, relative=relative)
        else:
            return get_module_doc_link(self.module.name, relative=relative)

    def get_resource_paths(self, resource):
        """
        Retrieve paths associated with a resource.

        :param resource: The resource object.
        :return: List of paths if available, None otherwise.
        """
        paths = []
        if hasattr(resource, "versions"):
            paths = [r.path for r in resource.versions.values()]
        elif hasattr(resource, "path"):
            paths = [resource.path]
        elif isinstance(resource, str):
            paths = [resource]
        else:
            return None

        return [p for p in paths if
                check_file_exists_raise_error(p, error_if_exists=False)]

    def add_resource_to_flowchart(self, flowchart, resource, bind_kwargs=None):
        """
        Add a resource and its related paths to a flowchart.

        :param flowchart: The flowchart object to add the resource to.
        :param resource: The resource function or object.
        :param bind_kwargs: Optional keyword arguments for resource function.
        """
        paths = self.get_resource_paths(resource)
        if paths is None:
            err_msg = f" with {bind_kwargs}" if bind_kwargs else ""
            raise ValueError(
                f"Resource {self.obj}{err_msg} did not return a valid path.")

        is_raw_resource = self.name in EXTERNAL_RESOURCES
        for p in paths:
            flowchart.add_node(
                p,
                node_type=self.type,
                node_id=self.get_resource_id(bind_kwargs),
                label=self.get_resource_display_name(bind_kwargs),
                module_name=self.module.name,
                href=self.get_doc_link(),
                is_raw_resource=is_raw_resource
            )

    def add_to_flowchart(self, flowchart):
        """
        Iterate over function call combinations and add the resulting resources to the flowchart.

        :param flowchart: The flowchart object to update.
        """
        for bind_kwargs, func_call in self.iter_func_call_combos():
            self.add_resource_to_flowchart(flowchart, func_call,
                                           bind_kwargs=bind_kwargs)


########################################################################################
# Resource and PipelineStepNode parsers to help construct the flowchart.
########################################################################################
class ResourceParser(BaseQCParser):
    """
    Specialized parser for handling resource functions within a quality control pipeline.

    This parser is designed to manage how resource functions are parsed, visualized,
    and added to a QC Flowchart.
    """

    def __init__(self, res_name, res_obj, module):
        """
        Initialize a ResourceParser object with a resource function.

        :param res_name: The name of the resource function.
        :param res_obj: The resource function object.
        :param module: The module where the resource function is defined.
        """
        super().__init__(res_name, res_obj, module)
        self.type = "resource"

    @staticmethod
    def remove_excluded_parameters(bind_kwargs):
        """
        Filter out excluded parameters from the resource function keyword arguments.

        :param bind_kwargs: Dictionary of resource function keyword arguments.
        :return: Dictionary with excluded parameters removed.
        """
        if bind_kwargs is None:
            bind_kwargs = {}
        return {
            k: v for k, v in bind_kwargs.items() if k not in EXCLUDE_RESOURCE_PARAMS
        }

    def get_resource_kwargs_combos(self, kwargs_options=None):
        """
        Generate all combinations of keyword arguments for the resource function,
        tailored for specific resources.

        :param kwargs_options: Predefined options for keyword arguments.
        :return: List of all combinations of keyword arguments.
        """
        if kwargs_options is None:
            kwargs_options = KWARG_OPTIONS_BY_RESOURCE.get(self.name, {})

        return super().get_resource_kwargs_combos(kwargs_options)

    def add_to_flowchart(self, flowchart):
        """
        Add the resource or function to a flowchart, considering its type and specifications.

        :param flowchart: The QCFlowchart object to which the resource will be added.
        """
        if self.is_module_func and self.name != RAW_RESOURCE_FUNC:
            super().add_to_flowchart(flowchart)
        else:
            # Handle cases where the resource is a raw resource or not a module function
            super().add_resource_to_flowchart(flowchart, self.obj, bind_kwargs=None)


class PipelineResourceParser(BaseQCParser):
    """
    Parser for handling and integrating pipeline steps into a quality control pipeline flowchart.

    This parser is specifically designed for extracting and organizing the workflow of
    pipeline steps based on their input and output resources, and visualizing them
    within a flowchart.
    """

    def __init__(self, name, obj, module):
        """
        Initialize a PipelineResourceParser object.

        :param name: The name of the pipeline step.
        :param obj: The object representing the pipeline step.
        :param module: The module where the pipeline step is defined.
        """
        super().__init__(name, obj, module)
        self.script_name = module.name.split(".")[-1] + ".py"
        self.type = "step"

    def get_resource_kwargs(self, kwargs_options=None):
        """
        Retrieve a dictionary of keyword arguments for configuring the pipeline step function.

        :param kwargs_options: Predefined options for keyword arguments.
        :return: Dictionary with keyword arguments for the pipeline step function.
        """
        if kwargs_options is None:
            kwargs_options = KWARG_OPTIONS_BY_MODULE.get(self.module.name, {})

        return super().get_resource_kwargs(kwargs_options)

    def add_step_to_flowchart(self, flowchart, step_obj):
        """
        Add a pipeline step and its associated resources to the flowchart.

        :param flowchart: The flowchart object to update.
        :param step_obj: The pipeline step object containing input and output resources.
        """
        step_name = step_obj.pipeline_step
        step_id = step_name.replace("--", "").replace("-", "_")
        output_paths = self._gather_output_paths(step_obj)
        edges = self._prepare_edges(flowchart, step_obj, output_paths, step_id)

        if edges:
            self._add_step_node(flowchart, step_id, step_name, edges)

    def _gather_output_paths(self, step_obj):
        output_paths = []
        for resources in step_obj.output_resources.values():
            for resource in resources:
                output_paths.extend(self.get_resource_paths(resource) or [])
        return output_paths

    def _prepare_edges(self, flowchart, step_obj, output_paths, step_id):
        edges = []
        step_created = {}
        steps_using = defaultdict(set)

        for resources in step_obj.input_resources.values():
            for resource in resources:
                paths = self.get_resource_paths(resource) or []
                for p in paths:
                    if not flowchart.has_edge(p, step_id, "resource", "step", to_script_name=self.script_name):
                        steps_using[f"resource::{p}"].add(f"step::{self.script_name}::{step_id}")
                        edges.append((p, step_id, "resource", "step", None, self.script_name))

        for p in output_paths:
            step_created[f"resource::{p}"] = f"step::{self.script_name}::{step_id}"
            if not flowchart.has_edge(step_id, p, "step", "resource", from_script_name=self.script_name):
                edges.append((step_id, p, "step", "resource", self.script_name, None))

        return edges

    def _add_step_node(self, flowchart, step_id, step_name, edges):
        flowchart.add_node(
            step_id,
            node_type=self.type,
            node_id=f"step_{step_id}",
            label=step_name,
            module_name=self.module.name,
            script_name=self.script_name,
        )

        for e in edges:
            flowchart.add_edge(*e)

    def add_to_flowchart(self, flowchart):
        """
        Iterate over function call combinations and integrate each pipeline step into the flowchart.

        :param flowchart: The flowchart object to update.
        """
        for bind_kwargs, func_call in self.iter_func_call_combos():
            for step_obj in func_call.pipeline_steps.values():
                self.add_step_to_flowchart(flowchart, step_obj)


########################################################################################
# ScriptMainParser and CallVisitor classes to parse the main function of a script.
########################################################################################
class CallVisitor(ast.NodeVisitor):
    """
    Custom AST NodeVisitor that captures the name and value from calls in the AST.

    It builds the fully qualified name of functions or attributes accessed in the AST and
    captures constant values if encountered.
    """

    def __init__(self):
        """
        Initialize the CallVisitor with empty name components and no constant value.
        """
        self._name = []
        self._const_value = None

    @property
    def name(self):
        """
        Return the fully qualified name built from visited nodes, removing any 'args.' substrings.
        """
        return '.'.join(reversed(self._name)).replace("args.", "")

    @name.deleter
    def name(self):
        """
        Reset the name components.
        """
        self._name = []

    @property
    def const_value(self):
        """
        Return the constant value captured from visited nodes.
        """
        return self._const_value

    @const_value.deleter
    def const_value(self):
        """
        Reset the constant value.
        """
        self._const_value = None

    def visit_Name(self, node):
        """
        Visit a Name node and append its identifier to the name components.
        """
        self._name.append(node.id)

    def visit_Attribute(self, node):
        """
        Visit an Attribute node and append its identifier and the value identifier to the name components.
        """
        self._name.append(node.attr)
        self.generic_visit(node)  # Recursively visit the value to handle nested attributes.

    def visit_Constant(self, node):
        """
        Visit a Constant node and capture its value.
        """
        self._const_value = node.value

    def visit_keyword(self, node):
        """
        Visit a keyword node; this method is necessary to ensure all parts of an expression are visited.
        """
        self.generic_visit(node)


class ScriptMainParser(BaseQCParser):
    """
    Parser for extracting and processing function calls within the main function of a script,
    integrating these calls into a QC flowchart.
    """

    def __init__(self, name, obj, module):
        super().__init__(name, obj, module)
        self.script_name = module.name.split(".")[-1] + ".py"
        self.type = "func"
        self.steps = []
        self.step_func_calls = defaultdict(list)
        self.step_raw_input = {}
        self.function_to_module = self.get_module_function_call_dict()

    def get_module_function_call_dict(self):
        """
        Build a dictionary mapping function names to their module paths and names.

        :return: Dictionary of function names to (module path, function name).
        """
        module_f = importlib.import_module(self.module.name)
        return {name: (obj.__module__, obj.__name__) for name, obj in inspect.getmembers(module_f, inspect.isfunction)}

    def parse_ast_node(self, node, flowchart):
        """
        Recursively parse AST nodes, identifying function calls and structuring them for flowchart integration.

        :param node: The root AST node to start parsing from.
        :param flowchart: The QCFlowchart object to which the parsed information will be added.
        """
        for child in ast.iter_child_nodes(node):
            self._handle_ast_node(child, flowchart)

    def _handle_ast_node(self, child, flowchart):
        """
        Handle an individual AST node based on its type.

        :param child: The AST node being handled.
        :param flowchart: The QCFlowchart object being updated.
        """
        if isinstance(child, ast.If):
            self._handle_if_node(child, flowchart)
        elif isinstance(child, ast.Call):
            self._handle_call_node(child, flowchart)
        self.parse_ast_node(child, flowchart)

    def _handle_if_node(self, node, flowchart):
        """
        Handle 'if' statements in the AST, particularly for determining control flow steps.

        :param node: The 'if' node in the AST.
        :param flowchart: The flowchart object being updated.
        """
        callvisitor = CallVisitor()
        callvisitor.visit(node.test)
        if flowchart.has_node(callvisitor.name, "step", self.script_name):
            node.step = callvisitor.name
            self.steps.append(callvisitor.name)

    def _handle_call_node(self, node, flowchart):
        """
        Handle 'call' nodes in the AST, extracting and processing function calls.

        :param node: The 'call' node in the AST.
        :param flowchart: The flowchart object being updated.
        """
        callvisitor = CallVisitor()
        callvisitor.visit(node.func)
        func_module, func_name = self.function_to_module.get(callvisitor.name, ("", ""))
        pkgs = func_module.split(".")
        in_resource_pkg = len(pkgs) > 2 and pkgs[2] == "resources"

        self._process_function_call(node, func_name, func_module, in_resource_pkg, flowchart)

    def _process_function_call(self, node, func_name, func_module, in_resource_pkg, flowchart):
        """
        Process a function call found in the AST, updating the flowchart as necessary.

        :param node: The node representing the function call.
        :param func_name: The name of the function being called.
        :param func_module: The module where the function is located.
        :param in_resource_pkg: Boolean indicating if the function is in a resource package.
        :param flowchart: The flowchart object being updated.
        """
        # Process function call details here (similar to existing complex logic in original method)
        pass

    def parse_main_ast(self, flowchart):
        """
        Parse the main function of the script and update the flowchart based on function calls.

        :param flowchart: The flowchart object being updated.
        """
        main_ast_node = ast.parse(inspect.getsource(self.obj))
        self.parse_ast_node(main_ast_node, flowchart)

    def add_to_flowchart(self, flowchart):
        """
        Add the main function's details and associated function calls to the flowchart.

        :param flowchart: The flowchart object to update.
        """
        self.parse_main_ast(flowchart)
        # Additional processing based on parsed AST nodes


########################################################################################
# Main function to parse the repo and construct the flowchart.
########################################################################################
def parse_qc_package(flowchart, package, full_package_name, resource_pkg=False):
    """
    Parse a resource package and add its resource nodes to the QCFlowchart.

    :param flowchart: QCFlowchart object.
    :param package: Package object for resources.
    :param full_package_name: Full package name.
    :param resource_pkg: Boolean indicating if the package is a resource package.
    :return: None.
    """
    resource_pkg = is_resource_package(full_package_name, resource_pkg)
    for module in iter_package_modules(package, full_package_name):
        parse_module(flowchart, module, resource_pkg)


def is_resource_package(full_package_name, resource_pkg):
    """
    Determine if the given package is a resource package.

    :param full_package_name: Full name of the package.
    :param resource_pkg: Current state if it's a resource package.
    :return: True if it's a resource package, else False.
    """
    return full_package_name in RESOURCE_PACKAGES or resource_pkg


def iter_package_modules(package, full_package_name):
    """
    Generate an iterator over the modules in the given package.

    :param package: Package object.
    :param full_package_name: Full name of the package, used as prefix.
    :return: Iterator over package modules.
    """
    return pkgutil.iter_modules(package.__path__, prefix=full_package_name + '.')


def parse_module(flowchart, module, resource_pkg):
    """
    Parse a single module and manage the parsing based on package type.

    :param flowchart: QCFlowchart object.
    :param module: Module object to parse.
    :param resource_pkg: Boolean indicating if it is a resource package.
    :return: None.
    """
    module_f = importlib.import_module(module.name)
    if module.ispkg:
        parse_qc_package(flowchart, module_f, module.name, resource_pkg=resource_pkg)
        return

    rel_module_name = module.name.split(".")[-1]
    if not resource_pkg and rel_module_name in MODULES_TO_SKIP:
        return

    process_module_members(flowchart, module_f, module, resource_pkg)


def process_module_members(flowchart, module_f, module, resource_pkg):
    """
    Process all members of a given module to extract resources, pipeline resources, or main scripts.

    :param flowchart: QCFlowchart object.
    :param module_f: Module from which members are inspected.
    :param module: Module object, needed for context in parsers.
    :param resource_pkg: Boolean indicating if it is a resource package.
    :return: None.
    """
    resource_parsers = []
    pipeline_parsers = []
    main_parsers = []

    for name, obj in inspect.getmembers(module_f):
        if resource_pkg:
            res_parse = handle_resource_package_member(module, name, obj)
            resource_parsers.extend([] if res_parse is None else [res_parse])
        elif name == "get_pipeline_resources":
            pl_parse = handle_pipeline_resource(module, name, obj)
            pipeline_parsers.extend([] if pl_parse is None else [pl_parse])
        elif name == "main":
            main_parsers.append(ScriptMainParser(name, obj, module))

    for p in resource_parsers + pipeline_parsers + main_parsers:
        p.add_to_flowchart(flowchart)


def handle_resource_package_member(module, name, obj):
    """
    Handle parsing of a resource member within a resource package.

    :param module: Module containing the member.
    :param name: Name of the member.
    :param obj: Object (typically a function or class) associated with the name.
    :return: None.
    """
    if not (name.startswith("_") or name in EXCLUDE_RESOURCES):
        res_parser = ResourceParser(name, obj, module)
        if res_parser.is_resource or res_parser.is_module_func:
            return res_parser


def handle_pipeline_resource(module, name, obj):
    """
    Handle parsing of pipeline resource functions within a module.

    :param module: Module containing the member.
    :param name: Name of the pipeline resource function.
    :param obj: Object (function) associated with the name.
    :return: None.
    """
    pipeline_parser = PipelineResourceParser(name, obj, module)
    if pipeline_parser.is_module_func:
        return pipeline_parser


def iter_top_level_package_modules():
    """
    Generate an iterator over the top-level package names in the gnomAD QC pipeline.

    This function finds all namespace packages under the specified repository root path
    that match the pattern "gnomad_qc.*" and filters them to include only those specified
    in INCLUDE_TOP_LEVEL_PACKAGES.

    :yield: Package names that are considered top-level for the QC pipeline.
    """
    # Find all potential top-level packages that match the "gnomad_qc.*" pattern.
    candidate_packages = setuptools.find_namespace_packages(
        where=REPOSITORY_ROOT_PATH, include=["gnomad_qc.*"]
    )

    # Filter the packages to include only those defined as top-level.
    top_level_packages = [pkg for pkg in candidate_packages if pkg in INCLUDE_TOP_LEVEL_PACKAGES]

    for package_name in top_level_packages:
        yield package_name


########################################################################################
# Functions to format the flowchart nodes and connections for mermaid.js.
########################################################################################
def format_node(node_identifier, node_label, node_type, href=None):
    """
    Format a node string for visual representation in a graph, optionally adding a hyperlink.

    :param node_identifier: The identifier of the node, used in the graph.
    :param node_label: The label of the node to be displayed.
    :param node_type: The type of the node which determines the node's visual style.
    :param href: Optional hyperlink to be embedded in the node label.
    :return: A formatted node string.
    """
    # Clean up node identifier
    node_identifier = node_identifier.strip("_")

    # Format the hyperlink if provided
    if href:
        link_type = "a href" if href.startswith("http") else "a class='reference internal' href"
        node_label = f"<{link_type}='{href}'>{node_label}</{link_type.split()[0]}>"

    # Format the node based on its type
    if node_type == "step":
        formatted_node = f"  {node_identifier}{{{{\"{node_label}\"}}}}:::{node_type}_color;"
    elif node_type in {"func", "gnomad_methods"}:
        formatted_node = f"  {node_identifier}[[\"{node_label}\"]]:::{node_type}_color;"
    elif node_type in {"resource", "main_resource", "raw_resource"}:
        formatted_node = f"  {node_identifier}[/\"{node_label}\"/]:::{node_type}_color;"
    else:
        raise ValueError(f"Node type '{node_type}' not recognized.")

    return formatted_node


def format_connection(node1, node2, invisible=False):
    """
    Format a graph connection between two nodes with an option to make the connection invisible.

    :param node1: The identifier of the first node.
    :param node2: The identifier of the second node.
    :param invisible: A boolean indicating if the connection should be rendered as invisible.
    :return: A formatted string representing the connection between two nodes.
    """
    # Remove any leading/trailing underscores from node identifiers for clean formatting
    node1 = node1.strip("_")
    node2 = node2.strip("_")

    # Choose the connector based on the visibility desired
    connector = "~~~" if invisible else "-->"

    # Construct the formatted connection string
    return f"  {node1} {connector} {node2};"


def networkx_sg_to_mermaid(graph, subgraph, subgraph_name, href_relative_to=None):
    """
    Converts a NetworkX subgraph into groups of Mermaid-formatted strings based on node and edge types.

    :param graph: The complete NetworkX graph.
    :param subgraph: The subgraph to convert.
    :param subgraph_name: The name used for subgraph identity in diagrams.
    :param href_relative_to: Base path to calculate relative hrefs for nodes.
    :return: A dictionary grouping Mermaid-formatted strings by their types (e.g., nodes, edges).
    """
    mermaid_by_group = defaultdict(list)
    resource_types = {"raw_resource", "main_resource", "resource"}

    for node, data in subgraph.nodes(data=True):
        if "id" not in data or subgraph.degree(node) == 0:
            continue

        href = data.get("href")
        if href and not href.startswith("http") and href_relative_to:
            base_path = href_relative_to.split("v4/")[-1]
            href = os.path.relpath(href, base_path)

        node_format = format_node(data["id"], data["label"], data["type"], href)

        if data["type"] in resource_types:
            group_key = "raw_input_resources" if data["type"] == "raw_resource" else "output_resources"
            if data["type"] == "main_resource" and len(list(graph.neighbors(node))) == 0:
                mermaid_by_group[group_key].append(node_format)
        else:
            mermaid_by_group["nodes"].append(node_format)

    for from_node, to_node in subgraph.edges(data=False):
        from_node_data, to_node_data = subgraph.nodes[from_node], subgraph.nodes[to_node]

        if "id" not in from_node_data or "id" not in to_node_data:
            continue

        if from_node_data["type"] in resource_types:
            con_format = format_connection(from_node_data["id"], subgraph_name)
            mermaid_by_group["resource_edges"].append(con_format)
        elif to_node_data["type"] not in resource_types:
            con_format = format_connection(from_node_data["id"], to_node_data["id"])
            mermaid_by_group["edges"].append(con_format)

    return mermaid_by_group


def networkx_to_mermaid(graph, href_relative_to=None, add_module_subgraphs=False,
                        step_sub_graphs=None):
    """
    Converts a NetworkX graph into a dictionary of Mermaid formatted strings grouped by node and edge types.

    :param graph: NetworkX graph to convert.
    :param href_relative_to: Base directory for relative hyperlinks.
    :param add_module_subgraphs: Boolean indicating whether to add module subgraphs.
    :param step_sub_graphs: Dictionary of subgraphs for steps.
    :return: Dictionary of lists containing Mermaid formatted strings.
    """
    mermaid_by_group = defaultdict(list)
    resource_types = {"raw_resource", "main_resource", "resource"}

    # Process nodes
    for node, data in graph.nodes(data=True):
        if "id" not in data or graph.degree(node) == 0:
            continue

        href = data.get("href")
        if href and not href.startswith("http") and href_relative_to:
            base_path = href_relative_to.split("v4/")[-1]
            href = os.path.relpath(href, base_path)

        node_format = format_node(data["id"], data["label"], data["type"], href)
        if data["type"] in resource_types:
            group_key = f"{data['type']}_resources"
            if data["type"] == "main_resource" and len(
                    list(graph.predecessors(node))) == 0:
                group_key = "input_resources"
            mermaid_by_group[group_key].append(node_format)
        elif data["type"] == "step" and step_sub_graphs and data[
            "id"] in step_sub_graphs:
            sg = step_sub_graphs[data["id"]]
            step_by_group = networkx_sg_to_mermaid(graph, sg, data["id"],
                                                   href_relative_to=href_relative_to)
            mermaid_by_group[data["id"]].append(
                format_node("n_" + data["id"], data["label"], data["type"], href))
            mermaid_by_group[data["id"]].extend(
                step_by_group["nodes"] + step_by_group["edges"])
            mermaid_by_group["edges"].extend(step_by_group["resource_edges"])
            for res_type in resource_types:
                res_group_key = f"{res_type}_resources"
                mermaid_by_group[res_group_key].extend(step_by_group[res_group_key])
        else:
            mermaid_by_group["nodes"].append(node_format)

    # Process edges
    for from_node, to_node in graph.edges(data=False):
        from_data = graph.nodes[from_node]
        to_data = graph.nodes[to_node]

        if "id" not in from_data or "id" not in to_data:
            continue
        if from_data["type"] == "step" and step_sub_graphs and from_data[
            "id"] in step_sub_graphs:
            continue

        if add_module_subgraphs and from_data["type"] == "step" and from_data.get(
                "module"):
            module = from_data["module"].split(".")[-1]
            mermaid_by_group[module].append(
                format_connection(from_data["id"], to_data["id"]))
        elif add_module_subgraphs and to_data["type"] == "step" and to_data.get(
                "module"):
            module = to_data["module"].split(".")[-1]
            mermaid_by_group[module].append(
                format_connection(from_data["id"], to_data["id"]))
        else:
            mermaid_by_group["edges"].append(
                format_connection(from_data["id"], to_data["id"]))

    return mermaid_by_group


def traversal(G, n, start=None, seen=set([]), int_nodes=defaultdict(list), pred=None):
    if start is None:
        start = n
    next_pred = None
    for c in G.neighbors(n):
        if c == start or (c not in seen and not c.startswith("resource::")):
            if c == start:
                next_pred = n
            else:
                if pred is not None and pred != n:
                    int_nodes[pred].append(c)
            seen.add(c)
            traversal(G, c, start=start, seen=seen, int_nodes=int_nodes, pred=next_pred)

    return int_nodes


def write_mermaid_flowchart(flowchart, path, step_sub_graphs=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    edges_to_remove = []

    for n, data in dict(flowchart.nodes(data=True)).items():
        if data.get("type", "") == "main_resource":
            try:
                nx.find_cycle(flowchart, n, orientation="original")
                from_nodes = set(flowchart.predecessors(n))
                to_nodes = set(flowchart.neighbors(n))
                all_nodes = from_nodes | to_nodes

                all_cycles = nx.recursive_simple_cycles(flowchart)
                in_all_cycles = True
                for c in all_cycles:
                    if n in c:
                        all_nodes -= set(c)
                    else:
                        in_all_cycles = False
                if in_all_cycles:
                    edges_to_change = traversal(flowchart, all_nodes.pop(), n)
                    for i, (pred, descs) in enumerate(edges_to_change.items()):
                        new_data = data.copy()
                        new_data["type"] = "resource"
                        new_data["id"] = f"{new_data['id']}_int{i}"
                        new_data["label"] = f"intermediate {new_data['label']}"
                        new_n = f"{n}::intermediate{i}"
                        flowchart.add_node(new_n, **new_data)
                        flowchart.remove_edge(pred, n)
                        flowchart.add_edge(pred, new_n)
                        for desc in descs:
                            flowchart.remove_edge(n, desc)
                            flowchart.add_edge(new_n, desc)

            except nx.exception.NetworkXNoCycle:
                pass

    resource_types = {"raw_resource", "main_resource", "resource"}

    if step_sub_graphs is not None:
        flowchart = QCFlowchart(flowchart)
        for step, sg in step_sub_graphs.items():
            for from_node, to_node in sg.edges(data=False):
                from_node_data = sg.nodes[from_node]
                to_node_data = sg.nodes[to_node]
                if "id" not in from_node_data or "id" not in to_node_data:
                    continue

                if from_node_data.get("type") not in {"step"} and to_node_data.get("type") not in resource_types:
                    edges_to_remove.append((from_node, to_node))

        for from_node, to_node in edges_to_remove:
            flowchart.remove_edge(from_node, to_node)

        flowchart = flowchart.graph

    with open(path, "w") as f:
        f.write("flowchart-elk TD;\n")
        for node_type, color in NODE_COLORS.items():
            f.write(
                f"""  classDef {node_type}_color fill:{color},color:#000000\n"""
            )
        mermaid_by_group = networkx_to_mermaid(
            flowchart,
            href_relative_to=os.path.relpath(path, FLOWCHART_DIRECTORY),
            add_module_subgraphs=False,
            step_sub_graphs=step_sub_graphs,
        )

        raw_inputs = mermaid_by_group.pop("raw_input_resources", [])
        inputs = mermaid_by_group.pop("input_resources", [])
        outputs = mermaid_by_group.pop("output_resources", [])
        nodes = mermaid_by_group.pop("nodes", [])
        edges = mermaid_by_group.pop("edges", [])

        if raw_inputs and inputs:
            f.write("  subgraph AllInputResources [Input resources]\n")
            raw_in_str = "External/raw data"
            in_str = "From other modules"
        else:
            raw_in_str = "External/raw data input resources"
            in_str = "Input resources from other modules"

        if raw_inputs:
            f.write(f"  subgraph RawInputResources [{raw_in_str}]\n")
            for n in raw_inputs:
                f.write(f"  {n}\n")
            f.write("  end\n")

        if inputs:
            f.write(f"  subgraph InputResources [{in_str}]\n")
            for n in inputs:
                f.write(f"  {n}\n")
            f.write("  end\n")

        if raw_inputs and inputs:
            f.write("  end\n")

        if outputs:
            f.write("  subgraph OutputResources [Output resources]\n")
            for n in outputs:
                f.write(f"  {n}\n")
            f.write("  end\n")

        for l in nodes:
            f.write(l + "\n")

        for module, mod_nodes in mermaid_by_group.items():
            f.write(f'  subgraph {module} [" "]\n')
            for n in mod_nodes:
                f.write(f"  {n}\n")
            f.write("  end\n")

        for l in edges:
            f.write(l + "\n")


def divide_graph_by_module(graph):
    """
    Divides a graph into subgraphs based on modules, removing edges between resources and steps
    of different modules, and grouping nodes by their associated steps and modules.

    :param graph: A directed graph where nodes represent resources or steps and edges the connections.
    :return: A dictionary mapping each module to its corresponding subgraphs and their step-based divisions.
    """
    # Definitions
    resource_types = {"raw_resource", "main_resource", "resource"}

    # Prepare the graph by removing certain edges
    g_no_input_edges = graph.copy()
    resource_creation_step = dict(
        g_no_input_edges.nodes(data="step_created", default=""))
    node_module = dict(g_no_input_edges.nodes(data="module", default=""))

    # Remove edges between resources and steps from different modules
    g_no_input_edges.remove_edges_from([
        e for e in g_no_input_edges.edges if (
                e[0].startswith("resource::") and e[1].startswith("step::") and
                node_module.get(resource_creation_step[e[0]], "") != node_module.get(
            e[1], "")
        )
    ])

    # Find all connected components
    by_module = defaultdict(list)
    all_sub_graphs = list(nx.connected_components(g_no_input_edges.to_undirected()))

    # Process each component
    for s in all_sub_graphs:
        s_g = g_no_input_edges.subgraph(s)
        n_types = Counter(s_g.nodes(data="type", default=""))

        # Ignore single resource nodes
        if len(s) == 1 and n_types.get("resource", 0) == 1:
            continue

        # Determine modules not categorized under "resources"
        non_resource_modules = [
            m for k, m in s_g.nodes(data="module", default="")
            if m != "" and ((len(m.split(".")) < 3) or (m.split(".")[2] != "resources"))
        ]

        # Group subgraphs by module
        if non_resource_modules:
            by_module[set(non_resource_modules).pop()].append(s_g)

    # Assemble subgraphs by module
    sub_graphs_by_module = {}
    for module, s_g_list in by_module.items():
        module_nodes = set()
        for g in s_g_list:
            module_nodes.update(g.nodes())
            module_nodes.update(n2 for n in g.nodes() for n2 in graph.predecessors(n) if
                                graph.nodes[n2]["type"] in resource_types)
            module_nodes.update(n2 for n in g.nodes() for n2 in graph.neighbors(n) if
                                graph.nodes[n2]["type"] in resource_types)

        # Create and clean subgraph for the module
        module_sg = graph.subgraph(module_nodes).copy()
        sub_graphs_by_module[module] = {"module_subgraph": module_sg}
        module_sg.remove_edges_from([
            e for e in module_sg.edges if (
                    "resource::" in e[0] or "resource::" in e[1] or
                    (not any(
                        prefix in e[0] for prefix in ("func::", "gnomad_methods::")) and
                     not any(
                         prefix in e[1] for prefix in ("func::", "gnomad_methods::")))
            )
        ])

        # Extract connected components by steps
        sub_graph_by_step = {}
        for c in nx.connected_components(module_sg.to_undirected()):
            if len(c) > 1:
                s_g2 = module_sg.subgraph(c)
                step_name = next((s_g2_n for s_g2_n, step_id in s_g2.nodes(data="id") if
                                  s_g2_n.startswith("step::")), None)
                if step_name:
                    sub_graph_by_step[step_name] = s_g2

        sub_graphs_by_module[module]["steps"] = sub_graph_by_step

    return sub_graphs_by_module


def write_all_mermaid_flowcharts(flowchart, full_package_name, sub_graphs_by_module=None):
    if sub_graphs_by_module is None:
        sub_graphs_by_module = divide_graph_by_module(flowchart)

    package = importlib.import_module(full_package_name)
    modules = pkgutil.iter_modules(package.__path__, prefix=full_package_name + '.')
    for module in modules:
        if module.name.split(".")[-1] in MODULES_TO_SKIP:
            continue

        module_f = importlib.import_module(module.name)
        if module.ispkg:
            path = package_flowchart_path(module_f, local=True)
            write_all_mermaid_flowcharts(
                flowchart, module.name, sub_graphs_by_module=sub_graphs_by_module
            )
        else:
            path = module_flowchart_path(module_f, local=True)
            module_sub_graphs = sub_graphs_by_module[module.name]

            sub_fc = module_sub_graphs["module_subgraph"]
            write_mermaid_flowchart(sub_fc, path, step_sub_graphs=module_sub_graphs["steps"])


def find_resource_node_creator_module(flowchart, resource_node):
    """
    Search recursively through the predecessors of a resource node to find the module
    where the resource node is created, specifically looking for a 'step' type node.

    :param flowchart: A graph object representing the flowchart.
    :param resource_node: The resource node for which the creator module is sought.
    :return: The name of the module that creates the resource node, or None if not found.
    """
    # Iterate through each predecessor of the resource node
    for predecessor in flowchart.predecessors(resource_node):
        # Fetch node data
        node_data = flowchart.nodes[predecessor]

        # Check if the current node is a 'step' type
        if node_data.get("type") == "step":
            return node_data.get("module")

        # Recursive call to check the predecessors of the current node
        creator_module = find_resource_node_creator_module(flowchart, predecessor)
        if creator_module is not None:
            return creator_module

    # Return None if no creator module found
    return None


def find_input_modules(flowchart, resource_node, prev_nodes=None):
    """
    Recursively collect the modules that serve as inputs to a given resource node within the flowchart.
    This function identifies "step" type nodes that directly connect to the resource node and recursively
    collects modules for any predecessors not yet visited.

    :param flowchart: A graph object representing the flowchart.
    :param resource_node: The node for which input modules are sought.
    :param prev_nodes: A set of nodes that have been visited to avoid cycles.
    :return: A list of modules that serve as inputs to the specified resource node.
    """
    if prev_nodes is None:
        prev_nodes = set()

    input_steps = []
    # Iterate over all neighbors of the resource node
    for neighbor in flowchart.neighbors(resource_node):
        node_data = flowchart.nodes[neighbor]

        # Collect modules from "step" type nodes directly
        if node_data.get("type") == "step":
            input_steps.append(node_data.get("module"))
        # Recurse into neighbors not previously visited to avoid cycles
        elif neighbor not in prev_nodes:
            prev_nodes.add(resource_node)  # Mark the current node as visited
            input_steps.extend(find_input_modules(flowchart, neighbor, prev_nodes))

    return input_steps


def main(args):
    hl.init()
    sys.path.insert(0, REPOSITORY_ROOT_PATH)

    if args.create_full_qc_graph:
        fc = QCFlowchart()
        for full_module_name in iter_top_level_package_modules():
            module = importlib.import_module(full_module_name)
            parse_qc_package(fc, module, full_module_name)

        with open('QC_flowchart.gpickle', 'wb') as f:
            pickle.dump(fc.graph, f, pickle.HIGHEST_PROTOCOL)

    if args.write_mermaid_flowchart:
        with open('QC_flowchart.gpickle', 'rb') as f:
            fc = pickle.load(f)

        node_attribute_change = {}
        for n, d in fc.nodes(data=True):
            if d.get("type") not in {"resource"}:
                continue
            if d.get("is_raw_resource", ""):
                node_attribute_change[n] = {'type': "raw_resource"}

            creator_step = find_resource_node_creator_module(fc, n)
            input_steps = find_input_modules(fc, n)

            if creator_step is not None and len(input_steps) > 0 and {creator_step} != set(input_steps):
                node_attribute_change[n] = {'type': "main_resource"}

        nx.set_node_attributes(fc, node_attribute_change)

        for full_module_name in iter_top_level_package_modules():
            if full_module_name not in MODULES_TO_SKIP:
                write_all_mermaid_flowcharts(fc, full_module_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--create-full-qc-graph", help="Create the full QC graph", action="store_true"
    )
    parser.add_argument(
        "--write-mermaid-flowchart",
        help="Write the mermaid flowchart to a file",
        action="store_true",
    )
    args = parser.parse_args()
    main(args)
