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
from collections import defaultdict

import hail as hl

from gnomad.resources.resource_utils import (
    BaseResource,
    DataException,
    BaseVersionedResource,
)
from gnomad.utils.file_utils import check_file_exists_raise_error

from generate_api_reference import module_flowchart_path, package_flowchart_path

# Need to add the repository root to the path before importing any gnomad_qc modules.
REPOSITORY_ROOT_PATH = str(pathlib.Path(os.path.abspath(__file__)).parent.parent)
sys.path.insert(0, REPOSITORY_ROOT_PATH)

from gnomad_qc.v4.resources.variant_qc import TRUTH_SAMPLES
from gnomad_qc.v4.resources.basics import gnomad_v4_genotypes

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
    Class to represent a directed graph (flowchart) of a quality control pipeline.
    """
    def __init__(self, graph=None):
        """
        Initialize a QCFlowchart object.
        """
        if graph is not None:
            self.graph = nx.DiGraph(graph)
        else:
            self.graph = nx.DiGraph()

    def add_node(
        self,
        node,
        node_type,
        node_id,
        label,
        module_name,
        href=None,
        script_name=None,
        **kwargs,
    ):
        """
        Add a node to the graph.
        """
        full_node_name = node
        if script_name is not None:
            full_node_name = f"{script_name}::{full_node_name}"

        node_attr = {
            "name": node,
            "type": node_type,
            "id": node_id,
            "label": label,
            "module": module_name,
            "href": href if href is not None else "",
            "script_name": script_name,
            **kwargs,
        }
        self.graph.add_node(f"{node_type}::{full_node_name}", **node_attr)

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
        Add an edge to the graph.
        """
        from_node_new = from_node
        if from_script_name is not None:
            from_node_new = f"{from_script_name}::{from_node_new}"

        if from_node_type is not None:
            if not self.has_node(from_node, from_node_type, from_script_name):
                self.add_node(
                    from_node,
                    node_type=from_node_type,
                    node_id=f"{from_node_type}_{from_node}",
                    label=from_node,
                    module_name="",
                    script_name=from_script_name,
                )
            from_node_new = f"{from_node_type}::{from_node_new}"

        to_node_new = to_node
        if to_script_name is not None:
            to_node_new = f"{to_script_name}::{to_node_new}"

        if to_node_type is not None:
            if not self.has_node(to_node, to_node_type, to_script_name):
                self.add_node(
                    to_node,
                    node_type=to_node_type,
                    node_id=f"{to_node_type}_{to_node}",
                    label=to_node,
                    module_name="",
                    script_name=to_script_name,
                )
            to_node_new = f"{to_node_type}::{to_node_new}"

        self.graph.add_edge(from_node_new, to_node_new, **kwargs)

    def has_node(self, node, node_type=None, script_name=None):
        """
        Check if the graph has a node.
        """
        if script_name is not None:
            node = f"{script_name}::{node}"

        if node_type is not None:
            node = f"{node_type}::{node}"

        return self.graph.has_node(node)

    def has_edge(
        self,
        from_node,
        to_node,
        from_node_type=None,
        to_node_type=None,
        from_script_name=None,
        to_script_name=None,
    ):
        """
        Check if the graph has an edge.
        """
        if from_script_name is not None:
            from_node = f"{from_script_name}::{from_node}"
        if from_node_type is not None:
            from_node = f"{from_node_type}::{from_node}"

        if to_script_name is not None:
            to_node = f"{to_script_name}::{to_node}"
        if to_node_type is not None:
            to_node = f"{to_node_type}::{to_node}"

        return self.graph.has_edge(from_node, to_node)

    def remove_edge(
        self,
        from_node,
        to_node,
        from_node_type=None,
        to_node_type=None,
        from_script_name=None,
        to_script_name=None,
    ):
        """
        Remove an edge from the graph.
        """
        if from_script_name is not None:
            from_node = f"{from_script_name}::{from_node}"
        if from_node_type is not None:
            from_node = f"{from_node_type}::{from_node}"

        if to_script_name is not None:
            to_node = f"{to_script_name}::{to_node}"
        if to_node_type is not None:
            to_node = f"{to_node_type}::{to_node}"

        self.graph.remove_edge(from_node, to_node)

    def remove_node(self, node, node_type=None, script_name=None):
        """
        Remove a node from the graph, and resolve edges.
        """
        if script_name is not None:
            node = f"{script_name}::{node}"
        if node_type is not None:
            node = f"{node_type}::{node}"

        predecessors = list(self.graph.predecessors(node))
        children = list(self.neighbors(node))

        self.graph.remove_node(node)
        for p in predecessors:
            for c in children:
                self.add_edge(p, c)

    def neighbors(self, node, node_type=None, script_name=None):
        """
        Get the neighbors of a node.
        """
        if script_name is not None:
            node = f"{script_name}::{node}"
        if node_type is not None:
            node = f"{node_type}::{node}"

        return list(self.graph.neighbors(node))

    def get_graph(self):
        """
        Get the graph.
        """
        return self.graph

    def insert_nodes(self, parent, parent_type, parent_script, nodes):
        """
        Insert nodes into the graph with edges from parent.
        """
        if not self.has_node(parent, parent_type, parent_script):
            children = []
        else:
            children = list(self.neighbors(parent, parent_type, parent_script))

        last_node = parent
        last_node_type = parent_type
        last_node_script = parent_script
        for node, kwargs in nodes:
            self.add_node(node, **kwargs)
            node_type = kwargs["node_type"]
            script_name = kwargs["script_name"]
            self.add_edge(last_node, node, last_node_type, node_type, last_node_script, script_name)
            last_node = node
            last_node_type = node_type
            last_node_script = script_name

        for node in children:
            self.add_edge(
                last_node,
                self.graph.nodes[node]["name"],
                last_node_type,
                self.graph.nodes[node]["type"],
                last_node_script,
                self.graph.nodes[node]["script_name"],
            )
            self.remove_edge(
                parent, node, from_node_type=parent_type, from_script_name=parent_script
            )

    def set_node_attributes(self, attrs, values=None):
        """
        Set attributes for a node.
        """
        nx.set_node_attributes(self.graph, attrs, values)


class BaseQCParser:
    def __init__(self, name, obj, module):
        self.name = name
        self.obj = obj
        self.module = module
        self.is_resource = self.get_is_resource()
        self.is_module_func = self.get_is_module_func()
        self.type = None

    def get_is_resource(self):
        """
        Check if the resource object is a resource.
        """
        return (
            issubclass(type(self.obj), BaseResource)
            or issubclass(type(self.obj), BaseVersionedResource)
        )

    def get_is_module_func(self):
        """
        Check if the resource object is a function defined in the current module.
        """
        return (
            inspect.isfunction(self.obj)
            and self.obj.__module__ == self.module.name
        )

    def get_signature(self):
        """
        Get the signature of the resource function.
        """
        return inspect.signature(self.obj)

    def get_resource_kwargs(self, kwargs_options):
        """
        Get Dictionary of kwargs for the PipelineStepNode function.

        :return: Dictionary of kwargs for the PipelineStepNode function.
        """
        for k, v in self.get_signature().parameters.items():
            if (k == "test") or (k == "overwrite"):
                kwargs_options[k] = [False]
            elif k in kwargs_options:
                continue
            elif v.annotation.__name__ == "bool":
                kwargs_options[k] = [True, False]
            elif k == "data_type":
                kwargs_options[k] = ["exomes", "genomes", "joint"]

        return kwargs_options

    def get_resource_kwargs_combos(self, kwargs_options):
        """
        Get all combinations of kwargs for the ResourceNode function.

        :return: List of all combinations of kwargs for the ResourceNode function.
        """
        all_options = []
        for k, v in self.get_resource_kwargs(kwargs_options).items():
            all_options.append([(k, val) for val in v])

        return list(itertools.product(*all_options))

    @staticmethod
    def remove_excluded_parameters(bind_kwargs):
        """
        Remove excluded parameters from the resource function kwargs.

        :return: Dictionary of resource function kwargs with excluded parameters removed.
        """
        if bind_kwargs is None:
            bind_kwargs = {}
        return {k: v for k, v in bind_kwargs}

    def get_resource_display_name(self, bind_kwargs=None):
        """
        Get the display name for the ResourceNode.

        :param bind_kwargs: Optional resource function kwargs to run the function with.
        :return: Display name for the ResourceNode.
        """
        if bind_kwargs is None:
            res_display = ""
        else:
            bind_kwargs = self.remove_excluded_parameters(bind_kwargs)
            res_display = ",\n".join("{}={}".format(k, v) for k, v in bind_kwargs.items())

            if len(bind_kwargs) > 1:
                res_display = "\n" + res_display

        return f"{self.name}({res_display})"

    def get_resource_id(self, bind_kwargs=None):
        """
        Get the ID for the ResourceNode.

        :param bind_kwargs: Optional resource function kwargs to run the function with.
        :return: ID for the ResourceNode.
        """
        bind_kwargs = self.remove_excluded_parameters(bind_kwargs)
        name_reformat = self.name.replace('--', '').replace('-', '_')
        if bind_kwargs is None or len(bind_kwargs) == 0:
            return f"{self.type}_{name_reformat}"
        else:
            res_arg_str = "_".join("{}_{}".format(k, v) for k, v in bind_kwargs.items())

            return f"{self.type}_{name_reformat}_{res_arg_str}"

    def get_doc_link(self, relative=True, add_function=True):
        """
        Get the full link to the resource's documentation.
        """
        if add_function:
            return get_function_doc_link(self.module.name, self.name, relative=relative)
        else:
            return get_module_doc_link(self.module.name, relative=relative)

    @staticmethod
    def get_resource_paths(resource):
        """
        Run a resource function and add its resource nodes to the flowchart.

        """
        if hasattr(resource, "__name__") and resource.__name__ == RAW_RESOURCE_FUNC:
            resource = gnomad_v4_genotypes
        if hasattr(resource, "versions"):
            paths = [r.path for r in resource.versions.values()]
        elif hasattr(resource, "path"):
            paths = [resource.path]
        elif isinstance(resource, str):
            paths = [resource]
        else:
            return None

        paths = [
            p for p in paths if check_file_exists_raise_error(p, error_if_exists=False)
        ]

        return paths

    def add_resource_to_flowchart(self, flowchart, resource, bind_kwargs=None):
        """
        Run a resource function and add its resource nodes to the flowchart.

        :param flowchart: QCFlowchart object.
        """
        paths = self.get_resource_paths(resource)
        if paths is None:
            err_msg = "" if bind_kwargs is None else f"with {bind_kwargs}"
            raise ValueError(
                f"Resource {self.obj} {err_msg}did not return a valid resource."
            )

        is_raw_resource = True if self.name in EXTERNAL_RESOURCES else False
        for p in paths:
            flowchart.add_node(
                p,
                node_type=self.type,
                node_id=self.get_resource_id(bind_kwargs),
                label=self.get_resource_display_name(bind_kwargs),
                module_name=self.module.name,
                href=self.get_doc_link(),
                is_raw_resource=is_raw_resource,
                #parser=self,
            )

    def iter_func_call_combos(self, kwargs_options=None):
        if kwargs_options is None:
            kwargs_options = self.get_resource_kwargs_combos(kwargs_options)
        for bind_kwargs in kwargs_options:
            res_bind = self.get_signature().bind(**dict(bind_kwargs))
            try:
                yield bind_kwargs, self.obj(*res_bind.args, **res_bind.kwargs)
            except (ValueError, DataException, KeyError):
                continue

    def add_to_flowchart(self, flowchart):
        """
        Run a resource function and add its resource nodes to the flowchart.

        :param flowchart: QCFlowchart object.
        """
        for bind_kwargs, func_call in self.iter_func_call_combos():
            self.add_resource_to_flowchart(
                flowchart, func_call, bind_kwargs=bind_kwargs
            )


########################################################################################
# Resource and PipelineStepNode parsers to help construct the flowchart.
########################################################################################
class ResourceParser(BaseQCParser):

    def __init__(self, res_name, res_obj, module):
        """
        Initialize a ResourceNode object.

        :param res_name: Function name.
        :param res_obj: Resource function object.
        """
        super().__init__(res_name, res_obj, module)
        self.type = "resource"

    @staticmethod
    def remove_excluded_parameters(bind_kwargs):
        """
        Remove excluded parameters from the resource function kwargs.

        :return: Dictionary of resource function kwargs with excluded parameters removed.
        """
        if bind_kwargs is None:
            return {}
        return {
            k: v for k, v in bind_kwargs if k not in EXCLUDE_RESOURCE_PARAMS
        }

    def get_resource_kwargs_combos(self, kwargs_options=None):
        """
        Get all combinations of kwargs for the ResourceNode function.

        :return: List of all combinations of kwargs for the ResourceNode function.
        """
        if kwargs_options is None:
            kwargs_options = KWARG_OPTIONS_BY_RESOURCE.get(self.name, {})

        return super().get_resource_kwargs_combos(kwargs_options)

    def add_to_flowchart(self, flowchart):
        """
        Run a resource function and add its resource nodes to the flowchart.

        :param flowchart: QCFlowchart object.
        """
        if self.is_module_func and self.name != RAW_RESOURCE_FUNC:
            super().add_to_flowchart(flowchart)
        else:
            super().add_resource_to_flowchart(flowchart, self.obj, bind_kwargs=None)


class PipelineResourceParser(BaseQCParser):

    def __init__(self, name, obj, module):
        """
        Initialize a PipelineStepNode object.

        :param name: Pipeline step name.
        :param obj: Pipeline step object.
        :param module: Module object.
        """
        super().__init__(name, obj, module)
        self.script_name = module.name.split(".")[-1] + ".py"
        self.type = "step"

    def get_resource_kwargs(self, kwargs_options=None):
        """
        Get Dictionary of kwargs for the PipelineStepNode function.

        :return: Dictionary of kwargs for the PipelineStepNode function.
        """
        if kwargs_options is None:
            kwargs_options = KWARG_OPTIONS_BY_MODULE.get(self.module.name, {})

        return super().get_resource_kwargs(kwargs_options)

    def add_step_to_flowchart(self, flowchart, step_obj):
        step_name = step_obj.pipeline_step
        step_id = step_name.replace("--", "").replace("-", "_")

        output_paths = []
        for n, r in step_obj.output_resources.items():
            for r2 in r:
                output_paths.extend(self.get_resource_paths(r2) or [])

        edges = []
        step_created = {}
        steps_using = defaultdict(set)
        if len(output_paths) > 0:
            for n, r in step_obj.input_resources.items():
                for r2 in r:
                    paths = self.get_resource_paths(r2) or []
                    for p in paths:
                        if not flowchart.has_edge(p, step_id, "resource", "step", to_script_name=self.script_name):
                            steps_using[f"resource::{p}"].add(f"step::{self.script_name}::{step_id}")
                            edges.append((p, step_id, "resource", "step", None, self.script_name))
            for p in output_paths:
                step_created[f"resource::{p}"] = f"step::{self.script_name}::{step_id}"
                if not flowchart.has_edge(step_id, p, "step", "resource", from_script_name=self.script_name):
                    edges.append((step_id, p, "step", "resource", self.script_name, None))

        if len(edges) > 0:
            flowchart.add_node(
                step_id,
                node_type=self.type,
                node_id=f"step_{step_id}",
                label=step_name,
                module_name=self.module.name,
                script_name=self.script_name,
                #step_node=step_obj,
            )

            for e in edges:
                flowchart.add_edge(*e)

            if step_created:
                flowchart.set_node_attributes(step_created, values="step_created")

            if steps_using:
                flowchart.set_node_attributes(steps_using, values="steps_using")

    def add_to_flowchart(self, flowchart):
        for bind_kwargs, func_call in self.iter_func_call_combos():
            for step_obj in func_call.pipeline_steps.values():
                self.add_step_to_flowchart(flowchart, step_obj)


########################################################################################
# ScriptMainParser and CallVisitor classes to parse the main function of a script.
########################################################################################
class CallVisitor(ast.NodeVisitor):
    def __init__(self):
        self._name = []
        self._const_value = None

    @property
    def name(self):
        return '.'.join(reversed(self._name)).replace("args.", "")

    @name.deleter
    def name(self):
        self._name = []

    @property
    def const_value(self):
        return self._const_value

    @const_value.deleter
    def const_value(self):
        self._const_value = None

    def visit_Name(self, node):
        self._name.append(node.id)

    def visit_Attribute(self, node):
        try:
            self._name.append(node.attr)
            self._name.append(node.value.id)
        except AttributeError:
            self.generic_visit(node)

    def visit_Constant(self, node):
        self._const_value = node.value

    def visit_keyword(self, node):
        return self.generic_visit(node)


class ScriptMainParser(BaseQCParser):

    def __init__(self, name, obj, module):
        super().__init__(name, obj, module)
        self.rel_module_name = module.name.split(".")[-1]
        self.script_name = self.rel_module_name + ".py"
        self.steps = []
        self.step_func_calls = defaultdict(list)
        self.step_raw_input = {}
        self.type = "func"
        self.function_to_module = self.get_module_function_call_dict()

    def get_module_function_call_dict(self):
        module_f = importlib.import_module(self.module.name)
        function_to_module = {}
        for name, obj in inspect.getmembers(module_f, inspect.isfunction):
            function_to_module[name] = (obj.__module__, obj.__name__)

        return function_to_module

    def parse_ast_node(self, node, flowchart, temp=""):
        temp += "\t"
        for child in ast.iter_child_nodes(node):
            child.parent = node
            child.step = node.step or None
            child.caller = node.caller or []

            if isinstance(child, ast.If):
                callvisitor = CallVisitor()
                callvisitor.visit(child.test)
                if flowchart.has_node(callvisitor.name, "step", self.script_name):
                    child.step = callvisitor.name
                    self.steps.append(callvisitor.name)
            elif isinstance(child, ast.Call):
                callvisitor = CallVisitor()
                callvisitor.visit(child.func)
                func_module, func_name = self.function_to_module.get(
                    callvisitor.name, ("", "")
                )
                pkgs = func_module.split(".")
                in_resource_pkg = len(pkgs) > 2 and pkgs[2] == "resources"
                if in_resource_pkg and func_module and func_name == RAW_RESOURCE_FUNC:
                    kwargs = {}
                    for k in child.keywords:
                        if k.arg in EXCLUDE_RESOURCE_PARAMS:
                            continue
                        callvisitor.visit(k)
                        kwargs[k.arg] = callvisitor.const_value

                    self.step_raw_input[self.script_name] = kwargs
                    if child.step is not None:
                        self.step_raw_input[child.step] = kwargs

                if func_module and func_name not in EXCLUDE_FUNCTIONS and not in_resource_pkg:
                    script_calls = self.step_func_calls[self.script_name]
                    if not script_calls or script_calls[-1] != (func_name, func_module):
                        self.step_func_calls[self.script_name].append(
                            (func_name, func_module)
                        )
                    step_calls = self.step_func_calls.get(child.step, [])
                    if child.step is not None and (not step_calls or step_calls[-1] != (func_name, func_module)):
                        self.step_func_calls[child.step].append(
                            (func_name, func_module)
                        )

                    pkgs = func_module.split(".")
                    # Prevents infinite recursion if a function calls itself.
                    if pkgs[0] == "gnomad_qc" and func_name not in child.caller:
                        func_parse = ast.parse(
                            inspect.getsource(
                                getattr(importlib.import_module(func_module), func_name)
                            )
                        )

                        func_parse.step = child.step
                        func_parse.parent = child
                        func_parse.caller = child.caller + [func_name]
                        self.parse_ast_node(func_parse, flowchart, temp)
                        continue

            self.parse_ast_node(child, flowchart, temp)

    def parse_main_ast(self, flowchart):
        # Get all AST nodes for each step in the script that we can use to extract
        # all function calls within that step.
        main_ast_node = ast.parse(inspect.getsource(self.obj))
        main_ast_node.parent = None
        main_ast_node.step = None
        main_ast_node.caller = []
        self.parse_ast_node(main_ast_node, flowchart)

    def add_to_flowchart(self, flowchart):
        """
        Run a resource function and add its resource nodes to the flowchart.

        :param flowchart: QCFlowchart object.
        """
        self.parse_main_ast(flowchart)
        steps = self.steps or [self.script_name]
        for step in steps:
            step_func_nodes = []
            for step_func in self.step_func_calls[step]:
                pkgs = step_func[1].split(".")
                if pkgs[0] not in {"gnomad", "gnomad_qc"}:
                    if pkgs[0] == "hail":
                        print("Using hail function:", step_func)
                    continue
                if pkgs[0] == "gnomad_qc" and len(pkgs) > 2 and pkgs[2] == "resources":
                    continue
                if pkgs[0] == "gnomad_qc":
                    p = f"{step_func[0]}::{step}"
                    node_type = "func"
                    node_id = f"{node_type}_{step_func[0]}_{step}"
                    href = get_function_doc_link(
                        step_func[1], step_func[0], relative=True
                    )
                else:
                    p = f"{step_func[0]}::{step}"
                    node_type = "gnomad_methods"
                    node_id = f"{node_type}_{step_func[0]}_{step}"
                    href = get_function_doc_link(
                        step_func[1], step_func[0], gnomad_methods=True
                    )

                step_func_nodes.append(
                    (
                        p,
                        {
                            "node_type": node_type,
                            "node_id": node_id,
                            "label": f"{step_func[0]}()",
                            "href": href,
                            "module_name": self.module.name,
                            "script_name": self.script_name,
                            "step_name": step,
                            #"script_main": self,
                        },
                    )
                )

            if step_func_nodes:
                flowchart.insert_nodes(step, "step", self.script_name, step_func_nodes)

            if step in self.step_raw_input:
                bind_kwargs = self.step_raw_input[step]
                p = gnomad_v4_genotypes.path
                label = RAW_RESOURCE_FUNC
                node_id = f"resource_{label}"

                if bind_kwargs:
                    label_kwargs = ",\n".join("{}={}".format(k, v) for k, v in bind_kwargs.items())
                    label = f"{label}({label_kwargs})"
                    node_id = f"{node_id}_{'_'.join(f'{k}_{v}' for k, v in bind_kwargs.items())}"
                    p = f"{p}_{'_'.join(f'{k}_{v}' for k, v in bind_kwargs.items())}"

                m_name = "gnomad_qc.v4.resources.basics"

                if not flowchart.has_node(p, "resource"):
                    flowchart.add_node(
                        p,
                        node_type="resource",
                        node_id=node_id,
                        label=label,
                        module_name=m_name,
                        href=get_function_doc_link(m_name, RAW_RESOURCE_FUNC, relative=True),
                        is_raw_resource=True,
                    )

                flowchart.add_edge(
                    p,
                    step,
                    "resource",
                    "step",
                    None,
                    self.script_name,
                )


########################################################################################
# Main function to parse the repo and construct the flowchart.
########################################################################################
def parse_qc_package(flowchart, package, full_package_name, resource_pkg=False):
    """
    Parse a resource package and add its resource nodes to the QCFlowchart.

    :param flowchart: QCFlowchart object.
    :param package: Package object for resources.
    :param full_package_name: Full package name.
    :return: None.
    """
    if (full_package_name in RESOURCE_PACKAGES) or resource_pkg:
        resource_pkg = True

    modules = pkgutil.iter_modules(package.__path__, prefix=full_package_name + '.')
    for module in modules:
        rel_module_name = module.name.split(".")[-1]
        module_f = importlib.import_module(module.name)
        if module.ispkg:
            parse_qc_package(flowchart, module_f, module.name, resource_pkg=resource_pkg)
            continue

        # Skip modules that are not resources and are in the MODULES_TO_SKIP list for
        # the flowchart scripts.
        if not resource_pkg and rel_module_name in MODULES_TO_SKIP:
            continue

        resource_parsers = []
        pipeline_parsers = []
        main_parsers = []

        for name, obj in inspect.getmembers(module_f):
            # If the member is in a resource package, skip private functions and
            # excluded resources.
            if resource_pkg and not(name.startswith("_") or name in EXCLUDE_RESOURCES):
                res_parser = ResourceParser(name, obj, module)

                # Skip members that are not resources or functions defined in the
                # current module.
                if not res_parser.is_resource and not res_parser.is_module_func:
                    continue

                resource_parsers.append(res_parser)

            elif name == "get_pipeline_resources":
                pipeline_parser = PipelineResourceParser(name, obj, module)
                # Skip members that are not functions defined in the current module.
                if not pipeline_parser.is_module_func:
                    continue

                pipeline_parsers.append(pipeline_parser)

            elif name == "main":
                main_parsers.append(ScriptMainParser(name, obj, module))

        for parser in resource_parsers + pipeline_parsers + main_parsers:
            parser.add_to_flowchart(flowchart)


def iter_top_level_package_modules():
    """
    Create the flowchart for the gnomAD QC pipeline.
    """
    packages = setuptools.find_namespace_packages(
        where=REPOSITORY_ROOT_PATH, include=["gnomad_qc.*"]
    )
    top_level_packages = [pkg for pkg in packages if pkg in INCLUDE_TOP_LEVEL_PACKAGES]

    for package_name in top_level_packages:
        yield package_name
        #package = importlib.import_module(package_name)
        #for module in pkgutil.iter_modules(package.__path__):
            #if not module.ispkg:
            #    continue
        #    full_module_name = f"{package_name}.{module.name}"
        #    yield full_module_name


########################################################################################
# Functions to format the flowchart nodes and connections for mermaid.js.
########################################################################################
def format_node(node_name, node_str, node_type, href=None):
    node_name = node_name.strip("_")
    if href:
        if href.startswith("http"):
            node_str = r"""<a href='{href}'>{node_str}</a>""".format(href=href, node_str=node_str)
        else:
            node_str = r"""<a class="reference internal" href='{href}'>{node_str}</a>""".format(href=href, node_str=node_str)
    if node_type == "step":
        node = f"""  {node_name}{{{{"{node_str}"}}}}:::{node_type}_color;"""
    elif node_type in {"func", "gnomad_methods"}:
        node = f"""  {node_name}[["{node_str}"]]:::{node_type}_color;"""
    elif node_type in {"resource", "main_resource", "raw_resource"}:
        node = f"""  {node_name}[/"{node_str}"/]:::{node_type}_color;"""
    else:
        raise ValueError(f"Node type {node_type} not recognized.")

    return node


def format_connection(node1, node2, invisible=False):
    node1 = node1.strip("_")
    node2 = node2.strip("_")
    if invisible:
        connector = "~~~"
    else:
        connector = "-->"
    return f"""  {node1} {connector} {node2};"""


def networkx_sg_to_mermaid(graph, subgraph, subgraph_name, href_relative_to=None):
    mermaid_by_group = defaultdict(list)

    resource_types = {"raw_resource", "main_resource", "resource"}
    for node, data in subgraph.nodes(data=True):
        if "id" not in data or subgraph.degree(node) == 0:
            continue
        else:
            href = data.get("href")
            if href and not href.startswith("http") and href_relative_to:
                href_relative_to = href_relative_to.split("v4/")[-1]
                href = os.path.relpath(href, href_relative_to)

            node_format = format_node(data["id"], data["label"], data["type"], href)
            if data.get("type") == "raw_resource":
                mermaid_by_group["raw_input_resources"].append(node_format)
            elif data.get("type") == "main_resource":
                if len(list(graph.neighbors(node))) == 0:
                    mermaid_by_group["output_resources"].append(node_format)
                #elif len(list(graph.predecessors(node))) == 0:
                #    print("Adding as input resource in networkx_sg_to_mermaid:", node)
                #    mermaid_by_group["input_resources"].append(node_format)
            elif data.get("type") not in {"step", "resource"}:
                mermaid_by_group["nodes"].append(node_format)

    for from_node, to_node in subgraph.edges(data=False):
        from_node_data = subgraph.nodes[from_node]
        to_node_data = subgraph.nodes[to_node]
        if "id" not in from_node_data or "id" not in to_node_data:
            continue

        if from_node_data.get("type") in resource_types:
            con_format = format_connection(from_node_data["id"], subgraph_name)
            group_type = "resource_edges"
        elif to_node_data.get("type") in resource_types:
            continue
        #    con_format = format_connection(subgraph_name, to_node_data["id"])
        #    group_type = "resource_edges"
        elif from_node_data.get("id") == subgraph_name:
            con_format = format_connection("n_"+from_node_data["id"], to_node_data["id"], invisible=True)
            group_type = "edges"
        else:
            con_format = format_connection(from_node_data["id"], to_node_data["id"])
            group_type = "edges"

        mermaid_by_group[group_type].append(con_format)

    return mermaid_by_group


def networkx_to_mermaid(graph, href_relative_to=None, add_module_subgraphs=False, step_sub_graphs=None):
    mermaid_by_group = defaultdict(list)
    for node, data in graph.nodes(data=True):
        if "id" not in data:
            continue
            #print("Missing resource node:", node, data)
        elif graph.degree(node) == 0:
            continue
            #print("Node with no edges:", node, data)
            #print("Main resource node with no edges:", node, data)
        else:
            href = data.get("href")
            #print("\t", data)
            if href and not href.startswith("http") and href_relative_to:
                href_relative_to = href_relative_to.split("v4/")[-1]
                href = os.path.relpath(href, href_relative_to)

            node_format = format_node(data["id"], data["label"], data["type"], href)
            if data.get("type") in {"raw_resource"}:
                mermaid_by_group["raw_input_resources"].append(node_format)
                continue

            if data.get("type") in {"main_resource"}:
                if len(list(graph.predecessors(node))) == 0:
                    mermaid_by_group["input_resources"].append(node_format)
                    continue
                #if len(list(graph.neighbors(node))) == 0:
                mermaid_by_group["output_resources"].append(node_format)
                continue

            if data.get("type") == "step" and step_sub_graphs is not None and data.get("id") in step_sub_graphs:
                sg = step_sub_graphs[data.get("id")]
                step_by_group = networkx_sg_to_mermaid(graph, sg, data.get("id"), href_relative_to=href_relative_to)
                mermaid_by_group[data.get("id")].append(
                    format_node("n_" + data["id"], data["label"], data["type"], href)
                )
                mermaid_by_group[data.get("id")].extend(step_by_group["nodes"] + step_by_group["edges"])
                mermaid_by_group["edges"].extend(step_by_group["resource_edges"])
                for resource_type in {"output_resources", "input_resources", "raw_input_resources"}:
                    mermaid_by_group[resource_type].extend(step_by_group[resource_type])
                continue

            mermaid_by_group["nodes"].append(node_format)

    for from_node, to_node in graph.edges(data=False):
        from_node_data = graph.nodes[from_node]
        to_node_data = graph.nodes[to_node]
        if "id" not in from_node_data or "id" not in to_node_data:
            continue

        if from_node_data.get("type") == "step" and step_sub_graphs is not None and from_node_data.get("id") in step_sub_graphs:
            continue

        con_format = format_connection(from_node_data["id"], to_node_data["id"])
        if from_node_data.get("type") == "step" or to_node_data.get("type") == "step":
            module = from_node_data.get("module")
            if from_node_data.get("type") == "step" and module and add_module_subgraphs:
                mermaid_by_group[module.split(".")[-1]].append(con_format)
                continue

            module = to_node_data.get("module")
            if to_node_data.get("type") == "step" and module and add_module_subgraphs:
                mermaid_by_group[module.split(".")[-1]].append(con_format)
                continue

        mermaid_by_group["edges"].append(
            format_connection(from_node_data["id"], to_node_data["id"])
        )

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


def write_mermaid_flowchart(flowchart, path, step_sub_graphs=None):  # add_module_subgraphs=False):
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
            #add_module_subgraphs,
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
            #f.write(f"  subgraph {module} [Module: {module}]\n")
            for n in mod_nodes:
                f.write(f"  {n}\n")
            f.write("  end\n")

        for l in edges:
            f.write(l + "\n")


def divide_graph_by_module(graph):
    from collections import Counter
    resource_types = {"raw_resource", "main_resource", "resource"}

    g_no_input_edges = graph.copy()
    resource_creation_step = dict(
        g_no_input_edges.nodes(data="step_created", default="")
    )
    node_module = dict(g_no_input_edges.nodes(data="module", default=""))
    g_no_input_edges.remove_edges_from(
        [
            e for e in g_no_input_edges.edges
            if (
                e[0].startswith("resource::")
                and e[1].startswith("step::")
                and (
                    node_module.get(resource_creation_step[e[0]], "")
                    != node_module.get(e[1], "")
                )
            )
        ]
    )

    by_module = defaultdict(list)
    all_sub_graphs = list(nx.connected_components(g_no_input_edges.to_undirected()))
    for s in all_sub_graphs:
        s_g = g_no_input_edges.subgraph(s)
        n_types = Counter([t for k, t in s_g.nodes(data="type", default="")])
        if len(s) == 1 and n_types.get("resource", 0) == 1:
            continue
        non_resource_modules = [
            m for k, m in s_g.nodes(data="module", default="")
            if m != "" and ((len(m.split(".")) < 3) or (m.split(".")[2] != "resources"))
        ]
        if len(non_resource_modules) != 0:
            by_module[set(non_resource_modules).pop()].append(s_g)

    sub_graphs_by_module = {}
    for module, s_g in by_module.items():
        module_nodes = []
        for g in s_g:
            for n in g.nodes():
                module_nodes.append(n)
                for n2 in list(graph.predecessors(n)) + list(graph.neighbors(n)):
                    if graph.nodes[n2]["type"] in resource_types:
                        module_nodes.append(n2)

        module_sg = graph.subgraph(module_nodes).copy()
        sub_graphs_by_module[module] = {"module_subgraph": module_sg}
        module_sg = module_sg.copy()
        module_sg.remove_edges_from(
            [
                e for e in module_sg.edges
                if (
                    e[0].startswith("resource::")
                    or e[1].startswith("resource::")
                    or (
                        not (e[0].startswith("func::") or e[0].startswith("gnomad_methods::"))
                        and not (e[1].startswith("func::") or e[1].startswith("gnomad_methods::"))
                    )
                )
            ]
        )
        sub_graph_by_step = {}
        for c in list(nx.connected_components(module_sg.to_undirected())):
            if len(c) > 1:
                s_g2 = module_sg.subgraph(c)
                step_name = [
                    step_id for s_g2_n, step_id in s_g2.nodes(data="id", default="")
                    if s_g2_n.startswith("step::")
                ][0]
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
            #if module.name != "gnomad_qc.v4.sample_qc.assign_ancestry":
            #    continue
            write_mermaid_flowchart(sub_fc, path, step_sub_graphs=module_sub_graphs["steps"])  #, add_module_subgraphs=module.ispkg)


def write_all_mermaid_flowcharts_old(flowchart, full_package_name):
    package = importlib.import_module(full_package_name)
    modules = pkgutil.iter_modules(package.__path__, prefix=full_package_name + '.')
    for module in modules:
        if module.name.split(".")[-1] in MODULES_TO_SKIP:
            continue

        module_f = importlib.import_module(module.name)
        if module.ispkg:
            path = package_flowchart_path(module_f, local=True)
            write_all_mermaid_flowcharts(flowchart, module.name)
        else:
            path = module_flowchart_path(module_f, local=True)

        nodes_in_module = []
        for k, v in flowchart.nodes(data=True):
            if v.get("module", "").startswith(module.name) or v.get("type") in {"resource", "main_resource", "raw_resource"}:
                nodes_in_module.append((k, v))

        ids_in_module = []
        nodes_to_keep = []
        for n, v in nodes_in_module:
            if v.get("id") not in ids_in_module:
                nodes_to_keep.append(n)
                ids_in_module.append(v.get("id"))

            for x in list(flowchart.neighbors(n)):
                d = flowchart.nodes(data=True)[x]
                if d.get("id") in nodes_in_module:
                    continue

                if d.get("type") in {"resource", "main_resource", "raw_resource", "gnomad_methods"}:
                    nodes_to_keep.append(x)
                    ids_in_module.append(d.get("id"))

        sub_fc = flowchart.subgraph(nodes_to_keep)
        #if module.ispkg:
        #    pkg_fc = QCFlowchart(sub_fc)
        #    nodes_to_remove = []
        #    for n, v in sub_fc.nodes(data=True):
        #        if v.get("type") not in {"step", "resource", "main_resource"}:
        #            nodes_to_remove.append(n)

        #    for n in nodes_to_remove:
        #        print(f"Removing node {n}")
        #        pkg_fc.remove_node(n)

        #    sub_fc = pkg_fc.get_graph()
        #    print()

        write_mermaid_flowchart(sub_fc, path)  #, add_module_subgraphs=module.ispkg)


def find_resource_node_creator_module(flowchart, resource_node):
    for n in list(flowchart.predecessors(resource_node)):
        d = flowchart.nodes(data=True)[n]
        if d.get("type") == "step":
            return d.get("module")
        creator = find_resource_node_creator_module(flowchart, n)
        if creator is not None:
            return creator

    return None


def find_input_modules(flowchart, resource_node, prev_nodes=None):
    input_steps = []
    if prev_nodes is None:
        prev_nodes = []
    for n in list(flowchart.neighbors(resource_node)):
        d = flowchart.nodes(data=True)[n]
        if d.get("type") == "step":
            input_steps.append(d.get("module"))
        elif n not in prev_nodes:
            input_steps.extend(find_input_modules(flowchart, n, prev_nodes + [resource_node]))

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
