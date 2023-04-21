"""Generic gnomAD QC resource utils."""
from collections import defaultdict
from typing import Any, Dict, List, Optional, Union

from gnomad.utils.file_utils import check_file_exists_raise_error


def check_resource_existence(
    input_step_resources: Optional[Dict[str, List]] = None,
    output_step_resources: Optional[Dict[str, List]] = None,
    overwrite: bool = False,
) -> None:
    """
    Check the existence of all specified input and output resources.

    If any of the input resources (`input_step_resources` values) don't exist, an error
    will be raised indicating which input resources are missing.

    If any of the output resources (`output_step_resources` values) already exist and
    the `overwrite` parameter is not set to True, an error will be raised indicating
    which output resources already exist.

    If no parameters are passed to the function, nothing is done.

    :param input_step_resources: A dictionary with keys as pipeline steps that generate
        input files and the value as a list of the input files to check the existence
        of. Default is None.
    :param output_step_resources: A dictionary with keys as pipeline step that generate
        output files and the value as a list of the output files to check the existence
        of. Default is None.
    :param overwrite: The overwrite parameter used when writing the output files.
        Default is False.
    :return: None.
    """
    # Check if the input resources exist
    if input_step_resources:
        for step, input_resources in input_step_resources.items():
            check_file_exists_raise_error(
                [r if isinstance(r, str) else r.path for r in input_resources],
                error_if_not_exists=True,
                error_if_not_exists_msg=(
                    f"Not all input resources exist. Please add {step} to the command "
                    "line. The following files are missing: "
                ),
            )

    # Check if the output resources exist when `overwrite` is False
    if not overwrite and output_step_resources:
        for step, output_resources in output_step_resources.items():
            check_file_exists_raise_error(
                [r if isinstance(r, str) else r.path for r in output_resources],
                error_if_exists=True,
                error_if_exists_msg=(
                    f"Some of the output resources that will be created by {step} "
                    "already exist and the --overwrite argument was not set. Please "
                    f"rerun {step} with --overwrite. The following files already exist:"
                ),
            )


# TODO: This is a temporary solution to input/output checking for a workflow/pipeline
#  that was added for use in v4. We will need to revisit this class post gnomAD v4
#  MVP because it is fairly complex and there are some external libraries like luigi
#  that could do some of this dependency checking for us.
# TODO: Rethink `add_input_resources`. The current structure results in large keys and
#  there is likely a better way to provide these added inputs.
class PipelineStepResourceCollection:
    """
    A collection of resources used in a specific step of a gnomad_qc pipeline.

    Example usage:
        Assume a pipeline that has 4 pipeline steps (a1, a2, b, and c) and 6 resources
        used and/or created by the steps in the pipeline (r1, a2_r_in, a1_r_out,
        a2_r_out, b_r_out, and c_r_out).
            - Steps a1 and a2 can be run in parallel, and use input from different
              pipelines.
            - Step b uses, as input, the output from both steps a1 and a2.
            - Step c uses, as input, the output from step b, and the resources used by
              a1.

        r1_input = {"r1_pipeline.py --r1_step": {"r1": r1}
        a1 = PipelineStepResourceCollection(
            "--a1",
            output_resources={"a1_r_out": a1_r_out},
            input_resources=r1_input
        )
        a2 = PipelineStepResourceCollection(
            "--a2",
            output_resources={"a2_r_out": a2_r_out},
            input_resources={"a2_r_in_pipeline.py --a2_r_in_step": {"a2_r_in": a2_r_in}}
        )
        b = PipelineStepResourceCollection(
            "--b",
            pipeline_input_steps=[a1, a2],
            output_resources={"b_r_out": b_r_out},
        )
        c = PipelineStepResourceCollection(
            "--c",
            pipeline_input_steps=[b],
            add_input_resources=r1_input,
            output_resources={"b_r_out": b_r_out},
        )

        # Check that input exists and that output doesn't already exist for a1 and a2 .
        a1.check_resource_existence()
        a2.check_resource_existence()

        # Run steps a1 and a2 to produce a1_r_out and a2_r_out.
        out = a1_step_run(a1.r1)
        out.write(a1.a1_r_out.path)
        out = a2_step_run(a2.a2_r_in)
        out.write(a2.a2_r_out.path)

        # Check that input exists and that output doesn't already exist for b.
        b.check_resource_existence()

        # Run step b to produce b_r_out.
        out = b_step_run(b.a1_r_out, b.a2_r_out)
        out.write(b.b_r_out.path)

        # Check that input exists and that output doesn't already exist for c.
        c.check_resource_existence()

        # Run step c to produce c_r_out.
        out = c_step_run(c.r1, c.b_r_out)
        out.write(c_r_out.path)
    """

    def __init__(
        self,
        step_name: str,
        output_resources: Dict[str, Any],
        pipeline_input_steps: List = [],
        input_resources: Optional[Dict[str, Union[Dict[str, Any], List[Any]]]] = None,
        add_input_resources: Optional[
            Dict[str, Union[Dict[str, Any], List[Any]]]
        ] = None,
        pipeline_name: Optional[str] = None,
        overwrite: Optional[bool] = None,
    ) -> None:
        """
        Construct all the necessary attributes for the object.

        :param step_name: Name to be assigned to current step. The recommendation is to
            use the step's argparse call as the name.
        :param pipeline_input_steps: Previous steps of the pipeline that have input
            needed for the current step. Default is an empty List.
        :param output_resources: Dictionary of the pipeline step's output
            resources. Keyed by the name to use as an attribute name and with the
            associated resource as the value.
        :param input_resources: Optional Dictionary of the pipeline step's input
            resources. Keyed by a string used to describe where the resource comes
            from or is created for example "generate_qc_mt.py --generate-qc-meta". The
            Values can be either a list of resources, or a Dictionary keyed by each
            resource's name (use as an attribute name) and with the associated resource
            as the value. By default, the output from steps listed in
            `previous_pipeline_steps` will be used as input resources to this step, but
            if input_resources is supplied, then only those resources will be used.
        :param add_input_resources: Optional Dictionary of additional input resources
            to add for this step on top of either `input_resources` or the output from
            steps listed in `previous_pipeline_steps`. Keyed by a string used to
            describe where the resource comes from or is created for example:
            "generate_qc_mt.py --generate-qc-meta". The Values can be either a list of
            resources, or a Dictionary keyed by each resource's name (use as an
            attribute name) and with the associated resource as the value.
        :param pipeline_name: Optional name of the full pipeline this step belongs to.
        :param overwrite: Whether these resources can be overwritten. Used in
            `check_resource_existence`.
        :return: None.
        """
        self.pipeline_name = pipeline_name
        self.overwrite = overwrite
        self.pipeline_step = step_name
        self.previous_steps = pipeline_input_steps
        self._output_resource_dict = output_resources

        self._add_output_resource_attributes(output_resources)
        output_resources = {self.pipeline_step: list(output_resources.values())}

        self.output_resources = output_resources

        self.input_resources = defaultdict(list)
        if input_resources is None:
            for step in pipeline_input_steps:
                self.input_resources.update(step.output_resources)
                # Add the output of all the previous pipeline steps as attributes to the
                # current object. This allows the resource to be called from this
                # object. For example (from the class example above):
                # b = PipelineStepResourceCollection(
                #    "--b",
                #    pipeline_input_steps=[a1, a2],
                #    output_resources={"b_r_out": b_r_out},
                # )
                # c = PipelineStepResourceCollection(
                #    "--c",
                #    pipeline_input_steps=[b],
                #    add_input_resources=r1_input,
                #    output_resources={"b_r_out": b_r_out},
                # )
                # # Step c can easily be run using the attribute name 'b_r_out' on the
                # # 'c' class instance.
                # c_step_run(c.r1, c.b_r_out)
                self._add_output_resource_attributes(step._output_resource_dict)
        else:
            self.add_input_resources(input_resources)

        if add_input_resources is not None:
            self.add_input_resources(add_input_resources)

    def __getattr__(self, name: str) -> None:
        """Raise AttributeError with modified message if the attribute is not found."""
        AttributeError(
            f"Pipeline step {self.pipeline_step} has no resource attribute named"
            f" {name}!"
        )

    def _add_output_resource_attributes(
        self,
        output_resources: Dict[str, Any],
    ) -> None:
        """
        Add output resources as class attributes.

        :param output_resources: Dictionary of the pipeline step's output resources.
            Keyed by the name to use as an attribute name and with the associated
            resource as the value.
        :return: None.
        """
        for name, resource in output_resources.items():
            setattr(self, name, resource)

    def add_input_resources(
        self,
        input_resources: Dict[str, Union[Dict[str, Any], List[Any]]],
    ) -> None:
        """
        Add input resources to the pipeline step.

        :param input_resources: Dictionary of resources to add as input for this
            pipeline step. Keyed by a string used to describe where the resource comes
            from or is created for example: "generate_qc_mt.py --generate-qc-meta". The
            Values can be either a list of resources, or a Dictionary keyed by each
            resource's name (used as an attribute name) and with the associated resource
            as the value.
        :return: None.
        """
        for step, resources in input_resources.items():
            if isinstance(resources, dict):
                for name, resource in resources.items():
                    setattr(self, name, resource)
                resources = list(resources.values())
            self.input_resources[step].extend(resources)

    def check_resource_existence(self, overwrite: bool = None) -> None:
        """
        Check for existence of input and output resources.

        An error is returned if input resources don't exist, or if overwrite is False
        and an output resource does exist.

        :param overwrite: Whether these resources can be overwritten.
        :return: None.
        """
        if overwrite is None:
            if self.overwrite is None:
                overwrite = False
            else:
                overwrite = self.overwrite
        check_resource_existence(
            input_step_resources=self.input_resources,
            output_step_resources=self.output_resources,
            overwrite=overwrite,
        )


class PipelineResourceCollection:
    """
    A collection of PipelineStepResourceCollections used in a gnomad_qc pipeline.

    Example usage extending the example in `PipelineStepResourceCollection`:

        # Initialize pipeline resource collection
        pipeline = PipelineResourceCollection(
            pipeline_name="example",
            pipeline_resources={"r1", r1},
            overwrite=False,
        )

        # Add all steps to the example pipeline resource collection.
        pipeline.add_steps({"a1": a1, "a2": a2, "b": b, "c": c})

        # Check that input exists and that output doesn't already exist for c.
        pipeline.check_resource_existence("c")

        # Run step c to produce c_r_out.
        out = c_step_run(pipeline.c.r1, pipeline.c.b_r_out)
        out.write(pipeline.c.c_r_out.path)
    """

    def __init__(
        self,
        pipeline_name: str,
        pipeline_steps: Dict[str, PipelineStepResourceCollection] = {},
        pipeline_resources: Optional[Dict[str, Any]] = None,
        overwrite: bool = None,
    ) -> None:
        """
        Construct all the necessary attributes for the object.

        :param pipeline_name: Name of the pipeline.
        :param pipeline_steps: Dictionary of all steps of the pipeline Keyed by a name
            to use for the step (added as an attribute of the object) and with the
            associated PipelineStepResourceCollection as the Value. Default is an empty
            Dictionary.
        :param pipeline_resources: Optional Dictionary of resources used by multiple
            steps of the pipeline. Keyed by the name to use as an attribute name and
            with the associated resource as the value.
        :param overwrite: Whether these resources can be overwritten. Used in
            `check_resource_existence`.
        :return: None.
        """
        self.pipeline_name = pipeline_name
        self.pipeline_steps = {}
        self.add_steps(pipeline_steps)
        self.overwrite = overwrite

        if pipeline_resources is not None:
            for name, resource in pipeline_resources.items():
                setattr(self, name, resource)

    def __getattr__(self, name: str) -> None:
        """Raise AttributeError with modified message if the attribute is not found."""
        AttributeError(
            f"Pipeline {self.pipeline_name} has no resource or pipeline step attribute"
            f" named {name}!"
        )

    def add_steps(self, steps: Dict[str, PipelineStepResourceCollection]) -> None:
        """
        Add PipelineStepResourceCollections to the object.

        :param steps: Dictionary of steps to add to the pipeline Keyed by a name
            to use for the step (added as an attribute of the object) and with the
            associated PipelineStepResourceCollection as the Value.
        :return: None.
        """
        for step_name, step in steps.items():
            step.pipeline_name = self.pipeline_name
            step.overwrite = self.overwrite
            setattr(self, step_name, step)
        self.pipeline_steps.update(steps)

    def check_resource_existence(
        self, step: str, overwrite: Optional[bool] = None
    ) -> None:
        """
        Check for existence of input and output resources.

        An error is returned if input resources don't exist, or if overwrite is False
        and an output resource does exist.

        :param step: Name of step to check resource existence of.
        :param overwrite: Whether these resources can be overwritten. By default, the
            class attribute of overwrite will be used if it exists, otherwise overwrite
            will be set to False.
        :return: None.
        """
        if overwrite is None:
            if self.overwrite is None:
                overwrite = False
            else:
                overwrite = self.overwrite

        # TODO: When revisiting, if this is kept, throw an error message here that
        #  tells you the available steps in the pipeline collection object.
        self.pipeline_steps[step].check_resource_existence(overwrite=overwrite)
