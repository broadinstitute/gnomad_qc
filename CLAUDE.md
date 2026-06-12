# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project Overview

gnomad_qc contains the QC and variant annotation pipelines for the Genome
Aggregation Database (gnomAD). Each major release has its own versioned package
(`v2`–`v5`); `v5` (All of Us integration) is under active development, `v4`
still receives fixes, and older versions are historical record — match the
conventions of the version you're editing. Pipeline scripts are standalone
`argparse` CLI tools that read/write Hail Tables, MatrixTables, and VDS objects
on GCS.

Two tightly coupled dependencies:

- [gnomad_methods](https://github.com/broadinstitute/gnomad_methods) — shared
  resource classes and utility functions. Installed from git `main` (not PyPI)
  and co-developed with this repo; a fix may belong there rather than here.
- [Hail](https://hail.is/) — the compute framework. Behavior can differ across
  Hail versions, so apparent bugs in pipeline code are sometimes version
  mismatches in the local environment.

## Finding Code and APIs

Function names, signatures, and locations change between releases — discover
them at point of use rather than relying on memory or this file:

- **Dataset paths**: every dataset is referenced through a `get_*()` resource
  function in the version's `resources/` package, organized by domain
  (`basics.py`, `sample_qc.py`, `variant_qc.py`, `annotations.py`, `meta.py`,
  `release.py`, `constants.py`). Grep `def get_` in the relevant module. Never
  hardcode `gs://` paths in pipeline code.
- **Shared utilities**: check gnomad_methods (`gnomad.utils.*`,
  `gnomad.sample_qc.*`, `gnomad.variant_qc.*`, `gnomad.resources.*`) before
  writing a new helper — it probably exists. A local checkout often sits next
  to this repo at `../gnomad_methods` and may be the editable install; verify
  which copy imports resolve to before debugging it.
- **Precedent**: the previous version package usually contains the pattern the
  next one needs (v4 → v5). When adding a pipeline step, find and read the
  closest analog in the current and prior version first.
- **Docs**: https://broadinstitute.github.io/gnomad_qc/ and the gnomad_methods
  docs are generated from docstrings and searchable.

## Commands

### Setup
```
pip install -r requirements.txt
pip install -r requirements-dev.txt
pre-commit install
```

### Formatting (auto-applied by pre-commit hooks)
```
black gnomad_qc
isort gnomad_qc
autopep8 --in-place gnomad_qc
```

### Linting (matches CI)
```
black --check gnomad_qc
isort --check-only gnomad_qc
autopep8 --exit-code --diff gnomad_qc
./lint --disable=W                                    # pylint, warnings disabled
pydocstyle --match-dir='(?!v2|cuKING).*' gnomad_qc    # skips v2 and cuKING
```

### Docs
```
./docs/build.sh    # sphinx with -W (warnings are errors); CI publishes to gh-pages
```

## Architecture

### Resource pattern

Resource modules define `get_*()` functions returning typed resource objects
(`TableResource`, `MatrixTableResource`, `VariantDatasetResource`, and their
`Versioned*` wrappers from `gnomad.resources.resource_utils`). Getters
typically take `version`, `test`, and (in v5) `environment` parameters; `test`
routes output to temp paths with limited retention. Read data via the resource
(`.ht()`, `.mt()`, `.vds()`), write to its `.path`.

### Environments and execution (v5)

v5 runs in three environments — `"rwb"` (All of Us Researcher Workbench),
`"batch"` (Hail Query-on-Batch, the primary backend), and `"dataproc"` —
selected by the `--environment` flag. `v5/resources/basics.py` owns
environment validation, GCS bucket routing, and Hail initialization; scripts
call its init helper before any Hail work. Some resource-existence checks are
intentionally skipped under `"batch"` (users lack read access to the read-only
bucket). Earlier versions submit jobs with `hailctl dataproc`.

### Pipeline script anatomy

Scripts follow a common skeleton — copy it from a neighboring script in the
same version rather than from memory:

- Module docstring describing the pipeline step.
- Module-level logger: `logging.getLogger(<script name>)` with
  `logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")`.
- Helper functions (typed, docstringed) that transform Tables; `main(args)`
  dispatches steps gated by boolean `--step-name` flags and owns all
  read/write/checkpoint I/O.
- Standard flags: `--overwrite`, `--test`, and in v5 `--environment`.
- A `finally` block that copies the Hail log to a logging path.
- Parser built in a `get_script_argument_parser()`-style function, invoked
  under `if __name__ == "__main__":`.

v4 scripts manage inter-step dependencies with `PipelineStepResourceCollection`
/ `check_resource_existence` from `gnomad_qc/resource_utils.py`, and optionally
wrap `main` in `slack_notifications(...)` (token from `gnomad_qc/slack_creds.py`,
channel from `--slack-channel`). v5 favors simpler `main()` bodies.

## Code Style

- black (`preview = true`, py39–py311), isort (`black` profile, `gnomad` is
  third-party), autopep8. `cuKING` is excluded from all formatting/linting;
  `v2` is excluded from pydocstyle.
- All functions have type annotations. Use `Optional[...]` when the default is
  `None`; never use mutable defaults. Hail types: `hl.Table`,
  `hl.MatrixTable`, `hl.vds.VariantDataset`, `hl.expr.*Expression`.
- Sphinx-style docstrings (pydocstyle, pep257 convention) on every module,
  class, and function:

```python
def my_function(ht: hl.Table, max_af: Optional[float] = None) -> hl.Table:
    """Short imperative summary line.

    Extended description if needed. Use ``double backticks`` for field/code
    references and ``.. note::`` blocks for caveats.

    :param ht: Input Table with ``freq`` field.
    :param max_af: Maximum allele frequency threshold. Default is None.
    :return: Table with ``filters`` field added.
    """
```

- Reusable transformations should be Table/expression in, Table/expression out
  with no file I/O — and if genuinely general, belong in gnomad_methods, not
  here.

## Hail Notes

- `checkpoint(new_temp_file(...))` after expensive computations feeding
  multiple downstream uses; it breaks the query plan. After a checkpoint,
  `.count()` is cheap — otherwise never call `.count()` just for logging on
  large tables.
- After aggressively filtering a large table, `.naive_coalesce(N)` to avoid
  empty-partition skew in downstream aggregations.
- Check optional numeric params with `is not None` — truthiness drops valid
  `0`/`0.0` values.
- Check field existence with `field_name in ht.row`; Tables have no `.get()`.
- `ht.aggregate(expr, _localize=False)` keeps results as Hail expressions for
  downstream use without a Python round-trip.

## Testing and Validation

There is no unit test suite; CI runs linting and the docs build only.
Validation is empirical: run the affected pipeline step with `--test` (small
test datasets, temp output paths) and inspect the result. Don't claim a change
is verified on the basis of lint passing.

## Repo Gotchas

- `tmp_*.py` and other `tmp_*` files at the repo root are scratch/debugging
  artifacts, not pipeline code — don't pattern-match off them or import them.
- `cuKING/` (in versions that have it) is a vendored CUDA relatedness tool,
  exempt from all style tooling.
- Machine-specific environment notes (Python paths, Hail version pairing) live
  in the gitignored `CLAUDE.local.md`.

## Maintaining CLAUDE.md

Proactively add gotchas, non-obvious behavior, and conventions discovered
during development. Keep entries concise and durable: document patterns and
how to find things, not inventories of function names or dataset paths — those
churn between releases and rot here.
