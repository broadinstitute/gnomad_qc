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

**This is a pipeline repo, not a library or an installable package.** Nothing
here is published to PyPI and nothing outside this repo imports it. Code is run
as scripts, not consumed as an API. The practical consequence: internal
refactors are cheap here (there are no downstream consumers to break), and
reusable logic should not accumulate here at all — if it's genuinely general, it
belongs in [gnomad_methods](https://github.com/broadinstitute/gnomad_methods).

**Do not over-optimize for style.** Correctness and efficiency are what matter
in this repo; pipeline code is allowed to be uglier than library code. Don't
spend effort on cosmetic restructuring, don't refactor working pipeline code for
elegance, and don't let style polish displace attention from whether the step
produces the right output at an acceptable cost. Formatting is enforced by
tooling (below) — beyond that, leave it alone.

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

### Per-release versioning

Code is split by release: `gnomad_qc/v2/`, `v3/`, `v4/`, `v5/`. Each release
directory repeats the same structure, laid out roughly in pipeline order:

| Subpackage | Purpose |
|------------|---------|
| `resources/` | Versioned input/output paths and resource objects. Single source of truth for where data lives in GCS. |
| `sample_qc/` | Hard filters, sex/genetic-ancestry inference, relatedness, outlier filtering, metadata HT. |
| `annotations/` | Info HT, coverage, frequencies, Ensembl VEP, in-silico predictors, variant QC features. |
| `variant_qc/` | Random forest / VQSR training, evaluation, final filter. |
| `create_release/` | Combine annotations, prepare the VCF release, validity checks. |

Not every version has every subpackage (v5 has no `variant_qc/` yet; v4 also has
`assessment/`, `analyses/`, `plots/`), so check before assuming a path exists.
Because the structure repeats, **the equivalent file in the previous version is
almost always the best starting point** for new work.

`v3/README.md` documents the v3 pipeline step by step and is the most complete
prose description of the overall flow. `v4/README.md` is minimal;
`v4/sample_qc/README.md` covers cuKING setup.

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

### Script invocation pattern

Pipeline scripts expose argparse step flags (e.g. `--create-info-ht`,
`--generate-trio-stats`, `--run-pc-relate`). A single script typically implements
many discrete steps, and `main(args)` dispatches by flag — you run one step at a
time, not the whole script.

When adding a step, follow the existing pattern:

1. Add a CLI flag for it.
2. Add a corresponding output resource in `resources/<area>.py`.
3. Gate I/O with `check_resource_existence` (`gnomad_qc/resource_utils.py`) so
   reruns fail fast on missing inputs and refuse to clobber existing outputs
   unless `--overwrite` is passed.

**Always support a test mode** on a new step, and run it that way — small or
mock data, temp output paths — before anything runs genome-wide. Most scripts use
`--test`; some also or instead take `--test-n-partitions`, `--test-partitions`,
`--test-gene`, or `--use-test-dataset`. Match the neighboring script.

### How code is executed

Steps are submitted to a cluster or batch, not run locally on real data:

```bash
hailctl dataproc submit <cluster> gnomad_qc/v4/annotations/generate_freq.py --run-freq-and-dense-annotations --use-test-dataset
```

```bash
python gnomad_qc/v5/annotations/generate_variant_qc_annotations.py --environment batch --create-info-ht --test
```

Exploratory work happens in Jupyter notebooks on a Dataproc cluster
(`hailctl dataproc connect <cluster> nb`). Small checks can run against Hail's
local Spark backend, streaming data from GCS via the storage connector.

Cost conventions: read from the public gnomAD buckets where possible to avoid
requester-pays and egress charges, use **us-central1-b**, prefer autoscaling
clusters, and prefer preemptible workers for initial runs.

## Code Style

- black (`preview = true`, py39–py311), isort (`black` profile, `gnomad` is
  third-party), autopep8. `cuKING` is excluded from all formatting/linting;
  `v2` is excluded from pydocstyle.
- All functions have type annotations. Use `Optional[...]` when the default is
  `None`; never use mutable defaults. Hail types: `hl.Table`,
  `hl.MatrixTable`, `hl.vds.VariantDataset`, `hl.expr.*Expression`.
- Sphinx-style docstrings (pydocstyle, pep257 convention) on every module,
  class, and function. The summary line goes on the line *after* the opening
  `"""` — both styles exist in the repo, but that one dominates. Document
  `:param:` and `:return:`; don't document `:raises:`.

```python
def my_function(ht: hl.Table, max_af: Optional[float] = None) -> hl.Table:
    """
    Short imperative summary line.

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

- **Lazy evaluation drives most surprises.** Anything that pulls Hail data into
  Python (`.count()`, `.collect()`, `hl.eval()`, `.show()`, `.aggregate()`)
  forces the whole upstream plan to execute. Filtering, then `.count()` for a log
  line, then writing does the work twice — **write or checkpoint first, then log**
  against the materialized result.
- `checkpoint(new_temp_file(...))` after expensive computations feeding
  multiple downstream uses; it breaks the query plan. After a checkpoint,
  `.count()` is cheap — otherwise never call `.count()` just for logging on
  large tables. Checkpoints cost compute and storage at gnomAD scale, so add them
  where they actually pay for themselves (joins, repeated large computations),
  not by default.
- **Avoid shuffles.** They're expensive, and shuffle failures rarely name the
  operation that caused them in the stack trace. Re-keying a table shuffles it.
- **Partitioning:** large datasets need many partitions, small ones don't, and
  joins benefit from more. Joining a many-partition table to a few-partition
  table is pathologically slow. After aggressively filtering a large table,
  `.naive_coalesce(N)` to avoid empty-partition skew in downstream aggregations.
  **Never use `repartition()`** — use `naive_coalesce()` to reduce partitions, and
  repartition on read to increase them.
- **Aggregations are costly.** If a step needs many, build one aggregator
  expression (a single `hl.struct` of `hl.agg.*` fields) and make one pass.
- **Hail defaults to GRCh37.** Set the reference genome explicitly; gnomAD v3+ is
  GRCh38.
- Check optional numeric params with `is not None` — truthiness drops valid
  `0`/`0.0` values.
- Check field existence with `field_name in ht.row`; Tables have no `.get()`.
- `ht.aggregate(expr, _localize=False)` keeps results as Hail expressions for
  downstream use without a Python round-trip. It's a private argument and the
  Hail team has advised against relying on it — confirm with them before adding
  new uses.

### Hail versions

- **Hail data is backwards-compatible, not forwards-compatible**: files written
  by a newer Hail cannot be read by an older one. This bites when a cluster and a
  local environment are on different versions.
- **Performance differs between versions** — there was a significant regression
  from 0.2.130 to 0.2.131. Newer isn't automatically better; weigh it against the
  features you need.
- **GCP/GCS argument conventions change between versions and can change costs.**
  A past change in how requester-pays buckets were specified caused egress
  charges on *all* writes. Treat cloud config changes as cost-affecting until
  shown otherwise.

## Testing and Validation

There is no unit test suite; CI runs linting and the docs build only.
Validation is empirical: run the affected pipeline step with `--test` (small
test datasets, temp output paths) and inspect the result. Don't claim a change
is verified on the basis of lint passing.

A passing `--test` run is necessary but not sufficient. Out-of-memory failures,
shuffle failures, and partition skew appear only at full scale, and a step that
works on a test gene can fail genome-wide. Say what you actually ran and at what
scale, rather than calling a change verified.

## Never Do This

- **Never overwrite an existing HT/MT/VDS when testing.** Write a new copy to a
  temp path instead. `--overwrite` on a production output path is destructive and
  the data is expensive to regenerate.
- **Never delete data outside temp buckets**, and never delete files that are
  referenced or built by the repo.
- **Never remove code without asking first.**
- **Never use `repartition()`** — see the Hail notes above.
- **Never hardcode a `gs://` path** — go through a `get_*()` resource function.

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
