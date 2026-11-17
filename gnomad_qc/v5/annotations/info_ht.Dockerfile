# Dependencies-only image for the v5 info HT runs (generate_variant_qc_annotations.py
# --generate-ac-info-ht) on Hail 0.2.140 (official release,
# https://hail.is/docs/0.2/change_log.html). gnomad_qc code itself is NOT
# installed -- runs mount it with -v .../gnomad_qc:/pkg/gnomad_qc
# --env PYTHONPATH=/pkg.
#
# Built FROM the official hailgenetics/hail:0.2.140-py3.12 image and pushed as
# gnomad-qc-deps:20261117-hail0.2.140-py312. Replaces
# gnomad-qc-deps:20260710-hail882835b-py312, the dev-wheel image some
# production runs used before 0.2.140 shipped.
#
# gnomad_methods IS installed, pinned at build time (@main via requirements.txt
# here). To run local gnomad_methods changes without a rebuild, mount over the
# install:
#
#   -v /path/to/gnomad_methods/gnomad:/pkg/gnomad \
#   --env PYTHONPATH=/pkg
#
# PYTHONPATH precedes site-packages on sys.path, so the checkout wins. How
# completely depends on the checkout: with a gnomad/__init__.py it is a regular
# package and hides the install outright; without one it is a PEP 420
# namespace package, __path__ holds both directories, and only the modules
# present in the checkout are shadowed. Either way the installed dist-info
# stays, so pip keeps reporting the pinned commit.
#
# Tables written by runs on this image use this Hail's table format; anything
# that reads them (local envs, other images) must run the same or a newer
# Hail.
#
# Build context must be the REPO ROOT (requirements.txt lives there), and
# pushes go to Artifact Registry. No local docker/podman needed, from the repo
# root:
#
#   gcloud builds submit \
#   --project broad-mpg-gnomad \
#   --config gnomad_qc/v5/annotations/info_ht_cloudbuild.yaml \
#   --substitutions _TAG=us-central1-docker.pkg.dev/broad-mpg-gnomad/images/gnomad-qc-deps:YYYYMMDD-hail0.2.140-py312 \
#   .
FROM hailgenetics/hail:0.2.140-py3.12

# annoy (a gnomad_methods dep) has no py3.12 wheel and compiles from source, so
# it needs a C++ compiler at build time.
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends g++ \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/requirements.txt
# aiohttp 3.14 breaks Hail's GCS uploads over 8 MiB.
RUN pip install --no-cache-dir -r /tmp/requirements.txt hail==0.2.140 aiohttp==3.13.2

# Assert the Hail version so a base/pin mismatch fails the build, not a run.
RUN python -c "import hail, gnomad, gnomad.utils.annotations; assert hail.__version__.startswith('0.2.140'), hail.__version__; print('gnomad_methods OK:', gnomad.__file__)"
