# Dependencies-only image for the v5 info HT runs (generate_variant_qc_annotations.py
# --generate-ac-info-ht): Hail + gnomad_methods, via the repo root requirements.txt.
# gnomad_qc code itself is NOT installed -- runs mount it with
# -v .../gnomad_qc:/pkg/gnomad_qc --env PYTHONPATH=/pkg.
#
# gnomad_methods IS installed, pinned at build time (@main via requirements.txt here;
# an exact commit via info_ht_hail_dev_venv.pins in the hail-dev overlay). To run local
# gnomad_methods changes without a rebuild, mount over the install the same way:
#
#   -v /path/to/gnomad_methods/gnomad:/pkg/gnomad \
#   --env PYTHONPATH=/pkg
#
# PYTHONPATH precedes site-packages on sys.path, so the checkout wins. How completely
# depends on the checkout: with a gnomad/__init__.py it is a regular package and hides
# the install outright; without one it is a PEP 420 namespace package, __path__ holds
# both directories, and only the modules present in the checkout are shadowed. Either
# way the installed dist-info stays, so pip keeps reporting the pinned commit.
#
# Tables written by runs on this image use this Hail's table format; anything that
# reads them (local envs, other images) must run the same or a newer Hail.
#
# Build context must be the REPO ROOT (requirements.txt lives there), and pushes go to
# Artifact Registry. Easiest path (no local docker needed), from the repo root:
#
#   gcloud builds submit \
#   --project broad-mpg-gnomad \
#   --config gnomad_qc/v5/annotations/info_ht_cloudbuild.yaml \
#   --substitutions _TAG=us-central1-docker.pkg.dev/broad-mpg-gnomad/images/gnomad-qc-deps:YYYYMMDD \
#   .
#
# With podman locally instead, again from the repo root:
#
#   IMAGE=us-central1-docker.pkg.dev/broad-mpg-gnomad/images/gnomad-qc-deps:YYYYMMDD
#   podman build \
#   --platform linux/amd64 \
#   -f gnomad_qc/v5/annotations/info_ht.Dockerfile \
#   -t $IMAGE \
#   .
#   gcloud auth print-access-token | podman login -u oauth2accesstoken --password-stdin us-central1-docker.pkg.dev
#   podman push $IMAGE
FROM hailgenetics/hail:0.2.139-py3.11

# annoy (a gnomad_methods dep) has no py3.11 wheel and compiles from source, so
# it needs a C++ compiler at build time.
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends g++ \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt \
    && pip install --no-cache-dir hail==0.2.139

# Assert the Hail version so a base/pin mismatch fails the build, not a run.
RUN python -c "import hail, gnomad, gnomad.utils.annotations; assert hail.__version__.startswith('0.2.139'), hail.__version__; print('gnomad_methods OK:', gnomad.__file__)"
