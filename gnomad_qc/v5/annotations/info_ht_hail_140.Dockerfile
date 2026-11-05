# Overlay on gnomad-qc-deps:20260710 that runs Hail 0.2.140 (official release,
# https://hail.is/docs/0.2/change_log.html) for the v5 info HT runs. Supersedes
# info_ht_hail_dev.Dockerfile once 0.2.140 is confirmed to carry the fix that
# build existed for -- retire that file and info_ht_hail_dev_venv.pins at that
# point.
#
# 0.2.140 requires Python 3.12+ (hail-is/hail#15663) and the base's /opt/venv is
# 3.11, so this builds a second venv at /opt/venv-312 from
# info_ht_hail_dev_venv.pins (shared with the hail-dev overlay -- deps besides
# hail itself are unchanged) and puts that venv first on PATH. The 3.11 venv
# stays intact but shadowed.
#
# Because this venv is the one on PATH, its pins are what actually run --
# including gnomad_methods, frozen to an exact commit rather than the base
# image's @main. Mount over the install to run local gnomad_methods changes
# without a rebuild:
#
#   -v /path/to/gnomad_methods/gnomad:/pkg/gnomad \
#   --env PYTHONPATH=/pkg
#
# See info_ht.Dockerfile's header for why that works.
#
# Build (no local docker needed), from the repo root:
#
#   CTX=$(mktemp -d)
#   cp gnomad_qc/v5/annotations/info_ht_hail_140.Dockerfile $CTX/Dockerfile
#   cp gnomad_qc/v5/annotations/info_ht_hail_dev_venv.pins $CTX/
#   gcloud builds submit \
#   --tag us-central1-docker.pkg.dev/broad-mpg-gnomad/images/gnomad-qc-deps:YYYYMMDD-hail0.2.140-py312 \
#   --project broad-mpg-gnomad \
#   $CTX
FROM us-central1-docker.pkg.dev/broad-mpg-gnomad/images/gnomad-qc-deps:20260710

COPY info_ht_hail_dev_venv.pins /tmp/info_ht_hail_dev_venv.pins

RUN apt-get update \
  && apt-get install -y --no-install-recommends python3.12 python3.12-venv python3.12-dev \
  && rm -rf /var/lib/apt/lists/* \
  && python3.12 -m venv /opt/venv-312 \
  && /opt/venv-312/bin/pip install --no-cache-dir -r /tmp/info_ht_hail_dev_venv.pins \
  && /opt/venv-312/bin/pip install --no-cache-dir hail==0.2.140 \
  && rm /tmp/info_ht_hail_dev_venv.pins

ENV VIRTUAL_ENV=/opt/venv-312
ENV PATH=/opt/venv-312/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Catch a broken pin set or the wrong interpreter at build time, not in a run.
RUN python -c "import sys, hail, gnomad.utils.annotations; assert sys.version_info >= (3, 12), sys.version; assert hail.__version__.startswith('0.2.140'), hail.__version__; print('hail', hail.__version__)"
