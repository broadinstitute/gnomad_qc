# Overlay on gnomad-qc-deps:20260710 that runs a development Hail build for the v5
# info HT runs: Hail main at 882835b33794a5e05e8b8812b268400b687195ef, built
# 2026-10-15, versioned 0.2.139 but 84 commits ahead of the 0.2.139 release.
#
# 882835b itself is only a CI commit ([batch] db testing framework, hail-is/hail#15794)
# -- it is where main sat when the build was cut, not the reason for it. The fix this
# image exists for is somewhere in 0.2.139..882835b; `git log 0.2.139..882835b` in
# hail-is/hail is the list. Retire this file and info_ht_hail_dev_venv.pins once a
# Hail release carries it, and rebuild info_ht.Dockerfile on that release instead.
#
# That build requires Python >= 3.12 (hail-is/hail#15663 shifted supported versions
# from 3.10-3.13 to 3.12-3.14; the wheel METADATA says Requires-Python: >=3.12) and
# the base's /opt/venv is 3.11. So this builds a second venv at /opt/venv-312 from
# info_ht_hail_dev_venv.pins, adds the wheel with --no-deps, and puts that venv first
# on PATH. The 3.11 venv stays intact but shadowed.
#
# The wheel is a build artifact and is deliberately not in the repo. It was published
# to a 30-day Hail bucket that expires around 2026-11-14, so a copy is mirrored to:
#
#   gs://gnomad-resources/hail-dev-wheels/0.2.139-882835b33794/hail-0.2.139-py3-none-any.whl
#   (97002008 bytes, md5 UGHnomAsq7//TywWfowBLA==)
#
# The QoB JAR is inside that wheel at hail/backend/hail-all-spark.jar; it does not
# need to be mirrored separately. Failing both, rebuild Hail at 882835b from source.
#
# The build context has to be staged, since the wheel is not in the repo. From the
# repo root:
#
#   CTX=$(mktemp -d)
#   cp gnomad_qc/v5/annotations/info_ht_hail_dev.Dockerfile $CTX/Dockerfile
#   cp gnomad_qc/v5/annotations/info_ht_hail_dev_venv.pins $CTX/
#   gsutil cp gs://gnomad-resources/hail-dev-wheels/0.2.139-882835b33794/hail-0.2.139-py3-none-any.whl $CTX/
#   gcloud builds submit \
#   --tag us-central1-docker.pkg.dev/broad-mpg-gnomad/images/gnomad-qc-deps:YYYYMMDD-hail<short-sha>-py312 \
#   --project broad-mpg-gnomad \
#   $CTX
#
# --no-deps means the wheel's own dependency metadata is never consulted: if a dev
# build adds or bumps a dependency, pin it in info_ht_hail_dev_venv.pins as well or it
# will be silently missing at runtime. Regenerate those pins by resolving the repo root
# requirements.txt against the target Python, so the 3.12 venv matches what the 3.11
# base image resolved to. The file is named .pins rather than requirements*.txt so
# GitHub's dependency graph does not read it as a manifest for the whole repo; pip
# install -r takes any filename.
#
# Because this venv is the one on PATH, its pins are what actually run -- including
# gnomad_methods, frozen to an exact commit rather than the base image's @main. Local
# gnomad_methods changes therefore do NOT reach a run until the pins are regenerated
# and this image rebuilt, unless you mount over the install at submit time:
#
#   -v /path/to/gnomad_methods/gnomad:/pkg/gnomad \
#   --env PYTHONPATH=/pkg
#
# See info_ht.Dockerfile's header for why that works.
#
# :20260710 was built by info_ht.Dockerfile with 0.2.137 in place of 0.2.139 (FROM
# hailgenetics/hail:0.2.137-py3.11, hail==0.2.137); building that file as committed
# produces the equivalent of :20260826. It makes no difference to what runs here --
# the 3.11 venv is shadowed on PATH, and the base only supplies the OS, g++, git and
# the OpenJDK 11 JRE -- but the tag is not reproducible from the current Dockerfile.
FROM us-central1-docker.pkg.dev/broad-mpg-gnomad/images/gnomad-qc-deps:20260710

COPY hail-0.2.139-py3-none-any.whl /tmp/hail-0.2.139-py3-none-any.whl
COPY info_ht_hail_dev_venv.pins /tmp/info_ht_hail_dev_venv.pins

RUN apt-get update \
  && apt-get install -y --no-install-recommends python3.12 python3.12-venv python3.12-dev \
  && rm -rf /var/lib/apt/lists/* \
  && python3.12 -m venv /opt/venv-312 \
  && /opt/venv-312/bin/pip install --no-cache-dir -r /tmp/info_ht_hail_dev_venv.pins \
  && /opt/venv-312/bin/pip install --no-cache-dir --no-deps /tmp/hail-0.2.139-py3-none-any.whl \
  && rm /tmp/hail-0.2.139-py3-none-any.whl /tmp/info_ht_hail_dev_venv.pins

ENV VIRTUAL_ENV=/opt/venv-312
ENV PATH=/opt/venv-312/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Catch a broken pin set or the wrong interpreter at build time, not in a run.
RUN python -c "import sys, hail, gnomad.utils.annotations; assert sys.version_info >= (3, 12), sys.version; print('hail', hail.__version__)"
