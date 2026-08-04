FROM hailgenetics/hail:0.2.128-py3.11

# Install build tools needed by annoy (a gnomad_methods dependency).
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ curl && \
    rm -rf /var/lib/apt/lists/*

# Install gnomad_methods and its transitive deps into the same env as hail
# (system pip at /usr/local/bin/pip -> the python3.10 dist-packages that carry
# hail 0.2.128 in this image; the "-py3.11" tag is actually python 3.10).
RUN pip install --quiet \
    https://github.com/broadinstitute/gnomad_methods/archive/main.tar.gz

# gnomad_methods' deps pull setuptools >= 81, which dropped pkg_resources -- and
# hail 0.2.128 imports pkg_resources at startup. Pin setuptools back below 81 (still
# ships pkg_resources) so `import hail` works in the container.
RUN pip install --quiet "setuptools<81"
