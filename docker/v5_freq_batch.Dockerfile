FROM hailgenetics/hail:0.2.137-py3.11

# g++ is only needed at build time (annoy, a gnomad_methods dependency, compiles
# from source); install, build, and remove it in one layer to keep the image small.
# curl stays: the relay setup command uses it to pull gnomad_methods / gnomad_qc.
USER root
RUN apt-get update && apt-get install -y --no-install-recommends g++ curl && \
    pip install --quiet --no-cache-dir \
    https://github.com/broadinstitute/gnomad_methods/archive/main.tar.gz && \
    apt-get purge -y g++ && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

# Fail the build (not a runtime job) if hail/gnomad don't import in the container
# or the baked Hail drifts from the version this image is tagged with.
RUN python3 -c "import hail; assert hail.__version__.startswith('0.2.137'), hail.__version__; print('hail', hail.__version__); import gnomad; print('gnomad import ok')"
