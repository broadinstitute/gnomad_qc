FROM us.gcr.io/broad-gatk/gatk:4.2.6.1

LABEL maintainer="gnomAD Team"
LABEL description="GATK 4.2.6.1 with BCFtools 1.10.2 for VQSR ApplyRecalibration"

# Fix expired Google Cloud SDK GPG key, then install dependencies for BCFtools
RUN rm -f /etc/apt/sources.list.d/google-cloud-sdk.list \
    && apt-get update && apt-get install -y \
    autoconf \
    automake \
    make \
    gcc \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install BCFtools 1.10.2
RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 \
    && tar -xjf bcftools-1.10.2.tar.bz2 \
    && cd bcftools-1.10.2 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf bcftools-1.10.2 bcftools-1.10.2.tar.bz2

# Install tabix (part of htslib, needed for indexing VCFs)
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 \
    && tar -xjf htslib-1.10.2.tar.bz2 \
    && cd htslib-1.10.2 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf htslib-1.10.2 htslib-1.10.2.tar.bz2

# Verify installations
RUN gatk --version && bcftools --version && tabix --version
