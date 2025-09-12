FROM ubuntu:25.04

ARG BCFTOOLS_VERSION="1.22"
ARG BCFTOOLS_URI="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

ARG HTSLIB_VERSION="1.22.1"
ARG HTSLIB_URI="https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"

ARG BEDTOOLS_VERSION="2.31.1"
ARG BEDTOOLS_URI="https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz"

ARG DUCKDB_VERSION="1.3.2"
ARG DUCKDB_URI="https://github.com/duckdb/duckdb/releases/download/v${DUCKDB_VERSION}/duckdb_cli-linux-amd64.zip"

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    apt-transport-https \
    build-essential \
    bzip2 \
    ca-certificates \
    curl \
    gawk \
    gnupg \
    jq \
    libbz2-dev \
    libcurl4-gnutls-dev \
    libdeflate-dev \
    liblzma-dev \
    libopenblas0-pthread \
    libssl-dev \
    python3-minimal \
    unzip \
    xz-utils \
    zlib1g-dev \
    zstd \
  && rm -rf /var/lib/apt/lists/*

RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
  | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
  | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
  && apt-get update -y \
  && apt-get install -y --no-install-recommends google-cloud-cli \
  && rm -rf /var/lib/apt/lists/*

RUN curl -L -o htslib.tar.bz2 "${HTSLIB_URI}" \
  && tar -jxf htslib.tar.bz2 \
  && cd "htslib-${HTSLIB_VERSION}" \
  && ./configure --prefix=/usr/local \
    --enable-libcurl \
    --enable-gcs \
    --with-libdeflate \
  && make \
  && make install \
  && cd .. \
  && rm -rf "htslib-${HTSLIB_VERSION}" htslib.tar.bz2

RUN ldconfig

RUN curl -L -o bcftools.tar.bz2 "${BCFTOOLS_URI}" \
  && tar -jxf bcftools.tar.bz2 \
  && cd "bcftools-${BCFTOOLS_VERSION}" \
  && ./configure --prefix=/usr/local \
    --with-htslib=/usr/local \
  && make \
  && make install \
  && cd .. \
  && rm -rf "bcftools-${BCFTOOLS_VERSION}" bcftools.tar.bz2

RUN curl -L -o bedtools.tar.gz "${BEDTOOLS_URI}" \
  && tar -zxf bedtools.tar.gz \
  && cd bedtools2 \
  && make prefix=/usr/local \
  && make install \
  && cd .. \
  && rm -rf "bedtools-${BEDTOOLS_VERSION}" bedtools.tar.gz

RUN curl -L -o duckdb_cli-linux-amd64.zip "${DUCKDB_URI}" \
  && unzip duckdb_cli-linux-amd64.zip -d /usr/local/bin \
  && rm duckdb_cli-linux-amd64.zip

RUN mkdir -p /opt/gatk-sv-utils/scripts
COPY src/scripts/base /opt/gatk-sv-utils/scripts

ENV AWKPATH='.:/usr/local/share/awk/pongo'
RUN mkdir -p /usr/local/share/awk/pongo
COPY src/pongo /usr/local/share/awk/pongo

RUN --mount=type=bind,source=src/bonobo,target=/tmp/bonobo,rw \
  cd /tmp/bonobo \
  make all && make install

CMD ["bash"]
