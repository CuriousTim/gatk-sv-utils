FROM ubuntu:24.04

ARG BCFTOOLS_VERSION="1.21"
ARG BCFTOOLS_URI="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		curl \
		ca-certificates \
		bzip2 \
		xz-utils \
		build-essential \
		zlib1g-dev \
		libbz2-dev \
		liblzma-dev \
		libcurl4-gnutls-dev \
		libssl-dev \
	&& rm -rf /var/lib/apt/lists/*

RUN curl -L -o bcftools.tar.bz2 "${BCFTOOLS_URI}" \
  && tar -jxf bcftools.tar.bz2 \
  && cd "bcftools-${BCFTOOLS_VERSION}" \
  && ./configure --prefix=/usr/local \
  && make \
  && make install \
  && cp "bcftools-${BCFTOOLS_VERSION}/htslib"*/{bgzip,tabix} /usr/local/bin \
  && cd .. \
  && rm -fr "bcftools-${BCFTOOLS_VERSION}" bcftools.tar.bz2

RUN mkdir /opt/task-scripts
COPY task-scripts /opt/task-scripts

CMD ["bash"]
