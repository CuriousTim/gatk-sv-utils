FROM ubuntu:24.04

ARG BCFTOOLS_VERSION="1.21"
ARG BCFTOOLS_URI="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

ARG HTSLIB_VERSION="1.21"
ARG HTSLIB_URI="https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"

ARG BEDTOOLS_VERSION="2.31.1"
ARG BEDTOOLS_URI="https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz"

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
		libdeflate-dev \
                gawk \
                python3-minimal \
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

RUN mkdir /opt/task_scripts
COPY --chmod=755 task_scripts /opt/task_scripts

RUN mkdir -p /opt/gatk-sv-utils/scripts
COPY src/scripts /opt/gatk-sv-utils/scripts

ENV AWKPATH='.:/usr/local/share/awk/pongo'
RUN mkdir -p /usr/local/share/awk/pongo
COPY src/pongo /usr/local/share/awk/pongo

CMD ["bash"]
