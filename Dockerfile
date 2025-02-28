FROM ubuntu:24.04

ARG BCFTOOLS_VERSION="1.21"
ARG BCFTOOLS_URI="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

RUN curl -L -o bcftools.tar.bz2 "${BCFTOOLS_URI}" \
  && tar -jxf bcftools.tar.bz2 \
  && cd "bcftools-${BCFTOOLS_VERSION}" \
  && ./configure --prefix=/usr/local \
  && make \
  && make install \
  && cd .. \
  && rm -fr "bcftools-${BCFTOOLS_VERSION}" bcftools.tar.bz2

COPY task-scripts /opt/

CMD ["bash"]
