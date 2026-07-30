FROM rocker/r-ver:4.5.1

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    autoconf \
    automake \
    git \
    libcurl4-openssl-dev \
    libbz2-dev \
    libdeflate-dev \
    liblzma-dev \
    libssl-dev \
    libuv1-dev \
    libxml2-dev \
    pkgconf \
    zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

RUN cd /opt \
  && git clone --recurse-submodules 'https://github.com/samtools/htslib' \
  && cd htslib \
  && git checkout 'd94071f3dd9613ab66f177e1237660855a56f0ed' \
  && autoreconf -i \
  && ./configure --prefix=/usr/local \
    --enable-libcurl \
    --enable-gcs \
    --with-libdeflate \
  && make \
  && make install \
  && cd /opt \
  && rm -rf htslib \
  && apt-get remove -y git autoconf automake

RUN R -e 'withCallingHandlers(install.packages("remotes"), warning = function(w) stop(w))'

RUN --mount=type=bind,source=docker/r/install.R,target=/tmp/install.R \
  Rscript /tmp/install.R

RUN R -e 'remotes::install_github("talkowski-lab/gorilla", build_opts = c("--without-keep.source", "--no-help"))'

RUN R -e 'remove.packages("remotes")'

RUN mkdir -p /opt/gatk-sv-utils/scripts
COPY src/scripts/r /opt/gatk-sv-utils/scripts

CMD ["bash"]
