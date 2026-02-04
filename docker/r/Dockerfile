FROM rocker/r-ver:4.5.1

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

RUN R -e 'withCallingHandlers(install.packages("remotes"), warning = function(w) stop(w))'

RUN --mount=type=bind,source=docker/r/install.R,target=/tmp/install.R \
  Rscript /tmp/install.R

RUN --mount=type=bind,source=src/gorilla,target=/tmp/gorilla \
  cd /tmp \
  && R CMD build gorilla \
  && R CMD INSTALL --no-docs --no-html --no-help --without-keep.source \
       gorilla_0.0.1.tar.gz \
  && rm gorilla_0.0.1.tar.gz

RUN R -e 'remove.packages("remotes")'

RUN mkdir -p /opt/gatk-sv-utils/scripts
COPY src/scripts/r /opt/gatk-sv-utils/scripts

CMD ["bash"]
