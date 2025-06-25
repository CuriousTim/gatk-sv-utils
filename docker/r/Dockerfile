FROM r-base:4.5.1

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    openssl-dev \
  && rm -rf /var/lib/apt/lists/*

RUN R -e 'withCallingHandlers(install.packages("remotes"), warning = function(w) stop(w))'

RUN --mount=type=bind,source=docker/r/install.R,target=/tmp/install.R \
  Rscript /tmp/install.R

RUN R -e 'remove.packages("remotes")'

RUN mkdir -p /opt/gatk-sv-utils/scripts
COPY src/scripts/r /opt/gatk-sv-utils/scripts

CMD ["bash"]
