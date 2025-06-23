FROM r-base:4.5.1

RUN R -e 'install.packages(c("data.table", "optparse", "BiocManager"), INSTALL_opts = c("--no-docs", "--no-help", "--no-html", "--no-data"))' \
  && R -e 'BiocManager::install(c("Rsamtools", "GenomicRanges"), INSTALL_opts = c("--no-docs", "--no-help", "--no-html", "--no-data"))'

RUN mkdir -p /opt/gatk-sv-utils/scripts
COPY src/scripts/r /opt/gatk-sv-utils/scripts

CMD ["bash"]
