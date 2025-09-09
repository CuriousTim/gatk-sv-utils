FROM python:3.13.7-trixie

ARG PYSAM_VERSION="0.23.3"

RUN pip3 install "pysam==${PYSAM_VERSION}"

RUN mkdir -p /opt/gatk-sv-utils/scripts
COPY src/scripts/python /opt/gatk-sv-utils/scripts

CMD ["bash"]
