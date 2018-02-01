FROM continuumio/miniconda3

# docker build -t vanessa/snakemake.scif .

ADD snakemake_tutorial.scif /

ENV PATH /opt/conda/bin:$PATH

RUN apt-get -y install build-essential git

# Install scif, snakemake

RUN /opt/conda/bin/pip install scif && \
    /opt/conda/bin/pip install snakemake==4.4.0 && \
    /opt/conda/bin/scif install /snakemake_tutorial.scif

ENTRYPOINT ["scif"]
