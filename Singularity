Bootstrap: docker
From: continuumio/miniconda3

# sudo singularity build snakemake Singularity

%files
    snakemake_tutorial.scif

%environment
    PATH=/opt/conda/bin:$PATH
    export PATH

%post
    apt-get -y install build-essential

    # Install scif, snakemake
    /opt/conda/bin/pip install scif 
    /opt/conda/bin/pip install snakemake==4.4.0
    /opt/conda/bin/scif install /snakemake_tutorial.scif

%runscript
    exec scif "$@"
