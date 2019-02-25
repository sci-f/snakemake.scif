Bootstrap: docker
From: continuumio/miniconda3

# sudo singularity build snakemake.simg Singularity

%files
    snakemake_tutorial.scif
    Snakefile
    config.yaml

%environment
    PATH=/opt/conda/bin:$PATH
    export PATH

%post
    apt-get update && apt-get -y install build-essential valgrind time

    # Install snakemake
    /opt/conda/bin/conda install --yes -c bioconda -c conda-forge snakemake==4.4.0
    /opt/conda/bin/pip install docutils==0.14

    # Install scif and scif-apps
    /opt/conda/bin/pip install scif 
    /opt/conda/bin/scif install /snakemake_tutorial.scif
    
%runscript
    exec scif "$@"
