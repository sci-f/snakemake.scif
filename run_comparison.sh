#!/bin/sh

# This script assumes that we start in the same directory as this script,
# and will run and save results for different pipelines (to compare later)

# If running interactively, use HERE=$PWD
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
DATA=$HERE/data
RESULT=$HERE/results

mkdir -p $RESULT

function clean () {
    echo "Cleaning up previous results..."
    rm -rf $HERE/data/sorted_reads
    rm -rf $HERE/data/mapped_reads
    rm -rf $HERE/data/calls
    rm -rf $HERE/data/genome.fa.*
    rm -rf $HERE/data/report.html
    rm -rf $HERE/data/dag.svg
}

function saveResult() {
    container=$1
    output=${RESULT}/${container}
    mkdir -p $output
    cp $DATA/calls/all.vcf $output/calls.vcf
    cp $DATA/report.html $output/report.html
    cp $DATA/genome.fa.* $output/
    cp $DATA/dag.svg $output/
}

# Run for Docker first
clean
echo "Running docker pipeline"
docker build --no-cache -t vanessa/snakemake.scif .
docker run -v $PWD/data:/scif/data -v $PWD:/scif/data/snakemake -it vanessa/snakemake.scif run snakemake all
docker run -v $PWD/data:/scif/data -v $PWD:/scif/data/graphviz_create_dag -it vanessa/snakemake.scif run graphviz_create_dag /scif/data/graphviz_create_dag report.html dag.svg
saveResult docker
sudo chown -R vanessa data

# Singularity
clean
echo "Running singularity pipeline"
sudo singularity build snakemake.simg Singularity
singularity run --bind data/:/scif/data snakemake.simg run snakemake all
singularity run --bind data/:/scif/data snakemake.simg run graphviz_create_dag $PWD report.html dag.svg
saveResult singularity

# Charliecloud
clean
echo "Running charliecloud pipeline"
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif run snakemake all
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif run graphviz_create_dag $PWD report.html dag.svg
saveResult charliecloud

# Runc
clean
echo "Running runc pipeline"
mkdir -p /tmp/runc/rootfs
sudo singularity build --sandbox /tmp/runc/rootfs snakemake.simg
cp configs/runc-config-noninteractive.json /tmp/runc/config.json
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
cp configs/runc-config-dag.json /tmp/runc/config.json
sudo chown -R $USER /tmp/runc/rootfs
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
saveResult runc

# Shifter
clean
echo "Shifter pipeline run on cluster"
echo "shifter --image=vanessa/snakemake.scif --volume `pwd`/data:/scif/data --workdir /scif/data /opt/conda/bin/scif run snakemake all"
