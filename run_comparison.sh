#!/bin/bash

# This script assumes that we start in the same directory as this script,
# and will run and save results for different pipelines (to compare later)

# usageL source run_comparison.sh

HERE=$PWD
DATA=$HERE/data
RESULT=$HERE/results
DOBUILD=0

mkdir -p $RESULT

function clean () {
    echo
    echo "Cleaning up previous results..."
    rm -rf $HERE/data/sorted_reads
    rm -rf $HERE/data/mapped_reads
    rm -rf $HERE/data/calls
    rm -rf $HERE/data/genome.fa.*
    rm -rf $HERE/data/report.html
    rm -rf $HERE/data/dag.svg
}

function doBuild() {
   if [ ${DOBUILD} -eq "1" ]
       then 
       $1
   fi 
}

function saveResult() {
    container=$1
    echo "Saving $container results..."
    output=${RESULT}/${container}
    mkdir -p $output
    cp $DATA/calls/all.vcf $output/calls.vcf
    cp $DATA/report.html $output/report.html
    cp $DATA/genome.fa.* $output/
    cp $DATA/dag.svg $output/
    echo
}

function introduce() {
    container=$1
    clean
    cat $HERE/asciinema/$container.txt
    echo
    echo "Running $container pipeline"
    sleep 2
}

cat asciinema/intro.txt
sleep 3

# Run for Docker first
introduce docker
doBuild "docker build --no-cache -t vanessa/snakemake.scif ."
docker run -v $PWD/data:/scif/data -v $PWD:/scif/data/snakemake -it vanessa/snakemake.scif run snakemake all
docker run -v $PWD/data:/scif/data -v $PWD:/scif/data/graphviz_create_dag -it vanessa/snakemake.scif run graphviz_create_dag /scif/data/graphviz_create_dag report.html dag.svg
sudo chown -R vanessa data/
saveResult docker

# Singularity
introduce singularity
doBuild "sudo singularity build snakemake.simg Singularity"
singularity run --bind data/:/scif/data snakemake.simg run snakemake all
singularity run --bind data/:/scif/data snakemake.simg run graphviz_create_dag $PWD report.html dag.svg
saveResult singularity

# Charliecloud
introduce charliecloud
doBuild "ch-build -t snakemake.ch ." 
doBuild "ch-docker2tar snakemake.ch /var/tmp"
doBuild "ch-tar2dir /var/tmp/snakemake.ch.tar.gz /var/tmp"
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif run snakemake all
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif run graphviz_create_dag $PWD report.html dag.svg
saveResult charliecloud

# Runc
introduce runc
mkdir -p /tmp/runc/rootfs
doBuild "sudo singularity build --sandbox /tmp/runc/rootfs snakemake.simg"
cp configs/runc-config-noninteractive.json /tmp/runc/config.json
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
cp configs/runc-config-dag.json /tmp/runc/config.json
sudo chown -R $USER /tmp/runc/rootfs
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
saveResult runc

cat asciinema/goals.txt
sleep 2
cat asciinema/goalsfor.txt
sleep 2
cat asciinema/reproducibility.txt
sleep 2
cat asciinema/portability.txt
sleep 2
cat asciinema/friends.txt
sleep 2

# Shifter
# Shifter must be done manually to use shifter in a box, see run_shifter.sh
