#!/bin/bash

# This script assumes that we start in the same directory as this script,
# and will run and save results for different pipelines (to compare later)

# sudo apt-get install -y valgrind
# usage source measure_comparison.sh

HERE=$PWD
DATA=$HERE/data
RESULT=$HERE/results

mkdir -p $RESULT

function clean () {
    echo
    echo "Cleaning up previous results..."
    rm -rf $HERE/data/sorted_reads
    rm -rf $HERE/data/mapped_reads
    rm -rf $HERE/data/calls
    rm -rf $HERE/data/genome.fa.*
    rm -rf $HERE/data/report.html
}

function doBuild() {
   if [ ${DOBUILD} -eq "1" ]
       then 
       $1
   fi 
}

# Docker
clean
doBuild "docker build --no-cache -t vanessa/snakemake.scif ."
docker run -v $PWD/data:/scif/data -it vanessa/snakemake.scif run valgrind
docker run -v $PWD/data:/scif/data -it vanessa/snakemake.scif run timeit
mv data/massif.* $RESULT/docker
mv data/snakemake-times.log $RESULT/docker
sudo chown -R vanessa data/

# Singularity
clean
doBuild "sudo singularity build snakemake.simg Singularity"
singularity run --bind data/:/scif/data snakemake.simg run valgrind
singularity run --bind data/:/scif/data snakemake.simg run timeit
mv data/massif.* $RESULT/singularity
mv data/snakemake-times.log $RESULT/singularity

# Charliecloud
clean
doBuild "ch-build -t snakemake.ch ." 
doBuild "ch-docker2tar snakemake.ch /var/tmp"
doBuild "ch-tar2dir /var/tmp/snakemake.ch.tar.gz /var/tmp"
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif run valgrind
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif run timeit
mv data/massif.* $RESULT/charliecloud
mv data/snakemake-times.log $RESULT/charliecloud

# Runc
clean
mkdir -p /tmp/runc/rootfs
mkdir -p /tmp/runc/tmp
doBuild "sudo singularity build --sandbox /tmp/runc/rootfs snakemake.simg"
cp configs/runc-config-valgrind.json /tmp/runc/config.json
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
mv data/massif.* $RESULT/runc
cp configs/runc-config-timeit.json /tmp/runc/config.json
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
mv data/snakemake-times.log $RESULT/runc

# Shifter
clean
#echo "Shifter pipeline will be run on cluster!"
#echo "shifter --image=vanessa/snakemake.scif --volume `pwd`/data:/scif/data --workdir /scif/data /opt/conda/bin/scif run snakemake all"
