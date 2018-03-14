#!/bin/bash

# here we will create a shifter image and calculate metrics using the Shifter
# Box! This would usually be installed on a cluster.

docker run -it --rm --privileged scanon/shifterbox:2

# We need git to clone the repo!
yum install -y git

/src/start.sh
su - auser

shifterimg pull vanessa/snakemake.scif

# If you get an error about contacting the gateway, it's a race condition likely, 
# and you should exit the user, as the superuser remove the lock for MongoDB,
# and try again:
# exit
# rm -rf /data/db/
# /src/start.sh

git clone https://github.com/sci-f/snakemake.scif
cd snakemake.scif/

shifter --image=vanessa/snakemake.scif --volume `pwd`/data:/scif/data --workdir /scif/data /opt/conda/bin/scif run valgrind
shifter --image=vanessa/snakemake.scif --volume `pwd`/data:/scif/data --workdir /scif/data /opt/conda/bin/scif run timeit

mv data/massif.* $RESULT/charliecloud
mv data/snakemake-times.log $RESULT/charliecloud
