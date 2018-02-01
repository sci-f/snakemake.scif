# Singularity

## build container

```
sudo singularity build snakemake Singularity
```

or use the Makefile, if you have `build-essential` installed:

```
make
```

This builds a squashfs filesystem, meaning that it is read only. We will write data by way of mapping a directory to the host.

## Map with bwa mem 

### Inside container

```
# shell into the container
singularity shell --bind data/:/scif/data snakemake
$ mkdir -p /scif/data/mapped_reads
```
```
# environment variables for scif are referenced with [e]
$ scif run bwa mem -o [e]SCIF_DATA/mapped_reads/r1_subset.sam [e]SCIF_DATA/index/ref.fa [e]SCIF_DATA/raw_reads/r1_subset.fq
[bwa] executing /scif/apps/bwa/bin/bwa mem -o $SCIF_DATA/mapped_reads/r1_subset.sam $SCIF_DATA/index/ref.fa $SCIF_DATA/raw_reads/r1_subset.fq
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 10000 sequences (2064809 bp)...
[M::mem_process_seqs] Processed 10000 reads in 0.679 CPU sec, 0.679 real sec
[main] Version: 0.7.17-r1188
[main] CMD: /scif/apps/bwa/bin/bwa mem -o /scif/data/mapped_reads/r1_subset.sam /scif/data/index/ref.fa /scif/data/raw_reads/r1_subset.fq
[main] Real time: 0.728 sec; CPU: 0.700 sec
```

### Outside the container

```
mkdir -p data/mapped_reads
```
```
singularity run --bind data:/scif/data snakemake run bwa mem -o [e]SCIF_DATA/mapped_reads/r1_subset.sam [e]SCIF_DATA/index/ref.fa [e]SCIF_DATA/raw_reads/r1_subset.fq
```

## Sam -> Bam

### Inside the container

```
singularity shell --bind data/:/scif/data snakemake

# Note the use of [out] as a substitute for >. If you wanted to use > you could put the entire thing in quotes.
scif run samtools view -Sb /scif/data/mapped_reads/r1_subset.sam [out] /scif/data/mapped_reads/r1_subset.bam
[samtools] executing /bin/bash /scif/apps/samtools/scif/runscript view -Sb /scif/data/mapped_reads/r1_subset.sam > /scif/data/mapped_reads/r1_subset.bam
```

### Outside the container
```
singularity run --bind data:/scif/data snakemake run samtools view -Sb [e]SCIF_DATA/mapped_reads/r1_subset.sam [out] [e]SCIF_DATA/mapped_reads/r1_subset.bam
[samtools] executing /bin/bash /scif/apps/samtools/scif/runscript view -Sb $SCIF_DATA/mapped_reads/r1_subset.sam > $SCIF_DATA/mapped_reads/r1_subset.bam
```
Notice that we are using `[e]` instead of the traditional `$` for an environment variable, and `[out]` instead of a traditional pipe so the entire thing is passed into the container with SCIF to run. Otherwise, it would be evaluated on the host.

# Docker

## build container

```
docker build -t vanessa/snakemake.scif .
```

## Map with bwa mem 
For these examples, we will bind data to the host as a volume at `/data`.

### Inside container

```
# shell into the container
docker run -it --entrypoint /bin/bash vanessa/snakemake.scif
mkdir -p /scif/data/mapped_reads

# if you want to map data to the host, you need a volume
docker run -v $PWD/data:/scif/data -it --entrypoint /bin/bash vanessa/snakemake.scif
```
```
# environment variables for scif are referenced with [e]
$ scif run bwa mem -o [e]SCIF_DATA/mapped_reads/r1_subset.sam [e]SCIF_DATA/index/ref.fa [e]SCIF_DATA/raw_reads/r1_subset.fq
[bwa] executing /scif/apps/bwa/bin/bwa mem -o $SCIF_DATA/mapped_reads/r1_subset.sam $SCIF_DATA/index/ref.fa $SCIF_DATA/raw_reads/r1_subset.fq
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 10000 sequences (2064809 bp)...
[M::mem_process_seqs] Processed 10000 reads in 0.679 CPU sec, 0.679 real sec
[main] Version: 0.7.17-r1188
[main] CMD: /scif/apps/bwa/bin/bwa mem -o /scif/data/mapped_reads/r1_subset.sam /scif/data/index/ref.fa /scif/data/raw_reads/r1_subset.fq
[main] Real time: 0.728 sec; CPU: 0.700 sec
```

### Outside the container

```
mkdir -p data/mapped_reads
```
```
docker run -v $PWD/data:/scif/data vanessa/snakemake.scif run bwa mem -o [e]SCIF_DATA/mapped_reads/r1_subset.sam [e]SCIF_DATA/index/ref.fa [e]SCIF_DATA/raw_reads/r1_subset.fq
```

## Sam -> Bam

### Inside the container

```
docker run -v $PWD/data:/scif/data -it --entrypoint /bin/bash vanessa/snakemake.scif

# Note the use of [out] as a substitute for >
scif run samtools view -Sb /scif/data/mapped_reads/r1_subset.sam [out] /scif/data/mapped_reads/r1_subset.bam
[samtools] executing /bin/bash /scif/apps/samtools/scif/runscript view -Sb /scif/data/mapped_reads/r1_subset.sam > /scif/data/mapped_reads/r1_subset.bam
```

You will notice in the above that the `[out]` is parsed internally as the correct `>`.

### Outside the container
```
docker run -v $PWD/data:/scif/data vanessa/snakemake.scif run samtools view -Sb /scif/data/mapped_reads/r1_subset.sam [out] /scif/data/mapped_reads/r1_subset.bam

singularity run --bind data:/scif/data snakemake run samtools "view -Sb [e]SCIF_DATA/mapped_reads/r1_subset.sam > [e]SCIF_DATA/mapped_reads/r1_subset.bam"
[samtools] executing /bin/bash /scif/apps/samtools/scif/runscript view -Sb $SCIF_DATA/mapped_reads/r1_subset.sam > $SCIF_DATA/mapped_reads/r1_subset.bam
```
Notice that we are using `[e]` instead of the traditional `$` for an environment variable, and the entire command is in quotes so it gets passed into the container. This is to ensure that the environment variable is not evaluated on the host.


## Interactive development
This can be done for Docker or Singularity, just with different commands to shell into the container!

```
docker run -it -v $PWD/data:/scif/data vanessa/snakemake.scif pyshell
singularity run --bind data/:/scif/data snakemake pyshell

# If you do need to write data (and make the bind)
Found configurations for 3 scif apps
bwa
samtools
snakemake
[scif] /scif bwa | samtools | snakemake
Python 3.6.2 |Anaconda, Inc.| (default, Sep 22 2017, 02:03:08) 
[GCC 7.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
(InteractiveConsole)
>>> 

# Activate bwa 
client.activate('samtools')

# Environment variables active
client.environment

# Run bwa interactively
args = ['mem', '-o', '[e]SCIF_DATA/mapped_reads/r1_subset.sam', '[e]SCIF_DATA/index/ref.fa', '[e]SCIF_DATA/raw_reads/r1_subset.fq']
client.run('bwa', args=args)

# Run sam--bam interactively
args = ["view", "-Sb", "/scif/data/mapped_reads/r1_subset.sam", ">", "/scif/data/mapped_reads/r1_subset.bam"]
client.run('samtools', args=args)
```
