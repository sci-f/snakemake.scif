# Snakemake tutorial workflow with SCIF

**and container friends!**

This repository implements the SnakeMake tutorial workflow and uses the Scientific Filesystem (SCIF) to provide a reproducible research environment. It has been extended to not just include Docker and Singularity, but (as many container technologies as the author can figure out) to demonstrate the variable use for SCIF. To get started, clone this repository and build the containers per the instructions below.

## Installation

 - [Docker](https://docs.docker.com/install/)
 - [Singularity](http://singularity.lbl.gov/install-linux)
 - [Charliecloud](https://hpc.github.io/charliecloud/install.html#manual-build-and-install)


We will take the following format in each section, to show the same command for each
container technology, followed by the output that each results in.

## Building containers
```
sudo singularity build snakemake.simg Singularity
```
```
docker build -t vanessa/snakemake.scif .
```
```
ch-build -t snakemake.ch .
ch-docker2tar snakemake.ch /var/tmp
ch-tar2dir /var/tmp/snakemake.ch.tar.gz /var/tmp
```

The interesting thing about the above is that we could actually build all containers starting with Docker. Like this:

```
docker build -t vanessa/snakemake.scif .
```

Then build on Docker hub, and generate a Singularity container from it.
```
sudo singularity build snakemake.simg docker://vanessa/snakemake.scif
```
```
ch-docker2tar vanessa/snakemake.scif /var/tmp
ch-tar2dir /var/tmp/snakemake.ch.tar.gz /var/tmp
```

This reinforces the idea that the underlying guts of a container are the file contents. Docker delivers the contents nicely in layers, so it's a starting point to dump them into another base, and then interact. Do all roads lead to Docker layers? I'll let you decided.

## Running Containers
For each of the examples below, we will show examples running commands **inside** and **outside** of the containers.

### Run the workflow 

The `Snakefile` specifies rules that should be executed. State a target file and Snakemake will build a directed acyclic graph (DAG) from the rules until the target file is reached. The rules are then executed to create the target file. The Snakefile hides the details of the environment provided by SCIF in the container as you can see from the commands below. All commands can be executed from outside the container or from within the container.

In all cases we bind the example dataset to the data directory provided by SCIF. SCIF provides the environment in which the workflow steps are executed. The workflow steps themself and therefore the interaction with SCIF is implemented in the Snakefile.

**Inside container, Singularity**

```
singularity shell --bind data/:/scif/data snakemake.simg
```


**Inside container, Docker**

To run the whole workflow we need the Snakefile which is not part of the container. To access the Snakefile from within the Docker container, we need an additional mount.

```
docker run -v $PWD/data:/scif/data:z -v $PWD:/working_dir -it -w /working_dir --entrypoint /bin/bash vanessa/snakemake.scif
snakemake all
```

**Inside container, Charliecloud**

```
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /bin/bash
export PATH=/opt/conda/bin:$PATH
snakemake all
```

<hr>

**Outside Container, Singularity**

```
singularity exec --bind data/:/scif/data snakemake.simg snakemake all
```

**Outside Container, Docker**

Because Docker does not mount the current working directory, we mount it to a the SCIF data directory of the SCIF snakemake app. The snakemake app will execute snakemake in `/scif/data/snakemake`, thus snakemake can find the Snakefile.

```
docker run -v $PWD/data:/scif/data:z -v $PWD:/scif/data/snakemake:z -it vanessa/snakemake.scif run snakemake all
```

**Outside Container, Charliecloud**

```
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/snakemake all
#TODO: need to test updated recipe with Docker and Singularity
```

## Generate graphical representation of the workflow

The whole workflow should finish within 5 seconds ... too fast to realize what is happening. A visualization of the whole workflow we just ran would be nice ...

The SCIF app `graphviz_create_dag` generates a directed acyclic graph of the workflow execution plan. This is a feature of Snakemake, but needs Graphviz as additional dependency. Therefore Graphviz is installed in the `appinstall` section of this app.

```
rm -r data/calls/ /data/mapped_reads/ data/sorted_reads/ data/report.html
```

Now we can generate the graph of the execution plan. The plan is computed from the Snakefile in our current working directory and the target file we specify. In this case we want to create the final report file
`report.html`. The directed acyclic graph of the workflow should be saved at`data/dag.svg`. The commands for Singularity and Docker are as follows:

**Singularity**

```
singularity run --bind data/:/scif/data snakemake.simg run graphviz_create_dag $PWD report.html dag.svg
```

**Docker**

In contrast to Singularity, Docker cannot access the current working directory. Therefore we mount the current working directory to the SCIF data folder of the `graphviz_create_dag` app.

```
docker run -v $PWD/data:/scif/data:z -v $PWD:/scif/data/graphviz_create_dag:z -it vanessa/snakemake.scif run graphviz_create_dag /scif/data/graphviz_create_dag report.html dag.svg
```

**Charliecloud**

```
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif run graphviz_create_dag $PWD report.html dag.svg
```

```
> [graphviz_create_dag] executing /bin/bash /scif/apps/graphviz_create_dag/scif/runscript /home/fbartusch/github/snakemake_tutorial report.html dag.svg
> Building DAG of jobs...
```

## Map with bwa mem 
Instead of running the whole workflow at once with Snakemake, we can also use the SCIF environment to execute computations. The following examples also illustrate special SCIF syntax for working with environment variables or piping output to a file. First, we map the reads with bwa mem to the reference genome.

**Inside container, Singularity** 

Note the use of `[e]` so we can easily pass the environment variable to the SCIF. Otherwise, it would be evaluated on the host.

```
singularity shell --bind data/:/scif/data snakemake.simg
mkdir /scif/data/mapped_reads
scif run bwa mem -o [e]SCIF_DATA/mapped_reads/A.sam [e]SCIF_DATA/genome.fa [e]SCIF_DATA/samples/A.fastq
```

**Inside container, Docker**

```
docker run -v $PWD/data:/scif/data:z -it --entrypoint /bin/bash vanessa/snakemake.scif
mkdir -p /scif/data/mapped_reads
scif run bwa mem -o [e]SCIF_DATA/mapped_reads/A.sam [e]SCIF_DATA/genome.fa [e]SCIF_DATA/samples/A.fastq
```

**Inside container, Charliecloud**

```
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif shell
/opt/conda/bin/scif run bwa mem -o [e]SCIF_DATA/mapped_reads/A.sam [e]SCIF_DATA/genome.fa [e]SCIF_DATA/samples/A.fastq
```


<hr>

```
mkdir -p data/mapped_reads
```

**Outside container, Singularity**

```
singularity run --bind data:/scif/data snakemake.simg run bwa mem -o [e]SCIF_DATA/mapped_reads/A.sam [e]SCIF_DATA/genome.fa [e]SCIF_DATA/samples/A.fastq
```

**Outside the container, Docker**

```
docker run -v $PWD/data:/scif/data:z vanessa/snakemake.scif run bwa mem -o [e]SCIF_DATA/mapped_reads/A.sam [e]SCIF_DATA/genome.fa [e]SCIF_DATA/samples/A.fastq
```

**Outside the container, Charliecloud**

```
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif run bwa mem -o [e]SCIF_DATA/mapped_reads/A.sam [e]SCIF_DATA/genome.fa [e]SCIF_DATA/samples/A.fastq
```

## Sam to Bam Conversion

**Inside the container, Singularity**

Note the use of `[out]` as a substitute for `>`. If you wanted to use `>` you could put the entire thing in quotes.

```
singularity shell --bind data/:/scif/data snakemake.simg
scif run samtools view -Sb [e]SCIF_DATA/mapped_reads/A.sam [out] [e]SCIF_DATA/mapped_reads/A.bam    # or
scif run samtools 'view -Sb $SCIF_DATA/mapped_reads/A.sam > $SCIF_DATA/mapped_reads/A.bam'
> [samtools] executing /bin/bash /scif/apps/samtools/scif/runscript view -Sb $SCIF_DATA/mapped_reads/A.sam > $SCIF_DATA/mapped_reads/A.bam
```

**Inside the container, Docker**

```
docker run -v $PWD/data:/scif/data:z -it --entrypoint /bin/bash vanessa/snakemake.scif
scif run samtools view -Sb [e]SCIF_DATA/mapped_reads/A.sam [out] [e]SCIF_DATA/mapped_reads/A.bam    # or
scif run samtools 'view -Sb $SCIF_DATA/mapped_reads/A.sam > $SCIF_DATA/mapped_reads/A.bam'
> [samtools] executing /bin/bash /scif/apps/samtools/scif/runscript view -Sb $SCIF_DATA/mapped_reads/A.sam > $SCIF_DATA/mapped_reads/A.bam
```

**Inside the container, Charliecloud**

```
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif shell
/opt/conda/bin/scif run samtools view -Sb [e]SCIF_DATA/mapped_reads/A.sam [out] [e]SCIF_DATA/mapped_reads/A.bam
```

<hr>

**Outside the container, Singularity**

```
singularity run --bind data:/scif/data snakemake.simg run samtools view -Sb [e]SCIF_DATA/mapped_reads/A.sam [out] [e]SCIF_DATA/mapped_reads/A.bam   # or
singularity run --bind data:/scif/data snakemake.simg run samtools 'view -Sb $SCIF_DATA/mapped_reads/A.sam > $SCIF_DATA/mapped_reads/A.bam'
```

**Outside the container, Docker**

```
docker run -v $PWD/data:/scif/data vanessa/snakemake.scif run samtools view -Sb [e]SCIF_DATA/mapped_reads/A.sam [out] [e]SCIF_DATA/mapped_reads/A.bam   # or
docker run -v $PWD/data:/scif/data vanessa/snakemake.scif run samtools 'view -Sb $SCIF_DATA/mapped_reads/A.sam > $SCIF_DATA/mapped_reads/A.bam'
```

**Outside the container, Charliecloud**

```
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif run samtools view -Sb [e]SCIF_DATA/mapped_reads/A.sam [out] [e]SCIF_DATA/mapped_reads/A.bam
```

## Interactive development
This can be done for Docker or Singularity, just with different commands to shell into the container!

```
docker run -it -v $PWD/data:/scif/data:z vanessa/snakemake.scif pyshell
singularity run --bind data/:/scif/data snakemake.simg pyshell
ch-run /var/tmp/snakemake.ch -- /opt/conda/bin/scif pyshell
```
```
Found configurations for 4 scif apps
bwa
graphviz_create_dag
samtools
snakemake
[scif] /scif bwa | graphviz_create_dag | samtools | snakemake
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
args = ['mem', '-o', '[e]SCIF_DATA/mapped_reads/a.sam', '[e]SCIF_DATA/genome.fa', '[e]SCIF_DATA/samples/A.fastq']
client.run('bwa', args=args)

# Run sam--bam interactively
args = ["view", "-Sb", "/scif/data/mapped_reads/r1_subset.sam", ">", "/scif/data/mapped_reads/r1_subset.bam"]
client.run('samtools', args=args)
```

This repo is also used as example on the SCIF GitHub repo, so take a [look:](https://github.com/sci-f/snakemake.scif)
