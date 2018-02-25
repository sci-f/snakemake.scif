# Snakemake tutorial workflow with SCIF

**and container friends!**

This repository implements the SnakeMake tutorial workflow and uses the Scientific Filesystem (SCIF) to provide a reproducible research environment. It has been extended to not just include Docker and Singularity, but (as many container technologies as @vsoch can figure out) to demonstrate the variable use for SCIF. To get started, clone this repository and build the containers per the instructions below.

## Installation

 - [Docker](https://docs.docker.com/install/)
 - [Singularity](http://singularity.lbl.gov/install-linux)
 - [Charliecloud](https://hpc.github.io/charliecloud/install.html#manual-build-and-install)
 - [runc](https://github.com/opencontainers/runc)

The following container technologies don't have user friendly build clients available (but are used by admins successfully in their respective locations!).

 - [Inception](https://github.com/NCAR/Inception): I couldn't figure out how to use it, my [issue is posted here](https://github.com/NCAR/Inception/issues/9).

and the following containers had various bugs, or I just couldn't figure it out from the documentation available.

 - [shifter](https://github.com/NERSC/Shifter-Tutorial/tree/master/examples/shifter-in-a-box)

We will take the following format in each section, to show the same command for each
container technology, followed by the output that each results in.

## Building containers

**Singularity**

```
sudo singularity build snakemake.simg Singularity
```

**Docker**
```
docker build -t vanessa/snakemake.scif .
```

**Charliecloud**

```
ch-build -t snakemake.ch .
ch-docker2tar snakemake.ch /var/tmp
ch-tar2dir /var/tmp/snakemake.ch.tar.gz /var/tmp
```

**runc**

```
mkdir -p /tmp/runc/rootfs
sudo singularity build --sandbox /tmp/runc/rootfs snakemake

# I'm not sure if this is needed, but I did it.
sudo chown -R vanessa /tmp/runc/rootfs
cd /tmp/runc
runc spec --rootless --bundle /tmp/runc
ls /tmp/runc
config.json  rootfs
```

At this point we need to customize the config.json to allow for specific mounts, 
and add the `/opt/conda/bin` to the `$PATH` to find the scif executable.

```
"env": ["PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
        "TERM=xterm"
       ],
...

	{
		"destination": "/scif/data",
		"type": "bind",
                "options": ["bind"],
		"source": "/home/vanessa/Documents/Dropbox/Code/scif/examples/snakemake.scif/data"
	},
```


**shifter**

The user should know this is the "shifter in a box" which is likely much harder to use than a cluster
installation. I don't have shifter at my insitution so this was the approach I chose to test.

```
docker run -it --rm --privileged -v $PWD/data:/scratch -v $PWD:/code scanon/shifterbox:2
```

This is a bit confusing, but we will need to again bind the folders with the snakefile ($PWD on the host)
somewhere inside the container, and the /scratch data directory to /scif/data.

```
chmod 1777 /scratch
chmod 1777 /code
```

Next, edit the config to change

```
siteFs=/home:/home
# to
siteFs=/home:/home;/scratch:/scif/data;/code:/code
```

The above is a bit weird - it says that we need to bind the data folder to /scif/data
in the container (we can't bind /scif because we lose our apps) and we want to bind
/code to /code, as this is the folder bound from the host with the Snakefile 
for the workflow.

You will need an editor and then to open the right file!

```
yum install -y vim
vim /etc/shifter/udiRoot.conf
```
and then start! Unfortunately this only works for me about half the time, and usually on the second go. But you can try again.

```
/src/start.sh
su - auser
shifterimg pull vanessa/snakemake.scif:container-friends
```
If you get an error about contacting the gateway, it's a race condition likely, and you should remove the lock for MongoDB and try again:

```
$ rm -rf /data/db/
$ /src/start.sh
```

If you intend to run / test the shifter box, leave it open! You need the pulled container.

The interesting thing about the above is that we are building many of the containers by way of starting with Docker layers. This reinforces the idea that the underlying guts of a container are the file contents. Docker delivers the contents nicely in layers, so it's a starting point to dump them into another base, and then interact. Do all roads lead to Docker layers? I'll let you decided.

## Running Containers
For each of the examples below, we will show examples running commands **inside** and **outside** of the containers.

### Run the workflow 

The `Snakefile` specifies rules that should be executed. State a target file and Snakemake will build a directed acyclic graph (DAG) from the rules until the target file is reached. The rules are then executed to create the target file. The Snakefile hides the details of the environment provided by SCIF in the container as you can see from the commands below. All commands can be executed from outside the container or from within the container.

In all cases we bind the example dataset to the data directory provided by SCIF. SCIF provides the environment in which the workflow steps are executed. The workflow steps themself and therefore the interaction with SCIF is implemented in the Snakefile.

**Inside container, Singularity**

```
singularity shell --bind data/:/scif/data snakemake.simg
$ scif run snakemake all
```


**Inside container, Docker**

To run the whole workflow we need the Snakefile which is not part of the container. To access the Snakefile from within the Docker container, we need an additional mount.

```
docker run -v $PWD/data:/scif/data:z -v $PWD:/working_dir -it -w /working_dir --entrypoint /bin/bash vanessa/snakemake.scif
$ scif run snakemake all
```

**Inside container, Charliecloud**

```
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /bin/bash
export PATH=/opt/conda/bin:$PATH
$ scif run snakemake all
```

**Inside container, runc**

The complete config file can be seen in [configs/runc-config-interactive.json].
This interactive version uses a shell entrypoint, and then runs scif. If the Snakemake and config.yaml
aren't found in the data directory, they are copied there from the container snakemake application.

```
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
# which scif
/opt/conda/bin/scif
# scif run snakemake all
[snakemake] executing /bin/bash /scif/apps/snakemake/scif/runscript all
...
Complete log: /scif/data/.snakemake/log/2018-02-25T180544.725734.snakemake.log
```

**Inside container, shifter**

Note that we are working in a container that has had the build step done (as shown above to pull)

```
$ shifterimg images
mycluster  docker     READY    e4e467178d   2018-02-23T21:29:58 vanessa/snakemake.scif:container-friends

# I would want to do this, because I have /scratch in the container. 
# The subfolder "data" of scratch should be mounted to /scif/data, and /scratch should be bound anywhere really.
# The base folder /scif cannot be bound as
$ shifter --image=vanessa/snakemake.scif:container-friends --volume /scratch:/scif/data /bin/bash
$ cd /
```

This doesn't work
```
shifter --image=vanessa/snakemake.scif:container-friends  /bin/bash
FAILED to create volume "to": /var/udiMount//scif/data, cannot create mount points in that location
FAILED to mount siteFs volumes
FAILED to properly setup site modifications
FAILED to mount image into UDI
FAILED to setup image.
```
```
# but I figured out pwd command...
shifter --image=vanessa/snakemake.scif:container-friends --volume /scratch:/scif/data --workdir /scratch /bin/bash
```

<hr>

**Outside Container, Singularity**

```
singularity run --bind data/:/scif/data snakemake.simg run snakemake all
```

**Outside Container, Docker**

Because Docker does not mount the current working directory, we mount it to a the SCIF data directory of the SCIF snakemake app. The snakemake app will execute snakemake in `/scif/data/snakemake`, thus snakemake can find the Snakefile.

```
docker run -v $PWD/data:/scif/data:z -v $PWD:/scif/data/snakemake:z -it vanessa/snakemake.scif run snakemake all
```

**Outside Container, Charliecloud**

```
ch-run --cd $PWD -b data:/scif/data /var/tmp/snakemake.ch -- /opt/conda/bin/scif run snakemake all
```

**Outside container, runc**

This is a very easy switch to do! All that is needed is to define the entrypoint to the
container via the configuration file to be the scif entrypoint, and we can even "hard code" the application
to run as snakemake. This modified config file can be found in [configs/runc-config-noninteractive.json].

```
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
[snakemake] executing /bin/bash /scif/apps/snakemake/scif/runscript all
...
Complete log: /scif/data/.snakemake/log/2018-02-25T180544.725734.snakemake.log
```

So it comes down to changing:

```
"args": [
	"sh"
],
```

to

```
"args": [
	"scif", "run", "snakemake", "all"
],
```

I can imagine that I would have some strategy to either manage a small 
collection of config files to use per rootfs, **or** a method to produce
the one I need on demand. I suspect this is what many popular container 
technologies are already doing and I just haven't realized it :)


**Outside container, shifter**

Note that we are working in a container that has had the build step done (as shown above to pull)

```
$ shifterimg images
mycluster  docker     READY    e4e467178d   2018-02-23T21:29:58 vanessa/snakemake.scif:container-friends
$ shifter --image=vanessa/snakemake.scif:container-friends --volume /scratch:/scif/data /opt/conda/bin/snakemake all
# STOPPED here - need $PWD or equivalent command.
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

**runc**
This is the same deal as before, we are just going to use a config to specify a different entry point and arguments to generate the dag. Since you've done this twice, I'll just show you the modified command:

```
"args": [
	"scif", "run", "graphviz_create_dag", "/scif/data", "report.html","dag.svg"
	],
```

and then run! Note that if you do this without changing permission, you would get an error:

```
process:10): Pango-WARNING **: error opening config file '/root/.config/pango/pangorc': Permission denied
```

and this would need to be resolved by specifying a different file, or changing the permission of the
file from the initial build.

```
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc[graphviz_create_dag] executing /bin/bash /scif/apps/graphviz_create_dag/scif/runscript /scif/data report.html dag.svg
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

**Inside container, runc**
Here we are shelling inside, so we are using the first config.json with args as a list with "sh".

```
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
scif run bwa mem -o [e]SCIF_DATA/mapped_reads/A.sam [e]SCIF_DATA/genome.fa [e]SCIF_DATA/samples/A.fastq
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

**Outside container, runc**

The modified "args" in the config.json is:

```
"args": [
	"scif", "run", "bwa", "mem", "-o", "[e]SCIF_DATA/mapped_reads/A.sam", "[e]SCIF_DATA/genome.fa", "[e]SCIF_DATA/samples/A.fastq"
	],
```

and then run! 

```
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
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

**Inside container, runc**
Here we are again shelling inside, so we are using the first config.json with args as a list with "sh".

```
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
scif run samtools view -Sb [e]SCIF_DATA/mapped_reads/A.sam [out] [e]SCIF_DATA/mapped_reads/A.bam
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

**Outside container, runc**
Change the config.json args to:

```
"args": [
	"scif", "run", "samtools", "view", "-Sb", "[e]SCIF_DATA/mapped_reads/A.sam", "[out]", "[e]SCIF_DATA/mapped_reads/A.bam"
	],
```

and then run
```
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
```

## Interactive development
This can be done for Docker, Singularity, Charliecloud, or runc, just with different commands to shell into the container. Note that for runc, the args in the config.json needs to be:

```
"args": [
	"scif", "pyshell"
	],
```
```
docker run -it -v $PWD/data:/scif/data:z vanessa/snakemake.scif pyshell
singularity run --bind data/:/scif/data snakemake.simg pyshell
ch-run /var/tmp/snakemake.ch -- /opt/conda/bin/scif pyshell
runc --root /tmp/runc run --bundle /tmp/runc snakmake.runc
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
