run: clean build

logfile := file_$(shell date +%FT%T%Z).log

clean:
	rm -f snakemake.simg

build: clean
	sudo singularity build snakemake.simg Singularity 2>&1 | tee build.log
