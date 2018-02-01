run: clean build

clean:
	rm -f snakemake

build: clean
	sudo singularity build snakemake Singularity
