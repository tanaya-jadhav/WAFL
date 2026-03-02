srun -J snakemake \
	--export=ALL \
	-t 120:00:00 \
	-n 1 \
	--mem-per-cpu=100G \
	-D $PWD \
	-o logs/snakemake.out \
	-e logs/snakemake.err \
	$PWD/snakemake.script.sh &
