srun -J snakemake \
	--export=ALL \
	-p defq \
	-t 120:00:00 \
	-n 1 \
	--mem-per-cpu=100G \
	-D $PWD \
	-o logs/snakemake_SV.out \
	-e logs/snakemake_SV.err \
	$PWD/snakemake.script_SV.sh &