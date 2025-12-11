#!/bin/sh
snakemake --configfile config.yaml --use-singularity --jobs 500 --singularity-args "--bind /mnt/isilon" --printshellcmds --latency-wait 20 --keep-going
