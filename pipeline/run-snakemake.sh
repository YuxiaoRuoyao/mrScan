#!/bin/bash

mkdir -p log
mkdir -p results
mkdir -p data
snakemake \
   --keep-going \
   --jobs 96 \
   --max-jobs-per-second 5 \
   --latency-wait 60 \
   --cluster-config cluster.yaml  \
   --cluster "sbatch \
              --output={cluster.log}_%j.out \
              --error={cluster.log}_%j.err \
              --account=jvmorr1 \
              --job-name={cluster.name} \
              --time={cluster.time}  \
              --cpus-per-task={cluster.cpus}  \
              --mem={cluster.mem}"



