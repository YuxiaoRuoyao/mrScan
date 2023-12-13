#!/bin/bash

mkdir -p log
mkdir -p results
snakemake \
   --keep-going \
   --jobs 96 \
   --max-jobs-per-second 5 \
   --latency-wait 30 \
   --cluster-config cluster.yaml  \
   --cluster "sbatch \
              --output={cluster.log}_%j.out \
              --error={cluster.log}_%j.err \
              --account=yuxiaow \ # You should change to your own account
              --job-name={cluster.name} \
              --time={cluster.time}  \
              --cpus-per-task={cluster.cpus}  \
              --mem={cluster.mem}"

