#!/bin/bash
# Loop to submit one job per chromosome (1 to 22)
for c in {1..22}
do
  sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=combine_chr$c
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=2:00:00
#SBATCH --account=jvmorr1
#SBATCH --partition=standard
#SBATCH --output=out_chr$c.log
#SBATCH --error=err_chr$c.log

module load R/4.2.0

Rscript combine_gwas.R $c
EOT
done
