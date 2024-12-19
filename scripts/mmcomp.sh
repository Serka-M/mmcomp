#!/usr/bin/bash -l
#SBATCH --job-name=mmcomp
#SBATCH --output=job_%j_%x.out
#SBATCH --error=job_%j_%x.err
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@email.com

# Activate conda environment with only snakemake
conda activate snakemake_8.14.0

# Set default paths
export TMPDIR=/path/to/tmp
export SLURM_TMPDIR=/path/to/tmp

# Main workflow
snakemake --profile custom --executor slurm --default-resources slurm_partition={resources.slurm_partition} mem_mb={resources.mem_mb} runtime={resources.runtime} --software-deployment-method conda --cores 250 --jobs 20 --retries 4 -s mmcomp.smk --configfile config.yaml --rerun-incomplete --nolock
