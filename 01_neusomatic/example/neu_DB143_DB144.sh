#!/bin/bash
#SBATCH -p compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -o /home/path/neu-bin/sbatch/out/neuDB143_DB144.%J
#SBATCH -e /home/path/neu-bin/sbatch/err/neuDB143_DB144.%J
#SBATCH --job-name=neu
#SBATCH --account=projectAccount
#SBATCH --mail-user=email@address.ac.uk
#SBATCH --mail-type=ALL
module load singularity
cd /scratch/path/
singularity exec /home/path/dockertools/neusomatic_latest.sif python /opt/neusomatic/neusomatic/python/call.py --candidates_tsv work_DB143_DB144/dataset/*/candidates*.tsv --reference hg38.fa --out work_DB143_DB144 --checkpoint /opt/neusomatic/neusomatic/models/NeuSomatic_v0.1.4_standalone_SEQC-WGS-Spike.pth --num_threads 20 --batch_size 100
