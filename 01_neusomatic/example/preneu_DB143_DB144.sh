#!/bin/bash
#SBATCH -p compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -o /home/path/neu-bin/sbatch/out/preneuDB143_DB144.%J
#SBATCH -e /home/path/neu-bin/sbatch/err/preneuDB143_DB144.%J
#SBATCH --job-name=preneu
#SBATCH --account=projectAccount
#SBATCH --mail-user=email@address.ac.uk
#SBATCH --mail-type=ALL
module load singularity
cd /scratch/path/
singularity exec /home/path/dockertools/neusomatic_latest.sif preprocess.py --mode call --reference hg38.fa --region_bed hg38.bed --tumor_bam aligns/DB144.bam --normal_bam aligns/DB143.bam --work work_DB143_DB144 --min_mapq 10 --num_threads 20 --scan_alignments_binary /opt/neusomatic/neusomatic/bin/scan_alignments
