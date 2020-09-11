#!/bin/bash
#SBATCH -p compute
#SBATCH --mem-per-cpu=100G
#SBATCH -o /home/path/neu-bin/sbatch/out/postneuDB143_DB144.%J
#SBATCH -e /home/path/neu-bin/sbatch/err/postneuDB143_DB144.%J
#SBATCH --job-name=postneu
#SBATCH --account=projectAccount
#SBATCH --mail-user=email@address.ac.uk
#SBATCH --mail-type=ALL
module load singularity
cd /scratch/path/
singularity exec /home/path/dockertools/neusomatic_latest.sif python /opt/neusomatic/neusomatic/python/postprocess.py --reference hg38.fa --tumor_bam aligns/DB144.bam --pred_vcf work_DB143_DB144/pred.vcf --candidates_vcf work_DB143_DB144/work_tumor/filtered_candidates.vcf --output_vcf work_DB143_DB144/NeuSomatic.vcf --work work_DB143_DB144
mkdir -p output/neusomatic/DB143_DB144
mv work_DB143_DB144/*.vcf output/neusomatic/DB143_DB144
mv work_DB143_DB144/*.bed output/neusomatic/DB143_DB144
rm -r work_DB143_DB144
