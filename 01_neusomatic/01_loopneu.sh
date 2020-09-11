#!/bin/bash

#mkdir -p ./sbatch/tmp
#mkdir -p ./sbatch/out
#mkdir -p ./sbatch/err

cwd=$(pwd)
ref="hg38.fa"
ref_dir="/scratch/path"
bam_dir="/scratch/path/aligns"
out_dir="/scratch/path/output/neusomatic"
samplepairsdir="/home/path"

NEUSOMATIC_BIN="/opt/neusomatic/neusomatic/bin"
NEUSOMATIC_MODELS="/opt/neusomatic/neusomatic/models"
NEUSOMATIC_SCAN_ALIGNMENTS="/opt/neusomatic/neusomatic/bin/scan_alignments"

userProject="projectName"

{
read
while IFS=, read -r normal tumour rest
do
	echo $normal
	echo $tumour
	sample="${normal}_${tumour}"
	echo $sample
	pretmp=${cwd}/sbatch/tmp/preneu_${sample}.sh
        rm -rf ${pretmp} || true
        touch ${pretmp}

        #set variables used by slurm management
        echo "#!/bin/bash" >> ${pretmp} 
        echo "#SBATCH -p compute" >> ${pretmp}
        #echo "#SBATCH --mem=${mem}" >> ${pretmp}
        #echo "#SBATCH --nodes=${nodes}" >> ${pretmp}
        echo "#SBATCH --ntasks=1" >> ${pretmp}
	echo "#SBATCH --cpus-per-task=20" >> ${pretmp}
        #echo "#SBATCH --tasks-per-node=${nodes}" >> ${pretmp}
        #echo "#SBATCH -t ${runTime}" >> ${pretmp}
        echo "#SBATCH -o ${cwd}/sbatch/out/preneu${sample}.%J" >> ${pretmp}
        echo "#SBATCH -e ${cwd}/sbatch/err/preneu${sample}.%J" >> ${pretmp}
        echo "#SBATCH --job-name=preneu" >> ${pretmp}
        echo "#SBATCH --account=${userProject}" >> ${pretmp}
	echo "#SBATCH --mail-user=email@address.ac.uk" >> ${pretmp}
	echo "#SBATCH --mail-type=ALL" >> ${pretmp}

	echo "module load singularity" >> ${pretmp}
	echo "cd /scratch/path/" >> ${pretmp}
	echo "singularity exec /home/path/dockertools/neusomatic_latest.sif preprocess.py --mode call --reference hg38.fa --region_bed hg38.bed --tumor_bam aligns/${tumour}.bam --normal_bam aligns/${normal}.bam --work work_${sample} --min_mapq 10 --num_threads 20 --scan_alignments_binary /opt/neusomatic/neusomatic/bin/scan_alignments" >> ${pretmp}
	chmod u+x ${pretmp}

	#jid1=$(sbatch ${pretmp})

	tmp=${cwd}/sbatch/tmp/neu_${sample}.sh
        rm -rf ${tmp} || true
        touch ${tmp}

        #set variables used by slurm management
        echo "#!/bin/bash" >> ${tmp} 
        echo "#SBATCH -p compute" >> ${tmp}
        #echo "#SBATCH --mem=${mem}" >> ${tmp}
        #echo "#SBATCH --nodes=${nodes}" >> ${tmp}
        echo "#SBATCH --ntasks=1" >> ${tmp}
	echo "#SBATCH --cpus-per-task=20" >> ${tmp}
        #echo "#SBATCH --tasks-per-node=${nodes}" >> ${tmp}
        #echo "#SBATCH -t ${runTime}" >> ${tmp}
        echo "#SBATCH -o ${cwd}/sbatch/out/neu${sample}.%J" >> ${tmp}
        echo "#SBATCH -e ${cwd}/sbatch/err/neu${sample}.%J" >> ${tmp}
        echo "#SBATCH --job-name=neu" >> ${tmp}
        echo "#SBATCH --account=${userProject}" >> ${tmp}
	echo "#SBATCH --mail-user=email@address.ac.uk" >> ${tmp}
	echo "#SBATCH --mail-type=ALL" >> ${tmp}

	echo "module load singularity" >> ${tmp}
	echo "cd /scratch/path/" >> ${tmp}
	echo "singularity exec /home/path/dockertools/neusomatic_latest.sif python /opt/neusomatic/neusomatic/python/call.py --candidates_tsv work_${sample}/dataset/*/candidates*.tsv --reference hg38.fa --out work_${sample} --checkpoint /opt/neusomatic/neusomatic/models/NeuSomatic_v0.1.4_standalone_SEQC-WGS-Spike.pth --num_threads 20 --batch_size 100" >> ${tmp}
	chmod u+x ${tmp}

	##jid2=$(sbatch --dependency=afterok:${jid1} ${tmp})

	posttmp=${cwd}/sbatch/tmp/postneu_${sample}.sh
        rm -rf ${posttmp} || true
        touch ${posttmp}

        #set variables used by slurm management
        echo "#!/bin/bash" >> ${posttmp} 
        echo "#SBATCH -p compute" >> ${posttmp}
        #echo "#SBATCH --mem=${mem}" >> ${posttmp}
        #echo "#SBATCH --nodes=${nodes}" >> ${posttmp}
        echo "#SBATCH --mem-per-cpu=100G" >> ${posttmp}
	#echo "#SBATCH --ntasks=1" >> ${posttmp}
	echo "#SBATCH --cpu-per-task=10" >> ${posttmp}
        #echo "#SBATCH --tasks-per-node=${nodes}" >> ${posttmp}
        #echo "#SBATCH -t ${runTime}" >> ${posttmp}
        echo "#SBATCH -o ${cwd}/sbatch/out/postneu${sample}.%J" >> ${posttmp}
        echo "#SBATCH -e ${cwd}/sbatch/err/postneu${sample}.%J" >> ${posttmp}
        echo "#SBATCH --job-name=postneu" >> ${posttmp}
        echo "#SBATCH --account=${userProject}" >> ${posttmp}
	echo "#SBATCH --mail-user=email@address.ac.uk" >> ${posttmp}
	echo "#SBATCH --mail-type=ALL" >> ${posttmp}

	echo "module load singularity" >> ${posttmp}
	echo "cd /scratch/path/" >> ${posttmp}
	echo "singularity exec /home/path/dockertools/neusomatic_latest.sif python /opt/neusomatic/neusomatic/python/postprocess.py --reference hg38.fa --tumor_bam aligns/${tumour}.bam --pred_vcf work_${sample}/pred.vcf --candidates_vcf work_${sample}/work_tumor/filtered_candidates.vcf --output_vcf work_${sample}/NeuSomatic.vcf --work work_${sample}"  >> ${posttmp}
	echo "mkdir -p output/neusomatic/${sample}" >> ${posttmp}
	echo "mv work_${sample}/*.vcf output/neusomatic/${sample}" >> ${posttmp}
	echo "mv work_${sample}/*.bed output/neusomatic/${sample}" >> ${posttmp}
	echo "rm -r work_${sample}" >> ${posttmp}

	chmod u+x ${posttmp}

	#jid3=$(sbatch --dependency=afterok:${jid2} ${posttmp})

done
} < ${samplepairsdir}/sample_pairs.csv
