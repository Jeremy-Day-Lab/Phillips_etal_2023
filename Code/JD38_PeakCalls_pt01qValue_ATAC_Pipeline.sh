#!/bin/bash
#
#SBATCH --job-name=JD38_MACS2_rn7
#SBATCH --output=JD38_MACS2_rn7.out
#SBATCH --error=JD38_MACS2_rn7.err
#
# Number of tasks needed for this job. Generally, used with MPI jobs
#SBATCH --ntasks=1
#SBATCH --partition=medium
#
#
# Number of CPUs allocated to each task.
#SBATCH --cpus-per-task=10
#
# Mimimum memory required per allocated  CPU  in  MegaBytes.
#SBATCH --mem-per-cpu=5000
#
# Send mail to the email address when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rphill3@uab.edu

#Set your environment here

#Run your commands here

#load MACS2
module load MACS2

macs2 --version


for x in `ls /data/project/daylab/2019-JD-0038/Rn7/NBI_ATAC_seq_PE_pipeline/BAM_outputs/*dedub_qc_filtered.bam`;
do
	y="$(basename "$x" | cut -d "." -f1)"
	macs2 callpeak \
                --treatment $x \
                --qvalue .01 \
                --gsize 2626580772 \
                --format BAMPE \
                --outdir /data/project/daylab/2019-JD-0038/Rn7/NBI_ATAC_seq_PE_pipeline/PeakCalling \
                --name "$y"_output
done
