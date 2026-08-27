#!/bin/bash
#SBATCH --job-name=cellranger_count        ## Job name
#SBATCH -A vswarup_lab                     ## Account name
#SBATCH -p standard                        ## Partition/queue name
#SBATCH --nodes=1                          ## Number of nodes
#SBATCH --ntasks=1                         ## Number of tasks
#SBATCH --cpus-per-task=8                  ## Number of CPU cores per task
#SBATCH --mem=228G                         ## Total memory for job
#SBATCH --error=err_%j.log                 ## Error log file (%j is the job ID)
#SBATCH --output=out_%j.log                ## Output log file (%j is the job ID)
#SBATCH --mail-user=msaddala@uci.edu       ## Email for notifications
#SBATCH --mail-type=END,FAIL               ## Email on job completion or failure
#SBATCH --time=96:00:00                    ## Time limit (4 days)

# Load Cell Ranger module (Modify this based on your HPC environment)
module load cellranger/8.0.1  # Adjust the version if necessary

# Define input/output directories
DATA_DIR="/dfs7/swaruplab/msaddala/test/ganji/TCR_seq/data"
OUTPUT_DIR="/dfs7/swaruplab/msaddala/test/ganji/TCR_seq/cellranget_out_transcriptome"
GENOME_REFERENCE="/dfs7/swaruplab/msaddala/test/ganji/TCR_seq/mouse_genome/refdata-gex-GRCm39-2024-A"  # Update with correct reference genome

# Create output directory if it does not exist
mkdir -p ${OUTPUT_DIR}

# Define sample names and their corresponding conditions
declare -A SAMPLES
SAMPLES["12801-SL-01"]="TCR1-Control"
SAMPLES["12801-SL-02"]="TCR1-Control"
SAMPLES["12801-SL-03"]="TCR2-Naive_T_cells"
SAMPLES["12801-SL-04"]="TCR2-Naive_T_cells"
SAMPLES["12801-SL-05"]="TCR3-Adaptive_T_cells"
SAMPLES["12801-SL-06"]="TCR3-Adaptive_T_cells"

# Loop through each sample and run Cell Ranger Count
for SAMPLE in "${!SAMPLES[@]}"; do
    echo "Processing ${SAMPLE} - ${SAMPLES[$SAMPLE]}"

    cellranger count \
        --id=${SAMPLE} \
        --fastqs=${DATA_DIR} \
        --sample=${SAMPLE} \
        --transcriptome=${GENOME_REFERENCE} \
        --localcores=8 \
        --localmem=228 \
        --expect-cells=5000 \
        --output-dir=${OUTPUT_DIR}/${SAMPLES[$SAMPLE]} \
        --create-bam=true   # Required argument to generate BAM files

    echo "Finished processing ${SAMPLE}"
done

echo "All samples have been processed successfully!"
