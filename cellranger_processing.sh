#!/usr/bin/bash

# cellranger_processing.sh
# Script to process 10X Genomics single-cell RNA-seq data using Cell Ranger

# SLURM job submission script for Cell Ranger count
cat > run_cellranger.sh << 'EOF'
#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=100G
#SBATCH --job-name="cellranger"

# Set these variables before running
FASTQ_DIR=""        # Path to directory containing FASTQ files
SAMPLE_ID=""        # Sample identifier
TRANSCRIPTOME=""    # Path to Cell Ranger reference transcriptome

# Run Cell Ranger count
cellranger count --id="run_count_${SAMPLE_ID}" \
    --fastqs="${FASTQ_DIR}" \
    --sample="${SAMPLE_ID}" \
    --transcriptome="${TRANSCRIPTOME}" \
    --localcores=10 \
    --localmem=95
EOF

chmod +x run_cellranger.sh

echo "Created Cell Ranger processing script: run_cellranger.sh"
echo "Before running, please set FASTQ_DIR, SAMPLE_ID, and TRANSCRIPTOME variables"
echo ""
echo "Example usage:"
echo "export FASTQ_DIR=/path/to/fastqs"
echo "export SAMPLE_ID=sample1"
echo "export TRANSCRIPTOME=/path/to/reference"
echo "sbatch run_cellranger.sh"
