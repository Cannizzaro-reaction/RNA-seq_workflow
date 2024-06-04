#!/bin/bash

set -e 
LOGFILE="setup_log.txt"

# Function to log messages
log() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') - $1" | tee -a $LOGFILE
}

log "Starting environment setup..."

# Check if conda is installed
if ! command -v conda &> /dev/null
then
    log "Conda could not be found. Please install Conda and rerun this script."
    exit 1
fi

# Create working directory
working_dir=$(pwd)
log "Setting up working directory at $working_dir/RNA-seq_analysis"
mkdir -p "$working_dir/RNA-seq_analysis/database"
cd RNA-seq_analysis

#################### Set up the rna_seq environment ####################
# Add conda mirror channels
log "Adding conda mirror channels..."
conda config --add channels https://mirror.sjtu.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirror.sjtu.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirror.sjtu.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirror.sjtu.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes

# Load Conda environment
log "Loading Conda environment..."
source $(conda info --base)/etc/profile.d/conda.sh

# Create and activate conda environment
if conda env list | grep -q "rna_seq"; then
    log "Environment 'rna_seq' already exists. Activating..."
else
    log "Creating and activating 'rna_seq' environment..."
    conda create --name rna_seq -y
fi
conda activate rna_seq

# Install necessary packages via conda
log "Installing necessary packages for rna_seq..."
conda install -c bioconda sra-tools fastqc trim-galore bwa samtools subread -y

#################### Set up the r_env environment ####################
if conda env list | grep -q "r_env"; then
    log "Environment 'r_env' already exists. Activating..."
else
    log "Creating and activating 'r_env' environment..."
    conda create --name r_env -y
fi
conda activate r_env
conda install -c r r-base -y

# Install R packages for differential expression analysis
log "Installing R packages..."
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" >> $LOGFILE 2>&1
Rscript -e "BiocManager::install('DESeq2')" >> $LOGFILE 2>&1
Rscript -e "BiocManager::install('EnhancedVolcano')" >> $LOGFILE 2>&1
Rscript -e "BiocManager::install('pheatmap')" >> $LOGFILE 2>&1
Rscript -e "if (!requireNamespace('ggplot2', quietly = TRUE)) install.packages('ggplot2')" >> $LOGFILE 2>&1

log "Environment setup and installation completed."
