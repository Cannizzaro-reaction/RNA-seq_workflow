#!/bin/bash

# create working directory
working_dir=$(pwd)
mkdir -p "$working_dir/RNA-seq_analysis/database"
cd RNA-seq_analysis

#################### Set up the rna_seq environment ####################
# Add conda mirror channels
conda config --add channels https://mirror.sjtu.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirror.sjtu.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirror.sjtu.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirror.sjtu.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes

# Create and activate conda environment
conda create --name rna_seq -y
source activate rna_seq

# Install necessary packages
sudo apt-get update
sudo apt-get install -y unzip wget fastqc bwa samtools

# Download and install SRA tools (used in prefetch)
SRA_VERSION="3.0.0"
INSTALL_DIR="$HOME/sratoolkit.$SRA_VERSION-ubuntu64"
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/$SRA_VERSION/sratoolkit.$SRA_VERSION-ubuntu64.tar.gz -O sratoolkit.tar.gz
tar -xvzf sratoolkit.tar.gz -C $HOME
rm sratoolkit.tar.gz
export PATH=$PATH:$INSTALL_DIR/bin

# Download and install trim_galore
mkdir software
cd software
wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.7.tar.gz
tar -xzvf 0.6.7.tar.gz
cd TrimGalore-0.6.7
chmod +x trim_galore
sudo cp trim_galore /usr/local/bin
cd ..
rm 0.6.7.tar.gz
cd ..

# Install subread for expression quantification analysis
conda install -c bioconda subread -y

#################### Set up the r_env environment ####################
conda create --name r_env -y
source activate r_env
conda install -c r r-base -y

# Install R packages for differential expression analysis
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
Rscript -e "BiocManager::install('DESeq2')"
Rscript -e "BiocManager::install('EnhancedVolcano')"
Rscript -e "BiocManager::install('pheatmap', force = TRUE)"
Rscript -e "if (!requireNamespace('ggplot2', quietly = TRUE)) install.packages('ggplot2')"

echo "Environment setup and installation completed."
