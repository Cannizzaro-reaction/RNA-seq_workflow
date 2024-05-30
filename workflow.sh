#!/bin/bash

#################### workflow for rna-seq analysis ####################
#################### Upstream Analysis ####################
# Activate conda environment
conda activate rna_seq

# Create necessary directories
mkdir -p database result/qc result/trim result/expression result/bam result/ref_compare

# Download SRA data
echo "Downloading SRA data..."
cd database
for sra in SRR15174659 SRR15174661 SRR15174663 SRR15174670 SRR15174671 SRR15174672; do
    prefetch $sra
    fasterq-dump $sra/$sra.sra -O .
    rm -r $sra
done
echo "SRA data downloaded."

# Download reference genome data (Klebsiella pneumoniae)
echo "Downloading reference genome data..."
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/694/555/GCA_000694555.1_Kleb_pneu_MGH_66_V1/GCA_000694555.1_Kleb_pneu_MGH_66_V1_genomic.fna.gz
gunzip -c GCA_000694555.1_Kleb_pneu_MGH_66_V1_genomic.fna.gz > reference.fna
rm GCA_000694555.1_Kleb_pneu_MGH_66_V1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/694/555/GCA_000694555.1_Kleb_pneu_MGH_66_V1/GCA_000694555.1_Kleb_pneu_MGH_66_V1_genomic.gtf.gz
gunzip -c GCA_000694555.1_Kleb_pneu_MGH_66_V1_genomic.gtf.gz > reference.gtf
rm GCA_000694555.1_Kleb_pneu_MGH_66_V1_genomic.gtf.gz
echo "Reference genome data downloaded."

# Convert .sra files to fastq format
echo "Converting .sra files to fastq format..."
for sra in SRR15174659 SRR15174661 SRR15174663 SRR15174670 SRR15174671 SRR15174672; do
    fastq-dump $sra.sra
done
echo ".sra files converted to fastq format."

# Use fastqc to evaluate the quality of each .fastq file (raw data), the result will be saved to result/qc
echo "Evaluating fastq file quality with fastqc..."
cd ..
fastqc -t 6 -o ./result/qc/ ./database/SRR*.fastq
echo "Fastq file quality evaluation completed."

# Trim the original fastq files using trim_galore
echo "Trimming fastq files with trim_galore..."
cd database
for sra in SRR15174659 SRR15174661 SRR15174663 SRR15174670 SRR15174671 SRR15174672; do
    trim_galore --output_dir ./ ${sra}.fastq
    mv ${sra}_trimmed.fq ../result/trim/${sra}_trimmed.fq
    mv ${sra}.fastq_trimming_report.txt ../result/trim/${sra}.fastq_trimming_report.txt
done
echo "Fastq file trimming completed."

# Align the trimmed fastq files with the reference genome data and generate .sam files
echo "Aligning trimmed fastq files with the reference genome..."
bwa index reference.fna
for sra in SRR15174659 SRR15174661 SRR15174663 SRR15174670 SRR15174671 SRR15174672; do
    bwa mem reference.fna ../result/trim/${sra}_trimmed.fq > ../result/ref_compare/${sra}_alignment.sam
done
echo "Alignment completed."

# Convert .sam files to .bam format, sort and index the bam files
echo "Converting .sam files to .bam format, sorting, and indexing..."
for sra in SRR15174659 SRR15174661 SRR15174663 SRR15174670 SRR15174671 SRR15174672; do
    samtools view -S -b ../result/ref_compare/${sra}_alignment.sam > ../result/bam/${sra}_alignment.bam
    samtools sort ../result/bam/${sra}_alignment.bam -o ../result/bam/${sra}_sorted.bam
    samtools index ../result/bam/${sra}_sorted.bam
done
echo ".bam file conversion, sorting, and indexing completed."

# Expression quantification using featureCounts
echo "Performing expression quantification with featureCounts..."
featureCounts -a reference.gtf -o ../result/expression/counts.txt -g ID ../result/bam/*_sorted.bam
echo "Expression quantification completed."

# Clean up intermediate files
echo "Cleaning up intermediate files..."
rm database/*.sra
rm database/*.fastq
rm result/trim/*.fq
rm -r result/ref_compare
echo "Clean up completed."

echo "RNA-seq analysis upstream workflow completed."

#################### Downstream Analysis and Data Visualization ####################

# Process the gene expression outcome
cd result/expression
awk 'NR > 1 { printf "%s", $1; for (i = 7; i <= NF; i++) printf "\t%s", $i; print "" }' counts.txt | sed -e 's/result\/bam\///g' -e 's/_sorted.bam//g' > raw_counts.txt

cd ../..
conda activate r_env
Rscript deseq_analysis.R
mv differential_expression_results.csv result/expression/differential_expression_results.csv
mv Rplots.pdf result/expression/Rplots.pdf
echo "Expression analysis and data visualization completed."

# Final structure confirmation
echo "Directory structure:"
tree