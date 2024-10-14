#!/usr/bin/env bash

# %% Make data directory

mkdir -p data

# %% Generate read quality control reports

mkdir -p data/fastqc
fastqc \
    -t 8 \
    -o data/fastqc \
    biohub-data/Peter_Colias/*.fastq.gz

# %% Trim adapters

mkdir -p data/trimmed
for name in \
    "WT_PJC_S1" \
    "CD55off_PJC_S2" \
    "WT_JPL_S3" \
    "CD55off_JPL_S4"
do
    echo "Running cutadapt on ${name}..."
    cutadapt \
        --cores=0 \
        -m 1 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o "data/trimmed/${name}_R1.fastq.gz" \
        -p "data/trimmed/${name}_R2.fastq.gz" \
        "biohub-data/Peter_Colias/${name}_R1_001.fastq.gz" \
        "biohub-data/Peter_Colias/${name}_R2_001.fastq.gz"
done

# %% Download and unzip human genome (hg38) with spike-in EM-seq controls

mkdir -p data/em-seq-genome
wget \
    https://neb-em-seq-sra.s3.amazonaws.com/grch38_core%2Bbs_controls.fa \
    -O data/em-seq-genome/hg38.fasta

# # Rename genome file for Bismark
#  mv \
#    data/hg38/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna \
#    data/hg38/data/GCF_000001405.40/genome.fasta

# %% Bismark genome preparation

bismark_genome_preparation \
    --verbose \
    data/em-seq-genome/

#%% Bismark alignment

mkdir -p data/aligned

for name in \
    "WT_PJC_S1" \
    "CD55off_PJC_S2" \
    "WT_JPL_S3" \
    "CD55off_JPL_S4"
do
    bismark \
       --bam \
       --parallel 20 \
       --genome data/em-seq-genome/ \
       -o data/aligned \
       -1 data/trimmed/${name}_R1.fastq.gz \
       -2 data/trimmed/${name}_R2.fastq.gz
done

# %% Bismark methylation extraction

mkdir -p data/methylation

for name in \
    "WT_PJC_S1" \
    "CD55off_PJC_S2" \
    "WT_JPL_S3" \
    "CD55off_JPL_S4"
do
    bismark_methylation_extractor \
        --parallel 20 \
        -o data/methylation \
        data/aligned/${name}_R1_bismark_bt2_pe.bam
done

# %% Convert Bismark files to bedGraph

mkdir -p data/methylation-bedGraph

for name in \
    "WT_PJC_S1" \
    "CD55off_PJC_S2" \
    "WT_JPL_S3" \
    "CD55off_JPL_S4"
do
    # Original top (OT) strand
    bismark2bedGraph \
        --dir data/methylation-bedGraph/ \
        -o ${name}_OT.txt \
        data/methylation/CpG_OT_${name}_R1_bismark_bt2_pe.txt.gz &

    # Original bottom (OB) strand
    bismark2bedGraph \
        --dir data/methylation-bedGraph/ \
        -o ${name}_OB.txt \
        data/methylation/CpG_OB_${name}_R1_bismark_bt2_pe.txt.gz &
done

wait

# %% Convert Bismark bedGraph files to DSS input files

mkdir -p data/methylation-dss-separate

for name in \
    "WT_PJC_S1" \
    "CD55off_PJC_S2" \
    "WT_JPL_S3" \
    "CD55off_JPL_S4"
do
    ./scripts/bismark_to_dss.py \
        data/methylation-bedGraph/${name}_OT.txt.gz.bismark.cov.gz \
        data/methylation-dss-separate/${name}_OT.txt &

    ./scripts/bismark_to_dss.py \
        data/methylation-bedGraph/${name}_OB.txt.gz.bismark.cov.gz \
        data/methylation-dss-separate/${name}_OB.txt &
done

wait

# %% Aggregate original top and bottom strands

mkdir -p data/methylation-dss-combined

for name in \
    "WT_PJC_S1" \
    "CD55off_PJC_S2" \
    "WT_JPL_S3" \
    "CD55off_JPL_S4"
do
    uv \
        --project scripts/ \
        run ./scripts/combine_dss.py \
        data/methylation-dss-separate/${name}_OT.txt \
        data/methylation-dss-separate/${name}_OB.txt \
        data/methylation-dss-combined/${name}.txt &
done

wait

# % % Filter low read counts

mkdir -p data/methylation-dss-combined-filtered

uv \
    --project scripts/ \
    run ./scripts/filter_dss.py \
    10 \
    data/methylation-dss-combined \
    data/methylation-dss-combined-filtered

# % % Call DSS to perform differential methylation analysis

mkdir -p data/dss-output/
./scripts/dss.r

# % % Plot results

mkdir -p graphs
uv \
    --project scripts/ \
    run scripts/plot.py \
    data/dss-output/dss.csv \
    graphs/
