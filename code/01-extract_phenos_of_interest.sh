#!/bin/bash

date

# Extract COPD phenotypes of interest from UKBB application file:
phenofile="/hpf/largeprojects/struglis/datasets/uk_biobank_40946/phenotypes/ukb24727.tab"

# Note: column number is 1-based (not zero-based); so, add 1 to provided column number in data dictionary
# eid (column 1)
# age at recruitment: 21022 (column 8149)
# standing height: 50 (columns 68-70)
# genetic sex: 22001 (column 8188)
# genetic ethnic grouping (Caucasian or not, by PCA): 22006 (column 8193)
# Spirometry method: 23 (columns 18-20)
# FEV1: 3063 (columns 1166-1174)
# FVC: 3062 (columns 1157-1165)
# PEF: 3064 (columns 1175-1183)
# blow curve time series: 3066 (columns 1193-1201)
# Vitalograph spirometer blow quality metrics: 20031 (columns 6345-6353)
# Recommended genomic analysis exclusions: 22010 (column 8236)
# Genetic relatedness exclusions: 22018 (column 8252)
# Genetic kinship to other participants: 22021 (column 8255)
# Ever smoked: 20160 (column 7450)

cols="1,8149,68-70,8188,8193,18-20,1166-1174,1157-1165,1175-1183,1193-1201,6345-6353,8236,8252,8255,7450"
echo "Extracting desired spirometry measure columns to assess GOLD-defined moderate-severe (2-4) lung function"
cut -f"$cols" "$phenofile" > data/intermediate_files/ukb24727_spirometry.tab
echo "Done"
date

# from previous analyses (36,100 samples):
cp /hpf/largeprojects/struglis/datasets/uk_biobank_40946/phenotypes/set_of_related_ind_to_rm.txt data/intermediate_files/

#qsub extract_phenos_ukbb.sh -l walltime=23:59:00 -l nodes=1:ppn=1 -l mem=8g -l vmem=8g -o ./jobout -e ./jobout -d `pwd` -N ukbb_pheno_extraction
