#!/bin/bash

# Extract COPD phenotypes of interest from UKBB application file:
phenofile="/hpf/largeprojects/struglis/datasets/uk_biobank_40946/phenotypes/ukb24727.tab"
# Columns of interest (looked up from ukb24727.html):
# 1: EID
# 8306: doctor-diagnosed COPD (1,768 cases, 119,516 controls, 381,259 NA)
# 8322: age COPD diagnosed by doctor
# 8338: recent medication for COPD
# 8193: genetic ethnic grouping. data coding 1002 (409,634 Caucasian (1), 92,902 NA)
# 8188: genetic sex. data-coding 9 (1=male=223,481; 0=female=264,814)
# 8236: recommended genomic analysis exclusions. data-coding 1101
# 8252: genetic relatedness exclusions. data-coding 947 (692 mixed ancestral background (1); 840 high heterozygosity or missing (2); 501,011 NA)
# 22: year of birth
# 8149: age at recruitment
# 8255: genetic kinship to other participants. data-coding 682 (339,604 no kinship found (0); etc.)
# 1166-1174: FEV1
# 7417: FEV1 best measure
# 7421: FEV1 predicted percentage
# 1157-1165: FVC
# 7418: FVC best measure
# 1175-1183: PEF
# 1148-1156: Acceptability of each blow result - decimal derived from hexadecimal sum
# 7419: Reproduciblity of spirometry measurement using ERS/ATS criteria - yes/no
# 6345-6353: Acceptability of each blow result (text)


cols="1,8306,8322,8338,8193,8188,8236,8252,22,8149,8255,1166-1174,7417,7421,1157-1165,7418,1175-1183,1148-1156,7419,6345-6353"
echo "Extracting desired COPD and spirometry measure columns"
cut -f"$cols" "$phenofile" > ukb24727_copd_and_spirometry2.tab # added field 3061 (columns 1148-1156) for acceptability of each blow result
echo "Done"

#qsub extract_phenos_of_interest.sh -l walltime=23:59:00 -l nodes=1:ppn=1 -l mem=8g -l vmem=8g -o ./jobout -e ./jobout -d `pwd` -N ukbb_pheno_extraction

# from previous analyses:
ln -s /hpf/largeprojects/struglis/datasets/uk_biobank_40946/phenotypes/set_of_related_ind_to_rm.txt # (36,100 samples)
