#!/bin/bash

pruned1KGplusUKBB="/hpf/largeprojects/struglis/cfcentre/strug/Naim/uk_biobank/Ethnic_PCA/12-unrelated_1KG_plus_UKBiobank_pruned_set"
pcadir="data/intermediate_files/pca"
copddata="data/intermediate_files/GOLD2-4_copd_ukbb_spirodata.csv"

[ ! -d "$pcadir" ] && mkdir "$pcadir"

# The 1000 Genomes and UKBB genotype arrays have been combined and pruned before, so create symbolic links to them
ln -s "${pruned1KGplusUKBB}.bed" "$pcadir"
ln -s "${pruned1KGplusUKBB}.bim" "$pcadir"
ln -s "${pruned1KGplusUKBB}.fam" "$pcadir"
ln -s "${pruned1KGplusUKBB}.log" "$pcadir"

# intermediate files:
keepfile="${pcadir}/ukbb_gold2-4_copd_cases.txt"
plinkout="${pcadir}/13-ukbb_copd_set"
kingout="${pcadir}/14-copdset_kinshipmat"
pcairout="${pcadir}/15-ukbb_copd_pcair"

sed '1d' "$copddata" |awk -F "," '{print $1"\t"$1}' > "$keepfile"
plink --bfile "$pruned1KGplusUKBB" --keep "$keepfile" --make-bed --out "$plinkout"
king -b "${plinkout}.bed" --kinship --prefix "$kingout"

# Run PCA for the set with COPD (i.e. GOLD2-4)
Rscript ~/scripts/R_scripts/PC-AiR_v1_empty_kinFile2.R "$plinkout" "$kingout" "$pcairout"



# Next, we want to run PCA for an analysis of all UKBB with good spirometry/genotype data:
ukbbdataQCd="data/clean/ukbb_spiro_and_geno_qc.csv"
keepfile="${pcadir}/ukbb_good_spirometry_and_genoQC_set.txt"
plinkout="${pcadir}/16-ukbb_spiro_set"
#kingout="${pcadir}/17-ukbbspiro_kinshipmat"
flashpcaout="${pcadir}/18-ukbb_ukbbspiro_flashpca2"

sed '1d' "$ukbbdataQCd" |awk -F "," '{print $1"\t"$1}' > "$keepfile"
plink --bfile "$pruned1KGplusUKBB" --keep "$keepfile" --make-bed --out "$plinkout"
#king -b "${plinkout}.bed" --kinship --prefix "$kingout" # very large matrix and not using it in PCA step below

# Run PCA for the set of UKBB individuals that underwent spirometry and genotyping QC
# Rscript ~/scripts/R_scripts/PC-AiR_v1_empty_kinFile2.R "$plinkout" "$kingout" "$pcairout"
# fails due to very large memory requirement; run flashPCA2 instead

Rscript ~/scripts/run_flashPCA2.R "$plinkout" "$flashpcaout"


# qsub code/03-pca.sh -l walltime=23:59:00 -l nodes=1:ppn=1 -l mem=110g -l vmem=110g -o ./pcajobout -e ./pcajobout -d `pwd` -N copd_pca

