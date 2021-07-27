#!/bin/bash

ln -s /hpf/largeprojects/struglis/cfcentre/strug/Naim/uk_biobank/Ethnic_PCA/12-unrelated_1KG_plus_UKBiobank_pruned_set.bed
ln -s /hpf/largeprojects/struglis/cfcentre/strug/Naim/uk_biobank/Ethnic_PCA/12-unrelated_1KG_plus_UKBiobank_pruned_set.bim
ln -s /hpf/largeprojects/struglis/cfcentre/strug/Naim/uk_biobank/Ethnic_PCA/12-unrelated_1KG_plus_UKBiobank_pruned_set.fam
ln -s /hpf/largeprojects/struglis/cfcentre/strug/Naim/uk_biobank/Ethnic_PCA/12-unrelated_1KG_plus_UKBiobank_pruned_set.log

sed '1d' ../ukbb24727_spirometry_in_copd_irnt.tsv |awk '{print $1"\t"$1}' > ukbb_copd_cases.txt
#plink --bfile 12-unrelated_1KG_plus_UKBiobank_pruned_set --keep ukbb_copd_cases.txt --make-bed --out 13-ukbb_copd_cases_only
#king -b 13-ukbb_copd_cases_only.bed --kinship --prefix 14-copd_kinshipmat 
#Rscript ~/scripts/R_scripts/PC-AiR_v1_empty_kinFile2.R 13-ukbb_copd_cases_only 14-copd_kinshipmat 15-ukbb_copd_pcair

keepfile="fev_copd_set.txt"
plinkout="13-ukbb_copd_fev_set"
kingout="14-fevset_kinshipmat"
pcairout="15-ukbb_copd_fevset_pcair"
plink --bfile 12-unrelated_1KG_plus_UKBiobank_pruned_set --keep "$keepfile" --make-bed --out "$plinkout"
king -b "${plinkout}.bed" --kinship --prefix "$kingout"
Rscript ~/scripts/R_scripts/PC-AiR_v1_empty_kinFile2.R "$plinkout" "$kingout" "$pcairout"

keepfile="bestpef_set.txt"
plinkout="13-ukbb_copd_bestpef_set"
kingout="14-bestpefset_kinshipmat"
pcairout="15-ukbb_copd_bestpefset_pcair"
plink --bfile 12-unrelated_1KG_plus_UKBiobank_pruned_set --keep "$keepfile" --make-bed --out "$plinkout"
king -b "${plinkout}.bed" --kinship --prefix "$kingout"
Rscript ~/scripts/R_scripts/PC-AiR_v1_empty_kinFile2.R "$plinkout" "$kingout" "$pcairout"

keepfile="fev1pp_set.txt"
plinkout="13-ukbb_copd_fev1pp_set"
kingout="14-fev1ppset_kinshipmat"
pcairout="15-ukbb_copd_fev1ppset_pcair"
plink --bfile 12-unrelated_1KG_plus_UKBiobank_pruned_set --keep "$keepfile" --make-bed --out "$plinkout"
king -b "${plinkout}.bed" --kinship --prefix "$kingout"
Rscript ~/scripts/R_scripts/PC-AiR_v1_empty_kinFile2.R "$plinkout" "$kingout" "$pcairout"

