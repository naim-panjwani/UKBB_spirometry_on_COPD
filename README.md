# Introduction

This project contains the steps that were undertaken to perform association analyses of spirometry measures in the UK Biobank (UKBB) with variants around _SLC26A9_ for patients with COPD.


# Steps in order

  1. extract_phenos_of_interest.sh - generates ukb24727_copd_and_spirometry.tab, a smaller file that may be more easily loaded in R  
  2. subset_qc_copd_individuals.R - generates ukbb24727_spirometry_in_copd_irnt.tsv; this step removes individuals that did not pass QC and inverse-rank normalizes the spirometry measures  
  3. pca/steps.sh - generates pca/15-ukbb_copd_pcair_eigenvectors.txt; this step conducts principal component analysis (PCA) for the individuals that will be analyzed  
  4. assoc_copd_spirometry.R - generates association files in assoc/ folder; this step carries out the association analyses for each spirometry measure of interest  
  5. locusfocus_prep/merge_and_convert_to_html.py - generates locusfocus_prep/ukbb_spirometry_assoc_in_copd.html for use in LocusFocus; this step prepares the association results for loading as secondary datasets into [LocusFocus](https://locusfocus.research.sickkids.ca)  


