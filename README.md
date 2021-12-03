# Introduction

This project contains the steps undertaken to perform association analyses of spirometry measures in the UK Biobank (UKBB) with variants around _SLC26A9_ gene for patients with spirometry-defined COPD with modified [GOLD criteria 2-4 of moderate to very severe lung function](https://pubmed.ncbi.nlm.nih.gov/17765526/). The definition here is relaxed in that measurements are not required to be post-bronchodilator measurements. According to [Mannino and Buist, 2007](https://pubmed.ncbi.nlm.nih.gov/17765526/), in instances where pre-bronchodilator lung function has been recorded, an overestimate of airflow obstruction may result.


# Analysis

  1. [01-extract_phenos_of_interest.sh](code/01-extract_phenos_of_interest.sh)  
      - input:  
          - `ukb24727.tab`, which contains all phenotypic information from UKBB        
      - output:  
          - `ukb24727_spirometry.tab`, a smaller file containing the required variables only  

  2. [02-subset_qc_copd_individuals.R](code/02-subset_qc_copd_individuals.R)  
      - input:  
          - `ukb24727_spirometry.tab`
      - output:  
          - `ukbb_spiro_and_geno_qc.csv`, which contains all the individuals passing spirometry and genotyping QC, and their spirometry measures (best FEV<sub>1</sub>, best FVC, and FEV1pp)
          - `GOLD2-4_copd_ukbb_spirodata.csv`, which is the subset of individuals from `ukbb_spiro_and_geno_qc.csv` that fit GOLD class 2-4 criteria for lung function (i.e. FEV<sub>1</sub>/FVC ratio &lt; 0.7 and FEV1pp &lt; 80&percnt;)  
          
      This step removes individuals that did not pass spirometry and genotyping QC, removes related and non-European individuals, and calculates FEV1pp using the [GLI calculator](http://gli-calculator.ersnet.org/index.html) and inverse-rank normalizes the spirometry measures  

  3. [03-pca.sh](code/03-pca.sh)  
      - input: genotype array data for:  
          - All individuals defined in `GOLD2-4_copd_ukbb_spirodata.csv`  
          - All individuals defined in `ukbb_spiro_and_geno_qc.csv`
      - output:  
          - `15-ukbb_copd_pcair_eigenvectors.txt` for individuals corresponding to `GOLD2-4_copd_ukbb_spirodata.csv`  
          - `18-ukbb_ukbbspiro_flashpca2_eigenvectors.txt` for individuals corresponding to `ukbb_spiro_and_geno_qc.csv`

  4. `04-assoc_spirometry.R`  
      - input:
          - d

  5. locusfocus_prep/merge_and_convert_to_html.py - generates locusfocus_prep/ukbb_spirometry_assoc_in_copd.html for use in LocusFocus; this step prepares the association results for loading as secondary datasets into [LocusFocus](https://locusfocus.research.sickkids.ca)  


