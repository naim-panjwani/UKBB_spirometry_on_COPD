library(data.table)
library(RNOmni)
library(httr)
library(jsonlite)
source(".apikey")

# Load data
dat <- fread(file.path("data","intermediate_files","ukb24727_spirometry.tab"), header=T, stringsAsFactors=F) # 502,543

# Take individuals with at least two FEV1 measures (Variable ID: 3063) and FVC (3062), complete info for spirometry method used (23), age (21022), genetic sex(22001), and standing height (50)
num_fvc_missing <- apply(dat[,grep("3062.0",colnames(dat)),with=F], 1, function(x) sum(is.na(x)))
num_fev_missing <- apply(dat[,grep("3063.0",colnames(dat)),with=F], 1, function(x) sum(is.na(x)))

dat <- dat[ (num_fev_missing < 2) & (num_fvc_missing < 2) ] # 453,454
dat <- dat[ f.23.0.0 %in% c(0,1) ] # 453,443 - 0=direct entry; 1=manual
dat <- dat[ !is.na(f.21022.0.0) ] # no age missing; min age: 37, max age: 73
dat <- dat[ !is.na(f.22001.0.0) ] # 444,620
dat <- dat[ !is.na(f.50.0.0) ] # 444,052

# Spirometry QC as per Shrine et al.
## Acceptability of blows
acceptable_blows <- c("ACCEPT", "BELOW6SEC ACCEPT", "BELOW6SEC")
acceptable_blows_mat_idx <- apply(dat[ ,grep("20031.0", colnames(dat)),with=F ], 2, function(x) ifelse(is.na(x) | (x %in% acceptable_blows), 1, 0)) # 891,784 acceptable blows from 444,052 participants

## Assess start of blow quality
source("code/back_ev_calculation.R")
blow_curves_mat <- dat[, grep("3066.0", colnames(dat)), with=F]
acceptable_blow_start <- matrix(rep(NA, nrow(blow_curves_mat)*ncol(blow_curves_mat)), ncol=ncol(blow_curves_mat))
for(i in 1:nrow(blow_curves_mat)) {
  if((i %% 10000) == 0) print(paste0(i,"/",nrow(acceptable_blow_start)," ",format(round(i/nrow(acceptable_blow_start)*100,digits = 2),nsmall=1),"% done"))
  for(j in 1:ncol(blow_curves_mat)) {
    acceptable_blow_start[i,j] <- blow_start_quality(blow_curves_mat[i,j,with=F])
  }
}
# 968,546 acceptable blow starts


x <- unname(acceptable_blows_mat_idx * acceptable_blow_start)
x <- ifelse(x==0,NA,x) # 703,772 acceptable blows and blow starts in 365,700 participants
participants_to_rm <- rowSums(is.na(x))==3 # 78,352 marked for removal

eid <- dat[,'f.eid']
fvc_mat <- dat[, grep("3062.0", colnames(dat)), with=F]
fev_mat <- dat[, grep("3063.0", colnames(dat)), with=F]

fvc_mat <- fvc_mat * x
fvc_mat <- fvc_mat[-which(participants_to_rm),]
fev_mat <- fev_mat * x
fev_mat <- fev_mat[-which(participants_to_rm),]
eid <- eid[-which(participants_to_rm),]

## Get best measures
fvc_best <- apply(fvc_mat,1,max,na.rm=T)
fev_best <- apply(fev_mat,1,max,na.rm=T)

## Assess reproducibility of measures
### best measures have to be within 250mL from any other blow (including unacceptable blows)
fvc_mat2 <- dat[-which(participants_to_rm), grep("3062.0", colnames(dat)), with=F]
fev_mat2 <- dat[-which(participants_to_rm), grep("3063.0", colnames(dat)), with=F]
fvc_reproducible <- abs(fvc_mat2 - fvc_best)
fvc_reproducible <- apply(fvc_reproducible, 1, function(x) !any(x>0.25,na.rm=T)) # 231,566
fev_reproducible <- abs(fev_mat2 - fev_best)
fev_reproducible <- apply(fev_reproducible, 1, function(x) !any(x>0.25,na.rm=T)) # 273,478

fev_and_fvc_reproducible <- fev_reproducible & fvc_reproducible # 214,928

fvc_best_reproducible <- fvc_best[fev_and_fvc_reproducible] # 214,928
fev_best_reproducible <- fev_best[fev_and_fvc_reproducible] # 214,928
fev_fvc_ratio <- fev_best_reproducible / fvc_best_reproducible # 214,928
eid <- eid[fev_and_fvc_reproducible,]

dat <- dat[-which(participants_to_rm),] # 365,700 with acceptable blows
dat <- dat[fev_and_fvc_reproducible,] # 214,928 with reproducible blows

spiro_qc_df <- cbind(eid, fvc_best_reproducible, fev_best_reproducible, fev_fvc_ratio) # 214,928
# 30,863 have FEV1/FVC ratio < 0.7

dat <- cbind(dat, spiro_qc_df[,-"f.eid",with=F])

#==================================================
# 1. Remove failed QC individuals (column 8236, UDI=22010-0.0; poor heterozygosity/missingness)
dat <- dat[!(f.22010.0.0 %in% 1)] # 214,717; 211 removed

# 2. Remove individuals with no genetic sex information
dat <- dat[!is.na(f.22001.0.0)] # 214,717;; none removed

# 3. Remove related individuals
relinds <- fread(file.path("data","intermediate_files", "set_of_related_ind_to_rm.txt"), header=F, stringsAsFactors=F) # 36,100
# check all fid==iid
# for(i in 1:nrow(relinds)) {
#   if(relinds[i,1] != relinds[i,2]) print(relinds[i,])
# }
#V1      V2
#1: VN061 HG02061 --> not in the dataset
relinds <- relinds[V2!="HG02061"]
relinds$V1 <- as.integer(relinds$V1)
relinds$V2 <- as.integer(relinds$V2) # 36,099
dat <- dat[!(f.eid %in% relinds$V2)] # 198,732; 15,985 removed

# 4. Remove non-Caucasian (column 8193; 22006-0.0)
dat <- dat[f.22006.0.0==1] # 168,340; 30,392 removed

# 5. Save the current dataset to calculate FEV1pp on GLI calculator
gli_query_data <- dat[, grep("eid|21022.0|50.0|22001.0|22006.0|fev_best_reproducible|fvc_best_reproducible", colnames(dat)), with=F]

setnames(gli_query_data,
         c("f.eid", "f.50.0.0", "f.21022.0.0", "f.22001.0.0", "f.22006.0.0", "fvc_best_reproducible", "fev_best_reproducible"),
         c("eid", "height","age","sex","ethnic","fvc","fev1"))
setcolorder(gli_query_data,
            c("eid","age","height","sex","ethnic","fev1","fvc"))
gli_query_data[sex := ifelse(sex==0,'F',ifelse(sex==1,'M',NA))]
fwrite(gli_query_data
       ,file.path("data", "intermediate_files","gli_calc_data.csv")
       ,quote=F, row.names=F, col.names=T, sep=",")

# 7. Use GLI Calculator API to calculate FEV1pp
# data retrieval limited to 100k per month, so subset into smaller set
gli_query_data_lowfunc <- gli_query_data[(fev1/fvc)<0.7] # 24,953
rest_url <- "https://gli-api.ersnet.org/public/"
fev1pp <- NULL
results <- NULL
for(i in 1:nrow(gli_query_data_lowfunc)) {
  print(paste0("Acquiring sample ",i," of ",nrow(gli_query_data_lowfunc),
               " ",format(round((i/nrow(gli_query_data_lowfunc))*100
                                ,digits = 2),nsmall=2), "% done"))
  response <- GET(paste0(rest_url,"type/spiro"
                         ,"/age/",gli_query_data_lowfunc$age[i]
                         ,"/height/",gli_query_data_lowfunc$height[i]
                         ,"/sex/",tolower(gli_query_data_lowfunc$sex[i])
                         ,"/ethnic/",gli_query_data_lowfunc$ethnic[i]
                         ,"/fev1/",gli_query_data_lowfunc$fev1[i]
                         ,"/fvc/",gli_query_data_lowfunc$fvc[i]
                         )
                  , add_headers("x-api-key" = apikey))
  if(response$status_code==200){
    spiro <- content(response, type="application/json")
    fev1pp <- as.numeric(spiro$fev1_pp)
    results <- rbind(results, c(eid=gli_query_data_lowfunc$eid[i], fev1pp=fev1pp))
  } else {
    message <- content(response, type="application/json")$message
    print(paste0("API status code ",response$status_code,": ",message))
    break
  }
}
x <- data.table(results)
gli_query_data_lowfunc <- merge(gli_query_data_lowfunc, x, by='eid')
# remove one extreme outlier result with fev1pp returned as 5789855.000
gli_query_data_lowfunc <- gli_query_data_lowfunc[fev1pp<200] #24,952
fwrite(gli_query_data_lowfunc, "data/intermediate_files/gli_lowfunc_fev1pp.csv", quote=F, row.names=F, col.names=T, sep=",")


# try to get the rest - later (exceeded limit)
# gli_query_data <- fread("data/intermediate_files/gli_calc_data.csv")
gli_query_data_normfunc <- gli_query_data[(fev1/fvc)>=0.7]
#fwrite(gli_query_data_normfunc, "data/intermediate_files/gli_normfunc.csv", quote=F, row.names=F, col.names=T, sep=",") # save temporarily then add fev1pp after queries complete
#gli_query_data_normfunc <- fread("data/intermediate_files/gli_normfunc.csv")
#gli_query_data_normfunc_fev1pp <- fread("data/intermediate_files/gli_normfunc_fev1pp.csv")
#gli_query_data_normfunc <- gli_query_data_normfunc[-which(eid %in% gli_query_data_normfunc_fev1pp$eid)]
rest_url <- "https://gli-api.ersnet.org/public/"
fev1pp <- NULL
results <- NULL
for(i in 1:nrow(gli_query_data_normfunc)) {
  print(paste0("Acquiring sample ",i," of ",nrow(gli_query_data_normfunc),
               " ",format(round((i/nrow(gli_query_data_normfunc))*100
                                ,digits = 2),nsmall=2), "% done"))
  response <- GET(paste0(rest_url,"type/spiro"
                         ,"/age/",gli_query_data_normfunc$age[i]
                         ,"/height/",gli_query_data_normfunc$height[i]
                         ,"/sex/",tolower(gli_query_data_normfunc$sex[i])
                         ,"/ethnic/",gli_query_data_normfunc$ethnic[i]
                         ,"/fev1/",gli_query_data_normfunc$fev1[i]
                         ,"/fvc/",gli_query_data_normfunc$fvc[i]
  )
  , add_headers("x-api-key" = apikey))
  if(response$status_code==200){
    spiro <- content(response, type="application/json")
    fev1pp <- as.numeric(spiro$fev1_pp)
    results <- rbind(results, c(eid=gli_query_data_normfunc$eid[i], fev1pp=fev1pp))
  } else {
    message <- content(response, type="application/json")$message
    print(paste0("API status code ",response$status_code,": ",message))
    break
  }
}
x <- data.table(results)
gli_query_data_normfunc2 <- gli_query_data[(fev1/fvc)>=0.7]
gli_query_data_normfunc2 <- merge(gli_query_data_normfunc2, x, by='eid')
fwrite(gli_query_data_normfunc2, "data/intermediate_files/gli_normfunc_fev1pp2.csv", quote=F, row.names=F, col.names=T, sep=",")


# Combine all API queries and save:
gli_query_data_normfunc1 <- fread("data/intermediate_files/gli_normfunc_fev1pp1.csv")
gli_query_data_normfunc2 <- fread("data/intermediate_files/gli_normfunc_fev1pp2.csv")
gli_query_data_normfunc <- rbind(gli_query_data_normfunc1, gli_query_data_normfunc2)
fwrite(gli_query_data_normfunc, "data/intermediate_files/gli_normfunc_fev1pp.csv", quote=F, row.names=F, col.names=T, sep=",")
gli_query_data <- fread("data/intermediate_files/gli_calc_data.csv")
gli_query_data_lowfunc <- fread("data/intermediate_files/gli_lowfunc_fev1pp.csv")
x <- rbind(gli_query_data_lowfunc, gli_query_data_normfunc) # one outlier was removed due to extreme fev1pp number
gli_query_data <- merge(gli_query_data, x[,c('eid','fev1pp')], by='eid')
hasCOPD <- with(gli_query_data, (fev1/fvc)<0.7 & fev1pp<80) # as per GOLD2-4 definitions
GOLDlevel <- with(gli_query_data, ifelse((fev1/fvc)<0.7, ifelse(fev1pp>=80,1,
                                                               ifelse(fev1pp>=50,2,
                                                                      ifelse(fev1pp>=30,3,
                                                                             ifelse(fev1pp<30,4,NA)))),NA))
gli_query_data <- cbind(gli_query_data,hasCOPD,GOLDlevel)
samplefile <- fread("/hpf/largeprojects/struglis/datasets/uk_biobank_40946/imputation/sample_files/ukb40946_imp_chr1_v3_s487324.sample")
samplefile <- samplefile[-1,] # 487,409
samplefile$index <- 1:nrow(samplefile)
samplefile$samplename <- paste0("(anonymous_sample_", samplefile$index, ")")
samplefile <- samplefile[,-c("missing","sex","ID_2")]
setnames(samplefile, "ID_1","eid")
i <- which(samplefile$eid %in% gli_query_data$eid) # 168,104 out of 168,339
gli_query_data_samples <- paste0("(anonymous_sample_", i, ")")

x <- merge(gli_query_data, samplefile, by="eid")
gli_query_data <- x # 168,104
fwrite(gli_query_data, file.path("data","clean","ukbb_spiro_and_geno_qc.csv")
       , quote=F, row.names=F, col.names=T, sep=",")

# 6. Want spirometrically-defined COPD cases as per GOLD 2-4 class definitions and their spirometry measures for association analysis with variants.
# Definition GOLD 2-4 (moderate to very severe lung function): FEV1/FVC < 0.7 and FEV1pp < 0.8
# i.e. hypothesis = SLC26A9 locus variants influence the severity of lung disease in spirometrically-defined COPD patients
# get GOLD 2-4 level lung function (moderate-very severe)
gold24 <- gli_query_data_lowfunc[fev1pp<80] # 14,297

# Inverse rank normal transform (IRNT) the data:
spirodat_irnt <- gold24
columns_of_interest <- c("fev1","fvc","fev1pp")
for(coli in columns_of_interest) {
  colvector <- gold24[[coli]]
  vec <- colvector[[coli]]
  narows <- which(is.na(vec))
  if(length(narows)>0) vec <- vec[-narows]
  vecRankNorm <- rep(NA, nrow(spirodat_irnt))
  vecRankNorm <- rankNorm(vec)
  if(length(narows)>0) vecRankNorm[-narows] <- rankNorm(vec)
  spirodat_irnt[[coli]] <- vecRankNorm
  newcolname <- paste0(coli,".irnt")
  colnames(spirodat_irnt)[which(colnames(spirodat_irnt) %in% coli)] <- newcolname
}

spirodat_irnt <- cbind(spirodat_irnt
                       ,fev1=gold24$fev1
                       ,fvc=gold24$fvc
                       ,fev1pp=gold24$fev1pp
)


# Next, determine PCA subsets for each phenotype that will be analyzed
samplefile <- fread("/hpf/largeprojects/struglis/datasets/uk_biobank_40946/imputation/sample_files/ukb40946_imp_chr1_v3_s487324.sample")
samplefile <- samplefile[-1,] # 487,409
samplefile$index <- 1:nrow(samplefile)
samplefile$samplename <- paste0("(anonymous_sample_", samplefile$index, ")")
samplefile <- samplefile[,-c("missing","sex","ID_2")]
setnames(samplefile, "ID_1","eid")
i <- which(samplefile$ID_1 %in% spirodat_irnt$eid) # 14,281 (out of 14,297)
copd_samples <- paste0("(anonymous_sample_", i, ")")

x <- merge(spirodat_irnt, samplefile, by="eid")
spirodat_irnt <- x

fwrite(spirodat_irnt
       ,file.path("data","intermediate_files","GOLD2-4_copd_ukbb_spirodata.csv")
       ,quote=F, row.names=F, col.names=T, sep=","
       )
