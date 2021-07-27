library(data.table)
library(RNOmni)

# Load data
dat <- fread("ukb24727_copd_and_spirometry2.tab", header=T, stringsAsFactors=F) # 502,543

# Note: column numbers 1-based; add +1 to ukb24727.html table (e.g. column 8236 is 8235 in ukb24727.html table)
# temp <- subset(dat, dat$f.22130.0.0 %in% 1) # 1,768 COPD cases

# 1. Remove failed QC individuals (column 8236, UDI=22010-0.0; poor heterozygosity/missingness)
dat1 <- subset(dat, !(dat$f.22010.0.0 %in% 1)) # 502,063; 480 removed
dat <- dat1

# 2. Remove individuals with no genetic sex information
dat2 <- subset(dat, !is.na(dat$f.22001.0.0)) # 487,826; 14,237 removed
dat <- dat2

# 3. Remove related individuals
relinds <- fread("set_of_related_ind_to_rm.txt", header=F, stringsAsFactors=F) # 36,100
# check all fid==iid
# for(i in 1:nrow(relinds)) {
#   if(relinds[i,1] != relinds[i,2]) print(relinds[i,])
# }
#V1      V2
#1: VN061 HG02061 --> not in the dataset
relinds <- relinds[V2!="HG02061"]
relinds$V1 <- as.integer(relinds$V1)
relinds$V2 <- as.integer(relinds$V2)
dat3 <- subset(dat, !(dat$f.eid %in% relinds$V2)) # 451,796; 36,030 removed
dat <- dat3

# 4. Remove non-Caucasian (column 8193; 22006-0.0)
dat4 <- subset(dat, dat$f.22006.0.0==1) # 377,668; 74,128 removed (all Caucasians in column were present)
dat <- dat4

# 5. Get doctor-diagnosed COPD cases/controls
dat5 <- subset(dat, !is.na(dat$f.22130.0.0)) # 95,585 (94,166 controls (0), 1419 cases (1))
dat <- dat5

# 6. Want COPD cases only and their spirometry measures for association analysis with variants:
# i.e. hypothesis = SLC26A9 locus variants influence the severity of lung disease in COPD patients
dat6 <- subset(dat, dat$f.22130.0.0 %in% 1) # 1,419
dat <- dat6

# 7. Now get a dataframe with 
#   FEV1 (3063-0.0 to -2.2), 
#   FEV1 best measure (20150-0.0)
#   predicted percentage (20154-0.0)
#   PEF (3064-0.0 to -2.2)
#   FVC (3062-0.0 to -2.2)
#   FVC  best measure (20151-0.0)

spiro_cols <- c( colnames(dat)[grep("306",colnames(dat))], colnames(dat)[grep("2015",colnames(dat))] )
desired_cols <- c("f.eid", spiro_cols)
spirodat <- dat[, ..desired_cols]
# We want the best measure for the *-0.* data point (highest N--this is the main recruitment phase/instance)
desired_cols <- desired_cols[c(1,grep("0.0", desired_cols))]
spirodat <- dat[, ..desired_cols]
#fev_fvc <- dat[,"f.3063.0.0"] / dat[,"f.3062.0.0"]
fev_fvc_best <- dat[,"f.20150.0.0"] / dat[,"f.20151.0.0"] # best fev1 / best fvc
#colnames(fev_fvc) <- "fev_fvc"
colnames(fev_fvc_best) <- "fev_fvc_best"
#fev_fvc <- fev_fvc$fev_fvc
fev_fvc_best <- fev_fvc_best$fev_fvc_best
spirodat <- cbind(spirodat, fev_fvc_best)

# Inverse rank normal transform (IRNT) the data:
spirodat_irnt <- data.table(matrix(ncol=ncol(spirodat),nrow=nrow(spirodat)))
spirodat_irnt[,1] <- spirodat[,"f.eid"]
for(coli in 2:ncol(spirodat)) {
  colvector <- spirodat[,..coli]
  narows <- which(is.na(colvector))
  vec <- c(colvector)[[1]][-narows]
  vecRankNorm <- rep(NA, nrow(spirodat))
  vecRankNorm[-narows] <- rankNorm(vec)
  spirodat_irnt[,coli] <- vecRankNorm
}
colnames(spirodat_irnt) <- colnames(spirodat)

# add genetic sex (column 8188. UID=22001-0.0. 223,481 males (1) and 264,814 females (0)) and age at recruitment (column 8149. UID=21022-0.0):
spirodat_irnt <- cbind(spirodat_irnt, dat[,c("f.22001.0.0","f.21022.0.0")]) # 568 females (0) and 851 males (1)



# Next, determine PCA subsets for each phenotype that will be analyzed
samplefile <- fread("/hpf/largeprojects/struglis/datasets/uk_biobank_40946/imputation/sample_files/ukb40946_imp_chr1_v3_s487324.sample")
samplefile <- samplefile[-1,]
i <- which(samplefile$ID_1 %in% spirodat_irnt$f.eid)
copd_samples <- paste0("(anonymous_sample_", i, ")") # 1,417 samples
covar <- spirodat_irnt[,c("f.eid","f.22001.0.0","f.21022.0.0")]
colnames(covar)<-c("f.eid","sex","age")
covar <- cbind(covar, age2=covar$age^2)

spiropass <- spirodat[,c("f.eid", "f.20152.0.0")]
# get "best" PEF:
bestpef <- NULL
desired_cols <- c("f.eid", colnames(dat)[grep("3064.0",colnames(dat))] )
pefdat <- dat[,..desired_cols]
pefdat <- cbind(pefdat, miss=apply(pefdat,1,function(x) sum(is.na(x))))
pefdat_nmiss <- subset(pefdat, pefdat$miss != 3) # 1125
x <- merge(pefdat_nmiss, spiropass, by="f.eid")
pefdat_pass <- subset(x, x$f.20152.0.0 %in% 1) # 848
# select best PEF according to sum(best fev + best fvc)
desired_cols <- c("f.eid", colnames(dat)[grep("3062.0",colnames(dat))], colnames(dat)[grep("3063.0",colnames(dat))] )
fevfvcdat <- dat[,..desired_cols]
fevfvcdat <- cbind(fevfvcdat, miss=apply(fevfvcdat,1,function(x) sum(is.na(x))))
fevfvcdat <- subset(fevfvcdat, fevfvcdat$miss != 6) # 1125
sum0 <- fevfvcdat$f.3062.0.0 + fevfvcdat$f.3063.0.0
sum1 <- fevfvcdat$f.3062.0.1 + fevfvcdat$f.3063.0.1
sum2 <- fevfvcdat$f.3062.0.2 + fevfvcdat$f.3063.0.2
blowsums <- cbind(sum0,sum1,sum2)
bestpefblow <- unlist(apply(blowsums, 1, function(x) which(x %in% max(x,na.rm=T))[1]))
fevfvcdat <- cbind(fevfvcdat, bestpefblow)
pefdat_pass <- merge(pefdat_pass, fevfvcdat[,c("f.eid","bestpefblow")], by="f.eid")
bestpef <- apply(pefdat_pass, 1, function(x) x[x['bestpefblow']+1])
pefdat_pass <- cbind(pefdat_pass, bestpef)

spirodat <- merge(spirodat, pefdat_pass[,c("f.eid","bestpef")], by="f.eid", all.x=T)
coli <- which(colnames(spirodat) %in% "bestpef")
colvector <- spirodat[,..coli]
narows <- which(is.na(colvector))
vec <- c(colvector)[[1]][-narows]
vecRankNorm <- rep(NA, nrow(spirodat))
vecRankNorm[-narows] <- rankNorm(vec)
  
spirodat_irnt <- cbind(spirodat_irnt, bestpef=vecRankNorm)


phenos_to_analyze <- c("f.20150.0.0", "f.20151.0.0", "fev_fvc_best", "bestpef", "f.20154.0.0") # best FEV1, best FVC, FEV1/FVC best, Best PEF, FEV1pp
for(j in 1:length(phenos_to_analyze)) {
  phenoj <- spirodat_irnt[,c("f.eid", phenos_to_analyze[j]),with=F]
  colnames(phenoj) <- c("f.eid","pheno")
  phenoname <- phenos_to_analyze[j]
  phenoname <- gsub("f.(\\d+).0.0","\\1",phenoname)
  phenoname <- paste0(phenoname,".irnt.copd")
  print(paste0("Phenotype: ", phenoname))
  x <- merge(phenoj, covar, by="f.eid")
  N <- nrow(na.omit(x))
  print(paste0("NMISS: ", N))
}

# Save the IRNT spirometry measures for association analysis:
cols_to_save <- c("f.eid", phenos_to_analyze, "f.22001.0.0","f.21022.0.0") # 22001, 21022 = sex, age; 
fwrite(spirodat_irnt[,..cols_to_save], "ukbb24727_spirometry_in_copd_irnt2.tsv", quote=F, row.names=F, col.names=T, sep="\t")



# PCA subset for FEV1, FVC, and FEV1/FVC ratio (N=1008); 7 not genotyped, so 1001
phenoj <- spirodat_irnt[,c('f.eid', 'fev_fvc_best')]
phenoj <- na.omit(phenoj)
fev_set <- cbind(phenoj$f.eid, phenoj$f.eid)
fwrite(fev_set, "pca/fev_copd_set.txt", quote=F, row.names=F, col.names=F, sep="\t")

# PCA subset for best PEF (N=848); 6 not genotyped, so 842
phenoj <- spirodat_irnt[,c('f.eid', 'bestpef')]
phenoj <- na.omit(phenoj)
bestpef_set <- cbind(phenoj$f.eid, phenoj$f.eid)
fwrite(bestpef_set, "pca/bestpef_set.txt", quote=F, row.names=F, col.names=F, sep="\t")

# PCA subset for FEV1pp (N=543); 6 not genotyped, so 537
phenoj <- spirodat_irnt[,c('f.eid', 'f.20154.0.0')]
phenoj <- na.omit(phenoj)
fev1pp_set <- cbind(phenoj$f.eid, phenoj$f.eid)
fwrite(fev1pp_set, "pca/fev1pp_set.txt", quote=F, row.names=F, col.names=F, sep="\t")


