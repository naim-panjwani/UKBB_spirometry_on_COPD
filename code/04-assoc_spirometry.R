# Association analysis of inverse rank normal transformed (IRNT)
# spirometry measured in GOLD2-4 spirometry-defined COPD patients in the UKBB around SLC26A9 locus

library(data.table)
library(rbgen)
library(RNOmni)
library(bigmemory)

get_irnt <- function(vec) {
  narows <- which(is.na(vec))
  if(length(narows)>0) vec <- vec[-narows]
  vecRankNorm <- rep(NA, length(vec))
  vecRankNorm <- rankNorm(vec)
  if(length(narows)>0) vecRankNorm[-narows] <- rankNorm(vec)
  return(vecRankNorm)
}

calculate_dr2 <- function(y1,y2) {
  if(length(y1)!=length(y2)) stop("vector lengths differ")
  u <- y1 + 2*y2
  z <- u
  w <- y1 + 4*y2
  n <- length(y1)
  numerator <- ( sum(z*u) - (1/n)*(sum(u)*sum(z)) )^2
  denominator <- (sum(w) - (1/n)*(sum(u))^2) * (sum(z^2) - (1/n)*(sum(z))^2)
  r2 <- numerator / denominator
  return(r2)
}

convert_to_eid <- function(anonymous_sample_name) {
  if(length(anonymous_sample_name)!=1 & class(anonymous_sample_name)!="character") {
    stop("input must be of length 1 and character type")
  }
  return(pheno$eid[which(pheno$samplename %in% anonymous_sample_name)])
}

get_dosagemat <- function(bgen_filename, ranges, samples, chunksize, fileout) {
  rangestr <- with(ranges, paste0(as.character(chromosome),":",as.character(start),"-",as.character(end)))
  print(paste0("Getting number of SNPs in ", rangestr))
  tmp <- bgen.load(bgen_filename, ranges, samples = samples[1:10])
  mapdf <- data.table(tmp$variants)
  fwrite(mapdf, paste0(fileout,".map"), sep="\t", quote=F)
  numsnps <- dim(tmp$data)[1]
  numsamples <- length(samples)
  print(paste0("Initializing big matrix with dimension: ",numsamples," x ",numsnps))
  dat <- big.matrix(nrow=numsamples, ncol=numsnps, type="double", init=-9.0)
  #dat <- data.table(matrix(rep(-9.0,numsnps*numsamples), nrow=numsnps, ncol=numsamples))
  
  samplecounter <- 1
  while(samplecounter <= numsamples) {
    starti <- samplecounter
    endi <- starti + chunksize - 1
    if(endi>numsamples) endi <- numsamples
    ivcfmat <- matrix(rep("",(endi-starti+1)*numsnps), nrow=endi-starti+1, ncol=numsnps) # inverted vcf data contents (sample x snp)
    print(paste0("Getting samples ",starti," to ",endi, " out of ", numsamples))
    currdata <- bgen.load(bgen_filename, ranges, samples = samples[starti:endi])
    currsamplenames <- dimnames(currdata$data)[[2]]
    currsamplenames <- unname(sapply(currsamplenames, convert_to_eid))
    print(paste0("Writing sample eid's to: ", paste0(fileout,".sample")))
    fwrite(data.table(eid=currsamplenames), paste0(fileout,".sample"), sep="\t", quote=F, append=TRUE)
    
    for(j in 1:numsnps){
      if(j==1 | (j %% 100)==0) {
        print(paste0("SNP ",j," of ",numsnps,". "
                     ,format(round((j/numsnps)*100,2),nsmall=2)
                     ,"% done"))
      }
      snpi <- currdata$data[j,,]
      dosage <- round(snpi[,"g=1"] + 2*snpi[,"g=2"], 3)
      dat[starti:endi,j] <- dosage
      
      maxGP <- unname(unlist(apply(snpi, 1, function(x) which(x==max(x))[1]-1)))
      names(maxGP) <- rownames(snpi)
      GT <- ifelse(maxGP==0, "0/0", ifelse(maxGP==1, "0/1", ifelse(maxGP==2, "1/1", NA)))
      DS <- dosage
      GP <- apply(snpi, 1, function(x) paste0(round(x,3),collapse=","))
      ivcfmat[,j] <- paste0(GT,":",DS,":",GP)
      #currsnpnames <- dimnames(currdata$data)[[1]]
      #colnames(ivcfmat) <- currsnpnames
      #rownames(ivcfmat) <- currsamplenames
    }
    
    print(paste0("Writing imputation data to: ", paste0(fileout,".dat")))
    fwrite(ivcfmat, paste0(fileout,".dat"), sep="\t", quote=F, append=TRUE)
    
    samplecounter <- endi+1
  }
  return(dat)
}
fileout <- "data/intermediate_files/ukbb_imputation_slc26a9/ukbb_spiroqc"
outdir <- "data/clean/assoc/"

print("Loading spirometry phenotypes")
pheno <- fread("data/clean/ukbb_spiro_and_geno_qc.csv") # 168,104
pheno[,ratio := fev1/fvc]
pheno[,eid := as.character(eid)]

print("Loading PCA")
pca <- fread("data/intermediate_files/pca/18-ukbb_ukbbspiro_flashpca2_eigenvectors.txt", header=F, stringsAsFactors = F) # 167,655
pca <- pca[,c("V1","V2","V3","V4"),with=F]  # adjust with 3 PCs
pca[,V1 := gsub("(.*):(.*)","\\2",V1)]
colnames(pca) <- c("eid","PC1","PC2","PC3")

covar <- pheno[,c("eid","sex","age")]
covar <- cbind(covar, age2=covar$age^2)
covar[,eid := as.character(eid)]

print("Merging phenotypes and PCA")
pheno <- merge(pheno,pca,by="eid") # 167,655 - some eid's (449) were missing from genotype arrays
pheno[,ratio.irnt := get_irnt(ratio)]
pheno[,fev1pp.irnt := get_irnt(fev1pp)]


print(paste0("Loading imputation data for ",nrow(pheno)," samples"))
#data = bgen.load( "/hpf/largeprojects/struglis/datasets/uk_biobank_40946/imputation/ukb_imp_chr1_v3.bgen", ranges, samples = pheno$samplename) # too big to load
bgen_filename <- "/hpf/largeprojects/struglis/datasets/uk_biobank_40946/imputation/ukb_imp_chr1_v3.bgen"
ranges = data.frame(
  chromosome = "01",
  start = 205780000,
  end = 205940000
)
samples <- pheno$samplename
chunksize <- 15000
dosage <- get_dosagemat(bgen_filename, ranges, samples, chunksize, fileout)
write.big.matrix(dosage, paste0(fileout,"_dosage.bigmat"), sep="\t")

#dosage <- read.big.matrix(paste0(fileout,"_dosage.bigmat"), sep="\t", type="double")
map <- fread(paste0(fileout,".map")) # generated above with get_dosagemat function
sample <- fread(paste0(fileout,".sample")) # generated above with get_dosagemat function


# Association with all 167,655 individuals
print(paste0("Starting association analysis for ",nrow(pheno)," samples and ",nrow(map)," variants"))
phenocols <- c("ratio.irnt","fev1pp.irnt", "hasCOPD")
results <- NULL
for(j in 1:length(phenocols)) {
  phenoj <- pheno[,c("eid", phenocols[j]),with=F]
  phenoname <- phenocols[j]
  print(paste0("Running association analysis for ", phenoname))
  phenoj_assoc <- NULL
  # for each snp
  for(i in 1:ncol(dosage)) {
    if(i==1 | (i %% 100)==0) {
      print(paste0("SNP ",i," of ",ncol(dosage),". "
                   ,format(round((i/ncol(dosage))*100,2),nsmall=2), "% done"))
    }
    dosagei <- dosage[,i]
    af <- sum(dosagei,na.rm=T) / (2*length(dosagei))
    maf <- af
    if(maf>0.5) maf <- 1-maf
    if(maf<(10/length(dosagei))) next
    snpi <- data.table(cbind(sample, dosagei))
    snpi[,eid := as.character(eid)]
    
    x <- merge(phenoj, covar, by="eid")
    x <- merge(x, pca, by="eid")
    x <- merge(x, snpi[,c("eid","dosagei")], by="eid")
    #x <- x[,-c("eid"),with=F]
    x <- na.omit(x)
    x[,dosagei := as.numeric(dosagei)]
    x[,sex := as.factor(sex)]
    realmaf <- sum(x$dosagei,na.rm=T)/(2*nrow(x))
    if(realmaf>0.5) maf <- 1-maf
    if(realmaf<(10/nrow(x))) next
    pcol <- 'Pr(>|t|)'
    if(phenoname %in% "hasCOPD") {
      model <- glm(formula(paste0(phenoname," ~ dosagei + PC1 + PC2 + PC3 + sex + age + age2 + sex:age + sex:age2")), family="binomial", data=x)
      pcol <- 'Pr(>|z|)'
    } else {
      model <- lm(formula(paste0(phenoname," ~ dosagei + PC1 + PC2 + PC3 + sex + age + age2 + sex:age + sex:age2")), data=x)
    }
    snp <- map[i,]
    coefmat <- coef(summary(model))
    assoc <- data.frame(chrom=as.integer(snp$chromosome)
                        ,pos=snp$position
                        ,rsid=as.character(snp$rsid)
                        ,allele0=snp$allele0
                        ,allele1=snp$allele1
                        ,beta=coefmat['dosagei','Estimate']
                        ,std.err=coefmat['dosagei','Std. Error']
                        ,P=coefmat['dosagei',pcol]
                        ,N=nrow(x)
                        ,allele_freq=af
                        ,phenotype=phenoname
                        )
    phenoj_assoc <- rbind(phenoj_assoc, assoc)
    results <- rbind(results, assoc)
  }
  outname <- paste0(outdir,phenoname,".assoc.csv")
  print(paste0("Saving results to: ", outname))
  fwrite(phenoj_assoc, outname, quote=F, row.names=F, col.names=T)
}
outname <- paste0(outdir,"spirometry_assoc_ukbb_26a9.csv")
print(paste0("Saving all spirometry association results to: ", outname))
fwrite(results, outname, quote=F, row.names=F, col.names=T)



# Association with GOLD2-4 COPD individuals
print("")
print("Starting association for GOLD2-4 COPD subset")
pheno <- pheno[,-c("ratio.irnt","fev1pp.irnt","PC1","PC2","PC3"),with=F]
pheno <- pheno[hasCOPD==TRUE]

print("Loading PCA")
pca <- fread("data/intermediate_files/pca/15-ukbb_copd_pcair_eigenvectors.txt", header=F, stringsAsFactors = F)
pca <- pca[,c("V1","V2","V3","V4"),with=F]  # adjust with 3 PCs
pca[,V1 := gsub("(.*):(.*)","\\2",V1)]
colnames(pca) <- c("eid","PC1","PC2","PC3")

print("Merging phenotypes with PCA")
pheno <- merge(pheno,pca,by="eid")
pheno[,ratio.irnt := get_irnt(ratio)]
pheno[,fev1pp.irnt := get_irnt(fev1pp)]
ranges = data.frame(
  chromosome = "01",
  start = 205780000,
  end = 205940000
)
print(paste0("Loading imputation data for ",nrow(pheno)," samples"))
data = bgen.load( "/hpf/largeprojects/struglis/datasets/uk_biobank_40946/imputation/ukb_imp_chr1_v3.bgen", ranges, samples = pheno$samplename)


covar <- pheno[,c("eid","sex","age")]
covar <- cbind(covar, age2=covar$age^2)


# Dosage sample file:
print("Loading imputation sample file")
samplefile <- fread("/hpf/largeprojects/struglis/datasets/uk_biobank_40946/imputation/sample_files/ukb40946_imp_chr1_v3_s487324.sample")
samplefile <- samplefile[-1,] # 487,409
samplefile$index <- 1:nrow(samplefile)
samplefile$samplename <- paste0("(anonymous_sample_", samplefile$index, ")")
samplefile <- samplefile[,-c("missing","sex","ID_2")]
setnames(samplefile, "ID_1","eid")
#i <- which(samplefile$eid %in% pheno$eid) # 14,196
copd_samples <- paste0("(anonymous_sample_", i, ")")


print(paste0("Starting association analysis for ",nrow(pheno)," samples and ",nrow(data$variants)," variants"))
phenocols <- c("ratio.irnt","fev1pp.irnt")
results <- NULL
for(j in 1:length(phenocols)) {
  phenoj <- pheno[,c("eid", phenocols[j]),with=F]
  phenoname <- phenocols[j]
  print(paste0("Running association analysis for ", phenoname))
  phenoj_assoc <- NULL
  # for each snp
  for(i in 1:nrow(data$variants)) {
    if(i==1 | (i %% 100)==0) {
      print(paste0("SNP ",i," of ",nrow(data$variants),". "
                   ,format(round((i/nrow(data$variants))*100,2),nsmall=2), "% done"))
    }
    snpi <- data$data[i,,]
    dosage <- snpi[,"g=1"] + 2*snpi[,"g=2"]
    af <- sum(dosage,na.rm=T) / (2*length(dosage))
    maf <- af
    if(maf>0.5) maf <- 1-maf
    if(maf<(10/length(dosage))) next
    snpi <- cbind(snpi,dosage)
    sample_indices <- as.integer(gsub("\\(anonymous_sample_(\\d+)\\)","\\1",names(dosage)))
    sample_ids <- as.character(samplefile$eid[sample_indices])
    snpi <- data.table(cbind(snpi, sample_ids))
    covar[,eid := as.character(eid)]
    x <- merge(phenoj, covar, by="eid")
    x <- merge(x, pca, by="eid")
    x <- merge(x, snpi[,c("sample_ids","dosage")], by.x="eid", by.y="sample_ids")
    x <- x[,-c("eid"),with=F]
    x <- na.omit(x)
    x[,dosage := as.numeric(dosage)]
    x[,sex := as.factor(sex)]
    realmaf <- sum(x$dosage,na.rm=T)/(2*nrow(x))
    if(realmaf>0.5) maf <- 1-maf
    if(realmaf<(10/nrow(x))) next
    model <- glm(formula(paste0(phenoname,"~ dosage + PC1 + PC2 + PC3 + sex + age + age2 + sex:age + sex:age2")), data=x)
    snp <- data$variants[i,]
    coefmat <- coef(summary(model))
    assoc <- data.frame(chrom=as.integer(snp$chromosome)
                        ,pos=snp$position
                        ,rsid=as.character(snp$rsid)
                        ,REF=snp$allele0
                        ,ALT=snp$allele1
                        ,beta=coefmat['dosage','Estimate']
                        ,std.err=coefmat['dosage','Std. Error']
                        ,P=coefmat['dosage','Pr(>|t|)']
                        ,N=nrow(x)
                        ,allele_freq=af
                        ,phenotype=phenoname
                        )
    phenoj_assoc <- rbind(phenoj_assoc, assoc)
    results <- rbind(results, assoc)
  }
  outname <- paste0(outdir,phenoname,".assoc_copd_only.csv")
  print(paste0("Saving results to: ", outname))
  fwrite(phenoj_assoc, outname, quote=F, row.names=F, col.names=T)
}
outname <- paste0(outdir,"spirometry_assoc_ukbb_copd_only_26a9.csv")
print(paste0("Saving all results of COPD association to: ", outname))
fwrite(results, outname, quote=F, row.names=F, col.names=T)

