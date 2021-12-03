# Association analysis of inverse rank normal transformed (IRNT)
# spirometry measured in COPD patients in the UKBB around SLC26A9 locus

library(data.table)
library(rbgen)

samples_not_in_pca <- c("1285771", "1424258", "1605913", "2133505", "4062544", "4664239", "5591128") # these are not genotyped

pheno <- fread("ukbb24727_spirometry_in_copd_irnt2.tsv", header=T, stringsAsFactors = F, sep="\t")
pheno <- pheno[!(f.eid %in% samples_not_in_pca)]
samplefile <- fread("/hpf/largeprojects/struglis/datasets/uk_biobank_40946/imputation/sample_files/ukb40946_imp_chr1_v3_s487324.sample")
samplefile <- samplefile[-1,]
i <- which(samplefile$ID_1 %in% pheno$f.eid)
copd_samples <- paste0("(anonymous_sample_", i, ")") # 1,410 samples
ranges = data.frame(
  chromosome = "01",
  start = 205780000,
  end = 205940000
)
data = bgen.load( "/hpf/largeprojects/struglis/datasets/uk_biobank_40946/imputation/ukb_imp_chr1_v3.bgen", ranges, samples = copd_samples)
#pca <- fread("pca/15-ukbb_copd_pcair_eigenvectors.txt", header=F, stringsAsFactors = F)
#pca <- pca[,c("V1","V2","V3"),with=F]  # adjust with 2 PCs
#colnames(pca) <- c("f.eid","PC1","PC2")

covar <- pheno[,c("f.eid","f.22001.0.0","f.21022.0.0")]
colnames(covar)<-c("f.eid","sex","age")
covar <- cbind(covar, age2=covar$age^2)

#phenocols <- c("f.3062.0.0", "f.3063.0.0", "f.3064.0.0", "f.20150.0.0", "f.20151.0.0", "f.20154.0.0", "fev_fvc", "fev_fvc_best")
phenocols <- c("f.20150.0.0", "f.20151.0.0", "f.20154.0.0", "fev_fvc_best", "bestpef")
results <- NULL
for(j in 1:length(phenocols)) {
  phenoj <- pheno[,c("f.eid", phenocols[j]),with=F]
  colnames(phenoj) <- c("f.eid","pheno")
  phenoname <- phenocols[j]
  phenoname <- gsub("f.(\\d+).0.0","\\1",phenoname)
  phenoname <- paste0(phenoname,".irnt.copd")
  print(paste0("Association analysis for ", phenoname))
  pcafile <- ifelse(phenocols[j] %in% c("f.20150.0.0", "f.20151.0.0", "fev_fvc_best"), "15-ukbb_copd_fevset_pcair_eigenvectors.txt",
                    ifelse(phenocols[j] %in% "bestpef", "15-ukbb_copd_bestpefset_pcair_eigenvectors.txt", "15-ukbb_copd_fev1ppset_pcair_eigenvectors.txt"))
  pcafile <- paste0("pca/",pcafile)
  pca <- fread(pcafile, header=F, stringsAsFactors = F)
  pca <- pca[,c("V1","V2","V3"),with=F]  # adjust with 2 PCs (could probably get away with just one PC...)
  colnames(pca) <- c("f.eid","PC1","PC2")
  phenoj_assoc <- NULL
  # for each snp
  for(i in 1:nrow(data$variants)) {
    snpi <- data$data[i,,]
    dosage <- snpi[,"g=1"] + 2*snpi[,"g=2"]
    af <- sum(dosage,na.rm=T) / (2*length(dosage))
    maf <- af
    if(maf>0.5) maf <- 1-maf
    if(maf<(10/length(dosage))) next
    snpi <- cbind(snpi,dosage)
    sample_indices <- as.integer(gsub("\\(anonymous_sample_(\\d+)\\)","\\1",names(dosage)))
    sample_ids <- samplefile$ID_1[sample_indices]
    snpi <- cbind(snpi, sample_ids)
    x <- merge(phenoj, snpi[,c("sample_ids","dosage")], by.x="f.eid", by.y="sample_ids")
    x <- merge(x, covar, by="f.eid")
    x <- merge(x, pca, by="f.eid")
    x <- x[,-c("f.eid"),with=F]
    x <- na.omit(x)
    model <- glm(pheno ~ dosage + PC1 + PC2 + sex + age + age2 + sex*age + sex*age2, data=x)
    snp <- data$variants[i,]
    coefmat <- coef(summary(model))
    assoc <- data.frame(chrom=as.integer(snp$chromosome), pos=snp$position, rsid=as.character(snp$rsid), REF=snp$allele0, ALT=snp$allele1,
                        beta=coefmat['dosage','Estimate'], std.err=coefmat['dosage','Std. Error'], P=coefmat['dosage','Pr(>|t|)'], N=nrow(x),
                        allele_freq=af, phenotype=phenoname)
    phenoj_assoc <- rbind(phenoj_assoc, assoc)
    results <- rbind(results, assoc)
  }
  fwrite(phenoj_assoc, paste0("assoc/",phenoname,".tsv"), quote=F, row.names=F, col.names=T, sep="\t")
}
fwrite(results, "assoc/spirometry_assoc_in_ukbb_copd_around_26a9_2.tsv", quote=F, row.names=F, col.names=T, sep="\t")
