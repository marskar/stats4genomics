library(GEOquery)
#source("http://bioconductor.org/biocLite.R"); biocLite("lumi")
#library(lumi)
library(limma)
#source("https://bioconductor.org/biocLite.R"); biocLite("illuminaHumanv4.db")
library(illuminaHumanv4.db)
#source("https://bioconductor.org/biocLite.R"); biocLite("massiR")
library(massiR)
library(dplyr)

# download the T cell files: GSE56580
files <- getGEOSuppFiles('GSE56580')
untar('GSE56580/GSE56580_RAW.tar', exdir="./GSE56580")
# gunzip file : gunzip GSE56580_non_normalized.txt.gz

## read expression data 
# this works
test<-read.ilmn("./GSE56580/GSE56580_non_normalized.txt.gz", probeid= "ID_REF", expr='intensity')

## background correction and quantile normalization
test.norm <- neqc(test)
# select probe_ids
ids  <- as.character(rownames(test.norm))
#ids2 <- unlist(mget(ids,revmap(illuminaHumanv4ARRAYADDRESS),ifnotfound = NA )) #unnecessary
qual <- unlist(mget(ids, illuminaHumanv4PROBEQUALITY, ifnotfound = NA ))
# table quiality of probes
table(qual)
# get average signal from the $E = expression data
AveSignal = rowMeans(test.norm$E)

# plot average expression by quailty 
pdf('Plot1.pdf')
boxplot(AveSignal~qual) #plot average signals
dev.off()

# remove probes based on quality
rem <- qual=="No match "|qual=="Bad"
test.norm.filt <- test.norm[!rem ,]
# compare 
dim(test.norm )
dim(test.norm.filt )


# so we remove
nrow(test.norm) - nrow(test.norm.filt) # [1] 12375 probes 

## load Y probles from MassiR
data(y.probes)

# illumina_humanht_12

# extract the correct chip
illyprobe <- data.frame(y.probes["illumina_humanht_12"])

head(y.probes)

# extract y probes from eset 

temp <- as.data.frame(test.norm.filt$E)
y.eset <- massi_y(temp, illyprobe)
class(temp)
class(y.eset)
# how many y probes?
length(y.eset$id) # [1] 92

# Expression variation plot
pdf('test.Plot2.pdf')
massi_y_plot(y.eset)
dev.off()

gset <- getGEO("GSE56580", GSEMatrix =TRUE, AnnotGPL=TRUE)

if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno <- Biobase::pData(gset)
edata <- Biobase::exprs(gset)
class(gset)
# selecting only probes with CV on the top 25%
topprobes <- massi_select(gset, illyprobe, threshold=4)
# how many probes?
nrow(topprobes) #[1] 28
# look at first 5
head(topprobes)[,1:5]

## predict sex
results <- massi_cluster(topprobes)

s.results <- data.frame(results[[2]])

head(s.results)

table(s.results$sex)

meta <- pData(gset) # ok this is a data.frame

# subset meta select: 
pheno <- meta[,c(2,12,14:17)]

# add predicted sex

temp <- s.results[,c(1,5)]
temp <- merge(pheno,temp, by.x='geo_accession',by.y='ID')

pheno <- temp
head(pheno)

# fix variables
colnames(pheno)[2] <- 'Age_years'
#pheno$Age_years <- pheno$characteristics_ch1.2
pheno$Age_years <- gsub('.*:', '', pheno$Age_years)
pheno$Age_years <- as.numeric(as.character(pheno$Age_years))

summary(pheno$Age_years)
# plot Age 

pdf('plot3.pdf')
hist(pheno$Age_years)
dev.off()
