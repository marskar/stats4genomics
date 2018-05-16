library(GEOquery)
#source("http://bioconductor.org/biocLite.R"); biocLite("lumi")
#library(lumi)
library(limma)
#source("http://bioconductor.org/biocLite.R"); biocLite("samr")
library(samr)
#source("https://bioconductor.org/biocLite.R"); biocLite("illuminaHumanv4.db")
library(illuminaHumanv4.db)
#source("https://bioconductor.org/biocLite.R"); biocLite("massiR")
library(massiR)
#install.packages("DT")
library(DT)

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
head(gset)
pheno2 <- Biobase::pData(gset)
edata <- Biobase::exprs(gset)

dim(gset)
slotNames(gset)
dim(phenoData(gset))
head(pheno)
head(pheno2)
head(edata)


head(pheno)
head(edata)
enorm <- normalizeQuantiles(edata)

genenames <- rownames(edata)

data <- list(x = enorm, y = pheno$Age_years, geneid = genenames, genenames = genenames, logged2 = TRUE)
dim(enorm)
head(data)
length(pheno$Age_years)
### Create samr object
samr.obj <- samr(data, resp.type = "Quantitative", nperms=100)
### Look at structure of the samr object
str(samr.obj)

delta.table <- samr.compute.delta.table(samr.obj, min.foldchange = 1.5)  # Compute thresholds for different deltas
datatable(delta.table)  # Look at the whole range

delta.table[delta.table[, "median FDR"] < 0.1, ][1, ]  # Check delta corresponding to median FDR ~0.1

delta <- 0.49
samr.plot(samr.obj, delta)

hist(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar))
head(samr.obj$tt)
siggenes.table <- samr.compute.siggenes.table(samr.obj, delta, data, delta.table, min.foldchange = 1.5)  # Summarize significant genes
names(siggenes.table)  # What data we have in the summary list

nrow(siggenes.table$genes.up)  # How many upregulated genes

nrow(siggenes.table$genes.lo)  # How many downregulated

datatable(siggenes.table$genes.up)  # Check how table with the results look like

# model 1: age as continuos:
dmat1 <- model.matrix(~pheno$Age_years)
colnames(dmat1) <- c('Intercept', 'Age')

# contrast matrix to set comparisons:
# cmat <- makeContrasts(levels=colnames(dmat),Age_yearsvsCTRL=(Age_years - control))

# Fit the model 1
fit.ls <- lmFit(test.norm.filt,dmat1,method="ls")

#fit.ls <- contrasts.fit(fit.ls,cmat)

# moderation of standard erros using empirical Bayes
eb.ls <- eBayes(fit.ls)

## plot p values model1

pdf('plot4.pvalueMod1.pdf')
hist(eb.ls$p.value[,1], main='')
dev.off()
