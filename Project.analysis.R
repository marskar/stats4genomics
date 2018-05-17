# Project.analysis.R

library(GEOquery)
#source("http://bioconductor.org/biocLite.R"); biocLite("lumi")
#library(lumi)
library(limma)
#source("https://bioconductor.org/biocLite.R"); biocLite("illuminaHumanv4.db")
library(illuminaHumanv4.db)
#source("https://bioconductor.org/biocLite.R"); biocLite("massiR")
library(massiR)
library(scales)
library(preprocessCore)
library(calibrate)


# download the T cell files: GSE56580
files <- getGEOSuppFiles('GSE56580')
untar('GSE56580/GSE56580_RAW.tar', exdir="./GSE56580")
# gunzip file : gunzip GSE56580_non_normalized.txt.gz

## read expression data 
# this works
test<-read.ilmn("/users/cvalenci/GSE56580/GSE56580_non_normalized.txt", probeid= "ID_REF", expr='intensity')

## background correction and quantile normalization
test.norm <- neqc(test)
# select probe_ids
ids  <- as.character(rownames(test.norm))
#ids2 <- unlist(mget(ids,revmap(illuminaHumanv4ARRAYADDRESS),ifnotfound = NA )) #unnecessary
qual <- unlist(mget(ids,illuminaHumanv4PROBEQUALITY , ifnotfound = NA ))
# table quiality of probes
table(qual)
# get average signal from the $E = expression data
AveSignal = rowMeans(test.norm$E)

# plot average expression by quailty 
pdf('/users/cvalenci/GSE56580/Plot1.pdf')
boxplot(AveSignal~qual) #plot average signals
dev.off()

# remove probes based on quality
rem <- qual=="No match "|qual=="Bad"
test.norm.filt <- test.norm[!rem ,]
# compare 
dim(test.norm )
dim(test.norm.filt)


# so we remove
nrow(test.norm) - nrow(test.norm.filt) # [1] 12375 probes 

## load Y probles from MassiR
data(y.probes)

# illumina_humanht_12

# extract the correct chip
illyprobe <- data.frame(y.probes["illumina_humanht_12"])

# extract y probes from eset 

temp <- as.data.frame(test.norm.filt$E)
y.eset <- massi_y(temp, illyprobe)

# how many y probes?
length(y.eset$id) # [1] 92

# Expression variation plot
pdf('test.Plot2.pdf')
massi_y_plot(y.eset)
dev.off()

# selecting only probles with CV on the top 25%

topprobes <- massi_select(temp, illyprobe, threshold=4)
# how many probes?
nrow(topprobes) #[1] 28 ; 23?
# look at first 5
head(topprobes)[,1:5]

## predict sex

results <- massi_cluster(topprobes)

s.results <- data.frame(results[[2]])

head(s.results)

table(s.results$sex)

## save results :
save(test.norm.filt, s.results, file='/users/cvalenci/GSE56580/Tcell.data.Rdata')


## pheno data
tset <- getGEO("GSE56580", GSEMatrix =TRUE, AnnotGPL=TRUE)

if (length(tset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
tset <- tset[[idx]]

meta <- pData(tset) # ok this is a data.frame

# subset meta select: 
pheno <- meta[,c(2,12,14:17,)]

# add predicted sex

temp <- s.results[,c(1,5)]
temp <- merge(pheno,temp, by.x='geo_accession',by.y='ID')

pheno <- temp

# fix variables
colnames(pheno)[2] <- 'Age_years'
#pheno$Age_years <- pheno$characteristics_ch1.2
pheno$Age_years <- gsub('.*:', '', pheno$Age_years)
pheno$Age_years <- as.numeric(as.character(pheno$Age_years))

# plot Age 

pdf('/users/cvalenci/GSE56580/plot3.pdf')
hist(pheno$Age_years)
dev.off()

## Diff expression analysis
# test.norm.filt : expression normalize data
# pheno: all pheno data including predicted sex



# Prepare models

# model 1: age as continuos:
dmat1 <- model.matrix(~pheno$Age_years)
colnames(dmat1) <- c('Intercept', 'Age')

# contrast matrix to set comparisons:
# cmat <- makeContrasts(levels=colnames(dmat),Age_yearsvsCTRL=(Age_years - control))

# Fit the model 1 : Looking at Age 
fit.ls <- lmFit(test.norm.filt,dmat1,method="ls")

#fit.ls <- contrasts.fit(fit.ls,cmat)

# moderation of standard erros using empirical Bayes
eb.ls <- eBayes(fit.ls)

## plot p values model1

pdf('/users/cvalenci/GSE56580/plot4.pvalueMod1.pdf')
hist(eb.ls$p.value[,2], main='Model1: Age')
dev.off()

##
tg <- topTable(eb.ls,coef=2,number=Inf,adjust="BH",genelist = rownames(eb.ls))

nrow(subset(tg,adj.P.Val < 0.05 & abs(logFC)>1 ))

# Volcano plot
pdf('/users/cvalenci/GSE56580/Vplot.Mod1.pdf')
with(tg, plot(logFC, -log10(adj.P.Val),pch=20,cex=0.4,main="Model 1: Age" ) )
with(subset(tg,adj.P.Val < 0.05 ), points(logFC,-log10(adj.P.Val),pch=20,col='red' ) )
dev.off()

## model 2: Age+sex

# model 2: age as continuos:

dmat2 <- model.matrix(~pheno$Age_years+pheno$sex)
colnames(dmat2) <- c('Intercept', 'Age', 'sex')


# contrast matrix to set comparisons:
# cmat <- makeContrasts(levels=colnames(dmat),Age_yearsvsCTRL=(Age_years - control))

# Fit the model 2
fit.ls <- lmFit(test.norm.filt,dmat2,method="ls")

#fit.ls <- contrasts.fit(fit.ls,cmat)

# moderation of standard erros using empirical Bayes
eb.ls <- eBayes(fit.ls)

## plot p values model1

pdf('/users/cvalenci/GSE56580/plot.pvalueMod2.pdf')
hist(eb.ls$p.value[,2], main='Model 2: Age + Sex')
dev.off()

##
tg <- topTable(eb.ls,coef=2,number=Inf,adjust="BH",genelist = rownames(eb.ls))

nrow(subset(tg,adj.P.Val < 0.05 & abs(logFC)>1 ))

# Volcano plot
pdf('/users/cvalenci/GSE56580/Vplot.Mod2.pdf')
with(tg, plot(logFC, -log10(adj.P.Val),pch=20,cex=0.4,main="Model 2: Age + sex" ) )
with(subset(tg,adj.P.Val < 0.05 ), points(logFC,-log10(adj.P.Val),pch=20,col='red' ) )
dev.off()


## Model 3: Age continuous + sex 
# categorize by the median
pheno$Age_cat <- ifelse(pheno$Age_years >= 58,1,0)

# model 1: age as continuos:
dmat3 <- model.matrix(~pheno$Age_cat+pheno$sex)
colnames(dmat3) <- c('Intercept', 'Age_cat','sex')

# contrast matrix to set comparisons:
# cmat <- makeContrasts(levels=colnames(dmat),Age_yearsvsCTRL=(Age_years - control))

# Fit the model 3
fit.ls <- lmFit(test.norm.filt,dmat3,method="ls")

#fit.ls <- contrasts.fit(fit.ls,cmat)

# moderation of standard erros using empirical Bayes
eb.ls <- eBayes(fit.ls)

## plot p values model1

pdf('/users/cvalenci/GSE56580/plot.pvalueMod3.pdf')
hist(eb.ls$p.value[,2], main='')
dev.off()

##
tg <- topTable(eb.ls,coef=2,number=Inf,adjust="BH",genelist = rownames(eb.ls))

nrow(subset(tg,adj.P.Val < 0.05 & abs(logFC)>1 ))

# Volcano plot
pdf('/users/cvalenci/GSE56580/Vplot.Mod3.pdf')
with(tg, plot(logFC, -log10(adj.P.Val),pch=20,cex=0.4,main="Model 3: Median Age + Sex" ) )
with(subset(tg,adj.P.Val < 0.05 ), points(logFC,-log10(adj.P.Val),pch=20,col='red' ) )
dev.off()


#################################################################################################
# this part below was not run
##
# model 4: age_cat + sex. Now I care about sex coefficient
# sex will be 1= male vs 0 = female
dmat1 <- model.matrix(~pheno$Age_cat+pheno$sex)
colnames(dmat1) <- c('Intercept', 'Age_cat', 'sex')

# contrast matrix to set comparisons:
# cmat <- makeContrasts(levels=colnames(dmat),Age_yearsvsCTRL=(Age_years - control))

# Fit the model
fit.ls <- lmFit(test.norm.filt,dmat1,method="ls")

#fit.ls <- contrasts.fit(fit.ls,cmat)

# moderation of standard erros using empirical Bayes
eb.ls <- eBayes(fit.ls)

## plot p values model this time lookig at diff gene expression by sex

pdf('/users/cvalenci/GSE56580/plot4.pvalueMod1.pdf')
hist(eb.ls$p.value[,3], main='')
dev.off()
 
## selecting values for sex
tg <- topTable(eb.ls,coef=3,number=Inf,adjust="BH",genelist = rownames(eb.ls))

## anotate genes
x <- illuminaHumanv4SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
annotation <- data.frame(names(xx), unlist(xx))

tg <- merge(tg,annotation, by.x='ID',by.y='names.xx.',all.x=TRUE)

# how many genes look diff express?
nrow(subset(tg,adj.P.Val < 0.05 & abs(logFC)>1 )) # 14

# check if these probes correspond to top Y probes
y <- as.data.frame(rownames(topprobes))
colnames(y) <- 'probe_id'

nrow(tg[tg$ID %in% y$probe_id, ]) # 24 ; 23?
# remove genes from Y chrm (those used to predict sex)
tg <- tg[! tg$ID %in% y$probe_id,  ]

nrow(subset(tg,adj.P.Val < 0.05 & abs(logFC)>1 )) # 4
nrow(subset(tg,adj.P.Val < 0.05)) # 105


# Volcano plot
with(tg, plot(logFC, -log10(adj.P.Val),pch=20,cex=0.4,main="Volcano Plot" ) )
with(subset(tg,adj.P.Val < 0.05 & abs(logFC)>1 ), points(logFC,-log10(adj.P.Val),pch=20,col='red' ) )

## save results :
save(test.norm.filt,s.results,meta,pheno, file='/users/cvalenci/GSE56580/Tcell.data.Rdata')

df <- df[unique(df$chi_id), ]

