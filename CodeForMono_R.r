for (i in 1:NROW(t))
{
  subsetMono = rbind(subsetMono,sample.rows(subset(mono,X.Sample_characteristics_ch1 == t$X.Sample_characteristics_ch1[i] & X.Sample_characteristics_ch1.1 == t$X.Sample_characteristics_ch1.1[i] ),1))
}

files <- getGEOSuppFiles('GSE56045')
nonNormalizedExpression = read.delim("GSE56045/GSE56045_non_normalized.txt")
test<-read.ilmn("nonNormalizedExpressioncopy.txt",probeid= "ID_REF", expr='intensity')

test.norm <- neqc(test )
ids  <- as.character(rownames(test.norm))
qual <- unlist(mget(ids, illuminaHumanv4PROBEQUALITY, ifnotfound = NA ))
table(qual)
AveSignal = rowMeans(test.norm$E)

pdf('Plot1.pdf')
boxplot(AveSignal~qual) #plot average signals
dev.off()

rem <- qual=="No match "|qual=="Bad"
test.norm.filt <- test.norm[!rem ,]
dim(test.norm )
dim(test.norm.filt )

nrow(test.norm) - nrow(test.norm.filt)


data(y.probes)
illyprobe <- data.frame(y.probes["illumina_humanht_12"])
head(y.probes)

temp <- as.data.frame(test.norm.filt$E)
y.eset <- massi_y(temp, illyprobe)
class(temp)
class(y.eset)

length(y.eset$id)

pdf('test.Plot2.pdf')
massi_y_plot(y.eset)
dev.off()


gset <- getGEO("GSE56045", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pheno <- Biobase::pData(gset)
edata <- Biobase::exprs(gset)
topprobes <- massi_select(gset, illyprobe, threshold=4)
results <- massi_cluster(topprobes)

s.results <- data.frame(results[[2]])

head(s.results)

table(s.results$sex)

meta <- pData(gset)
pheno <- meta[,c(1,2,12,14:17)]

temp <- s.results[,c(1,5)]
temp <- merge(pheno,temp, by.x='geo_accession',by.y='ID')

pheno <- temp
head(pheno)
colnames(pheno)[3] <- 'Age_years'
pheno$Age_years <- gsub('.*:', '', pheno$Age_years)
pheno$Age_years <- as.numeric(as.character(pheno$Age_years))
pheno$chip = gsub('_.*','',pheno$title)

pheno_obs = pheno
pheno = merge(pheno_obs,subsetMono, by.y = 'Sample_geo_accession', by.x = 'geo_accession' )
pheno = pheno[,c(1,2,7)]


dmat1 <- model.matrix(~pheno_trial$Age_years)
colnames(dmat1) <- c('Intercept', 'Age')
dmat1

fit.ls <- lmFit(test.norm.filt,dmat1,method="ls")

eb.ls <- eBayes(fit.ls)
names(eb.ls)
class(eb.ls)

pdf('plot4.pvalueMod1.pdf')
hist(eb.ls$p.value[,2], xlab = "p-value", main = "")
dev.off()

table <- topTable(eb.ls, coef = 2, adjust.method = "BH", genelist = rownames(eb.ls), number=Inf)
head(table)

pdf('Vplot.Mod1.pdf')
with(table, plot(logFC, -log10(adj.P.Val),pch=20,cex=0.4,main="Volcano Plot" ) )
with(subset(table,adj.P.Val < 0.05 ), points(logFC,-log10(adj.P.Val),pch=20,col='red' ) )
dev.off()



dmat2 <- model.matrix(~pheno_trial$Age_years + pheno_trial$sex)
colnames(dmat2) <- c('Intercept', 'Age', 'Sex')
dmat2

fit.ls2 <- lmFit(test.norm.filt,dmat2,method="ls")

eb.ls2 <- eBayes(fit.ls2)
names(eb.ls2)
class(eb.ls2)

pdf('plot6.pvalueMod1.pdf')
hist(eb.ls2$p.value[,2], xlab = "p-value", main = "")
dev.off()


table2 <- topTable(eb.ls2, coef = 2, adjust.method = "BH", genelist = rownames(eb.ls2), number=Inf)
head(table2)

pdf('Vplot.Mod2.pdf')
with(table2, plot(logFC, -log10(adj.P.Val),pch=20,cex=0.4,main="Volcano Plot" ) )
with(subset(table2,adj.P.Val < 0.05), points(logFC,-log10(adj.P.Val),pch=20,col='red' ) )
dev.off()



pheno_trial$Age_cat <- ifelse(pheno_trial$Age_years >= 58,1,0)


dmat3 <- model.matrix(~pheno_trial$Age_cat + pheno_trial$sex)
colnames(dmat3) <- c('Intercept', 'Age', 'Sex')
dmat3

fit.ls3 <- lmFit(test.norm.filt,dmat3,method="ls")

eb.ls3 <- eBayes(fit.ls3)
names(eb.ls3)
class(eb.ls3)

pdf('plot7.pvalueMod1.pdf')
hist(eb.ls3$p.value[,2], xlab = "p-value", main = "")
dev.off()


table3 <- topTable(eb.ls3, coef = 2, adjust.method = "BH", genelist = rownames(eb.ls3), number=Inf)
head(table3)

pdf('Vplot.Mod3.pdf')
with(table3, plot(logFC, -log10(adj.P.Val),pch=20,cex=0.4,main="Volcano Plot" ) )
with(subset(table3,adj.P.Val < 0.05), points(logFC,-log10(adj.P.Val),pch=20,col='red' ) )
dev.off()