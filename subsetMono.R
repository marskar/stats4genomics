#Read in mono.txt
#Read in t.txt
#These are curated files from GSE56045_series_matrix.txt and GSE56580_series_matrix.txt

mono = read.delim("~/mono.txt",sep = "\t")
t = read.delim("~/t.txt",sep = "\t")
subsetMono = mono[FALSE,]

for (i in 1:NROW(t))
{
  subsetMono = rbind(subsetMono,sample.rows(subset(mono,X.Sample_characteristics_ch1 == t$X.Sample_characteristics_ch1[i] & X.Sample_characteristics_ch1.1 == t$X.Sample_characteristics_ch1.1[i] ),1))
}
