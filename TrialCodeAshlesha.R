library(oligo)
library(affy)
library(ggplot2)
library(gcrma)
library(Biostrings)
celpath="/Users/ashleshagogate/Documents/StoriciLab/Analysis"
list = list.files(celpath,full.names=TRUE)
data = oligo::read.celfiles(list)
oligo::pm(data, "probesetname")
probeNrs = rep(rownames(pm(data, "probesetname")),7)
ints=c(oligo::pm(data, "probesetname")[,1],oligo::pm(data, "probesetname")[,2],oligo::pm(data, "probesetname")
       [,3],oligo::pm(data, "probesetname")[,4],oligo::pm(data, "probesetname")[,5],oligo::pm(data, "probesetname")[,6],oligo::pm(data, "probesetname")[,7])
pset = data.frame(probeNrs=probeNrs,ints=ints)
pset$PNs = factor(pset$probeNrs,levels=pset$probeNrs)
scatter = ggplot(pset,aes(PNs,ints))
scatter + geom_point() + labs(x="Probe Number",y="Intensity of PM probe")

# annotation

ph = data@phenoData
ph@data
colnames(ph)<-c("sample")
ph@data[ ,1] = c("HEK293_1","HEK293_2","H9_1","H9_2","H9_3","HEK293T_1", "HEK293T_2")

feat = data@featureData
feat@data
exp = data@experimentData
exp


#Plots to assess the quality of the data
color=c('green','green','red','red','red','blue','blue')
hist(data[,1:7],lwd=2,which='all',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram of raw data',target="core")

name = "boxplot1.jpg"
jpeg(name)
boxplot(data,which='all',col='red',names=ph@data$sample, target="core") 
dev.off()

df<-oligo::exprs(data)
color=c('green','green','red','yellow','black','blue','blue')
data.PC = prcomp(t(df),scale.=TRUE)
plot(data.PC$x[,1],data.PC$x[,2],col=color)

#Normalisation using RMA method
data.rma = oligo::rma(data, target="core", normalize=FALSE)
data.matrix = oligo::exprs(data.rma)

#Plots to assess the quality of the data after normalisation

name = "boxplotnorm.jpg"
jpeg(name)
boxplot(data.matrix,col='red',names=ph@data$sample)
dev.off()

color=c('green','green','red','red','red','blue','blue')
data.PC = prcomp(t(data.matrix),scale.=TRUE)
plot(data.PC$x[,1],data.PC$x[,2],col=color)

pmSeq<- oligo::pmSequence(data)
pmsLog2 <- log2(oligo::pm(data))
coefs <- oligo::getAffinitySplineCoefficients(pmsLog2, pmSeq)

counts <- Biostrings::alphabetFrequency(pmSeq, baseOnly=TRUE)
GCcontent <- ordered(counts[, "G"]+counts[, "C"])

