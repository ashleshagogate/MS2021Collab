library(Matrix)
library(lattice)
library(limma)
library(fdrtool)
library(rpart)
library(ggplot2)
library(limma)
library(Biobase)
library(Biostrings)
library(affy)
library(oligo)
celpath = "/Users/kirti/Downloads/CELfiles"
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)
expr = exprs(data)
ph = data@phenoData
colnames(ph)<-c("samples")
ph@data[ ,1] = c("HEk293_1", "HEK293_2", "H9_1", "H9_2", "H9_3", "HEK293T_1", "HEK293T_2")
ph@data
expr[1:5,]
pm(data)[1:5,]
expr[51,]
pm(data)[5,]

#individualhistograms
for (i in 1:7){
  name = paste("histogram",i,".jpg",sep="")
  jpeg(name)
  hist(data[,i],lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main=ph@data$sample[i],target="core")
  dev.off()
}

#combinedhistogram
color=c('green','green','red','red','red','blue','blue')
hist(data[,1:7],lwd=2,which='pm',col=color,ylab='Density',xlab='Log2 intensities', main='Histogram Before Normalization', target="core")

#boxplot before normalization
name = "boxplot.jpg"
jpeg(name, width = 720, height = 480)
boxplot(data,which='pm',col='red',names=ph@data$sample,target="core")
dev.off()

#PCA before normalization
color=c('green','green','red','red','red','blue','blue')
exp_data <-(Biobase::exprs(data))
data.PC = prcomp(t(exp_data),scale.=TRUE)
plot (data.PC$x[,1],data.PC$x[,2],col=color, xlab="PCA1", ylab="PCA2")

#RMA
data.rma = rma(data,target='core',normalize=FALSE)
data.matrix = exprs(data.rma)

#Boxplot after normalization
name="boxplotnormalized.jpg"
boxplot(data.matrix,col='red',names=ph@data$sample, target='core')
jpeg(name, width = 720, height = 480)
dev.off()

#PCA after normalization
color=c('green','green','red','red','red','blue','blue')
data.PC = prcomp(t(data.matrix),scale.=TRUE)
plot (data.PC$x[,1],data.PC$x[,2],col=color, xlab="PCA1", ylab="PCA2")

#Histogram after normalization(something is wrong with this)
color=c('green','green','red','red','red','blue','blue')
hist(data.matrix[,1:7],lwd=2,which='pm',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram After Nomalization', target="core")


#Starting with GC
pmSeq <- pmSequence(data)
pmSeq
pmsLog2 <- log2(pm(data))

coefs <- getAffinitySplineCoefficients(pmsLog2, pmSeq)

counts <- Biostrings::alphabetFrequency(pmSeq, baseOnly=TRUE)
GCcontent <- ordered(counts[, "G"]+counts[, "C"])
