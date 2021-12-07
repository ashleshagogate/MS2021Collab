#Combined Expression analysis

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(fdrtool))
suppressPackageStartupMessages(library(rpart))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(Biobase))             
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(gcrma))


celpath = "C:/Users/deepa/Dropbox (GaTech)/scripts/GitHub/MS2021Collab/CELfiles"
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)

ph = data@phenoData
ph@data
colnames(ph)<-c("samples")
ph@data[ ,1] = c("HEk293_1", "HEK293_2", "H9_1", "H9_2", "H9_3", "HEK293T_1", "HEK293T_2")
ph@data
expr = exprs(data)

#feat = data@featureData
#feat@data
#exp = data@experimentData
#exp
#individualhistograms
#for (i in 1:7){
#name = paste("histogram",i,".jpg",sep="")
#jpeg(name)
#hist(data[,i],lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main=ph@data$sample[i],target="core")
#dev.off()
#}

#Plots to assess the quality of the data <---------------------- needs to have a legend to show which cell line is which
#all
color=c('green','green','red','red','red','blue','blue')
hist(data[,1:7],lwd=2,which='all',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram of raw data',target="core")
#PM
hist(data[,1:7],lwd=2,which='pm',col=color,ylab='Density',xlab='Log2 intensities', main='Histogram Before Normalization', target="core")

#boxplot before normalization #Colours consistent with density/histogram plots
#<-------------------------Add labels to axes

name = "boxplot-all.jpg"
jpeg(name, width = 720, height = 480)
boxplot(data,which='all',col=color,names=ph@data$sample,target="core")
dev.off()
name = "boxplot-pm.jpg"
jpeg(name, width = 720, height = 480)
boxplot(data,which='pm',col=color,names=ph@data$sample,target="core")
dev.off()

#PCA before normalization <------------ Add legend or labels to the dots
color=c('green','green','red','red','red','blue','blue')
color=c('green','green','red','yellow','black','blue','blue') #This might not be correct for your presentation
exp_data <-(Biobase::exprs(data))
data.PC = prcomp(t(exp_data),scale.=FALSE)
summary(data.PC)

#autoplot(data.PC, data = ph@data, colour = 'samples',label = TRUE, label.size = 3)
plot(data.PC$x[,1],data.PC$x[,2],col=color, xlab="PCA1", ylab="PCA2", label = c("HEk293_1", "HEK293_2", "H9_1", "H9_2", "H9_3", "HEK293T_1", "HEK293T_2"))
plot(data.PC$x[,2],data.PC$x[,3],col=color, xlab="PCA2", ylab="PCA3")
plot(data.PC$x[,3],data.PC$x[,4],col=color, xlab="PCA3", ylab="PCA4")

#RMA # we want to normalize here!!
data.rma = oligo::rma(data, target="core", normalize=TRUE)
data.matrix = exprs(data.rma)

#Boxplot after normalization
name="boxplotnormalized.jpg"
boxplot(data.matrix,col=color,names=ph@data$sample, target='core')
jpeg(name, width = 720, height = 480)
dev.off()

#PCA after normalization
color=c('green','green','red','red','red','blue','blue')
data.PC = prcomp(t(data.matrix),scale.=FALSE)
summary(data.PC)
plot (data.PC$x[,1],data.PC$x[,2],col=color, xlab="PCA1", ylab="PCA2")
plot (data.PC$x[,2],data.PC$x[,3],col=color, xlab="PCA1", ylab="PCA2")
#plot (data.PC$x[,3],data.PC$x[,4],col=color, xlab="PCA1", ylab="PCA2")

#Histogram after normalization(something is wrong with this) 
#previously you used Exonfeatureset with this code, you might have to modify the chart or view it with a dfferent code. 
color=c('green','green','red','red','red','blue','blue')
hist(data.rma,lwd=2,which='pm',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram After Nomalization', target="core")

#Starting with GC
pmSeq<- oligo::pmSequence(data)
pmsLog2 <- log2(oligo::pm(data))
coefs <- oligo::getAffinitySplineCoefficients(pmsLog2, pmSeq)
counts <- Biostrings::alphabetFrequency(pmSeq, baseOnly=TRUE)
GCcontent <- (counts[, "G"]+counts[, "C"])/rowSums(counts)
GCcontent<-round(GCcontent, digits = 1)
boxplot(pmsLog2[,2]~GCcontent)
boxplot(data.matrix[,2]~GCcontent)

#limma
groups = ph@data$sample
f<- factor(groups,levels=c("HEK293","HEK293T","H9"))
design = model.matrix(~ 0 + f)
colnames(design) = c("HEK293","HEK293T","H9")
data.fit = lmFit(data.matrix,design)
data.fit$coefficients[1:10,]
contrast.matrix = makeContrasts(HEK293T-H9,levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)
names(data.fit.eb)
data.fit.eb$coefficients[1:10,]
r<-rownames(data.fit.eb)
df1<-as.data.frame(data.fit.eb, row.names=r)
df<-topTable(data.fit.eb, n=Inf)
up<-subset(df,(df$logFC>1&df$P.Value<0.05))
down<-subset(df,(df$logFC<(-1)&df$P.Value<0.05))
df$status<-"none"
df$status[df$logFC>1&df$P.Value<0.05]<-"up"
df$status[df$logFC<(-1)&df$P.Value<0.05]<-"down"
mycolors=c("blue","black","red")
names(mycolors)=c("down","none","up")
plot<- ggplot(data=df, aes(x=logFC, y=-log10(P.Value), col=status)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
plot

#Getting probe annotations from a database
library(huex10sttranscriptcluster.db)
library(annotate)
results<-select(huex10sttranscriptcluster.db, probes, c("SYMBOL","ENTREZID", "GENENAME"))
write.csv(results, file="AnnotationForAllCellLines.csv")









