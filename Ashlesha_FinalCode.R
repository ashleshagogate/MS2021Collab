library(oligo)
library(affy)
library(ggplot2)
library(gcrma)
library(Biostrings)
library(limma)
library(Matrix)
library(lattice)
library(limma)
library(fdrtool)
library(rpart)
library(Biobase)

celpath="/Users/ashleshagogate/Documents/GT/StoriciLab/Analysis"
list = list.files(celpath,full.names=TRUE)
data = oligo::read.celfiles(list)

# annotation
expr = exprs(data)
ph = data@phenoData
colnames(ph)<-c("samples")
ph@data[ ,1] = c("HEK293", "HEK293", "HCT","HCT","H9", "H9", "HEK293T", "HEK293T","HELA","HELA","HELA","HEPG2","HEPG2","HEPG2","HUVEC","HUVEC","K562","K562","K562")
ph@data$Replicates<-c("1","2","1","2","1","2","1","2","1","2","3","1","2","3","1","2","1","2","3")
ph@data


#Plots to assess the quality of the data
color=c('red','red','blue','blue','green','green','yellow','yellow','pink','pink','pink','grey','grey','grey','orange','orange','purple','purple','purple')
hist(data[,1:19],lwd=2,which='pm',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram of raw data',target="core")

boxplot(data,which='pm',col=color, names=ph@data$sample,target="core",main='Boxplot before Normalization')


df<-oligo::exprs(data)
color=c('red','red','blue','blue','green','green','yellow','yellow','pink','pink','pink','grey','grey','grey','orange','orange','purple','purple','purple')
data.PC = prcomp(t(df),scale.=FALSE)
percentVar <- round(100*data.PC$sdev^2/sum(data.PC$sdev^2),1)
plot(data.PC$x[,1],data.PC$x[,2],col=color)

#Normalisation using RMA method
data.rma = oligo::rma(data, target="core", normalize=TRUE)
data.matrix = oligo::exprs(data.rma)
write.csv(data.matrix,file="ExpressionWithAllCellLines.csv")

#Plots to assess the quality of the data after normalisation
#Dataframe to Density plot
library(tidyr)


norm_data<-as.data.frame(data.matrix)
colnames(norm_data)<-paste(ph@data[,1],ph@data[,2],sep="_")
norm_data$Probe<-row.names(norm_data)
norm_data_gather<-gather(norm_data, key="Celllines", value="Log2RMA", -Probe)
#color=c('red','red','green','green','blue','blue')
names(color)=c("HEK293_1", "HEK293_2", "HCT_1","HCT_2","H9_1", "H9_2", "HEK293T_1", "HEK293T_2","HELA_1","HELA_2","HELA_3","HEPG2_1","HEPG2_2","HEPG2_3","HUVEC_1","HUVEC_2","K562_1","K562_2","K562_3")
ggplot(norm_data_gather, aes(x=Log2RMA, col=Celllines, main="Histogram after Normalization"))+
  geom_density(alpha=0)+
  theme_classic()+
  scale_color_manual(values=color)+
  scale_y_continuous(expand=c(0,0))

boxplot(data.matrix,col=color,names=ph@data$sample, target='core')

#color=c('red','red','green','green','green','blue','blue')
data.PC = prcomp(t(data.matrix),scale.=FALSE)
percentVar <- round(100*data.PC$sdev^2/sum(data.PC$sdev^2),1)
plot(data.PC$x[,1],data.PC$x[,2],col=color)
plot(data.PC$x[,2],data.PC$x[,3],col=color)
plot(data.PC$x[,3],data.PC$x[,4],col=color)


#limma
groups = ph@data$sample
f<- factor(groups,levels=c("HEK293", "HCT","H9","HEK293T","HELA","HEPG2","HUVEC","K562"))
design = model.matrix(~ 0 + f)
colnames(design) = c("HEK293", "HCT","H9","HEK293T","HELA","HEPG2","HUVEC","K562")
data.fit = lmFit(data.matrix,design)
data.fit$coefficients[1:10,]
avgexp<-as.data.frame(data.fit$coefficients)
write.csv(avgexp,file="LimmaExpressionAllCellLines.csv")
#baseline-test
contrast.matrix = makeContrasts(H9-HCT,levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)
names(data.fit.eb)
data.fit.eb$coefficients[1:10,]
r<-rownames(data.fit.eb)
df1<-as.data.frame(data.fit.eb, row.names=r)
df<-topTable(data.fit.eb, n=Inf)
up<-subset(df,(df$logFC>1.5&df$P.Value<0.05))
down<-subset(df,(df$logFC<(-1.5)&df$P.Value<0.05))
write.csv(up,file="UpregulatedinHCTvsH9.csv")
write.csv(down,file="DownregulatedinHCTvsH9.csv")
df$status<-"none"
df$status[df$logFC>1.5&df$P.Value<0.05]<-"up"
df$status[df$logFC<(-1.5)&df$P.Value<0.05]<-"down"
mycolors=c("blue","black","red")
names(mycolors)=c("down","none","up")
plot<- ggplot(data=df, aes(x=logFC, y=-log10(P.Value), col=status)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
plot
write.csv(df,file="HCTvsH9.csv")

library(huex10sttranscriptcluster.db)
library(annotate)
results<-select(huex10sttranscriptcluster.db, probes, c("SYMBOL","ENTREZID", "GENENAME"))
write.csv(results, file="AnnotationForAllCellLines.csv")

