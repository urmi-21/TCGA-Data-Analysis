library(readr)
library(dplyr)
library(ggplot2)
library(reshape)
#Read expression data
brca_nontumor <- read_delim("brca-rsem-fpkm-tcga.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
brca_tumor <- read_delim("brca-rsem-fpkm-tcga-t.txt", "\t", escape_double = FALSE, trim_ws = TRUE)


#formatting
rn<-brca_nontumor$Hugo_Symbol
brca_nontumor<-brca_nontumor[,3:dim(brca_nontumor)[2]]
rownames(brca_nontumor)<-rn

rn<-brca_tumor$Hugo_Symbol
brca_tumor<-brca_tumor[,3:dim(brca_tumor)[2]]
rownames(brca_tumor)<-rn

#find correlation of tp53
tp53data<-brca_nontumor[which(rownames(brca_nontumor)=="TP53"),]

brca_nontumorCor<-cor(t(brca_nontumor))

x<-brca_nontumorCor[which(rownames(brca_nontumorCor)=="TP53"),]
sort(x, decreasing = T)

#plot box plots of top mutated genes
topGenes<-c("PIK3CA","TP53","TTN","GATA3","CDH1","MAP3K1","MUC16","KMT2C","MUC4","PTEN")
nt<-as.data.frame(t(brca_nontumor[topGenes,]))
colnames(nt)<-topGenes


ggplot(melt(nt), aes(x=factor(variable),y=value,fill=factor(variable)))+geom_boxplot()+scale_y_log10()+theme(legend.position = "top")+
  theme(axis.text.x = element_text(angle=45,size = 15,face = "bold"),axis.text.y = element_text(size = 10,face = "bold"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ scale_fill_discrete(name = "Gene")+ylab("log(expression)")+xlab("")

#tumordata
tdata<-as.data.frame(t(brca_tumor[topGenes,]))
colnames(tdata)<-topGenes

ggplot(melt(tdata), aes(x=factor(variable),y=value,fill=factor(variable)))+geom_boxplot()+scale_y_log10()+theme(legend.position = "top",legend.text = element_text(size = 10,face = "bold"))+
  theme(axis.text.x = element_text(angle=45,size = 15,face = "bold"),axis.text.y = element_text(size = 10,face = "bold"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ scale_fill_discrete(name = "Gene")+ylab("log(expression)")+xlab("")

nt<-nt%>%mutate(source="Non-Tumor")
tdata<-tdata%>%mutate(source="Tumor")

topCombined<-bind_rows(nt,tdata)
ggplot(melt(topCombined), aes(x=factor(variable),y=value,fill=factor(source)))+geom_boxplot(outlier.colour=NA)+scale_y_log10()+theme(legend.position = "top",legend.text = element_text(size = 15,face = "bold"))+
  theme(axis.text.x = element_text(angle=45,size = 15,face = "bold"),axis.text.y = element_text(size = 10,face = "bold"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ scale_fill_discrete(name = "Gene")+ylab("log(expression)")+xlab("")+
  scale_fill_brewer(palette="Set1")


#ref for data https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5470989/
upGenes<-c("E2F1", "EZH2", "FOXM1", "MYBL2", "PLK1")
dwnGenes<-c("SCARA5", "MYOM1", "NKAPL", "PEG3", "USP2")


ntUP<-as.data.frame(t(brca_nontumor[upGenes,]))
colnames(ntUP)<-upGenes
ntUP<-ntUP%>%mutate(source="Non-Tumor")
tdataUP<-as.data.frame(t(brca_tumor[upGenes,]))
colnames(tdataUP)<-upGenes
tdataUP<-tdataUP%>%mutate(source="Tumor")

diffExpgenes<-bind_rows(ntUP,tdataUP)
dodge <- position_dodge(width = 0.5)

ggplot(melt(diffExpgenes), aes(x=factor(variable),y=value,fill=factor(source)))+geom_violin(position = dodge)+geom_boxplot(width=.1, outlier.colour=NA, position = dodge)+scale_y_log10()+theme(legend.position = "top")+
  theme(axis.text.x = element_text(angle=45,size = 15,face = "bold"),axis.text.y = element_text(size = 10,face = "bold"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ scale_fill_discrete(name = "Gene")+ylab("log(expression)")+xlab("")+scale_fill_brewer(palette="Set1")

#for down regulated in cancer
ntDWN<-as.data.frame(t(brca_nontumor[dwnGenes,]))
colnames(ntDWN)<-dwnGenes
ntDWN<-ntDWN%>%mutate(source="Non-Tumor")
tdataDWN<-as.data.frame(t(brca_tumor[dwnGenes,]))
colnames(tdataDWN)<-dwnGenes
tdataDWN<-tdataDWN%>%mutate(source="Tumor")
diffExpgenesDWN<-bind_rows(ntDWN,tdataDWN)


ggplot(melt(diffExpgenesDWN), aes(x=factor(variable),y=value,fill=factor(source)))+geom_violin(position = dodge)+geom_boxplot(width=.1, outlier.colour=NA, position = dodge)+scale_y_log10()+theme(legend.position = "top")+
  theme(axis.text.x = element_text(angle=45,size = 15,face = "bold"),axis.text.y = element_text(size = 10,face = "bold"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ scale_fill_discrete(name = "Gene")+ylab("log(expression)")+xlab("")+scale_fill_brewer(palette="Set1")


