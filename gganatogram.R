#gganatogram
#devtools::install_github("jespermaag/gganatogram")
library(gganatogram)
library(dplyr)
library(viridis)
library(gridExtra)


#example
organPlot <- data.frame(organ = c("heart", "leukocyte", "nerve", "brain", "liver", "stomach", "colon"), 
                        type = c("circulation", "circulation",  "nervous system", "nervous system", "digestion", "digestion", "digestion"), 
                        colour = c("red", "red", "purple", "purple", "orange", "orange", "orange"), 
                        value = c(100, 5, 100, 8000, 2000, 5, 5), 
                        stringsAsFactors=F)

head(organPlot)

gganatogram(data=organPlot, fillOutline='#a6bddb', organism='human', sex='male', fill="colour")+ 
  theme_void()+  scale_fill_gradient(low = "white", high = "red")


hgMale_key$organ

gganatogram(data=hgMale_key, fillOutline='#a6bddb', organism='human', sex='male', fill="colour") +theme_void()

gganatogram(data=organPlot, fillOutline='#a6bddb', organism='human', sex='female', fill="value") + 
  theme_void() +
  scale_fill_gradient(low = "white", high = "red")

#read tumor data for brca,liver,stomach,colon,Lung
brca <- read_delim("tumorExp/brca-rsem-fpkm-tcga-t.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
brcaRownames<-brca$Hugo_Symbol
brca<-brca[,3:dim(brca)[2]]
rownames(brca)<-brcaRownames

coad<-read_delim("tumorExp/coad-rsem-fpkm-tcga-t.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
coadRownames<-coad$Hugo_Symbol
coad<-coad[,3:dim(coad)[2]]
rownames(coad)<-coadRownames

lihc<-read_delim("tumorExp/lihc-rsem-fpkm-tcga-t.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
lihcRownames<-lihc$Hugo_Symbol
lihc<-lihc[,3:dim(lihc)[2]]
rownames(lihc)<-lihcRownames

luad<-read_delim("tumorExp/luad-rsem-fpkm-tcga-t.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
luadRownames<-luad$Hugo_Symbol
luad<-luad[,3:dim(luad)[2]]
rownames(luad)<-luadRownames

stad<-read_delim("tumorExp/stad-rsem-fpkm-tcga-t.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
stadRownames<-stad$Hugo_Symbol
stad<-stad[,3:dim(stad)[2]]
rownames(stad)<-stadRownames

#get mean exp of top10 BRCA genes in all samps
topGenes<-c("PIK3CA","TP53","TTN","GATA3","CDH1","MAP3K1","MUC16","KMT2C","MUC4","PTEN")
expVals<-c()
temp<-t(brca[topGenes,])
colnames(temp)<-topGenes
meanbrca<-colMeans(temp)

temp<-t(coad[topGenes,])
colnames(temp)<-topGenes
meancoad<-colMeans(temp)

temp<-t(lihc[topGenes,])
colnames(temp)<-topGenes
meanlihc<-colMeans(temp)

temp<-t(luad[topGenes,])
colnames(temp)<-topGenes
meanluad<-colMeans(temp)

temp<-t(stad[topGenes,])
colnames(temp)<-topGenes
meanstad<-colMeans(temp)

#expr of PIK#ca
organs<-c("breast","colon","liver","lung","stomach")
type<-c( "other","digestion","digestion","respiratory","digestion")
colour<-c("#41ab5d","orange","orange","steelblue","orange")
i<-10
vals<-c(meanbrca[i],meancoad[i],meanlihc[i],meanluad[i],meanstad[i])
gganatogramData<-data.frame(organ=organs,type=type,colour=colour,value=as.numeric(vals),stringsAsFactors=F)

gganatogram(data=gganatogramData, fillOutline='#a6bddb', organism='human', sex='female', fill="value")+ 
  theme_void()+  scale_fill_gradient(low = "yellow", high = "red",name= paste(topGenes[i],"(fpkm)"))+theme(legend.text = element_text(size=15,face = "bold"),legend.position = c(0.75, 0.2))

plots<-list() 
for(i in 1:10){
  vals<-c(meanbrca[i],meancoad[i],meanlihc[i],meanluad[i],meanstad[i])
  gganatogramData<-data.frame(organ=organs,type=type,colour=colour,value=as.numeric(vals),stringsAsFactors=F)
  
  p<-gganatogram(data=gganatogramData, fillOutline='#a6bddb', organism='human', sex='female', fill="value")+ 
    theme_void()+  scale_fill_gradient(low = "yellow", high = "red",name= paste(topGenes[i],"(fpkm)"))+theme(legend.text = element_text(size=15,face = "bold"),legend.position = c(0.75, 0.2))
  
  plots[[i]]<-p
}

grid.arrange(p, p, nrow = 1)
n <- length(plots)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(plots, ncol=5))
