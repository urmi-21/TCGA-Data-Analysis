#process mutation data
library(dplyr)
library(DT)
library(plyr)
library(data.table)
library(maftools)
library("readr")
library(TCGAbiolinks)
library(ggplot2)
packageVersion("TCGAbiolinks")

###########################################define functions##############################################################
#split maf file into categories by a given variable in the clinical metadata file
splitMafby<-function(clinicalData,by,mafData){
  #get tumor samps fo different values of by
  myList<-list()
  uniqVals<-clinicalData%>%select(by)%>%unique
  for(i in 1:nrow(uniqVals)){
    #i<-1
    s<-as.character(uniqVals[i,1])
    print(s)
    tList<-clinicalData %>% filter(clinicalData[,by]==s) %>% select(portions.analytes.aliquots.submitter_id)%>%unique
    thisData<-mafData%>%filter(Tumor_Sample_Barcode %in% tList$portions.analytes.aliquots.submitter_id)
    myList[[s]]<-thisData
  }
  
  return(myList)
}

plotSummaryTofile<-function(mafList,fname){
  l1<-mafList
  #for each item in list do calculations
  lnames<-names(l1)
  plotList<-list()
  k<-0
  pdf(fname)
  for(i in lnames){
    print((i))
    print(dim(l1[[i]]))
    if(dim(l1[[i]])[1]<1){
      next
    }
    #plot summary and save to pdf
    maf<-read.maf(l1[[i]],isTCGA = T)
    plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = T, titvRaw = FALSE,showBarcodes=F, top = 20)
    mtext(paste("Mutation summary by",i), outer=T,  cex=NA, line=-1.5,side = 3)
  }
  dev.off()
}


##############################################################################################################################
#read clinical metadata
TCGAbrcaMetadata_reduced <- read_csv("TCGAbrcaMetadata_reduced.csv")
colnames(TCGAbrcaMetadata_reduced)[which(colnames(TCGAbrcaMetadata_reduced)=="portions.analytes.aliquots.submitter_id")]="Tumor_Sample_Barcode"
##download gene mutation metadata
brcaMAF <- GDCquery_Maf("BRCA", pipelines = "varscan2")

brcaMAF <- GDCquery_Maf("LIHC", pipelines = "varscan2")

write_tsv(brcaMAF,"LIHC_mutations.tsv")
#extract mutated genes list
brcaMAFFreq<-brcaMAF %>% select(Hugo_Symbol) %>% group_by(Hugo_Symbol) %>% count %>% arrange(desc(freq))
write_tsv(brcaMAFFreq,"LIHC_mutations_freq.tsv")

par(mfrow=c(2,2))
plotmafSummary(maf = read.maf(brcaMAF), rmOutlier = TRUE, addStat = 'median', dashboard = T, titvRaw = FALSE,showBarcodes=F, top = 30)
par(mfrow=c(1,1)) 
brcaMAF_MD<-join(brcaMAF,TCGAbrcaMetadata_reduced)
ggplot(data=brcaMAF_MD,aes(x=clinical.race,fill=VARIANT_CLASS))+geom_bar()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

#what percent of mutations are in BRCA1 BRCA 2
brcaMAF_MD_geneGroupBRCA<-brcaMAF_MD%>%select(Hugo_Symbol,Variant_Classification)%>%mutate(geneName=ifelse(Hugo_Symbol=="BRCA1","BRCA1",ifelse(Hugo_Symbol=="BRCA2","BRCA2","Other"))) %>%select(geneName)%>% group_by(geneName) %>% count %>% mutate(logFreq=log(freq))

p1<-ggplot(data=brcaMAF_MD_geneGroupBRCA, aes(x=geneName,y=logFreq,fill=geneName))+geom_bar(stat = "identity")+theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 15,face = "bold"),axis.text.y = element_text(size = 10,face = "bold"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ scale_fill_brewer(palette="Paired")



p2<-ggplot(data=brcaMAF_MD%>%filter(Hugo_Symbol %in% c("BRCA1","BRCA2")),aes(x=Hugo_Symbol,fill=Variant_Classification))+geom_bar()+
  theme(axis.text.x = element_text(size = 15,face = "bold"),axis.text.y = element_text(size = 10,face = "bold"),legend.text=element_text(size = 10,face = "bold"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ scale_fill_brewer(palette="Dark2")

grid.arrange(p1, p2, nrow = 1)
#count mutations in BRCA1/2
brcaMAF_brca <- brcaMAF %>% filter(Hugo_Symbol %in% c("BRCA1","BRCA2"))
brcaMAF_brca<-join(brcaMAF_brca,TCGAbrcaMetadata_reduced,type="left")

ggplot(data=brcaMAF_brca,aes(x=Hugo_Symbol,fill=VARIANT_CLASS))+geom_bar()
ggplot(data=brcaMAF_brca,aes(x=Hugo_Symbol,fill=clinical.race))+geom_bar()





tempList<-list()
tempList[["BRCADataset"]]<-brcaMAF
plotSummaryTofile(tempList,"mutationSummary20.pdf")

l1<-splitMafby(TCGAbrcaMetadata_reduced,"clinical.race",brcaMAF)
plotSummaryTofile(l1,"mutationByRace.pdf")

l2<-splitMafby(TCGAbrcaMetadata_reduced,"clinical.primary_diagnosis",brcaMAF)
plotSummaryTofile(l2,"mutationByprimarydiagnosis.pdf")

l3<-splitMafby(TCGAbrcaMetadata_reduced,"clinical.gender",brcaMAF)
plotSummaryTofile(l3,"mutationBygender.pdf")

primaryDiagnosisFreq<-TCGAbrcaMetadata_reduced%>%select(clinical.primary_diagnosis)%>%group_by(clinical.primary_diagnosis)%>%count()%>%arrange(desc(freq))%>%mutate(logfreq=log10(freq))
primaryDiagnosisFreq$clinical.primary_diagnosis<-factor(primaryDiagnosisFreq$clinical.primary_diagnosis, levels = primaryDiagnosisFreq$clinical.primary_diagnosis)
ggplot(data=primaryDiagnosisFreq,aes(x=clinical.primary_diagnosis,y=logfreq,fill=freq))+geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

#group by diagnosis and gender
primaryDiagnosisFreq<-TCGAbrcaMetadata_reduced%>%select(clinical.gender,clinical.primary_diagnosis)%>%group_by(clinical.gender,clinical.primary_diagnosis)%>%count()%>%arrange(desc(freq))%>%mutate(logfreq=log10(freq))
ggplot(data=primaryDiagnosisFreq,aes(x=clinical.primary_diagnosis,y=logfreq,fill=clinical.gender))+geom_bar(stat = "identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#plot dnarna
ggplot(data=TCGAbrcaMetadata_reduced,aes(x=sample_type,fill=portions.analytes.analyte_type))+geom_bar()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

maf<-read.maf(l1[["white"]],isTCGA = T)
drugInteractions(maf = maf, fontSize = 0.75)
drugInteractions(genes = "DNMT3A", drugs = TRUE)

whiteMaf<-l1[["white"]]
#maftest<-whiteMaf %>% select(Hugo_Symbol,Variant_Classification) %>% group_by(Variant_Classification) %>% count %>% group_by(Hugo_Symbol) %>% count %>% filter(Hugo_Symbol =="USO1")

maftest<-brcaMAF %>% select(Hugo_Symbol,Variant_Classification) %>% group_by(Hugo_Symbol) %>% count %>% arrange(desc(freq)) %>% top_n(n=20)

#find top mutated genes
topGenes<-brcaMAF %>% filter(Variant_Classification != "Silent") %>% select(Hugo_Symbol) %>% group_by(Hugo_Symbol) %>% count %>% arrange(desc(freq)) %>% top_n(n=10)

maftest<-brcaMAF %>% filter(Hugo_Symbol %in% topGenes$Hugo_Symbol & Variant_Classification != "Silent") %>% select(Hugo_Symbol,Variant_Classification) 

#for fill colors
library("RColorBrewer", lib.loc="~/R/win-library/3.5")
colourCount = length(unique(maftest$Variant_Classification))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# assign colors manually to be consistent with missing data
variants<-c(
  "Missense_Mutation",
  "Silent",
  "3'UTR",
  "Nonsense_Mutation",
  "5'Flank",
  "Intron",
  "Splice_Region",
  "RNA",
  "5'UTR",
  "Splice_Site",
  "In_Frame_Del",
  "Frame_Shift_Ins",
  "Frame_Shift_Del",
  "In_Frame_Ins",
  "3'Flank",
  "Nonstop_Mutation",
  "Translation_Start_Site",
  "IGR"
)
colourCount= length(unique(variants))
palette <- getPalette(colourCount)
names(palette)<-variants

ggplot(data=maftest,aes(x=Hugo_Symbol,fill=Variant_Classification))+geom_bar(stat = "count")+ scale_x_discrete(limits = rev(topGenes$Hugo_Symbol))+coord_flip()+
  theme(legend.position = "top",legend.text = element_text(size = 11),legend.title = element_blank())+
  theme(axis.text.x = element_text(angle=45,size = 11,face = "bold"),axis.text.y = element_text(size = 11,face = "bold"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ 
    ylab("")+xlab("")+scale_fill_manual(values = palette)+ ggtitle("Petal and sepal length of iris")


plotGeneVarFreq<-function(mafList){
  #for each item in list 
  lnames<-names(mafList)
  plotList<-list()
  k<-1
  plist<-list()
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  # assign colors manually to be consistent with missing data
  variants<-c(
    "Missense_Mutation",
    "Silent",
    "3'UTR",
    "Nonsense_Mutation",
    "5'Flank",
    "Intron",
    "Splice_Region",
    "RNA",
    "5'UTR",
    "Splice_Site",
    "In_Frame_Del",
    "Frame_Shift_Ins",
    "Frame_Shift_Del",
    "In_Frame_Ins",
    "3'Flank",
    "Nonstop_Mutation",
    "Translation_Start_Site",
    "IGR"
  )
  colourCount= length(unique(variants))
  palette <- getPalette(colourCount)
  names(palette)<-variants
  
  for(i in lnames){
    print(i)
    print(dim(mafList[[i]]))
    if(dim(mafList[[i]])[1]<1){
      next
    }
    
    #find top mutated genes
    thisMaf<-mafList[[i]]
    topGenes<-thisMaf %>% filter(Variant_Classification != "Silent") %>% select(Hugo_Symbol) %>% group_by(Hugo_Symbol) %>% count %>% arrange(desc(freq)) %>% top_n(n=10)
    print(topGenes)
    maftest<-thisMaf %>% filter(Hugo_Symbol %in% topGenes$Hugo_Symbol & Variant_Classification != "Silent") %>% select(Hugo_Symbol,Variant_Classification) 
    
    
    p<-ggplot(data=maftest,aes(x=Hugo_Symbol,fill=Variant_Classification))+geom_bar(stat = "count")+ scale_x_discrete(limits = rev(topGenes$Hugo_Symbol))+coord_flip()+
      theme(legend.position = "right")+
      theme(axis.text.x = element_text(size = 11,face = "bold"),axis.text.y = element_text(size = 11,face = "bold"),panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ 
      ylab("")+xlab("")+scale_fill_manual(values = palette)+ ggtitle(i)
    
    plist[[k]]<-p
    k=k+1
    
  }
  #arrange plots on a grid
  do.call("grid_arrange_shared_legend", c(plist))
  
  
}

#plot on a grid with single legend
#function reference https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

names(l1)
l1<-l1[c("white","black or african american","asian" )]
l1[[1]]
plotGeneVarFreq(l1)

filterList<- brcaMAF_MD %>% select(clinical.primary_diagnosis) %>% count %>% arrange(desc(freq)) %>% top_n(n=2)
l2<-splitMafby(TCGAbrcaMetadata_reduced%>%filter(clinical.primary_diagnosis %in% filterList$clinical.primary_diagnosis),"clinical.primary_diagnosis",brcaMAF)
names(l2)
plotGeneVarFreq(l2)

l3<-splitMafby(TCGAbrcaMetadata_reduced,"clinical.gender",brcaMAF)
plotGeneVarFreq(l3)
