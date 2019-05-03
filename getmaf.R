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

#download and combine maf files
getmaf<-function(projName){
  thisMaf<-GDCquery_Maf(projName, pipelines = "varscan2")
  return(thisMaf)
}

getmafs<-function(projList){
  allMaf<-data.frame()
  for(p in projList){
    temp<-GDCquery_Maf(p, pipelines = "varscan2")
    allMaf<-rbind(allMaf,temp)
  }
  return(allMaf)
}

geneToTranscript<-function(mafFile){
  return(distinct(mafFile[,c("Transcript_ID","Gene","Hugo_Symbol","ENSP","SWISSPROT","TREMBL","BIOTYPE","RefSeq")]))
}

############################################################################
tcgamafProjList<-c("BLCA","BRCA","CESC","UCEC","UCS","READ","COAD","LIHC","HNSC","ESCA","PRAD","STAD","THCA","LUAD","LUSC","KIRC","KIRP","KICH")
ucsmaf<-getmaf("UCS")
ucecmaf<-getmaf("UCEC")
uvmmaf<-getmaf("UVM")
mafs<-getmafs(c("UCS","UVM","BLCA"))

maf<-read.maf(mafs,isTCGA = T)
plotmafSummary(maf = maf)


colnames(brcaMAF)
#geneDataCols<-c( "Hugo_Symbol","Entrez_Gene_Id","NCBI_Build", "Chromosome","Start_Position","End_Position","Strand","Gene"  )

brcaTids<-geneToTranscript(brcaMAF)

#create a list of tids

