#R script process mutation data
library(dplyr)
library(maftools)
library(TCGAbiolinks)



#function to download maf file given a project name
getmaf<-function(projName){
  thisMaf<-GDCquery_Maf(projName, pipelines = "varscan2")
  return(thisMaf)
}

#function to download maf files given a list project names. Returned files is combined into one bigger file.
getmafs<-function(projList){
  allMaf<-data.frame()
  for(p in projList){
    temp<-GDCquery_Maf(p, pipelines = "varscan2")
    allMaf<-rbind(allMaf,temp)
  }
  return(allMaf)
}

#function to extract transcript information from maf file
geneToTranscript<-function(mafFile){
  return(distinct(mafFile[,c("Transcript_ID","Gene","Hugo_Symbol","ENSP","SWISSPROT","TREMBL","BIOTYPE","RefSeq")]))
}

###################################Example usage############################
tcgamafProjList<-c("BLCA","BRCA","CESC","UCEC","UCS","READ","COAD","LIHC","HNSC","ESCA","PRAD","STAD","THCA","LUAD","LUSC","KIRC","KIRP","KICH")
ucsmaf<-getmaf("UCS")
ucecmaf<-getmaf("UCEC")
uvmmaf<-getmaf("UVM")
mafs<-getmafs(c("UCS","UVM","BLCA"))

#read maf file with read.maf function
maf<-read.maf(mafs,isTCGA = T)
plotmafSummary(maf = maf)


#extract transcripts from maf
ucsTids<-geneToTranscript(ucsmaf)



