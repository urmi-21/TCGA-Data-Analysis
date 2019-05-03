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

#split maf file into categories by a given variable in the TCGA metadata file
#e.g., maf_race<-splitMafby(brcaMetadata,"clinical.race",brcaMAF)
#This will return a list of maf files each coressponding to a "clinical.race" category
#For information of how to download metadata see "tcgaMetadataFunctions.R"
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


#This function takes a list of maf files as input and saves the summary to a .pdf file using "plotmafSummary {maftools}"
#e.g., plotSummaryTofile(splitMafby(brcaMetadata_reduced,"clinical.race",brcaMAF),"mutationSummary20.pdf")
#this will produce variant summary of BRCA data separated by clinical.race category
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



