#a collection of functions/code to use TCGA metadata

library(plyr)
library("readr")
library(TCGAbiolinks)


#function to expand columns of TCGA metadata from list to df
expand<-function(df,colName){
  res<-data.frame()
  #for each row
  for(i in 1: dim(df)[1]){
    thisRow<-df[i, ! (colnames(df) %in% c(colName))]
    tempdf<-as.data.frame(df[i, c(colName)])
    #if list is empty skip that row
    if(dim(tempdf)[1]<1){
      next
    }
    #change colnames so they are unique
    colnames(tempdf)<-paste(paste(colName,".",sep = ""),colnames(tempdf),sep = "")
    #print(paste(i,colnames(tempdf)))
    newRow<-cbind(thisRow,tempdf,row.names = NULL)
    res<-bind_rows(res,newRow)
    
  }
  #print(res)
  return(res)
}

#function to download and combine TCGA clinical and Biospecimen metadata given a project name
#example usage: getTCGAMetadata("TCGA-BRCA")
getTCGAMetadata<-function(projName){
  print(paste("Downloading",projName))
  clinicalBRCA <- GDCquery_clinic(project = projName, type = "clinical")
  biospecimenBRCA <- GDCquery_clinic(project = projName, type = "Biospecimen")
  
  #rename all cols from clinical table with suffix clinical
  colnames(clinicalBRCA)<- paste0("clinical.",colnames(clinicalBRCA))
  
  #expand biospecimen data in the order portions, portions.analytes, portions.analytes.aliquots
  toUnpack<-c("portions", "portions.analytes", "portions.analytes.aliquots")
  for(s in toUnpack){
    biospecimenBRCA<-expand(biospecimenBRCA,s)
  }
  #add patient barcode to biospecimen data
  biospecimenBRCA<- biospecimenBRCA %>% mutate(clinical.bcr_patient_barcode=substr(submitter_id,1,nchar(as.character(submitter_id))-4))
  #join clinical and biospecimen
  finalJoined<-join(clinicalBRCA,biospecimenBRCA,by="clinical.bcr_patient_barcode")
  return(finalJoined)
}

#function to examine dimentions of TCGA metadata
printDatadim<-function(projName){
  print(paste("Downloading",projName))
  clinicalBRCA <- GDCquery_clinic(project = projName, type = "clinical")
  biospecimenBRCA <- GDCquery_clinic(project = projName, type = "Biospecimen")
  print(paste("clinical dim:",dim(clinicalBRCA)))
  print(paste("biospc dim:",dim(biospecimenBRCA)))
}

#a list of useful columns to retain from the TCGA metadata (to reduce the dimensions)
colsToKeep<-c("clinical.submitter_id",
              "clinical.classification_of_tumor",
              "clinical.primary_diagnosis",
              "clinical.tumor_stage",
              "clinical.age_at_diagnosis",
              "clinical.vital_status",
              "clinical.days_to_death",
              "clinical.tissue_or_organ_of_origin",
              "clinical.days_to_birth",
              "clinical.site_of_resection_or_biopsy",
              "clinical.days_to_last_follow_up",
              "clinical.cigarettes_per_day",
              "clinical.weight",
              "clinical.alcohol_history",
              "clinical.bmi",
              "clinical.years_smoked",
              "clinical.height",
              "clinical.gender",
              "clinical.year_of_birth",
              "clinical.race",
              "clinical.ethnicity",
              "clinical.year_of_death",
              "clinical.bcr_patient_barcode",
              "clinical.disease",
              "submitter_id",
              "sample_type",
              "tissue_type",
              "portions.submitter_id",
              "portions.analytes.analyte_type",
              "portions.analytes.submitter_id",
              "portions.analytes.analyte_type_id",
              "portions.analytes.aliquots.analyte_type",
              "portions.analytes.aliquots.submitter_id")

#################Examples#######################

#download BRCA metadata
brcaMetadata<-getTCGAMetadata("TCGA-BRCA")
#keep only selected columns
brcaMetadata<-brcaMetadata[,colsToKeep]

#download metadata of following projects into a single dataframe 
tcgaProjList<-c("TCGA-BLCA","TCGA-HNSC","TCGA-ESCA","TCGA-PRAD")
#mdList will have all metadata for tcgaProjList
mdListDF<-data.frame()
for(s in tcgaProjList){
  #mdList<-c(mdList,getjoinedBiospcClinc(s))
  if(dim(mdListDF)[1]<1){
    mdListDF<-getjoinedBiospcClinc(s)
  }else{
    print("joining")
    temp<-getjoinedBiospcClinc(s)
    mdListDF<-bind_rows(mdListDF,temp)  
  }
}


