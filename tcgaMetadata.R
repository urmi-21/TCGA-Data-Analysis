
library(dplyr)
library(DT)
library(plyr)
library(data.table)
library(maftools)
library("readr")
library(TCGAbiolinks)
packageVersion("TCGAbiolinks")

#######################################################################################
colsToKeep<-c("clinical.submitter_id","clinical.classification_of_tumor","clinical.primary_diagnosis","clinical.tumor_stage","clinical.age_at_diagnosis","clinical.vital_status","clinical.days_to_death","clinical.tissue_or_organ_of_origin","clinical.days_to_birth","clinical.site_of_resection_or_biopsy","clinical.days_to_last_follow_up","clinical.cigarettes_per_day","clinical.weight","clinical.alcohol_history","clinical.bmi","clinical.years_smoked","clinical.height","clinical.gender","clinical.year_of_birth","clinical.race","clinical.ethnicity","clinical.year_of_death","clinical.bcr_patient_barcode","clinical.disease","submitter_id","sample_type","tissue_type","portions.submitter_id","portions.analytes.analyte_type","portions.analytes.submitter_id","portions.analytes.analyte_type_id","portions.analytes.aliquots.analyte_type","portions.analytes.aliquots.submitter_id")
#Function takes a df and expands it by unlisting elements at a column
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

getjoinedBiospcClinc<-function(projName){
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

printDatadim<-function(projName){
  print(paste("Downloading",projName))
  clinicalBRCA <- GDCquery_clinic(project = projName, type = "clinical")
  biospecimenBRCA <- GDCquery_clinic(project = projName, type = "Biospecimen")
  print(paste("clinical dim:",dim(clinicalBRCA)))
  print(paste("biospc dim:",dim(biospecimenBRCA)))
}

##########################End Functions##########################################

#download and merge BRCA metadata
BRCAMetadata<-getjoinedBiospcClinc("TCGA-BRCA")
clinical <- GDCquery_clinic(project = "TCGA-UCS", type = "clinical")

tcgaProjList<-c("TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-UCEC","TCGA-UCS","TCGA-READ","TCGA-COAD","TCGA-LIHC","TCGA-HNSC","TCGA-ESCA","TCGA-PRAD","TCGA-STAD","TCGA-THCA","TCGA-LUAD","TCGA-LUSC","TCGA-KIRC","TCGA-KIRP","TCGA-KICH")
#tcgaProjList<-c("TCGA-BLCA","TCGA-HNSC","TCGA-ESCA","TCGA-PRAD")

#mdList will have all data frames for rcgaProjList
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

#remove cols with all NA values
naCols<-colnames(mdListDF)[sapply(mdListDF, function(x)all(is.na(x)))]
mdListDFNONA<-mdListDF[,!(colnames(mdListDF) %in% naCols)]
#keep rows with RNA samples only
mdListDFRNA<-mdListDF%>%filter(portions.analytes.analyte_type_id == "R")

ulMD<-unlist(mdList)
mdJoined<-rbindlist(unlist(mdList))
n1<-colnames(t)
n2<-colnames(BRCAMetadata)
n3<-colnames(mdListDF)
setdiff(n2,n1)
#"updated_datetime" "submitter_id"     "created_datetime" "state" colnames are repeated


biospecimentest<- GDCquery_clinic(project = "TCGA-UCS", type = "Biospecimen")
clinicaltest<- GDCquery_clinic(project = "TCGA-UCS", type = "Clinical")

length(colnames(biospecimentest))
length(unique(colnames(biospecimentest)))

length(colnames(clinicaltest))
length(unique(colnames(clinicaltest)))

intersect(colnames(biospecimentest),colnames(clinicaltest))


#check data dim for all projects
start.time <- Sys.time()
for(s in tcgaProjList){
  printDatadim(s)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

biospecimenUCS<- GDCquery_clinic(project = "TCGA-UCS", type = "Biospecimen")
biospecimenCESC<- GDCquery_clinic(project = "TCGA-CESC", type = "Biospecimen")
clinicalUCS<- GDCquery_clinic(project = "TCGA-UCS", type = "Clinical")
clinicalCESC<- GDCquery_clinic(project = "TCGA-CESC", type = "Clinical")

cnClUCS<-colnames(clinicalUCS)
cnClCESC<-colnames(clinicalCESC)
cnBSUCS<-colnames(biospecimenUCS)
cnBSCESC<-colnames(biospecimenCESC)

setdiff(cnClUCS,cnClCESC)
all.equal(cnClCESC,cnClUCS)
setdiff(cnBSCESC,cnBSUCS)
head(biospecimenCESC[,setdiff(cnBSCESC,cnBSUCS)])

#join two bstables
bsJ<-bind_rows(biospecimenCESC,biospecimenUCS)

ucsJ<-getjoinedBiospcClinc("TCGA-UCS")
cescJ<-getjoinedBiospcClinc("TCGA-CESC")

cnCESC<-colnames(cescJ)
cnUCS<-colnames(ucsJ)
extraCols<-setdiff(cnCESC,cnUCS)

temp<-cescJ[,setdiff(cnCESC,cnUCS)]

joined<-bind_rows(ucsJ,cescJ)


##Download only BRCA metadata
brcaDF<-getjoinedBiospcClinc("TCGA-BRCA")
brcaDF<-brcaDF[,colsToKeep]
#remove cols with all NA values
naCols<-colnames(brcaDF)[sapply(brcaDF, function(x)all(is.na(x)))]
brcaDFNONA<-brcaDF[,!(colnames(brcaDF) %in% naCols)]
#keep rows with RNA samples only
brcaDFRNA<-brcaDFNONA%>%filter(portions.analytes.analyte_type_id == "R")

write.csv(brcaDFNONA,"TCGAbrcaMetadata_reduced.csv",row.names = F)

##download gene mutation metadata
brcaMAF <- GDCquery_Maf("BRCA", pipelines = "varscan2")
normalSampList<-brcaDF%>% filter(sample_type=="Solid Tissue Normal") %>% select(portions.analytes.aliquots.submitter_id)

maf_norm<-brcaMAF %>% filter(Tumor_Sample_Barcode %in% normalSampList$portions.analytes.aliquots.submitter_id)
#visualise mutation



brcaClinic<-GDCquery_clinic(project = "TCGA-BRCA", type = "Clinical")
#maf<-read.maf(brcaMAF,clinicalData =brcaClinic,isTCGA = T)
maf<-read.maf(brcaMAF,isTCGA = T)

plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,showBarcodes=F)
#We will draw oncoplots for top ten mutated genes.
oncoplot(maf = maf, top = 10, fontSize = 12)

geneCloud(input = maf, minMut = 30)

#identify cols to keep
write.csv(colnames(brcaDF),"TCGAcolnames.csv")



