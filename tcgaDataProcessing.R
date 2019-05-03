library(TCGAbiolinks)
library(dplyr)
library(DT)
library(data.table)
library(plyr)
packageVersion("TCGAbiolinks")


query <- GDCquery(project = "TCGA-CHOL",  data.category = "Clinical", file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

queryB <- GDCquery(project = "TCGA-COAD",  data.category = "Biospecimen", file.type = "xml")
queryB <- GDCquery(project = "TCGA-BRCA",  data.category = "Biospecimen", file.type = "xml")
GDCdownload(queryB,method = "client")

aliquot <- GDCprepare_clinic(queryB, clinical.info = c("aliquot"))
sample <- GDCprepare_clinic(queryB, clinical.info = c("sample"))
bio_patient <- GDCprepare_clinic(queryB, clinical.info = c("bio_patient"))
analyte <- GDCprepare_clinic(queryB, clinical.info = c("analyte"))
portion <- GDCprepare_clinic(queryB, clinical.info = c("portion"))
protocol <- GDCprepare_clinic(queryB, clinical.info = c("protocol"))
slide <- GDCprepare_clinic(queryB, clinical.info = c("slide"))

hbp<-as.data.frame(names(bio_patient))
hanalyte<-as.data.frame(names(analyte))
hportion<-as.data.frame(names(portion))
hprot<-as.data.frame(names(protocol))
hslid<-as.data.frame(names(slide))


join(aliquot,bio_patient)

length(unique(aliquot$bcr_patient_barcode))
length(unique(sample$bcr_patient_barcode))
length(unique(bio_patient$bcr_patient_barcode))
length(unique(analyte$bcr_patient_barcode))
length(unique(portion$bcr_patient_barcode))
length(unique(protocol$bcr_patient_barcode))
length(unique(slide$bcr_patient_barcode))
length(Reduce(intersect,list(aliquot$bcr_patient_barcode,sample$bcr_patient_barcode,bio_patient$bcr_patient_barcode,analyte$bcr_patient_barcode,portion$bcr_patient_barcode,protocol$bcr_patient_barcode,slide$bcr_patient_barcode)))

#remove duplicated rows from all tables
bio_patient_nr<-distinct(bio_patient)
aliquot_nr<-distinct(aliquot)
samp_nr<-distinct(sample)
analyte_nr<-distinct(analyte)
portion_nr<-distinct(portion)
protocol_nr<-distinct(protocol)
slide_nr<-distinct(slide)

#1 join aliquot with analyte, add analyte barcode in aliquot. for aliquot barcode TCGA-3L-AA1B-01A-01D-YYYY-23, analyte barcode is TCGA-3L-AA1B-01A-01D
#more info on barcodes https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
#tcga center codes https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes

aliquot_nr<-aliquot_nr%>% mutate(bcr_analyte_barcode=substr(bcr_aliquot_barcode,1,nchar(as.character(bcr_aliquot_barcode))-8))
j1<-join(analyte_nr,aliquot_nr,by="bcr_analyte_barcode")

#2 add bioportion barcode to j1 and join with bioportion
#portion: TCGA-3L-AA1B-01A-11 ; analyte: TCGA-3L-AA1B-01A-11
j1<-j1 %>% mutate(bcr_portion_barcode=substr(bcr_analyte_barcode,1,nchar(as.character(bcr_analyte_barcode))-1))
j2<-join(j1,portion_nr,by="bcr_portion_barcode")

#3 add biosample barcode to j2
#sample: TCGA-3L-AA1B-01A ; portion: TCGA-3L-AA1B-01A-11
j2<-j2 %>% mutate(bcr_sample_barcode=substr(bcr_portion_barcode,1,nchar(as.character(bcr_portion_barcode))-3))
j3<-join(j2,samp_nr,by="bcr_sample_barcode")

#finally join by biopatient
j4<-join(j3,bio_patient_nr,by="bcr_patient_barcode")


#download clinical data
queryC <- GDCquery(project = "TCGA-COAD",  data.category = "Clinical", file.type = "xml")
GDCdownload(queryC)
drug<-GDCprepare_clinic(queryC, clinical.info = "drug")
admin<-GDCprepare_clinic(queryC, clinical.info = "admin")
follow_up<-GDCprepare_clinic(queryC, clinical.info = "follow_up")
radiation<-GDCprepare_clinic(queryC, clinical.info = "radiation")
patient<-GDCprepare_clinic(queryC, clinical.info = "patient")
stage_event<-GDCprepare_clinic(queryC, clinical.info = "stage_event")
new_tumor_event<-GDCprepare_clinic(queryC, clinical.info = "new_tumor_event")


admin_nr<-distinct(admin)
patient_nr<-distinct(patient)

j5<-join(j4,patient_nr,by="bcr_patient_barcode")




#keep only selected columns



#other method
clinicalBRCA <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
biospecimenBRCA <- GDCquery_clinic(project = "TCGA-BRCA", type = "Biospecimen")
biospecimenCOAD <- GDCquery_clinic(project = "TCGA-COAD", type = "Biospecimen")

sampdf<-head(biospecimenCOAD,10)[,c("sample_type_id","tumor_code_id","sample_id","submitter_id","portions")]
sampPor<-as.data.frame(sampdf[2,  c("portions")])

#use apply to unlist columns
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
    
    newRow<-cbind(thisRow,tempdf)
    res<-bind_rows(res,newRow)
    #for(j in 1: dim(tempdf)[1]){
      #convert to dataframe in case there is only single column
     #res<-bind_rows(res,newRow)
    #}
  }
  #print(res)
  return(res)
  
}
res<-NULL
sampdfExbrnew<-expand(biospecimenBRCA,"portions")

sampdfExanalyte<-expand(sampdfExbrnew,"portions.analytes")

sampdfExaliquot<-expand(sampdfExanalyte,"portions.analytes.aliquots")

#add patient barcode to biospecimen data
brcaTabRNA<- brcaTabRNA %>% mutate(bcr_patient_barcode=substr(submitter_id,1,nchar(as.character(submitter_id))-4))

#join clinical and biospecimen

brcaJoinedRNA<-join(clinicalBRCA,brcaTabRNA,by="bcr_patient_barcode")


##Download only BRCA metadata
brcaDF<-getjoinedBiospcClinc("TCGA-BRCA")

#remove cols with all NA values
naCols<-colnames(brcaDF)[sapply(brcaDF, function(x)all(is.na(x)))]
brcaDFNONA<-brcaDF[,!(colnames(brcaDF) %in% naCols)]
#keep rows with RNA samples only
brcaDFRNA<-brcaDFNONA%>%filter(portions.analytes.analyte_type_id == "R")

##download gene mutation metadata
brcaMAF <- GDCquery_Maf("BRCA", pipelines = "varscan2")

#visualise mutation
library(maftools)
brcaClinic<-GDCquery_clinic(project = "TCGA-BRCA", type = "Clinical")
#maf<-read.maf(brcaMAF,clinicalData =brcaClinic,isTCGA = T)
maf<-read.maf(brcaMAF,isTCGA = T)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#We will draw oncoplots for top ten mutated genes.
oncoplot(maf = maf, top = 10, fontSize = 12)

geneCloud(input = maf, minMut = 30)
