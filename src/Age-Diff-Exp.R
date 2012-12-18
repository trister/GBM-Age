### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20120606

### Start looking at the differential gene expression in the TCGA GBM data based on the age differences
###


require(synapseClient)
require(ggplot2)
require(gplots)
require(limma)
require(survival)
require(HDclassif)
library(org.Hs.eg.db)
source("./ggkm.R")

synapseLogin()

## PULL IN GBM DATA BY SOURCING "intermediate data a" ENTITY FROM SYNAPSE - output from populateGBMdata function
print("Loading data from Synapse")
dataReturn <- loadEntity(299091)
gbmPat <- dataReturn$objects$gbmPat
gbmClin <- dataReturn$objects$gbmClin
gbmMat <- dataReturn$objects$gbmMat
tmpSurv <- dataReturn$objects$tmpSurv



#Now let's get an object to have geneIDs from Entrez
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the SYMBOL for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}


#####
## ASSOCIATION OF EXPRESSION WITH SURVIVAL
#####
print("Building Cox models and assessing transcript significance")
plot(survfit(tmpSurv ~ 1))

coxRes <- apply(gbmMat, 1, function(x){
  summary(coxph(tmpSurv ~ x))$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
})

## P-VALUE HISTOGRAMS
print("Generating p-value histograms")
pvalHist <- qplot(coxRes[2, ], geom = "histogram") + 
  opts(title = "GBM Transcripts and Survival Unadjusted p-values")
adjPvals <- p.adjust(coxRes[2, ], method = "BH")
adjPvalHist <- qplot(adjPvals, geom = "histogram") + 
  opts(title = "GBM Transcripts and Survival B-H Adjusted p-values")

## VOLCANO PLOT
print("Generating volcano plot")
vplotDF <- as.data.frame(t(rbind(log2(coxRes[1, ]), -1*log10(coxRes[2, ]))))
colnames(vplotDF) <- (c("Column1", "Column2"))

volcanoPlot <- ggplot(vplotDF, aes(Column1, Column2)) + geom_point() +
  opts(title = "Volcano Plot GBM Transcripts") +
  ylab("- log10 P-values") +
  xlab("Transcripts") +
  opts(plot.title = theme_text(size = 14))





####
#Look at a linear model fit with limma for age
#

#First get rid of the patients with NA for the age

id <- which(!is.na(gbmPat$age_at_initial_pathologic_diagnosis))
gbmMat.culled <- gbmMat[,id]
gbmPat.culled <- gbmPat[id,]

temp <- quantile(gbmPat.culled$age_at_initial_pathologic_diagnosis,na.rm=T)

age.group <- sapply(gbmPat.culled$age_at_initial_pathologic_diagnosis, function(x) {
  if(!is.na(x)) {
  if (x < temp[2] & x > 25) 1
  else if (x > temp[4]) 3
  else 2
  }
  else NA
})



color.map <- function(x) {
  if (x==1) "#FF0000"
  else if(x==2) "#0000FF"
  else if(x==3) "#FF00FF"
  else "#00FF00"
}



design <- model.matrix(~0+factor(age.group))
#design <- model.matrix(~age.group)
colnames(design) = c("low", "middle", "high")

fit <- lmFit(gbmMat.culled, design)
fit <- eBayes(fit)
topTable(fit, coef=2)

contrast.matrix <- makeContrasts(low-high, low-middle, levels=design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=2, adjust="BH",lfc=0.5)
results <- decideTests(fit2)
table(results)



temp <- length(gbmPat.culled$age_at_initial_pathologic_diagnosis)
colors <- colorpanel(temp,"#0000FF","#FF0000")
patientcolors <- colors[rank(gbmPat.culled$age_at_initial_pathologic_diagnosis,ties.method="random")]

#gbmPat.selected <- gbmPat.culled[which(age.group==1|age.group==3),]
gbmMat.selected <- gbmMat.culled[,which(age.group==1|age.group==3)]
patientcolors <- patientcolors[which(age.group==1|age.group==3)]

##
#reorder all of these by age
temp <- order(gbmPat.selected$age_at_initial_pathologic_diagnosis)
gbmMat.selected <- gbmMat.selected[,temp]
gbmPat.selected <- gbmPat.selected[temp,]
patientcolors <- patientcolors[temp]

vennDiagram(results)

#here we will try to use ROAST
ind <- as.numeric(row.names(topTable(fit2,coef=1,number=50,p.value=0.05)))
ind <- which(results[,1]!=0)
roast(ind,gbmMat.culled,design,c(1,0,0))



esetSel <- gbmMat.selected[ind,]
esetSel <- gbmMat.selected[results[,1]!=0,]
esetSel <- gbmMat.selected[results[,1]==-1,]
esetSel <- gbmMat.selected[results[,1]!=0&results[,2]!=0,]




gene.names <- lapply(rownames(esetSel),function(x){
  return(xx[strsplit(x,"_mt")[[1]]])
})
paste(unlist(gene.names),collapse=" ")



##Survival with clustering based on these genes
full.cluster <- gbmMat[results[,1]!=0,]
bigClust <- hddc(t(full.cluster))
patientcolors <- unlist(lapply(bigClust$class, color.map))

heatmap(full.cluster, col=topo.colors(100),scale="none",ColSideColors=patientcolors,main="Full expression with clusters")
gbmSurvFit <- survfit(tmpSurv ~ bigClust$class)
#ggkm(gbmSurvFit)

gbmStrata <- bigClust$class
gbmKMPlot <- ggkm(gbmSurvFit, 
                  ystratalabs = (c("Mid", "Young", "Old")), 
                  timeby = 365,
                  main = "GBM K-M Plot By Class")

#which classes do each of the different groups fall into

smallClust <- unlist(lapply(colnames(gbmMat.selected),function(x){
  bigClust$class[grep(x,colnames(gbmMat))]
}))





ggkm(survfit(tmpSurv~bigClust$class=="1"|bigClust$class=="4"),timeby=365,ystratalabs=c("2 or 3","1 or 4"),main="GBM K-M plot by Class")
ggkm(survfit(tmpSurv ~ gbmMat["800_mt",]>median(gbmMat["800_mt",])),timeby=365)


patientcolors <- unlist(lapply(smallClust, color.map))

heatmap(esetSel, col=topo.colors(100), scale = "none",main="Differential Expression by Age")
heatmap(esetSel, col=topo.colors(100), scale = "none", Colv=NA,ColSideColors=patientcolors,main="Differential Expression by Age")




####
####
#### Make some nice graphs for presentations
####
####
p <- ggplot(gbmPat.culled)
p + geom_histogram(aes(x=age_at_initial_pathologic_diagnosis,fill= ..count..),binwidth=3) +
    opts(title="TCGA Age at Diagnosis",labels=list(x= "Age (years)")) +
    scale_fill_gradient("Count",low="red",high="green")

n <- ggplot(gbmPat)
n+geom_density(aes(x=age_at_initial_pathologic_diagnosis, fill=factor(bigClust$class)),alpha=0.2) +
  scale_fill_hue(name="Class",
                 breaks=c(2, 1, 3),
                 labels=c("Young", "Mid", "Old")) +
  opts(title="Age by Cluster",labels=list(x= "Age (years)"))




temp <- quantile(gbmPat$age_at_initial_pathologic_diagnosis,na.rm=T)
full.ages <- sapply(gbmPat$age_at_initial_pathologic_diagnosis, function(x) {
  if(!is.na(x)) {
    if (x < temp[2]) 1
    else if (x < temp[3]) 2
    else if (x < temp[4]) 3
    else 4
  }
  else NA
})


gbmSurvFit <- survfit(tmpSurv ~ full.ages)

gbmStrata <- full.ages
gbmKMPlot <- ggkm(gbmSurvFit, 
                  ystratalabs = (c("25-49", "50-59", "60-68",">69")), 
                  timeby = 365,
                  main = "GBM K-M Plot By Age")






###
# Now let's do something with the Rembrandt Data
###

## PULL IN REMBRANDT data from synapse
print("Loading data from Synapse")
rembrandt.dataReturn <- loadEntity(376920)
rembrandt.matrix <- rembrandt.dataReturn$objects$rembrandt.matrix.matched
rembrandt.phenotype <- rembrandt.dataReturn$objects$phenotype.data.matched

rembrandt.phenotype$Survival..months.[which(rembrandt.phenotype$Survival..months.=="--")] <- NA

rembrandt.select <- rembrandt.matrix[which(rownames(rembrandt.matrix) %in% rownames(esetSel)),]

rembrandt.patients <- sapply(rembrandt.phenotype$Disease, function(x) {
  if (x==" ASTROCYTOMA") 1
  else if (x==" GBM") 3
  else if (x==" OLIGODENDROGLIOMA") 2
  else 4
})

rembrandt.patients <- unlist(lapply(rembrandt.patients, color.map))
  
heatmap(rembrandt.select[,which(rembrandt.phenotype$Disease==" GBM")], col=topo.colors(100),scale = "none", main="REMBRANDT Expression")

rembrandt.cluster <- hddc(t(rembrandt.select[,which(rembrandt.phenotype$Disease==" GBM")]))
  
rembrandtSurvFit <- survfit(Surv(as.numeric(rembrandt.phenotype$Survival..months.[which(rembrandt.phenotype$Disease==" GBM")])) ~ rembrandt.cluster$class)
rembrandtStrata <- rembrandt.cluster$class
rembrandtKMPlot <- ggkm(rembrandtSurvFit, 
                  ystratalabs = (c("Mid", "Old", "Young")), 
                  timeby = 12,
                  main = "GBM K-M Plot By Class")









######
######
######
### Now let's look at normals from the Harvard Brain Registry
#####
#####
#####









#let's look at MDK - just one of a few different genes that could be interesting
plot(gbmPat.culled$age_at_initial_pathologic_diagnosis,gbmMat.culled[5043,])
boxplot(gbmMat.culled[5043,which(gbmPat.culled$age_at_initial_pathologic_diagnosis<50)],
        gbmMat.culled[5043,which(gbmPat.culled$age_at_initial_pathologic_diagnosis>65)],
        xlab="age-group",ylab="expression",main="Expression of MDK")

t.test(gbmMat.culled[5043,which(gbmPat.culled$age_at_initial_pathologic_diagnosis<50)],
       gbmMat.culled[5043,which(gbmPat.culled$age_at_initial_pathologic_diagnosis>65)])
coxph(tmpSurv~gbmMat[5043,])

mfit <- survfit(tmpSurv~gbmMat[5043,]>median(gbmMat[5043,]))
plot(mfit,col=c("black","grey"),lty=1:2,main="KM of MDK expression",xlab="Days",ylab="Probability of Survival")
survdiff(tmpSurv~gbmMat[5043,]>median(gbmMat[5043,]))





### Load in the Level 3 CNV data
TCGA.directory <- "~/TCGA/Copy_Number_Results/MSKCC__HG-CGH-244A/Level_3/"
TCGA.files <- list.files(TCGA.directory,full.names=TRUE)
















populateGBMdata <- function(){
  print("Loading require packages")
  require(synapseClient)
  require(Biobase)
  require(survival)
  
  
  #####
  ## EXPRESSION DATA
  #####
  print("Loading molecular data from Synapse")
  ge <- loadEntity("274865")
  gbmExprSet <- ge$objects$coherentEset
  gbmMat <- exprs(gbmExprSet)
  
  ## GET THE ROW NAMES FROM ANOTHER LAYER
  ne <- loadEntity("274544")
  tmpNames <- featureNames(ne$objects$eset)
  
  ## SUBSET TO THOSE WITHOUT NAS
  print("Processing molecular data")
  id <- which(colSums(is.na(gbmMat)) != nrow(gbmMat))
  gbmMat <- gbmMat <- gbmMat[, id]
  rownames(gbmMat) <- tmpNames
  
  ## GET RID OF NON-TUMOR SAMPLES
  theseTissues <- sapply(strsplit(colnames(gbmMat), "-", fixed=T), function(x){
    substr(x[4], 1, 2)
  })
  gbmMat <- gbmMat[ , theseTissues == "01" ]
  
  ## LOOK AT PATIENT LEVEL AND REMOVE DUPS
  thesePats <- sapply(strsplit(colnames(gbmMat), "-", fixed=T), function(x){
    paste(x[1:3], collapse="-")
  })
  idd <- duplicated(thesePats)
  gbmMat <- gbmMat[ , !idd]
  thesePats <- thesePats[!idd]
  
  
  #####
  ## CLINICAL DATA
  #####
  print("Loading clinical data from Synapse")
  gc <- loadEntity("274426")
  gbmClin <- gc$objects$clinAll
  gbmPat <- gbmClin$clinical_patient_public_gbm
  rownames(gbmPat) <- gbmPat$bcr_patient_barcode
  gbmPat <- gbmPat[ thesePats, ]
  
  print("Creating survival object")
  gbmPat$vital_status[ gbmPat$vital_status == "[Not Available]" ] <- NA
  gbmPat$survTime <- as.numeric(gbmPat$days_to_death)
  gbmPat$surv <- ifelse( gbmPat$vital_status=="DECEASED", 1, 0)
  gbmPat$survTime[ which(gbmPat$vital_status == "LIVING") ] <- as.numeric(gbmPat$days_to_last_followup)[ which(gbmPat$vital_status == "LIVING") ]
  tmpSurv <- Surv(gbmPat$survTime, gbmPat$surv)
  
  
  print("Generating output list")
  rm(list=setdiff(ls(), c("gbmPat", "gbmClin", "gbmMat", "tmpSurv")))
  
  return(list("gbmPat" = gbmPat,
              "gbmClin" = gbmClin,
              "gbmMat" = gbmMat,
              "tmpSurv" = tmpSurv))
}

