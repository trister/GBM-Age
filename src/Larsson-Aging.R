require(synapseClient)
require(affy)
require(mg.hgu133a.db)

#########################
# Install dependencies  #
#########################
#source('http://sage.fhcrc.org/CRAN.R')
#pkgInstall("syanpseClient")
#pkgInstall("mg.hgu133plus2.db", type="source")   #This has to be done better somehow
#pkgInstall("mg.hgu133a.db", type="source")

loadData <- function(exprEntity, metaEntity){
  # Extracts an expression entitity and a TCGA metadata entity from Synapse
  #
  # Args:
  #   exprEntity: a synapse ID pointing to a expression data set
  #   metaEntity: a synapse ID pointing to a metadata object as defined under TCGA
  #
  # Return:
  #   list of expresion set and filtered metadata
  ent <- loadEntity(exprEntity)     #load expression entity
  eset <- exprs(ent$objects$eset)   #Extract expression object and convert to matrix
  
  metadata <- loadEntity(metaEntity)
  metadata <- metadata$objects$metadata
  
  # Find the column containing hthgu133a file names
  columnMatchTotals <- apply(metadata,2,function(x){ sum(!is.na(match(colnames(eset), x)))})
  columnId <- which(columnMatchTotals > 0)
  
  #Extract rows corresponding to hthgu133a
  metadata <- metadata[which(!is.na(metadata[,columnId])),]
  rownames(metadata) <- metadata[,columnId]
  #Extract columns containing data for this subdata type.
  metadata <- metadata[,which(colSums(!is.na(metadata)) > 0)]
  
  #Filter out un-annotated expression sets.
  inCommon <- intersect(colnames(eset), rownames(metadata))
  return(list(eset=eset[,inCommon], metadata=metadata[inCommon,]))
}

diffExpression <- function(y, eset){
  # Given response variable returns linear association between each
  # Args:
  #   y: Phenotypic variable representing age.
  #   eset: Expression data matrix
  
  removeIdx= is.na(y) | (y>50 & y<65) | y<25
  y<-y[!removeIdx]
  eset<-eset[,!removeIdx]
  
  vals<-apply(eset, 1, function(eset) {
    t.test(eset[y<=50], eset[y>=65])$p.value })
  return(sort(vals))
}

hif1OverlapScore <- function(exprEntity, metaEntity, hifGenes) {
  #Loads data, performs ttest in data versus age, and looks for enrichment
  # compared to hif1 targets.  Returns p-value for enrichment
  
  data <- loadData(exprEntity, metaEntity)
  age <- data$metadata["age_at_initial_pathologic_diagnosis"]
  pvals <- diffExpression(age, data$eset)
  nSigP <- sum(p.adjust(pvals, method="BH")<0.05)
  
  #Load R annotation library for hgu133a platform
  symbols <- toTable(mg.hgu133aSYMBOL[names(pvals[1:nSigP])])
  
  #Merge with symbols and find intersections with hif1 targets
  overlapping <- intersect(symbols[,2], hif1Genes)
  p <- fisher.test(matrix(c(length(overlapping), length(hif1Genes)-length(overlapping),
                            nSigP-length(overlapping), dim(data$eset)[1]-length(hif1Genes)), 2,2))$p.value
  print("Overlapping genes:")
  print(overlapping)
  print(sprintf("Enrichment P value = %0.2g", p))
  
  # #Plot top genes
  # layout(matrix(1:6, 2,3, byrow=TRUE))
  # for (i in 1:6) {
  #   plot(cbind(age, eset[names(x)[i],]), ylab=names(x)[i])
  # }
  return(p)
}


################################
#Main method
################################
 synapseLogin()
setwd("~/ROSTO/")
ids <-matrix(c('syn313583', 'syn372761',  #GBM_hghgu133a, GBM_mergedClinical
               'syn313579', 'syn376054',  #LUSC_hthgu133a
               'syn313584', 'syn376070'),  #OV_hthgu133a
             ncol=2, byrow=TRUE)

load("hif1Genes.Rda")
for (i in 1:3) {
  hif1OverlapScore(ids[i,1], ids[i,2], hif1Genes)
}