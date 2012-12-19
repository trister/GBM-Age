### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20121112

### Build an elasticnet model trained on POSTN class in TCGA to learn what other genes may be implicated
### and verify in REMBRANDT


require(synapseClient)
require(glmnet)
require(randomForest)
require(caret)
require(affy)
require(survival)
require(pROC)
require(pls)
require(gplots)
require(mixtools)
require(ggplot2)
source("./src/ggkm.R")
library(org.Hs.eg.db)






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




## let's make a survival object
metadata$vital_status[metadata$vital_status=="[Not Available]"] <- NA
gbmPat <- c()
gbmPat$survTime <- as.numeric(metadata$days_to_death)
gbmPat$surv <- ifelse(metadata$vital_status=="DECEASED", 1,0)
gbmPat$survTime[ which(metadata$vital_status=="LIVING")] <- as.numeric(metadata$days_to_last_followup)[ which(metadata$vital_status=="LIVING")]

tmpSurv <- Surv(gbmPat$survTime,gbmPat$surv)


## let's make a survival object for REMBRANDT
gbmPat.rembrandt <- c()
gbmPat.rembrandt$survTime <- as.numeric(rembrandtPat$Survival..months.)
gbmPat.rembrandt$surv <- ifelse(rembrandtPat$Survival..months.=="--", 0,1)
gbmPat.rembrandt$survTime[ which(gbmPat.rembrandt$surv==0)] <- as.numeric(rembrandtPat$Followup.Month)[ which(gbmPat.rembrandt$surv==0)]

tmpSurv.rembrandt <- Surv(gbmPat.rembrandt$survTime,gbmPat.rembrandt$surv)

#####
## ASSOCIATION OF EXPRESSION WITH SURVIVAL
#####
print("Building Cox models and assessing transcript significance")
plot(survfit(tmpSurv ~ 1))

coxRes <- apply(eset, 1, function(x){
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


adjPvalHist
volcanoPlot





#cull the REMBRANDT U133Aplus2 to U133A genesets
temp <- unlist(lapply(rownames(rembrandtEset),function(x){
  return(gsub("_mt","_eg",x))
}))

rembrandtEset.culled <- rembrandtEset
rownames(rembrandtEset.culled) <- temp

temp <- intersect(temp,rownames(eset))

eset <- eset[temp,]
rembrandtEset.culled <- rembrandtEset.culled[temp,]

eset.culled <- eset[-grep("10631_eg",rownames(eset)),]
rembrandtEset.culled <- rembrandtEset.culled[-grep("10631_eg",rownames(rembrandtEset.culled)),]




#First get rid of the patients with NA for the age

id <- which(!is.na(metadata$age_at_initial_pathologic_diagnosis))
eset.culled <- eset.culled[,id]
metadata.culled <- metadata[id,]



# first rescale the data (chemores) to have the same mean and variance than the LUAD
# Justin's function to rescale the VS to get the same mean/var than the TS
normalize_to_X <- function(mean.x, sd.x, Y){
  m.y <- rowMeans(Y)
  sd.y <- apply(Y, 1, sd)
  Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
  Y.adj[sd.y == 0] <- mean.x[sd.y==0]
  Y.adj
}

tmp <- apply(eset.culled,1,sd)
rembrandtEset.scaled <- normalize_to_X(rowMeans(eset.culled),tmp,rembrandtEset.culled)

#get rid of the least variant probes
tmp1 <- which(tmp>quantile(tmp,probs=0.2))
eset.culled2 <- eset.culled[tmp1,]
rembrandtEset.scaled <- rembrandtEset.scaled[tmp1,]
rm(tmp,tmp1)



## begin by looking at a model of age with elastic net
cv.fit <- cv.glmnet(x=t(eset.culled2),y=metadata.culled$age_at_initial_pathologic_diagnosis,nfolds=10,alpha=0.1,family="gaussian")
plot(cv.fit)
fitEnet <- glmnet(x=t(eset.culled2), y=metadata.culled$age_at_initial_pathologic_diagnosis, family="gaussian", alpha=.1, lambda=cv.fit$lambda.min)


###############################################################################
#  Look now at all patients in REMBRANDT                                      #
###############################################################################


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled), type="response", s="lambda.min")  #all of REMBRANDT

boxplot(yhatEnet ~ as.factor(rembrandtPat$Age.at.Dx..years.), ylab="Predicted Age", xlab="age group", main="elastic net validation")
stripchart(yhatEnet ~ as.factor(rembrandtPat$Age.at.Dx..years.), ,pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.factor(rembrandtPat$Age.at.Dx..years.),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)


plot(survfit(tmpSurv.rembrandt  ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.rembrandt  ~ riskEnet, rho=0)









###############################################################################
#  Look now at only the patients without an age                               #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Age.at.Dx..years.==" --")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Age.at.Dx..years.==" --")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Age.at.Dx..years.==" --")]), type="response", s="lambda.min")

#boxplot(yhatEnet ~ as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Age.at.Dx..years.==" --")]), ylab="Age", xlab="predicted Age group", main="elastic net validation")
#stripchart(yhatEnet ~ as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Age.at.Dx..years.==" --")]),pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

#rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Age.at.Dx..years.==" --")]),ci=TRUE)
#plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)



#ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
#summary(survfit(tmpSurv.lowgrade  ~ riskEnet))
tmp1 <- quantile(yhatEnet)
tmp <- unlist(lapply(yhatEnet,function(x){
  if(x<tmp1[2]) 1
  else if(x<tmp1[3]) 1
  else if(x<tmp1[4]) 1
  else 2
}))

plot(survfit(tmpSurv.lowgrade ~ tmp), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
ggkm(survfit(tmpSurv.lowgrade  ~ tmp),timeby=12,main="KM of REMBRANDT",ystratalabs=c("Younger","Older"))
survdiff(tmpSurv.lowgrade ~ tmp, rho=0)








###############################################################################
#  Look now at only the Grade IV                                              #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" IV")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" IV")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Grade==" IV")]), type="response", s="lambda.min")

boxplot(yhatEnet ~ as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Grade==" IV")]), ylab="Age", xlab="predicted Age group", main="elastic net validation")
stripchart(yhatEnet ~ as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Grade==" IV")]),pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Grade==" IV")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)



#ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
#summary(survfit(tmpSurv.lowgrade  ~ riskEnet))

plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)



###############################################################################
#  Look now at only the Grade III                                             #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" III")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" III")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Grade==" III")]), type="response", s="lambda.min")

boxplot(yhatEnet ~ as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Grade==" III")]), ylab="Age", xlab="predicted Age group", main="elastic net validation")
stripchart(yhatEnet ~ as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Grade==" III")]),pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Grade==" III")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)



#ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
#summary(survfit(tmpSurv.lowgrade  ~ riskEnet))

plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)










###############################################################################
#  Look now at only the Grade II                                              #
###############################################################################

lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" II")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" II")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.scaled[,which(rembrandtPat$Grade==" II")]), type="response", s="lambda.min")

boxplot(yhatEnet ~ as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Grade==" II")]), ylab="Age", xlab="predicted Age group", main="elastic net validation")
stripchart(yhatEnet ~ as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Grade==" II")]),pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.factor(rembrandtPat$Age.at.Dx..years.[which(rembrandtPat$Grade==" II")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)



#ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
#summary(survfit(tmpSurv.lowgrade  ~ riskEnet))

plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)













heatmap(x=rembrandtEset.culled[which(abs(fitEnet$beta)>0),],
        ColSideColors=c("red","green","blue","yellow")[factor(rembrandtPat$Disease)],
        scale="row",col=redgreen(50),main="Elastic Net")



gene.names <- lapply(rownames(eset)[which(abs(fitEnet$beta)>0)],function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names),collapse=" ")







### Build a more robust classifer using bootstrapping




# alpha = 0.1 optimize lambda nfolds=10
# bootstrap 100 times
N <- 100
fit <- c()
selected <- rep(0,length(rownames(eset.culled2)))
features <- c()
yhat_REMBRANDT <- c()
models <- 0
i <- 0
for(i in 1:N) {
  
  j <- sample(1:length(metadata.culled$age_at_initial_pathologic_diagnosis),replace=TRUE)
  
  cv.fit <- cv.glmnet(x=t(eset.culled2[,j]),y=metadata.culled$age_at_initial_pathologic_diagnosis[j],nfolds=10,alpha=0.1,family="gaussian")
  fit <- glmnet(x=t(eset.culled2[,j]), y=metadata.culled$age_at_initial_pathologic_diagnosis[j], family="gaussian", alpha=.1, lambda=cv.fit$lambda.min)
  features <- c(features,fit$beta)
  print(i)
}





# weighted aggregation
#function from In Sock to help with the determination of the best features to select
weightAggregation<-function(resultsModel){
  
  ResultBS<-c()
  for(k in 1:length(resultsModel)){  
    #     a<-abs(resultsModel[[k]][-1])
    #     A<-sort(a,decreasing = T,index.return=T)
    ResultBS<-cbind(ResultBS,rank(abs(as.numeric(resultsModel[[k]][-1])),ties.method="min")/length(resultsModel[[k]]))
  }
  rownames(ResultBS)<-rownames(resultsModel[[1]])[-1]  
  reference <- apply(ResultBS,1,sum) 
  return(reference)
}


selected <- weightAggregation(features)
plot(sort(selected))

#use the weighted aggregation to grab the top 100
select_features <- sort(selected,decreasing=TRUE)[1:100]


heatmap(x=rembrandtEset.scaled[names(sort(select_features)[90:100]),],
        ColSideColors=c("red","green","blue","yellow")[factor(rembrandtPat$Disease)],
        scale="none",col=redgreen(50),main="Elastic Net")

gene.names <- lapply(names(sort(select_features)[90:100]),function(x){
  return(xx[strsplit(x,"_eg")[[1]]])
})

paste(unlist(gene.names),collapse=" ")


###############################################################################
#  Look now at only the Grade III                                             #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" III")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" III")]


tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.culled[,which(rembrandtPat$Grade==" III")]), type="response", s="lambda.min")

boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")], ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade==" III")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= median(yhatEnet), 1, 0))
names(riskEnet) <- rownames(yhatEnet)



ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
summary(survfit(tmpSurv.lowgrade  ~ riskEnet))

plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)










###############################################################################
#  Look now at only the Grade II                                              #
###############################################################################


lowgradePat <- c()
lowgradePat$surv <- gbmPat.rembrandt$surv[which(rembrandtPat$Grade==" II")]
lowgradePat$survTime <- gbmPat.rembrandt$survTime[which(rembrandtPat$Grade==" II")]

tmpSurv.lowgrade <- Surv(lowgradePat$survTime,lowgradePat$surv)


yhatEnet <- predict(fitEnet, t(rembrandtEset.culled[,which(rembrandtPat$Grade==" II")]), type="response", s="lambda.min") 

boxplot(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")], ylab="3-year OS prediction (%)", xlab="predicted POSTN group", main="elastic net validation")
stripchart(yhatEnet ~ predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")],pch=20, col="royalblue", vertical=TRUE, add=TRUE, cex=.6)

rocEnet <- roc(predictor=as.numeric(yhatEnet), response=as.numeric(predictedClasses.rembrandt[which(rembrandtPat$Grade==" II")]),ci=TRUE)
plot.roc(rocEnet,col="red")

riskEnet <- as.vector(ifelse(yhatEnet >= 0.15, 1, 0))
names(riskEnet) <- rownames(yhatEnet)



ggkm(survfit(tmpSurv.lowgrade  ~ riskEnet),timeby=12,main="KM of REMBRANDT")
summary(survfit(tmpSurv.lowgrade  ~ riskEnet))


plot(survfit(tmpSurv.lowgrade ~ riskEnet), main="elastic net model", xlab="months",ylab="probability of OS",col= c("blue","magenta"),lwd=3)
survdiff(tmpSurv.lowgrade ~ riskEnet, rho=0)




