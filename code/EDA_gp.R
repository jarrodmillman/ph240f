####Group_Project####
setwd("/Users/lktroszak/src/ph240f/data")
load("eset.rda")
#Check to make sure everything is as it should be
sum(eset$type=="KIRC")
sum(eset$type=="KIRP")
sum(eset$train)
sum(eset$train & eset$type=="KIRC")
sum(eset$train & eset$type=="KIRP")

#total rsem counts for each subject
colSums(exprs(eset))
#proportion of 0 counts for each subject
prop0=colSums(exprs(eset)==0)/dim(eset)[1]
max(prop0)
min(prop0)


#Filter out some more low count genes
a=rowMeans(exprs(eset)[,eset$type=="KIRC"])>=log(10+1)
b=rowMeans(exprs(eset)[,eset$type=="KIRP"])>=log(10+1)
f=(a|b)
new_eset=eset[f]
hist(exprs(eset))
hist(exprs(new_eset))

boxplot(exprs(new_eset),main="Boxplots of Log-Expression for KIRC & KIRP Subjects")
#look at Mean-Difference Plots from subjects within cell-types and between cell-types
#might need to untransform eset from log since I think MDPlot takes the log for you
MDPlot(new_eset,c(1,5))
MDPlot(new_eset,c(48,60))
MDPlot(new_eset,c(12,50))
MDPlot(new_eset,c(12,70))
#looks like we may have an issue deciphering between biological and technical effects
#potentially a batch effect?


#boxplot(eset,main="Boxplots of Log-Expression for KIRC & KIRP Subjects")

#### NORMALIZATION ####


#### Attempt to figure out genes which have the most variance ####
#Take genes with top 500 variances (do we want to do this in just the training set?)
KIRC.genes=names(sort(rowVars(exprs(eset[eset$type=="KIRC"])),decreasing=T))[1:500]
KIRP.genes=names(sort(rowVars(exprs(eset[eset$type=="KIRP"])),decreasing=T))[1:500]

#### KNN Classification ####
require("class")
#use function knn

#### Support Vector Machines Classification ####
require("e1071")
#use functions svm, predict.svm

#### Random Forest Classification ####
require("randomForest")
#use function randomForest