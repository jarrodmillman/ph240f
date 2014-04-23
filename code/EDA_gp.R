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
a=rowMeans(exprs(eset)[,eset$type=="KIRC"])>=10
b=rowMeans(exprs(eset)[,eset$type=="KIRP"])>=10
f=(a|b)
table(f)
new_eset=eset[f]
hist(log(exprs(eset)+1))
hist(log(exprs(new_eset)+1))

boxplot(new_eset,main="Boxplots of Log-Expression for KIRC & KIRP Subjects")
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
norm_eset=betweenLaneNormalization(new_eset,which="upper",offset=F,round=F)
lognorm_eset=log(exprs(norm_eset)+1)

#check normalization 
par(mfrow=c(2,2))
MDPlot(norm_eset,c(1,5),main="within")
MDPlot(norm_eset,c(48,60),main="within")
MDPlot(norm_eset,c(12,50),main="bw")
MDPlot(norm_eset,c(12,70),main="bw")

hist(log(exprs(norm_eset)+1))
#### Attempt to figure out genes which have the most variance ####
#Take genes with top 500 variances (do we want to do this in just the training set?)
KIRC.genes=names(sort(rowVars(exprs(lognorm_eset[norm_eset$type=="KIRC"])),decreasing=T))[1:500]
KIRP.genes=names(sort(rowVars(exprs(lognorm_eset[norm_eset$type=="KIRP"])),decreasing=T))[1:500]
top.genes=names(sort(rowVars(exprs(lognorm_eset)),decreasing=T))[1:500]


KIRC.var=(rowVars(lognorm_eset[,norm_eset$type=="KIRC"]))
KIRP.var=(rowVars(lognorm_eset[,norm_eset$type=="KIRP"]))
KIRC.mean=(rowMeans(lognorm_eset[,norm_eset$type=="KIRC"]))
KIRP.mean=(rowMeans(lognorm_eset[,norm_eset$type=="KIRP"]))
KIRC.n=sum(norm_eset$type=="KIRC")
KIRP.n=sum(norm_eset$type=="KIRP")
#### t-tests for difference of mean expression ####
#Using wiki t-test for unequal sample sizes/variances
#work in progress
se=sqrt(((KIRC.var^2)/KIRC.n)+((KIRP.var^2)/KIRP.n))
t=(KIRC.mean-KIRP.mean)/se

denomt1=(((KIRC.var^2)/(KIRC.n))^2)/(KIRC.n-1)
denomt2=(((KIRP.var^2)/(KIRP.n))^2)/(KIRP.n-1)
df=(se^4)/(denomt1+denomt2)
pvals=2*pt(-abs(t),df)

#### KNN Classification ####
require("class")
#use function knn

#### Support Vector Machines Classification ####
require("e1071")
#use functions svm, predict.svm

#### Random Forest Classification ####
require("randomForest")
#use function randomForest