####Group_Project####
setwd("/Users/lktroszak/src/ph240f/data")
load("eset.rda")
require("Biobase")
require("EDASeq")
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

par(mfrow=c(1,1))
hist(log(exprs(norm_eset)+1))
#### Attempt to figure out genes which have the most variance ####
#Take genes with top 500 variances (do we want to do this in just the training set?)
kirc.genes=names(sort(rowVars(lognorm_eset[,norm_eset$type=="KIRC"]),decreasing=T))[1:500]
kirp.genes=names(sort(rowVars(lognorm_eset[,norm_eset$type=="KIRP"]),decreasing=T))[1:500]
top.genes=names(sort(rowVars(lognorm_eset),decreasing=T))[1:500]


kirc.var=(rowVars(lognorm_eset[,norm_eset$type=="KIRC"]))
kirp.var=(rowVars(lognorm_eset[,norm_eset$type=="KIRP"]))
kirc.mean=(rowMeans(lognorm_eset[,norm_eset$type=="KIRC"]))
kirp.mean=(rowMeans(lognorm_eset[,norm_eset$type=="KIRP"]))
kirc.n=sum(norm_eset$type=="KIRC")
kirp.n=sum(norm_eset$type=="KIRP")
#### t-tests for difference of mean expression ####
#Using wiki t-test for unequal sample sizes/variances
#work in progress
se=sqrt(((kirc.var^2)/kirc.n)+((kirp.var^2)/kirp.n))
t=(kirc.mean-kirp.mean)/se

denomt1=(((kirc.var^2)/(kirc.n))^2)/(kirc.n-1)
denomt2=(((kirp.var^2)/(kirp.n))^2)/(kirp.n-1)
df=(se^4)/(denomt1+denomt2)
pvals=2*pt(-abs(t),df)
par(mfrow=c(1,2))
hist(pvals,main="Gene P-Values")
hist(log(pvals),main="Gene Log P-Values")

##### P-values using FDR ####
fdr.pvals=p.adjust(pvals,method="fdr")
par(mfrow=c(1,2))
hist(log(pvals),main="Log Regular P")
hist(log(fdr.pvals),main="Log FDR P")


top.pvals=(log(pvals)<=-80)
pval.ind=which(log(pvals)<=-80)
filt.lognorm_eset=lognorm_eset[top.pvals,]

fit=princomp(t(filt.lognorm_eset),cor=T)

loadings(fit)
par(mfrow=c(1,1))
plot(fit)
score1=fit$scores[,1]
score2=fit$scores[,2]
plot(score1,score2,col=ifelse((norm_eset$type=="KIRC"),"red","blue"),main="PC1 vs. PC2, Colored by Cancer")
legend("bottomleft",c("KIRC","KIRP"),text.col=c("red","blue"),cex=.8)

fit$loadings[,1]

#### KNN Classification ####
require("class")
#use function knn

#### Support Vector Machines Classification ####
require("e1071")
#use functions svm, predict.svm
SVM=svm(norm_eset$type~fit$scores,cross=10)

#### Random Forest Classification ####
require("randomForest")
#use function randomForest