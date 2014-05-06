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
a=rowMeans(exprs(eset)[,eset$type=="KIRC" & eset$train])>=10
b=rowMeans(exprs(eset)[,eset$type=="KIRP" & eset$train])>=10
f=(a|b)
table(f)
new_eset=eset[f]
#17,136 genes after filter
par(mfrow=c(1,2))
hist(log(exprs(eset)+1),main="Pre-Filter Log Exprs.")
hist(log(exprs(new_eset)+1),main="Filtered Log Exprs")

boxplot(new_eset,main="Boxplots of Log-Expression for KIRC & KIRP Subjects")
#look at Mean-Difference Plots from subjects within cell-types and between cell-types
#might need to untransform eset from log since I think MDPlot takes the log for you
par(mfrow=c(2,2))
MDPlot(new_eset,c(1,5),main="Within KIRC")
MDPlot(new_eset,c(48,60),main="Within KIRP")
MDPlot(new_eset,c(12,50),main="KIRC vs. KIRP")
MDPlot(new_eset,c(30,70),main="KIRC vs. KIRP")
#looks like we may have an issue deciphering between biological and technical effects
#potentially a batch effect?


#boxplot(eset,main="Boxplots of Log-Expression for KIRC & KIRP Subjects")
#### Specify training and test sets to use ####

train.eset=new_eset[,eset$train]
test.eset=new_eset[,!eset$train]

#### NORMALIZATION ####
#Normalize Training Set
train.norm_eset=betweenLaneNormalization(train.eset,which="upper",offset=F,round=F)
type=eset$type[eset$train]
train.lognorm_exprs=log(exprs(train.norm_eset)+1)
upper=fivenum(exprs(train.norm_eset))[4]

#Normalize the Test Set based on the upper quantile of the train set
#find upper quantile of each gene
test.norm_exprs=matrix(data=NA,nrow=dim(exprs(new_eset))[1],ncol=dim(test.eset)[2])
up.test=vector()
for(i in 1:dim(test.eset)[2]){
  up.test[i]=quantile(exprs(test.eset)[,i],.75)
  test.norm_exprs[,i]=upper*exprs(test.eset)[,i]/up.test[i] 
}
par(mfrow=c(1,1))
boxplot(log(test.norm_exprs+1))

#check normalization 
par(mfrow=c(2,2))
MDPlot(train.norm_eset,c(1,5),main="Within KIRC")
MDPlot(train.norm_eset,c(27,35),main="Within KIRP")
MDPlot(train.norm_eset,c(12,40),main="KIRC vs. KIRP")
MDPlot(train.norm_eset,c(3,29),main="KIRC vs. KIRP")

par(mfrow=c(1,1))
hist(log(exprs(train.norm_eset)+1))
#### Attempt to figure out genes which have the most variance ####
#Take genes with top 500 variances (do we want to do this in just the training set?)
kirc.genes=names(sort(rowVars(train.lognorm_exprs[,type=="KIRC"]),decreasing=T))[1:500]
kirp.genes=names(sort(rowVars(train.lognorm_exprs[,type=="KIRP"]),decreasing=T))[1:500]
top.genes=names(sort(rowVars(train.lognorm_exprs),decreasing=T))[1:500]


kirc.var=(rowVars(train.lognorm_exprs[,type=="KIRC"]))
kirp.var=(rowVars(train.lognorm_exprs[,type=="KIRP"]))
kirc.mean=(rowMeans(train.lognorm_exprs[,type=="KIRC"]))
kirp.mean=(rowMeans(train.lognorm_exprs[,type=="KIRP"]))
kirc.n=sum(type=="KIRC")
kirp.n=sum(type=="KIRP")
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
sort.pvals=sort(pvals)

ngenes=40
top.pvals=sort.pvals[1:ngenes]

pval.ind=vector()
for(i in 1:ngenes){
  pval.ind[i]=which(pvals==top.pvals[i])
}

filt.train=exprs(train.norm_eset)[pval.ind,]
filt.test=test.norm_exprs[pval.ind,]
rownames(filt.test)=rownames(filt.train)
fit=princomp(t(filt.train),cor=T)

loadings(fit)
par(mfrow=c(1,1))
plot(fit)
score1=fit$scores[,1]
score2=fit$scores[,2]
plot(score1,score2,col=ifelse((train.norm_eset$type=="KIRC"),"red","blue"),main="PC1 vs. PC2, Colored by Cancer")
legend("topright",c("KIRC","KIRP"),text.col=c("red","blue"),cex=.8)

fit$loadings[,1]


#### Classifications ####
trainy=eset$type[which(eset$train)]
testy=eset$type[which(!eset$train)]
trainx=as.data.frame(t(filt.train))
testx=as.data.frame(t(filt.test))
#### LDA Classification ####
require("MASS")
#use function lda
LDA=lda(trainy~.,data=trainx)
lda.pred=predict(LDA,testx)
table(pred=lda.pred$class,true=testy)

#### Support Vector Machines Classification ####
require("e1071")
#use functions svm, predict.svm

SVM=svm(trainy~.,data=trainx)
svm.pred=predict(SVM,testx)
table(pred=svm.pred,true=testy)

#### Random Forest Classification ####
require("randomForest")
#I tried throwing all the genes into random forest but it was taking a very long time
#trainx=t(lognorm_eset[,which(eset$train)])
#testx=t(lognorm_eset[,which(!eset$train)])

#use function randomForest
ngenes.rf=17136
top.pvals.rf=sort.pvals[1:ngenes.rf]

pval.ind.rf=vector()
for(i in 1:ngenes.rf){
  pval.ind.rf[i]=which(pvals==top.pvals.rf[i])
}

filt.train.rf=exprs(train.norm_eset)[pval.ind.rf,]
filt.test.rf=test.norm_exprs[pval.ind.rf,]
trainx.rf=as.data.frame(t(filt.train.rf))
testx.rf=as.data.frame(t(filt.test.rf))
colnames(testx.rf)=colnames(trainx.rf)
rf.fit=randomForest(x=trainx.rf,y=trainy,proximity=T,importance=T)
rf.pred=predict(rf.fit,testx.rf)
table(pred=rf.pred,true=testy)

#Interestingly random forest does worse when I feed it more genes...

