require("Biobase")
require("EDASeq")
a.all <- read.csv('kirc.csv', row.names=1)
b.all <- read.csv('kirp.csv', row.names=1)
set.seed(1)
a.in <- as.logical(rbinom(n=ncol(a.all), size=1, prob=.5))
b.in <- as.logical(rbinom(n=ncol(b.all), size=1, prob=.5))
a <- a.all[,a.in]
b <- b.all[,b.in]
c.all <- cbind(a,b)
c <- c.all[rowSums(c.all)>100,]

pData <- data.frame(
    type=factor(colnames(c)==colnames(a), levels=c(TRUE,FALSE), labels=c("KIRC","KIRP")),
    train=as.logical(c(rbinom(length(a), 1, 2/3), rbinom(length(b), 1, 2/3))),
    row.names=names(c)
)

phenoData <- new("AnnotatedDataFrame", data = pData)
#eset <- ExpressionSet(assayData=as.matrix(log(c+1)), phenoData=phenoData)

#put data in a SeqExpressionSet
eset=newSeqExpressionSet(exprs=as.matrix(log(c+1)), phenoData=phenoData)

save(file='eset.rda', eset)
#boxplot(log(a+1), xlab='Subjects (n=73)', ylab="log(E(gene counts)+1)")
