require("Biobase")

load("../data/eset.rda")
summary(matrix(exprs(eset)))
sum(matrix(exprs(eset))==0)

boxplot(exprs(eset))

