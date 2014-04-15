a = read.csv('data.csv', row.names=1)
boxplot(log(a+1))
