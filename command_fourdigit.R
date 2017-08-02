args <- commandArgs(trailingOnly = TRUE)

x<-unlist(strsplit(args[2],split=","))
x<-as.numeric(x)
#args[3]
paril<-1-pnorm(as.numeric(args[1]),mean(x),sd(x))
poutlier<-pbinom(0,length(x),paril)
1-poutlier

