args <- commandArgs(trailingOnly = TRUE)

n <- as.numeric(args[1])/3
offset <- 0
for(i in 1:n)
{
	x <- unlist(strsplit(args[i*2+1+offset], split=","))
      	x <- as.numeric(x)
      	print(args[i*2 + 2 + offset])
     	paril <- 1 - pnorm(as.numeric(args[i*2 + offset]), mean(x), sd(x))
	poutlier <- pbinom(0, length(x), paril)
     	print(1 - poutlier)
      	offset <- offset + 1
}
