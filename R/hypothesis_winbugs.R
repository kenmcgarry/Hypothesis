# hypothesis_winbugs.R

library(R2WinBUGS)

setwd("C:/common_laptop/R-files/hypothesis")  

k <- 5
n <- 10

data1 <- list("k","n")
myinits <- list(list(theta = 0.1),
                list(theta = 0.9))

bugsdir <- "C:\\WinBUGS14"
parameters <- c("theta")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
mysamples <- bugs(data1, inits=myinits, parameters,
                    model.file ="Rate_1.txt",
                    n.chains=2, n.iter=20000, n.burnin=1, n.thin=1,
                    DIC=TRUE, bugs.directory=bugsdir,
                    codaPkg=FALSE, debug=FALSE)

# The commands below are useful for a quick overview:
print(mysamples)  # a rough summary
plot(mysamples)   # a visual representation
names(mysamples)  # summarizes the variables
mysamples$summary # more detailed summary
mysamples$sims.array[1:15,,2]# array: sample, chain, parameter 

# Collect posterior samples across all chains:
theta <- mysamples$sims.list$theta 

# Now let's plot a histogram for theta. NB. Some the plots will not look good in RStudio.
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), xlab="Rate", ylab="Posterior Density") 



