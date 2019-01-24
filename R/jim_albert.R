# jim_albert.R
# https://www.wiley.com/en-gb/Causal+Inference+in+Statistics%3A+A+Primer-p-9781119186847
# http://www.statisticsviews.com/details/video/7941841/Teaching-Causality---Part-II-with-Judea-Pearl-and-Rob-Gould.html


library(LearnBayes)
library(eigenmodel)

P <- matrix(c( 0.5, 0.5, 0, 0, 0, 0, 
              0.25, 0.5, 0.25, 0, 0, 0,
              0, 0.25, 0.5, 0.25, 0, 0,
              0, 0, 0.25, 0.5, 0.25, 0,
              0, 0, 0, 0.25, 0.5, 0.25,
              0, 0, 0, 0, 0.5, 0.5), nrow=6,byrow=TRUE)

s <- array(0,c(50000, 1)) 

s[1] <- 3 # starting location is 3

for(j in 2:50000){
  s[j] <- sample(1:6, size=1, prob = P[s[j-1],])}

m <- c(500, 2000, 8000, 50000)

for(i in 1:4){
  print(table(s[1:m[i]])/m[i])}


