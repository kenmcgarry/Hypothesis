# hypothesis_main.R
# Hypothesis (theory) revision using Bayesian techniques.
# https://www.semanticscholar.org/paper/Necessity%2C-possibility-and-belief%3A-a-study-of-Evans-Handley/e5bd73cbab1477d95719fc3fc8209759bf3fc1eb?tab=abstract&citingPapersSort=is-influential&citingPapersLimit=10&citingPapersOffset=10&citedPapersSort=is-influential&citedPapersLimit=10&citedPapersOffset=0
# http://www.bnlearn.com/documentation/man/
# https://cran.r-project.org/web/packages/OneR/vignettes/OneR.html
# https://cran.r-project.org/web/packages/LaplacesDemon/vignettes/LaplacesDemonTutorial.pdf
# started 19/5/18

#memory.limit(1010241024*1024) # allocate RAM memory (10 GBs)
setwd("C:/common_laptop/R-files/hypothesis")  # now point to where the new code lives
source("hypothesis_functions.R")  # load in the functions required for this work. 

# using bnlearn model2network() function
# NB: all incoming nodes must be done in one command not peicemeal. 
bn1 <- model2network("[A][C][B|A][D|C][F|A:B:C][E|F:D]")
plot(bn1)

bn2 <- model2network("[C][F][A|C][E|A:C][D|C:E:F][B|A:C:D:F]")
plot(bn2)

bn3 <- model2network("[C][F][B|A][D|A:C:E][E|B:F][A|F]",debug = FALSE)
plot(bn3)

#
# load the data.
data(asia)
# create and plot the network structure.
dag = model2network("[Asia][smoking][tuberculosis|Asia][lungCancer|smoking][bronchitis|smoking][dyspnoea|bronchitis:either][either|tuberculosis:lungCancer][Xray|either]")
graphviz.plot(dag)


