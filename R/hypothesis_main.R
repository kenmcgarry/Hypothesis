# hypothesis_main.R
# Hypothesis (theory) revision using Bayesian techniques.
# # https://github.com/vonjd/OneR
# https://www.rdocumentation.org/packages/LaplacesDemon/versions/16.1.0
# https://cran.r-project.org/web/packages/LaplacesDemon/vignettes/LaplacesDemonTutorial.pdf
# started 19/5/18

memory.limit(1010241024*1024) # allocate RAM memory (10 GBs)
setwd("C:/R-files/hypothesis")  # now point to where the new code lives
source("hypothesis_functions.R")  # load in the functions required for this work. 

# Stuff below is simply to provide a webbased front-end 9eventually)
#rsconnect::setAccountInfo(name='rulebase',token='DD8F9D725C3798523EAD8601AC33C0AD',secret='2uXoDpo+Yg/1CKYnzacibYdFpLfP02NEKbWTBLme')
#rsconnect::deployApp('C:/R-files/hypothesis/frontend.Rmd')
# https://www.shinyapps.io/admin/#/dashboard
# https://rmarkdown.rstudio.com/flexdashboard/shiny.html#getting_started
# https://rulebase.shinyapps.io/frontend/