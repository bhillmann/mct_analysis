options(repos = c(CRAN = "https://cran.revolutionanalytics.com"))

install.packages(c("bnlearn", "igraph"))

source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "Rgraphviz", "RBGL", "phyloseq"))
install.packages("gRain")
install.packages("devtools")
install.packages("pcalg")

library(devtools)
install_github("zdk123/SpiecEasi")
