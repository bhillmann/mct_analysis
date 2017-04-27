options(repos = c(CRAN = "https://cran.revolutionanalytics.com"))

install.packages(c("bnlearn", "igraph"))

source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "Rgraphviz", "RBGL"))
install.packages("gRain")

library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)

library(graph)
library(bnlearn)
library(pcalg)
library(igraph)
library(Rgraphviz)

load('data/amgut1.filt.rda')
data = amgut1.filt
norm_data = t(clr(data+1, 1))

# Number of nodes
n = dim(norm_data)[1]
p = dim(norm_data)[2]
# Inclusion prob

## use predefined test for conditional independence on gaussian data
indepTest <- gaussCItest

## the functin gaussCItest needs as input the correlation matrix C and
## the sample size n
suffStat <- list(C = cor(norm_data), n = n)
## estimate the causal structure
pc.fit <- pc(suffStat, indepTest, p = p, alpha = 0.001)

se.mb.amgut <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, icov.select.params=list(rep.num=50))
se.gl.amgut <- spiec.easi(amgut1.filt, method='glasso', lambda.min.ratio=1e-2,
                          nlambda=20, icov.select.params=list(rep.num=50))
sparcc.amgut <- sparcc(amgut1.filt)

## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
diag(sparcc.graph) <- 0
sparcc.graph <- as.matrix(sparcc.graph, sparse=TRUE)

## Create igraph objects
ig.mb <- adj2igraph(se.mb.amgut$refit)
ig.gl <- adj2igraph(se.gl.amgut$refit)
ig.sparcc <- adj2igraph(sparcc.graph)

ig.pc = graph_from_adjacency_matrix(as(pc.fit@graph, 'matrix'))

## set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(amgut1.filt, 1))+6
am.coord <- layout.fruchterman.reingold(ig.pc)

plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
plot(ig.pc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="PC Algo (.001)")

dd.gl <- degree.distribution(ig.pc)


plot(0:(length(dd.gl)-1), dd.gl, ylim=c(0,.35), type='b',
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

