library(SpiecEasi)
library(graph)
library(bnlearn)
library(pcalg)
library(igraph)
library(Rgraphviz)

## Summarize taxonomy to level 5 (family)
summarize_taxonomy = function(data, level=3) {
  split <- strsplit(data$taxonomy, "; ")
  drop_rows <- sapply(split, function(x) length(x) >= level)
  split <- split[drop_rows]
  data <- data[drop_rows,]
  data$taxonomy <- sapply(split,function(x) paste(x[1:level], collapse=";"))

  sample_no <- ncol(data)-1    # No. samples = no. columns (- 1 for taxonomy col)

  data <- aggregate(data[,1:sample_no], by=list(data$taxonomy), FUN=sum)
  rownames(data) <- data[,1]    # Set the rownames to the taxonomy column
  data[,-1]
}

clr_norm <- function(data) {
  data = t(data)
  t(clr(data+1, 1))
}

fit_pc <- function(data) {
  norm_data = clr_norm(data)

  ## use predefined test for conditional independence on gaussian data
  indepTest <- gaussCItest

  # Number of nodes
  n = dim(norm_data)[1]
  p = dim(norm_data)[2]

  suffStat <- list(C = cor(norm_data), n = n)

  pc.fit <- pc(suffStat, indepTest, p = p, alpha = 0.05)
}

foodtable_subset = read.delim("data/mct-v3/foodtable_for_taxonomy-subset.txt", sep="\t", skip = 1, row = 1, as.is=T)
pre_subset = read.delim("data/mct-v3/otutable-subset-n50000-s10-norm-no-soylent-no-transition-pre.txt", sep="\t", skip = 1, row = 1, as.is=T)

pre_subset_l1 = summarize_taxonomy(pre_subset, level=1)
pre_subset_l7 = summarize_taxonomy(pre_subset, level=7)

#colnames(data)

#data = summarize_taxonomy(data)

# norm_pre_subset_l7 <- clr_norm(pre_subset_l7)

#library(bnlearn)
#boot.strength(as.data.frame(norm_pre_subset_l7), algorithm = "hc")

pc.fit = fit_pc(pre_subset_l7)
se.mb.mct <- spiec.easi(pre_subset_l7, method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, icov.select.params=list(rep.num=50))
se.gl.mct <- spiec.easi(pre_subset_l7, method='glasso', lambda.min.ratio=1e-2,
                          nlambda=20, icov.select.params=list(rep.num=50))
sparcc.mct <- sparcc(pre_subset_l7)


## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc.mct$Cor) >= 0.3
diag(sparcc.graph) <- 0
sparcc.graph <- as.matrix(sparcc.graph, sparse=TRUE)

## Create igraph objects
ig.mb <- adj2igraph(se.mb.mct$refit)
ig.gl <- adj2igraph(se.gl.mct$refit)
ig.sparcc <- adj2igraph(sparcc.graph)
ig.pc = graph_from_adjacency_matrix(as(pc.fit@graph, 'matrix'))

## set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(pre_subset_l7, 1))+6
am.coord <- layout.fruchterman.reingold(ig.pc)

plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
plot(ig.pc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="PC Algo (.05)")
