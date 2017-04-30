library(phyloseq)
library(igraph)
library(SpiecEasi)
library(Matrix)
library(pcalg)
library(graph)

## Load Pre and Post Datasets
mct.otu.pre <- read.delim("data/mct-v3/otutable-subset-n50000-s10-norm-no-soylent-no-transition-pre.txt", sep="\t", row = 1, as.is=T, skip = 1)
mct.otu.post <- read.delim("data/mct-v3/otutable-subset-n50000-s10-norm-no-soylent-no-transition-post.txt", sep="\t", row = 1, as.is=T, skip = 1)
mct.food <- read.delim("results/taxatable-foodtable_for_taxonomy-subset_L3.txt", sep="\t", row = 1, as.is=T, skip = 1)


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

## Collapse the Datasets
mct.otu.pre.l7 = summarize_taxonomy(mct.otu.pre, level=7)
mct.otu.post.l7 = summarize_taxonomy(mct.otu.post, level=7)
# mct.food.l3 = summarize_taxonomy(mct.food, level=3)

spiec_easi_analysis = function(data, tax) {
  pc.fit = fit_pc(data)
  ig.pc = graph_from_adjacency_matrix(as(pc.fit@graph, 'matrix'))

  se.mb.amgut <- spiec.easi(data, method='mb', lambda.min.ratio=1e-2,
                            nlambda=20, icov.select.params=list(rep.num=50))
  se.gl.amgut <- spiec.easi(data, method='glasso', lambda.min.ratio=1e-2,
                            nlambda=20, icov.select.params=list(rep.num=50))
  sparcc.amgut <- sparcc(data)

  ## Define arbitrary threshold for SparCC correlation matrix for the graph
  sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
  diag(sparcc.graph) <- 0
  sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
  ## Create igraph objects
  ig.mb <- adj2igraph(se.mb.amgut$refit)
  ig.gl <- adj2igraph(se.gl.amgut$refit)
  ig.sparcc <- adj2igraph(sparcc.graph)


  ## set size of vertex proportional to clr-mean
  vsize <- rowMeans(clr(data, 1))+6
  am.coord <- layout.fruchterman.reingold(ig.mb)

  plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
  plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
  plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
  plot(ig.pc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="PC Algo (.05)")


  elist.gl <- summary(triu(cov2cor(se.gl.amgut$opt.cov)*se.gl.amgut$refit, k=1))
  elist.mb <- summary(symBeta(getOptBeta(se.mb.amgut), mode='maxabs'))
  elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)


  hist(elist.sparcc[,3], main="", xlab="edge weights")
  hist(elist.mb[,3], add=TRUE, col='forestgreen')
  hist(elist.gl[,3], add=TRUE, col='red')

  dd.gl <- degree.distribution(ig.gl)
  dd.mb <- degree.distribution(ig.mb)
  dd.sparcc <- degree.distribution(ig.sparcc)

  plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,.35), type='b',
       ylab="Frequency", xlab="Degree", main="Degree Distributions")
  points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
  points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
  legend("topright", c("MB", "glasso", "sparcc"),
         col=c("forestgreen", "red", "black"), pch=1, lty=1)

  ig2.mb <- adj2igraph(se.mb.amgut$refit,  vertex.attr=list(name=taxa_names(mct.phy)))
  plot_network(ig2.mb, tax, type='taxa', color="Kingdom", label=NULL)
  plot_network(ig2.mb, tax, type='taxa', color="Phylum", label=NULL)
  plot_network(ig2.mb, tax, type='taxa', color="Class", label=NULL)
  plot_network(ig2.mb, tax, type='taxa', color="Order", label=NULL)
  plot_network(ig2.mb, tax, type='taxa', color="Family", label=NULL)
  plot_network(ig2.mb, tax, type='taxa', color="Genus", label=NULL)


  ig2.mb <- adj2igraph(se.gl.amgut$refit,  vertex.attr=list(name=taxa_names(mct.phy)))
  plot_network(ig2.mb, tax, type='taxa', color="Order", label=NULL)
}

mct.otu.tax <- as.matrix(read.delim("results/taxatable-subset.txt", row.names =  1,  sep="\t", as.is=T))
#mct.food.tax <- as.matrix(read.delim("results/taxatable-foodtable.txt", row.names =  1,  sep="\t", as.is=T))

mct.otu <- otu_table(mct.otu.pre.l7, taxa_are_rows = T)
mct.tax <- tax_table(mct.otu.tax)


mct.phy <- phyloseq(mct.otu, mct.tax)

spiec_easi_analysis(mct.otu.pre.l7, mct.phy)
