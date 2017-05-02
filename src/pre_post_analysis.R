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

rbind.match.columns <- function(input1, input2) {
  n.input1 <- ncol(input1)
  n.input2 <- ncol(input2)

  if (n.input2 < n.input1) {
    TF.names <- which(names(input2) %in% names(input1))
    column.names <- names(input2[, TF.names])
  } else {
    TF.names <- which(names(input1) %in% names(input2))
    column.names <- names(input1[, TF.names])
  }

  return(rbind(input1[, column.names], input2[, column.names]))
}

# feature x samples
clr_norm <- function(data) {
  data = t(data)
  clr(data+1, 1)
}

fit_pc <- function(data) {
  norm_data = t(clr_norm(data))

  ## use predefined test for conditional independence on gaussian data
  indepTest <- gaussCItest

  # Number of nodes
  n = dim(norm_data)[1]
  p = dim(norm_data)[2]

  print(dim(norm_data))
  suffStat <- list(C = cor(norm_data), n = n)

  pc(suffStat, indepTest, labels = colnames(norm_data), alpha = 0.05, numCores=4)
}

# feature x samples
prevalence_filt <- function(data, thresh=.1) {
  rowMeans(data > 0) >= thresh
}

spiec_easi_analysis = function(data, tax, basename) {
  pc.fit = fit_pc(data)
  ig.pc = adj2igraph(as(pc.fit@graph, 'matrix'), vertex.attr=list(name=taxa_names(tax)))

  data <- t(data)
  dim(data)
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
  ig.mb <- adj2igraph(se.mb.amgut$refit, vertex.attr=list(name=taxa_names(tax)))
  ig.gl <- adj2igraph(se.gl.amgut$refit, vertex.attr=list(name=taxa_names(tax)))
  ig.sparcc <- adj2igraph(sparcc.graph, vertex.attr=list(name=taxa_names(tax)))


  ## set size of vertex proportional to clr-mean
  vsize <- rowMeans(clr(data, 1))+6
  am.coord <- layout.fruchterman.reingold(ig.mb)

  png(filename=sprintf("%s-gl.png", basename))
  plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
  dev.off()

  png(filename=sprintf("%s-mb.png", basename))
  plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
  dev.off()

  png(filename=sprintf("%s-sparcc.png", basename))
  plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
  dev.off()

  png(filename=sprintf("%s-pc.png", basename))
  plot(ig.pc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="PC Algo (.05)")
  dev.off()


  elist.gl <- summary(triu(cov2cor(se.gl.amgut$opt.cov)*se.gl.amgut$refit, k=1))
  elist.mb <- summary(symBeta(getOptBeta(se.mb.amgut), mode='maxabs'))
  elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)


  png(filename=sprintf("-dd-hist.png", basename))
  hist(elist.sparcc[,3], main="", xlab="edge weights")
  hist(elist.mb[,3], add=TRUE, col='forestgreen')
  hist(elist.gl[,3], add=TRUE, col='red')
  dev.off()

  dd.gl <- degree.distribution(ig.gl)
  dd.mb <- degree.distribution(ig.mb)
  dd.sparcc <- degree.distribution(ig.sparcc)
  dd.pc <- degree.distribution(ig.pc)

  png(filename=sprintf("%s-dd.png", basename))
  plot(0:(length(dd.mb)-1), dd.mb, col="red" , type='b', ylab="Frequency", xlab="Degree", main="Degree Distributions", ylim=c(0, .40))
  points(0:(length(dd.sparcc)-1), dd.sparcc, type='b')
  points(0:(length(dd.gl)-1), dd.gl, col="forestgreen", type='b')
  points(0:(length(dd.pc)-1), dd.pc, col="blue", type='b')
  legend("topright", c("MB", "glasso", "sparcc", "PC"),
         col=c("forestgreen", "red", "black", "blue"), pch=1, lty=1)
  dev.off()

  list(gl=ig.gl, mb=ig.mb, pc=ig.pc, sparcc=ig.sparcc, pc.fit=pc.fit)
}

## Collapse the Datasets by Species
mct.otu.pre.l7 = summarize_taxonomy(mct.otu.pre, level=7)
mct.otu.post.l7 = summarize_taxonomy(mct.otu.post, level=7)

dim(mct.otu.pre.l7)
dim(mct.otu.post.l7)


## Remove Low Prevalence Taxa
## Microbes Must Be Present in Greater than 10% of Samples
filter.taxa <- prevalence_filt(cbind(mct.otu.pre.l7, mct.otu.post.l7))

mct.otu.pre.filt.l7 <- mct.otu.pre.l7[filter.taxa,]
mct.otu.post.filt.l7 <- mct.otu.post.l7[filter.taxa,]

dim(mct.otu.pre.filt.l7)
dim(mct.otu.post.filt.l7)

clr_norm(mct.otu.pre.filt.l7)

## Counts per Genome
png(filename="results/plots-hist-mct-otu-pre-filt-l7.png")
otu.counts <- rowSums(mct.otu.pre.filt.l7 > 0)
hist(otu.counts, breaks=30, main="Histogram of Taxa Counts (Pre)")
dev.off()

png(filename="results/plots-hist-mct-otu-post-filt-l7.png")
otu.counts <- rowSums(mct.otu.post.filt.l7 > 0)
hist(otu.counts, breaks=30, main="Histogram of Taxa Counts (Post)")
dev.off()

## Load taxonomy Table
mct.otu.tax <- as.matrix(read.delim("results/taxatable-subset.txt", row.names =  1,  sep="\t", as.is=T))

## Create the Phylogenetics Tables
mct.otu.pre <- otu_table(mct.otu.pre.l7, taxa_are_rows = T)
mct.otu.post <- otu_table(mct.otu.post.l7, taxa_are_rows = T)
mct.otu.pre.filt <- otu_table(mct.otu.pre.filt.l7, taxa_are_rows = T)
mct.otu.post.filt <- otu_table(mct.otu.post.filt.l7, taxa_are_rows = T)

## Analysis Pipeline
mct.tax <- tax_table(mct.otu.tax)

## Phylo Objects
mct.phy.pre <- phyloseq(mct.otu.pre.l7, mct.tax)
mct.phy.post <- phyloseq(mct.otu.post.l7, mct.tax)
mct.phy.pre.filt <- phyloseq(mct.otu.post.filt.l7, mct.tax)
mct.phy.post.filt <- phyloseq(mct.otu.post.filt.l7, mct.tax)

## Create the IG2
# Calculate the number of cores
# no_cores <- detectCores() - 1

# Initiate cluster
# cl <- makeCluster(4)

jobs = list(list(mct.otu.pre.l7, "results/mct-otu-pre"), list(mct.otu.post.l7, "results/mct-otu-post"), list(mct.otu.post.filt.l7, "results/mct-otu-pre-filt"), list(mct.otu.post.filt.l7, "results/mct-otu-post-filt"))

# clusterExport(cl, list("spiec_easi_analysis", "clr_norm", "fit_pc", "mct.tax"))

mapping <- function(x) {
  spiec_easi_analysis(x[[1]], mct.tax, x[[2]])
}

ig3 <- lapply(jobs, mapping)

# ig2.pre <- spiec_easi_analysis(mct.otu.pre.l7, mct.tax, "results/mct-otu-pre")
# ig2.post <- spiec_easi_analysis(mct.otu.post.l7, mct.tax, "results/mct-otu-post")
# ig2.pre.filt <- spiec_easi_analysis(mct.otu.post.filt.l7, mct.tax, "results/mct-otu-pre-filt")
# ig2.post.filt <- spiec_easi_analysis(mct.otu.post.filt.l7, mct.tax, "results/mct-otu-post-filt")

