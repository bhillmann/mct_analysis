library(phyloseq)
library(igraph)
library(SpiecEasi)
library(Matrix)
library(pcalg)
library(graph)

## Load Pre and Post Datasets
mct.otu.pre <- read.delim("data/mct-v3/otutable-subset-n50000-s10-norm-no-soylent-no-transition-pre.txt", sep="\t", row = 1, as.is=T, skip = 1)
mct.otu.post <- read.delim("data/mct-v3/otutable-subset-n50000-s10-norm-no-soylent-no-transition-post.txt", sep="\t", row = 1, as.is=T, skip = 1)

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

# feature x samples
prevalence_filt <- function(data, thresh=.1) {
  rowMeans(data > 0) >= thresh
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

norm.otu <- clr_norm(cbind(mct.otu.pre.filt.l7, mct.otu.post.filt.l7))
time_col <- as.matrix(rep("post", dim(norm.otu)[1]))
time_col[1:dim(mct.otu.pre.filt.l7)[2]] <- "pre"
colnames(time_col) <- "time"
norm.otu <- cbind(norm.otu, time_col)

mapping <- read.delim("data/mct-v3/map-subset.txt", sep="\t", row = 1, as.is=T)
norm.otu[match(rownames(mapping), rownames(norm.otu)),]

dim(merge(mapping, norm.otu, by=intersect(rownames(mapping), rownames(norm.otu))))


write.table(norm.otu, file = "results/mct.otu.filt.l7.clr.txt", sep="\t")

# write.table(, file = "results/mct.otu.post.filt.l7.clr.txt", sep="\t")
