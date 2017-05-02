library(phyloseq)
library(igraph)
library(Matrix)
library(pcalg)
library(graph)
library(plyr)
library(bnlearn)

source("src/normalization.R")

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

mapping <- read.delim("data/mct-v3/map-subset.txt", sep="\t", row = 1, as.is=T)

mapping <- mapping[c("Supplement", "Treatment")]

norm.mapping.otu <- cbind(norm.otu, mapping[match(rownames(norm.otu), rownames(mapping)),])

write.table(norm.mapping.otu, file = "results/mct.otu.filt.l7.clr.txt", sep="\t")

library(bnlearn)
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

norm.mapping.otu$Treatment <- as.factor(norm.mapping.otu$Treatment)
norm.mapping.otu$Supplement <- as.factor(norm.mapping.otu$Supplement)

s1 <- norm.mapping.otu[norm.mapping.otu$Supplement == 1,]
s1 <- s1[, !names(norm.mapping.otu) %in% c("Supplement")]
s2 <- norm.mapping.otu[norm.mapping.otu$Supplement == 2,]
s2 <- s2[, !names(norm.mapping.otu) %in% c("Supplement")]

pre <- norm.mapping.otu[norm.mapping.otu$Treatment == "Pre",]
pre <- pre[, !names(norm.mapping.otu) %in% c("Treatment")]
post <- norm.mapping.otu[norm.mapping.otu$Treatment == "Post",]
post <- post[, !names(norm.mapping.otu) %in% c("Treatment")]


blacklist_supplement <- tiers2blacklist(list("Supplement", colnames(norm.mapping.otu[, !names(norm.mapping.otu) %in% c("Supplement", "Treatment")])))
blacklist_treatment <- tiers2blacklist(list("Treatment", colnames(norm.mapping.otu[, !names(norm.mapping.otu) %in% c("Supplement", "Treatment")])))

norm.mapping.otu$Treatment <- as.factor(norm.mapping.otu$Treatment)
s1.val = bn.cv(s1, cluster=cl, bn = "si.hiton.pc", loss = "pred-lw-cg", algorithm.args = list(blacklist = blacklist_treatment), loss.args = list(target = "Treatment"))
save(s1.val, file="results/s1.l5.val.RData")

s2.val = bn.cv(s2, cluster=cl, bn = "si.hiton.pc", loss = "pred-lw-cg", algorithm.args = list(blacklist = blacklist_treatment), loss.args = list(target = "Treatment"))
save(s2.val, file="results/s2.l5.val.RData")

pre.val = bn.cv(pre, cluster=cl, bn = "si.hiton.pc", loss = "pred-lw-cg", algorithm.args = list(blacklist = blacklist_supplement), loss.args = list(target = "Supplement"))
save(pre.val, file="results/pre.l5.val.RData")

post.val = bn.cv(post, cluster=cl, bn = "si.hiton.pc", loss = "pred-lw-cg", algorithm.args = list(blacklist = blacklist_supplement), loss.args = list(target = "Supplement"))
save(post.val, file="results/post.l5.val.RData")

stopCluster(cl)


arclist = list()
relevant.nodes = list()

for (i in seq_along(xval)) {
  run = xval[[i]]$fitted

  arclist[[length(arclist) + 1]] = arcs(run)
  relevant.nodes[[length(relevant.nodes) + 1]] = list(children(run, "Supplement"), "Supplement")
}

nodes = unique(unlist(arclist))
strength = custom.strength(arclist, nodes = nodes)
averaged = averaged.network(strength)

# relevant.nodes = nodes(averaged)[sapply(nodes, degree, object = averaged) > 0]
relevant.nodes = unique(unlist(relevant.nodes))

averaged2 = subgraph(averaged, relevant.nodes)
strength2 = strength[(strength$from %in% relevant.nodes) & (strength$to %in% relevant.nodes), ]
png(filename="results/plots-prediction.png")
gR = strength.plot(averaged2, strength2, shape = "rectangle", layout = "fdp")
dev.off()
