source("src/normalization.R")

library(bnlearn)
library(parallel)

mapping.pre <- read.delim("data/mct-v3/map-subset-no-soylent-no-transition-pre.txt", sep="\t", as.is=T)
mapping.post <- read.delim("data/mct-v3/map-subset-no-soylent-no-transition-post.txt", sep="\t", as.is=T)

mapping.post <- mapping.post[mapping.post$Supplement == 1, ]

mapping <- rbind.data.frame(mapping.pre, mapping.post)

user.names <- unique(mapping$UserName)

sample_names = list()
for (user.name in user.names) {
  subset <- mapping[mapping$UserName == user.name, ]
  days <- as.numeric(lapply(strsplit(subset$StudyDayNo, "[.]"), function(x) x[2]))
  for (day in days) {
    if ((day + 1) %in% days) {
      from <- subset[days == day, ]
      to <- subset[days == (day + 1), ]
      colnames(from) <- as.character(sapply(colnames(from), function(x) sprintf("%s.from", x)))
      colnames(to) <- as.character(sapply(colnames(to), function(x) sprintf("%s.to", x)))
      sample_names[[length(sample_names) + 1]] <- cbind(from, to)
    }
  }
}

mapping.subset <- do.call(rbind.data.frame, sample_names)
write.table(mapping.subset, file="results/mapping-subset-no-soylent-long.txt")

# Load the Food and Phylogenetic Tables
mct.food <- read.delim("data/mct-v3/foodtable_for_taxonomy-subset.txt", sep="\t", row = 1, as.is=T, skip = 1)
mct.otu <- read.delim("data/mct-v3/otutable-subset-n50000-s10-norm-no-soylent-no-transition.txt", sep="\t", row = 1, as.is=T, skip = 1)

summarize_taxonomy = function(data, level=7) {
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

mct.food.l2 <- summarize_taxonomy(mct.food, level = 2)
mct.otu.l7 <- summarize_taxonomy(mct.otu, level = 6)

mapping.subset[c("X.SampleID.from", "X.SampleID.to")]

mct.food.from <- t(mct.food.l2[mapping.subset$X.SampleID.from])
colnames(mct.food.from) <- as.character(sapply(colnames(mct.food.from), function(x) sprintf("%s.from", x)))

mct.otu.from <- t(mct.otu.l7[mapping.subset$X.SampleID.from])
colnames(mct.otu.from) <- as.character(sapply(colnames(mct.otu.from), function(x) sprintf("%s.from", x)))

mct.otu.to <- t(mct.otu.l7[mapping.subset$X.SampleID.to])
colnames(mct.otu.to) <- as.character(sapply(colnames(mct.otu.to), function(x) sprintf("%s.to", x)))

# samples x features
prevalence_filt <- function(data, thresh=.1) {
  colMeans(data > 0) >= thresh
}

mct.food.from <- mct.food.from[ ,prevalence_filt(mct.food.from)]
mct.otu.from <- mct.otu.from[ ,prevalence_filt(mct.otu.from)]
mct.otu.to <- mct.otu.to[ ,prevalence_filt(mct.otu.to)]

mct.food.otu.otu <- cbind(mct.food.from, mct.otu.from, mct.otu.to)
mct.otu.otu <- cbind(mct.otu.from, mct.otu.to)
mct.food.otu <- cbind(mct.food.from, mct.otu.to)


# samples x features
clr_norm <- function(data) {
  clr(data+1, 1)
}

mct.food.otu.otu <- clr_norm(mct.food.otu.otu)
mct.food.otu <- clr_norm(mct.food.otu)
mct.otu.otu <- clr_norm(mct.otu.otu)

# From Cannot Depend on From
otu.from.otu.to <- tiers2blacklist(list(colnames(mct.otu.from), colnames(mct.otu.to)))
food.from.otu.to <- tiers2blacklist(list(colnames(mct.food.from), colnames(mct.otu.to)))

n <- dim(mct.food.from)[2]
rows <- 2*((n * n) - n)
res <- matrix("", rows, 2)
pos <- 1
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      x <- colnames(mct.food.from)[i]
      y <- colnames(mct.food.from)[j]
      res[pos, 1] <- x
      res[pos + 1, 1] <- y
      res[pos, 2] <- y
      res[pos + 1, 2] <- x
      pos <- pos + 2
    }
  }
}

food.from.food.from <- as.data.frame(res)
colnames(food.from.food.from) <- c("from", "to")

n <- dim(mct.otu.from)[2]
rows <- 2*((n * n) - n)
res <- matrix("", rows, 2)
pos <- 1
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      x <- colnames(mct.otu.from)[i]
      y <- colnames(mct.otu.from)[j]
      res[pos, 1] <- x
      res[pos + 1, 1] <- y
      res[pos, 2] <- y
      res[pos + 1, 2] <- x
      pos <- pos + 2
    }
  }
}

otu.from.otu.from <- as.data.frame(res)
colnames(otu.from.otu.from) <- c("from", "to")


n <- dim(mct.otu.from)[2]
m <- dim(mct.food.from)[2]
rows <- 2*(m * n)
res <- matrix("", rows, 2)
pos <- 1
for (i in 1:m) {
  for (j in 1:n) {
    x <- colnames(mct.food.from)[i]
    y <- colnames(mct.otu.from)[j]
    res[pos, 1] <- x
    res[pos + 1, 1] <- y
    res[pos, 2] <- y
    res[pos + 1, 2] <- x
    pos <- pos + 2
  }
}

food.from.otu.from <- as.data.frame(res)
colnames(food.from.otu.from) <- c("from", "to")


# Create the blacklists
blacklist.food.otu.otu <- rbind(food.from.otu.to, food.from.food.from, otu.from.otu.to, otu.from.otu.from, food.from.otu.from)
dim(blacklist.food.otu.otu)
blacklist.food.otu <- rbind(food.from.otu.to, food.from.food.from)
dim(blacklist.food.otu)
blacklist.otu.otu <- rbind(otu.from.otu.to, otu.from.otu.from)
dim(blacklist.otu.otu)


# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

traits <- names(tail(sort(colSums(mct.otu.to)), 5))

# norm.mapping.otu$Treatment <- as.factor(norm.mapping.otu$Treatment)
food.val = bn.cv(as.data.frame(t(mct.food.otu)), cluster=cl, bn = "hc", loss = "mse", algorithm.args = list(blacklist = blacklist.food.otu), loss.args = list(target = traits[1]))
save(food.val, file="results/mct.food.l6.val.RData")

otu.val = bn.cv(as.data.frame(t(mct.otu.otu)), cluster=cl, bn = "hc", loss = "mse", algorithm.args = list(blacklist = blacklist.otu.otu), loss.args = list(target = traits[1]))
save(otu.val, file="results/otu.l6.val.RData")

# mct.food.otu.otu -> as.data.frame(t(mct.food.otu.otu))
# for (name in blacklist.food.otu.otu$from) {
#   if (!name %in% rownames(mct.food.otu.otu)) {
#     print(name)
#   }
# }

food.otu.val = bn.cv(as.data.frame(t(mct.food.otu.otu)), cluster=cl, bn = "hc", loss = "mse", algorithm.args = list(blacklist = blacklist.food.otu.otu), loss.args = list(target = traits[1]))
save(food.otu.val, file="results/food.otu.l6.val.RData")

stopCluster(cl)

