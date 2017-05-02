load("results/ig3.RData")

library(igraph)

results <- list()
for (i in 3:4) {
  ig3[[i]]$pc <- ig3[[i]]$pc
  inter <- list(mb <- as_adj(ig3[[i]]$mb, sparse=F),
  gl <- as_adj(ig3[[i]]$gl, sparse=F),
  pc <- as_adj(ig3[[i]]$pc, sparse=F),
  sparcc <- as_adj(ig3[[i]]$sparcc, sparse=F))
  v <- lapply(inter, function(x) {x[upper.tri(x)]})
  results[[length(results) + 1]] <- v
}

x <- matrix(unlist(results), ncol = 8)

custom.dist <- function(my.list, my.function) {
  n <- ncol(my.list)
  mat <- matrix(0, ncol = n, nrow = n)
  for(i in 1:n) {
    for(j in 1:n) {
      mat[j, i] <- my.function(my.list[,i], my.list[,j])
    }}
  return(mat)
}

r <- custom.dist(x, function(a, b) sum(a != b))
