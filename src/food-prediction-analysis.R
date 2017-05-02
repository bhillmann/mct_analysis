load("results/mct.food.l5.val.RData")
load("results/otu.l5.val.RData")
load("results/food.otu.l5.val.RData")
load("results/food.otu.val.RData")

library("bnlearn")

plot_results <- function(xval) {
  arclist = list()

  for (i in seq_along(xval)) {
    run = xval[[i]]$fitted
    print(arcs(run))
    arclist[[length(arclist) + 1]] = arcs(run)
  }

  nodes = unique(unlist(arclist))
  strength = custom.strength(arclist, nodes = nodes)
  averaged = averaged.network(strength)

  relevant.nodes = nodes(averaged)[sapply(nodes, degree, object = averaged) > 0]
  relevant.nodes = unique(unlist(relevant.nodes))

  averaged2 = subgraph(averaged, relevant.nodes)
  strength2 = strength[(strength$from %in% relevant.nodes) & (strength$to %in% relevant.nodes), ]
  # png(filename="results/plots-prediction.png")
  gR = strength.plot(averaged2, strength2, shape = "rectangle", layout = "fdp")
  # dev.off()
}

plot_results(food.otu.val)

