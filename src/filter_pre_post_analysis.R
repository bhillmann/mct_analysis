load("results/s1.val.RData")
load("results/s2.val.RData")

s1.val
s2.val

plot_results <- function(xval) {
  arclist = list()

  for (i in seq_along(xval)) {
    run = xval[[i]]$fitted
    arcs(run)

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

plot_results(s2.val)
