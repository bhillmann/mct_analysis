library(phyloseq)
library(igraph)
library(SpiecEasi)


## Load round 2 of American gut project

mct.otu <- read.delim("results/otutable-subset-n50000-s10-norm-no-soylent-no-transition-pre-amg2-order.txt", sep="\t", row = 1, as.is=T)
mct.tax <- as.matrix(read.delim("results/taxatable-subset-n50000-s10-norm-no-soylent-no-transition-pre-amg2-order.txt", row.names =  1,  sep="\t", as.is=T))
mct.otu <- otu_table(mct.otu, taxa_are_rows = T)
mct.tax <- tax_table(mct.tax)
mct.phy <- phyloseq(mct.otu, mct.tax)

se.mb.mct <- spiec.easi(mct.phy, method='mb', lambda.min.ratio=1e-2, nlambda=20, icov.select.params=list(rep.num=50))
ig2.mb <- adj2igraph(se.mb.mct$refit,  vertex.attr=list(name=taxa_names(mct.phy)))

plot_network(ig2.mb, mct.phy, type='taxa', color="Order")
