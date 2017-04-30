library(phyloseq)

data('amgut2.filt.phy')

write.table(tax_table(amgut2.filt.phy), "results/taxtable-amgut2-filt.txt", sep = "\t")
