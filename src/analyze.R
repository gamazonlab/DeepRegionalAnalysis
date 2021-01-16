a = read.table('../prediction_performance/r2.txt', header=T)

b = read.table('../data/EvoStats_full.txt', header=T)

b$NEAN = b$GeneName %in% a$genename

b.1 = subset(b, b$Type=="protein_coding")
dim(b.1)
# [1] 19644    21

wilcox.test(subset(b, b$NEAN)$tau, subset(b, !b$NEAN)$tau)
# W = 2346800, p-value = 1.053e-10  # this includes all non-coding and protein-coding genes

wilcox.test(subset(b.1, b.1$NEAN)$tau, subset(b.1, !b.1$NEAN)$tau)
# W = 979330, p-value = 0.006321  # protein-coding genes only



