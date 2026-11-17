# ============================================================
# COMMIT 3: EXPLORATORY ANALYSIS
# Use the GEO expression values as supplied.
# No additional transformation or normalisation is applied.
# ============================================================

dir.create("figures", showWarnings = FALSE)

# Confirm sample alignment and valid values
stopifnot(
  identical(colnames(expr), sample_table$Sample_ID),
  all(is.finite(expr))
)

# Use consistent colours for the two groups
group_colours <- c(
  "Emmetropia" = "#2878B5",
  "High myopia" = "#D9534F"
)

sample_colours <- unname(
  group_colours[as.character(sample_table$Group)]
)

# 1. Compare expression distributions across samples_______________________________
png(
  "figures/01_Expression_Boxplot.png",
  width = 1800, height = 1200, res = 180
)

par(mar = c(7, 5, 4, 2))

boxplot(
  expr,
  names = sample_table$Sample_ID,
  col = sample_colours,
  las = 2,
  outline = FALSE,
  main = "Expression distributions across samples",
  ylab = "Expression value (as supplied by GEO)"
)

legend(
  "topright",
  legend = names(group_colours),
  fill = group_colours,
  bty = "n"
)

dev.off()

# 2. PCA: samples are observations; probes are variables
# Centre each probe without scaling it to unit variance.
pca <- prcomp(t(expr), center = TRUE, scale. = FALSE)

variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)

png(
  "figures/02_PCA_Plot.png",
  width = 1800, height = 1200, res = 180
)

plot(
  pca$x[, 1], pca$x[, 2],
  col = sample_colours,
  pch = 19,
  cex = 1.5,
  xlim = extendrange(pca$x[, 1], f = 0.15),
  ylim = extendrange(pca$x[, 2], f = 0.15),
  xlab = sprintf("PC1 (%.1f%%)", variance_percent[1]),
  ylab = sprintf("PC2 (%.1f%%)", variance_percent[2]),
  main = "PCA of lens epithelial samples"
)




text(
  pca$x[, 1], pca$x[, 2],
  labels = sample_table$Sample_ID,
  pos = 3,
  cex = 0.7
)

legend(
  "topright",
  legend = names(group_colours),
  col = group_colours,
  pch = 19,
  bty = "n"
)

axis(
  side = 2,
  at = pretty(c(0, max(sample_clusters$height))),
  las = 1
)

dev.off()

print(round(variance_percent[1:2], 2))





# 3. Hierarchical clustering using all supplied probes___________________________
sample_distances <- dist(t(expr), method = "euclidean")

sample_clusters <- hclust(
  sample_distances,
  method = "complete"
)

png(
  "figures/03_Sample_Clustering.png",
  width = 1800, height = 1200, res = 180
)

par(mar = c(10, 5, 4, 2))

plot(
  sample_clusters,
  labels = paste(
    sample_table$Sample_ID,
    sample_table$Group,
    sep = " - "
  ),
  main = "Hierarchical clustering of lens epithelial samples",
  xlab = "",
  sub = "",
  ylab = "Euclidean distance (complete linkage)",
  hang = -1,
  cex = 0.8
)

dev.off()