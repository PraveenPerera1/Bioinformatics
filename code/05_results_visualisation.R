# COMMIT 5: RESULTS VISUALISATION
#_______________________________________________________________________________

# Read the results saved during differential analysis
results <- read.csv(
  "results/All_Results.csv",
  stringsAsFactors = FALSE
)

dir.create("figures", showWarnings = FALSE)

# Retain valid values for plotting
plot_data <- results[
  is.finite(results$Expression_Difference) &
    is.finite(results$adj.P.Val),
]

stopifnot(nrow(plot_data) > 0)

# Use adjusted p-values to show the multiple-testing result
plot_data$Significant <- plot_data$adj.P.Val < 0.05

plot_data$MinusLog10FDR <- -log10(
  pmax(plot_data$adj.P.Val, .Machine$double.xmin)
)

point_colours <- ifelse(
  plot_data$Significant, "#D9534F", "#777777"
)

threshold <- -log10(0.05)

png(
  "figures/04_Effect_Significance_Plot.png",
  width = 1800, height = 1200, res = 180
)

par(mar = c(5, 5, 5, 2))

plot(
  plot_data$Expression_Difference,
  plot_data$MinusLog10FDR,
  pch = 16,
  cex = 0.5,
  col = point_colours,
  ylim = c(
    0,
    max(threshold, plot_data$MinusLog10FDR) * 1.2
  ),
  xlab = "Expression difference (High myopia - Emmetropia)",
  ylab = "-log10(BH-adjusted p-value)",
  main = "Expression differences and statistical evidence"
)

abline(h = threshold, col = "#D9534F", lty = 2)
abline(v = 0, col = "#444444", lty = 3)

mtext(
  sprintf(
    "%d of %d probes have adjusted p < 0.05",
    sum(plot_data$Significant),
    nrow(plot_data)
  ),
  side = 3,
  line = 0.4,
  cex = 0.85
)

dev.off()

#_______________________________________________________________________________

# Exploratory heatmap of the 30 top-ranked probes
library(pheatmap)

stopifnot(
  identical(colnames(expr), sample_table$Sample_ID)
)

# Rank by unadjusted p-value, with probe ID breaking ties
ranked_results <- results[
  order(results$P.Value, results$Probe_ID, na.last = NA),
]

top_probes <- head(ranked_results$Probe_ID, 30)

stopifnot(
  length(top_probes) == 30,
  all(top_probes %in% rownames(expr))
)

heatmap_values <- expr[top_probes, , drop = FALSE]

# Standardise each probe across samples for visual comparison
probe_sd <- apply(heatmap_values, 1, sd)
stopifnot(all(is.finite(probe_sd) & probe_sd > 0))

heatmap_z <- t(scale(t(heatmap_values)))

annotation <- data.frame(
  Group = sample_table$Group,
  row.names = sample_table$Sample_ID
)

annotation_colours <- list(
  Group = c(
    "Emmetropia" = "#2878B5",
    "High myopia" = "#D9534F"
  )
)

pheatmap(
  heatmap_z,
  scale = "none",
  annotation_col = annotation,
  annotation_colors = annotation_colours,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(100),
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 9,
  main = "Top 30 probes: exploratory, none significant at FDR < 0.05",
  filename = "figures/05_Top30_Exploratory_Heatmap.png",
  width = 11,
  height = 10
)

