
# COMMIT 4: EXPLORATORY DIFFERENTIAL EXPRESSION
# Contrast: High myopia minus Emmetropia
# Use GEO values unchanged; exact transformation is undocumented.
#_______________________________________________________________________________

library(limma)

dir.create("results", showWarnings = FALSE)

# Confirm alignment and valid expression values
stopifnot(
  identical(colnames(expr), sample_table$Sample_ID),
  all(is.finite(expr))
)

group <- factor(
  sample_table$Group,
  levels = c("Emmetropia", "High myopia")
)

stopifnot(
  !anyNA(group),
  all(table(group) == 3)
)

# Create separate mean estimates for each group
design <- model.matrix(~ 0 + group)
colnames(design) <- c("Emmetropia", "HighMyopia")
rownames(design) <- sample_table$Sample_ID

print(design)

# Positive effect means higher expression in high myopia
contrast_matrix <- makeContrasts(
  HighMyopia_vs_Emmetropia = HighMyopia - Emmetropia,
  levels = design
)

# Fit probe-wise linear models and moderate their variances
fit <- lmFit(expr, design)
fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_final <- eBayes(fit_contrast)

# Correct for multiple testing using Benjamini-Hochberg
results <- topTable(
  fit_final,
  coef = "HighMyopia_vs_Emmetropia",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

results$Probe_ID <- rownames(results)

# limma calls this column logFC, but that interpretation
# requires a confirmed log2 input scale.
names(results)[names(results) == "logFC"] <-
  "Expression_Difference"

results <- results[, c(
  "Probe_ID", "Expression_Difference", "AveExpr",
  "t", "P.Value", "adj.P.Val", "B"
)]

# Select probes using adjusted p-values only
significant <- results[
  !is.na(results$adj.P.Val) & results$adj.P.Val < 0.05,
]

write.csv(
  results,
  "results/All_Results.csv",
  row.names = FALSE
)

write.csv(
  significant,
  "results/Significant_Probes.csv",
  row.names = FALSE
)

cat("\nProbes tested:", nrow(results))
cat("\nProbes with adjusted p < 0.05:", nrow(significant), "\n")

print(head(results, 10))