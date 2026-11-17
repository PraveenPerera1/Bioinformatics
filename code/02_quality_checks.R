# Sample groups and quality checks
library(Biobase)

# Confirm sample alignment
stopifnot(identical(colnames(expr), rownames(metadata)))

# Extract diagnosis from the verified metadata column
diagnosis <- trimws(sub(
  "^diagnosis:\\s*", "",
  metadata$characteristics_ch1.3
))

# Check that all diagnoses match the expected groups
stopifnot(all(diagnosis %in% c("Emmetropia", "High myopia")))

group <- factor(
  diagnosis,
  levels = c("Emmetropia", "High myopia")
)

# Create a sample information table
sample_table <- data.frame(
  Sample_ID = rownames(metadata),
  Group = group,
  Age = as.numeric(sub(
    "^age: ([0-9]+).*", "\\1",
    metadata$characteristics_ch1.2
  ))
)

print(sample_table)
print(table(group))

# Verify expression values
stopifnot(
  ncol(expr) == nrow(sample_table),
  all(is.finite(expr))
)

# Save verified sample information
dir.create("results", showWarnings = FALSE)

write.csv(
  sample_table,
  "results/Sample_Information.csv",
  row.names = FALSE
)

# Inspect the processing documented by the original researchers
processing_columns <- grep(
  "data_processing",
  colnames(metadata),
  value = TRUE
)

print(unique(metadata[, processing_columns, drop = FALSE]))