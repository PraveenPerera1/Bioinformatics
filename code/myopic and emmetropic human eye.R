library(GEOquery)
library(Biobase)

options(timeout = 600)

# Expected GEO Series Matrix address
matrix_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE136nnn/GSE136701/matrix/",
  "GSE136701_series_matrix.txt.gz"
)

# Download to a temporary local file
matrix_file <- tempfile(fileext = ".txt.gz")

download.file(
  url = matrix_url,
  destfile = matrix_file,
  method = "libcurl",
  mode = "wb"
)



# Read the downloaded file
eset <- getGEO(
  filename = matrix_file,
  getGPL = FALSE
)

# Extract and inspect expression values
expr <- exprs(eset)

dim(expr)
annotation(eset)

# Inspect sample labels
metadata <- pData(eset)

metadata[, grepl(
  "title|source_name|characteristics",
  colnames(metadata)
), drop = FALSE]


#________________________________________________________________________________

file.exists(matrix_file)
file.info(matrix_file)$size
#___________________reading the dataset_________________________________________

# Check expression matrix dimensions and platform
dim(expr)
annotation(eset)

# Inspect sample groups
sample_info <- metadata[, grepl(
  "title|source_name|characteristics",
  colnames(metadata)
), drop = FALSE]

print(sample_info)

# Check expression values
sum(is.na(expr))
sum(is.infinite(expr))
summary(as.vector(expr))

# Confirm samples are in the same order
identical(colnames(expr), rownames(metadata))



