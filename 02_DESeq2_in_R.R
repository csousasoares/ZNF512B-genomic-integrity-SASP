
# Load Packages -----------------------------------------------------------

library(tidyverse)
library(DESeq2)

# Load Data ---------------------------------------------------------------


raw_counts <- read.csv(
  "output_data/feature_counts/raw_counts.tsv",
  sep = "\t"
)

raw_counts |> head()

raw_counts <- raw_counts |> 
  column_to_rownames(var = "gene_symbol")

raw_counts |> head()


# Create Metadata ---------------------------------------------------------


metadata <- data.frame(
  "condition" = c("siZNF512B", "siZNF512B", "Control", "Control")
)

row.names(metadata) <- colnames(raw_counts)

