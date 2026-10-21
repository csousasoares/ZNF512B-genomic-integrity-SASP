
# Load Packages -----------------------------------------------------------

renv::restore()

library(Rsubread)
library(Rsamtools)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readr)

# Check Data Structure ----------------------------------------------------

# Use any bam file:

bamfile <- "input_data/2_hisat2/SRR19647642_GSM6241345_siZNF512B_R2_Homo_sapiens_RNA-Seq.sorted.bam"

# Get header:

hdr <- scanBamHeader(bamfile)

# View chromosome names:

head(names(hdr[[1]]$targets), 20)

transcripts <- unique(names(hdr[[1]]$targets))
transcripts <- as.data.frame(transcripts)

transcripts |> head(20)


# Start FeatureCounts -----------------------------------------------------

annotations <- getInBuiltAnnotation(annotation = "hg38")

# Gets the built-in genome annotations. No need for GTF files.

# If needed, replace the annotations$Chr with "" instead of "chr":

unique(annotations$Chr)

# Now FeatureCounts:

bam_files <- list.files("input_data/2_hisat2/", # BAM folder.
                        pattern = "bam$", # Ends with .bam.
                        full.names = T) # Full directory name.
bam_files

fc <- featureCounts(
  files = bam_files,
  annot.ext = annotations,
  isGTFAnnotationFile = F,
  isPairedEnd = F,
  strandSpecific = 0,
  nthreads = 8
)

fc_df <- as.data.frame(fc$counts)

# Simplify Sample Names:

colnames_data <- stringr::str_split(
  colnames(fc_df), pattern = "_", simplify = T
) |> as.data.frame()

colnames_data$name <- paste0(colnames_data$V3, "_", colnames_data$V4)

colnames(fc_df) <- colnames_data$name

# ENTREZID TO Gene Symbol:

genes_names <- row.names(fc_df)
genes_names |> head()

annots <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys = genes_names, 
  columns = "SYMBOL", 
  keytype = "ENTREZID") 

annots <- na.omit(annots)

# Keep only genes mapped and replace with symbols:

fc_df_filtered <- fc_df[annots$ENTREZID,]

fc_df_filtered <- as.data.frame(fc_df_filtered)

row.names(fc_df_filtered) <- annots$SYMBOL

fc_df_filtered <- fc_df_filtered |> 
  tibble::rownames_to_column(var = "gene_symbol")

# Export Data -------------------------------------------------------------

dir.create("output_data/feature_counts", recursive = T)

readr::write_tsv(fc_df_filtered, 
                 "output_data/feature_counts/raw_counts.tsv", )

dir.create("output_data/session_info_logs")

writeLines(capture.output(sessionInfo()), "output_data/session_info_logs/01_featureCounts_sessionInfo.txt")
