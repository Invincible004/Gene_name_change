# Required libraries
library(data.table)

# Function to parse gene_info file and create mapping of symbols to Entrez IDs
parse_gene_info <- function(gene_info_file) {
  gene_info <- fread(gene_info_file, sep="\t", header=TRUE, data.table=FALSE, na.strings="", quote="", fill=TRUE)
  symbol_to_entrez <- list()

  for (i in 1:nrow(gene_info)) {
    entrez_id <- gene_info[i, "GeneID"]
    symbol <- gene_info[i, "Symbol"]
    synonyms <- unlist(strsplit(gene_info[i, "Synonyms"], "\\|"))

    symbol_to_entrez[[symbol]] <- entrez_id
    for (synonym in synonyms) {
      symbol_to_entrez[[synonym]] <- entrez_id
    }
  }

  return(symbol_to_entrez)
}

# Function to replace gene names in the GMT file with Entrez IDs
replace_gene_names <- function(gmt_file, symbol_to_entrez) {
  gmt_data <- fread(gmt_file, sep="\t", header=FALSE, data.table=FALSE, na.strings="", quote="", fill=TRUE)
  new_gmt_data <- list()

  for (i in 1:nrow(gmt_data)) {
    pathway_name <- gmt_data[i, 1]
    pathway_description <- gmt_data[i, 2]
    genes <- gmt_data[i, -c(1, 2)]

    new_genes <- sapply(genes, function(gene) {
      if (gene %in% names(symbol_to_entrez)) {
        return(symbol_to_entrez[[gene]])
      } else {
        return(gene)
      }
    })

    new_row <- c(pathway_name, pathway_description, new_genes)
    new_gmt_data[[i]] <- new_row
  }

  return(data.table::as.data.table(do.call(rbind, new_gmt_data)))
}

# Function to write the new GMT file with Entrez IDs
write_new_gmt <- function(new_gmt_data, output_file) {
  fwrite(new_gmt_data, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
  library(data.table)
}

# Input file names
gene_info_file <- "Homo_sapiens.gene_info.gz"
gmt_file <- "h.all.v2023.1.Hs.symbols.gmt"
output_file <- "output.gmt"

# Step 1: Parse gene_info file and create mapping of symbols to Entrez IDs
symbol_to_entrez <- parse_gene_info(gene_info_file)

# Step 2: Replace gene names in the GMT file with Entrez IDs
new_gmt_data <- replace_gene_names(gmt_file, symbol_to_entrez)

# Step 3: Write the new GMT file with Entrez IDs
write_new_gmt(new_gmt_data, output_file)
