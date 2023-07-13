#partial saga for converting victors database
#problem being legacy references to now removed records and inavailability of metadata
#13/07/23
library(tidyverse)
library(stringr)

#input file containing headers from victors_dna.fa
headers <- readLines("dna_heads.txt")
df <- data.frame(headers)

#generate new columns with extracted nucleotide accession, start and end pos
df <- df %>%
  mutate(accession = str_extract(headers, "\\|([^:]+)"),
         start_pos = str_extract(headers, "(?<=:)[0-9]+"),
         end_pos = str_extract(headers, "(?<=-)[0-9]+"))
#remove errant pipe symbols
df$accession <- str_replace_all(df$accession, "\\|", "")
# Function to retrieve accession number from nucleotide database
retrieve_accession <- function(accession_number) {
  query <- rentrez::entrez_fetch(db = "nucleotide", id = accession_number, rettype = "acc")
  query <- gsub("\n", "", query)
  return(query)
}
#check: retrieve_accession(110804074)
# Apply the function to the accession column in the dataframe
df$retrieved_accession <- lapply(df$accession, retrieve_accession)

#remove any na
df %>% dplyr::filter(grepl("Error",retrieved_accession)) -> missing
df %>% dplyr::filter(!grepl("Error", retrieved_accession)) -> db
db %>% dplyr::filter(is.na(start_pos)) -> na_values
db %>% dplyr::filter(!is.na(start_pos)) -> db
#flatten db and write as output
db <- apply(db,2,as.character)
write.csv(db,"accession_table.csv")

########################
# try to access genome #
########################
db <- read.csv("accession_table.csv")

# # Specify the accession, genomic range, and report type
# accession <- "NC_008438.1"
# from <- 2617
# to <- 3153
# report_type <- "genbank"
# record <- entrez_fetch(db = "nucleotide", id = accession, rettype = report_type, seq_start = from, seq_stop = to)
# file_string <- paste(record, collapse = " ")
# # Find the amino acid sequence using regular expressions
# protid <- str_extract(file_string, "(?<=/protein_id=\")[^\"/]+")

# Function to fetch protein ID and add it as a new column
fetch_protein_id <- function(df, accession, start_pos, end_pos) {
  # Specify the report type
  report_type <- "genbank"
  # Fetch protein ID for each row in the dataframe
  protein_ids <- apply(df, 1, function(row) {
    # Extract the accession, start position, and end position from the row
    accession <- row[[accession]]
    from <- row[[start_pos]]
    to <- row[[end_pos]]
    # Fetch the record
    record <- entrez_fetch(db = "nucleotide", id = accession, rettype = report_type, seq_start = from, seq_stop = to)
    # Convert the record to a string
    file_string <- paste(record, collapse = " ")
    # Find the protein ID using regular expressions
    protein_id <- str_extract(file_string, "(?<=/protein_id=\")[^\"/]+")
    # Return the protein ID
    protein_id
  })
  # Add the protein IDs as a new column in the dataframe
  df$protein_id <- protein_ids
  # Return the modified dataframe
  df
}
#apply function
db2 <- fetch_protein_id(db, "accession", "start_pos", "end_pos")
db <- apply(db,2,as.character)
write.csv(db2,"accession_table2.csv")
