LCA <- function(
    BLASTOUTPUT = BLASTOUTPUT, #path to tabular blast output in format "6 sscinames scomnames qseqid sseqid pident length mismatch gapopen qcovus qstart qend sstart send evalue bitscore staxids qlen qcovs"
    FASTA = FASTA, #from which to get the actual sequences, because for some reason they aren't in the blast results table
    DB_PATH_IN = NULL, # existing database to which to add these sequence annotations, if desired
    DB_PATH_OUT = NULL #if desired, name of database file to write out; may be the same as `DB_PATH_IN`, if appending to existing db
    ){

  suppressMessages(require(here))
  suppressMessages(require(tidyverse))
  suppressMessages(require(sys))
  suppressMessages(require(ShortRead))

  taxonkit_path <- "/opt/miniconda3/bin/taxonkit" # laptop
  
 blastoutput <- read_delim(BLASTOUTPUT, delim = "\t", col_names = F, show_col_types = FALSE)
  fasta <- readFasta(FASTA)

  # Initialize `db` as an empty dataframe with correct columns if no database exists
  if (!is.null(DB_PATH_IN) && file.exists(DB_PATH_IN)){
    db <- read.csv(DB_PATH_IN, row.names = 1) %>%
      mutate(BlastTaxIDs = as.character(BlastTaxIDs))
  } else {
    db <- data.frame(
      Hash = character(),
      BestTaxon = character(),
      HaplotypeNumber = numeric(),
      Tax_list = character(),
      Mismatches = character(),
      BlastTaxIDs = character(),
      LCA_taxid = character(),
      Kingdom = character(),
      Phylum = character(),
      Class = character(),
      Order = character(),
      Family = character(),
      Genus = character(),
      Species = character(),
      Sequence = character(),
      DateAdded = character(),
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(DB_PATH_OUT)){
    db_outfile <- DB_PATH_OUT
  }

  # Process blastoutput
  names(blastoutput) <- c("Taxon", "CommonName", "Hash", "Accession",
                          "pident", "length", "mismatch", "gapopen",
                          "qcovus", "qstart", "qend", "sstart",
                          "send", "evalue", "bitscore", "staxids", "qlen", "qcovs")

  result_list <- list()

  # Apply decision tree for each unique Hash
  for (hash_id in unique(blastoutput$Hash)) {
    hash_hits <- blastoutput %>%
      filter(Hash == hash_id)

    top_pident_hits <- hash_hits %>%
      filter(pident == 100)

    if (nrow(top_pident_hits) > 0) {
      bitscore_50th <- quantile(top_pident_hits$bitscore, 0.5)
      top_hit <- top_pident_hits[1, ]

      if (top_hit$bitscore >= bitscore_50th) {
        selected_hits <- top_pident_hits
      } else {
        selected_hits <- hash_hits %>%
          filter(pident > 99)
      }
    } else {
      selected_hits <- hash_hits %>%
        filter(pident > 99)
    }

    if (nrow(selected_hits) == 0) {
      selected_hits <- hash_hits %>%
        filter(pident > 97)
    }

    if (nrow(selected_hits) == 0) {
      selected_hits <- hash_hits %>%
        arrange(desc(bitscore)) %>%
        head(3)
    }

    result_list[[hash_id]] <- selected_hits
  }

  selected_data <- do.call(rbind, result_list)

  a <- selected_data %>%
    select(Hash, staxids) %>%
    drop_na() %>%
    unique() %>%
    group_by(Hash) %>%
    summarise(tax_list = paste(staxids, collapse = " ")) %>%
    mutate(tax_list = gsub(";", " ", tax_list))

  # Filter out hashes already in the db
  if (!is.null(DB_PATH_IN) && nrow(db) > 0){
    a <- a %>%
      filter(!Hash %in% db$Hash)
  }

  if (dir.exists(here("tmp")) == FALSE){
    system2("mkdir", shQuote(here("tmp")))
  }
  filepath <- here("tmp/taxids.txt")
  fileout <- here("tmp/taxids_lca.txt")
  lineagesout <- here("tmp/taxids_lineages.txt")

  write.table(a, file = filepath, row.names = F, col.names = F, quote = F, sep = "\t")
  exec_internal(taxonkit_path, args = c("lca", filepath, "-i", "2", "-o", fileout))
  exec_internal(taxonkit_path, args = c("reformat", fileout, "-I", "3", "-o", lineagesout))

  b <- read.table(lineagesout, sep ="\t") %>%
    separate(V4, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";")
  names(b)[1:3] <- c("Hash", "BlastTaxIDs", "LCA_taxid")

  d <- blastoutput %>%
    select(Hash, Taxon) %>%
    drop_na() %>%
    unique() %>%
    group_by(Hash) %>%
    summarise(Tax_list = paste(Taxon, collapse = ", "))

  e <- blastoutput %>%
    select(Hash, mismatch) %>%
    drop_na() %>%
    unique() %>%
    group_by(Hash) %>%
    summarise(Mismatches = paste(mismatch, collapse = ", "))

  f <- d %>%
    left_join(e, by = "Hash") %>%
    left_join(b, by = "Hash") %>%
    filter(Family != "")  %>%
    mutate(BestTaxon = case_when(
      Species != "" ~ Species,
      Genus != "" & Species == "" ~ Genus,
      Family != "" & Genus == "" & Species == "" ~ Family,
      TRUE ~ NA_character_
    )) %>%
    group_by(BestTaxon) %>%
    mutate(
      MaxHaplotypeNumber = if_else(
        BestTaxon %in% db$BestTaxon,
        coalesce(
          as.numeric(max(db$HaplotypeNumber[db$BestTaxon == BestTaxon], na.rm = TRUE)),
          0
        ),
        0
      ),
      MaxHaplotypeNumber = replace_na(MaxHaplotypeNumber, 0),
      HaplotypeNumber = row_number() + MaxHaplotypeNumber
    ) %>%
    mutate(HaplotypeNumber = if_else(MaxHaplotypeNumber == 0, row_number(), HaplotypeNumber)) %>%
    select(-MaxHaplotypeNumber) %>%
    relocate(Hash, BestTaxon, HaplotypeNumber, Tax_list, Mismatches)

  f <- f %>%
    mutate(
      HaplotypeNumber = as.numeric(HaplotypeNumber),
      HaplotypeNumber = if_else(is.na(HaplotypeNumber) | HaplotypeNumber == "#NAME?", NA_real_, HaplotypeNumber)
    ) %>%
    filter(
      !is.na(HaplotypeNumber) &
        !is.na(BestTaxon) &
        !is.na(Hash) &
        !is.na(BlastTaxIDs)
    )

  f$Sequence <- sread(fasta[match(f$Hash, ShortRead::id(fasta))]) %>% as.character()

  f$DateAdded <- as.character(Sys.Date())

  if (!is.null(DB_PATH_IN) && nrow(db) > 0){
    g <- db %>%
      bind_rows(f)
  } else {
    g <- f
  }

  if (!is.null(DB_PATH_OUT)){
    write.csv(g, file = DB_PATH_OUT)
  } else {
    return(f)
  }
}
