##############################################################################
##                         Welcome to The Kelly Lab's                       ##
##              Metabarcoding Taxonomic Assignment Pipeline                 ##
## This pipeline will take you from de-multiplexed fastq files to a matrix  ##
## with taxon names and counts. It should work for linux or mac. PC users,  ##
## look into WSL2.                                                          ##
##############################################################################
##    Before you start, ensure you have access to the following files.      ##
##            Many can be found in the OneDrive under:                      ##
##             "2. KellyLab/Projects/AnnotationDatabases/code/"             ##
##                                                                          ##
## Required Files:                                                          ##
## 1. Your Fastqs, of course. Their names must start with the primer name   ##
## 2. CEG_BLAST_function.R - Function to run BLAST on the CEG cluster       ##
## 3. LCA_function.R - Function to run LCA                                  ##
## 4. Pipeline_code.R - This very code                                      ##
## 5. MURIblast_new_seqs.sh                                                 ##
## 6. MURIblast_*_template.sh - Code for taxonomic assignment               ##
##    depending on delimitation                                             ##
## 7. make_primer_shell_script.R - A script that writes another script      ##
## 8. primer_data.csv - Sheet with primer information like sequence, etc    ##
##                                                                          ##
## Required Programs:                                                       ##
## 1. Taxonkit: https://bioinf.shenwei.me/taxonkit/                         ##
## 2. Cutadapt: https://cutadapt.readthedocs.io/en/stable/installation.html ##
##############################################################################

#Restart your Rstudio session before running each time.
#You may also clean your environment if you want to. Be careful when running this
#rm(list = ls())
#.rs.restartR()

##put the path to this very code here.
here::i_am("pedro_pipeline_code-NVP-intercal-MFU.R")

#load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(dada2))
suppressMessages(library(digest))
suppressMessages(library(seqinr))
suppressMessages(library(sys))
suppressMessages(library(ShortRead)) 
suppressMessages(library(here))

#Your R working directory. This will change again downstream 
setwd(here("/Users/nastassiapatin/OneDrive - SCCWRP/Experiments/Intercalibration/MFU"))

#the folder in which a subfolder, called `Fastq` lives, and which contains only fastq or fastq.gz files
PARENT_LOCATION <- "/Volumes/CalCOFI_OMICS/Sequence_data/UCSD_NovaSeq_20240403/intercal-MFU/"
  if (strsplit(PARENT_LOCATION,"")[[1]][nchar(PARENT_LOCATION)] != "/"){PARENT_LOCATION = paste0(PARENT_LOCATION, "/")} #ensures a trailing "/" as the last character on the path here

RUN_NAME <- "Intercal-MFU" #basename(PARENT_LOCATION) #

PRIMERNAME <- "MFU" #used to look up primer information, so needs to be one of a set number of known primers
#if there are multiple primers, you run the pipeline multiple times, once for each primer. This is why fastqs have to be named a certain way

CLEANUP <- FALSE #logical; delete all intermediate files and just keep logs and final outputs?
BLASTSCOPE <- "vertebrate" #"eukaryote" ##either "vertebrate" or "eukaryote" .  Default is 97% identity for vertebrate, 90% for eukaryote

#where the hash database is. This is to check if sequences have been seen before and avoid double effort in taxon assignment
DATABASE_LOCATION <- "/Volumes/CalCOFI_OMICS/Databases/MURI hash DBs"
  if (strsplit(DATABASE_LOCATION,"")[[1]][nchar(DATABASE_LOCATION)] != "/"){DATABASE_LOCATION = paste0(DATABASE_LOCATION, "/")} #ensures a trailing "/" as the last character on the path here

#then, ultimately, where do you want the small handful of relevant output files to live?
PROCESSED_LOCATION <- "/Users/nastassiapatin/OneDrive - SCCWRP/Experiments/Intercalibration/MFU/"

#system2("mkdir", shQuote(paste0(PROCESSED_LOCATION, "/", PRIMERNAME, "_", RUN_NAME))) #create folder for code, etc
#FASTQ_LOCATION <- paste0(PARENT_LOCATION, "fastq") #folder within Parent Location
CODE_LOCATION <- paste0(PROCESSED_LOCATION, "code_etc") #folder within Parent Location
  system2("mkdir", shQuote(CODE_LOCATION)) #create folder for code, etc
TRIMMED_LOCATION <- paste0(PARENT_LOCATION, "cutadapt-trimmed-v2") #folder within Parent Location
FILTERED_LOCATION <- paste0(PARENT_LOCATION, "dada2_filtered-v2") #folder within Parent Location
OUTPUT_LOCATION <- paste0(PROCESSED_LOCATION, "outputs/") #folder within Parent Location
  system2("mkdir", args=shQuote(FILTERED_LOCATION)) #make folder for processed reads
  system2("mkdir", args=shQuote(OUTPUT_LOCATION)) #make folder for pipeline outputs

#check dependencies and get path for taxonkit
TAXONKIT_PATH = "/opt/miniconda3/bin/taxonkit" # office

#check if taxonkit is working
system2(TAXONKIT_PATH, args=shQuote("version"))

#calling function for LCA
source("/Users/nastassiapatin/OneDrive - SCCWRP/Bioinformatics/MURI metabarcoding/Accessory_scripts/LCA_function_pedro_1.3-NVP.R")

#read in primer data
primer.data <- read.csv("/Users/nastassiapatin/OneDrive - SCCWRP/Bioinformatics/MURI metabarcoding/Accessory_scripts/primer.data.csv")

#-----------------------------------------------------------------------------------------------
#Downstream of here, no further user input is required. 
#Therefore, you may select all lines up to "ANNOTATION" and run all at once if you wish to. 
#Good luck!-------------------------------------------------------------------------------------

#set up relevant info for primer trimming
PRIMERSEQ_F <- primer.data %>% filter(name == PRIMERNAME) %>% pull(seq_f)
PRIMERSEQ_R <- primer.data %>% filter(name == PRIMERNAME) %>% pull(seq_r)
PRIMERLENGTH_F <- primer.data %>% filter(name == PRIMERNAME) %>% pull(primer_length_f)
PRIMERLENGTH_R <- primer.data %>% filter(name == PRIMERNAME) %>% pull(primer_length_r)
MAX_AMPLICON_LENGTH <- primer.data %>% filter(name == PRIMERNAME) %>% pull(max_amplicon_length)
MIN_AMPLICON_LENGTH <- primer.data %>% filter(name == PRIMERNAME) %>% pull(min_amplicon_length)
OVERLAP <- primer.data %>% filter(name == PRIMERNAME) %>% pull(overlap)

###COPY FILES TO SCRIPTS FOLDER and write primer-trimming script via make_primer_shell_script.R

    write.csv(primer.data, paste0(CODE_LOCATION, "/primer.data.csv"))
    
    if (BLASTSCOPE == "vertebrate"){
      system2("cp", args = c(shQuote("/Users/nastassiapatin/OneDrive - SCCWRP/Bioinformatics/MURI metabarcoding/Accessory_scripts/MURIblast_vertebrate_template-NVP-localnteuk.sh"), 
                             shQuote(CODE_LOCATION))) #template shell script for blast  
    }
    if (BLASTSCOPE == "eukaryote"){
      system2("cp", args = c(shQuote(here("MURIblast_eukaryote_template.sh")), shQuote(CODE_LOCATION))) #template shell script for blast  
    }
    

###RUN DADA2---------------------------------------------------------------------

filelist <- system2("ls", args = shQuote(TRIMMED_LOCATION), stdout = TRUE)

fnFs <- filelist[str_detect(filelist, "_R1-v2")] #filelist[grep(pattern="_R1_001.fastq", filelist)]
fnRs <- filelist[str_detect(filelist, "_R2-v2")] #filelist[grep(pattern="_R2_001.fastq", filelist)]
  
sample.names <- sapply(strsplit(basename(fnFs), "_12S"), `[`, 1)

### Name filtered files in filtered/subdirectory ----------------------------------

filtFs <- file.path(FILTERED_LOCATION, 
                    paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(FILTERED_LOCATION, 
                    paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### Filter out Empty Samples -----------------------------------------------------
# if we don't filter out empty samples, an error happens during finding qual trimming length

setwd(TRIMMED_LOCATION)
file.empty <- function(filenames) file.info(filenames)$size == 20
empty_files <- file.empty(fnFs) | file.empty(fnRs)
fnFs <- fnFs[!empty_files]
fnRs <- fnRs[!empty_files]
filtFs <- filtFs[!empty_files]
filtRs <- filtRs[!empty_files]
sample.names <- sample.names[!empty_files]

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     trimRight = c(PRIMERLENGTH_R,PRIMERLENGTH_F),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE)

### Dereplicate ---------------------------------------------------------------

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

# Name the derep-class objects by the sample names

sample.names <- sample.names[exists]

### Learn Error Rates ---------------------------------------------------------------

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

##It is possible it will not converge here. If that is the case, check if it was at least close:
#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)
#dada2:::checkConvergence(dadaRs[[1]])

### Sample Inference -----------------------------------------------------------------------

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE) 

### Merge Paired Reads ---------------------------------------------------------------------

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, 
                      minOverlap = OVERLAP, verbose=TRUE, 
                      trimOverhang = TRUE)

### Construct sequence table ---------------------------------------------------------------

seqtab <- makeSequenceTable(mergers)

### Remove chimeras ------------------------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)

freq.nochim <- sum(seqtab.nochim)/sum(seqtab)

### Filter by Size -------------------------------------------------------------------------

indexes.to.keep <- which((nchar(colnames(seqtab.nochim)) <= MAX_AMPLICON_LENGTH) & 
                           ((nchar(colnames(seqtab.nochim))) >= MIN_AMPLICON_LENGTH))
cleaned.seqtab.nochim <- seqtab.nochim[, indexes.to.keep]
filteredout.seqtab.nochim <- seqtab.nochim[, !indexes.to.keep]
write.csv(filteredout.seqtab.nochim, paste0(OUTPUT_LOCATION,
                                           "filtered_out_asv.csv"))

### Track reads through pipeline -----------------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out[exists,],
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim),
               rowSums(cleaned.seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", 
                     "denoisedR", "merged", "nonchim", "length_filter")
rownames(track) <- sample.names
head(track)
write.csv(track, paste0(OUTPUT_LOCATION, "tracking_reads.csv"))

### Create Hashing  ------------------------------------------------------------------------

# define output files

conv_file <- paste0(OUTPUT_LOCATION, 
                    paste0(RUN_NAME,"_hash_key.csv"))
conv_file.fasta <- file.path(OUTPUT_LOCATION, 
                             paste0(RUN_NAME,"_hash_key.fasta"))
ASV_file <-  file.path(OUTPUT_LOCATION, 
                       paste0(RUN_NAME,"_ASV_table.csv"))
taxonomy_file <- file.path(OUTPUT_LOCATION, 
                           paste0(RUN_NAME,"_taxonomy_output.csv"))
bootstrap_file <- file.path(OUTPUT_LOCATION, 
                            paste0(RUN_NAME,"_tax_bootstrap.csv"))

# create ASV table and hash key 
print(paste0("creating ASV table and hash key...", Sys.time()))
seqtab.nochim.df <- as.data.frame(cleaned.seqtab.nochim)
conv_table <- tibble( Hash = "", Sequence ="")
Hashes <- map_chr (colnames(seqtab.nochim.df), 
                   ~ digest(.x, algo = "sha1", serialize = F, 
                            skip = "auto"))
conv_table <- tibble (Hash = Hashes,
                      Sequence = colnames(seqtab.nochim.df))

write_csv(conv_table, conv_file) # write hash key into a csv
write.fasta(sequences = as.list(conv_table$Sequence), # write hash key into a fasta
            names     = as.list(conv_table$Hash),
            file.out = conv_file.fasta)
sample.df <- tibble::rownames_to_column(seqtab.nochim.df,"Sample_name")
sample.df <- data.frame(append(sample.df,c(Label=PRIMERNAME), after = 1))
current_asv <- bind_cols(sample.df %>%
                           dplyr::select(Sample_name, Label),
                         seqtab.nochim.df)
current_asv <- current_asv %>%
  pivot_longer(cols = c(- Sample_name, - Label),
               names_to = "Sequence",
               values_to = "nReads") %>%
  filter(nReads > 0)
current_asv <- merge(current_asv, conv_table, by="Sequence") %>%
  select(-Sequence) %>%
  relocate(Hash, .after=Label)

write_csv(current_asv, ASV_file) # write asv table into a csv

###ANNOTATION------------------------------------------------------------------------------------------

##Looking at the existing hash database to check if any observed sequences have been seen before.
if (file.exists(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))){
  db <- read.csv(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"), row.names = 1)
  seen <- which(conv_table$Hash %in% db$Hash) #hashes in this run that already occur in the database 
  notseen <- which(!conv_table$Hash %in% db$Hash) #hashes in this run that have not previously been annotated
  
  #write new (as-yet-unannotated sequences) to fasta
  write.fasta(sequences = as.list(conv_table$Sequence[notseen]), # write hash key into a fasta
              names     = as.list(conv_table$Hash[notseen]),
              file.out = paste0(OUTPUT_LOCATION, "seqs_to_annotate.fasta"))
#if you want to blast everything instead, only run these lines below
  } else {
  write.fasta(sequences = as.list(conv_table$Sequence), # write hash key into a fasta
              names     = as.list(conv_table$Hash),
              file.out = paste0(OUTPUT_LOCATION, "seqs_to_annotate.fasta"))
}

### Local BLAST against nt_euk database ###
system2("sh", args = shQuote(paste0(CODE_LOCATION, 
                                    "/MURIblast_vertebrate_template-NVP-localnteuk.sh")))

# Import BLAST outputs and sequence fasta
BLASTOUTPUT = "/Users/nastassiapatin/OneDrive - SCCWRP/Experiments/Intercalibration/MFU/outputs/BLAST_output.txt"
FASTA = "/Users/nastassiapatin/OneDrive - SCCWRP/Experiments/Intercalibration/MFU/outputs/seqs_to_annotate.fasta"
#TAX = paste0(PROCESSED_LOCATION, "/", PRIMERNAME, "_", 
 #            RUN_NAME, "/outputs/taxids_taxonomy.txt")

##Taxonomic assignment

if (file.size(paste0(OUTPUT_LOCATION, "BLAST_output.txt")) == 0L){
  db <- read.csv(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"), 
                 row.names = 1)
} else if (file.exists(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))){
  #LCA table for output, and add to database
  LCA(BLASTOUTPUT = paste0(OUTPUT_LOCATION, "/", "BLAST_output.txt"),
      FASTA = paste0(OUTPUT_LOCATION, "/", "seqs_to_annotate.fasta"),
      DB_PATH_IN = paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"),
      DB_PATH_OUT = paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))
  
  # re-load newly updated db
  db <- read.csv(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"), row.names = 1)
  #if you have blasted everything again instead, only run these lines below  
} else {
  db <- LCA(BLASTOUTPUT = BLASTOUTPUT,
            FASTA = FASTA)
  write.csv(db %>% distinct(), paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))
}

tax_table <- current_asv %>%
  left_join(db %>% dplyr::select(Hash, BestTaxon)) %>%
  group_by(Sample_name, BestTaxon) %>%
  summarise(nReads = sum(nReads))

write_csv(tax_table, file = paste0(OUTPUT_LOCATION,
                                   "taxon_table.csv"))

tax_table %>%
  pivot_wider(id_cols = BestTaxon, names_from = Sample_name, 
              values_from = nReads, values_fill = 0) %>% 
  write_csv(file = paste0(OUTPUT_LOCATION, "taxon_table_wide.csv"))

message("Pipeline complete. Outputs are now available in the Processed Location.")

#OPTIONAL SECTION-----------------------------------------------------------------------

### (Optional) 
#A taxon table but that lists haplotypes separately

hash_database <- read_csv(file.path(paste0(PROCESSED_LOCATION, "/", 
                                           PRIMERNAME, "_", RUN_NAME, 
                                           "/outputs/", PRIMERNAME,
                                           "_database.csv")))
hash_info <- hash_database %>%
  select(Hash, BestTaxon, Tax_list, LCA_taxid, Kingdom, Phylum,
         Class, Order, Family, Genus, Species)
#str(hash_database)
#str(current_asv)
merged_data <- current_asv %>%
  left_join(hash_info, by = "Hash")
#str(merged_data)

# Reshape the data to have Taxons as columns and samples as rows
reshaped_data <- merged_data %>%
  pivot_wider(names_from = Sample_name, 
              values_from = nReads,
              values_fill = 0) # id_cols = BestTaxon, values_fill = 0

# Export the reshaped data to a file if needed
write.csv(reshaped_data, paste0(PROCESSED_LOCATION, "/", 
                                PRIMERNAME, "_", RUN_NAME, 
                                "/outputs/ASV_table_with_tax.csv"), 
          row.names = FALSE)

###CLEANUP
#Note this will throw a shell-init error after completing, because it deletes the directory in which some active code is sitting (?)
if (CLEANUP){
  for (j in c("code_etc", "outputs", "filtered", "for_dada2", "logs")){
    system2("rm", args = c("-r", shQuote(paste0(PARENT_LOCATION, j))))
  }
}
