
# DADA2 Course script (v1.8)

## Input: Illumina paired-end FASTQs (demultiplexed), adapters/primers already removed (Cutadapt done)
#
# Output: Quality control plots, filtering summaries, read-tracking table, chimera-free ASV table + sequences

## Optional: taxonomy assignment, phyloseq object





### load required packages

library(tidyverse)
library(dada2)
library(ggplot2)
library(phyloseq)
library(Biostrings)
library(kableExtra)
library(readxl)
library(ShortRead)
library(stringr)
library(vegan)
library(tibble)
library(dplyr)
library(writexl)
library(gridExtra)

#These are essential: 
library(tidyverse)
library(dada2)
library(ggplot2)
library(dplyr)
library(Biostrings)
library(ShortRead)
library(stringr)
library(tibble)



# This block installs any missing packages automatically.
# Safe to re-run anytime.

#--- Define helper function ---#
install_load_package <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  library(pkg, character.only = TRUE)
}

#--- CRAN packages ---#
cran_packages <- c(
  "tidyverse", "ggplot2", "kableExtra", "readxl", "stringr",
  "vegan", "tibble", "dplyr", "writexl", "gridExtra"
)

#--- Bioconductor packages ---#
bioc_packages <- c("dada2", "ShortRead", "Biostrings", "phyloseq")

#--- Install & load all ---#
sapply(cran_packages, install_load_package)
sapply(bioc_packages, install_load_package, bioc = TRUE)

###########################################################################

#### Prepare Directories ####################################################
# BEFORE RUNNING THIS SCRIPT:
# 1. In your working directory, create a folder called "DADA2_output"
# 2. Inside that folder, create another folder called "fastq"
# 3. Put all your FASTQ.GZ files in:
#       DADA2_output/fastq/
#############################################################################


# Check current working directory:
getwd()
 # if you need to fix - setwd("C:/Users/Documents/DADA2_project")  -> edit this to path in your PC

# Define the name of directories to use.
# These will be created in your current working directory
fastq_dir <- "DADA2_output/fastq"  # fastq directory with the samples that will be used
filtered_dir <- "DADA2_output/fastq_filtered/"  # for the fastq-files after filtering
qual_dir <- "DADA2_output/qual_pdf/"  # quality scores plots
dada2_dir <- "DADA2_output/dada2_results/"  # final dada2 results


# Create the directories
dir.create(filtered_dir)
dir.create(qual_dir)
dir.create(dada2_dir)



#### 1) First, we will check that primers/adapters have been already removed from our reads. Run all code below until step 2)

# List the used Primer sequences
fwd_primer <- "ATGCGATACTTGGTGTGAAT"
rev_primer <- "TCCTCCGCTTATTGATATGC"

# Path to fastq files
fastq_path <- "DADA2_output/fastq" 

# List forward and reverse reads
fwd_files <- list.files(fastq_path, pattern = "_R1.*\\.fastq\\.gz$", full.names = TRUE)
rev_files <- list.files(fastq_path, pattern = "_R2.*\\.fastq\\.gz$", full.names = TRUE)

# Sort to align pairs
fwd_files <- sort(fwd_files)
rev_files <- sort(rev_files)

# Function to count primer hits in a file
count_primer_hits <- function(file, primer) {
  fq <- readFastq(file)
  seqs <- sread(fq)
  sum(vcountPattern(primer, seqs, fixed = TRUE) > 0)
}

# Create a summary data frame
primer_summary <- data.frame(
  Sample = character(),
  Fwd_Primer_Hits = integer(),
  Rev_Primer_Hits = integer(),
  stringsAsFactors = FALSE
)

# Loop through file pairs
for (i in seq_along(fwd_files)) {
  fwd_file <- fwd_files[i]
  rev_file <- rev_files[i]
  
  sample_name <- str_replace(basename(fwd_file), "_R1.*", "")
  
  fwd_hits <- count_primer_hits(fwd_file, fwd_primer)
  rev_hits <- count_primer_hits(rev_file, rev_primer)
  
  primer_summary <- rbind(primer_summary, data.frame(
    Sample = sample_name,
    Fwd_Primer_Hits = fwd_hits,
    Rev_Primer_Hits = rev_hits
  ))
}

# View results 
## should show 0 - our data is ready to go (primers removed)
print(primer_summary)


#####  
# 2) Examining files, setting correct file names and tyding up

#### Examine fastq files ####
# get a list of all fastq files in the fastq" directory and separate R1 and R2
fns <- sort(list.files(fastq_dir, full.names = T))
fns <- fns[str_detect(basename(fns), ".fastq.gz")]
fns_R1 <- fns[str_detect(basename(fns), "R1")]
fns_R2 <- fns[str_detect(basename(fns), "R2")]
list.files(fastq_dir)

# Extract sample names, assuming filenames have format: RL-XXXX_SXXX_L001_RX_trimmed.fastq.gz
sample.names <- str_split(basename(fns_R1), pattern = "_", simplify = TRUE)
sample.names <- sample.names[, 1] 

#### Prepare data: Make a dataframe with the number of sequences in each file ####
# Empty files (i.e files without sequences) will be removed
df <- data.frame()
for (i in seq_along(fns_R1)) {
  # skip empty files
  if (file.info(fns_R1[i])$size == 0) {
    warning(paste("Empty file skipped:", fns_R1[i]))
    next
  }
  # use the Biostrings function fastq.geometry
  geom <- fastq.geometry(fns_R1[i])
  # extract number of sequences and file name
  df_one_row <- data.frame(n_seq = geom[1], file_name = basename(fns_R1[i]))
  df <- bind_rows(df, df_one_row)
}

View(df)

# Optional: write the table to the working directory:
#write.table(df, file = 'DADA2_output/n_seq.txt', sep='\t', row.names = FALSE, na='',quote=FALSE)

# plot the histogram with number of sequences:
ggplot(df, aes(x = n_seq)) +
  geom_histogram(alpha = 0.5, position = "identity", binwidth = 15000)
hist(df$n_seq, breaks = 10)

# plot boxplot with number of sequences: 
ggplot(df, aes(x = "", y = n_seq)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.5) 

#### Plot Quality for each fastq file #### done once so dont do it again ##
for (i in seq_along(fns)) {
  
  # Plot quality profile
  p1 <- plotQualityProfile(fns[i])
  
  # Show on screen only for first two samples (optional)
  if (i <= 2) {
    print(p1)
  }
  
  # Define output file name
  p1_file <- file.path(qual_dir, paste0(basename(fns[i]), ".qual.pdf"))
  
  # Save as PDF
  ggsave(
    filename = p1_file,
    plot = p1,
    device = "pdf",
    width = 15,
    height = 15,
    units = "cm"
  )
}


#### DADA2 ####
# Running filtering steps to clean data
# Each parameter below controls how strict the filtering is.
# Changing them will alter how many reads survive vs. how clean they are.
# The goal is to find a good balance between read quality and retention.

out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, 
                     maxN = 0, ## Maximum number of ambiguous bases ("N") allowed per read.
                     # DADA2 requires 0 to function properly.
                     
                     maxEE = c(2, 2), # Maximum "expected errors" (EE) allowed per read: forward, reverse.
                     # A better quality metric than average Q-scores.
                     # Typical range: c(2,2) (high quality) to c(2,5) (lower-quality reverse reads) .
                     # - Smaller numbers = stricter (higher quality, fewer reads kept)
                     # - Larger numbers = more lenient (lower quality, more reads kept)
                     
                     truncQ = 2,  #Truncate a read at the first base where quality score ≤ this value.
                     # Removes poor-quality tails dynamically.
                     # Default = 2 (very lenient, trims at Q2).
                     # You can make this stricter (e.g., truncQ = 10 or 15) to remove bad tails faster,
                     # but beware: overly short reads may fail to merge later.
                     
                     minLen = 50, # Minimum read length to keep (after trimming).
                     # Reads shorter than this are discarded.
                     # For ITS or variable-length markers, set a lower threshold (e.g. 50).
                     # For fixed-length amplicons like 16S V4, truncLen usually controls this anyway.
                     
                     rm.phix = TRUE, # Remove reads matching the PhiX genome (control spike-in used by Illumina). Keep TRUE
                     compress = FALSE, multithread = FALSE)

#### STEP 2. Learn Errors
# # Uses the quality scores in your filtered FASTQs to estimate per-cycle error rates.

# The defualt error funciton assumes Illumina data (HiSeq, MiSeq, NextSeq But *NOT* NovaSeq). 

err_R1 <- learnErrors(filt_R1, multithread = FALSE, errorEstimationFunction = loessErrfun)
plotErrors(err_R1, nominalQ = TRUE)
ggsave("DADA2_output/error_model_FWD.pdf")
#saveRDS(err_R1, "DADA2_output/err_R1.rds")

err_R2 <- learnErrors(filt_R2, multithread = FALSE, errorEstimationFunction = loessErrfun)
plotErrors(err_R2, nominalQ = TRUE)
ggsave("DADA2_output/error_model_REV.pdf")
#saveRDS(err_R2, "DADA2_output/err_R2.rds")


#### STEP 3. Dereplicate the reads 
#### this will group identical sequences together 
derep_R1 <- derepFastq(filt_R1, verbose = TRUE)
derep_R2 <- derepFastq(filt_R2, verbose = TRUE)


# Name the derep-class objects by the sample names
names(derep_R1) <- sample.names
names(derep_R2) <- sample.names


#### STEP 4. DADA2 main inference
# Goal: Use the learned error models (err_R1/err_R2) to infer exact ASVs per sample.
# What it does:
#   - Models sequencing errors explicitly and "denoises" dereplicated reads
#   - Distinguishes true variants from errors

#Key parameters:
#   pool         = FALSE / PSEUDO / TRUE
#                 FALSE (default)       → fastest; per-sample inference
#                 PSEUDO              → good sensitivity with modest cost 
#                 TRUE                  → most sensitive, slowest, needs RAM
#   multithread  = TRUE to speed up on multi-core machines (set FALSE on Windows if it misbehaves)

dada_R1 <- dada(derep_R1, err = err_R1, multithread = TRUE, pool = FALSE)
dada_R2 <- dada(derep_R2, err = err_R2, multithread = TRUE, pool = FALSE)

# Quick peek to confirm denoising ran and how many variants were inferred in a sample: 
dada_R1[[3]]
dada_R2[[1]]

#### STEP 5. Merge paired reads ################################################
# Goal: Align denoised forward and reverse reads; only keep pairs that overlap
#       exactly (within rules) to reconstruct the full amplicon.
#
# Why it matters:
#   - If reads don’t sufficiently overlap after trimming, merge rate drops.
#   - For 2x250 V4 16S, overlaps are big; for longer/variable regions (e.g., ITS),
#     be careful with truncation and overlap.
#
# Key parameters in mergePairs():
#   minOverlap       = 12 (default). Increase if you want stricter merges; decrease if overlap is tight.
#   maxMismatch      = 0 (default). Set to >0 to allow mismatches in the overlap (rarely needed).
#   trimOverhang     = TRUE trims hanging overhangs (useful when primers not fully removed).
#   justConcatenate  = TRUE allows non-overlapping merges (NOT recommended unless you know what you’re doing).
#
# Troubleshooting low merge rates:
#   - Ensure forward/reverse truncLen choices still leave ≥20–30 bp overlap.
#   - Relax truncation (keep reads longer) or reduce minOverlap slightly.


mergers <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, minOverlap = 12, maxMismatch = 0,  verbose = TRUE, )

# Inspect one sample’s merged table
head(mergers[[2]])

#### STEP 6. Build the ASV table & basic length check ##########################
# Goal: Create a samples × ASVs count matrix from merged reads.
#       Optionally inspect/limit sequence lengths to expected ranges.

seqtab <- makeSequenceTable(mergers) # rows = samples, cols = unique ASV sequences
dim(seqtab)  # quick size check (n_samples × n_ASVs)

# Length diagnostics: helpful to spot off-targets or residual junk.
len_vec <- nchar(getSequences(seqtab))
len_tab <- table(len_vec)

#Check sequence lenght of ASVs
print(len_tab)



#### STEP 7. Remove chimeras ####
# Goal: Remove sequences that are recombinations of more abundant “parent” ASVs.
# Why: Chimeras are PCR artifacts; keeping them inflates richness and biases results.
#
# Key knobs:
#   method = "consensus"   → robust across samples (good default)
#            "pooled"      → slightly more sensitive; can be stricter
#   multithread            → TRUE to speed up if your machine allows

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)

## Save the clean sequence table to a file ### 
saveRDS(seqtab.nochim, "DADA2_output/seqtab.nochim.rds")

# Get some stats:
# Compute % of non chimeras
paste0("% of non chimeras : ", sum(seqtab.nochim)/sum(seqtab) * 100)
paste0("total number of sequences : ", sum(seqtab.nochim))

#### STEP 8. Track reads across steps #########################################
# Goal: Summarize how many reads survived each step, per sample.
# Columns:
#   input → filtered → denoised → merged → tabled → nonchim
# Use this to spot bottlenecks (e.g., low merging from too-aggressive truncation).


getN <- function(x) sum(getUniques(x))  # helper: sum of unique reads

track <- cbind(
  out,                                   # input, filtered
  sapply(dada_R1, getN),                 # denoised (forward)
  sapply(mergers, getN),                 # merged
  rowSums(seqtab),                       # tabled (merged counts)
  rowSums(seqtab.nochim)                 # nonchim (final counts)
)

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names

print(track)
# view(track)






#### Transforming and saving the ASV sequences
seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% rownames_to_column(var = "sequence") %>%
  rowid_to_column(var = "ASVNumber") %>% mutate(ASVNumber = sprintf("ASV_%05d",
                                                                    ASVNumber)) %>% mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))

#### Optional: Extract the sequences and export them in a fasta file:
df <- seqtab.nochim_trans
seq_out <- Biostrings::DNAStringSet(df$sequence)
names(seq_out) <- df$ASVNumber
seq_out

Biostrings::writeXStringSet(seq_out, str_c(dada2_dir, "ASV_no_taxonomy.fasta"),
                            compress = FALSE, width = 20000)
