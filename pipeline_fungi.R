# Make sure you have access to a unix system
# install conda (or miniconda etc)
# create an environment in conda, call it a name and install:

# cutadapt
conda install -c bioconda cutadapt
# sabre
conda install -c bioconda sabre
# dos2unix
conda install -c conda-forge dos2unix
# R
conda install -c r r

#create data folder in my data directory
mkdir data
cd data 
mikdir seqfiles

# Open the terminal and navigate to your sequencing raw_data folder. For me it is here 
cd /data/bigexpansion/michadm/seqdata/2024-06-10_Novogene_NovaSeq_250PE_Pernettya_fungi_metabarcoding_DanielAnguloSerrano/X204SC24050751-Z01-F001/01.RawData/
cd /data/bigexpansion/michadm/seqdata/2024-06-20_Novogene_NovaSeq_250PE_Pernettya_fungi_metabarcoding_DanielAnguloSerrano_batch2/X204SC24050751-Z01-F002/01.RawData/
  
  
  # Copy raw data folder to my data directory 
  cp -r /data/bigexpansion/michadm/seqdata/2024-06-10_Novogene_NovaSeq_250PE_Pernettya_fungi_metabarcoding_DanielAnguloSerrano/X204SC24050751-Z01-F001/01.RawData /data/lastexpansion/danieang/data/
  
  #Create folders for each plate and split samples by plate
  
  cd /data/lastexpansion/danieang/data/rawdata/
  mkdir P1 P2 P3 P4
cp -r *_P1 /data/lastexpansion/danieang/data/newrawdata01/P1/
  cp -r *_P2 /data/lastexpansion/danieang/data/newrawdata01/P2/
  cp -r *_P3 /data/lastexpansion/danieang/data/newrawdata01/P3/
  cp -r *_P4 /data/lastexpansion/danieang/data/newrawdata01/P4/
  
  #Since I received the raw data in two separated batches, need to concatenate the files
  #Concatenating the two forward file togheter and the two reverse files togheter -> this script is for P1
  
  cd /data/lastexpansion/danieang/data/newrawdata01/P1/
  mkdir seqfiles
cp **/*.fq.gz /data/lastexpansion/danieang/data/newrawdata01/P1/seqfiles


cd /data/lastexpansion/danieang/data/newrawdata01/P1/seqfiles

# Get a list of unique sample prefixes (up to the FKDL part of naming) THIS WILL WORK ONLY FOR P1-> CHANGE SCRIPT FOR OTHER PLATES

for sample in $(ls *_1.fq.gz | sed -E 's/(.+_P1)_.*_L1_[12]\.fq\.gz/\1/' | sort | uniq); do

# Concatenate forward reads for each sample
cat ${sample}*_L1_1.fq.gz > ${sample}_concat_1.fq.gz

# Concatenate reverse reads for each sample
cat ${sample}*_L1_2.fq.gz > ${sample}_concat_2.fq.gz

#Verify the concatenation
echo "Forward reads for $sample combined into ${sample}_concat_1.fq.gz"
echo "Reverse reads for $sample combined into ${sample}_concat_2.fq.gz"
done

#Concatenation for P2:
for sample in $(ls *_1.fq.gz | sed -E 's/(.+_P2)_.*_[12]\.fq\.gz/\1/' | sort | uniq); do


cat ${sample}*_1.fq.gz > ${sample}_concat_1.fq.gz
cat ${sample}*_2.fq.gz > ${sample}_concat_2.fq.gz

echo "Forward reads for $sample combined into ${sample}_concat_1.fq.gz"
echo "Reverse reads for $sample combined into ${sample}_concat_2.fq.gz"

done

#Move concatenated files to /data/lastexpansion/danieang/data/newrawdata01/P1/concatP1

mv *_concat_1.fq.gz /data/lastexpansion/danieang/data/newrawdata01/P1/concatP1
mv *_concat_2.fq.gz /data/lastexpansion/danieang/data/newrawdata01/P1/concatP1

#I will start demultiplexing P1


## based on the forward barcode read for plate1 that is designated in another text file in the same directory (plate1_barcode_data.txt)
## make shell script to run sabre 

## if you made the barcode file with a tex editor - use dos2unix
source ~/miniforge3/bin/activate 
conda activate myenv 

dos2unix euka01_repbarcodes.txt

echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /data/lastexpansion/danieang/data/euka01_repbarcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
mv rep1f.fq ${bn}_rep1f.fq;
mv rep1r.fq ${bn}_rep1r.fq;
mv rep2f.fq ${bn}_rep2f.fq;
mv rep2r.fq ${bn}_rep2r.fq;
mv rep3f.fq ${bn}_rep3f.fq;
mv rep3r.fq ${bn}_rep3r.fq;
mv rep4f.fq ${bn}_rep4f.fq;
mv rep4r.fq ${bn}_rep4r.fq;  done' > eukasabre.sh

# run shell script
bash eukasabre.sh


## Using cut adapt to remove primers - Doing replicate by replicate and make reverse primer an exact match to the corresponding reverse barcode that was not removed by sabre
## obviously you will need to adapt these sequences for the specific index used for each replicate of each plate..  rc = reverse compliment

#### base primers (i.e. tag IDs etc) are :
# ITS4 Fwd: TCCTCCGCTTATTGATATGC
# ITS86F Rv: GTGAATCATCGAATCTTTGAA

#Plate 1 rep 1 reverse barcode is AGGAA
# ITS4 (fwd): TCCTCCGCTTATTGATATGC
# ITS86F (rcRv): TTCAAAGATTCGATGATTCACTTCCT
# ITS4 (rcfwd): GCATATCAATAAGCGGAGGA
# ITS86F (Rv): AGGAAGTGAATCATCGAATCTTTGAA

echo 'for i in *rep1f.fq; do bn=${i/rep1f.fq}; cutadapt -a TCCTCCGCTTATTGATATGC...TTCAAAGATTCGATGATTCACTTCCT -A AGGAAGTGAATCATCGAATCTTTGAA...GCATATCAATAAGCGGAGGA --untrimmed-output ${bn}.rep1out1.fq.gz --untrimmed-paired-output ${bn}.rep1out2.fq.gz -o ${bn}.rep1.trim1.fq.gz -p ${bn}.rep1.trim2.fq.gz ${bn}rep1f.fq ${bn}rep1r.fq; done' > ITSrep1.sh
bash ITSrep1.sh


#Plate 1 rep 2 reverse barcode is GAGTGG
# ITS4 (fwd): TCCTCCGCTTATTGATATGC
# ITS86F (rcRv): TTCAAAGATTCGATGATTCACCCACTC
# ITS4 (rcfwd): GCATATCAATAAGCGGAGGA
# ITS86F (rv): GAGTGGGTGAATCATCGAATCTTTGAA

echo 'for i in *rep2f.fq; do bn=${i/rep2f.fq}; cutadapt -a TCCTCCGCTTATTGATATGC...TTCAAAGATTCGATGATTCACCCACTC -A GAGTGGGTGAATCATCGAATCTTTGAA...GCATATCAATAAGCGGAGGA --untrimmed-output ${bn}.rep2out1.fq.gz --untrimmed-paired-output ${bn}.rep2out2.fq.gz -o ${bn}.rep2.trim1.fq.gz -p ${bn}.rep2.trim2.fq.gz ${bn}rep2f.fq ${bn}rep2r.fq; done' > ITSrep2.sh
bash ITSrep2.sh

#Plate 1 rep 3 reverse barcode is CCACGTC
# ITS4 (fwd): TCCTCCGCTTATTGATATGC
# ITS86F (rcRv): TTCAAAGATTCGATGATTCACGACGTGG
# ITS4 (rcfwd): GCATATCAATAAGCGGAGGA
# ITS86F (Rv): CCACGTCGTGAATCATCGAATCTTTGAA

echo 'for i in *rep3f.fq; do bn=${i/rep3f.fq}; cutadapt -a TCCTCCGCTTATTGATATGC...TTCAAAGATTCGATGATTCACGACGTGG -A CCACGTCGTGAATCATCGAATCTTTGAA...GCATATCAATAAGCGGAGGA --untrimmed-output ${bn}.rep3out1.fq.gz --untrimmed-paired-output ${bn}.rep3out2.fq.gz -o ${bn}.rep3.trim1.fq.gz -p ${bn}.rep3.trim2.fq.gz ${bn}rep3f.fq ${bn}rep3r.fq; done' > ITSrep3.sh
bash ITSrep3.sh

#Plate 1 rep 4 reverse barcode is TTCTCAGC
# ITS4 (fwd): TCCTCCGCTTATTGATATGC
# ITS86F (rcRv): TTCAAAGATTCGATGATTCACGCTGAGAA
# ITS4 (rcfwd): GCATATCAATAAGCGGAGGA
# ITS86F (Rv): TTCTCAGCGTGAATCATCGAATCTTTGAA

echo 'for i in *rep4f.fq; do bn=${i/rep4f.fq}; cutadapt -a TCCTCCGCTTATTGATATGC...TTCAAAGATTCGATGATTCACGCTGAGAA -A TTCTCAGCGTGAATCATCGAATCTTTGAA...GCATATCAATAAGCGGAGGA --untrimmed-output ${bn}.rep4out1.fq.gz --untrimmed-paired-output ${bn}.rep4out2.fq.gz -o ${bn}.rep4.trim1.fq.gz -p ${bn}.rep4.trim2.fq.gz ${bn}rep4f.fq ${bn}rep4r.fq; done' > ITSrep4.sh
bash ITSrep4.sh

## moving trimmed files to their own directory for dada2 analysis
## make a new folder called p1trimmed
mkdir P1trimmed
## move all your demultiplexed, oligo-trimmed files there.
mv *trim1.fq.gz /data/lastexpansion/danieang/data/trimmed/P1trimmed/
mv *trim2.fq.gz /data/lastexpansion/danieang/data/trimmed/P1trimmed/

  
  #Open R inside conda env
  
  conda activate my_r_env

R

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

## Making filepath based on where the trimmed files are. Now is a good point to follow the DADA2 tutorial while trying these steps
setwd("/data/lastexpansion/danieang/data/trimmed/P1trimmed/")
path <- "/data/lastexpansion/danieang/data/trimmed/P1trimmed/"
list.files(path)
fnFs <- sort(list.files(path, pattern=".trim1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".trim2.fq.gz", full.names = TRUE))
primerHits), REV.ForwardReads = sapply(REV.orients, primerHits), REV.ReverseReads = sapply(REV.orients, primerHits))

###tidying up sample names and replicates
sample.names <-sapply(strsplit(basename(fnFs), "_"), function(x){paste(x[[1]], x[[2]], x[[3]], x[[4]], sep="_")})
reps <- rep(c("r1", "r2", "r3", "r4"))

#Some samples still contain "concat" in their sample.names. doing this to remove "concat" from sample.names:
sample.names <- gsub("concat", "", sample.names)

sample.names <-  paste0(sample.names, reps)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## checking some quality plots 
plotQualityProfile(fnFs[1:20])  
plotQualityProfile(fnRs[1:20]) 

#Quality plots seems to be okay, no quality issues on the tails

## Filtering with standard parameters (did not truncate reads as they are already very short - no quality issues on the tails)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, # truncLen=c(45,45),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out, n=150)

## Learning error rates on these data
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## error - some files did not pass trim and need to be removed from lists - which files??
df.fe <-  data.frame(theref = file.exists(filtFs), therer = file.exists(filtRs), filef = filtRs, filer = filtRs)
subset(df.fe, theref == "FALSE") 
subset(df.fe, therer == "FALSE")

#For Plate 4 there are two files that did not pass trim. Rest of plates had 0 files that did not pass trim.

#Description of process for removing files from subset in case it is needed. Example from Euka script:
//
  ## IF there are any files flagged from this subset, they can be removed by row number like so:
  ## In EUKA example, in line 224 there were no sequences remaining after filtering
  ##redefining 'out' matrix so these files are not included
  #out <- out[file.exists(filtFs),]
  
  ### dropping samples that are empty - these numbers will change plate by plate - BE CAREFUL THESE ARE FOR PLATE 1 ONLY
  ### to tell which numbers to add, look at the output from the above subset commands.
  
  #FOR P4 ONLY
  filtFs <- filtFs[-c(138, 140)]
filtRs <- filtRs[-c(138, 140)]

sample.names <- sample.names[-c(138, 140)]

## Learning error rates on these data (again)

#errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
#errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)



## Having a look at error plots
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

#Seems like there is some issue with this error function in my ITS data. Applying same solution as Euka script: 
#### There is a problem with these error functions with NOVASEQ data. Thus we need to hack the function fitter
## Ok we have an issue here with the error models due to NOvaSeq's binning of quality scores. Implementing the solution of 
## JacobRPrice, Here: https://github.com/benjjneb/dada2/issues/1307
## by modifying the error model, his trial 1: alter loess arguments (weights and span) & enforce monotonicity
library(magrittr)
library(dplyr)

loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

errF <- learnErrors(
  filtFs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

errR <- learnErrors(
  filtRs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

## Having a look at new error plots with altered learn errors function
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

#Looks like the new error plots are aligning better with the black points (actual error rates) and there is consistent decrease in error rates with increasing quality scores. 

## dereplicating our sequence data set to remove redundency before applying dada denoising algorithm
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names - use sample.names1 if you needed to remove samples after filtering
names(derepFs) <- sample.names
names(derepRs) <- sample.names
### applying dada2 core inference algorithm - using default of all libraries processed seperately - no pooling
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


###  also applying core inference algorith - but using pooling (all samples pooled together for sample inference

dadaFPPs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRPPs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)

# For P1: 332 samples were pooled: 15699511 reads in 1263234 unique sequences.

###  also applying core inference algorith - but using pseudo pooling (samples processed independent after sharing info
### between samples)

dadaFpsPPs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRpsPPs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")


## Merging the paired-ends - no pooling
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Merging the paired-ends - true pooling
mergersPP <- mergePairs(dadaFPPs, filtFs, dadaRPPs, filtRs, verbose=TRUE)

### Merging the paired-ends - pseudo-pooling
mergers_psPP <- mergePairs(dadaFpsPPs, filtFs, dadaRpsPPs, filtRs, verbose=TRUE)

## Making sequence table
P1seqtab <- makeSequenceTable(mergers)
P1seqtabPP <- makeSequenceTable(mergersPP)
P1seqtabpsPP <- makeSequenceTable(mergers_psPP)

dim(Pseqtab)
dim(P1seqtabPP)
dim(P1seqtabpsPP)


table(nchar(getSequences(P1seqtab)))
table(nchar(getSequences(P1seqtabPP)))
table(nchar(getSequences(P1seqtabpsPP)))

######## ### 

#Moving the sequence table of each plate into new working directory /data/lastexpansion/danieang/data/trimmed/mergedPlates/


saveRDS(P1seqtab, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P1seqtab.rds")
saveRDS(P1seqtabPP, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P1seqtabPP.rds")
saveRDS(P1seqtabpsPP, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P1seqtabpsPP.rds")

saveRDS(P2seqtab, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P2seqtab.rds")
saveRDS(P2seqtabPP, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P2seqtabPP.rds")
saveRDS(P2seqtabpsPP, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P2seqtabpsPP.rds")

saveRDS(P3seqtab, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P3seqtab.rds")
saveRDS(P3seqtabPP, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P3seqtabPP.rds")
saveRDS(P3seqtabpsPP, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P3seqtabpsPP.rds")

saveRDS(P4seqtab, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P4seqtab.rds")
saveRDS(P4seqtabPP, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P4seqtabPP.rds")
saveRDS(P4seqtabpsPP, "/data/lastexpansion/danieang/data/trimmed/mergedPlates/P4seqtabpsPP.rds")



------------------------------------------------------------------------------------------

library(stringr)
library(abind)
library(tidyverse)
library(dada2) 
library(devtools)
install_github("tobiasgf/lulu")
library(lulu)
library(dplyr)
library(data.table)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library(phyloseq)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(ggplot2)

##loading objects in new working directory
setwd("/data/lastexpansion/danieang/data/trimmed/mergedPlates/")

P1seqtab <- readRDS("P1seqtab.rds")
P1seqtabPP <- readRDS("P1seqtabPP.rds")
P1seqtabpsPP <- readRDS("P1seqtabpsPP.rds")

P2seqtab <- readRDS("P2seqtab.rds")
P2seqtabPP <- readRDS("P2seqtabPP.rds")
P2seqtabpsPP <- readRDS("P2seqtabpsPP.rds")

P3seqtab <- readRDS("P3seqtab.rds")
P3seqtabPP <- readRDS("P3seqtabPP.rds")
P3seqtabpsPP <- readRDS("P3seqtabpsPP.rds")

P4seqtab <- readRDS("P4seqtab.rds")
P4seqtabPP <- readRDS("P4seqtabPP.rds")
P4seqtabpsPP <- readRDS("P4seqtabpsPP.rds")

# Merge sequence tables 
merged_seqtab <- mergeSequenceTables(P1seqtab, P2seqtab, P3seqtab, P4seqtab)
merged_seqtabPP <- mergeSequenceTables(P1seqtabPP, P2seqtabPP, P3seqtabPP, P4seqtabPP)
merged_seqtabpsPP <- mergeSequenceTables(P1seqtabpsPP, P2seqtabpsPP, P3seqtabpsPP, P4seqtabpsPP)

# Remove chimeras from merged sequence tables
merged_seqtab.nochim <- removeBimeraDenovo(merged_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(merged_seqtab.nochim)

merged_seqtabPP.nochim <- removeBimeraDenovo(merged_seqtabPP, method="consensus", multithread=TRUE, verbose=TRUE)
dim(merged_seqtabPP.nochim)

merged_seqtabpsPP.nochim <- removeBimeraDenovo(merged_seqtabpsPP, method="consensus", multithread=TRUE, verbose=TRUE)
dim(merged_seqtabpsPP.nochim)


## Moving to LULU curation

## getting sequences

uniquesToFasta(merged_seqtab.nochim, fout="/data/lastexpansion/danieang/data/trimmed/mergedPlates/merged_seqtab.nochim.fasta", ids=paste0("OTU", seq(length(getSequences(merged_seqtab.nochim)))))
uniquesToFasta(merged_seqtabPP.nochim, fout="/data/lastexpansion/danieang/data/trimmed/mergedPlates/merged_seqtabPP.nochim.fasta", ids=paste0("OTU", seq(length(getSequences(merged_seqtabPP.nochim)))))
uniquesToFasta(merged_seqtabpsPP.nochim, fout="/data/lastexpansion/danieang/data/trimmed/mergedPlates/merged_seqtabpsPP.nochim.fasta", ids=paste0("OTU", seq(length(getSequences(merged_seqtabpsPP.nochim)))))

## Make LULU OTU tables (OTUs: rows, samples: columns)
npool.lulu <- merged_seqtab.nochim
colnames(npool.lulu) <- paste0("OTU", seq(length(getSequences(merged_seqtab.nochim))))
npool.lulu <- t(npool.lulu)

pool.lulu <- merged_seqtabPP.nochim
colnames(pool.lulu) <- paste0("OTU", seq(length(getSequences(merged_seqtabPP.nochim))))
pool.lulu <- t(pool.lulu)

pspool.lulu <- merged_seqtabpsPP.nochim
colnames(pspool.lulu) <- paste0("OTU", seq(length(getSequences(merged_seqtabpsPP.nochim))))
pspool.lulu <- t(pspool.lulu)

#Save workspace 
save.image(file = "ITS.RData")

# Summarize sequence lengths
sequence_lengths1 <- nchar(getSequences(merged_seqtab.nochim))
summary(sequence_lengths1)

sequence_lengths2 <- nchar(getSequences(merged_seqtabPP.nochim))
summary(sequence_lengths2)


########### THIS NEXT PART IN THE BASH TERMINAL WITH BLAST INSTALLED IN PATH (or conda environment)
################### 

#First produce a blast databases with the OTUs
makeblastdb -in merged_seqtab.nochim.fasta -parse_seqids -dbtype nucl
makeblastdb -in merged_seqtabPP.nochim.fasta -parse_seqids -dbtype nucl
makeblastdb -in merged_seqtabpsPP.nochim.fasta -parse_seqids -dbtype nucl

# Then blast the OTUs against the database
blastn -db merged_seqtab.nochim.fasta -outfmt '6 qseqid sseqid pident' -out NoPool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query merged_seqtab.nochim.fasta
blastn -db merged_seqtabPP.nochim.fasta -outfmt '6 qseqid sseqid pident' -out Pool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query merged_seqtabPP.nochim.fasta
blastn -db merged_seqtabpsPP.nochim.fasta -outfmt '6 qseqid sseqid pident' -out psPool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query merged_seqtabpsPP.nochim.fasta

### Now running LULU algorithm in R 
setwd("/data/lastexpansion/danieang/data/trimmed/mergedPlates/")
NoPool_match_list.txt <- read.table("NoPool_match_list.txt")
str(NoPool_match_list.txt)
str(npool.lulu)
Pool_match_list.txt <- read.table("Pool_match_list.txt")
psPool_match_list.txt <- read.table("psPool_match_list.txt")
nopool.nochim.curated_result <- lulu(as.data.frame(npool.lulu), NoPool_match_list.txt)
pool.nochim.curated_result <- lulu(as.data.frame(pool.lulu), Pool_match_list.txt)
pspool.nochim.curated_result <- lulu(as.data.frame(pspool.lulu), psPool_match_list.txt)

## Check out how many OTUs were collapsed:
print(paste0("Not Pooled: ", "OTUs after Lulu: ", nopool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(nopool.nochim.curated_result$original_table)))

print(paste0("Pooled: ", "OTUs after Lulu: ", pool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(pool.nochim.curated_result$original_table)))

print(paste0("PsudoPooled: ", "OTUs after Lulu: ", pspool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(pspool.nochim.curated_result$original_table)))


## Making sequence tables compatible with summary routines below...
## Making vector of row numbers of OTUs kept after curation - correspond to column numbers from the precurated data
rownames(nopool.nochim.curated_result$curated_table)

nopool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(nopool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))
pool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(pool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))
pspool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(pspool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))

nopool.lulu <- t(nopool.nochim.curated_result$curated_table)
colnames(nopool.lulu) <- colnames(merged_seqtab.nochim[, nopool.kept.otus])
pool.lulu <- t(pool.nochim.curated_result$curated_table)
colnames(pool.lulu) <- colnames(merged_seqtabPP.nochim[, pool.kept.otus])
pspool.lulu <- t(pspool.nochim.curated_result$curated_table)
colnames(pspool.lulu) <- colnames(merged_seqtabpsPP.nochim[, pspool.kept.otus])

## now we have LULU curated dataframes. We will subtract the max # sequences in No Template Controls for each OTU from sample OTU

## Tidying up names (can be modified as needed to tidy any rownames)
name.change <- function(x) {
  rownames(x) <- gsub("_P4_r", "_P4r", rownames(x), ignore.case = FALSE, perl = TRUE)
  rownames(x) <- gsub("_P3_r", "_P3r", rownames(x), ignore.case = FALSE, perl = TRUE)
  rownames(x) <- gsub("_P2_r", "_P2r", rownames(x), ignore.case = FALSE, perl = TRUE)
  rownames(x) <- gsub("_P1_r", "_P1r", rownames(x), ignore.case = FALSE, perl = TRUE)
  return(x)
}


# Apply the name change function to the LULU curated data frames
nopool.lulu <- name.change(nopool.lulu)
pool.lulu <- name.change(pool.lulu)
pspool.lulu <- name.change(pspool.lulu)

#Check current rownames
print(rownames(nopool.lulu))
print(rownames(pool.lulu))
print(rownames(pspool.lulu))


library(stringr)

## Some summary info for the data to use later - can be modified as needed according to naming convention

index.info <- function(x){
  y <- data.frame(matrix(NA,    # Create empty data frame
                         nrow = nrow(x),
                         ncol = 0))
  y$rep <- str_sub(rownames(x), start= -1)
  y$sample <- str_sub(rownames(x),1,nchar(rownames(x))-2)
  #y$sample <- ifelse(startsWith(y$sample, "L2")==TRUE, substring(y$sample, 2), y$sample)
  y$full <- rownames(x)
  y$totseq <- rowSums(x)
  y$OTUs <- rowSums(x>0)
  return(y)
}

## getting summary info
nopool.lulu.index <- index.info(nopool.lulu)
pool.lulu.index <- index.info(pool.lulu)
pspool.lulu.index <- index.info(pspool.lulu)

## Stopping here. saving at: /data/lastexpansion/danieang/data/trimmed/mergedPlates/

save.image(file = "merged_ITS.RData")

#/data/lastexpansion/danieang/data/trimmed/mergedPlates/merged_ITS.RDATA

## need to subtract equivalent OTU sequences from NTCs, extraction blanks, field blanks etc.

## subtracting largest NTC value from each OTU from all samples for that OTU
ntc.change <- function(x){
  mind <- apply(x[grep("PCRB", rownames(x)), ], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

##### removing NTC counts from each OTU using the highest sequencing count per OTU in the NTCs

nopool.lulu.ntc <-  ntc.change(nopool.lulu)
pool.lulu.ntc <-  ntc.change(pool.lulu)
pspool.lulu.ntc <-  ntc.change(pspool.lulu)

### Controlling for extraction blanks

## need experimental data
exdata <- read.csv("ITS_metadata_with_habitat.csv")
str(exdata)

## function to convert data to list of matrices split by blank extract batch
#Modify "sample.info.nopool" and "nopool.lulu.index" for each dataset!! 
#I am hardcoring each sample.info in each function, so each function must be used only with the correct object

ntc.to.blankcontrol <-  function (x){
  z <- x
  sample.info.nopool <- dplyr::bind_rows(nopool.lulu.index, .id = "column_label")
  row.names(sample.info.nopool) = sample.info.nopool$full
  z = merge(nopool.lulu.ntc, sample.info.nopool, by =  'row.names', all.x=TRUE)
  z = merge(z, exdata, by.x="sample", by.y="sample", all.x=TRUE)
  rownames(z) <- z$Row.names 
  z = split(z, z$extract_blank)
  dropnames <- colnames(z[[2]][, c(which(nchar(colnames(z[[2]]))< 30))])
  z <- lapply(z, function(x) x[!(names(x) %in% dropnames)])
  z = lapply(z, as.matrix)
  return(z)
}

nopool.lulu.ntc.blankC <- ntc.to.blankcontrol(nopool.lulu.ntc)
pool.lulu.ntc.blankC <- ntc.to.blankcontrol(pool.lulu.ntc)
pspool.lulu.ntc.blankC <- ntc.to.blankcontrol(pspool.lulu.ntc)

## getting rid of experimental controls with no blank samples
nopool.lulu.ntc.blankC1 <- nopool.lulu.ntc.blankC[-1]
pool.lulu.ntc.blankC1 <- pool.lulu.ntc.blankC[-1]
pspool.lulu.ntc.blankC1 <- pspool.lulu.ntc.blankC[-1]

unique(exdata$extract_blank) # copy in blank extract names
blank.change <- function(x){
  mind <- apply(x[grep("E_1_7|E_2_12|E_3_12|E_4_12|E_5_12|E_6_12|E_7_12|E_8_12|E_9_12|E_10_12|E_11_12|E_12_6|E_S1_8|E_S2_6|E_BC|E_W4", rownames(x)), ], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

nopool.lulu.ntc.blankC2 <- lapply(nopool.lulu.ntc.blankC1, blank.change)
pool.lulu.ntc.blankC2 <- lapply(pool.lulu.ntc.blankC1, blank.change)
pspool.lulu.ntc.blankC2 <- lapply(pspool.lulu.ntc.blankC1, blank.change)

## inf values because trying to subtract 0 from 0 .. converting back 

nopool.lulu.ntc.blankC2 <- lapply(nopool.lulu.ntc.blankC2,function(x) replace(x, !is.finite(x), 0))
pool.lulu.ntc.blankC2 <- lapply(pool.lulu.ntc.blankC2, function(x) replace(x, !is.finite(x), 0))
pspool.lulu.ntc.blankC2 <- lapply(pspool.lulu.ntc.blankC2, function(x) replace(x, !is.finite(x), 0))

########## 


##### The same can occur for field blanks - merging dataset again
nopool.lulu.ntc.blankC2 <- lapply(nopool.lulu.ntc.blankC2,function(x) as.data.frame(x))
pool.lulu.ntc.blankC2 <- lapply(pool.lulu.ntc.blankC2,function(x) as.data.frame(x))
pspool.lulu.ntc.blankC2 <- lapply(pspool.lulu.ntc.blankC2,function(x) as.data.frame(x))

nopool.lulu.ntc.blank.fieldb <- dplyr::bind_rows(nopool.lulu.ntc.blankC2, .id = "column_label")
pool.lulu.ntc.blank.fieldb <- dplyr::bind_rows(pool.lulu.ntc.blankC2, .id = "column_label")
pspool.lulu.ntc.blank.fieldb <- dplyr::bind_rows(pspool.lulu.ntc.blankC2, .id = "column_label")


## function to group samples into lists via field blanks
#Modify "sample.info.nopool" and "nopool.lulu.index" for each dataset!!

exblank.to.fieldblankcontrol <-  function (x){
  z <- x
  sample.info.nopool <- dplyr::bind_rows(nopool.lulu.index, .id = "column_label")
  row.names(sample.info.nopool) = sample.info.nopool$full
  z = merge(nopool.lulu.ntc, sample.info.nopool , by =  'row.names', all.x=TRUE)
  z = merge(z, exdata, by.x="sample", by.y="sample", all.x=TRUE)
  rownames(z) <- z$Row.names 
  z = split(z, z$Fieldcontrol)
  dropnames <- colnames(z[[2]][, c(which(nchar(colnames(z[[2]]))< 30))])
  z <- lapply(z, function(x) x[!(names(x) %in% dropnames)])
  z = lapply(z, as.matrix)
  return(z)
}


nopool.lulu.ntc.blank.fieldblist <- exblank.to.fieldblankcontrol(nopool.lulu.ntc.blank.fieldb)
pool.lulu.ntc.blank.fieldblist <- exblank.to.fieldblankcontrol(pool.lulu.ntc.blank.fieldb)
pspool.lulu.ntc.blank.fieldblist <- exblank.to.fieldblankcontrol(pspool.lulu.ntc.blank.fieldb)


## getting rid non-samples experimental controls with no blank samples ->  eliminate blank extraction controls

nopool.lulu.ntc.blank.fieldblist <- nopool.lulu.ntc.blank.fieldblist[-1]
pool.lulu.ntc.blank.fieldblist <- pool.lulu.ntc.blank.fieldblist[-1]
pspool.lulu.ntc.blank.fieldblist <- pspool.lulu.ntc.blank.fieldblist[-1]

unique(exdata$Fieldcontrol) # copy in blank extract names

#Applying this function to remove OTUs (OTU count) present in soil controls by sampling location. This function is only applied for the root-only, dataset. 

unique(exdata$Fieldcontrol) # copy in blank extract names

fb.blank.change <- function(x){
  mind <- apply(x[grep("S_S1_1|S_S1_2|S_S1_3|S_S1_4|S_S1_5|S_S1_6|S_S1_7|S_S2_1|S_S2_2|S_S2_3|S_S2_4|S_S2_5", rownames(x)), ], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

nopool.lulu.ntc.blank.fieldblist1 <- lapply(nopool.lulu.ntc.blank.fieldblist, fb.blank.change)
pool.lulu.ntc.blank.fieldblist1 <- lapply(pool.lulu.ntc.blank.fieldblist, fb.blank.change)
pspool.lulu.ntc.blank.fieldblist1 <- lapply(pspool.lulu.ntc.blank.fieldblist, fb.blank.change)

nopool.lulu.ntc.blank.fieldblist1 <- lapply(nopool.lulu.ntc.blank.fieldblist1,function(x) replace(x, !is.finite(x), 0))
pool.lulu.ntc.blank.fieldblist1 <- lapply(pool.lulu.ntc.blank.fieldblist1, function(x) replace(x, !is.finite(x), 0))
pspool.lulu.ntc.blank.fieldblist1 <- lapply(pspool.lulu.ntc.blank.fieldblist1, function(x) replace(x, !is.finite(x), 0))


## Merging data back to just biological samples
##### The same can occur for field blanks - merging dataset again
nopool.lulu.ntc.blank.fieldblist1 <- lapply(nopool.lulu.ntc.blank.fieldblist1,function(x) as.data.frame(x))
pool.lulu.ntc.blank.fieldblist1 <- lapply(pool.lulu.ntc.blank.fieldblist1,function(x) as.data.frame(x))
pspool.lulu.ntc.blank.fieldblist1 <- lapply(pspool.lulu.ntc.blank.fieldblist1,function(x) as.data.frame(x))


nopool.lulu.ntc.blank.fielddone<- dplyr::bind_rows(nopool.lulu.ntc.blank.fieldblist1, .id = "column_label")

pool.lulu.ntc.blank.fielddone<- dplyr::bind_rows(pool.lulu.ntc.blank.fieldblist1, .id = "column_label")

pspool.lulu.ntc.blank.fielddone<- dplyr::bind_rows(pspool.lulu.ntc.blank.fieldblist1, .id = "column_label")

#removing column id labels
nopool.lulu.ntc.blank.fielddone<- nopool.lulu.ntc.blank.fielddone[,-1]
pool.lulu.ntc.blank.fielddone<- pool.lulu.ntc.blank.fielddone[,-1]
pspool.lulu.ntc.blank.fielddone<-  pspool.lulu.ntc.blank.fielddone[,-1]

## removing soil samples from dataframes and split into lists of individual samples for replicate control and combination etc
#Again, this function is only valid for the root-only (and soil OTU removed) dataset. 
unique(rownames(nopool.lulu.ntc.blank.fielddone))

controlledblanks.to.samplelist <- function(x) {
  minusNTCBLANK.sample <- x
  
  # Using pattern matching (grep) instead of exact matching
  samples_to_remove <- c("S_S1_1", "S_S1_2", "S_S1_3", "S_S1_4", "S_S1_5", 
                         "S_S1_6", "S_S1_7", "S_S2_1", "S_S2_2", "S_S2_3", 
                         "S_S2_4", "S_S2_5")
  
  # Use grep to match patterns correctly and remove rows
  rows_to_remove <- grep(paste(samples_to_remove, collapse = "|"), rownames(minusNTCBLANK.sample))
  minusNTCBLANK.sample <- minusNTCBLANK.sample[-rows_to_remove, , drop = FALSE]
  
  # Hardcoded metadata for nopool
  sample.info.nopool <- dplyr::bind_rows(nopool.lulu.index, .id = "column_label")
  row.names(sample.info.nopool) <- sample.info.nopool$full
  
  # Merge and process
  minusNTCBLANK.sample <- merge(minusNTCBLANK.sample, sample.info.pspool, by = 'row.names', all.x = TRUE)
  rownames(minusNTCBLANK.sample) <- minusNTCBLANK.sample$Row.names
  minusNTCBLANK.sample <- split(minusNTCBLANK.sample, minusNTCBLANK.sample$sample)
  
  # Drop columns with names shorter than 30 characters
  dropnames <- colnames(minusNTCBLANK.sample[[2]][, c(which(nchar(colnames(minusNTCBLANK.sample[[2]])) < 30))])
  minusNTCBLANK.sample <- lapply(minusNTCBLANK.sample, function(x) x[!(names(x) %in% dropnames)])
  
  return(minusNTCBLANK.sample)
}

nopool.lulu.controlled_1 <- controlledblanks.to.samplelist(nopool.lulu.ntc.blank.fielddone)
pool.lulu.controlled_1 <- controlledblanks.to.samplelist(pool.lulu.ntc.blank.fielddone)
pspool.lulu.controlled_1 <-  controlledblanks.to.samplelist(pspool.lulu.ntc.blank.fielddone)

str(nopool.lulu.controlled_1)
head(rownames(nopool.lulu.controlled_1))

View(sample.info.nopool)


## if you view just these raw data, some samples and some replicates did not sequence well  

View()

## so there is the "raw" data that can be combined (no PCR replicate filtering) - no filtering at all required.
## and then also filtering to rg

# Making 4 types of phyloseq objects based on OTUs presence in PCR replicates:
# Raw data
# Raw data without soil sample's sequences -> see script 3_nosoil_phyloseq
# OTU dropped if only in single PCR replicate (i.e. OTU in at least 2/4 replicated) w/soil samples
# OTU dropped if only in two or fewer PCR replicates (i.e. OTU in at least 3/4 replicated) w/soil samples
# Only OTUs in all four PCR replicates w/soil samples

## Raw data including soil sample's sequences 
nopool.lulu.controlled_1
pool.lulu.controlled_1
pspool.lulu.controlled_1


## Raw data without soil sample's sequences -> see script 3_nosoil_phyloseq


#OTU dropped if only in single PCR replicate (i.e. OTU in at least 2/4 replicated)
rep.groups2 <- function(x){
  r2 <- apply(x, 2, function(c) replace(c, sum(c!=0)<2, 0))
  return(r2)
}

rg2.nopool.lulu.controlled <- lapply(nopool.lulu.controlled_1 , rep.groups2)
rg2.pool.lulu.controlled <- lapply(pool.lulu.controlled_1 , rep.groups2)
rg2.pspool.lulu.controlled <- lapply(pspool.lulu.controlled_1 , rep.groups2)

#OTU dropped if only in two or fewer PCR replicates (i.e. OTU in at least 3/4 replicated)
rep.groups3 <- function(x){
  r3 <- apply(x, 2, function(c) replace(c, sum(c!=0)<3, 0))
  return(r3)
}

rg3.nopool.lulu.controlled <- lapply(nopool.lulu.controlled_1 , rep.groups3)
rg3.pool.lulu.controlled <- lapply(pool.lulu.controlled_1 , rep.groups3)
rg3.pspool.lulu.controlled <- lapply(pspool.lulu.controlled_1 , rep.groups3)

#Only OTUs in all four PCR replicates 

rep.groups4 <- function(x){
  r4 <- apply(x, 2, function(c) replace(c, sum(c!=0)<4, 0))
  return(r4)
}


rg4.nopool.lulu.controlled <- lapply(nopool.lulu.controlled_1 , rep.groups4)
rg4.pool.lulu.controlled <- lapply(pool.lulu.controlled_1 , rep.groups4)
rg4.pspool.lulu.controlled <- lapply(pspool.lulu.controlled_1 , rep.groups4)

## now we can extract the sequence data and make a phylopseq object - taxonomy can be added later

### extracting sequence data
library(dada2)
nopool.lulu.controlled_1 <- lapply(nopool.lulu.controlled_1, as.matrix)
pool.lulu.controlled_1 <- lapply(pool.lulu.controlled_1, as.matrix)
pspool.lulu.controlled_1 <- lapply(pspool.lulu.controlled_1, as.matrix)

uniquesToFasta(as.matrix(nopool.lulu.controlled_1[[1]]), fout="/data/lastexpansion/danieang/data/trimmed/mergedPlates/nopool_wSoil.fasta", ids=paste0("OTU", seq(length(getSequences(nopool.lulu.controlled_1[[1]])))))
uniquesToFasta(as.matrix(pool.lulu.controlled_1[[1]]), fout="/data/lastexpansion/danieang/data/trimmed/mergedPlates/pool_wSoil.fasta", ids=paste0("OTU", seq(length(getSequences(pool.lulu.controlled_1[[1]])))))
uniquesToFasta(as.matrix(pspool.lulu.controlled_1[[1]]), fout="/data/lastexpansion/danieang/data/trimmed/mergedPlates/pspool_wSoil.fasta", ids=paste0("OTU", seq(length(getSequences(pspool.lulu.controlled_1[[1]])))))

sequence_lengths3 <- nchar(getSequences("/data/lastexpansion/danieang/data/trimmed/mergedPlates/nopool_wSoil.fasta"))
summary(sequence_lengths3)


library(phyloseq)
library(seqinr)
#reading back in sequence files
nopoolseqs <- read.fasta("nopool_wSoil.fasta")
row.names(nopoolseqs) <- nopoolseqs$id
poolseqs  <- read.fasta("pool_wSoil.fasta")
row.names(poolseqs) <- poolseqs$seq.name
pspoolseqs  <- read.fasta("pspool_wSoil.fasta")
row.names(pspoolseqs) <- pspoolseqs$seq.name

## Function to make phylosequ object from OTU, tax data and sample data

to.one.matrix <- function(x){
  lah <- do.call(rbind.data.frame, x)
  rownames(lah) <- names(x)
  colnames(lah) <- names(x[[1]])
  lah <- as.matrix(lah)
  return(lah)
}

row.names(exdata) <- exdata$sample

make.phylo <- function (x, z, k){ # x = OTU data to this point, z = sample data, k = sequence data
  ## Summing all the replicates to get a single row for each sample
  single <- lapply(x, function(w) colSums(w) )
  ### getting a single matrix from list of samples (with single row)
  test <- to.one.matrix(single)
  ## Giving the columns the OTU name instead of the sequence 
  colnames(test) <-  names(k)
  ### Making phyloseq initial object
  wanted <- phyloseq(otu_table(test, taxa_are_rows = FALSE), sample_data(z))
  ## getting DNA data
  dna <- Biostrings::DNAStringSet(colnames(x[[1]]))
  ## naming sequences as OTUs
  names(dna) <- names(k)
  ##Finishinf phyloseq object
  wanted <- merge_phyloseq(wanted, dna)
  return(wanted)
}



## making phyoseq data from raw data w/ soil samples
nopoolps_wSoil <- make.phylo(nopool.lulu.controlled_1 , exdata, nopoolseqs)
poolps_wSoil <- make.phylo(pool.lulu.controlled_1 , exdata, poolseqs)
pspoolps_wSoil <- make.phylo(pspool.lulu.controlled_1 , exdata, pspoolseqs)


## making phyloseq data from rg2 (OTU in at least 2 replicates data) w/soil
rg2.nopoolps <- make.phylo(rg2.nopool.lulu.controlled , exdata, nopoolseqs)
rg2.poolps <- make.phylo(rg2.pool.lulu.controlled , exdata, poolseqs)
rg2.pspoolps <- make.phylo(rg2.pspool.lulu.controlled , exdata, pspoolseqs)

## making phyloseq data from rg3 (OTU in at least 3 replicates data) w/soil
rg3.nopoolps <- make.phylo(rg3.nopool.lulu.controlled , exdata, nopoolseqs)
rg3.poolps <- make.phylo(rg3.pool.lulu.controlled , exdata, poolseqs)
rg3.pspoolps <- make.phylo(rg3.pspool.lulu.controlled , exdata, pspoolseqs)

## making phyloseq data from rg4 (OTU in at least 4 replicates data) w/soil
rg4.nopoolps <- make.phylo(rg4.nopool.lulu.controlled , exdata, nopoolseqs)
rg4.poolps <- make.phylo(rg4.pool.lulu.controlled , exdata, poolseqs)
rg4.pspoolps <- make.phylo(rg4.pspool.lulu.controlled , exdata, pspoolseqs)




saveRDS(nopoolps_wSoil, "nopool_phyloseq_soil.rds")
saveRDS(poolps_wSoil, "pool_phyloseq_soil.rds")
saveRDS(pspoolps_wSoil, "pspool_phyloseq_soil.rds")

saveRDS(rg2.nopoolps, "rg2.nopoolps.rds")
saveRDS(rg2.poolps, "rg2.poolps.rds")
saveRDS(rg2.pspoolps, "rg2.pspoolps.rds")

saveRDS(rg3.nopoolps, "rg3.nopoolps.rds")
saveRDS(rg3.poolps, "rg3.poolps.rds")
saveRDS(rg3.pspoolps, "rg3.pspoolps.rds")

saveRDS(rg4.nopoolps, "rg4.nopoolps.rds")
saveRDS(rg4.poolps, "rg4.poolps.rds")
saveRDS(rg4.pspoolps, "rg4.pspoolps.rds")


nopoolps_wSoil <- "/data/lastexpansion/danieang/data/trimmed/mergedPlates/nopool_phyloseq_soil.rds"
poolps_wSoil <- "/data/lastexpansion/danieang/data/trimmed/mergedPlates/pool_phyloseq_soil.rds" 
pspoolps_wSoil <- "/data/lastexpansion/danieang/data/trimmed/mergedPlates/pspool_phyloseq_soil.rds"


total_reads <- sum(phyloseq::otu_table(nopoolps_wSoil))
cat("total reads:", total_reads, "\n")


### All data is in phyloseq objects now

## look into microviz for ways to examine your data..

save.image("merged_ITS2.RData")

################

#Moving to taxonomic assignment with DADA2

#Download the general release FASTA file for fungi and point to the path
unite <- "/data/lastexpansion/danieang/data/database/UNITE_small.fasta"

# Path to curated sequences FASTA file
nopool_fasta_wSoil <- "/data/lastexpansion/danieang/data/trimmed/mergedPlates/nopool_wSoil.fasta"
pool_fasta_wSoil <- "/data/lastexpansion/danieang/data/trimmed/mergedPlates/pool_wSoil.fasta"
pspool_fasta_wSoil <- "/data/lastexpansion/danieang/data/trimmed/mergedPlates/pspool_wSoil.fasta"

library(Biostrings)

# Reading sequences from FASTA file
nopool_fasta_wSoil <- readDNAStringSet("/data/lastexpansion/danieang/data/trimmed/mergedPlates/nopool_wSoil.fasta")
pool_fasta_wSoil <- readDNAStringSet("/data/lastexpansion/danieang/data/trimmed/mergedPlates/pool_wSoil.fasta")
pspool_fasta_wSoil <- readDNAStringSet("/data/lastexpansion/danieang/data/trimmed/mergedPlates/pspool_wSoil.fasta")


sequence_lengths <- nchar(getSequences(nopool_fasta_wSoil))
summary(sequence_lengths)

# Assign taxonomy to With/soil dataset
taxa_nopool_wSoil_80 <- assignTaxonomy(nopool_fasta_wSoil, unite, multithread=50, tryRC=TRUE, minBoot=80)
taxa_pool_wSoil_80 <- assignTaxonomy(pool_fasta_wSoil, unite, multithread=50, tryRC=TRUE, minBoot=80)
taxa_pspool_wSoil_80 <- assignTaxonomy(pspool_fasta_wSoil, unite, multithread=50, tryRC=TRUE, minBoot=80)

taxa_nopool_wSoil_60 <- assignTaxonomy(nopool_fasta_wSoil, unite, multithread=50, tryRC=TRUE, minBoot=60)
taxa_pool_wSoil_60 <- assignTaxonomy(pool_fasta_wSoil, unite, multithread=50, tryRC=TRUE, minBoot=60)
taxa_pspool_wSoil_60 <- assignTaxonomy(pspool_fasta_wSoil, unite, multithread=50, tryRC=TRUE, minBoot=60)


save.image(file = "dada2_taxa.RData")
load("dada2_taxa.RData")

# Function to summarize taxonomic assignment
taxa_summary <- function(taxa) {
  sapply(1:ncol(taxa), function(i) sum(!is.na(taxa[, i])))
}

# Summarize assignments for With/soil dataset at minBoot = 80
summary_nopool_wSoil_80 <- taxa_summary(taxa_nopool_wSoil_80)
summary_pool_wSoil_80 <- taxa_summary(taxa_pool_wSoil_80)
summary_pspool_wSoil_80 <- taxa_summary(taxa_pspool_wSoil_80)

# Summarize assignments for With/soil dataset at minBoot = 60
summary_nopool_wSoil_60 <- taxa_summary(taxa_nopool_wSoil_60)
summary_pool_wSoil_60 <- taxa_summary(taxa_pool_wSoil_60)
summary_pspool_wSoil_60 <- taxa_summary(taxa_pspool_wSoil_60)

# Combine results for With/soil dataset into a data frame
summary_wSoil_df <- data.frame(
  Rank = colnames(taxa_nopool_wSoil_80),
  Nopool_80 = summary_nopool_wSoil_80,
  Pool_80 = summary_pool_wSoil_80,
  PsPool_80 = summary_pspool_wSoil_80,
  Nopool_60 = summary_nopool_wSoil_60,
  Pool_60 = summary_pool_wSoil_60,
  PsPool_60 = summary_pspool_wSoil_60
)

print(summary_wSoil_df)

#Now I will proceed to add taxonomic information into phyloseq objects


# Load the phyloseq objects 

nopoolps_wSoil <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/nopool_phyloseq_soil.rds")
poolps_wSoil <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/pool_phyloseq_soil.rds")
pspoolps_wSoil <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/pspool_phyloseq_soil.rds")

rg2.nopoolps <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/rg2.nopoolps.rds")
rg2.poolps <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/rg2.poolps.rds")
rg2.pspoolps <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/rg2.pspoolps.rds")

rg3.nopoolps <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/rg3.nopoolps.rds")
rg3.poolps <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/rg3.poolps.rds")
rg3.pspoolps <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/rg3.pspoolps.rds")

rg4.nopoolps <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/rg4.nopoolps.rds")
rg4.poolps <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/rg4.poolps.rds")
rg4.pspoolps <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/rg4.pspoolps.rds")

#Convert taxonomy table into phyloseq compatible format (matrix)
tax_table_nopoolps <- tax_table(taxa_nopool_wSoil_80)
tax_table_poolps <- tax_table(taxa_pool_wSoil_80)
tax_table_pspoolps <- tax_table(taxa_pspool_wSoil_80)

# Check if the taxa IDs match between OTU table and taxonomy table
all(taxa_names(nopoolps_wSoil) == rownames(tax_table_nopoolps))
all(taxa_names(poolps_wSoil) == rownames(tax_table_poolps))
all(taxa_names(pspoolps_wSoil) == rownames(tax_table_pspoolps))

#Taxa names do not match between phyloseq objects and taxonomy tables... 

# Sanity checks for consistent OTU names
head(taxa_names(nopoolps_wSoil), 2)
head(taxa_names(tax_table_nopoolps), 2)

# Extract the sequences from tax_table_nopoolps
seqs_tax_nopool <- rownames(tax_table_nopoolps)
seqs_tax_pool <- rownames(tax_table_poolps)
seqs_tax_pspool <- rownames(tax_table_pspoolps)

# Create a mapping between sequences and OTU identifiers
seq_to_otu_nopool <- setNames(taxa_names(nopoolps_wSoil), seqs_tax_nopool)
seq_to_otu_pool <- setNames(taxa_names(poolps_wSoil), seqs_tax_pool)
seq_to_otu_pspool <- setNames(taxa_names(pspoolps_wSoil), seqs_tax_pspool)

# Rename the rows of tax_table_nopoolps using this mapping
rownames(tax_table_nopoolps) <- seq_to_otu_nopool[rownames(tax_table_nopoolps)]
rownames(tax_table_poolps) <- seq_to_otu_pool[rownames(tax_table_poolps)]
rownames(tax_table_pspoolps) <- seq_to_otu_pspool[rownames(tax_table_pspoolps)]

# Now checking that the names match
all(rownames(tax_table_nopoolps) == taxa_names(nopoolps_wSoil))
all(rownames(tax_table_poolps) == taxa_names(poolps_wSoil))
all(rownames(tax_table_pspoolps) == taxa_names(pspoolps_wSoil))

# Combining the taxonomy with existing phyloseq object
nopoolps.dada2.soil <- merge_phyloseq(nopoolps_wSoil, tax_table_nopoolps)
poolps.dada2.soil <- merge_phyloseq(poolps_wSoil, tax_table_poolps)
pspoolps.dada2.soil <- merge_phyloseq(pspoolps_wSoil, tax_table_pspoolps)

saveRDS(nopoolps.dada2.soil, "nopoolps.dada2.soil.rds")
saveRDS(poolps.dada2.soil, "poolps.dada2.soil.rds")
saveRDS(pspoolps.dada2.soil, "pspoolps.dada2.soil.rds")

# Check the merged object
print(nopoolps.dada2.S)
print(poolps.dada2.S)
print(pspoolps.dada2.S)

tax_table(nopoolps.dada2.S)

# Checking the number of taxa before and after the update
ntaxa(nopoolps_wSoil) == ntaxa(nopoolps.dada2.S)
ntaxa(poolps_wSoil) == ntaxa(poolps.dada2.S)
ntaxa(pspoolps_wSoil) == ntaxa(pspoolps.dada2.S)

#Integrating taxonomy into PCR replicate filtered phyloseq objects

# Integrate taxonomy tables into phyloseq objects
rg2.nopoolps.Soil <- rg2.nopoolps
rg2.poolps.Soil <- rg2.poolps
rg2.pspoolps.Soil <- rg2.pspoolps

rg3.nopoolps.Soil <- rg3.nopoolps
rg3.poolps.Soil <- rg3.poolps
rg3.pspoolps.Soil <- rg3.pspoolps

rg4.nopoolps.Soil <- rg4.nopoolps
rg4.poolps.Soil <- rg4.poolps
rg4.pspoolps.Soil <- rg4.pspoolps

tax_table(rg2.nopoolps.Soil) <- tax_table_nopoolps
tax_table(rg2.poolps.Soil) <- tax_table_poolps
tax_table(rg2.pspoolps.Soil) <- tax_table_pspoolps

tax_table(rg3.nopoolps.Soil) <- tax_table_nopoolps
tax_table(rg3.poolps.Soil) <- tax_table_poolps
tax_table(rg3.pspoolps.Soil) <- tax_table_pspoolps

tax_table(rg4.nopoolps.Soil) <- tax_table_nopoolps
tax_table(rg4.poolps.Soil) <- tax_table_poolps
tax_table(rg4.pspoolps.Soil) <- tax_table_pspoolps


tax_table_summary <- tax_table(rg2.nopoolps.Soil)
head(tax_table_summary)


