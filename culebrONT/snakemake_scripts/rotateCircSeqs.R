
## Read Fasta file of a circular sequence, rotate it,
## and save the resulting sequences in a new fasta file.

#### Loading required packages #####
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(Biostrings))

#### COMMAND LINE ARGUMENTS PARSING ######
option_list <- list(
  make_option(c("-f", "--seqFile"),
              type = "character",
              default = NULL,
              help = "Path of fasta file of seqences to be rotated."),
  make_option(c("-o", "--outFilePath"),
              type = "character",
              default = NULL,
              help = "Path where output sequence file will be written."),
  make_option(c("-d", "--circlatorLog"),
              type = "character",
              default = NULL,
              help = "Path to a circlator '04.merge.circularise.log' file to be use to rotate only circularized sequences.")
)

##### RETRIEVEING PARAMS from optparse #####
myArgs <- parse_args(
  OptionParser(usage = "%prog [options]", option_list = option_list,
               description = "Rotate sequences in an input fasta file. Sequences will be rotated if the sequence title contains a 'circular' flag, they will be rotated and their title will be appended with a 'rotated' suffix. If a circlator.log file is provided and no info is available in the titles, the sequences described as circular in the log file will be rotated and their title will be appended with a 'circular|rotated' suffix. Sequences that do not meet these criteria are left untouched.")
)

###########################################
##### FOR INTERACTIVE ARGUMENTS PASSING TESTS
# testArgs  <-  c(
#   "--seqFile=/home/cunnac/Lab-Related/MyScripts/pacbacpipe/tests/tagCirSeq/talCmut_canu/06.fixstartcircFlagged.fasta",
#   "--outFilePath=/home/cunnac/rotateSeqTest.fasta",
#   "--circlatorLog=/home/cunnac/Lab-Related/MyScripts/pacbacpipe/tests/tagCirSeq/talCmut_canu/04.merge.circularise.log"
# )
# myArgs <- parse_args(OptionParser(usage = "%prog", option_list = option_list), args = testArgs)
###########################################

# Assign parameter values
seqFile <- myArgs$seqFile
outSeqFile <- myArgs$outFilePath
if(is.null(outSeqFile)) outSeqFile <- paste0(gsub("(^.*)\\..*$", "\\1", seqFile), "_rotated.fasta")
circlatorLog <- myArgs$circlatorLog
circLogProvided <- ifelse(is.null(circlatorLog) || circlatorLog == "", FALSE, TRUE)

# This is the number by which seq lenght is devided to determine length of offset
n <- 5

# Loading sequences and determining if any is flagged as circular
seq <- readDNAStringSet(seqFile) #seq <- DNAStringSet(c("AAATTTGGGCCCNNN", "AAAATTTTCCCC"))
circularInName <- grepl("circular", names(seq), ignore.case = TRUE)

if (circLogProvided && any(circularInName)) { # if log and flagged seq, sanity check on names
  circlatorDf <- read.delim(circlatorLog, stringsAsFactors = FALSE)
  seqNamePrefixes <- sapply(names(seq), function(x) {strsplit(x, split = "[_ ]")[[1]][1]})
  circNamePrefixes <- sapply(circlatorDf$X.Contig, function(x) {strsplit(x, split = "[_ ]")[[1]][1]})
  if(!identical(length(setdiff(seqNamePrefixes, circNamePrefixes)), 0L)) {
    stop("Sequence names in sequence file and circlatorLog do not match:\n",
         "Names in fasta file: ", names(seqNamePrefixes), "\n",
         "Names in circlator log: ", names(circNamePrefixes), "\n")
  } else {
    cat("## Provided a circlator log file and sequences in fasta file are also otherwise flagged as circular...\n",
        "## Circular flags in fasta sequence titles will prevail to select sequences that are rotated.", sep = "")
  }
}

# Deciding on what will be rotated
if (any(circularInName)) { # if some seqs have 'circular' in their name, flagged seq prevails regardless of log
  circSeq <- seq[circularInName]
  linearSeq <- seq[!circularInName]
} else if (circLogProvided && !any(circularInName)) { # if log and no flagged seq, log prevails
  circlatorDf <- read.delim(circlatorLog, stringsAsFactors = FALSE)
  circSeq <- seq[circlatorDf[circlatorDf$circularised == 1, "X.Contig"]]
  linearSeq <- seq[circlatorDf[circlatorDf$circularised == 0, "X.Contig"]]
}  else { # No seqs will be rotated!!!
  warning("No 'circular' hit in sequence names and no circlatorLog file provided either, NONE of the sequences will be rotated!!!!.")
  circSeq <- DNAStringSet()
  linearSeq <- seq
}

# Rotating circular sequences
newCircSeq <- DNAStringSet(
  sapply(X = circSeq, FUN = function(x) {
    seqLength <- nchar(x)
    newStartPos <- round(seqLength/n)
    c(subseq(x, start = newStartPos, end = seqLength), subseq(x, end = newStartPos-1, width = newStartPos-1))
  })
)
# Modify names of rotated circular sequences if necessary
if (length(circSeq) != 0) {
  names(newCircSeq) <- if (any(circularInName)) {
    paste(names(newCircSeq), "rotated", sep = "_")
  } else if (circLogProvided && !any(circularInName)) {
    paste(names(newCircSeq), "circular_rotated", sep = "_")
  }
}
# Building output seq set
newSeq <- c(linearSeq, newCircSeq)

# Sanity check
if(!identical(length(setdiff(nchar(seq), nchar(newSeq))), 0L)) stop("Rotated sequence has a different length than the original one!")
# Save output
writeXStringSet(x = newSeq, filepath = outSeqFile, compress = FALSE)

# Log messages
summaryDf <- data.frame(
  seqNames = names(newSeq),
  seqLength = nchar(newSeq),
  StartAtFormerPosition = round(nchar(newSeq)/n)
)
cat("\n")
cat("***********************************************************************************************\n")
if (!circLogProvided & !any(circularInName)) {
  cat("##", date(), ": No 'circular' hit in sequence names and no circlatorLog file provided either, NOTHING has been rotated!!!!")
  print(knitr::kable(summaryDf[, c("seqNames", "seqLength")]))
} else if (length(circSeq) != 0) {
  cat("##", date(), ": Summary of circular sequences rotation in '", seqFile, "' :")
  print(knitr::kable(summaryDf))
} else {
  cat("##", date(), ": No 'circular' hit in sequence names and none of the sequences in '", seqFile, "' were specified as circular in circlatorLog file and nothing has been rotated!")
  print(knitr::kable(summaryDf[, c("seqNames", "seqLength")]))
}
cat("\n")
cat("All sequences saved in:", outSeqFile, "\n")
cat("***********************************************************************************************\n")

quit(save = "no")
