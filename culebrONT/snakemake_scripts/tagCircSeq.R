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
              help = "Path of fasta file of sequences to be flagged"),
  make_option(c("-o", "--outFilePath"),
              type = "character",
              default = NULL,
              help = "Path where output sequence file will be written."),
  make_option(c("-d", "--logFile"),
              type = "character",
              default = NULL,
              help = "Path to a log file that will be parsed to extract circular status info for tagging."),
  make_option(c("-a", "--assembler"),
              type = "character",
              default = NULL,
              help = "Wildcard assembler.")
)

##### RETRIEVEING PARAMS from optparse #####
myArgs <- parse_args(
  OptionParser(usage = "%prog [options]", option_list = option_list,
               description = "Add a 'circular' tag to the title of sequences in a fasta file. Tagging is based on information from a log file provided by the assembler or Circlator.")
)

# Assign parameter values
seqFile <- myArgs$seqFile
outSeqFile <- myArgs$outFilePath
if(is.null(outSeqFile)) outSeqFile <- paste0(gsub("(^.*)\\..*$", "\\1", seqFile), "_circFlagged.fasta")
logFile <- myArgs$logFile
if (is.null(logFile)) stop("Log file must be provided. Run program with '--help' for more info on usage.")
assembler <- myArgs$assembler

#######################################################################################################################
#######################################################################################################################

### POSSIBLY NEED TO CHANGE BEHAVIOR HERE AND IN rotateCircSeqs to accomodate canu titles: suggestCircular=<yes|no>
### THIS WILL DEPEND ON WHETHER CIRCLATOR STRIPS THAT FROM THE TITLE OR NOT

#######################################################################################################################
#######################################################################################################################

circularTag <- "circular"

# Load data
seq <- readDNAStringSet(seqFile)

# Inspect type of log file by loading first line
# flye: "seq_name	length	cov.	circ.	repeat	mult.	alt_group	graph_path"
# circlator: "[merge circularised]	#Contig	repetitive_deleted	circl_using_nucmer	circl_using_spades	circularised"
# Miniasm: "S       utg000001c      "
firstLine <- readLines(con = logFile, n=1)

# Perform tagging accordingly
#if (grepl("cov.\\tcirc.", firstLine)) {
if (assembler=="FLYE") {
  logFileType <- "Flye assembler log"
  logInfo <- read.delim(logFile, stringsAsFactors = FALSE)
  newseq <- seq
  logInfo$title <- ifelse(logInfo$circ. %in% c("Y", "+"), paste(logInfo[[1]], circularTag, sep = "_"), logInfo[[1]])
  names(newseq) <- logInfo$title[match(names(seq),logInfo[[1]])]

#} else if (grepl("circl_using_spades\\tcircularised", firstLine)) {
}else if (assembler=="CANU" || assembler=="SMARTDENOVO") {
  logFileType <- "Circlator circularisation log"
  logInfo <- read.delim(logFile, stringsAsFactors = FALSE)
  newseq <- seq
  logInfo$title <- ifelse(logInfo$circularised == 1, paste(logInfo$X.Contig, circularTag, sep = "_"), logInfo$X.Contig)
  names(newseq) <- logInfo$title[match(names(seq),logInfo$X.Contig)]

} else if (assembler=="MINIASM") {
#} else if (grepl("^S\\t.*\\t[ATGC]", firstLine)) {
  logFileType <- "MINIASM GFA file"
  circDiagLines <- grep("^L\\t(.*?)\\t\\1.*0M$",
                        readLines(con = logFile),
                        value = TRUE)
  newseq <- seq
  if (length(circDiagLines) != 0L) {
    circNames <-unique(sub("^L\\t(.*?)\\t.*", "\\1", circDiagLines))
    names(newseq) <- ifelse(names(newseq) %in% circNames,
                            paste(names(newseq), circularTag, sep = "_"),
                            names(newseq)
                            )
  }

#} else if (grepl("^S\\t.*\\t[ATGC]", firstLine) || grepl("^H\\t.*", firstLine)) {
} else if (assembler=="RAVEN" || assembler=="SHASTA") {
  logFileType <- "RAVEN or SHASTA file"
  circDiagLines <- grep("^L\\t(.*?)\\t\\1.*",
                        readLines(con = logFile),
                        value = TRUE)
  newseq <- seq
  if (length(circDiagLines) != 0L) {
    circNames <-unique(sub("^L\\t(.*?)\\t.*", "\\1", circDiagLines))
    names(newseq) <- ifelse(names(newseq) %in% circNames,
                            paste(names(newseq), circularTag, sep = "_"),
                            names(newseq)
                            )
  }

}else {
  stop("The log file provided to extract circular status of sequences has an unrecognized type!")
}

cat("\n")
cat("***********************************************************************************************\n")
cat("##", date(), assembler, "\n")
cat("o Invoking command :", commandArgs(), fill = TRUE, labels = "## ")
cat("##  o Provided log file was processed as a '", logFileType, "'.\n", sep = "")
cat("##  o Sequence titles after tagging for circularity:\n")
cat(paste("##     ->", names(newseq), sep = "    ", collapse = "\n"), "\n")

writeXStringSet(x = newseq, filepath = outSeqFile, compress = FALSE)

cat("o Tagged sequences saved in:", outSeqFile, fill = TRUE, labels = "## ")
cat("***********************************************************************************************\n")



quit(save = "no", status = 0, runLast = FALSE)
