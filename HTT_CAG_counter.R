# R script for CAG repeat counting in human Huntingtin (HTT) gene from nanopore reads
# Developed by Harness Therapeutics Limited

# load libraries
suppressPackageStartupMessages(library("Biostrings"))
suppressPackageStartupMessages(library("ShortRead"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("openxlsx"))
source("helper_functions.R")

# Define command-line options
option_list <- list(
  make_option(c("-m", "--mapped_fq"), type = "character", 
              help = "Path to the input FASTQ file containing reads mapping to human HTT gene"),
  make_option(c("-n", "--name"), type = "character", default = "sample", 
              help = "Name to append to output files"),
  make_option(c("-p", "--primer_seq_file"), type = "character", 
              help = "Path to the input csv file containing the primer sequences used for PCR amplification of human HTT gene")
)

# Parse the arguments
parser <- OptionParser(option_list = option_list, 
                       description = "Harness Therapeutics CAG Repeat Counter - Analyze HTT CAG repeat length from nanopore reads.",
                       usage = "Rscript HTT_CAG_counter.R --mapped_fq <path-to-mapped-fq> 
                       --name <name> 
                       --primer_seq_file <path-to-primer-csv>")
args <- parse_args(parser)

# Validate required input
if (is.null(args$mapped_fq)) {
  print_help(parser)
  stop("Error: --mapped_fq is required.\n", call. = FALSE)
} else if (is.null(args$primer_seq_file)) {
  print_help(parser)
  stop("Error: --primer_seq_file is required.\n", call. = FALSE)
}

mapped_fq <- args$mapped_fq
name <- args$name
primer_seq_file <- args$primer_seq_file

# Check inputs in correct format, and exist
if (!grepl("\\.f.*q\\.gz$", args$mapped_fq)) {
  stop("--mapped_fq must be a .fastq.gz file.")
}
if (!file.exists(args$mapped_fq)) {
  stop("Input file does not exist: ", args$mapped_fq)
}


if (!grepl("\\.csv$", args$primer_seq_file)) {
  stop("--primer_seq_file must be a .csv file.")
}
if (!file.exists(args$primer_seq_file)) {
  stop("Input file does not exist: ", args$primer_seq_file)
}
# Load primer sequence file
message("\n")
message("===================")
message("Reading primer sequences")
primer_seq_df <- read.csv(primer_seq_file)
expected_cols <- c("Primer", "Sequence")

if (!all(expected_cols %in% colnames(primer_seq_df))) {
  stop("Input primer csv is missing one or more expected columns.")
}

####
# Sanity check for the input primers
####

# (1) check primer sequence to avoid space
primer_seq_df$Sequence <- gsub(" ", "", primer_seq_df$Sequence)

# (2) setup list_primers object
## a. set up DNAStringSet obj
primer_seq_dss <- DNAStringSet(primer_seq_df$Sequence)
names(primer_seq_dss) <- primer_seq_df$Primer
## b. Reverse complement primer seqs
primer_seq_rc_dss <- reverseComplement(primer_seq_dss)
names(primer_seq_rc_dss) <- paste0(names(primer_seq_dss),"_RC")

list_primers <- list("F" = primer_seq_dss$`F`, "R" = primer_seq_dss$`R`,
                     "F_RC" = primer_seq_rc_dss$F_RC, "R_RC" = primer_seq_rc_dss$R_RC)

# (3) check primer sequences were designed for HTT gene
## a. load HTT gene fasta file
htt_fasta <- "example_data/HTT_gene.fasta"
seq_HTTgene <- readDNAStringSet(htt_fasta, format = "fasta")
seq_HTTgene <- unlist(seq_HTTgene)

## b. check primers perfectly match to HTT gene only once
# for check the location of primers in HTT gene without error
message("===================")
message("Check primer information")
message("===================")
message("* Are the primers designed for the human HTT gene?")
count_HTT_map2_5primer <- countPattern(unlist(list_primers$`F`), seq_HTTgene)
count_HTT_map2_3primer <- countPattern(unlist(list_primers$R_RC), seq_HTTgene)
if (count_HTT_map2_5primer == 1 & count_HTT_map2_3primer == 1){
  message("\tPrimer check passed...")
}else{
  message("\tPrimer check failed. Program stopped.")
  message("===================")
  stop("Please make sure primers are designed for human HTT gene and the primer information is correct.")
}
HTT_map2_5primer <- matchPattern(unlist(list_primers$`F`), seq_HTTgene)
HTT_map2_3primer <- matchPattern(unlist(list_primers$R_RC), seq_HTTgene)

# check the primer locations can identity CAG track and interrupting motif
## get the HTT sequence between primers
HTT_seq_between <- seq_HTTgene[c((Biostrings::end(HTT_map2_5primer) + 1) : (Biostrings::start(HTT_map2_3primer) - 1))]
## make sure targeted HTT gene contains CAG_19 and CAACAGCCGCCA
# We use number of CAG = 19 as the number of CAG from the reference allele
ref_cag_seq <- DNAString(paste(rep("CAG", 19), collapse = ""))
ref_int_motif_seq <- DNAString("CAACAGCCGCCA")

count_map2_wt_cag <- countPattern(ref_cag_seq, HTT_seq_between)
count_map2_int_motif <- countPattern(ref_int_motif_seq, HTT_seq_between)
message("===================")
message("* Does the targeted PCR product contain the CAG tract and interrupting motif from the HTT gene?")
if (count_map2_wt_cag == 1 & count_map2_int_motif == 1){
  message("\tPCR product check passed...")
}else{
  message("\tPCR product check failed. Program stopped.")
  message("===================")
  stop("Please make sure primers are designed for human HTT gene and the primer information is correct.")
}

# check the distance between the end of the 5' primer and the start of the CAG tract on the reference allele.
wt.cag_map2_pcr.product <- matchPattern(ref_cag_seq, HTT_seq_between)
dis.between_5prime.end_cag.tract <- start(wt.cag_map2_pcr.product) - 1

# Define key DNA sequences
## Trinucleotides
cag_dna <- DNAString("CAG")
cag_dna_rc <- reverseComplement(cag_dna)
caa_dna <- DNAString("CAA")
car_dna <- DNAString("CAR")
car_dna_rc <- reverseComplement(car_dna)
## Interrupting motif
intrr_pattern <- "CAACAGCCGCCA"
intrr_pattern_dna <- DNAString(intrr_pattern)
intrr_pattern_dna_rc <- reverseComplement(intrr_pattern_dna)
intrr_pattern_iupac <- "CARCARCCRCCR"
intrr_pattern_iupac_dna <- DNAString(intrr_pattern_iupac)
intrr_pattern_iupac_dna_rc <- reverseComplement(intrr_pattern_iupac_dna)
intrr_pattern_poly_p <- "CCGCCA"
intrr_pattern_poly_p_dna <- DNAString(intrr_pattern_poly_p)
intrr_pattern_poly_p_dna_rc <- reverseComplement(intrr_pattern_poly_p_dna)

# Input fastq files
inputread_list <- list("mapped_fq" = mapped_fq)

interrupt_motif_reads <- NULL

for (m in seq_along(inputread_list)) {
  indfq <- inputread_list[[m]]
  show_read_cat <- names(inputread_list)[m]
  message("===================")
  message("Processing ", show_read_cat)
  message("===================")
  # Parse fastq file
  fastq <- readFastq(indfq)
  num_read <- length(fastq)
  read_len <- width(fastq)
  read_id <- as.character(id(fastq))
  names(read_len) <- read_id
  reads_dna <- sread(fastq)
  names(reads_dna) <- read_id
  
  if (length(reads_dna) > 0) {
    message("\tCounting reads containing primer sequences")
    # Identify reads containing PCR primers
    all_map_primer_fr_5prime <- vcountPattern(list_primers$`F`, reads_dna, max.mismatch = 1)
    all_map_primer_fr_3prime <- vcountPattern(list_primers$R_RC, reads_dna, max.mismatch = 1)
    
    all_map_primer_rc_5prime <- vcountPattern(list_primers$`R`, reads_dna, max.mismatch = 1)
    all_map_primer_rc_3prime <- vcountPattern(list_primers$F_RC, reads_dna, max.mismatch = 1)
    
    # Identify reads containing interrupting pattern
    message("\tCounting reads containing interrupting motif")
    all_int_pattern_1 <- vcountPattern(intrr_pattern_iupac_dna, reads_dna, fixed = F)
    all_int_pattern_2 <- vcountPattern(intrr_pattern_iupac_dna_rc, reads_dna, fixed = F)
    primer_df <- data.frame(read_name = names(reads_dna),
                            read_len = read_len[names(reads_dna)],
                            FR_5primeFR = all_map_primer_fr_5prime,
                            FR_3primerc = all_map_primer_fr_3prime,
                            RC_5primeFR = all_map_primer_rc_5prime,
                            RC_3primerc = all_map_primer_rc_3prime,
                            InterruptMotif_FR = all_int_pattern_1,
                            InterruptMotif_RC = all_int_pattern_2)
    # Add interrupting motif info to master table
    primer_dfchecked <- Fun_check_Primer_InterruptMotif(primer_df)
    # Interrupt motif loci
    int_pattern_loci_df <- Fun_Primer_InterruptMotif_Loci(primer_dfchecked, intrr_pattern_iupac_dna, reads_dna, list_primers)
    if (!is.null(int_pattern_loci_df)) {
      int_pattern_loci_df$QualityControl <- show_read_cat
      interrupt_motif_reads <- rbind(interrupt_motif_reads, int_pattern_loci_df)
    }
  }
}
message("===================")

# interrupt_motif_reads
# If not enough reads pass filter, interrupt_motif_reads will be null
if (is.null(interrupt_motif_reads)) {
  interrupt_motif_reads <- data.frame()
} else {
  interrupt_motif_2_primer_NumTriNucl_df <- Fun_Calculate_Trinucleotide(interrupt_motif_reads, dis.between_5prime.end_cag.tract)
}

obj2save <- paste0(name, "_HTTcagEstimationTable.xlsx")
message(paste0("Writing reads containing interrupting motif to ", obj2save))
write.xlsx(interrupt_motif_2_primer_NumTriNucl_df, file = obj2save)

message("===================")
message("Finished, bye!")
message("===================")
