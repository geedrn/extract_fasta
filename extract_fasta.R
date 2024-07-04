#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Rsamtools)
  library(GenomicRanges)
  library(Biostrings)
})

# Command line argument setup
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input BAM file", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output FASTA file", metavar="FILE"),
  make_option(c("-r", "--region"), type="character", default=NULL, 
              help="Genomic region (format: chr:start-end)", metavar="STRING")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Argument check
if (is.null(opt$input) || is.null(opt$output) || is.null(opt$region)) {
  stop("All arguments (input, output, and region) are required.")
}

# Region analysis
region_parts <- strsplit(opt$region, ":")[[1]]
chrom <- region_parts[1]
pos_parts <- as.integer(gsub(",", "", strsplit(region_parts[2], "-")[[1]]))
start <- pos_parts[1]
end <- pos_parts[2]

cat("Analyzing region:", chrom, format(start, big.mark=","), format(end, big.mark=","), "\n")

# BAM file reading and region extraction
param <- ScanBamParam(what=c("qname", "flag", "rname", "pos", "cigar", "seq"),
                      which=GRanges(chrom, IRanges(start, end)))
bam <- scanBam(opt$input, param=param)[[1]]

cat("Number of reads extracted:", length(bam$qname), "\n")

if (length(bam$qname) == 0) {
  cat("No reads found in the specified region.\n")
  file.create(opt$output)
  quit(save = "no", status = 0)
}

# Function to extract sequence from specified region
extract_region <- function(seq, pos, cigar, start, end) {
  cigar_ops <- strsplit(cigar, "(?<=\\d)(?=\\D)|(?<=\\D)(?=\\d)", perl=TRUE)[[1]]
  ops <- cigar_ops[seq(2, length(cigar_ops), 2)]
  lengths <- as.integer(cigar_ops[seq(1, length(cigar_ops), 2)])
  
  ref_pos <- pos
  seq_pos <- 1
  extracted <- ""
  start_in_read <- NULL
  end_in_read <- NULL
  
  for (i in seq_along(ops)) {
    if (ops[i] %in% c("M", "=", "X")) {
      if (ref_pos + lengths[i] - 1 >= start && ref_pos <= end) {
        start_offset <- max(0, start - ref_pos)
        end_offset <- min(lengths[i], end - ref_pos + 1)
        if (is.null(start_in_read)) start_in_read <- seq_pos + start_offset
        end_in_read <- seq_pos + end_offset - 1
        extracted <- paste0(extracted, substr(seq, seq_pos + start_offset, seq_pos + end_offset - 1))
      }
      ref_pos <- ref_pos + lengths[i]
      seq_pos <- seq_pos + lengths[i]
    } else if (ops[i] == "D" || ops[i] == "N") {
      ref_pos <- ref_pos + lengths[i]
    } else if (ops[i] == "I" || ops[i] == "S") {
      seq_pos <- seq_pos + lengths[i]
    }
    if (ref_pos > end) break
  }
  return(list(seq = extracted, start = start_in_read, end = end_in_read))
}

# Extract sequences from each read
extracted_info <- mapply(extract_region, bam$seq, bam$pos, bam$cigar, 
                         MoreArgs = list(start = start, end = end), SIMPLIFY = FALSE)

# Keep only non-empty sequences
extracted_info <- extracted_info[sapply(extracted_info, function(x) nchar(x$seq) > 0)]

# Create FASTA file
if (length(extracted_info) > 0) {
  sequences <- DNAStringSet(sapply(extracted_info, function(x) x$seq))
  valid_reads <- which(sapply(extracted_info, function(x) nchar(x$seq) > 0))
  
  names(sequences) <- paste0(bam$qname[valid_reads], "_", 
                             bam$rname[valid_reads], ":", format(start, big.mark=","), "-", format(end, big.mark=","), 
                             " (Read:", sapply(extracted_info[valid_reads], function(x) x$start), 
                             "-", sapply(extracted_info[valid_reads], function(x) x$end), ")")
  
  writeXStringSet(sequences, opt$output, format="fasta")
  cat(sprintf("Extracted %d sequences to %s\n", length(sequences), opt$output))
  
  # Output debug information
  cat("Specified region length:", end - start + 1, "\n")
} else {
  cat("No sequences found in the specified region.\n")
  file.create(opt$output)
}