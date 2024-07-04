# 20240704 Extract Fasta
Gabriel (my lab mate) and I were struggling with extracting read sequences from BAM files. While samtools view extracts reads with the full length of the original reads, we wanted to have the reads from an exact region. So we created this script to address that need.
# How to Use

Ensure you have R installed with the following packages:

```r=
optparse
Rsamtools
GenomicRanges
Biostrings
```

Run the script using this command:

```r=
Rscript extract_fasta.R -i <input_bam> -o <output_fasta> -r <region>

<input_bam> is your input BAM file
<output_fasta> is the name for your output FASTA file
<region> is the genomic region of interest (format: chr:start-end)
```

Example:
CopyRscript extract_fasta.R -i sample.bam -o extracted_sequences.fasta -r chr1:100000-100100

The script will extract the sequences from the specified region and save them in your output FASTA file.

# Output
The output FASTA file will contain sequences that match your region of interest. The FASTA headers will include information about the original read and its position.
