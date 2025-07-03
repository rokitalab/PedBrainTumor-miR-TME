# Creaete 3' UTR fasta from Gencode GTF

library(GenomicFeatures) # for creating TxDb from GTF
library(Biostrings) # for handling and writing DNAStringSet
library(rtracklayer) # for importing GTF
library(BSgenome.Hsapiens.UCSC.hg38) # pre-built genome for hg38
library(optparse)

parser <- OptionParser()

option_list <- list(
  make_option(c("-g", "--gencode"),
    type = "character",
    default = "../../data/gencode.v39.primary_assembly.annotation.gtf.gz",
    help = "Input genocde file [default %default]", metavar = "FILE"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "results/gencode.v39.3utr.fa",
    help = "Output 3' utr fasta file name [default %default]", metavar = "FILE"
  )
)

parser <- OptionParser(option_list = option_list)

opt <- parse_args(parser)

gtf_file <- opt$gencode # Access the input filename
output_filename <- opt$output # Access the output filename

# Create a transcript database (TxDb)

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Extract 3′ UTR ranges by transcript

three_utrs <- threeUTRsByTranscript(txdb, use.names = TRUE)

# Load the genome reference (hg38)

genome <- BSgenome.Hsapiens.UCSC.hg38
valid_chroms <- intersect(seqlevels(three_utrs), seqnames(genome))

# Extract 3′ UTR sequences using transcript annotations

three_utrs_filtered <- keepSeqlevels(three_utrs, valid_chroms, pruning.mode = "coarse")
utr_seqs <- extractTranscriptSeqs(genome, three_utrs_filtered)

utr_fasta_path <- file.path(output_filename)

# Save the extracted 3′ UTR sequences to a FASTA file

writeXStringSet(utr_seqs, filepath = utr_fasta_path)
