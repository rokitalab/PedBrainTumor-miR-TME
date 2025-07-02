# Creaete 3' UTR fasta from Gencode GTF

library(GenomicFeatures) # for creating TxDb from GTF
library(Biostrings) # for handling and writing DNAStringSet
library(rtracklayer) # for importing GTF
library(BSgenome.Hsapiens.UCSC.hg38)  # pre-built genome for hg38


# Set directory paths

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "human-novel-mirna-target-prediction")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

Set GTF file path 

gtf_file <- file.path(data_dir, "gencode.v39.primary_assembly.annotation.gtf.gz")

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

utr_fasta_path <- file.path(results_dir, "gencode.v39.3utr.fa")

# Save the extracted 3′ UTR sequences to a FASTA file

writeXStringSet(utr_seqs, filepath = utr_fasta_path)
