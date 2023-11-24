library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <output_directory> <input_file>")
}

output_dir <- args[1]
input_file <- args[2]

# Check if input file exists
if (!file.exists(input_file)) {
  stop("Input file does not exist.")
}

# Annotate peaks
annotateMacsPeaks <- function(peakFile, output_dir) {
  peaks <- readPeakFile(peakFile)
  peakAnno <- annotatePeak(peaks, tssRegion = c(-3000, 3000),
                           TxDb = txdb, annoDb = "org.Hs.eg.db")

  # Save plot directly to the output directory
  jpeg(file = file.path(output_dir, "piechart.jpeg"))
  plotAnnoPie(peakAnno)
  dev.off()
}

# Annotate peaks with the given input file
annotateMacsPeaks(input_file, output_dir)
