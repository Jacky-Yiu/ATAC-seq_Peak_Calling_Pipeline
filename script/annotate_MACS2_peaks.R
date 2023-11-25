library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <output_directory> <input_file>")
}

output_dir <- args[1]
peakFile <- args[2]

# Check if input file exists
if (!file.exists(peakFile)) {
  stop("Input file does not exist.")
}

# Annotate peaks

peaks <- readPeakFile(peakFile)
peakAnno <- annotatePeak(peaks, tssRegion = c(-3000, 3000),
                           TxDb = txdb, annoDb = "org.Hs.eg.db")

# Save plot directly to the output directory
pdf(file = file.path(output_dir, "peak_plot.pdf"))
plotAnnoPie(peakAnno)
upsetplot(peakAnno)
dev.off()

