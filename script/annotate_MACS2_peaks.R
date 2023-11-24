library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Annotate genomic intervals in bed format using ChIPpeakAnno 
# This script was designed for Arabidopsis, but can be easily changed for
# any other organism available through biomaRt
args <- commandArgs(trailingOnly = TRUE)

setwd(args[1])
list.files()

annotateMacsPeaks <- function(peakFile) 
{
  peaks <- readPeakFile(peakFile)
  
  peakAnno <- annotatePeak(peaks, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
  # Read in results table
  

  #cov_plot <- covplot(peaks, weightCol="V5")

  jpeg(file="piechart.jpeg")

  pie_chart <- plotAnnoPie(peakAnno)

  dev.off()
    
}

lapply( args[2], function(x) annotateMacsPeaks(x) )