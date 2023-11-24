# ATAC-seq Peak-calling Pipeline

### Background and Rationale
This Standard Operating Procedure outlines the steps for a ChIP-seq data analysis pipeline using Nextflow with Singularity containerization. The pipeline aims to process ChIP-seq raw data, perform quality control, alignment, peak calling, and downstream analysis. It leverages Singularity containers for encapsulating dependencies, ensuring reproducibility across different computing environments.

**Dependencies:**
- Nextflow
- Singularity
- Containers:
  - FastQC (`staphb/fastqc:latest`)
  - MultiQC (`ewels/multiqc:latest`)
  - Trimmomatic (`quay.io/biocontainers/trimmomatic:0.35--6`)
  - Minimap2 and Samtools (`niemasd/minimap2_samtools:latest`)
  - Picard (`broadinstitute/picard:latest`)
  - MACS2 (`dukegcb/macs2:latest`)
  - CHIPSEEKER (`thugenomefacility/chipseeker:0.1`)

### Directed Acyclic Graph (DAG)
The following is a simplified representation of the workflow's Directed Acyclic Graph (DAG):
![DAG](dag.png)

### Usage

#### Installation

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/Jacky-Yiu/ATAC-seq_Peak_Calling_Pipeline.git
   cd ATAC-seq_Peak_Calling_Pipeline
# ChIP-seq Data Analysis Pipeline

## Standard Operating Procedure (SOP)

This Standard Operating Procedure outlines the steps for a ChIP-seq data analysis pipeline. The pipeline processes raw ChIP-seq data through quality control, alignment, peak calling, and downstream analysis. It leverages Singularity containers to encapsulate dependencies and ensure reproducibility.

## Background and Rationale



## Usage

