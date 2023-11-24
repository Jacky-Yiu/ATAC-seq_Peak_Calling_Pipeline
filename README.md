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

The ChIP-seq pipeline aims to:

- Perform quality control on raw data.
- Trim and preprocess reads.
- Align reads to a reference genome.
- Mark duplicate reads.
- Call peaks using MACS2.
- Perform downstream analysis using CHIPSEEKER.

**Dependencies:**

- Nextflow
- Singularity
- FastQC
- Trimmomatic
- Minimap2
- Samtools
- Picard
- MACS2
- CHIPSEEKER
- R (for CHIPSEEKER)

## Usage

### Installation

1. **Clone the Repository:**

    ```bash
    git clone https://github.com/your/repository.git
    cd repository
    ```

2. **Create Conda Environment (if applicable):**

    ```bash
    conda env create -f environment.yml
    conda activate your_env_name
    ```

3. **Download Reference Genome:**

    ```bash
    nextflow run your_pipeline.nf --genome
    ```

### Input

- Paired-end FASTQ files with proper naming conventions (e.g., sample1_1.fq.gz, sample1_2.fq.gz).
- Specify the input data path in the `params.reads` parameter in your Nextflow script.

### Exact Step-by-Step Usage

1. **Quality Control and Pre-processing:**

    ```bash
    nextflow run your_pipeline.nf --preprocess
    ```

2. **Alignment:**

    ```bash
    nextflow run your_pipeline.nf --align
    ```

3. **Mark Duplicates:**

    ```bash
    nextflow run your_pipeline.nf --mark_duplicates
    ```

4. **Peak Calling:**

    ```bash
    nextflow run your_pipeline.nf --peak_calling
    ```

5. **Downstream Analysis (CHIPSEEKER):**

    ```bash
    nextflow run your_pipeline.nf --chipseeker
    ```

### Singularity Configuration

```groovy
process {
    withName:PRE_FASTQC {
        container = 'staphb/fastqc:latest'
    }
    // ... (Other processes with container specifications)
}

singularity {
    enabled = true
    autoMounts = true
}
