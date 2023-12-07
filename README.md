# BIOF501 Term Project: ATAC-seq Peak-calling and Annotation Pipeline

### *By Jacky Yiu*

------------------------------------------------------------------------

## Repository Contents

### Directories

-   `Inputs` : Contains inputs for the pipeline, which are pair-end reads in fastq format.

-   `script` : Contains the R scripts used in the pipeline for [ChIPseeker](https://bioconductor.org/packages/ChIPseeker/) peaks annotations.

### Files

-   dag.png: Representation of the overall workflow in PNGformat.

-   nextflow.config: Config file for Nextflow, contain all Docker container info and setting needed for using Singularity for this pipeline

-   workflow.nf: The main pipeline file that contain the instructions for running the pipeline through Nextflow.

------------------------------------------------------------------------

## Introduction

### Background

### Purpose

### Rationale

------------------------------------------------------------------------

## Usage

### Installation

Installing this pipeline requires `git` , `nextflow`, and `singularity`

Click [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) for instruction to install `git`

Click [here](https://www.nextflow.io/docs/latest/getstarted.html) for instruction to install `nextflow`

Click [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) for instruction to install `singularity`

After installing all three, clone the repository by running the following command in a terminal:

```         
git clone https://github.com/Jacky-Yiu/ATAC-seq_Peak_Calling_Pipeline.git
```

### Running the Pipeline

Th command above should download the pipeline in your local directory, enter the pipeline directory using the following command:

```         
cd ./ATAC-seq_Peak_Calling_Pipeline
```

The pipeline can then be run by typing in:

```         
nextflow run workflow.nf
```

### Pipeline Overview

The following is a simplified representation of the workflow in a Directed Acyclic Graph (DAG):

![DAG](dag.png)

This pipeline will take in pair-end Fastq generatedby ATAC-Seq and go through 4 phases

1.  Pre-alignment QC

2.  Alignment

3.  Post-alignment QC

4.  Peak Calling

5.  Peak Annotation

### Data Used

Two files are used as inputs; *WNV_Genomes.fasta* and *Reference.fasta*. These can be found in the Viral_Sequence_Data directory. *WNV_Genomes.fasta* includes all complete WNV genomes sampled in North America that are listed on [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:11082). Reference.fasta is an in-frame coding sequence for the genome.

The Genomes file can be altered to include or exclude different genomes to change the search space. For example, if I wanted to search for selected sites in European or Asian WNV variants, I would use a genomes file with only those sequences.

Changing the reference file can alter the scope of the analysis; e.g., if you wanted to determine if only a specific segment of the genome is under selection, like the envelope protein, you could change the reference file to only include that sequence.

Note: The reference file *must* be in frame for the pipeline to function correctly.

------------------------------------------------------------------------

## Outputs and Results

### Pipeline outputs

The four outputs of this pipeline are as follows:

1.  Pre-trim FastQC Report

2.  Post-trim FastQC Report

3.  Peak File

4.  Peak Annotation Plot

### Results

### Future directions

------------------------------------------------------------------------

## References {#references}
