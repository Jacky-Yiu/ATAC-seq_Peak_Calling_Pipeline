#!/usr/bin/env nextflow

// Define input parameters
params.reads = "$projectDir/input/*{1,2}.fastq.gz"
params.outputDir = "output"

log.info """\
    P I P E L I N E
    ===================================
    reads        : ${params.reads}
    outdir       : ${params.outputDir}
    """
    .stripIndent()


// Conda environment for FastQC
process DOWNLOAD_REFERENCE {
    tag "Download Reference"
    output:
    path "hg38.fa.gz"

    script:
    """
    curl  https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz > hg38.fa.gz
    """
}

process PRE_FASTQC {
    tag "PRE_FASTQC on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("pre_align_fastqc/${sample_id}")

    script:
    """
    mkdir -p pre_align_fastqc
    mkdir -p pre_align_fastqc/${sample_id}
    fastqc -o pre_align_fastqc/${sample_id} -f fastq -q $reads
    """ 
}

process PRE_MULTIQC {
    tag "PRE_MULTIQC on $sample_id"
    publishDir params.outputDir, mode:'copy'

    input:
    tuple val(sample_id), path(fastqc)

    output:
    path "${sample_id}_pre_trimming_multiQC_report.html"

    script:
    """
    multiqc --filename ${sample_id}_pre_trimming_multiQC_report $fastqc
    """
}

process POST_FASTQC {
    tag "POST_FASTQC on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("post_align_fastqc/${sample_id}")

    script:
    """
    mkdir -p post_align_fastqc
    mkdir -p post_align_fastqc/${sample_id}
    fastqc -o post_align_fastqc/${sample_id} -f fastq -q $reads
    """ 
}

process POST_MULTIQC {
    tag "POST_MULTIQC on $sample_id"
    publishDir params.outputDir, mode:'copy'

    input:
    tuple val(sample_id), path(fastqc)

    output:
    path "${sample_id}_post_trimming_multiQC_report.html"

    script:
    """
    multiqc --filename ${sample_id}_post_trimming_multiQC_report $fastqc
    """
}

// Trimmomatic process
process TRIMMOMATIC {
    tag "TRIMMOMATIC on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("trimmed_output/${sample_id}_{forward,reverse}_trimmed.fastq")

    script:
    """
    mkdir trimmed_output
    trimmomatic PE -threads 2 $reads \
    trimmed_output/${sample_id}_forward_trimmed.fastq \
    trimmed_output/${sample_id}_forward_unpaired.fastq \
    trimmed_output/${sample_id}_reverse_trimmed.fastq \
    trimmed_output/${sample_id}_reverse_unpaired.fastq \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// minimap2 process
process MINIMAP2_SAMTOOLS {
    tag "MINIMAP2 on $sample_id"
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    path(reference)

    output:
    tuple val(sample_id), path("minimap2_output/${sample_id}_sorted.bam")

    script:
    """
    mkdir -p minimap2_output
    minimap2 -a -x sr -Y -K 100M $reference $trimmed_reads | samtools view -bS | samtools sort -o minimap2_output/${sample_id}_sorted.bam
    """
}

// MARKDUP process
process MARKDUP {
    tag "MARKDUP on $sample_id"
    
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("markdup_output/${sample_id}_marked_duplicates.bam")
    

    script:
    """
    mkdir -p markdup_output
    java "-Xmx60g" -jar /usr/picard/picard.jar MarkDuplicates \
      I=$bam \
      O=markdup_output/${sample_id}_marked_duplicates.bam \
      M=markdup_output/marked_dup_metrics.txt
    """
}

process MACS2 {
    tag "MACS2 on $sample_id"
    publishDir params.outputDir, mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), 
    path("MACS2_output/${sample_id}/${sample_id}_peaks.narrowPeak"), 
    path("MACS2_output/${sample_id}/${sample_id}_peaks.xls"), 
    path("MACS2_output/${sample_id}/${sample_id}_summits.bed")
    

    script:
    """
    mkdir -p MACS2_output
    macs2 callpeak \
    --treatment $bam \
    --name "${sample_id}" \
    --format BAMPE \
    --nomodel \
    --keep-dup all \
    --gsize hs \
    --qvalue 0.05 \
    --outdir MACS2_output/${sample_id}\
    
    """
}

process CHIPSEEKER {
    tag "CHIPSEEKER on $sample_id"
    publishDir params.outputDir, mode: 'copy'
    
    input:
    tuple val(sample_id), 
    path(narrowPeak), 
    path(peaks_xls), 
    path(summits_bed)

    path(script)

    output:
    tuple val(sample_id), path("chipseeker_output/${sample_id}/")
    

    script:
    """
    mkdir -p chipseeker_output
    mkdir -p chipseeker_output/${sample_id}
    Rscript $script chipseeker_output/${sample_id} $narrowPeak
    
    """
}

// Define the workflow
workflow {
    reference_ch = DOWNLOAD_REFERENCE()

    script_ch = file("$projectDir/script/annotate_MACS2_peaks.R")
    
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    pre_fastqc_ch = PRE_FASTQC(read_pairs_ch)
    PRE_MULTIQC(pre_fastqc_ch)

    trimmomatic_ch = TRIMMOMATIC(read_pairs_ch)

    post_fastqc_ch = POST_FASTQC(trimmomatic_ch)
    
    POST_MULTIQC(post_fastqc_ch)

    minimap2_ch = MINIMAP2_SAMTOOLS(trimmomatic_ch, reference_ch)

    markdup_ch = MARKDUP(minimap2_ch)

    macs2_ch = MACS2(markdup_ch)

    CHIPSEEKER(macs2_ch, script_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! The ouput is in --> $projectDir/$params.outputDir\n" : "Oops .. something went wrong" )
}