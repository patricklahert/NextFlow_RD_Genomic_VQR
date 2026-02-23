/*
 * Run fastq on the read fastq files
 */
process fastp {
    
    label 'process_single'

    container 'staphb/fastp:0.23.4'

    // Add a tag to identify the process
    tag "$sample_id"

    // Specify the output directory for the fastp results
    publishDir("$params.outdir/fastp", mode: "copy")

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastp_${sample_id}_trimmed/*"

    script:
    """
    echo "Running fastp"
    mkdir -p fastp_${sample_id}_trimmed  

    fastp -i ${reads[0]} \
    -I ${reads[1]} \
    -o fastp_${sample_id}_trimmed/${sample_id}_R1_trimmed.fastq.gz \
    -O fastp_${sample_id}_trimmed/${sample_id}_R2_trimmed.fastq.gz \
    -h fastp_${sample_id}_trimmed/${sample_id}_fastp_report.html

    echo "fastp Complete"
    """
}


