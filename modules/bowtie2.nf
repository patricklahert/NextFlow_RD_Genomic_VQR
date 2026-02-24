/*
Build a Bowtie2 genome index from a FASTA file
*/
process indexGenomeBowtie2 {

    if (params.platform == 'local') {
        label 'process_medium'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }
    container 'biocontainers/bowtie2:v2.4.1_cv1'

    publishDir("$params.outdir/GENOME_IDX", mode: "copy")

    input:
    path genomeFasta

    output:
    path "*.bt2l"

    script:
    def basename = genomeFasta.getSimpleName()
    """
    echo "Running Bowtie2 Index"
    bowtie2-build --large-index ${genomeFasta} ${basename}
    echo "Bowtie2 Indexing complete"
    """
}

/*
Align reads to the indexed genome using Bowtie2
*/
process alignReadsBowtie2 {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }
    container 'biocontainers/bowtie2:v2.4.1_cv1'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path bt2_index_files
    val genome_basename

    output:
    tuple val(sample_id), file("${sample_id}.sam")

    script:
    """
    echo "Running Align Reads with Bowtie2"
    bowtie2 -x ${genome_basename} -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam
    echo "Alignment complete"
    """
}

/*
Convert SAM to BAM using samtools
*/
process samToBam {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/indexgenome:1.1.0'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(samFile)

    output:
    tuple val(sample_id), file("${sample_id}.bam")

    script:
    """
    echo "Converting SAM to BAM and adding read groups"
    samtools view -bS ${samFile} | \
        samtools addreplacerg \
            -r "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:illumina\tLB:${sample_id}" \
            - > ${sample_id}.bam
    echo "SAM to BAM conversion complete"
    """
}