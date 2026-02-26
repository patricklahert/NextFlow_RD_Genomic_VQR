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

    # Validate input FASTA
    if [[ ! -f "${genomeFasta}" ]]; then
        echo "ERROR: Genome FASTA not found: ${genomeFasta}" >&2
        exit 1
    fi
    if [[ ! -s "${genomeFasta}" ]]; then
        echo "ERROR: Genome FASTA is empty: ${genomeFasta}" >&2
        exit 1
    fi

    bowtie2-build --large-index ${genomeFasta} ${basename}

    # Verify at least one index file was produced
    if ! ls ${basename}*.bt2l 1>/dev/null 2>&1; then
        echo "ERROR: Bowtie2 index files not created for basename: ${basename}" >&2
        exit 1
    fi

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

    # Validate input reads
    if [[ ! -f "${reads[0]}" ]]; then
        echo "ERROR: R1 reads file not found: ${reads[0]}" >&2
        exit 1
    fi
    if [[ ! -f "${reads[1]}" ]]; then
        echo "ERROR: R2 reads file not found: ${reads[1]}" >&2
        exit 1
    fi

    # Validate index files are present (support both small .bt2 and large .bt2l index formats)
    if ! ls ${genome_basename}*.bt2* 1>/dev/null 2>&1; then
        echo "ERROR: Bowtie2 index files not found for basename: ${genome_basename}" >&2
        exit 1
    fi

    bowtie2 -x ${genome_basename} -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam

    # Verify output SAM was produced and is non-empty
    if [[ ! -s "${sample_id}.sam" ]]; then
        echo "ERROR: Output SAM file is missing or empty for sample: ${sample_id}" >&2
        exit 1
    fi

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

    # Validate input SAM
    if [[ ! -f "${samFile}" ]]; then
        echo "ERROR: SAM file not found: ${samFile}" >&2
        exit 1
    fi
    if [[ ! -s "${samFile}" ]]; then
        echo "ERROR: SAM file is empty: ${samFile}" >&2
        exit 1
    fi

    samtools view -bS ${samFile} | \
        samtools addreplacerg \
            -r "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:illumina\tLB:${sample_id}" \
            - > ${sample_id}.bam

    # Verify output BAM was produced and is non-empty
    if [[ ! -s "${sample_id}.bam" ]]; then
        echo "ERROR: Output BAM file is missing or empty for sample: ${sample_id}" >&2
        exit 1
    fi

    echo "SAM to BAM conversion complete"
    """
}