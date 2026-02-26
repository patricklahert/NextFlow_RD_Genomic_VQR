// Use newest nextflow dsl
nextflow.enable.dsl = 2

process freebayes {
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }
    container 'staphb/freebayes:1.3.10'

    tag "$bamFile"

    input:
    tuple val(sample_id), file(bamFile), file(bamIndex)
    path indexFiles

    output:
    tuple val(sample_id), file("*.vcf")

    script:
    """
    echo "Running FreeBayes for Sample: ${bamFile}"

    # Validate input BAM
    if [[ ! -f "${bamFile}" ]]; then
        echo "ERROR: BAM file not found: ${bamFile}" >&2
        exit 1
    fi
    if [[ ! -s "${bamFile}" ]]; then
        echo "ERROR: BAM file is empty: ${bamFile}" >&2
        exit 1
    fi

    # Validate BAM index
    if [[ ! -f "${bamIndex}" ]]; then
        echo "ERROR: BAM index not found: ${bamIndex}" >&2
        exit 1
    fi

    # Resolve genome FASTA
    if [[ -n "${params.genome_file}" ]]; then
        genomeFasta=\$(basename ${params.genome_file})
    else
        genomeFasta=\$(find -L . -name '*.fasta' -o -name '*.fa' | head -1)
    fi

    if [[ -z "\${genomeFasta}" || ! -f "\${genomeFasta}" ]]; then
        echo "ERROR: Reference genome FASTA could not be found in the work directory" >&2
        exit 1
    fi

    echo "Genome File: \${genomeFasta}"

    freebayes -f "\${genomeFasta}" ${bamFile} > "${sample_id}.vcf"

    # Verify VCF output was produced
    if [[ ! -s "${sample_id}.vcf" ]]; then
        echo "ERROR: Output VCF is missing or empty for sample: ${sample_id}" >&2
        exit 1
    fi

    echo "Variant Calling for Sample: ${sample_id} Complete"
    """
}