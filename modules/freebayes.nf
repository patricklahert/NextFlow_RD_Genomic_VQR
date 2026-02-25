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

    if [[ -n ${params.genome_file} ]]; then
        genomeFasta=\$(basename ${params.genome_file})
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi  

    freebayes -f "\${genomeFasta}" ${bamFile} > "${sample_id}.vcf"
    
    echo "Variant Calling for Sample: ${sample_id} Complete"
    """
}