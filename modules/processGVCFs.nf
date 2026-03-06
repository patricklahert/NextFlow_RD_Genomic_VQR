process combineGVCFs {
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/gatk4:4.3.0.0'
    tag "${sample_ids.join('_')}" // Add a tag based on the sample IDs

    input:
    tuple val(sample_ids), path(gvcf_files), path(gvcf_index_files)
    path indexFiles

    output:
    tuple val("${sample_ids.join('_')}"), file("*_combined.vcf"), file("*_combined.vcf.idx")

    script:
    def merged_sample_id = "${sample_ids.join('_')}"
    def gvcf_files_args = gvcf_files.collect { file -> "-V ${file}" }.join(' ')

    """
    echo "Combining GVCFs for samples: ${gvcf_files.collect { it.baseName }.join(', ')}"

    genomeFasta="\$(find -L . -name '*.fasta')"

    # Ensure dictionary exists
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    gatk CombineGVCFs -R "\${genomeFasta}"\
        ${gvcf_files_args} \
        -O ${merged_sample_id}_combined.vcf
    """
}

// Merge per-sample FreeBayes VCFs into a multi-sample VCF using bcftools merge.
process bcftoolsMergeVCFs {
    if (params.platform == 'local') {
        label 'process_medium'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }
    container 'staphb/bcftools:1.21'

    tag "${sample_ids.join('_')}"

    input:
    tuple val(sample_ids), path(vcf_files)

    output:
    tuple val("${sample_ids.join('_')}"), file("*_merged.vcf.gz"), file("*_merged.vcf.gz.tbi")

    script:
    def merged_sample_id = "${sample_ids.join('_')}"
    def vcf_list = vcf_files instanceof List ? vcf_files.collect { it.toString() }.join(' ') : vcf_files.toString()
    def num_samples = vcf_files instanceof List ? vcf_files.size() : 1

    """
    echo "Merging FreeBayes VCFs for: ${sample_ids.join(', ')}"

    # Compress and index each VCF using bcftools (avoids needing bgzip/tabix separately)
    for vcf in ${vcf_list}; do
        if [[ ! -s "\${vcf}" ]]; then
            echo "ERROR: VCF file is missing or empty: \${vcf}" >&2
            exit 1
        fi
        bcftools view "\${vcf}" -Oz -o "\${vcf}.gz"
        bcftools index -t "\${vcf}.gz"
    done

    if [[ "${num_samples}" -gt 1 ]]; then
        bcftools merge *.vcf.gz -Oz -o ${merged_sample_id}_merged.vcf.gz
    else
        bcftools view *.vcf.gz -Oz -o ${merged_sample_id}_merged.vcf.gz
    fi

    bcftools index -t ${merged_sample_id}_merged.vcf.gz

    if [[ ! -s "${merged_sample_id}_merged.vcf.gz" ]]; then
        echo "ERROR: Merged VCF is missing or empty" >&2
        exit 1
    fi

    echo "bcftools merge complete: ${merged_sample_id}_merged.vcf.gz"
    """
}

process genotypeGVCFs {
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/gatk4:4.3.0.0'
    tag "$combined_sample_id"

    input:
    tuple val(combined_sample_id), file(combined_gvcf), file(combined_gvcf_idx)
    path indexFiles

    output:
    tuple val(combined_sample_id), file("*_genotyped.vcf"), file("*_genotyped.vcf.idx")

    script:
    def merged_sample_id = combined_gvcf.baseName

    """
    echo "Genotyping combined GVCF: ${combined_gvcf.baseName}"

    if [[ -n ${params.genome_file} ]]; then
        genomeFasta=\$(basename ${params.genome_file})
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi

    echo "Genome File: \${genomeFasta}"

    # Rename the dictionary file to the expected name if it exists
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    gatk GenotypeGVCFs -R "\${genomeFasta}" \
        -V ${combined_gvcf} \
        -O ${merged_sample_id}_genotyped.vcf

    """
}