process variantRecalibrator {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'broadinstitute/gatk:4.6.1.0'

    tag "$vcf"

    publishDir("$params.outdir/VCF", mode: "copy")

    input:
    tuple val(sample_id), file(vcf), file(vcfIndex)
    val knownSitesArgs
    path genome
    path qsrc_vcf

    output:
    tuple val(sample_id), file("${vcf.baseName}.recalibrated.vcf")

    script:
    def knownSitesArgsStr = knownSitesArgs.join(' ')
    def degradedDna = params.degraded_dna == "true"
    def isFreeBayes = params.variant_caller == 'freebayes'
    // For HaplotypeCaller, choose annotation set based on coverage
    def hcSnpAnnotations   = degradedDna ? '-an QD -an FS -an SOR'
                                         : '-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR'
    def hcIndelAnnotations = degradedDna ? '-an QD -an FS -an SOR'
                                         : '-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum'

    """
    echo "Running VQSR"

    if [[ -n "${params.genome_file}" ]]; then
        genomeFasta=\$(basename ${params.genome_file})
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi

    echo "Genome File: \${genomeFasta}"

    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    if ${isFreeBayes}; then
        echo "FreeBayes mode: hard-filtering SNPs and INDELs (VQSR not suitable for FreeBayes VCFs)"
        # FreeBayes does not emit GATK INFO annotations (QD, MQ, FS, SOR, etc.)
        # and its INDELs don't match GATK training resource coordinates.
        # Hard-filter both variant types using QUAL and DP thresholds.

        # --- Extract and hard-filter SNPs ---
        gatk SelectVariants \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --select-type-to-include SNP \
            -O ${vcf.baseName}.raw_snps.vcf
        gatk VariantFiltration \
            -R "\${genomeFasta}" \
            -V ${vcf.baseName}.raw_snps.vcf \
            --filter-expression "QUAL < 30.0" --filter-name "LowQual" \
            --filter-expression "DP < 5"      --filter-name "LowDP" \
            -O ${vcf.baseName}.output_SNP.vcf

        # --- Extract and hard-filter INDELs ---
        gatk SelectVariants \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --select-type-to-include INDEL \
            -O ${vcf.baseName}.raw_indels.vcf
        gatk VariantFiltration \
            -R "\${genomeFasta}" \
            -V ${vcf.baseName}.raw_indels.vcf \
            --filter-expression "QUAL < 30.0" --filter-name "LowQual" \
            --filter-expression "DP < 5"      --filter-name "LowDP" \
            -O ${vcf.baseName}.output_INDEL.vcf

    elif ${degradedDna}; then
        echo "Running VQSR for degraded DNA (1x coverage)"
        # relaxed parameters for SNPs
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            ${hcSnpAnnotations} \
            -mode SNP \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -O ${vcf.baseName}.recalibrated_SNP.recal
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -recal-file ${vcf.baseName}.recalibrated_SNP.recal \
            -mode SNP \
            -O ${vcf.baseName}.output_SNP.vcf
        # relaxed parameters for INDELs
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            ${hcIndelAnnotations} \
            -mode INDEL \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -O ${vcf.baseName}.recalibrated_INDEL.recal
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -recal-file ${vcf.baseName}.recalibrated_INDEL.recal \
            -mode INDEL \
            -O ${vcf.baseName}.output_INDEL.vcf
    else
        echo "Running VQSR for standard DNA (10x+ coverage)"
        # stricter parameters for SNPs
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            ${hcSnpAnnotations} \
            -mode SNP \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -O ${vcf.baseName}.recalibrated_SNP.recal
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -recal-file ${vcf.baseName}.recalibrated_SNP.recal \
            -mode SNP \
            -O ${vcf.baseName}.output_SNP.vcf
        # stricter parameters for INDELs
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            ${hcIndelAnnotations} \
            -mode INDEL \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -O ${vcf.baseName}.recalibrated_INDEL.recal
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -recal-file ${vcf.baseName}.recalibrated_INDEL.recal \
            -mode INDEL \
            -O ${vcf.baseName}.output_INDEL.vcf
    fi

    gatk MergeVcfs \
        -I ${vcf.baseName}.output_SNP.vcf \
        -I ${vcf.baseName}.output_INDEL.vcf \
        -O ${vcf.baseName}.recalibrated.vcf

    echo "VQSR Complete"
    """
}
