// Declare syntax version
nextflow.enable.dsl=2

process DEEPVARIANT {

    container = "${projectDir}/../singularity-images/google-deepvariant-1.4.0.img"

    input:
    path fasta
    path fasta_index
    path cram 
    path cram_index
    val prefix

    output:
    path "*.vcf.gz"
    path "*.vcf.gz.tbi"
    path "*.g.vcf.gz"
    path "*.g.vcf.gz.tbi"

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --ref=${fasta} \\
        --reads=${cram} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        --model_type WGS \\
        --num_shards=${params.thread_num}
    cp ${prefix}.vcf.gz ${launchDir}/${params.outdir}/
    cp ${prefix}.vcf.gz.tbi ${launchDir}/${params.outdir}/
    cp ${prefix}.g.vcf.gz ${launchDir}/${params.outdir}/
    cp ${prefix}.g.vcf.gz.tbi ${launchDir}/${params.outdir}/
    """
}

workflow{
    fasta       = Channel.of(params.fasta)
    fasta_index = Channel.of(params.fasta_index)
    cram        = Channel.of(params.cram)
    cram_index  = Channel.of(params.cram_index)
    prefix      = Channel.of(params.prefix)
    DEEPVARIANT(fasta, fasta_index, cram, cram_index, prefix)
}

