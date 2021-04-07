#!/usr/bin/env nextflow

/*
 * Default pipeline parameters are placed at the top of the
 * script for now but will be extracted out to the config
 * file once the pipeline is working. This will be helpful
 * for when sharing pipelines as users will only change
 * the config file to match their environment and should
 * not need to change the main.nf file.
*/

adaptSeq = params.adaptseq

println """\
         M A P P I N G - N F   P I P E L I N E
         =====================================
         ref	      : ${params.ref}
	 reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
        .stripIndent()

Channel
    .fromFilePairs ( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs2_ch}
/*
process download_ref {

     tag "GET reference genome"
     publishDir "${params.outdir}/ref_genome", mode: 'copy'
     
     output:
     file "*" into ref_ch, ref2_ch

     """
     wget -O ref.gz https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Ureaplasma_parvum/representative/GCF_000019345.1_ASM1934v1/GCF_000019345.1_ASM1934v1_genomic.fna.gz
     gunzip -c ref.gz > ref.fa
     """
}

process bwa_index {
    
    tag "INDEX reference genome"
    publishDir "${params.outdir}/index_genome", mode: 'copy'

    input:
    file "*" from ref_ch
   
    output:
    file "*" into index_ch, index2_ch

    """
    bwa index ref.fa
    """
}

process fastqc_raw {
    tag "FASTQC on $sample_id of the raw data"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
   
    input:
    set sample_id, file(reads) from read_pairs_ch

    output:
    set sample_id, file("*") into fastqc_raw_ch

    """
    fastqc \
    ${reads[0]} ${reads[1]}
    """
}

process trim_galore {
    tag "TRIM_GALORE on $sample_id using $adaptSeq adapter sequence for trim"
    publishDir "${params.outdir}/trim-galore", mode: 'copy'

    input:
    set sample_id, file(reads) from read_pairs2_ch

    output:
    set sample_id, file("*fq") into trim_ch, trim2_ch, trim3_ch

    script:

    """
    trim_galore \
    --paired ${reads[0]} ${reads[1]} \
    --$adaptSeq  
    """
}

process fastqc_trim {
    tag "FASTQC on $sample_id trimmed data"
    publishDir "${params.outdir}/fastqc_trim", mode: 'copy'

    input:
    set sample_id, file(reads_trim) from trim_ch

    output:
    set sample_id, file("*") into fastqc_trim_ch

    """
    fastqc \
    ${reads_trim[0]} ${reads_trim[1]}
    """
}

/*
*process screen_fastq {
*    tag "FASTQ_SCREEN on $sample_id of the trimmed data"
*    publishDir "${params.outdir}/fastq_screen_trim", mode: 'copy'
*
*    input:
*    set sample_id, file(reads_trim) from trim2_ch
*
*    output:
*    set sample_id, file("*.html") into fastq_screen_h_ch
*    set sample_id, file("*.txt") into fastq_screen_t_ch
*
*    """
*    fastq_screen \
*    --aligner 'bwa' \
*    --conf /scicomp/home-pure/ptx4/ureaplasma/NF/db/fastq_screen_genomes \
*    ${reads_trim[0]} ${reads_trim[1]}
*    """
*}
*/

process bwa_map {
    tag "BWA_MEM on $sample_id of the trimmed data"
    publishDir "${params.outdir}/sam_files", mode: 'copy'

    input:
    set sample_id, file(reads_trim) from trim3_ch
    file(ref) from ref_ch
    file(index) from index_ch

    output:
    set sample_id, file("*.sam") into sam_ch

    """
    bwa mem ref.fa ${reads_trim[0]} ${reads_trim[1]} > ${sample_id}.sam
    """
}

process sam_view {
    tag "SAMTOOLS view on $sample_id after running bwa-mem"
    publishDir "${params.outdir}/sam_view", mode: 'copy'

    input:
    set sample_id, file(sam_file) from sam_ch

    output:
    set sample_id, file("*.bam") into sam_view_ch

   """
   samtools view -b ${sam_file} > ${sample_id}.bam
   """
}

process sam_sort {
   tag "SAMTOOLS sort on $sample_id after running samtools view"
   publishDir "${params.outdir}/sam_sort", mode: 'copy'

   input:
   set sample_id, file(sam_file2) from sam_view_ch

   output:
   set sample_id, file("*.sorted.bam") into bam_ch

   """
   samtools sort ${sam_file2} -o ${sample_id}.sorted.bam
   """
}

process bcftools_mpileup {
    tag "BCFTOOLS mpileup on $sample_id mapped data"
    publishDir "${params.outdir}/bcf_files", mode: 'copy'

    input:
    set sample_id, file(bam_file) from bam_ch
    file(ref) from ref2_ch
    file(index) from index2_ch

    output:
    set sample_id, file("*.bcf") into bcf_ch

    """
    bcftools \
    mpileup \
    -f ref.fa \
    ${bam_file} \
    -I \
    -O u \
    -o  ${sample_id}.bcf
    """
}

process bcftools_call {
    tag "BCFTOOLS call on $sample_id"
    publishDir "${params.outdir}/vcf_files", mode: 'copy'

    input:
    set sample_id, file(bcf_file) from bcf_ch

    output:
    set sample_id, file("*.vcf") into vcf_ch

    """
    bcftools call -c ${bcf_file} --ploidy 1 -O v -o ${sample_id}.vcf
    """
}

process zip_index {
   tag "BGZIP & TABIX on ${sample_id}"
   publishDir "${params.outdir}/bgzip_files", mode: 'copy'

   input:
   set sample_id, file(vcf_file) from vcf_ch

   output:
   set sample_id, file("*vcf.gz"), \
                  file("*vcf.gz.tbi") into bg_ch

   """
   bgzip ${vcf_file}
   tabix ${vcf_file}.gz
   """ 
}
