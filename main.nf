#!/usr/bin/env nextflow

/*
 * Default pipeline parameters are placed at the top of the
 * script for now but will be extracted out to the config
 * file once the pipeline is working. This will be helpful
 * for when sharing pipelines as users will only change
 * the config file to match their environment and should
 * not need to change the main.nf file.
*/

params.ref = "$baseDir/ref/u.parvumRef.fa"
params.reads = "$baseDir/data/*_R{1,2}_001.fastq.gz"
params.outdir = "$HOME/results"
params.thread = 1

println """\
         G E N O M E - A S S E M B L Y  - N F 
         ====================================
         ref	      : ${params.ref}
	 reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
        .stripIndent()

Channel
    .fromFilePairs ( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { reads_ch; reads2_ch }

process bwa_index {
    
    tag "INDEX reference genome"
    publishDir "${params.outdir}/index_genome", mode: 'copy'

    input:
    file ref_file from file(params.ref)
   
    output:
    file "reference.*" into index_ch, index2_ch

    """
    bwa index ${ref_file} -p reference 
    """
}

process fastqc_raw {
    tag "FASTQC on raw $sample_id reads"
    publishDir "${params.outdir}/raw_fastqc", mode: 'copy' 

    input:
    set sample_id, file(reads) from reads_ch

    output:
    file "*_fastqc.{zip,html}" into fastqc_raw_ch

    script:
    """
    fastqc -q ${reads[0]} ${reads[1]}
    """
}

process trim_galore {

    tag "TRIM_GALORE on $sample_id"
    publishDir "${params.outdir}/trim_files", mode: 'copy'

    input:
    set sample_id, file(reads) from reads2_ch
    
    output: 
    set sample_id, file("*.fq.gz") into trim_ch, trim2_ch, trim3_ch

    """
    trim_galore --paired ${reads[0]} ${reads[1]}
    """
}

process fastqc_trim {
    tag "FASTQC on trimmed $sample_id reads"
    publishDir "${params.outdir}/trim_fastqc", mode: 'copy'

    input:
    set sample_id, file(reads) from trim_ch

    output:
    file "*_fastqc.{zip,html}" into fastqc_trim_ch

    script:
    """
    fastqc -q ${reads[0]} ${reads[1]}
    """
}

process fastq_screen {

    tag "FASTQ_SCREEN on trimmed reads of $sample_id"
    publishDir "${params.outdir}/screen", mode: 'copy'

    input:
    set sample_id, file(reads_trimmed) from trim2_ch
 
    output:
    set sample_id, file("*.html") into screenHTMl_ch
    set sample_id, file("*.txt")  into screenTXT_ch

    """
    fastq_screen \
    --conf /scicomp/home-pure/ptx4/db/fastq_screen_genomes/fastq_screen.conf \
    --tag --filter 000000 \
    ${reads_trimmed[0]} ${reads_trimmed[1]}
    """
}

process unicycler {
    tag "UNICYCLER on trimmed reads of $sample_id"
    publishDir "${params.outdir}/unicycler", mode: 'copy'

    input:
    set sample_id, file(reads) from trim3_ch

    output:
    set sample_id, file("${sample_id}_assembly.fasta") into (quast_ch, prokka_ch)
    file("${sample_id}_assembly.gfa")
    file("${sample_id}_unicycler.log")
    file("unicycler.version.txt") into unicycler_ver_ch
    
    script:
    """
    unicycler -1 ${reads[0]} -2 ${reads[1]} --keep 0 -o .
    mv unicycler.log ${sample_id}_unicycler.log
    # rename so that quast can use the name 
     mv assembly.gfa ${sample_id}_assembly.gfa
    mv assembly.fasta ${sample_id}_assembly.fasta
    unicycler --version | sed -e "s/Unicycler v//g" > unicycler.version.txt
    """
}

/*
*process bwa_map {
* 
*    tag "BWA_MEM on $sample_id of data"
*    publishDir "${params.outdir}/sam_files", mode: 'copy'
* 
*    cpus params.thread
*
*    input:
*    set sample_id, file(reads) from reads2_ch
*    file(reference) from index_ch
*
*    output:
*    set sample_id, file("*.sam") into sam_ch
*
*    """
*    bwa mem reference ${reads[0]} ${reads[1]} -t ${params.thread} > ${sample_id}.sam
*    """
*}
*/
