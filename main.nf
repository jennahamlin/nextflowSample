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
params.reads = "$baseDir/data/*_R{1,2}_001.fastq"
params.outdir = "$HOME/results"
params.thread = 1

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

process bwa_index {
    
    tag "INDEX reference genome"
//    publishDir "${params.outdir}/index_genome", mode: 'copy'

    input:
    file ref_file from file(params.ref)
   
    output:
    file "reference.*" into index_ch, index2_ch

    """
    bwa index ${ref_file} -p reference 
    """
}

process bwa_map {
    tag "BWA_MEM on $sample_id of data"
    publishDir "${params.outdir}/sam_files", mode: 'copy'
 
    cpus params.thread

    input:
    set sample_id, file(reads) from read_pairs_ch
    file(reference) from index_ch

    output:
    set sample_id, file("*.sam") into sam_ch

    """
    bwa mem reference ${reads[0]} ${reads[1]} -t ${params.thread} > ${sample_id}.sam
    """
}
