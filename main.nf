#!/usr/bin/env nextflow
/*
========================================================================================
                         tsRNAP
========================================================================================
 tsRNAP Pipeline.
 #### Homepage / Documentation
 https://github.com/tsRNAP
----------------------------------------------------------------------------------------
*/

//log.info Headers.nf_core(workflow, params.monochrome_logs)
////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
// Print message for user
def helpMessage() {
    log.info """\
    
    ===========
    tsRNAP
    ===========
    Usage: Group comparison analysis:
    nextflow run tsRNAP -profile conda --input '*fq.gz' --outdir ./Results --genome GRCm38 --layout ./layout.csv

    Mandatory arguments:
    Input directory of files (/path/to/files): ${params.input}
    Output directory: ${params.output}
    Genome: ${params.genome}
    Layout: ${params.layout}
    
    Other arguments:
    Minimum read length: ${params.min_length}
    Trim_galore length: ${params.trim_galore_max_length}
    Print help" ${params.help}
    """
}

/*
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/smrnaseq --input '*_R{1,2}.fastq.gz' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}
*/

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
/*
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}
*/
////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// Check if genome exists in the config file
// if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
//     exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
// }

//CPU usage set to 50%
CPU_usage=Math.round((Runtime.runtime.availableProcessors()*50)/100)

// Genome options
params.bt_index = params.genome ? params.genomes[ params.genome ].bowtie ?: false : false
params.mirtrace_species = params.genome ? params.genomes[ params.genome ].mirtrace_species ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.mirna_gtf = params.mirtrace_species ? "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/${params.mirtrace_species}.gff3" : false

// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_adapter = params.three_prime_adapter
//protocol = params.protocol
// Presets
//if (params.protocol == "illumina"){
  //  clip_r1 = 0
    //three_prime_clip_r1 = 0
    //three_prime_adapter = "TGGAATTCTCGGGTGCCAAGG"
//} else if (params.protocol == "nextflex"){
  //  clip_r1 = 4
    //three_prime_clip_r1 = 4
    //three_prime_adapter = "TGGAATTCTCGGGTGCCAAGG"
//} else if (params.protocol == "qiaseq"){
  //  clip_r1 = 0
    //three_prime_clip_r1 = 0
    //three_prime_adapter = "AACTGTAGGCACCATCAAT"
//} else if (params.protocol == "cats"){
  //  clip_r1 = 3
    //three_prime_clip_r1 = 0
    // three_prime_adapter = "GATCGGAAGAGCACACGTCTG"
    //three_prime_adapter = "AAAAAAAA"
//} else {
    //custom protocol
  //  clip_r1 = params.clip_r1
    //three_prime_clip_r1 = params.three_prime_clip_r1
    //three_prime_adapter = params.three_prime_adapter
    //protocol = params.protocol
//}

// mirtrace protocol defaults to 'params.protocol' if not set
//if (!params.mirtrace_protocol){
  //  mirtrace_protocol = protocol
//} else {
  //  mirtrace_protocol = params.mirtrace_protocol
//}

if (params.mirna_gtf) {
    mirna_gtf = file(params.mirna_gtf, checkIfExists: true)
} else {
    mirna_gtf = false
}

// Validate inputs
if (params.skip_mirdeep){
    if (params.mature) { mature = file(params.mature, checkIfExists: true) } else { exit 1, "Mature file not found: ${params.mature}" }
    if (params.hairpin) { hairpin = file(params.hairpin, checkIfExists: true) } else { exit 1, "Hairpin file not found: ${params.hairpin}" }
    indices_mirdeep2 = Channel.empty()
    fasta = Channel.empty()
} else {
    if (params.references_parsed){
        fasta = file("$params.references_parsed/genome.fa", checkIfExists: true)
        hairpin = file("$params.references_parsed/hairpin.fa", checkIfExists: true)
        mature = file("$params.references_parsed/mature.fa", checkIfExists: true)
        indices_mirdeep2 = Channel.fromPath("$params.references_parsed/genome.*.ebwt", checkIfExists: true).ifEmpty { exit 1, "Reference parsed genome indices not found: ${references_parsed}"}
    } else {
        if (params.mature) { reference_mature = file(params.mature, checkIfExists: true) } else { exit 1, "Mature miRNA fasta file not found: ${params.mature}" }
        if (params.hairpin) { reference_hairpin = file(params.hairpin, checkIfExists: true) } else { exit 1, "Hairpin miRNA fasta file not found: ${params.hairpin}" }
        if (params.fasta) { reference_genome = file(params.fasta, checkIfExists: true) } else { exit 1, "Reference genome Fasta file not found: ${params.fasta}" }
    }
}

if( params.bt_index ){
    bt_indices = Channel.fromPath("${params.bt_index}*", checkIfExists: true).ifEmpty { exit 1, "Bowtie1 index directory not found: ${bt_dir}" }
} else if( params.bt_indices ){
    bt_indices = Channel.from(params.bt_indices).map{ file(it) }.toList()
} else {
    log.info "No Bowtie 1 index supplied - host reference genome analysis will be skipped."
}
if( !params.mirtrace_species ){
    exit 1, "Reference species for miRTrace is not defined."
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

/*
 * Create a channel for input read files
 */

if(params.input_paths){
    Channel
        .from(params.input_paths)
        .map { file(it) }
        .ifEmpty { exit 1, "params.input_paths was empty - no input files supplied" }
        .into { raw_reads_fastqc; raw_reads_trimgalore; raw_reads_mirtrace }
} else {
    Channel
        .fromPath( params.input )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}" }
        .into { raw_reads_fastqc; raw_reads_trimgalore; raw_reads_mirtrace }
}

/*
 * Create a path for adapter sequence
 */
params.adapter = "$baseDir/adapters_db/adapters.fa"

/*
 * Create a path for layout file
 */
layout_channel = Channel.fromPath(["$params.layout"])

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
//log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if(workflow.revision)          summary['Pipeline Release'] = workflow.revision
summary['Run Name']            = workflow.runName
summary['Input']               = params.input
summary['Genome']              = params.genome
summary['Min Trimmed Length']  = params.min_length
summary["Trim 5' R1"]          = clip_r1
summary["Trim 3' R1"]          = three_prime_clip_r1
summary['miRBase mature']      = params.mature
summary['miRBase hairpin']     = params.hairpin
summary['Reference Genome']    = params.fasta
if(params.bt_index)            summary['Bowtie Index for Ref'] = params.bt_index
summary['Save Reference']      = params.save_reference ? 'Yes' : 'No'
//summary['Protocol']            = params.protocol
//summary['Mirtrace Protocol']   = mirtrace_protocol
summary['miRTrace species']    = params.mirtrace_species
summary["3' adapter"]          = three_prime_adapter
summary['Output dir']          = params.outdir
summary['Launch dir']          = workflow.launchDir
summary['Working dir']         = workflow.workDir
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Script dir']          = workflow.projectDir
summary['Config Profile']      = (workflow.profile == 'standard' ? 'UPPMAX' : workflow.profile)
summary['Fasta Ref']        = params.fasta
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-smrnaseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/smrnaseq Workflow Summary'
    section_href: 'https://github.com/nf-core/smrnaseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }
/*
* Parse software version numbers
*/
process get_software_versions {
   publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
   saveAs: {filename ->
       if (filename.indexOf(".csv") > 0) filename
       else null
   }

   output:
   file 'software_versions_mqc.yaml' into software_versions_yaml
   file "software_versions.csv"

   script:
   java_mem = ''
   if(task.memory){
       tmem = task.memory.toBytes()
       java_mem = "-Xms${tmem} -Xmx${tmem}"
   }
   """
   export mirtracejar=\$(dirname \$(which mirtrace))
   echo $workflow.manifest.version > v_pipeline.txt
   echo $workflow.nextflow.version > v_nextflow.txt
   echo \$(R --version 2>&1) > v_R.txt
   fastqc --version > v_fastqc.txt
   trim_galore --version > v_trim_galore.txt
   bowtie --version > v_bowtie.txt
   samtools --version > v_samtools.txt
   htseq-count -h > v_htseq.txt
   fasta_formatter -h > v_fastx.txt
   java $java_mem -jar \$mirtracejar/mirtrace.jar --mirtrace-wrapper-name mirtrace --version > v_mirtrace.txt
   multiqc --version > v_multiqc.txt
   miRDeep2.pl -h > v_mirdeep2.txt

   scrape_software_versions.py > software_versions_mqc.yaml
   """
}

if (!params.references_parsed && !params.skip_mirdeep){
  /*
   * PREPROCESSING - Parse genome.fa and build Bowtie index for the reference genome
   */
  process bowtie_indices {
    cpus CPU_usage
    //label 'process_medium'
    publishDir path: { params.save_reference ? "${params.outdir}/references_parsed" : params.outdir },
               saveAs: { params.save_reference ? it : null }, mode: 'copy'

    input:
    file refgenome from reference_genome
    file mature from reference_mature
    file hairpin from reference_hairpin

    output:
    file 'genome.edited.fa' into fasta
    file 'genome.*.ebwt' into indices_mirdeep2
    file 'hairpin.fa' into hairpin
    file 'mature.fa' into mature

    script:
    """
    # Uncompress FASTA reference files if necessary
    MATURE="$mature"
    HAIRPIN="$hairpin"
    if [ \${MATURE: -3} == ".gz" ]; then
        gunzip -f \$MATURE
        MATURE=\${MATURE%%.gz}
    fi
    if [ \${HAIRPIN: -3} == ".gz" ]; then
        gunzip -f \$HAIRPIN
        HAIRPIN=\${HAIRPIN%%.gz}
    fi

    # Remove any special base characters from reference genome FASTA file
    sed '/^[^>]/s/[^ATGCatgc]/N/g' $refgenome > genome.edited.fa

    # Remove spaces from miRBase FASTA files
    sed -i 's, ,_,g' \$HAIRPIN
    sed -i 's, ,_,g' \$MATURE

    # Build bowtie index
    bowtie-build genome.edited.fa genome --threads ${task.cpus}
    """
  }
}
/*
 * PREPROCESSING - Build Bowtie index for mature and hairpin
 */
process make_bowtie_index {
    cpus CPU_usage
    //label 'process_medium'
    publishDir path: { params.save_reference ? "${params.outdir}/bowtie/reference" : params.outdir },
               saveAs: { params.save_reference ? it : null }, mode: 'copy'

    input:
    file mature from mature
    file hairpin from hairpin

    output:
    file 'mature_idx.*' into mature_index_bowtie
    file 'hairpin_idx.*' into hairpin_index_bowtie, hairpin_index_bowtie_2
    file 'hairpin_idx.fa' into hairpin_mirtop

    script:
    """
    seqkit grep -r --pattern \".*${params.mirtrace_species}-.*\" $mature > mature_sps.fa
    seqkit seq --rna2dna mature_sps.fa > mature_igenome.fa
    fasta_formatter -w 0 -i mature_igenome.fa -o mature_idx.fa
    # fasta_nucleotide_changer -d -i mature_igenome.fa -o mature_idx.fa
    bowtie-build mature_idx.fa mature_idx --threads ${task.cpus}

    seqkit grep -r --pattern \".*${params.mirtrace_species}-.*\" $hairpin > hairpin_sps.fa
    seqkit seq --rna2dna hairpin_sps.fa > hairpin_igenome.fa
    # fasta_nucleotide_changer -d -i hairpin_igenome.fa -o hairpin_idx.fa
    fasta_formatter -w 0 -i hairpin_igenome.fa -o hairpin_idx.fa
    bowtie-build hairpin_idx.fa hairpin_idx --threads ${task.cpus}
    """
}

/*
 * STEP 1 - FastQC
 */
process fastqc {
    //cpus CPU_usage
    //label 'process_low'
    tag "$reads"
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }
    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    file reads from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc --quiet --threads ${task.cpus} $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    cpus CPU_usage
    //label 'process_low'
    tag "$reads"
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    file reads from raw_reads_trimgalore
    path adapter from params.adapter
  
    output:
    file '*.gz' into trimmed_reads_bowtie, trimmed_reads_collapse, trimmed_reads_bowtie_ref, trimmed_reads_insertsize, trimmed_zipped_reads_mirdeep2
    file '*trimming_report.txt' into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
    file "*adapter_sequence.txt" into adapter_sequence

    script: 
    tg_length = "--length ${params.min_length}"
    c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    //tpa = (protocol == "qiaseq" | protocol == "cats") ? "--adapter ${three_prime_adapter}" : '--small_rna'
    """
    get_adapter.py $reads $adapter  > ${reads}.adapter_sequence.txt

    # get adapter as a string
    cat ${reads}.adapter_sequence.txt | while read -r F 
    do
        echo \$F
        trim_galore --adapter \$F $tg_length $c_r1 $tpc_r1 --max_length ${params.trim_galore_max_length} --gzip $reads --fastqc
    done
    """
}
/*
 * STEP 2.1 - Insertsize
 */
process insertsize {
    //label 'process_low'
    tag "$reads"
    publishDir "${params.outdir}/trim_galore/insertsize", mode: 'copy'

    input:
    file reads from trimmed_reads_insertsize

    output:
    file '*.insertsize' into insertsize_results

    script:
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    awk 'NR%4 == 2 {lengths[length(\$0)]++} END {for (l in lengths) {print l, lengths[l]}}' <(zcat $reads) >${prefix}.insertsize
    """
}

/*
 * STEP 2.2 - Collapse
 */
process collapse {
    //label 'process_medium'
    tag "$reads"

    input:
    file reads from trimmed_reads_collapse

    output:
    file 'final/*.fastq' into collapsed_fasta

    script:
    prefix = reads.toString() - '_trimmed.fq.gz'
    """
    seqcluster collapse -f $reads -m 1 --min_size 15 -o collapsed
    mkdir final
    mv collapsed/${prefix}_trimmed_trimmed.fastq final/${prefix}.fastq
    """
}


/*
 * STEP 3 - Bowtie miRBase mature miRNA
 */
process bowtie_miRBase_mature {
    cpus CPU_usage
    //label 'process_medium'
    tag "$reads"
    publishDir "${params.outdir}/bowtie/miRBase_mature", mode: params.publish_dir_mode, pattern: '*.mature_unmapped.fq.gz'

    input:
    file reads from trimmed_reads_bowtie
    file index from mature_index_bowtie

    output:
    file '*.mature.bam' into miRBase_mature_bam
    file '*.mature_unmapped.fq.gz' into mature_unmapped_reads

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    seq_center = params.seq_center ? "--sam-RG ID:${prefix} --sam-RG 'CN:${params.seq_center}'" : ''
    """
    bowtie \\
        $index_base \\
        -q <(zcat $reads) \\
        -p ${task.cpus} \\
        -t \\
        -k 50 \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        --un ${prefix}.mature_unmapped.fq \\
        -S $seq_center \\
        | samtools view -bS - > ${prefix}.mature.bam

    gzip ${prefix}.mature_unmapped.fq
    """
}

/*
 * STEP 4 - Bowtie against miRBase hairpin
 */
process bowtie_miRBase_hairpin {
    cpus CPU_usage
    //label 'process_medium'
    tag "$reads"
    publishDir "${params.outdir}/bowtie/miRBase_hairpin", mode: params.publish_dir_mode, pattern: '*.hairpin_unmapped.fq.gz'

    input:
    file reads from mature_unmapped_reads
    file index from hairpin_index_bowtie

    output:
    file '*.hairpin.bam' into miRBase_hairpin_bam, miRBase_hairpin_bam_mirtop
    file '*.hairpin_unmapped.fq.gz' into hairpin_unmapped_reads

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.toString() - '.mature_unmapped.fq.gz'
    seq_center = params.seq_center ? "--sam-RG ID:${prefix} --sam-RG 'CN:${params.seq_center}'" : ''
    """
    bowtie \\
        $index_base \\
        -p ${task.cpus} \\
        -t \\
        -a \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        -q <(zcat $reads) \\
        --un ${prefix}.hairpin_unmapped.fq \\
        -S $seq_center \\
        | samtools view -bS - > ${prefix}.hairpin.bam

    gzip ${prefix}.hairpin_unmapped.fq
    """
}

/*
 * STEP 4.1 - Bowtie against miRBase hairpin with collapsed reads
 */
process bowtie_miRBase_hairpin_collapsed {
    cpus CPU_usage
    //label 'process_medium'
    tag "$reads"

    input:
    file reads from collapsed_fasta
    file index from hairpin_index_bowtie_2

    output:
    file '*.bam' into miRBase_hairpin_collapse_bam

    script:
    index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    prefix = reads.baseName
    seq_center = params.seq_center ? "--sam-RG ID:${prefix} --sam-RG 'CN:${params.seq_center}'" : ''
    """
    bowtie \\
        $index_base \\
        -p ${task.cpus} \\
        -t \\
        -k 50 \\
        -a \\
        --best \\
        --strata \\
        -e 99999 \\
        --chunkmbs 2048 \\
        -q <(cat $reads) \\
        -S $seq_center \\
        | samtools view -bS - > ${prefix}.bam
    """
}

/*
 * STEP 5.1 - Post-alignment processing for miRBase mature and hairpin
 */
def wrap_mature_and_hairpin = { file ->
    if ( file.contains("mature") ) return "miRBase_mature/$file"
    if ( file.contains("hairpin") ) return "miRBase_hairpin/$file"
}

process mirna_post_alignment {
    //label 'process_medium'
    tag "$input"
    publishDir "${params.outdir}/bowtie", mode: params.publish_dir_mode, saveAs: wrap_mature_and_hairpin

    input:
    file input from miRBase_mature_bam.mix(miRBase_hairpin_bam)

    output:
    file "${input.baseName}.stats" into miRBase_counts
    file "*.{flagstat,idxstats,stats}" into ch_sort_bam_flagstat_mqc
    file "${input.baseName}.sorted.bam" into miRBase_bam
    file "${input.baseName}.sorted.bam.bai" into miRBase_bai

    script:
    """
    samtools sort ${input.baseName}.bam -o ${input.baseName}.sorted.bam
    samtools index ${input.baseName}.sorted.bam
    samtools idxstats ${input.baseName}.sorted.bam > ${input.baseName}.stats
    samtools flagstat ${input.baseName}.sorted.bam > ${input.baseName}.sorted.bam.flagstat
    samtools stats ${input.baseName}.sorted.bam > ${input.baseName}.sorted.bam.stats
    """
}


/*
 * STEP 5.2 - DESeq2 RNAseq count analysis 
 */
process DESEQ {
    publishDir "${params.outdir}/DeSEQ", mode: params.publish_dir_mode
    tag "$mature_counts"

    input:
    file input_files from miRBase_counts.toSortedList()
    path layout from layout_channel

    output:
    file '*.{pdf,csv,xlsx}' into DESEQ_results
    file "mature_counts.csv" into mature_counts_ch

    script:
    """
    DESeq2.r $layout $params.paired_samples $input_files 
    """
}

/*
 * STEP 5.3 - miRNA format conversion to mirGFF3
 */
process mirtop_bam_hairpin {
    //label 'process_medium'
    tag "$input"
    publishDir "${params.outdir}", mode: 'copy'

    when:
    mirna_gtf

    input:
    file input from miRBase_hairpin_collapse_bam.collect()
    file hairpin from hairpin_mirtop
    file gtf from mirna_gtf

    output:
    file "mirtop/mirtop.gff" into mirtop_gff
    file "mirtop/mirtop.tsv" into mirtop_tsv
    file "mirtop/mirna.tsv" into mirna_tsv
    file "mirtop/mirtop_rawData.tsv" into isomir_tsv

    script:
    """
    mirtop gff --hairpin $hairpin --gtf $gtf -o mirtop --sps $params.mirtrace_species $input
    mirtop counts --hairpin $hairpin --gtf $gtf -o mirtop --sps $params.mirtrace_species --add-extra --gff mirtop/mirtop.gff
    mirtop export --format isomir --hairpin $hairpin --gtf $gtf --sps $params.mirtrace_species -o mirtop mirtop/mirtop.gff
    collapse_mirtop.r mirtop/mirtop.tsv
    """
}


/*
 * STEP 6.1 and 6.2 IF A GENOME SPECIFIED ONLY!
 */
if( params.bt_index ) {

    /*
     * STEP 6.1 - Bowtie 1 against reference genome
     */
    process bowtie_ref {
        cpus CPU_usage
        //label 'process_high'
        tag "$reads"
        publishDir "${params.outdir}/bowtie_ref", mode: 'copy'

        input:
        file reads from trimmed_reads_bowtie_ref
        file indices from bt_indices.collect()

        output:
        file '*.genome.bam' into bowtie_bam, bowtie_bam_for_unmapped

        script:
        index_base = indices[0].toString() - ~/.rev.\d.ebwt?/ - ~/.\d.ebwt?/
        prefix = reads.toString() - ~/(.R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
        seq_center = params.seq_center ? "--sam-RG ID:${prefix} --sam-RG 'CN:${params.seq_center}'" : ''
        """
        bowtie \\
            $index_base \\
            -q <(zcat $reads) \\
            -p ${task.cpus} \\
            -t \\
            -k 50 \\
            --best \\
            --strata \\
            -e 99999 \\
            --chunkmbs 2048 \\
            -S $seq_center \\
            | samtools view -bS - > ${prefix}.genome.bam
        """
    }

    process genome_post_alignment  {
        //label 'process_low'
        tag "$input"
        publishDir "${params.outdir}/bowtie_ref", mode: 'copy'

        input:
        file input from bowtie_bam

        output:
        file "*.{flagstat,idxstats,stats}" into ch_genome_bam_flagstat_mqc

        script:
        """
        samtools sort ${input.baseName}.bam -o ${input.baseName}.sorted.bam
        samtools index ${input.baseName}.sorted.bam
        samtools idxstats ${input.baseName}.sorted.bam > ${input.baseName}.stats
        samtools flagstat ${input.baseName}.sorted.bam > ${input.baseName}.sorted.bam.flagstat
        samtools stats ${input.baseName}.sorted.bam > ${input.baseName}.sorted.bam.stats
        """
    }

    /*
     * STEP 6.2 - Statistics about unmapped reads against ref genome
     */

    process bowtie_unmapped {
        //label 'process_ignore'
        //label 'process_medium'
        tag "${input_files[0].baseName}"
        publishDir "${params.outdir}/bowtie_ref/unmapped", mode: 'copy'

        input:
        file input_files from bowtie_bam_for_unmapped.toSortedList()

        output:
        file 'unmapped_refgenome.txt' into bowtie_unmapped

        script:
        """
        for i in $input_files
        do
          printf "\${i}\t"
          samtools view -c -f0x4 \${i}
        done > unmapped_refgenome.txt
        """
    }

} else{
    ch_genome_bam_flagstat_mqc = Channel.empty()
}
/*
 * STEP 7.0 - unzip trim galore reads
 */

process uncompress_trimmed_reads{
    cpus CPU_usage
    //label 'process_low'

    when:
    !params.skip_mirdeep

    input:
    path reads from trimmed_zipped_reads_mirdeep2

    output:
    path "$unzip" into trimmed_reads_mirdeep2

    script:
    unzip = reads.toString() - '.gz'
    """
    pigz -f -d -p ${task.cpus} $reads
    """

}


/*
 * STEP 7.1 - miRDeep2 mapper
 */
process mapper_mirdeep2 {
    //label 'process_medium'
    tag "$reads"
    publishDir "${params.outdir}/mirdeep2/mapper", mode: 'copy'

    when:
    !params.skip_mirdeep

    input:
    file reads from trimmed_reads_mirdeep2
    file indices from indices_mirdeep2.collect()


    output:
    file '*_collapsed.fa' into mirdeep_reads_collapsed
    file '*reads_vs_refdb.arf' into reads_vs_refdb

    script:
    index_base = indices.toString().tokenize(' ')[0].tokenize('.')[0]

    """
    mapper.pl \\
    $reads \\
    -e \\
    -h \\
    -i \\
    -j \\
    -m \\
    -p $index_base \\
    -s ${reads.baseName}_collapsed.fa \\
    -t ${reads.baseName}_reads_vs_refdb.arf \\
    -o 4
    """
}

/*
 * STEP 7.2 - miRDeep2 novel mirnas
 */
process mirdeep2 {
    cpus CPU_usage
    //label 'process_medium'
    publishDir "${params.outdir}/mirdeep2/mirdeep", mode: 'copy'

    when:
    !params.skip_mirdeep

    input:
    file refgenome from fasta
    file reads_collapsed from mirdeep_reads_collapsed
    file reads_vs_refdb from reads_vs_refdb
    file mature from mature
    file hairpin from hairpin

    output:
    file 'result*.{bed,csv,html}'

    script:
    """
    perl -ane 's/[ybkmrsw]/N/ig;print;' $hairpin > hairpin_ok.fa
    sed 's/ .*//' $refgenome | awk '\$1 ~ /^>/ {gsub(/_/,"",\$1); print; next} {print}' > genome_nowhitespace.fa
    
    miRDeep2.pl \\
    $reads_collapsed \\
    genome_nowhitespace.fa \\
    $reads_vs_refdb \\
    $mature \\
    none \\
    hairpin_ok.fa \\
    -d \\
    -z _${reads_collapsed.simpleName}
    """
}


/*
 * STEP 8 - miRTrace
 */
process mirtrace {
    cpus CPU_usage
    //label 'process_medium'
    tag "$reads"
    publishDir "${params.outdir}/miRTrace", mode: 'copy'

    input:
    file reads from raw_reads_mirtrace.collect()
    path adapter from adapter_sequence.collect()

    output:
    file '*mirtrace' into mirtrace_results

    script:
    //primer = (protocol=="cats") ? " " : " --adapter $three_prime_adapter "
    //protocol_opt = (protocol=="custom") ? " " : " --protocol $protocol "
    java_mem = ''
    if(task.memory){
        tmem = task.memory.toBytes()
        java_mem = "-Xms${tmem} -Xmx${tmem}"
    }
    """
    export mirtracejar=\$(dirname \$(which mirtrace))
    for i in $reads
    do
        path=\$(realpath \${i})
        prefix=\$(echo \${i} | sed -e 's/.gz//' -e 's/.fastq//' -e 's/.fq//' -e 's/_val_1//' -e 's/_trimmed//' -e 's/_R1//' -e 's/.R1//')
        echo \$path","\$prefix
    done > mirtrace_config
    adapter=`cat $adapter | while read -r F; do echo \$F; break; done`

    java $java_mem -jar \$mirtracejar/mirtrace.jar --mirtrace-wrapper-name mirtrace qc \\
         --species $params.mirtrace_species \\
         --adapter \$adapter \\
         --config mirtrace_config \\
         --write-fasta \\
         --output-dir mirtrace \\
         --force
    """
}

/*
 * STEP 9 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

    when:
    !params.skip_qc && !params.skip_multiqc

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file ('fastqc/*') from fastqc_results.collect()
    file ('trim_galore/*') from trimgalore_results.collect()
    file ('mirtrace/*') from mirtrace_results.collect()
    file ('samtools/*') from ch_sort_bam_flagstat_mqc.collect()
    file ('samtools_genome/*') from ch_genome_bam_flagstat_mqc.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"

    script:
    rtitle = ''
    rfilename = ''
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        rtitle = "--title \"${workflow.runName}\""
        rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
    }
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc . -f $rtitle $rfilename $custom_config_file -m bowtie1 -m samtools -m cutadapt -m fastqc -m custom_content
    """
}


/*
 * STEP 10 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file 'results_description.html'

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}





/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/smrnaseq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/smrnaseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/smrnaseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/smrnaseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/smrnaseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/smrnaseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Switch the embedded MIME images with base64 encoded src
    smrnaseqlogo = new File("$baseDir/assets/smrnaseq_logo.png").bytes.encodeBase64().toString()
    email_html = email_html.replaceAll(~/cid:smrnaseqlogo/, "data:image/png;base64,$smrnaseqlogo")

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/smrnaseq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/smrnaseq]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "${c_red}====================================================${c_reset}\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "${c_red}====================================================${c_reset}\n"
                }
            }
        }
    }
}
