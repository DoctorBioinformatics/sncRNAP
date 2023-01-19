#!/usr/bin/env nextflow
/* 
sncRNAP is a pipeline that can be used for analyzing smallRNAseq datasets
More information can be found at https://github.com/sncRNAP
*/
////////////////////////////////////////////////////
//* --             Set defaults           -- *////
////////////////////////////////////////////////////
params.paired_samples = false
params.layout = null
params.input_dir = null
params.version = false
params.help = false
params.output_dir = "Results"
params.genome = 'human'
params.min_length = 16
params.max_length = 50
////////////////////////////////////////////////////
//* --    Print help, version messages     -- *////
////////////////////////////////////////////////////
def helpMessage() {
    log.info """\
    
    ===========
    sncRNAP
    ===========
    Usage: Group comparison analysis:
    nextflow run sncRNAP -profile conda --input_dir '*fq.gz' 
    --output_dir ./Results --genome GRCm38 --min_length 15 --trim_galore_max_length 50 
    --mature "https://mirbase.org/ftp/CURRENT/mature.fa.gz" 
    --hairpin "https://mirbase.org/ftp/CURRENT/hairpin.fa.gz" 
    --layout ./layout.csv

    Mandatory arguments:
    output_dir directory: ${params.output_dir}
    input_dir directory of files (/path/to/files): ${params.input_dir}
    Layout file: ${params.layout}
    
    Other arguments:
    genome (human/mouse/rat): ${params.genome}
    Minimum read length: ${params.min_read_length}
    Print help" ${params.help}
    """
}

// Print pipeline version
version = "Version:  sncRNAP 0.1"

// Print message for user
def versionMessage() {
    log.info """\
    ${version}
    """
}

// Show help message and quit
if (params.help) {
    helpMessage()
    versionMessage()
    exit 0
}

// Show version and quit
if (params.version) {
    versionMessage()
    exit 0
}
////////////////////////////////////////////////////
//* --         Validate parameters         -- *////
////////////////////////////////////////////////////
// input_dir parameter
if(!params.input_dir){
    exit 1, "Error: No input_dir provided. Provide --input_dir to pipeline"
}

// If output_dir is absolute path, do nothing. Otherwise, add launch dir working dir to path
if(params.output_dir.startsWith("/")){
    output_dir = "$params.output_dir"
} else {
    output_dir = "$launchDir/$params.output_dir"
}

// If layout is absolute path, do nothing. Otherwise, add launch dir working dir to path
if(params.layout){
    if(params.layout.startsWith("/")){
        layoutfile = "$params.layout"
    } else {
        layoutfile = "$launchDir/$params.layout"
    }
}

// Check if genome exists in the config file
if (!params.genome) {
    exit 1, "The provided genome '${params.genome}' is not available. Currently the available genomes are (human/mouse/rat)}"
}

//CPU usage set to 50%
CPU_usage = Math.round((Runtime.runtime.availableProcessors()*50)/100)
////////////////////////////////////////////////////
//* --         Channels and Paths          -- *////
////////////////////////////////////////////////////
// Create a channel for input_dir read files
reads_ch = Channel.fromPath( ["$params.input_dir/*.fastq.gz", "$params.input_dir/*.fq.gz"])

// Create a channel to read only a single fastq file from the input directory
import java.nio.file.*
all_files = Files.list(Paths.get("$params.input_dir")).toList()
random_index = Math.floor(Math.random() * all_files.size()).toInteger()
reads_single = all_files.get(random_index)

// Create a path for genomeFasta
params.genomeFasta = "$baseDir/DBs/${params.genome}_tRNAs-and-ncRNAs-and-lookalikes.fa_exc_miRNA.fa"

// Create a path for adapter sequence
params.adapter = "$baseDir/adapters_db/adapters.fa"

// Create a path for layout file
layout_channel = Channel.fromPath(["$params.layout"])
novel_layout_channel=Channel.fromPath(["$params.layout"])
////////////////////////////////////////////////////
//* --            Preprocessing             -- *////
////////////////////////////////////////////////////
/*
 * STEP 1 - Create genome_db
 */
process genome_db {
    cpus CPU_usage
    tag "$fasta"

    input:
    path fasta from params.genomeFasta

    output:
    path "star_db" into genome_db_ch

    script:
    def memory  = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    echo $memory
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_db/ \\
        --genomeFastaFiles $fasta \\
        --runThreadN $task.cpus \\
        --genomeSAsparseD 3 \\
        --genomeSAindexNbases 12 \\
        --genomeChrBinNbits 14 \\
        $memory
    STAR --version | sed -e "s/STAR_//g" > STAR.version.txt
    """
}

/*
 * STEP 2 - Get Adapter Sequence
 */
process get_adapter {
    cpus CPU_usage
    tag "$read"

    input:
    file read from reads_single
    path adapter from params.adapter

    output:
    file "adapter_sequence.txt" into adapter_sequence

    script: 
    """
    get_adapter.py $read $adapter  > adapter_sequence.txt
    """
}

/*
 * STEP 3 - Trim Galore
 */
process trim_galore {
    cpus CPU_usage
    tag "$reads"
    publishDir "${params.output_dir}/trim_galore", mode: 'copy'

    input:
    file reads from reads_ch
    file adapter from adapter_sequence

    output:
    file '*trimming_report.txt' into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
    file '*.gz' into trimmed_reads

    script: 
    """
    # read adapter file as a string and run trim_galore
    cat $adapter| while read -r F 
    do
        echo \$F
        trim_galore --adapter \$F \\
        --stringency 10 --length ${params.min_length} \\
        --max_length ${params.max_length} --gzip $reads --fastqc
    done
    """
}
////////////////////////////////////////////////////
//* --              Processes              -- *////
////////////////////////////////////////////////////
/*
 * STEP 4 - STAR Mapping
 */
process STAR_Mapping {
    cpus CPU_usage
    tag "$trimmed"

    input:
    file trimmed from trimmed_reads
    path star_db from genome_db_ch

    output:
    path("*Aligned.out.bam"), emit: bam
    path("*Log.final.out"), emit: log_final
    path("*_Stats.log"), emit: sam_stats

    script:
    """
    STAR \\
        --genomeDir ${star_db} \\
        --readFilesIn ${trimmed} \\
        --runThreadN ${task.cpus} \\
        --outSAMattributes AS nM HI NH \\
        --outFilterMultimapScoreRange 0 \\
        --outFilterMatchNmin ${params.min_length} \\
        --outFileNamePrefix ${trimmed}_ \\
        --readFilesCommand zcat
    grep 'Number of input reads' ${trimmed}_Log.final.out \\
        | sed -r 's/\\s+//g' \\
        | awk -F '|' '{print \$2}' \\
        > ${trimmed}_Stats.log
    echo ' reads; of these:' >> ${trimmed}_Stats.log
    samtools view -bS ${trimmed}_Aligned.out.sam > ${trimmed}_Aligned.out.bam
    rm ${trimmed}_Aligned.out.sam
    """
}

