#!/usr/bin/env nextflow
/* 
sncRNAP is a pipeline that can be used for analyzing smallRNAseq datasets
More information can be found at https://github.com/sncRNAP
*/

////////////////////////////////////////////////////
//* --             Set defaults           -- *////
////////////////////////////////////////////////////

// params.species = 'human'
// params.paired_samples=false
// params.all_plots = false
// params.layout = null
// params.input_dir = null
// params.output_dir = "Results"
// params.min_read_length = 16
// params.version = false
// params.help = false

////////////////////////////////////////////////////
//* --    Print help, version messages     -- *////
////////////////////////////////////////////////////

def helpMessage() {
    log.info """\
    
    ===========
    sncRNAP
    ===========
    Usage: Group comparison analysis:
    nextflow run sncRNAP -profile conda --input '*fq.gz' 
    --outdir ./Results --genome GRCm38 --min_length 15 --trim_galore_max_length 50 
    --mature "https://mirbase.org/ftp/CURRENT/mature.fa.gz" 
    --hairpin "https://mirbase.org/ftp/CURRENT/hairpin.fa.gz" 
    --layout ./layout.csv

    Mandatory arguments:
    Output directory: ${params.output_dir}
    Input directory of files (/path/to/files): ${params.input_dir}
    Layout file: ${params.layout}
    
    Other arguments:
    Species (human/mouse/rat): ${params.species}
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

// Input parameter
if(!params.input_dir){
    exit 1, "Error: No input provided. Provide --input_dir to pipeline"
}

// If output_dir is absolute path, do nothing. Otherwise, add launch dir working dir to path
if(params.output_dir.startsWith("/")){
    outputdir = "$params.output_dir"
} else {
    outputdir = "$launchDir/$params.output_dir"
}

// If layout is absolute path, do nothing. Otherwise, add launch dir working dir to path
if(params.layout.startsWith("/")){
    layoutfile = "$params.layout"
} else {
    layoutfile = "$launchDir/$params.layout"
}

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
}

//CPU usage set to 50%
CPU_usage=Math.round((Runtime.runtime.availableProcessors()*50)/100)

////////////////////////////////////////////////////
//* --         Channels and Paths          -- *////
////////////////////////////////////////////////////

// Create a channel for input read files
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

// Create a path for adapter sequence
params.adapter = "$baseDir/adapters_db/adapters.fa"

// Create a path for layout file
layout_channel = Channel.fromPath(["$params.layout"])
novel_layout_channel=Channel.fromPath(["$params.layout"])

////////////////////////////////////////////////////
//* --            Preprocessing             -- *////
////////////////////////////////////////////////////
