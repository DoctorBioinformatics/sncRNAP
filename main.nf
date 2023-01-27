#!/usr/bin/env nextflow
/* 
sncRNAP is a pipeline that can be used for analyzing smallRNAseq datasets
More information can be found at https://github.com/sncRNAP
*/
////////////////////////////////////////////////////
//* --             Set defaults           -- *////
////////////////////////////////////////////////////
params.output_dir = "Results"
params.min_length = 16
params.max_length = 50
params.genome = null
params.input_dir = null
params.mature = null
params.hairpin = null
params.paired_samples = false
params.version = false
params.help = false
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
 * STEP 1 - Create index for non_miRNA db
 */
process non_miRNA_db {
    cpus CPU_usage
    tag "$fasta"
    publishDir "${params.output_dir}/genome_ref", mode: 'copy'

    input:
    path fasta from params.genomeFasta

    output:
    path "non_miRNA_db" into non_miRNA_db_ch

    script:
    def memory  = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    echo $memory
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir non_miRNA_db/ \\
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
 * STEP 2 - Create index for mature miRNA
 */
process mature_idx {
    cpus CPU_usage
    publishDir "${params.output_dir}/mature_refs", mode: 'copy'

    input:
    path mature from params.mature

    output:
    path "mature_db" into mature_db_ch

    script:
    def memory  = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    # A) Uncompress FASTA files if files are zipped:
    MATURE="$mature"
    
    if [[ \$MATURE == *.gz ]]; then
        gunzip -f \$MATURE
        MATURE=\${MATURE%%.gz}
    fi

    # B) Remove spaces from files:
    sed -i 's, ,_,g' \$MATURE

    # C) Convert mRNA sequence into DNA:    
    seqkit grep -r --pattern ".*${params.genome == "human" || params.genome == "Human" ? "hsa" : params.genome == "mouse" || params.genome == "Mouse" ? "mmu" : params.genome == "rat" || params.genome == "Rat" ? "rno" : "unknown_species"}-.*" \$MATURE > mature_sps.fa
    seqkit seq --rna2dna mature_sps.fa > mature_genome.fa
    fasta_formatter -w 0 -i mature_genome.fa -o mature_idx.fa

    # D) Build genome idx for mature idx files:
    echo $memory
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir mature_db/ \\
        --genomeFastaFiles mature_idx.fa \\
        --runThreadN $task.cpus \\
        --genomeSAsparseD 3 \\
        --genomeSAindexNbases 6 \\
        --genomeChrBinNbits 14 \\
        $memory
    STAR --version | sed -e "s/STAR_//g" > STAR.version.txt
    """
}

/*
 * STEP 3 - Get Adapter Sequence
 */
process get_Adapter {
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
 * STEP 4 - Trim Galore
 */
process trim_Galore {
    cpus CPU_usage
    tag "$reads"
    publishDir "${params.output_dir}/trim_galore", mode: 'copy'

    input:
    file reads from reads_ch
    file adapter from adapter_sequence

    output:
    file '*trimming_report.txt' into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
    file '*.gz' into trimmed_reads_collapsed

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

/*
 * STEP 5 - read_Collapse
 */
process read_collapse {
    cpus CPU_usage
    tag "$reads"

    input:
    file reads from trimmed_reads_collapsed

    output:
    file 'final/*.fastq' into collapsed_non_miRNA_fasta, collapsed_miRNA_fasta

    script:
    prefix = reads.toString().replace("_trimmed.fq.gz", "")
    """
    seqcluster collapse -f $reads -m 1 --min_size 15 -o collapsed
    mkdir final
    mv collapsed/${prefix}_trimmed_trimmed.fastq final/${prefix}.fastq
    """
}
////////////////////////////////////////////////////
//* --              Processes              -- *////
////////////////////////////////////////////////////
/*
 * STEP 6 - STAR non_miRNA_Mapping
 */
process star_non_miRNA {
    cpus CPU_usage
    tag "$reads"

    input:
    file reads from collapsed_non_miRNA_fasta
    path non_miRNA_db from non_miRNA_db_ch

    output:
    file "*Aligned.out.bam" into star_non_miRNA_bam
    file "*Log.final.out" into star_non_miRNA_log_final
    file "*_Stats.log" into star_non_miRNA_sam_stats
    
    script:
    """
    STAR \\
        --genomeDir ${non_miRNA_db} \\
        --readFilesIn ${reads} \\
        --runThreadN ${task.cpus} \\
        --outSAMattributes AS nM HI NH \\
        --outFilterMultimapScoreRange 0 \\
        --outFilterMatchNmin ${params.min_length} \\
        --outFileNamePrefix ${reads}_ \\
    grep 'Number of input reads' ${reads}_Log.final.out \\
        | sed -r 's/\\s+//g' \\
        | awk -F '|' '{print \$2}' \\
        > ${reads}_Stats.log
    echo ' reads; of these:' >> ${reads}_Stats.log
    samtools view -bS ${reads}_Aligned.out.sam > ${reads}_Aligned.out.bam
    rm ${reads}_Aligned.out.sam
    """
}

/*
 * STEP 7 - Processing non_miRNA reads
 */
process non_mirna_processing {
    tag "$input"
    publishDir "${params.output_dir}/processed_reads/non_miRNA", mode: 'copy'

    input:
    file input from star_non_miRNA_bam

    output:
    file "${input.baseName}.stats" into non_mirna_counts
    file "*.{flagstat,idxstats,stats}" into non_mirna_ch_sort_bam_flagstat_mqc
    file "${input.baseName}.sorted.bam" into non_mirna_sorted_bam
    file "${input.baseName}.sorted.bam.bai" into non_mirna_bai

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
 * STEP 8 - DESeq2 RNAseq count analysis 
 */
process DESEQ {
    tag "$input_files"
    publishDir "${params.output_dir}/DESeq", mode: 'copy'
    
    input:
    file input_files from non_mirna_counts.toSortedList()
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
 * STEP 9 - miRNA format conversion to mirGFF3
 */
// process mirtop_bam_hairpin {
//     //label 'process_medium'
//     tag "$input"
//     publishDir "${params.outdir}", mode: 'copy'

//     when:
//     mirna_gtf

//     input:
//     file input from miRBase_hairpin_collapse_bam.collect()
//     file hairpin from hairpin_mirtop
//     file gtf from mirna_gtf

//     output:
//     file "mirtop/mirtop.gff" into mirtop_gff
//     file "mirtop/mirtop.tsv" into mirtop_tsv
//     file "mirtop/mirna.tsv" into mirna_tsv
//     file "mirtop/mirtop_rawData.tsv" into isomir_tsv

//     script:
//     """
//     mirtop gff --hairpin $hairpin --gtf $gtf -o mirtop --sps $params.mirtrace_species $input
//     mirtop counts --hairpin $hairpin --gtf $gtf -o mirtop --sps $params.mirtrace_species --add-extra --gff mirtop/mirtop.gff
//     mirtop export --format isomir --hairpin $hairpin --gtf $gtf --sps $params.mirtrace_species -o mirtop mirtop/mirtop.gff
//     collapse_mirtop.r mirtop/mirtop.tsv
//     """
// }

// /*
//  * STEP 7 - STAR mature_miRNA_Mapping
//  */
// process star_Mature {
//     cpus CPU_usage
//     tag "$reads"

//     input:
//     file reads from collapsed_miRNA_fasta
//     path mature_db from mature_db_ch

//     output:
//     path("*Aligned.out.bam"), emit: bam
//     path("*Log.final.out"), emit: log_final
//     path("*_Stats.log"), emit: sam_stats

//     script:
//     """
//     STAR \\
//         --genomeDir ${mature_db} \\
//         --readFilesIn ${reads} \\
//         --runThreadN ${task.cpus} \\
//         --outSAMattributes AS nM HI NH \\
//         --outFilterMultimapScoreRange 0 \\
//         --outFilterMatchNmin ${params.min_length} \\
//         --outFileNamePrefix ${reads}_ \\
//     grep 'Number of input reads' ${reads}_Log.final.out \\
//         | sed -r 's/\\s+//g' \\
//         | awk -F '|' '{print \$2}' \\
//         > ${reads}_Stats.log
//     echo ' reads; of these:' >> ${reads}_Stats.log
//     samtools view -bS ${reads}_Aligned.out.sam > ${reads}_Aligned.out.bam
//     rm ${reads}_Aligned.out.sam
//     """
// }

