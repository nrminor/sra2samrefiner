#!/usr/bin/env nextflow

workflow {

    assert params.accession_list :
    "No accession list text file provided with the `--accession_list` argument."
    assert file(params.accession_list).isFile() :
    "The provided accession list does not exist or cannot be read."
    assert params.ref_fasta :
    "No fasta reference file provided with the `--ref_fasta` argument."
    assert file(params.ref_fasta).isFile() :
    "The provided fasta reference file does not exist or cannot be read."
    assert params.ref_gbk :
    "No genbank reference file provided with the `--ref_gbk` argument."
    assert file(params.ref_gbk).isFile() :
    "The provided genbank reference file does not exist or cannot be read."

    // Channel for the SRA run accessions to download
    ch_sra_accession = Channel
        .fromPath(params.accession_list)
        .splitCsv(strip: true)
        .filter { it -> !it[0].startsWith("#") }
        .map { it -> it[0] }
        .unique()

    // The reference FASTA
    ch_ref_fasta = Channel.fromPath(params.ref_fasta)

    // The same reference in Genbank format
    ch_ref_gbk = Channel.fromPath(params.ref_gbk)


    FETCH_FASTQ(
        ch_sra_accession
    )

    // split the output FASTQs from the SRA download into two channels, where one
    // contains paired-end libraries with >=2 FASTQs, and the other contains the rest
    FETCH_FASTQ.out
        .branch { sample_id, fastq_files ->
            // Simple two-way branch: paired vs single
            // Check if we have 2 or more files (paired-end)
            paired: (fastq_files instanceof List || fastq_files instanceof Collection) && fastq_files.size() >= 2
                // Return all FASTQ files so Nextflow stages them in work dir
                return tuple(sample_id, fastq_files)
            
            // Everything else is treated as single-end
            single: true
                // For single-end, return the first (or only) FASTQ file
                def fastq_file = fastq_files instanceof List || fastq_files instanceof Tuple
                    ? fastq_files[0]
                    : fastq_files
                return tuple(sample_id, file(fastq_file))
        }
        .set { ch_sorted_fastqs }

    MERGE_PAIRS(
        ch_sorted_fastqs.paired
    )

    MAP_TO_REF(
        MERGE_PAIRS.out.mix(ch_sorted_fastqs.single).combine(ch_ref_fasta)
    )

    TRIM_ALIGNED_ENDS(
        MAP_TO_REF.out
            .filter { _id, count, _sam ->
                int count_int = count as Integer
                count_int  > 0 }
            .combine(ch_ref_fasta)
    )

    HANDLE_DUPLICATES(
        TRIM_ALIGNED_ENDS.out
        .filter { _id, count, _sam ->
            int count_int = count as Integer
            count_int  > 0 }
    )

    SAM_REFINER(
        HANDLE_DUPLICATES.out
        .filter { _id, count, _sam ->
            int count_int = count as Integer
            count_int  > 0 }
        .map { id, _count, sam -> tuple( id, file(sam) ) }
        .combine(ch_ref_gbk)
    )

    SORT_AND_CONVERT(
        HANDLE_DUPLICATES.out
        .map { id, _count, sam -> tuple( id, file(sam) ) }
        .combine(ch_ref_fasta)
    )
}

process FETCH_FASTQ {

    tag "${run_accession}"

    maxForks params.max_concurrent_downloads
    cpus 3

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    input:
    val run_accession

    output:
    tuple val(run_accession), path("${run_accession}*.fastq")

    script:
    """
    prefetch ${run_accession}
    fasterq-dump --split-files --skip-technical ${run_accession}
    """
}

process MERGE_PAIRS {

    tag "${run_accession}"
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(run_accession), path(fastqs)

    output:
    tuple val(run_accession), path("${run_accession}.merged.fastq.gz")

    script:
    """
    # Script auto-detects R1/R2 files from the accession name
    repair_merge_tag.sh \\
    -o ${run_accession} \\
    -t ${task.cpus}
	"""
}

process MAP_TO_REF {

    tag "${run_accession}"
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 4

    input:
    tuple val(run_accession), path(derep_fa_reads), path(ref_fasta)

    output:
    tuple val(run_accession), env('NUM_RECORDS'), path("${run_accession}.SARS2.wg.bam")

    script:
    """
    # run minimap2
    minimap2 -a ${ref_fasta} ${derep_fa_reads} \\
    --sam-hit-only --secondary=no \\
    | samtools view -@ ${task.cpus} -b -o ${run_accession}.SARS2.wg.bam

    # count the records in the SAM file
    NUM_RECORDS=\$(samtools view -@ ${task.cpus} -c ${run_accession}.SARS2.wg.bam)
    """
}

process TRIM_ALIGNED_ENDS {

    tag "${run_accession}"
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 1

    input:
    tuple val(run_accession), val(pre_trim_count), path(alignment), path(ref)

    output:
    tuple val(run_accession), env('NUM_RECORDS'), path("${run_accession}.SARS2.wg.trimmed.bam")

    script:
    """
    echo "Trimmimg ends on ${pre_trim_count} from ${alignment}..." >&2
    trim_aligned_reads.py \\
    --in ${alignment} \\
    --out "${run_accession}.SARS2.wg.trimmed.bam" \\
    --ref ${ref} \\
    --merged-left ${params.end_trim_bases} \\
    --merged-right ${params.end_trim_bases} \\
    --r1-left ${params.end_trim_bases} \\
    --r2-right ${params.end_trim_bases} \\
    --single-left ${params.end_trim_bases} \\
    --single-right ${params.end_trim_bases} \\
    --clipping-mode ${params.clipping_mode} \\
    ${params.strip_tags ? '--strip-tags' : ''} \\
    --min-len 20 \\
    -v

    # count the records in the SAM file
    NUM_RECORDS=\$(samtools view -c ${run_accession}.SARS2.wg.trimmed.bam)

    echo "End-trimming successful. \$NUM_RECORDS reads remain." >&2
    """
}

process HANDLE_DUPLICATES {

    tag "${run_accession}"
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 4

    input:
    tuple val(run_accession), val(pre_dedup_count), path(trimmed_cram)

    output:
    tuple val(run_accession), env('NUM_RECORDS'), path("${run_accession}.trimmed.deduped.bam")

    script:
    """
    set -euo pipefail

    samtools collate -@ ${task.cpus} -O -u ${trimmed_cram} \\
    | samtools fixmate -@ ${task.cpus} -m -u - -  \\
    | samtools sort -@ ${task.cpus} -u - \\
    | samtools markdup -@ ${task.cpus} -S -s -r --duplicate-count --include-fails - - \\
    | samtools view -h - \\
    | prefix_dup_count.awk \\
    | samtools view -h -b > ${run_accession}.trimmed.deduped.bam \\
    && seqkit bam -s ${run_accession}.trimmed.deduped.bam

    # count the records in the SAM file
    NUM_RECORDS=\$(samtools view -c ${run_accession}.trimmed.deduped.bam)

    echo "Deduplication successful. \$NUM_RECORDS reads remain." >&2
    """
}

process SAM_REFINER {

    tag "${run_accession}"
    publishDir params.results, mode: params.reporting_mode, overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus params.sam_refiner_procs

    input:
    tuple val(run_accession), path(sam), path(ref_gbk)

    output:
    tuple val(run_accession), path("*.tsv"), emit: tsv

    script:
    """
    SAM_Refiner -r ${ref_gbk} -S ${sam} \\
    --wgs 1 --collect 0 --seq 1 --indel 0 --covar 1 --max_covar 1 --max_dist 100 \\
    --AAcentered 0 --nt_call 1 --min_count 1 --min_samp_abund 0 --ntabund 0 \\
    --ntcover 1 --AAreport 1 --chim_rm 0 --deconv 0 --mp ${task.cpus}
    """
}

process SORT_AND_CONVERT {

    tag "${run_accession}"
    publishDir params.results, mode: params.reporting_mode, overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 8

    input:
    tuple val(run_accession), path(sam), path(ref_fasta)

    output:
    tuple val(run_accession), path("${run_accession}*.cram")

    script:
    """
    cat ${sam} \\
    | samtools sort \\
    | samtools view -T ${ref_fasta} -@${task.cpus} \\
    -o ${run_accession}.SARS2.wg.cram
    """
}
