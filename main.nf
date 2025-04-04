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
        .splitText()

    // The reference FASTA
    ch_ref_fasta = Channel.fromPath(params.ref_fasta)

    // The same reference in Genbank format
    ch_ref_gbk = Channel.fromPath(params.ref_gbk)


    FETCH_FASTQ(
        ch_sra_accession.combine(ch_ref_fasta).combine(ch_ref_gbk)
    )

    // split the output FASTQs from the SRA download into two channels, where one
    // contains paired-end libraries with >=2 FASTQs, and the other contains the rest
    FETCH_FASTQ.out
        .branch { sample_id, fastq_files ->
            single: fastq_files.size() == 1
                return tuple(sample_id, file(fastq_files[0]))
            paired: fastq_files.size() > 1
                return tuple(sample_id, file(fastq_files[0]), file(fastq_files[1]))
            other: false
        }
        .set { ch_fastq_cardinality }

    MERGE_PAIRS(
        ch_fastq_cardinality.paired
    )

    DEREPLICATE_READS(
        MERGE_PAIRS.out.mix(ch_fastq_cardinality.single)
    )

    TRIM_ENDS(
        FETCH_FASTQ.out
    )

    MAP_TO_REF(
        TRIM_ENDS.out.combine(ch_ref_fasta)
    )

    SAM_REFINER(
        TRIM_ENDS.out.combine(ch_ref_gbk)
    )

    SORT_AND_CONVERT(
        SAM_REFINER.out.combine(ch_ref_fasta)
    )
}


process FETCH_FASTQ {

    tag "${run_accession}"
    publishDir params.results

    maxForks params.max_concurrent_downloads

    input:
    val run_accession

    output:
    tuple val(run_accession), path("*.fastq")

    script:
    """
    prefetch ${run_accession}
    fasterq-dump ${run_accession} --split-3
    """
}

process MERGE_PAIRS {

    tag "${run_accession}"
    publishDir params.results

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    time { "${5 * task.attempt}minutes" }
    cpus 4

    input:
    tuple val(run_accession), path(reads1), path(reads2)

    output:
    tuple val(run_accession), path("${run_accession}.merged.fastq.gz")

    script:
    def Xmx = 4 * task.attempt
    """
	bbmerge.sh \
	-Xmx${Xmx}g \
	in1=`realpath ${reads1}` \
	in2=`realpath ${reads2}` \
	out=${run_accession}.merged.fastq.gz \
	outu=${run_accession}.unmerged.fastq.gz \
    qtrim=t \
	ihist=${run_accession}_ihist_merge.txt \
	threads=${task.cpus} \
	-eoom
	"""
}

process DEREPLICATE_READS {

    tag "${run_accession}"
    publishDir params.results

    cpus 1

    input:
    tuple val(run_accession), path(fastq)

    output:
    tuple val(run_accession), path("${run_accession}.collapsed.fasta")

    script:
    """
    derep.py ${fastq} ${run_accession}.collapsed.fasta 1
    """
}

process TRIM_ENDS {

    tag "${run_accession}"
    publishDir params.results

    input:
    tuple val(run_accession), path(fasta)

    output:
    tuple val(run_accession), path("${run_accession}*trimmed.fasta")

    script:
    if (params.end_trim_bases == 0 || params.end_trim_bases == false || params.end_trim_bases == null) {
        """
        cp ${fasta} ${run_accession}.collapsed.untrimmed.fasta
        """
    }
    else {
        """
        cat ${fasta} \
        | seqkit amplicon -r ${params.end_trim_bases}:-${params.end_trim_bases} \
        -o ${run_accession}.collapsed.${params.end_trim_bases}trimmed.fasta
        """
    }
}

process MAP_TO_REF {

    tag "${run_accession}"
    publishDir params.results

    cpus 4

    input:
    tuple val(run_accession), path(derep_fa_reads), path(ref_fasta)

    output:
    tuple val(run_accession), path("${run_accession}.SARS2.wg.sam")

    script:
    """
    minimap2 -a ${ref_fasta} ${derep_fa_reads} \
    --sam-hit-only --secondary=no \
    -o ${run_accession}.SARS2.wg.sam
    """
}

process SAM_REFINER {

    tag "${run_accession}"
    publishDir params.results

    cpus 8

    input:
    tuple val(run_accession), path(sam), path(ref_gbk)

    output:
    tuple val(run_accession), path("*.sam"), emit: sam
    tuple val(run_accession), path("*.tsv"), emit: tsv

    script:
    """
    SAM_Refiner.py -r ${ref_gbk} \
    --wgs 1 --collect 0 --seq 1 --indel 0 --covar 1 --nt_call 1 --min_count 1 \
    --min_samp_abund 0 --ntabund 0 --ntcover 1 --mp ${task.cpus} --chim_rm 0 --deconv 0 \
    -S ${run_accession}.SARS2.wg.sam
    """
}

process SORT_AND_CONVERT {

    tag "${run_accession}"
    publishDir params.results

    cpus 8

    input:
    tuple val(run_accession), path(sam), path(ref_fasta)

    output:
    tuple val(run_accession), path("${run_accession}*.cram")

    script:
    """
    cat ${sam} \
    | samtools sort \
    | samtools view -T ${ref_fasta} -@${task.cpus} \
    -o ${run_accession}.SARS2.wg.cram
    """
}

