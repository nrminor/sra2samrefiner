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

            // If there's just one FASTQ in an array, "unpack" it from the array returned by the glob
            single: fastq_files.size() == 1
                return tuple(sample_id, file(fastq_files[0]))

            // If there are two FASTQs, expect that the alphanumeric first will end with ".1.fastq.gz" and the second with ".2.fastq.gz",
            // which is a (mostly) reliable SRA convention.  Include empty file for merging process
            paired: fastq_files.size() == 2 && file(fastq_files[0]).getName().endsWith("1.fastq.gz") && file(fastq_files[1]).getName().endsWith("2.fastq.gz")
                return tuple(sample_id, file(fastq_files[0]), file(fastq_files[1]), file("assets/empty"))

            // There are a couple common cases of >2 FASTQs per accession that we can handle. The first is where the first two files
            // end with "1.fastq.gz" and ".2.fastq.gz" and the third ends with "3.fastq.gz". Assuming correct alphanumeric sorting, we handle
            // that in this branch.
            triple1: fastq_files.size() > 2 && file(fastq_files[0]).getName().endsWith("1.fastq.gz") && file(fastq_files[1]).getName().endsWith("2.fastq.gz")
                return tuple(sample_id, file(fastq_files[0]), file(fastq_files[1]), file(fastq_files[2]))

            // It's also possible that the third, non-R1/R2 reads are in a FASTQ that doesn't have a numbered suffix, e.g.,
            // SRR33146255.fastq.gz. When that's the case, that third FASTQ will end up first when the files are sorted by name.
            // We handle that case in this branch by indexing out the second and third FASTQ (groovy/nextflow are 0-indexed)
            triple2: fastq_files.size() > 2 && file(fastq_files[1]).getName().endsWith("1.fastq.gz") && file(fastq_files[2]).getName().endsWith("2.fastq.gz")
                return tuple(sample_id, file(fastq_files[1]), file(fastq_files[2]), file(fastq_files[0]))

            // If there's just one FASTQ and it isn't in an array.  Last so isFile() on array doesn't throw an error.
            single1: fastq_files.isFile()
                return tuple(sample_id, file(fastq_files))

            // Other cases are as-yet unsupported
            other: true
        }
        .set { ch_sorted_fastqs }
	// write out any fastqs that aren't sorted for troubleshooting code
	ch_sorted_fastqs.other.view()

    MERGE_PAIRS(
        ch_sorted_fastqs.paired
        .mix(
            ch_sorted_fastqs.triple1,
            ch_sorted_fastqs.triple2
        )
    )

    MAP_MERGED_TO_REF(
        MERGE_PAIRS.out.merged_fq
		.mix(
			ch_sorted_fastqs.single,
			ch_sorted_fastqs.single1,
		).combine(ch_ref_fasta)
    )

    MAP_UNMERGED_TO_REF(
        MERGE_PAIRS.out.unmerged_fq.combine(ch_ref_fasta)
    )

    // sort the joined outputs from the mapping to insure two files are passed to the trimming,
	// one each for merged and unmerged reads, with empty file taking the place a missing output
    MAP_MERGED_TO_REF.out.join(MAP_UNMERGED_TO_REF.out, remainder: true)
    .branch { sample_id, sam_file1, sam_file2 ->

            // If both files are present
            both: sam_file1 && sam_file2
                return tuple(sample_id, file(sam_file1), file(sam_file2))

            // If only the merged file is present, pass empty in place of unmerged
            merged: sam_file1.isFile()
                return tuple(sample_id, file(sam_file1), file("assets/empty"))

            // If only the unmerged file is present, pass empty in place of merged
            unmerged: sam_file2.isFile()
                return tuple(sample_id, file("assets/empty"), file(sam_file2))

            // Other cases are as-yet unsupported, which shouldn't exist
            other: true
        }
        .set { ch_sorted_sams }

	// write out any sams that aren't sorted for troubleshooting code
	ch_sorted_sams.other.view()

	TRIM_ENDS(
        ch_sorted_sams.both
		.mix(
			ch_sorted_sams.merged,
			ch_sorted_sams.unmerged
		)
    )

    SAM_REFINER(
        TRIM_ENDS.out
        .filter { _id, count, _sam ->
            int count_int = count as Integer
            count_int  > 0 }
        .map { id, _count, sam -> tuple( id, file(sam) ) }
        .combine(ch_ref_gbk)
    )

}

process FETCH_FASTQ {

    tag "${run_accession}"
    maxForks params.max_concurrent_downloads
    cpus 3

    maxRetries 2
    errorStrategy {
					task.exitStatus == 3 ? 'ignore' : //  ignore when no files for sample
					task.attempt <= maxRetries ? 'retry' : 'ignore'
					}

    input:
    val run_accession

    output:
    tuple val(run_accession), path("${run_accession}*.fastq.gz")

    script:
    """
    prefetch ${run_accession}
    fasterq-dump --split-3 ${run_accession}
	gzip ${run_accession}*.fastq
    """
}

process MERGE_PAIRS {

    tag "${run_accession}"
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    maxRetries 2
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }

    cpus 4

    input:
    tuple val(run_accession), path(reads1), path(reads2), path(reads3)

    output:
    tuple val(run_accession), path("${run_accession}.merged.fastq.gz"), emit: merged_fq
    tuple val(run_accession), path("${run_accession}.unmerged.fastq.gz"), emit: unmerged_fq

    script:
    """
	bbmerge.sh \
	in1=`realpath ${reads1}` \
	in2=`realpath ${reads2}` \
	out=${run_accession}.merged.fastq.gz \
	outu=${run_accession}.unmerged.fastq.gz \
    qtrim=t \
	ihist=${run_accession}_ihist_merge.txt \
	threads=${task.cpus} \
	-eoom

	# if orphaned reads file, cat to unmerged
	cat ${reads3} >> ${run_accession}.unmerged.fastq.gz

	"""
}

process MAP_MERGED_TO_REF {

    tag "${run_accession}"
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    maxRetries 2
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }

    cpus 4

    input:
    tuple val(run_accession), path(merged_file), path(ref_fasta)

    output:
    tuple val(run_accession), path("${run_accession}.merged.cram")

    script:
    """
    # run minimap2
    minimap2 -a ${ref_fasta} ${merged_file} \
    --sam-hit-only --secondary=no \
    | samtools view -T ${ref_fasta} -@${task.cpus} \
	-o ${run_accession}.merged.cram
    """
}

process MAP_UNMERGED_TO_REF {

    tag "${run_accession}"
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    maxRetries 2
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }

    cpus 4

    input:
    tuple val(run_accession), path(unmerged_file), path(ref_fasta)

    output:
    tuple val(run_accession), path("${run_accession}.unmerged.cram")

    script:
    """
    # run minimap2
    minimap2 -a ${ref_fasta} ${unmerged_file} \
    --sam-hit-only --secondary=no \
    | samtools view -T ${ref_fasta} -@${task.cpus} \
    -o ${run_accession}.unmerged.cram
    """
}

process  TRIM_ENDS {

    tag "${run_accession}"
    publishDir params.results, mode: params.reporting_mode, overwrite: true

    maxRetries 2
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }

    cpus 1

    input:
    tuple val(run_accession), path(merged_sam), path(unmerged_sam)

    output:
    tuple val(run_accession), env('NUM_RECORDS'), path("${run_accession}.cut.cram")

    script:
    """
    drtrimsam.py --in_file ${merged_sam} --in_file2 ${unmerged_sam} --out_file ${run_accession}.cut.cram --trim ${params.end_trim_bases}

	# count the records in the cram file
    NUM_RECORDS=\$(cat ${run_accession}.cut.cram.count)
    """
}


process SAM_REFINER {

    tag "${run_accession}"
    publishDir params.results, mode: params.reporting_mode, overwrite: true

    maxRetries 2
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }

    cpus params.sam_refiner_procs

    input:
    tuple val(run_accession), path(cram), path(ref_gbk)

    output:
    tuple val(run_accession), path("${run_accession}*.tsv.gz")

    script:
    """
    SAM_Refiner -r ${ref_gbk} -S ${cram} \
    --wgs 1 --collect 0 --seq 1 --indel 0 --covar 1 --max_covar 1 --max_dist 100 \
    --nt_call 1 --min_count 1 --min_samp_abund 0 --ntabund 0 \
    --ntcover 1 --AAreport 1 --chim_rm 0 --deconv 0 --mp ${task.cpus} &&
	gzip ${run_accession}*.tsv
    """
}

