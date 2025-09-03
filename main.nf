#!/usr/bin/env nextflow

workflow {

    // Validate required parameters
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
    
    // Validate dual-mode configuration
    if (params.retro_data && !params.retro_data.isEmpty()) {
        // Auto-enable retrospective mode if retro_data provided
        params.mode = 'retrospective'
        assert file(params.retro_data).isDirectory() :
        "The provided retro_data path '${params.retro_data}' is not a valid directory."
        log.info "Auto-enabled retrospective mode with historical data: ${params.retro_data}"
    }
    
    if (params.mode == 'retrospective' && (!params.retro_data || params.retro_data.isEmpty())) {
        error "Retrospective mode requires --retro_data parameter with path to historical SAMRefiner outputs directory"
    }
    
    log.info "Running in ${params.mode} mode"

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

            // If there's just one FASTQ, "unpack" it from the array returned by the glob
            single: fastq_files.size() == 1
                return tuple(sample_id, file(fastq_files[0]))

            // If there are two FASTQs, expect that the alphanumeric first will end with ".0.fq.gz" and the second with ".1.fq.gz",
            // which is a (mostly) reliable SRA convention
            paired: fastq_files.size() > 1 && file(fastq_files[0]).getName().endsWith("0.fq.gz") && file(fastq_files[1]).getName().endsWith("1.fq.gz")
                return tuple(sample_id, file(fastq_files[0]), file(fastq_files[1]))

            // There are a couple common cases of >2 FASTQs per accession that we can handle. The first is where the first two files
            // end with "0.fq.gz" and ".1.fq.gz" and the third ends with "3.fastq". Assuming correct alphanumeric sorting, we handle
            // that in this branch.
            triple1: fastq_files.size() > 2 && file(fastq_files[0]).getName().endsWith("0.fq.gz") && file(fastq_files[1]).getName().endsWith("1.fq.gz")
                return tuple(sample_id, file(fastq_files[0]), file(fastq_files[1]))

            // It's also possible that the third, non-R1/R2 reads are in a FASTQ that doesn't have a numbered suffix, e.g.,
            // SRR33146255.fastq. When that's the case, that third FASTQ will end up first when the files are sorted by name.
            // We handle that case in this branch by indexing out the second and third FASTQ (groovy/nextflow are 0-indexed)
            triple2: fastq_files.size() > 2 && file(fastq_files[1]).getName().endsWith("0.fq.gz") && file(fastq_files[2]).getName().endsWith("1.fq.gz")
                return tuple(sample_id, file(fastq_files[1]), file(fastq_files[2]))

            // Other cases are as-yet unsupported
            other: true
        }
        .set { ch_sorted_fastqs }

    MERGE_PAIRS(
        ch_sorted_fastqs.paired
        .mix(
            ch_sorted_fastqs.triple1,
            ch_sorted_fastqs.triple2
        )
    )

    DEREPLICATE_READS(
        MERGE_PAIRS.out.mix(ch_sorted_fastqs.single)
    )

    TRIM_ENDS(
        DEREPLICATE_READS.out
    )

    MAP_TO_REF(
        TRIM_ENDS.out.combine(ch_ref_fasta)
    )

    SAM_REFINER(
        MAP_TO_REF.out
        .filter { _id, count, _sam ->
            int count_int = count as Integer
            count_int  > 0 }
        .map { id, _count, sam -> tuple( id, file(sam) ) }
        .combine(ch_ref_gbk)
    )

    SORT_AND_CONVERT(
        MAP_TO_REF.out
        .map { id, _count, sam -> tuple( id, file(sam) ) }
        .combine(ch_ref_fasta)
    )

    // Create channels for surveillance parameters (will be empty if not provided)
    ch_start_date = Channel.value(params.surveillance_start_date ?: '')
    ch_end_date = Channel.value(params.surveillance_end_date ?: '')
    ch_target_lineage = Channel.value(params.surveillance_target_lineage ?: 'JN.1')
    
    // Download UShER database once and cache it
    DOWNLOAD_USHER_DB()
    
    // Prepare reference data for surveillance analysis using cached UShER data
    PREPARE_REFERENCE_DATA(
        ch_start_date,
        ch_end_date,
        ch_target_lineage,
        DOWNLOAD_USHER_DB.out.usher_db,
        DOWNLOAD_USHER_DB.out.lineage_definitions
    )
    
    // Handle dual-mode data collection
    if (params.mode == 'retrospective') {
        
        // Collect current run SAMRefiner outputs
        ch_current_samples = SAM_REFINER.out.tsv
            .map { _run_id, tsv_files -> tsv_files }
            .flatten()
        
        // Collect historical SAMRefiner outputs from retro_data directory
        ch_historical_samples = Channel
            .fromPath("${params.retro_data}/**/*_{covars,nt_calls}.tsv*")
            .ifEmpty { log.warn "No historical SAMRefiner files found in ${params.retro_data}" }
        
        // Mix current + historical for complete dataset
        ch_all_samrefiner = ch_current_samples
            .mix(ch_historical_samples)
            .collect()
            
        log.info "Retrospective mode: Combining current samples with historical data from ${params.retro_data}"
        
        // Run complete temporal analysis with historical context
        ANALYZE_TEMPORAL_VARIANTS(
            ch_all_samrefiner,
            PREPARE_REFERENCE_DATA.out.sample_metadata,
            PREPARE_REFERENCE_DATA.out.lineage_mutations,
            PREPARE_REFERENCE_DATA.out.variant_surveillance
        )
        
        // Generate surveillance reports
        GENERATE_SURVEILLANCE_REPORTS(
            ANALYZE_TEMPORAL_VARIANTS.out.temporal_variants,
            ANALYZE_TEMPORAL_VARIANTS.out.significant_changes,
            PREPARE_REFERENCE_DATA.out.sample_metadata,
            PREPARE_REFERENCE_DATA.out.variant_surveillance
        )
        
    } else {
        log.info "Daily mode: Processing current samples only, skipping surveillance analysis"
    }
}

process FETCH_FASTQ {

    tag "${run_accession}"

    maxForks params.max_concurrent_downloads
    cpus 8

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    input:
    val run_accession

    output:
    tuple val(run_accession), path("${run_accession}*.f*q*")

    script:
    """
    xsra prefetch ${run_accession}
    xsra dump \
    --prefix ${run_accession}_ --skip-technical --split \
    --compression g --threads ${task.cpus} --outdir . \
    ${run_accession}
    """
}

process MERGE_PAIRS {

    tag "${run_accession}"
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(run_accession), path(reads1), path(reads2)

    output:
    tuple val(run_accession), path("${run_accession}.merged.fastq")

    script:
    """
	bbmerge.sh \
	in1=`realpath ${reads1}` \
	in2=`realpath ${reads2}` \
	out=${run_accession}.merged.fastq \
	outu=${run_accession}.unmerged.fastq \
    qtrim=t \
	ihist=${run_accession}_ihist_merge.txt \
	threads=${task.cpus} \
	-eoom
	
    # Concat the unmerged reads to the end of the merged reads file
    cat ${run_accession}.unmerged.fastq >> ${run_accession}.merged.fastq
	"""
}

process DEREPLICATE_READS {

    tag "${run_accession}"
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

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
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 3

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
        | seqkit subseq -r ${params.end_trim_bases}:-${params.end_trim_bases} \
        -o ${run_accession}.collapsed.${params.end_trim_bases}trimmed.fasta
        """
    }
}

process MAP_TO_REF {

    tag "${run_accession}"
    // publishDir params.results, mode: params.reporting_mode, overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }

    cpus 4

    input:
    tuple val(run_accession), path(derep_fa_reads), path(ref_fasta)

    output:
    tuple val(run_accession), env('NUM_RECORDS'), path("${run_accession}.SARS2.wg.sam")

    script:
    """
    # run minimap2
    minimap2 -a ${ref_fasta} ${derep_fa_reads} \
    --sam-hit-only --secondary=no \
    -o ${run_accession}.SARS2.wg.sam

    # count the records in the SAM file
    NUM_RECORDS=\$(samtools view -c ${run_accession}.SARS2.wg.sam)
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
    SAM_Refiner -r ${ref_gbk} -S ${sam} \
    --wgs 1 --collect 0 --seq 1 --indel 0 --covar 1 --max_covar 1 --max_dist 100 \
    --AAcentered 0 --nt_call 1 --min_count 1 --min_samp_abund 0 --ntabund 0 \
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
    cat ${sam} \
    | samtools sort \
    | samtools view -T ${ref_fasta} -@${task.cpus} \
    -o ${run_accession}.SARS2.wg.cram
    """
}

// Split expensive operations for optimal caching
process DOWNLOAD_USHER_DB {
    
    tag "usher_db"
    
    // Cache UShER database permanently - only re-download if file doesn't exist
    storeDir "${params.results}/reference_cache/usher"
    
    cache 'deep' // Cache based on process definition, not inputs
    
    cpus 1
    memory '2 GB'
    
    output:
    path "public-latest.all.masked.pb.gz", emit: usher_db
    path "LineageDefinitions.tsv", emit: lineage_definitions
    
    when:
    params.surveillance_start_date && params.surveillance_end_date
    
    script:
    """
    # Download UShER database (large file, changes infrequently)
    curl -L -o public-latest.all.masked.pb.gz \\
        http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz
    
    # Extract lineage definitions
    matUtils extract -i public-latest.all.masked.pb.gz -C LineageDefinitions.tsv
    """
}

process PREPARE_REFERENCE_DATA {
    
    tag "reference_data_${start_date}_${end_date}"
    
    // Cache results based on date range and parameters
    storeDir "${params.results}/reference_cache/metadata"
    publishDir params.results, mode: params.reporting_mode, overwrite: true
    
    cache 'standard' // Cache based on inputs - will re-run if date range changes
    
    cpus 2
    memory '4 GB'
    
    input:
    val start_date
    val end_date
    val target_lineage
    path usher_db
    path lineage_definitions
    
    output:
    path "lineage_mutations_${target_lineage}_${start_date}_${end_date}.tsv", emit: lineage_mutations, optional: true
    path "variant_surveillance_${start_date}_${end_date}.tsv", emit: variant_surveillance, optional: true  
    path "sample_metadata_${start_date}_${end_date}.tsv", emit: sample_metadata, optional: true
    
    when:
    start_date && end_date && start_date != '' && end_date != ''
    
    script:
    """
    prepare_reference_data.py \\
        --start-date ${start_date} \\
        --end-date ${end_date} \\
        --target-lineage ${target_lineage} \\
        --usher-db ${usher_db} \\
        --lineage-definitions ${lineage_definitions} \\
        --lineage-mutations lineage_mutations_${target_lineage}_${start_date}_${end_date}.tsv \\
        --variant-surveillance variant_surveillance_${start_date}_${end_date}.tsv \\
        --sample-metadata sample_metadata_${start_date}_${end_date}.tsv \\
        --exclude-bioprojects ${params.surveillance_exclude_bioprojects ?: 'PRJNA748354 PRJEB44932'} \\
        --exclude-patterns ${params.surveillance_exclude_patterns ?: 'location:Maryland.*78365'} \\
        --min-reads ${params.surveillance_min_reads ?: 0}
    """
}

process ANALYZE_TEMPORAL_VARIANTS {
    
    tag "temporal_analysis_batch"
    
    publishDir params.results, mode: params.reporting_mode, overwrite: true
    
    cpus 8
    memory '16 GB'
    
    input:
    path "samrefiner_outputs/*"  // All SAMRefiner outputs as single directory
    path sample_metadata
    path lineage_mutations
    path variant_surveillance
    
    output:
    path "temporal_variants.tsv", emit: temporal_variants, optional: true
    path "weekly_summary.tsv", emit: weekly_summary, optional: true
    path "significant_changes.tsv", emit: significant_changes, optional: true
    
    when:
    params.surveillance_start_date && params.surveillance_end_date
    
    script:
    """
    analyze_temporal_variants.py \\
        --input-dir samrefiner_outputs \\
        --sample-metadata ${sample_metadata} \\
        --lineage-mutations ${lineage_mutations} \\
        --variant-surveillance ${variant_surveillance} \\
        --delta-threshold ${params.surveillance_delta_threshold ?: 0.02} \\
        --count-threshold ${params.surveillance_count_threshold ?: 100} \\
        --window-weeks ${params.surveillance_window_weeks ?: 3} \\
        --output-temporal temporal_variants.tsv \\
        --output-weekly weekly_summary.tsv \\
        --output-significant significant_changes.tsv
    """
}

process GENERATE_SURVEILLANCE_REPORTS {
    
    tag "surveillance_reports"
    
    publishDir params.results, mode: params.reporting_mode, overwrite: true
    
    cpus 4
    memory '8 GB'
    
    input:
    path temporal_variants
    path significant_changes
    path sample_metadata
    path variant_surveillance
    
    output:
    path "lineage_definitions.json", emit: json_lineages, optional: true
    path "regional_analysis.tsv", emit: regional_analysis, optional: true
    path "surveillance_summary.json", emit: summary_stats, optional: true
    
    when:
    params.surveillance_start_date && params.surveillance_end_date
    
    script:
    """
    generate_surveillance_reports.py \\
        --temporal-variants ${temporal_variants} \\
        --significant-changes ${significant_changes} \\
        --sample-metadata ${sample_metadata} \\
        --variant-surveillance ${variant_surveillance} \\
        --output-json lineage_definitions.json \\
        --output-regional regional_analysis.tsv \\
        --output-summary surveillance_summary.json
    """
}

