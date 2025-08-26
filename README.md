# SRA2SAMRefiner

A Nextflow pipeline for automated downloading, quality control, and mutational
haplotype analysis of NCBI SRA datasets using
[SAMRefiner](https://github.com/degregory/SAM_Refiner). The pipeline enables
high-throughput processing of viral genomic data, with particular focus on
SARS-CoV-2 wastewater surveillance and variant tracking.

## About SAMRefiner

SAMRefiner is a computational method designed for analyzing high-throughput
sequencing data to track viral populations, particularly in wastewater samples.
It processes SAM-formatted sequencing files from amplicon sequencing data and
generates detailed variant reports with chimera removal capabilities.

### Key Features

- **Variant Analysis**: Generates comprehensive reports including unique
  sequences, nucleotide calls, insertions/deletions, and covariant sequences
- **Chimera Removal**: Two independent algorithms to eliminate PCR-generated
  artifacts
- **Flexible Applications**: Suitable for monitoring viral variants, pathogenic
  bacteria virulence factors, and human disease alleles
- **Environmental Monitoring**: Optimized for wastewater surveillance
  applications

## Prerequisites

### System Requirements

- Unix-like operating system (Linux, macOS)
- At least 8GB RAM (16GB+ recommended for large datasets)
- ~10GB free disk space per sample

### Software Dependencies

This pipeline uses [Pixi](https://pixi.sh) for dependency management, which
handles both conda and PyPI packages automatically.

#### Installing Pixi

```bash
# Linux/macOS
curl -fsSL https://pixi.sh/install.sh | bash

# Windows (PowerShell)
powershell -ExecutionPolicy ByPass -c "irm -useb https://pixi.sh/install.ps1 | iex"

# Restart your shell or run:
source ~/.bashrc  # or ~/.zshrc
```

## Quick Start

1. **Clone and setup the environment**:

```bash
git clone https://github.com/nrminor/sra2samrefiner.git
cd sra2samrefiner

# Install all dependencies (conda and PyPI packages)
pixi install

# Activate the environment
pixi shell
```

2. **Prepare your accession list**:

```bash
# Create a file with one SRA accession per line
echo -e 'SRR32970560\nSRR32970561' > accessions.txt
```

3. **Run the pipeline**:

```bash
# Basic run with SARS-CoV-2 reference
nextflow run . \
  --ref_fasta ./assets/sars2.fasta \
  --ref_gbk ./assets/sars2.gbk \
  --accession_list accessions.txt \
  --end_trim_bases 30
```

## Configuration Options

### Required Parameters

- `--accession_list`: Text file with one SRA accession per line
- `--ref_fasta`: Reference genome in FASTA format
- `--ref_gbk`: Same reference in GenBank format

### Optional Parameters

- `--max_concurrent_downloads`: Maximum concurrent SRA downloads (default: 3)
- `--sam_refiner_procs`: Number of processes for SAMRefiner (default: 4)
- `--end_trim_bases`: Bases to trim from both ends of reads (default: 0)
- `--results`: Output directory (default: `./results`)
- `--no_cleanup`: Disable automatic cleanup of work directory

### Example with Custom Parameters

```bash
nextflow run . \
  --ref_fasta my_reference.fasta \
  --ref_gbk my_reference.gbk \
  --accession_list my_samples.txt \
  --max_concurrent_downloads 5 \
  --sam_refiner_procs 8 \
  --end_trim_bases 20 \
  --results ./my_results
```

## Execution Profiles

The pipeline supports multiple execution environments:

### Local Execution (Default)

```bash
nextflow run . -profile standard [parameters...]
```

### Container-based Execution

```bash
# Using Docker
nextflow run . -profile docker [parameters...]

# Using Singularity/Apptainer
nextflow run . -profile apptainer [parameters...]
```

### High-Performance Computing

```bash
# CHTC/SLURM clusters
nextflow run . -profile chtc_hpc [parameters...]
```

## Output Files

Results are organized by sample in the `results/` directory:

```
results/
├── {sample_id}.{reference}.wg.cram           # Aligned sequences
├── {sample_id}.{reference}.wg_covars.tsv     # Covariant sequences
├── {sample_id}.{reference}.wg_nt_calls.tsv   # Nucleotide calls
├── {sample_id}.{reference}.wg_unique_seqs.tsv # Unique sequences
└── workflow-visualization.png                # Pipeline diagram
```

### Output File Descriptions

- **CRAM files**: Compressed aligned sequences
- **covars.tsv**: Covariant sequence analysis for variant detection
- **nt_calls.tsv**: Per-position nucleotide call frequencies
- **unique_seqs.tsv**: Dereplicated unique sequences with abundance

## Pipeline Workflow

1. **SRA Download**: Fetches FASTQ files from NCBI SRA
2. **Quality Control**: Validates and processes sequencing data
3. **Read Trimming**: Optional end-trimming of low-quality bases
4. **Alignment**: Maps reads to reference using minimap2
5. **SAMRefiner Analysis**: Variant calling and chimera removal
6. **Report Generation**: Creates summary reports and visualizations

## Troubleshooting

### Common Issues

1. **SRA Download Failures**:
   - Check internet connection and NCBI SRA availability
   - Verify accession numbers are correct
   - Reduce `--max_concurrent_downloads` if hitting rate limits

2. **Memory Issues**:
   - Reduce `--sam_refiner_procs` for large datasets
   - Ensure sufficient RAM is available

3. **Pixi Environment Issues**:
   ```bash
   # Reset environment
   pixi clean
   pixi install
   ```

### Getting Help

- Check the [Nextflow documentation](https://nextflow.io/docs/latest/)
- Review [SAMRefiner manuscript](https://www.mdpi.com/1999-4915/13/8/1647)
- Open an issue on GitHub

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open
issues for bugs and feature requests.

## License

This project is licensed under the terms specified in the LICENSE file.
