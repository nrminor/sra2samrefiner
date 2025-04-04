# Nextflow wrapper for downloading SRA runs and running SAMRefiner on them

This is a simple Nextflow wrapper meant to parallelize a SRA run downloads, QC, and SAMRefiner executions on an arbitrary number of SRA run accessions. By default, three run accessions at most will be downloaded at once, though this number can be customized with `--max_concurrent_downloads`.

To use without any additional orchestrator support, run the following (this assumes you have [Pixi]() installed):

```bash
# set up conda and PyPI dependencies with Pixi
pixi shell --frozen

# put together a list of accessions to download
echo 'SRR32970560\nSRR32970561' > accessions.txt

# run with SARS2 reference while trimming 30 bases from each end of each read
nextflow run . \
--ref_fasta ./assets/sars2.fasta \
--ref_gbk ./assets/sars2.gbk \
--accession_list accessions.txt \
--end_trim_bases 30
```

