# Test Data Directory

This directory contains test data files for the SRA2SAMRefiner test suite.

## Files

- `small_reference.fasta` - Minimal reference sequence for testing
- `test_reads.sam` - Sample SAM file with various read types
- `test_reads.bam` - Binary version of the SAM file
- `empty.sam` - Empty SAM file with header only

## Usage

These files are generated programmatically by the test fixtures in `conftest.py`.
They are small, focused datasets designed for unit and integration testing.

## Data Sources

All test data is synthetic and created specifically for testing purposes.
No real biological samples are used in the test suite.