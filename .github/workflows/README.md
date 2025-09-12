# GitHub Actions CI/CD for SRA2SAMRefiner

This directory contains GitHub Actions workflows that provide comprehensive continuous integration for the SRA2SAMRefiner bioinformatics pipeline.

## Workflows

### 🔄 CI Workflow (`ci.yml`)

The main CI workflow runs on push/PR to main and develop branches, providing comprehensive testing across multiple dimensions:

#### Jobs Overview

1. **Pixi Environment Build** (`pixi-environment`)
   - Tests pixi environment builds on Ubuntu and macOS
   - Verifies all bioinformatics tools install correctly (samtools, minimap2, nextflow, seqkit)
   - Confirms Python dependencies import successfully
   - Uses pixi caching for faster builds

2. **UV Dependency Resolution** (`uv-dependencies`)
   - Tests that PyPI dependencies resolve cleanly with UV
   - Verifies dry-run dependency compilation
   - Creates UV virtual environment and tests imports
   - Ensures compatibility between pixi and UV package management

3. **Python Quality Checks** (`python-quality`)
   - Runs Ruff linting and formatting checks
   - Executes pytest test suite with coverage reporting
   - Uploads coverage to Codecov (requires `CODECOV_TOKEN` secret)
   - Tests type checking if basedpyright is available

4. **Pipeline Integration Test** (`pipeline-integration`)
   - End-to-end test of the complete Nextflow pipeline
   - Uses minimal sample data (first accession only for CI speed)
   - Tests both dry-run and actual execution
   - Verifies pipeline produces expected outputs
   - Gracefully handles network/SRA download issues
   - Uploads pipeline artifacts for debugging

5. **Python Script Tests** (`python-scripts`)
   - Tests standalone Python scripts (`test_slicing_bug.py`, `test_cigar_logic.py`)
   - Verifies module imports and CLI interface
   - Ensures Python code works outside of nextflow context

#### Configuration Features

- **Concurrency Control**: Cancels in-progress runs for the same PR
- **Matrix Testing**: Tests on Ubuntu and macOS where applicable
- **Caching**: Uses pixi caching to speed up dependency installation
- **Timeout Protection**: 60-minute timeout for integration tests
- **Artifact Upload**: Preserves test results and logs for 7 days

### 🛡️ Branch Protection Workflow (`branch-protection.yml`)

Provides a single required status check for branch protection rules:

- **Purpose**: GitHub branch protection requires specific named checks
- **Strategy**: Creates a single "All Checks Passed" job that depends on CI completion
- **Benefits**: Simplifies branch protection rules and provides clear merge status

## Repository Setup Instructions

### 1. Enable GitHub Actions
Ensure Actions are enabled in your repository settings.

### 2. Configure Branch Protection (Recommended)

Navigate to **Settings → Branches → Add rule** for `main` branch:

```
Branch name pattern: main
✅ Require status checks to pass before merging
✅ Require branches to be up to date before merging
Required status checks:
  - All Checks Passed
✅ Require review from code owners  
✅ Dismiss stale reviews when new commits are pushed
✅ Restrict pushes that create new files
```

### 3. Add Repository Secrets (Optional)

For enhanced functionality, add these secrets in **Settings → Secrets and variables → Actions**:

- `CODECOV_TOKEN`: For coverage reporting (get from codecov.io)

### 4. Asset Requirements

The pipeline tests require these files in the `assets/` directory:
- `test-accessions.txt`: SRA accession numbers (one per line)
- `sars2.fasta`: SARS-CoV-2 reference genome in FASTA format
- `sars2.gbk`: SARS-CoV-2 reference genome in GenBank format

## Workflow Behavior

### On Pull Requests
- All CI jobs run in parallel for maximum speed
- Integration test uses minimal data (1 accession) for fast feedback
- Results are reported as PR status checks
- Branch protection prevents merge until all checks pass

### On Push to Main
- Full CI validation runs
- Pixi cache is updated for future runs
- Results validate the main branch remains healthy

### Failure Handling
- **Network Issues**: Pipeline tests gracefully handle SRA download failures
- **Tool Missing**: Clear error messages if bioinformatics tools aren't available
- **Python Errors**: Comprehensive error reporting with coverage information

## Local Development

To run similar checks locally:

```bash
# Environment checks
pixi install --locked
pixi run samtools --version

# Python quality
pixi run ruff check .
pixi run pytest -v

# UV dependency check  
uv sync --dry-run

# Pipeline dry-run
pixi run nextflow run main.nf --accession_list assets/ci-test-accessions.txt --ref_fasta assets/sars2.fasta --ref_gbk assets/sars2.gbk -profile containerless -preview
```

## Troubleshooting

### Common Issues

1. **Pixi Install Fails**
   - Check `pixi.lock` is committed and up-to-date
   - Verify platform compatibility in `pyproject.toml`

2. **Pipeline Test Fails**
   - Network issues with SRA downloads are expected in CI
   - Check asset files exist and are readable

3. **Python Tests Fail**
   - Verify test dependencies in `pyproject.toml` dev group
   - Check test files don't have import path issues

### Debugging

- Check **Actions** tab for detailed logs
- Download pipeline artifacts to inspect results locally
- Review `.nextflow.log` in artifacts for pipeline issues

## Maintenance

### Regular Updates Needed

- Update action versions (e.g., `actions/checkout@v4` → `@v5`)
- Review timeout values based on actual runtime
- Update Python versions in matrix testing
- Monitor and update pinned tool versions

### Performance Optimization

The workflows are optimized for speed:
- Parallel job execution
- Aggressive caching 
- Minimal test data for integration tests
- Early failure detection

This setup provides comprehensive quality assurance while maintaining fast feedback cycles essential for productive development.