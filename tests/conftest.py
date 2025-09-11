# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pysam",
#     "pytest",
# ]
# ///
"""
Pytest fixtures and configuration for SRA2SAMRefiner testing.

This module provides shared fixtures for testing trim_aligned_reads.py and related
bioinformatics tools. It includes utilities for creating mock BAM files, SAM records,
and test data management.
"""

import sys
import tempfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pysam
import pytest

# Add bin directory to Python path so we can import the modules under test
BIN_DIR = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))

# Now we can import the modules we're testing
from trim_aligned_reads import TrimPolicy


@pytest.fixture
def temp_dir() -> Iterator[Path]:
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture
def reference_sequence() -> str:
    """Simple reference sequence for testing alignment."""
    return "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"


@pytest.fixture
def reference_fasta(temp_dir: Path, reference_sequence: str) -> Path:
    """Create a simple reference FASTA file for testing."""
    ref_path = temp_dir / "reference.fasta"
    with open(ref_path, "w") as f:
        f.write(">test_reference\n")
        f.write(f"{reference_sequence}\n")
    return ref_path


@pytest.fixture
def default_trim_policy() -> TrimPolicy:
    """Default trimming policy for testing."""
    return TrimPolicy(
        merged_left=10,
        merged_right=10,
        r1_left=10,
        r2_right=10,
        single_left=5,
        single_right=5,
        min_len=20,
    )


@pytest.fixture
def aggressive_trim_policy() -> TrimPolicy:
    """Aggressive trimming policy for edge case testing."""
    return TrimPolicy(
        merged_left=30,
        merged_right=30,
        r1_left=25,
        r2_right=25,
        single_left=20,
        single_right=20,
        min_len=10,
    )


@pytest.fixture
def minimal_trim_policy() -> TrimPolicy:
    """Minimal trimming policy for testing."""
    return TrimPolicy(
        merged_left=1,
        merged_right=1,
        r1_left=1,
        r2_right=1,
        single_left=1,
        single_right=1,
        min_len=5,
    )


class MockAlignedSegment:
    """Mock AlignedSegment for unit testing without pysam dependency."""

    def __init__(
        self,
        query_name: str = "test_read",
        query_sequence: str = "ATCGATCGATCG",
        query_qualities: list[int] | None = None,
        cigartuples: list[tuple[int, int]] | None = None,
        reference_start: int = 0,
        is_reverse: bool = False,
        is_unmapped: bool = False,
        is_secondary: bool = False,
        is_supplementary: bool = False,
    ) -> None:
        self.query_name = query_name
        self.query_sequence = query_sequence
        self.query_qualities = query_qualities or (
            [30] * len(query_sequence) if query_sequence else None
        )
        self.cigartuples = cigartuples
        self.reference_start = reference_start
        self.is_reverse = is_reverse
        self.is_unmapped = is_unmapped
        self.is_secondary = is_secondary
        self.is_supplementary = is_supplementary
        self.query_length = len(query_sequence) if query_sequence else 0


@pytest.fixture
def mock_merged_read() -> MockAlignedSegment:
    """Create a mock merged read for testing."""
    return MockAlignedSegment(
        query_name="MERGED_test_read_001",
        query_sequence="ATCGATCGATCGATCGATCGATCGATCGATCG",
        cigartuples=[(0, 32)],  # 32M
        reference_start=10,
    )


@pytest.fixture
def mock_unmerged_r1_read() -> MockAlignedSegment:
    """Create a mock unmerged R1 read for testing."""
    return MockAlignedSegment(
        query_name="UNMERGED_test_read_001/1",
        query_sequence="ATCGATCGATCGATCGATCG",
        cigartuples=[(0, 20)],  # 20M
        reference_start=15,
    )


@pytest.fixture
def mock_unmerged_r2_read() -> MockAlignedSegment:
    """Create a mock unmerged R2 read for testing."""
    return MockAlignedSegment(
        query_name="UNMERGED_test_read_001/2",
        query_sequence="CGATCGATCGATCGATCGAT",
        cigartuples=[(0, 20)],  # 20M
        reference_start=25,
        is_reverse=True,
    )


@pytest.fixture
def mock_other_read() -> MockAlignedSegment:
    """Create a mock 'other' (single-end/untagged) read for testing."""
    return MockAlignedSegment(
        query_name="nanopore_read_001",
        query_sequence="ATCGATCGATCGATCGATCGATCGATCG",
        cigartuples=[(0, 27)],  # 27M
        reference_start=5,
    )


@pytest.fixture
def mock_complex_cigar_read() -> MockAlignedSegment:
    """Create a mock read with complex CIGAR for testing."""
    return MockAlignedSegment(
        query_name="MERGED_complex_read",
        query_sequence="ATCGATCGATCGATCGATCGATCG",
        # 5M2I3M1D10M3S = 5 + 2 + 3 + 0 + 10 + 3 = 23 bases
        cigartuples=[(0, 5), (1, 2), (0, 3), (2, 1), (0, 10), (4, 3)],
        reference_start=20,
    )


def create_sam_header(reference_sequence: str) -> dict[str, Any]:
    """Create a minimal SAM header for testing."""
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "test_reference", "LN": len(reference_sequence)}],
        "PG": [{"ID": "test", "PN": "trim_aligned_reads_test", "VN": "0.1.0"}],
    }


@pytest.fixture
def empty_sam_file(temp_dir: Path, reference_sequence: str) -> Path:
    """Create an empty SAM file with header only."""
    sam_path = temp_dir / "empty.sam"
    header = create_sam_header(reference_sequence)

    with pysam.AlignmentFile(str(sam_path), "w", header=header):
        pass  # Just create the file with header

    return sam_path


@pytest.fixture
def sample_bam_file(temp_dir: Path, reference_sequence: str) -> Path:
    """Create a sample BAM file with a few test reads for integration testing."""
    bam_path = temp_dir / "sample.bam"
    header = create_sam_header(reference_sequence)

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_file:
        # Create some test reads programmatically
        reads_data = [
            ("MERGED_read_001", "ATCGATCGATCGATCGATCGATCG", [(0, 24)], 10, False),
            ("single_end_read_003", "ATCGATCGATCGATCGATCG", [(0, 20)], 15, False),
            ("UNMERGED_read_002/1", "ATCGATCGATCGATCG", [(0, 16)], 20, False),
            (
                "UNMERGED_read_002/2",
                "CGATCGATCGATCGAT",
                [(0, 16)],
                25,
                True,
            ),  # Changed from 30 to 25 for proper sorting
        ]

        for qname, seq, cigar, ref_start, is_rev in reads_data:
            read = pysam.AlignedSegment()
            read.query_name = qname
            read.query_sequence = seq
            read.query_qualities = [30] * len(seq)
            read.cigartuples = cigar
            read.reference_start = ref_start
            read.reference_id = 0  # First (and only) reference
            read.is_reverse = is_rev
            read.mapping_quality = 60
            read.flag = 16 if is_rev else 0  # Set reverse flag if needed

            bam_file.write(read)

    # Index the BAM file
    pysam.index(str(bam_path))
    return bam_path


@pytest.fixture
def cigar_test_cases() -> list[tuple[str, list[tuple[int, int]]]]:
    """Common CIGAR test cases for unit testing."""
    return [
        ("simple_match", [(0, 20)]),  # 20M
        ("with_insertion", [(0, 10), (1, 3), (0, 10)]),  # 10M3I10M
        ("with_deletion", [(0, 8), (2, 2), (0, 12)]),  # 8M2D12M
        ("with_soft_clips", [(4, 5), (0, 15), (4, 3)]),  # 5S15M3S
        ("complex", [(4, 2), (0, 5), (1, 1), (0, 8), (2, 1), (0, 6), (4, 2)]),  # 2S5M1I8M1D6M2S
        ("only_soft_clips", [(4, 20)]),  # 20S
        ("only_insertions", [(1, 15)]),  # 15I (unusual but valid)
    ]


@pytest.fixture(autouse=True)
def configure_logging_for_tests() -> None:
    """Configure logging for tests to reduce noise."""
    # Remove existing handlers and set to WARNING level for tests
    from loguru import logger

    logger.remove()
    logger.add(sys.stderr, level="WARNING")


# Test markers for organizing test runs
pytest_plugins = []  # Can add plugins here if needed
