# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pytest",
#     "unittest",
# ]
# ///
"""
Unit tests for trim_aligned_reads.py

This module tests the core functionality of the alignment-aware read trimming tool,
focusing on CIGAR string manipulation, read classification, and edge cases.
"""

import pytest
from trim_aligned_reads import (
    Cigar,
    CigarOp,
    ReadCategory,
    TrimPolicy,
    _consume_from_left,
    _consume_from_right,
    configure_logging,
    trim_alignment_in_place,
)

from .conftest import MockAlignedSegment


class TestTrimPolicy:
    """Test the TrimPolicy dataclass."""

    def test_default_construction(self):
        """Test TrimPolicy with default values."""
        policy = TrimPolicy()
        assert policy.merged_left == 30
        assert policy.merged_right == 30
        assert policy.r1_left == 30
        assert policy.r2_right == 30
        assert policy.single_left == 30
        assert policy.single_right == 30
        assert policy.min_len == 20

    def test_custom_construction(self):
        """Test TrimPolicy with custom values."""
        policy = TrimPolicy(
            merged_left=15,
            merged_right=10,
            r1_left=20,
            r2_right=25,
            single_left=5,
            single_right=3,
            min_len=15,
        )
        assert policy.merged_left == 15
        assert policy.merged_right == 10
        assert policy.r1_left == 20
        assert policy.r2_right == 25
        assert policy.single_left == 5
        assert policy.single_right == 3
        assert policy.min_len == 15

    def test_immutability(self):
        """Test that TrimPolicy is immutable (frozen dataclass)."""
        policy = TrimPolicy()
        with pytest.raises(AttributeError):
            policy.merged_left = 50


class TestReadCategory:
    """Test read classification and trimming extent calculation."""

    @pytest.mark.parametrize(
        "qname,expected_category",
        [
            ("MERGED_read_001", ReadCategory.MERGED),
            ("MERGED_SRR123456", ReadCategory.MERGED),
            ("UNMERGED_read_002/1", ReadCategory.UNMERGED_R1),
            ("UNMERGED_SRR123456/1", ReadCategory.UNMERGED_R1),
            ("UNMERGED_read_002/2", ReadCategory.UNMERGED_R2),
            ("UNMERGED_SRR123456/2", ReadCategory.UNMERGED_R2),
            ("nanopore_read", ReadCategory.OTHER),
            ("SRR123456", ReadCategory.OTHER),
            ("regular_read_name", ReadCategory.OTHER),
            ("MERGED_partial", ReadCategory.MERGED),  # Still starts with MERGED_
            ("UNMERGED_no_suffix", ReadCategory.OTHER),  # No /1 or /2
            ("", ReadCategory.OTHER),  # Empty name
        ],
    )
    def test_classify(self, qname: str, expected_category: ReadCategory):
        """Test read name classification."""
        assert ReadCategory.classify(qname) == expected_category

    def test_trim_extents_merged(self, default_trim_policy):
        """Test trim extents for merged reads."""
        left, right = ReadCategory.MERGED.trim_extents(default_trim_policy)
        assert left == 10
        assert right == 10

    def test_trim_extents_unmerged_r1(self, default_trim_policy):
        """Test trim extents for unmerged R1 reads."""
        left, right = ReadCategory.UNMERGED_R1.trim_extents(default_trim_policy)
        assert left == 10
        assert right == 0  # No right trimming for R1

    def test_trim_extents_unmerged_r2(self, default_trim_policy):
        """Test trim extents for unmerged R2 reads."""
        left, right = ReadCategory.UNMERGED_R2.trim_extents(default_trim_policy)
        assert left == 0  # No left trimming for R2
        assert right == 10

    def test_trim_extents_other(self, default_trim_policy):
        """Test trim extents for other reads."""
        left, right = ReadCategory.OTHER.trim_extents(default_trim_policy)
        assert left == 5  # single_left
        assert right == 5  # single_right

    def test_trim_extents_negative_values(self):
        """Test that negative trim values are handled correctly."""
        policy = TrimPolicy(
            merged_left=-5,  # Negative values should be max(0, value)
            merged_right=10,
            r1_left=-10,
            r2_right=15,
            single_left=-1,
            single_right=-3,
            min_len=20,
        )

        # All negative values should become 0
        left, right = ReadCategory.MERGED.trim_extents(policy)
        assert left == 0  # max(0, -5)
        assert right == 10

        left, right = ReadCategory.UNMERGED_R1.trim_extents(policy)
        assert left == 0  # max(0, -10)
        assert right == 0

        left, right = ReadCategory.OTHER.trim_extents(policy)
        assert left == 0  # max(0, -1)
        assert right == 0  # max(0, -3)


class TestCigarOp:
    """Test CIGAR operation handling."""

    def test_from_tuple(self):
        """Test CigarOp creation from tuple."""
        op = CigarOp.from_tuple((0, 10))
        assert op.op == 0
        assert op.length == 10

    def test_to_tuple(self):
        """Test CigarOp conversion to tuple."""
        op = CigarOp(2, 5)
        assert CigarOp.to_tuple(op) == (2, 5)

    def test_round_trip_conversion(self):
        """Test tuple -> CigarOp -> tuple conversion."""
        original = (1, 15)
        op = CigarOp.from_tuple(original)
        result = CigarOp.to_tuple(op)
        assert result == original


class TestCigar:
    """Test CIGAR string manipulation."""

    def test_from_pysam_none(self):
        """Test handling of None CIGAR."""
        result = Cigar.from_pysam(None)
        assert result is None

    def test_from_pysam_empty(self):
        """Test handling of empty CIGAR."""
        result = Cigar.from_pysam([])
        assert isinstance(result, Cigar)
        assert len(result) == 0

    def test_from_pysam_normal(self):
        """Test normal CIGAR conversion."""
        pysam_cigar = [(0, 10), (1, 2), (0, 5)]
        result = Cigar.from_pysam(pysam_cigar)
        assert len(result) == 3
        assert result[0] == CigarOp(0, 10)
        assert result[1] == CigarOp(1, 2)
        assert result[2] == CigarOp(0, 5)

    def test_to_pysam(self):
        """Test CIGAR conversion back to pysam format."""
        cigar = Cigar([CigarOp(0, 8), CigarOp(2, 3), CigarOp(0, 7)])
        result = cigar.to_pysam()
        expected = [(0, 8), (2, 3), (0, 7)]
        assert result == expected

    def test_push_compact_new_operation(self):
        """Test adding a new operation type."""
        cigar = Cigar([CigarOp(0, 5)])
        cigar.push_compact(1, 3)
        assert len(cigar) == 2
        assert cigar[-1] == CigarOp(1, 3)

    def test_push_compact_same_operation(self):
        """Test merging operations of the same type."""
        cigar = Cigar([CigarOp(0, 5)])
        cigar.push_compact(0, 7)
        assert len(cigar) == 1
        assert cigar[0] == CigarOp(0, 12)  # 5 + 7

    def test_push_compact_zero_length(self):
        """Test that zero-length operations are ignored."""
        cigar = Cigar([CigarOp(0, 5)])
        original_len = len(cigar)
        cigar.push_compact(1, 0)
        assert len(cigar) == original_len  # No change

    def test_push_compact_negative_length(self):
        """Test that negative-length operations are ignored."""
        cigar = Cigar([CigarOp(0, 5)])
        original_len = len(cigar)
        cigar.push_compact(1, -3)
        assert len(cigar) == original_len  # No change

    def test_push_compact_invalid_operation(self):
        """Test assertion for invalid CIGAR operation codes."""
        cigar = Cigar()
        with pytest.raises(AssertionError, match="Invalid CIGAR operation code"):
            cigar.push_compact(9, 5)  # Invalid op code
        with pytest.raises(AssertionError, match="Invalid CIGAR operation code"):
            cigar.push_compact(-1, 5)  # Invalid op code

    def test_push_compact_overflow_protection(self):
        """Test overflow protection in compaction."""
        cigar = Cigar([CigarOp(0, (1 << 28) - 2)])  # Near maximum
        with pytest.raises(AssertionError, match="CIGAR length overflow"):
            cigar.push_compact(0, 5)  # This would overflow


class TestConsumeFromLeft:
    """Test left-side CIGAR consumption."""

    def test_consume_zero(self):
        """Test consuming zero bases."""
        cigar = Cigar([CigarOp(0, 10)])
        result_cigar, ref_advance = _consume_from_left(cigar, 0)
        assert result_cigar == cigar
        assert ref_advance == 0

    def test_consume_empty_cigar(self):
        """Test consuming from empty CIGAR."""
        cigar = Cigar([])
        result_cigar, ref_advance = _consume_from_left(cigar, 5)
        assert len(result_cigar) == 0
        assert ref_advance == 0

    def test_consume_simple_match(self):
        """Test consuming from simple match operation."""
        cigar = Cigar([CigarOp(0, 20)])  # 20M
        result_cigar, ref_advance = _consume_from_left(cigar, 8)
        assert len(result_cigar) == 1
        assert result_cigar[0] == CigarOp(0, 12)  # 12M remaining
        assert ref_advance == 8  # 8 reference bases consumed

    def test_consume_with_insertion(self):
        """Test consuming across match and insertion."""
        # 10M5I10M - consuming 12 bases should take all of 10M, all of 5I, 2 of second 10M
        # Remaining: 3I from insertion + 10M from last match = 2 operations
        cigar = Cigar([CigarOp(0, 10), CigarOp(1, 5), CigarOp(0, 10)])
        result_cigar, ref_advance = _consume_from_left(cigar, 12)
        assert len(result_cigar) == 2
        assert result_cigar[0] == CigarOp(1, 3)  # 3I remaining from insertion
        assert result_cigar[1] == CigarOp(0, 10)  # 10M remaining (consumed 2 from beginning)
        assert ref_advance == 10  # Only M operations advance reference

    def test_consume_with_deletion(self):
        """Test consuming with deletion (D operations don't consume query)."""
        # 8M2D12M - deletions don't consume query bases
        # Consuming 10 query bases: all 8M + all 2D (no query) + 2M from last operation
        # Remaining: 2D + 10M from last operation = 2 operations
        cigar = Cigar([CigarOp(0, 8), CigarOp(2, 2), CigarOp(0, 12)])
        result_cigar, ref_advance = _consume_from_left(cigar, 10)
        assert len(result_cigar) == 2
        assert result_cigar[0] == CigarOp(2, 2)  # 2D remaining
        assert result_cigar[1] == CigarOp(0, 10)  # 10M remaining (consumed 2 from beginning)
        assert ref_advance == 10  # 8 from first M + 2 from D

    def test_consume_soft_clips(self):
        """Test consuming soft clips."""
        # 5S15M3S - soft clips consume query but not reference
        cigar = Cigar([CigarOp(4, 5), CigarOp(0, 15), CigarOp(4, 3)])
        result_cigar, ref_advance = _consume_from_left(cigar, 8)
        # Should consume all 5S + 3M
        assert len(result_cigar) == 2
        assert result_cigar[0] == CigarOp(0, 12)  # 12M remaining
        assert result_cigar[1] == CigarOp(4, 3)  # 3S remaining
        assert ref_advance == 3  # Only 3 bases from M operation

    def test_consume_more_than_available(self):
        """Test consuming more bases than available."""
        cigar = Cigar([CigarOp(0, 10)])  # 10M
        result_cigar, ref_advance = _consume_from_left(cigar, 15)
        # Should consume all available
        assert len(result_cigar) == 0 or all(op.length == 0 for op in result_cigar)
        assert ref_advance == 10

    def test_consume_negative_trim(self):
        """Test assertion for negative trim amount."""
        cigar = Cigar([CigarOp(0, 10)])
        with pytest.raises(AssertionError, match="Left trim amount must be non-negative"):
            _consume_from_left(cigar, -5)

    def test_consume_invalid_cigar(self):
        """Test assertion for invalid CIGAR operations."""
        cigar = Cigar([CigarOp(9, 10)])  # Invalid operation
        with pytest.raises(AssertionError, match="CIGAR contains invalid operations"):
            _consume_from_left(cigar, 5)


class TestConsumeFromRight:
    """Test right-side CIGAR consumption."""

    def test_consume_zero_right(self):
        """Test consuming zero bases from right."""
        cigar = Cigar([CigarOp(0, 10)])
        result_cigar = _consume_from_right(cigar, 0)
        assert result_cigar == cigar

    def test_consume_simple_match_right(self):
        """Test consuming from simple match operation from right."""
        cigar = Cigar([CigarOp(0, 20)])  # 20M
        result_cigar = _consume_from_right(cigar, 8)
        assert len(result_cigar) == 1
        assert result_cigar[0] == CigarOp(0, 12)  # 12M remaining

    def test_consume_complex_right(self):
        """Test consuming from complex CIGAR from right."""
        # 5M3I10M2S - consume 8 bases from right
        cigar = Cigar([CigarOp(0, 5), CigarOp(1, 3), CigarOp(0, 10), CigarOp(4, 2)])
        result_cigar = _consume_from_right(cigar, 8)
        # Should consume 2S + 6M, leaving 5M3I4M
        assert len(result_cigar) == 3
        assert result_cigar[0] == CigarOp(0, 5)
        assert result_cigar[1] == CigarOp(1, 3)
        assert result_cigar[2] == CigarOp(0, 4)  # 4M remaining

    def test_consume_with_non_query_ops_right(self):
        """Test consuming from right with deletions."""
        # 10M2D8M - deletions don't consume query
        cigar = Cigar([CigarOp(0, 10), CigarOp(2, 2), CigarOp(0, 8)])
        result_cigar = _consume_from_right(cigar, 5)
        # Should consume 5M from the last operation, leaving 10M2D3M
        assert len(result_cigar) == 3
        assert result_cigar[0] == CigarOp(0, 10)
        assert result_cigar[1] == CigarOp(2, 2)
        assert result_cigar[2] == CigarOp(0, 3)  # 3M remaining

    def test_consume_more_than_available_right(self):
        """Test consuming more bases than available from right."""
        cigar = Cigar([CigarOp(0, 5)])  # 5M
        result_cigar = _consume_from_right(cigar, 10)
        # Should consume everything
        assert len(result_cigar) == 0 or all(op.length == 0 for op in result_cigar)


class TestTrimAlignmentInPlace:
    """Test the main trimming function with mock objects."""

    def test_trim_no_sequence(self):
        """Test trimming read with no sequence."""
        mock_read = MockAlignedSegment(query_sequence=None)
        # Should return early without error
        trim_alignment_in_place(mock_read, 5, 5)
        assert mock_read.query_sequence is None

    def test_trim_no_cigar(self):
        """Test trimming read with no CIGAR (unmapped)."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCGATCGATCG",  # 16 bases
            cigartuples=None,
        )
        trim_alignment_in_place(mock_read, 3, 2)
        # Should trim sequence: remove 3 from left, 2 from right
        # 16 - 3 - 2 = 11 bases remaining
        expected_len = 16 - 3 - 2
        assert len(mock_read.query_sequence) == expected_len
        assert mock_read.query_sequence == "GATCGATCGAT"  # 11 bases remaining: [3:14]
        assert len(mock_read.query_qualities) == expected_len

    def test_trim_forward_read(self):
        """Test trimming forward strand read."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCGATCGATCG",  # 16 bases
            cigartuples=[(0, 16)],  # 16M
            reference_start=10,
            is_reverse=False,
        )
        trim_alignment_in_place(mock_read, 2, 3)
        # Forward read: left=2, right=3, so keep bases 2:(16-3) = 2:13 (11 bases)
        expected_len = 16 - 2 - 3
        assert len(mock_read.query_sequence) == expected_len
        assert mock_read.query_sequence == "CGATCGATCGA"  # [2:13]
        assert mock_read.reference_start == 12  # 10 + 2 (ref advance)
        assert len(mock_read.query_qualities) == expected_len

    def test_trim_reverse_read(self):
        """Test trimming reverse strand read."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCGATCGATCG",  # 16 bases
            cigartuples=[(0, 16)],  # 16M
            reference_start=20,
            is_reverse=True,
        )
        trim_alignment_in_place(mock_read, 2, 3)
        # Reverse read: CIGAR trim amounts are swapped
        # cig_left = 3 (right), cig_right = 2 (left)
        # Sequence trim: still left=2, right=3 in read orientation
        expected_len = 16 - 2 - 3
        assert len(mock_read.query_sequence) == expected_len
        assert mock_read.query_sequence == "CGATCGATCGA"  # Same sequence result as forward
        assert mock_read.reference_start == 23  # 20 + 3 (cig_left ref advance)
        assert len(mock_read.query_qualities) == expected_len

    def test_trim_negative_amounts(self):
        """Test that negative trim amounts are normalized to zero."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCGATCG",
            cigartuples=[(0, 12)],
        )
        # Negative amounts should be normalized to 0
        trim_alignment_in_place(mock_read, -5, -3)
        assert mock_read.query_sequence == "ATCGATCGATCG"  # Unchanged

    def test_trim_over_length(self):
        """Test trimming more than sequence length."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCG",  # 4 bases
            cigartuples=[(0, 4)],
        )
        trim_alignment_in_place(mock_read, 3, 3)
        # Should trim to empty or minimal sequence
        assert len(mock_read.query_sequence) == 0  # All trimmed away
        if mock_read.query_qualities is not None:
            assert len(mock_read.query_qualities) == 0

    def test_sequence_quality_consistency(self):
        """Test that sequence and quality arrays stay synchronized."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCGATCGATCG",
            query_qualities=[30, 25, 35, 40, 30, 25, 35, 40, 30, 25, 35, 40, 30, 25, 35, 40],
            cigartuples=[(0, 16)],
        )
        trim_alignment_in_place(mock_read, 4, 2)
        expected_seq_len = 10  # 16 - 4 - 2
        assert len(mock_read.query_sequence) == expected_seq_len
        assert len(mock_read.query_qualities) == expected_seq_len


class TestLogging:
    """Test logging configuration."""

    def test_configure_logging_default(self):
        """Test default logging configuration."""
        configure_logging(0, 0)
        # Should not raise any exceptions

    def test_configure_logging_verbose(self):
        """Test verbose logging configuration."""
        configure_logging(2, 0)  # DEBUG level
        # Should not raise any exceptions

    def test_configure_logging_quiet(self):
        """Test quiet logging configuration."""
        configure_logging(0, 2)  # ERROR level
        # Should not raise any exceptions


# Marker-based test organization
pytestmark = pytest.mark.unit
