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
    ClippingMode,
    ReadCategory,
    TagConfig,
    TrimPolicy,
    _consume_from_left,
    _consume_from_right,
    _convert_to_hard_clips,
    _convert_to_soft_clips,
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

    def test_clipping_mode_configuration(self):
        """Test that clipping mode can be configured properly."""
        # Default should be DELETE mode
        default_policy = TrimPolicy()
        assert default_policy.clipping_mode == ClippingMode.DELETE

        # Can configure soft clip mode
        soft_policy = TrimPolicy(clipping_mode=ClippingMode.SOFT_CLIP)
        assert soft_policy.clipping_mode == ClippingMode.SOFT_CLIP

        # Can configure hard clip mode
        hard_policy = TrimPolicy(clipping_mode=ClippingMode.HARD_CLIP)
        assert hard_policy.clipping_mode == ClippingMode.HARD_CLIP

        # Clipping mode should be immutable
        with pytest.raises(AttributeError):
            default_policy.clipping_mode = ClippingMode.SOFT_CLIP


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
        tag_config = TagConfig()  # Use default prefixes
        assert ReadCategory.classify(qname, tag_config) == expected_category

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
        assert result_cigar[1] == CigarOp(
            0, 10
        )  # 10M remaining (consumed 2 from beginning)
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
        assert result_cigar[1] == CigarOp(
            0, 10
        )  # 10M remaining (consumed 2 from beginning)
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
        with pytest.raises(
            AssertionError, match="Left trim amount must be non-negative"
        ):
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
        assert (
            mock_read.query_sequence == "CGATCGATCGA"
        )  # Same sequence result as forward
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
            query_qualities=[
                30,
                25,
                35,
                40,
                30,
                25,
                35,
                40,
                30,
                25,
                35,
                40,
                30,
                25,
                35,
                40,
            ],
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


class TestClippingFunctions:
    """Test the new clipping conversion functions."""

    def test_convert_to_soft_clips_basic(self):
        """Test basic soft clip conversion."""
        # 10M -> 3S10M5S with 3 left, 5 right trim
        input_cig = Cigar([CigarOp(0, 10)])  # 10M
        result = _convert_to_soft_clips(input_cig, 3, 5)
        expected = [(4, 3), (0, 10), (4, 5)]  # 3S10M5S
        assert result.to_pysam() == expected

    def test_convert_to_soft_clips_zero_trim(self):
        """Test soft clip conversion with zero trim amounts."""
        input_cig = Cigar([CigarOp(0, 10)])  # 10M
        result = _convert_to_soft_clips(input_cig, 0, 0)
        expected = [(0, 10)]  # Just 10M, no clips added
        assert result.to_pysam() == expected

    def test_convert_to_soft_clips_left_only(self):
        """Test soft clip conversion with only left trim."""
        input_cig = Cigar([CigarOp(0, 10)])  # 10M
        result = _convert_to_soft_clips(input_cig, 3, 0)
        expected = [(4, 3), (0, 10)]  # 3S10M
        assert result.to_pysam() == expected

    def test_convert_to_soft_clips_right_only(self):
        """Test soft clip conversion with only right trim."""
        input_cig = Cigar([CigarOp(0, 10)])  # 10M
        result = _convert_to_soft_clips(input_cig, 0, 5)
        expected = [(0, 10), (4, 5)]  # 10M5S
        assert result.to_pysam() == expected

    def test_convert_to_soft_clips_validation(self):
        """Test input validation for soft clip conversion."""
        input_cig = Cigar([CigarOp(0, 10)])

        # Negative left trim should raise assertion
        with pytest.raises(AssertionError, match="Left trim must be non-negative"):
            _convert_to_soft_clips(input_cig, -1, 5)

        # Negative right trim should raise assertion
        with pytest.raises(AssertionError, match="Right trim must be non-negative"):
            _convert_to_soft_clips(input_cig, 3, -1)

    def test_convert_to_hard_clips_basic(self):
        """Test basic hard clip conversion."""
        # 10M -> 2H10M3H with 2 left, 3 right trim
        input_cig = Cigar([CigarOp(0, 10)])  # 10M
        result = _convert_to_hard_clips(input_cig, 2, 3)
        expected = [(5, 2), (0, 10), (5, 3)]  # 2H10M3H
        assert result.to_pysam() == expected

    def test_convert_to_hard_clips_zero_trim(self):
        """Test hard clip conversion with zero trim amounts."""
        input_cig = Cigar([CigarOp(0, 10)])  # 10M
        result = _convert_to_hard_clips(input_cig, 0, 0)
        expected = [(0, 10)]  # Just 10M, no clips added
        assert result.to_pysam() == expected

    def test_convert_to_hard_clips_left_only(self):
        """Test hard clip conversion with only left trim."""
        input_cig = Cigar([CigarOp(0, 10)])  # 10M
        result = _convert_to_hard_clips(input_cig, 4, 0)
        expected = [(5, 4), (0, 10)]  # 4H10M
        assert result.to_pysam() == expected

    def test_convert_to_hard_clips_right_only(self):
        """Test hard clip conversion with only right trim."""
        input_cig = Cigar([CigarOp(0, 10)])  # 10M
        result = _convert_to_hard_clips(input_cig, 0, 6)
        expected = [(0, 10), (5, 6)]  # 10M6H
        assert result.to_pysam() == expected

    def test_convert_to_hard_clips_validation(self):
        """Test input validation for hard clip conversion."""
        input_cig = Cigar([CigarOp(0, 10)])

        # Negative left trim should raise assertion
        with pytest.raises(AssertionError, match="Left trim must be non-negative"):
            _convert_to_hard_clips(input_cig, -1, 5)

        # Negative right trim should raise assertion
        with pytest.raises(AssertionError, match="Right trim must be non-negative"):
            _convert_to_hard_clips(input_cig, 3, -1)

    def test_clipping_edge_cases(self):
        """Test edge cases for both clipping functions."""
        empty_cig = Cigar([])

        # Zero trim amounts should be no-op for both functions
        assert _convert_to_soft_clips(empty_cig, 0, 0).to_pysam() == []
        assert _convert_to_hard_clips(empty_cig, 0, 0).to_pysam() == []

        # Empty CIGAR with trim amounts should create only clips
        soft_result = _convert_to_soft_clips(empty_cig, 3, 5)
        hard_result = _convert_to_hard_clips(empty_cig, 3, 5)
        assert soft_result.to_pysam() == [(4, 8)]  # Compacted 3S + 5S = 8S
        assert hard_result.to_pysam() == [(5, 8)]  # Compacted 3H + 5H = 8H


class TestClippingModes:
    """Test the new clipping modes in trim_alignment_in_place."""

    def test_delete_mode_behavior(self):
        """Test DELETE mode preserves original behavior."""
        mock = MockAlignedSegment(
            query_name="delete_test",
            query_sequence="ATCGATCGATCG",  # 12 bases
            query_qualities=[30] * 12,
            cigartuples=[(0, 12)],  # 12M
            reference_start=100,
            is_reverse=False,
        )

        trim_alignment_in_place(mock, 3, 2, ClippingMode.DELETE)

        # Should trim sequence and advance reference
        assert mock.query_sequence == "GATCGAT"  # 7 bases (trimmed 3 left, 2 right)
        assert mock.cigartuples == [(0, 7)]  # 7M
        assert mock.reference_start == 103  # Advanced by 3
        assert len(mock.query_qualities) == 7

    def test_soft_clip_mode_behavior(self):
        """Test SOFT_CLIP mode preserves sequence and reference position."""
        mock = MockAlignedSegment(
            query_name="soft_test",
            query_sequence="ATCGATCGATCG",  # 12 bases
            query_qualities=[30] * 12,
            cigartuples=[(0, 12)],  # 12M
            reference_start=100,
            is_reverse=False,
        )

        trim_alignment_in_place(mock, 3, 2, ClippingMode.SOFT_CLIP)

        # Should preserve sequence and reference position, add soft clips to CIGAR
        assert mock.query_sequence == "ATCGATCGATCG"  # Unchanged
        assert mock.cigartuples == [(4, 3), (0, 12), (4, 2)]  # 3S12M2S
        assert mock.reference_start == 100  # Unchanged
        assert len(mock.query_qualities) == 12  # Unchanged

    def test_hard_clip_mode_behavior(self):
        """Test HARD_CLIP mode trims sequence but preserves reference position."""
        mock = MockAlignedSegment(
            query_name="hard_test",
            query_sequence="ATCGATCGATCG",  # 12 bases
            query_qualities=[30] * 12,
            cigartuples=[(0, 12)],  # 12M
            reference_start=100,
            is_reverse=False,
        )

        trim_alignment_in_place(mock, 3, 2, ClippingMode.HARD_CLIP)

        # Should trim sequence but preserve reference position, add hard clips to CIGAR
        assert mock.query_sequence == "GATCGAT"  # 7 bases (trimmed 3 left, 2 right)
        assert mock.cigartuples == [(5, 3), (0, 7), (5, 2)]  # 3H7M2H
        assert (
            mock.reference_start == 100
        )  # Unchanged (hard clips don't consume reference)
        assert len(mock.query_qualities) == 7

    def test_clipping_modes_with_reverse_strand(self):
        """Test clipping modes work correctly with reverse strand reads."""
        # For reverse strand: left/right are swapped for CIGAR operations
        mock_delete = MockAlignedSegment(
            query_name="reverse_test",
            query_sequence="ATCGATCGATCG",
            query_qualities=[30] * 12,
            cigartuples=[(0, 12)],
            reference_start=100,
            is_reverse=True,  # REVERSE STRAND
        )

        trim_alignment_in_place(mock_delete, 3, 2, ClippingMode.DELETE)

        # Should still trim sequence correctly in read orientation
        assert mock_delete.query_sequence == "GATCGAT"  # Same trimming pattern
        assert mock_delete.cigartuples == [(0, 7)]  # 7M
        # Reference start should advance by RIGHT amount (2) because strand is reversed
        assert mock_delete.reference_start == 102

    def test_clipping_modes_with_complex_cigar(self):
        """Test clipping modes work with complex CIGAR operations."""
        # Complex CIGAR: 3M2I4M1D3M (total 10 query bases, 8 reference bases)
        mock_soft = MockAlignedSegment(
            query_name="complex_test",
            query_sequence="ATCGATCGAT",  # 10 bases
            query_qualities=[30] * 10,
            cigartuples=[(0, 3), (1, 2), (0, 4), (2, 1), (0, 3)],  # 3M2I4M1D3M
            reference_start=100,
            is_reverse=False,
        )

        trim_alignment_in_place(mock_soft, 2, 1, ClippingMode.SOFT_CLIP)

        # Should preserve sequence and add soft clips around complex CIGAR
        assert mock_soft.query_sequence == "ATCGATCGAT"  # Unchanged
        expected_cigar = [
            (4, 2),
            (0, 3),
            (1, 2),
            (0, 4),
            (2, 1),
            (0, 3),
            (4, 1),
        ]  # 2S + original + 1S
        assert mock_soft.cigartuples == expected_cigar
        assert mock_soft.reference_start == 100  # Unchanged


class TestClippingModesComprehensive:
    """Comprehensive testing of all clipping modes with various CIGAR patterns."""

    def test_all_modes_with_simple_match(self):
        """Test all clipping modes with simple match CIGAR."""
        base_sequence = "ATCGATCGATCGATCG"  # 16 bases
        base_cigar = [(0, 16)]  # 16M

        for mode in [
            ClippingMode.DELETE,
            ClippingMode.SOFT_CLIP,
            ClippingMode.HARD_CLIP,
        ]:
            mock = MockAlignedSegment(
                query_name=f"simple_{mode.name.lower()}",
                query_sequence=base_sequence,
                query_qualities=[30] * 16,
                cigartuples=base_cigar,
                reference_start=100,
                is_reverse=False,
            )

            trim_alignment_in_place(mock, 3, 2, mode)

            match mode:
                case ClippingMode.DELETE:
                    assert mock.query_sequence == "GATCGATCGAT"  # 11 bases (16-3-2)
                    assert mock.cigartuples == [(0, 11)]  # 11M
                    assert mock.reference_start == 103  # Advanced
                case ClippingMode.SOFT_CLIP:
                    assert mock.query_sequence == base_sequence  # Unchanged
                    assert mock.cigartuples == [(4, 3), (0, 16), (4, 2)]  # 3S16M2S
                    assert mock.reference_start == 100  # Unchanged
                case ClippingMode.HARD_CLIP:
                    assert mock.query_sequence == "GATCGATCGAT"  # 11 bases (16-3-2)
                    assert mock.cigartuples == [(5, 3), (0, 11), (5, 2)]  # 3H11M2H
                    assert mock.reference_start == 100  # Unchanged

    def test_all_modes_with_insertions_deletions(self):
        """Test all clipping modes with complex CIGAR including insertions/deletions."""
        # CIGAR: 4M2I6M1D4M = 16 query bases, 15 reference bases
        base_sequence = "ATCGATCGATCGATCG"  # 16 bases (matches CIGAR query consumption)
        base_cigar = [(0, 4), (1, 2), (0, 6), (2, 1), (0, 4)]  # 4M2I6M1D4M

        for mode in [
            ClippingMode.DELETE,
            ClippingMode.SOFT_CLIP,
            ClippingMode.HARD_CLIP,
        ]:
            mock = MockAlignedSegment(
                query_name=f"complex_{mode.name.lower()}",
                query_sequence=base_sequence,
                query_qualities=[30] * 16,  # Match sequence length
                cigartuples=base_cigar,
                reference_start=200,
                is_reverse=False,
            )

            trim_alignment_in_place(mock, 2, 3, mode)

            match mode:
                case ClippingMode.DELETE:
                    # Should trim sequence and consume CIGAR
                    assert len(mock.query_sequence) == 11  # 16 - 2 - 3 = 11
                    # CIGAR consumption should handle the I/D operations correctly
                    assert mock.reference_start > 200  # Should advance

                case ClippingMode.SOFT_CLIP:
                    # Should preserve sequence and add soft clips around complex CIGAR
                    assert mock.query_sequence == base_sequence  # Unchanged
                    expected_cigar = (
                        [(4, 2)] + base_cigar + [(4, 3)]
                    )  # 2S + original + 3S
                    assert mock.cigartuples == expected_cigar
                    assert mock.reference_start == 200  # Unchanged

                case ClippingMode.HARD_CLIP:
                    # Should trim sequence and add hard clips around consumed CIGAR
                    assert (
                        len(mock.query_sequence) == 11
                    )  # Trimmed like DELETE (16-2-3)
                    # Should have hard clips around the consumed CIGAR
                    assert mock.cigartuples[0][0] == 5  # First op should be H
                    assert mock.cigartuples[-1][0] == 5  # Last op should be H
                    assert mock.reference_start == 200  # Unchanged

    def test_all_modes_with_existing_clips(self):
        """Test all clipping modes when reads already have soft/hard clips."""
        # CIGAR: 2S10M3S (already has soft clips) - must have 15 bases total (2+10+3)
        base_sequence = "ATCGATCGATCGATC"  # 15 bases total (matches CIGAR)
        base_cigar = [(4, 2), (0, 10), (4, 3)]  # 2S10M3S

        for mode in [
            ClippingMode.DELETE,
            ClippingMode.SOFT_CLIP,
            ClippingMode.HARD_CLIP,
        ]:
            mock = MockAlignedSegment(
                query_name=f"existing_clips_{mode.name.lower()}",
                query_sequence=base_sequence,
                query_qualities=[30] * 15,  # Match sequence length
                cigartuples=base_cigar,
                reference_start=150,
                is_reverse=False,
            )

            trim_alignment_in_place(mock, 1, 2, mode)

            match mode:
                case ClippingMode.SOFT_CLIP:
                    # Should preserve sequence and add more soft clips
                    assert mock.query_sequence == base_sequence  # Unchanged
                    # Should merge with existing soft clips: (1+2)S + 10M + (2+3)S = 3S10M5S
                    assert mock.cigartuples == [(4, 3), (0, 10), (4, 5)]
                    assert mock.reference_start == 150  # Unchanged

    def test_edge_cases_zero_trim(self):
        """Test all clipping modes with zero trim amounts (should be no-op)."""
        base_sequence = "ATCGATCGATCG"
        base_cigar = [(0, 12)]

        for mode in [
            ClippingMode.DELETE,
            ClippingMode.SOFT_CLIP,
            ClippingMode.HARD_CLIP,
        ]:
            mock = MockAlignedSegment(
                query_name=f"zero_trim_{mode.name.lower()}",
                query_sequence=base_sequence,
                query_qualities=[30] * 12,
                cigartuples=base_cigar,
                reference_start=100,
                is_reverse=False,
            )

            trim_alignment_in_place(mock, 0, 0, mode)

            # All modes should be no-op with zero trim
            assert mock.query_sequence == base_sequence
            assert mock.cigartuples == base_cigar
            assert mock.reference_start == 100

    def test_edge_cases_trim_entire_sequence(self):
        """Test clipping modes when trim amounts exceed sequence length."""
        base_sequence = "ATCG"  # 4 bases
        base_cigar = [(0, 4)]  # 4M

        for mode in [
            ClippingMode.DELETE,
            ClippingMode.SOFT_CLIP,
            ClippingMode.HARD_CLIP,
        ]:
            mock = MockAlignedSegment(
                query_name=f"over_trim_{mode.name.lower()}",
                query_sequence=base_sequence,
                query_qualities=[30] * 4,
                cigartuples=base_cigar,
                reference_start=100,
                is_reverse=False,
            )

            trim_alignment_in_place(
                mock, 2, 3, mode
            )  # Trim 5 bases from 4-base sequence

            match mode:
                case ClippingMode.DELETE:
                    assert mock.query_sequence == ""  # Empty after over-trimming
                    assert mock.cigartuples == []  # Empty CIGAR
                    assert mock.reference_start == 102  # Advanced by 2

                case ClippingMode.SOFT_CLIP:
                    assert mock.query_sequence == base_sequence  # Preserved
                    assert mock.cigartuples == [(4, 2), (0, 4), (4, 3)]  # 2S4M3S
                    assert mock.reference_start == 100  # Unchanged

                case ClippingMode.HARD_CLIP:
                    assert mock.query_sequence == ""  # Empty after trimming
                    assert mock.cigartuples == [
                        (5, 5)
                    ]  # Just hard clips (compacted 2H+3H=5H)
                    assert mock.reference_start == 100  # Unchanged

    def test_error_conditions_unmapped_reads(self):
        """Test clipping modes with unmapped reads (no CIGAR)."""
        base_sequence = "ATCGATCG"

        for mode in [
            ClippingMode.DELETE,
            ClippingMode.SOFT_CLIP,
            ClippingMode.HARD_CLIP,
        ]:
            mock = MockAlignedSegment(
                query_name=f"unmapped_{mode.name.lower()}",
                query_sequence=base_sequence,
                query_qualities=[30] * 8,
                cigartuples=None,  # No CIGAR (unmapped)
                reference_start=None,  # No reference position
                is_reverse=False,
            )

            # Should handle unmapped reads gracefully
            trim_alignment_in_place(mock, 2, 1, mode)

            # For unmapped reads, only sequence/quality trimming happens (no CIGAR)
            expected_seq = "CGATC"  # 8 - 2 - 1 = 5 bases
            assert mock.query_sequence == expected_seq
            assert mock.cigartuples is None  # Still no CIGAR
            assert mock.reference_start is None  # Still unmapped

    def test_strand_aware_behavior_all_modes(self):
        """Test that all clipping modes handle forward/reverse strand correctly."""
        base_sequence = "ATCGATCGATCG"  # 12 bases
        base_cigar = [(0, 12)]  # 12M

        # Test both forward and reverse strand
        for is_reverse in [False, True]:
            for mode in [
                ClippingMode.DELETE,
                ClippingMode.SOFT_CLIP,
                ClippingMode.HARD_CLIP,
            ]:
                mock = MockAlignedSegment(
                    query_name=f"strand_{is_reverse}_{mode.name.lower()}",
                    query_sequence=base_sequence,
                    query_qualities=[30] * 12,
                    cigartuples=base_cigar,
                    reference_start=100,
                    is_reverse=is_reverse,
                )

                trim_alignment_in_place(
                    mock, 3, 2, mode
                )  # 3 left, 2 right in read orientation

                # Sequence trimming should always be in read orientation (same result)
                if mode in [ClippingMode.DELETE, ClippingMode.HARD_CLIP]:
                    assert (
                        mock.query_sequence == "GATCGAT"
                    )  # Same trimming regardless of strand
                else:  # SOFT_CLIP
                    assert mock.query_sequence == base_sequence  # Preserved

                # Reference position behavior depends on strand and mode
                if mode == ClippingMode.DELETE:
                    if is_reverse:
                        # For reverse reads: right trim (2) affects reference position
                        assert mock.reference_start == 102  # Advanced by right trim (2)
                    else:
                        # For forward reads: left trim (3) affects reference position
                        assert mock.reference_start == 103  # Advanced by left trim (3)
                else:  # SOFT_CLIP or HARD_CLIP
                    assert mock.reference_start == 100  # Unchanged for both clip modes


# Marker-based test organization
pytestmark = pytest.mark.unit
