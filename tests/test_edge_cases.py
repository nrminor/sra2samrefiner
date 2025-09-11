# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "hypothesis",
#     "pysam",
#     "pytest",
#     "unittest",
# ]
# ///
"""
Edge case and property-based tests for trim_aligned_reads.py

This module focuses on boundary conditions, error cases, and property-based
testing to ensure robust behavior under unusual conditions.
"""

import pytest
import trim_aligned_reads
from hypothesis import given, settings
from hypothesis import strategies as st
from trim_aligned_reads import (
    Cigar,
    CigarOp,
    ReadCategory,
    TrimPolicy,
    _consume_from_left,
    _consume_from_right,
    trim_alignment_in_place,
)

from .conftest import MockAlignedSegment


class TestBoundaryConditions:
    """Test boundary conditions and edge cases."""

    def test_empty_sequence_trimming(self):
        """Test trimming reads with empty sequences."""
        mock_read = MockAlignedSegment(
            query_sequence="",
            query_qualities=[],
            cigartuples=[],
        )
        # Should handle gracefully without errors
        trim_alignment_in_place(mock_read, 5, 5)
        assert mock_read.query_sequence == ""

    def test_single_base_sequence(self):
        """Test trimming single-base sequences."""
        mock_read = MockAlignedSegment(
            query_sequence="A",
            query_qualities=[30],
            cigartuples=[(0, 1)],  # 1M
        )
        trim_alignment_in_place(mock_read, 0, 1)  # Trim right only
        assert mock_read.query_sequence == ""
        assert len(mock_read.query_qualities) == 0

    def test_trim_exactly_sequence_length(self):
        """Test trimming exactly the sequence length."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCG",  # 4 bases
            cigartuples=[(0, 4)],
        )
        trim_alignment_in_place(mock_read, 2, 2)  # Trim exactly all bases
        assert mock_read.query_sequence == ""

    def test_zero_trim_amounts(self):
        """Test that zero trim amounts leave sequence unchanged."""
        original_seq = "ATCGATCGATCG"
        mock_read = MockAlignedSegment(
            query_sequence=original_seq,
            cigartuples=[(0, len(original_seq))],
        )
        trim_alignment_in_place(mock_read, 0, 0)
        assert mock_read.query_sequence == original_seq

    def test_massive_trim_amounts(self):
        """Test trimming with amounts much larger than sequence."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCG",  # 4 bases
            cigartuples=[(0, 4)],
        )
        trim_alignment_in_place(mock_read, 1000, 1000)
        assert mock_read.query_sequence == ""

    def test_quality_array_none(self):
        """Test handling when quality array is explicitly None."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCG",
            query_qualities=None,
            cigartuples=[(0, 8)],
        )
        # Explicitly set to None after construction to override default behavior
        mock_read.query_qualities = None
        trim_alignment_in_place(mock_read, 2, 1)
        assert mock_read.query_sequence == "CGATC"  # Trimmed correctly
        assert mock_read.query_qualities is None  # Still None

    def test_reference_start_none(self):
        """Test handling when reference_start is None."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCG",
            cigartuples=[(0, 8)],
            reference_start=None,  # Unmapped read
        )
        # Should handle gracefully
        trim_alignment_in_place(mock_read, 1, 1)
        assert len(mock_read.query_sequence) == 6  # 8 - 1 - 1


class TestCigarEdgeCases:
    """Test CIGAR string edge cases."""

    def test_empty_cigar_list(self):
        """Test operations on empty CIGAR."""
        empty_cigar = Cigar([])
        result_left, ref_advance = _consume_from_left(empty_cigar, 5)
        assert len(result_left) == 0
        assert ref_advance == 0

        result_right = _consume_from_right(empty_cigar, 5)
        assert len(result_right) == 0

    def test_single_operation_cigars(self):
        """Test CIGAR with single operations."""
        test_cases = [
            (0, 20),  # 20M
            (1, 15),  # 15I
            (2, 10),  # 10D
            (4, 25),  # 25S
            (5, 30),  # 30H
        ]

        for op, length in test_cases:
            cigar = Cigar([CigarOp(op, length)])
            # Try consuming from left
            result, ref_advance = _consume_from_left(cigar, 5)
            assert isinstance(result, Cigar)
            assert ref_advance >= 0

            # Try consuming from right
            result = _consume_from_right(cigar, 5)
            assert isinstance(result, Cigar)

    def test_only_non_query_consuming_operations(self):
        """Test CIGAR with only deletions, introns, etc."""
        cigar = Cigar([CigarOp(2, 10), CigarOp(3, 5)])  # 10D5N
        result, ref_advance = _consume_from_left(cigar, 8)
        # Should consume 0 query bases since D and N don't consume query
        # Since no query bases are available to consume, no operations are processed
        assert ref_advance == 0  # No reference advance because no operations consumed
        assert len(result) == 2  # Operations preserved unchanged
        assert result[0] == CigarOp(2, 10)  # 10D unchanged
        assert result[1] == CigarOp(3, 5)   # 5N unchanged

    def test_mixed_consuming_non_consuming(self):
        """Test complex mix of query-consuming and non-consuming operations."""
        # 5M3D2I4N8M1D3S
        cigar = Cigar(
            [
                CigarOp(0, 5),  # 5M - consumes query and ref
                CigarOp(2, 3),  # 3D - consumes ref only
                CigarOp(1, 2),  # 2I - consumes query only
                CigarOp(3, 4),  # 4N - consumes ref only
                CigarOp(0, 8),  # 8M - consumes query and ref
                CigarOp(2, 1),  # 1D - consumes ref only
                CigarOp(4, 3),  # 3S - consumes query only
            ],
        )

        # Query-consuming bases: 5 + 2 + 8 + 3 = 18
        # Try consuming 10 query bases from left
        result, ref_advance = _consume_from_left(cigar, 10)
        assert ref_advance > 0  # Should advance some reference bases
        assert len(result) >= 0  # Should produce valid result

    def test_cigar_compaction_edge_cases(self):
        """Test CIGAR compaction with edge cases."""
        cigar = Cigar()

        # Add multiple operations of the same type
        cigar.push_compact(0, 5)
        cigar.push_compact(0, 7)
        cigar.push_compact(0, 3)
        assert len(cigar) == 1
        assert cigar[0] == CigarOp(0, 15)  # All merged

        # Add different operation
        cigar.push_compact(1, 4)
        assert len(cigar) == 2
        assert cigar[1] == CigarOp(1, 4)

        # Add zero length (should be ignored)
        cigar.push_compact(2, 0)
        assert len(cigar) == 2  # No change


class TestStrandOrientationEdgeCases:
    """Test edge cases in strand orientation handling."""

    def test_reverse_read_cigar_swap(self):
        """Test that reverse reads properly swap CIGAR trim amounts."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCGATCG",  # 12 bases
            cigartuples=[(0, 12)],  # 12M
            reference_start=100,
            is_reverse=True,
        )

        # For reverse read, left and right should be swapped for CIGAR operations
        trim_alignment_in_place(mock_read, 2, 3)

        # Should advance reference by 3 (the "right" trim becomes "left" CIGAR trim)
        assert mock_read.reference_start == 103  # 100 + 3

        # Sequence should still be trimmed in read orientation
        assert len(mock_read.query_sequence) == 7  # 12 - 2 - 3

    def test_forward_read_normal_behavior(self):
        """Test that forward reads behave normally."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCGATCG",  # 12 bases
            cigartuples=[(0, 12)],  # 12M
            reference_start=100,
            is_reverse=False,
        )

        trim_alignment_in_place(mock_read, 2, 3)

        # Should advance reference by 2 (the left trim)
        assert mock_read.reference_start == 102  # 100 + 2

        # Sequence should be trimmed normally
        assert len(mock_read.query_sequence) == 7  # 12 - 2 - 3


class TestErrorHandling:
    """Test error handling and validation."""

    def test_invalid_trim_policy_values(self):
        """Test handling of invalid TrimPolicy values in processing."""
        # The code should handle negative values by using max(0, value)
        policy = TrimPolicy(
            merged_left=-10,
            merged_right=-5,
            r1_left=5,
            r2_right=3,
            min_len=1,  # Very small but valid
        )

        # Should not raise exceptions
        left, right = ReadCategory.MERGED.trim_extents(policy)
        assert left == 0  # max(0, -10)
        assert right == 0  # max(0, -5)

    def test_sequence_quality_length_mismatch(self):
        """Test assertion when sequence and quality lengths don't match."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCG",  # 8 bases
            query_qualities=[30, 30, 30, 30, 30],  # 5 qualities - MISMATCH
        )

        with pytest.raises(AssertionError, match="Sequence/quality length mismatch"):
            trim_alignment_in_place(mock_read, 1, 1)

    def test_invalid_reference_start(self):
        """Test handling of invalid reference start values."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCG",
            cigartuples=[(0, 8)],
            reference_start=-5,  # Invalid negative position
        )

        with pytest.raises(AssertionError, match="Invalid reference_start.*must be non-negative"):
            trim_alignment_in_place(mock_read, 1, 1)

    def test_cigar_sequence_length_mismatch_detection(self):
        """Test detection of CIGAR/sequence length mismatches."""
        mock_read = MockAlignedSegment(
            query_sequence="ATCGATCGATCG",  # 12 bases
            cigartuples=[(0, 15)],  # 15M - MISMATCH with sequence length
        )

        # The current implementation validates this after trimming
        # This would likely cause an assertion error during trimming
        with pytest.raises(AssertionError, match="CIGAR/sequence mismatch"):
            trim_alignment_in_place(mock_read, 1, 1)


class TestPropertyBasedTesting:
    """Property-based tests using hypothesis."""

    @given(
        seq_length=st.integers(min_value=1, max_value=100),
        left_trim=st.integers(min_value=0, max_value=50),
        right_trim=st.integers(min_value=0, max_value=50),
    )
    @settings(max_examples=50, deadline=1000)
    def test_sequence_trimming_properties(self, seq_length, left_trim, right_trim):
        """Property-based test for sequence trimming invariants."""
        # Generate a test sequence
        sequence = "A" * seq_length

        mock_read = MockAlignedSegment(
            query_sequence=sequence,
            query_qualities=[30] * seq_length,
            cigartuples=[(0, seq_length)],  # All matches
        )

        trim_alignment_in_place(mock_read, left_trim, right_trim)

        # Property 1: Result length should be non-negative
        result_length = len(mock_read.query_sequence) if mock_read.query_sequence else 0
        assert result_length >= 0

        # Property 2: Result length should not exceed original length
        assert result_length <= seq_length

        # Property 3: If trim amounts are reasonable, check expected length
        expected_length = max(0, seq_length - left_trim - right_trim)
        assert result_length == expected_length

        # Property 4: Quality array should match sequence length (if present)
        if mock_read.query_qualities is not None:
            assert len(mock_read.query_qualities) == result_length

    @given(
        op_lengths=st.lists(
            st.integers(min_value=1, max_value=20),
            min_size=1,
            max_size=10,
        ),
        trim_amount=st.integers(min_value=0, max_value=50),
    )
    @settings(max_examples=30, deadline=1000)
    def test_cigar_consumption_properties(self, op_lengths, trim_amount):
        """Property-based test for CIGAR consumption invariants."""
        # Create a CIGAR with only match operations for simplicity
        cigar = Cigar([CigarOp(0, length) for length in op_lengths])
        original_query_len = sum(op_lengths)

        result_cigar, ref_advance = _consume_from_left(cigar, trim_amount)

        # Property 1: Reference advance should be non-negative and bounded
        assert 0 <= ref_advance <= original_query_len

        # Property 2: Result should be valid CIGAR
        assert isinstance(result_cigar, Cigar)

        # Property 3: No negative-length operations in result
        for op in result_cigar:
            assert op.length >= 0

        # Property 4: Total query consumption should be consistent
        consumed_query = min(trim_amount, original_query_len)
        remaining_query = sum(
            op.length for op in result_cigar if op.op in trim_aligned_reads.QRY_CONSUME
        )
        assert remaining_query + consumed_query == original_query_len

    @given(
        qname_prefix=st.sampled_from(["MERGED_", "UNMERGED_", "other_"]),
        qname_suffix=st.sampled_from(["/1", "/2", ""]),
        sequence_length=st.integers(min_value=10, max_value=100),
    )
    @settings(max_examples=30)
    def test_read_category_classification_properties(
        self,
        qname_prefix,
        qname_suffix,
        sequence_length,
    ):
        """Property-based test for read category classification."""
        qname = f"{qname_prefix}read_001{qname_suffix}"
        category = ReadCategory.classify(qname)

        # Property 1: Classification should be deterministic
        assert ReadCategory.classify(qname) == category

        # Property 2: Classification rules should be consistent
        if qname_prefix == "MERGED_":
            assert category == ReadCategory.MERGED
        elif qname_prefix == "UNMERGED_" and qname_suffix == "/1":
            assert category == ReadCategory.UNMERGED_R1
        elif qname_prefix == "UNMERGED_" and qname_suffix == "/2":
            assert category == ReadCategory.UNMERGED_R2
        else:
            assert category == ReadCategory.OTHER


# Marker for edge case tests
pytestmark = pytest.mark.unit
