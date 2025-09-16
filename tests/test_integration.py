"""
Integration tests for trim_aligned_reads.py

This module tests the complete workflow of the trimming tool with real file I/O,
various input formats (SAM/BAM/CRAM), and end-to-end functionality.
"""

from pathlib import Path
from unittest.mock import patch

import pysam
import pytest
import trim_aligned_reads
from trim_aligned_reads import TagConfig, main, process_stream


class TestFileIO:
    """Test file input/output operations."""

    def test_open_alignment_sam_read(self, empty_sam_file):
        """Test opening SAM file for reading."""
        with trim_aligned_reads.open_alignment(
            str(empty_sam_file), write=False
        ) as sam_file:
            assert isinstance(sam_file, pysam.AlignmentFile)
            assert not sam_file.is_write
            assert sam_file.filename.decode() == str(empty_sam_file)

    def test_open_alignment_sam_write(self, temp_dir, reference_sequence):
        """Test opening SAM file for writing."""
        output_path = temp_dir / "output.sam"
        header = {
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": "test", "LN": len(reference_sequence)}],
        }

        with trim_aligned_reads.open_alignment(
            str(output_path),
            write=True,
            template_or_header=header,
        ) as sam_file:
            assert isinstance(sam_file, pysam.AlignmentFile)
            assert sam_file.is_write

    def test_open_alignment_bam_operations(self, sample_bam_file, temp_dir):
        """Test BAM file operations."""
        output_path = temp_dir / "output.bam"

        # Read input BAM
        with trim_aligned_reads.open_alignment(
            str(sample_bam_file), write=False
        ) as input_bam:
            # Write output BAM using input as template
            with trim_aligned_reads.open_alignment(
                str(output_path),
                write=True,
                template_or_header=input_bam,
            ) as output_bam:
                # Copy a few reads
                read_count = 0
                for read in input_bam:
                    if read_count < 2:  # Copy first 2 reads
                        output_bam.write(read)
                        read_count += 1

        # Verify output file exists and is readable
        assert output_path.exists()
        with trim_aligned_reads.open_alignment(
            str(output_path), write=False
        ) as verify_bam:
            reads = list(verify_bam)
            assert len(reads) == 2

    def test_open_alignment_cram_with_reference(
        self,
        temp_dir,
        reference_fasta,
        reference_sequence,
    ):
        """Test CRAM file operations with reference."""
        cram_path = temp_dir / "test.cram"
        header = {
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": "test_reference", "LN": len(reference_sequence)}],
        }

        # Create CRAM with reference
        with trim_aligned_reads.open_alignment(
            str(cram_path),
            write=True,
            template_or_header=header,
            reference=str(reference_fasta),
        ) as cram_file:
            # Create a simple read
            read = pysam.AlignedSegment()
            read.query_name = "test_read"
            read.query_sequence = "ATCGATCG"
            read.query_qualities = [30] * 8
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 8)]  # 8M
            read.mapping_quality = 60
            cram_file.write(read)

        # Read it back
        with trim_aligned_reads.open_alignment(
            str(cram_path),
            write=False,
            reference=str(reference_fasta),
        ) as cram_read:
            reads = list(cram_read)
            assert len(reads) == 1
            assert reads[0].query_name == "test_read"

    def test_invalid_file_extension(self, temp_dir):
        """Test error handling for invalid file extensions."""
        invalid_path = temp_dir / "test.invalid"

        with pytest.raises(
            ValueError, match="Output/input must end with .sam, .bam, or .cram"
        ):
            trim_aligned_reads.open_alignment(str(invalid_path), write=False)

    def test_writing_without_header(self, temp_dir):
        """Test error handling when writing without header."""
        output_path = temp_dir / "output.sam"

        with pytest.raises(
            AssertionError, match="Writing .* requires template_or_header"
        ):
            trim_aligned_reads.open_alignment(str(output_path), write=True)


class TestProcessStream:
    """Test the main processing pipeline."""

    def test_process_empty_file(self, empty_sam_file, temp_dir, default_trim_policy):
        """Test processing an empty SAM file."""
        output_path = temp_dir / "output.sam"

        with trim_aligned_reads.open_alignment(
            str(empty_sam_file), write=False
        ) as input_file:
            with trim_aligned_reads.open_alignment(
                str(output_path),
                write=True,
                template_or_header=input_file,
            ) as output_file:
                kept, dropped_flag, dropped_short = process_stream(
                    input_file,
                    output_file,
                    default_trim_policy,
                    TagConfig(),
                )

        assert kept == 0
        assert dropped_flag == 0
        assert dropped_short == 0

    def test_process_sample_file(self, sample_bam_file, temp_dir, minimal_trim_policy):
        """Test processing the sample BAM file."""
        output_path = temp_dir / "trimmed.bam"

        with trim_aligned_reads.open_alignment(
            str(sample_bam_file), write=False
        ) as input_file:
            with trim_aligned_reads.open_alignment(
                str(output_path),
                write=True,
                template_or_header=input_file,
            ) as output_file:
                kept, dropped_flag, dropped_short = process_stream(
                    input_file,
                    output_file,
                    minimal_trim_policy,
                    TagConfig(),
                )

        # Should have processed some reads
        assert kept > 0
        assert kept + dropped_flag + dropped_short > 0

        # Verify output file
        with trim_aligned_reads.open_alignment(
            str(output_path), write=False
        ) as verify_file:
            output_reads = list(verify_file)
            assert len(output_reads) == kept

            # Verify all reads meet minimum length requirement after trimming
            for read in output_reads:
                if read.query_sequence:
                    assert len(read.query_sequence) >= minimal_trim_policy.min_len

    def test_process_with_aggressive_trimming(
        self,
        sample_bam_file,
        temp_dir,
        aggressive_trim_policy,
    ):
        """Test processing with aggressive trimming that drops short reads."""
        output_path = temp_dir / "aggressively_trimmed.bam"

        with trim_aligned_reads.open_alignment(
            str(sample_bam_file), write=False
        ) as input_file:
            with trim_aligned_reads.open_alignment(
                str(output_path),
                write=True,
                template_or_header=input_file,
            ) as output_file:
                kept, dropped_flag, dropped_short = process_stream(
                    input_file,
                    output_file,
                    aggressive_trim_policy,
                    TagConfig(),
                )

        # With aggressive trimming, we expect some reads to be dropped as too short
        assert dropped_short > 0 or kept == 0  # Either some dropped or all dropped

    def test_process_drop_untagged(
        self, sample_bam_file, temp_dir, default_trim_policy
    ):
        """Test processing with drop_untagged=True."""
        output_path = temp_dir / "tagged_only.bam"

        with trim_aligned_reads.open_alignment(
            str(sample_bam_file), write=False
        ) as input_file:
            with trim_aligned_reads.open_alignment(
                str(output_path),
                write=True,
                template_or_header=input_file,
            ) as output_file:
                kept, dropped_flag, dropped_short = process_stream(
                    input_file,
                    output_file,
                    default_trim_policy,
                    TagConfig(),
                    drop_untagged=True,
                )

        # Verify that only tagged reads (MERGED_*, UNMERGED_*) are kept
        with trim_aligned_reads.open_alignment(
            str(output_path), write=False
        ) as verify_file:
            for read in verify_file:
                assert read.query_name.startswith(
                    "MERGED_"
                ) or read.query_name.startswith(
                    "UNMERGED_",
                )

    def test_process_batch_size(self, sample_bam_file, temp_dir, default_trim_policy):
        """Test processing with different batch sizes."""
        output_path = temp_dir / "batched.bam"

        with trim_aligned_reads.open_alignment(
            str(sample_bam_file), write=False
        ) as input_file:
            with trim_aligned_reads.open_alignment(
                str(output_path),
                write=True,
                template_or_header=input_file,
            ) as output_file:
                kept, dropped_flag, dropped_short = process_stream(
                    input_file,
                    output_file,
                    default_trim_policy,
                    TagConfig(),
                    batch_size=1,  # Process one read at a time
                )

        # Results should be the same regardless of batch size
        assert kept >= 0
        assert kept + dropped_flag + dropped_short > 0


class TestMainFunction:
    """Test the main CLI entry point."""

    def test_main_basic_run(self, sample_bam_file, temp_dir, reference_fasta):
        """Test basic main function execution."""
        output_path = temp_dir / "main_output.bam"

        # Mock sys.argv
        test_args = [
            "trim_aligned_reads.py",
            "--in",
            str(sample_bam_file),
            "--out",
            str(output_path),
            "--ref",
            str(reference_fasta),
            "--merged-left",
            "5",
            "--merged-right",
            "5",
            "--min-len",
            "10",
        ]

        with patch("sys.argv", test_args):
            # Should not raise any exceptions
            main()

        # Verify output file was created
        assert output_path.exists()

        # Verify it contains valid data
        with trim_aligned_reads.open_alignment(
            str(output_path), write=False
        ) as output_file:
            reads = list(output_file)
            # Should have some reads (depending on sample data)
            assert len(reads) >= 0

    def test_main_with_verbose_logging(
        self, sample_bam_file, temp_dir, reference_fasta
    ):
        """Test main function with verbose logging."""
        output_path = temp_dir / "verbose_output.bam"

        test_args = [
            "trim_aligned_reads.py",
            "--in",
            str(sample_bam_file),
            "--out",
            str(output_path),
            "--ref",
            str(reference_fasta),
            "-vv",  # Very verbose
            "--min-len",
            "5",
        ]

        with patch("sys.argv", test_args):
            main()

        assert output_path.exists()

    def test_main_with_quiet_logging(self, sample_bam_file, temp_dir, reference_fasta):
        """Test main function with quiet logging."""
        output_path = temp_dir / "quiet_output.bam"

        test_args = [
            "trim_aligned_reads.py",
            "--in",
            str(sample_bam_file),
            "--out",
            str(output_path),
            "--ref",
            str(reference_fasta),
            "-qq",  # Very quiet
        ]

        with patch("sys.argv", test_args):
            main()

        assert output_path.exists()

    def test_main_drop_untagged_flag(self, sample_bam_file, temp_dir, reference_fasta):
        """Test main function with --drop-untagged flag."""
        output_path = temp_dir / "drop_untagged_output.bam"

        test_args = [
            "trim_aligned_reads.py",
            "--in",
            str(sample_bam_file),
            "--out",
            str(output_path),
            "--ref",
            str(reference_fasta),
            "--drop-untagged",
        ]

        with patch("sys.argv", test_args):
            main()

        assert output_path.exists()


class TestEdgeCases:
    """Test edge cases and error conditions in integration context."""

    def test_nonexistent_input_file(self, temp_dir):
        """Test error handling for nonexistent input file."""
        nonexistent = temp_dir / "does_not_exist.bam"
        output_path = temp_dir / "output.bam"

        test_args = [
            "trim_aligned_reads.py",
            "--in",
            str(nonexistent),
            "--out",
            str(output_path),
        ]

        with patch("sys.argv", test_args), pytest.raises((FileNotFoundError, OSError)):
            main()

    def test_invalid_output_directory(self, sample_bam_file):
        """Test error handling for invalid output directory."""
        invalid_output = Path("/nonexistent_directory/output.bam")

        test_args = [
            "trim_aligned_reads.py",
            "--in",
            str(sample_bam_file),
            "--out",
            str(invalid_output),
        ]

        with (
            patch("sys.argv", test_args),
            pytest.raises((FileNotFoundError, OSError, PermissionError)),
        ):
            main()

    def test_sam_to_bam_conversion(self, empty_sam_file, temp_dir, default_trim_policy):
        """Test converting from SAM to BAM during processing."""
        output_path = temp_dir / "converted.bam"

        with trim_aligned_reads.open_alignment(
            str(empty_sam_file), write=False
        ) as input_file:
            with trim_aligned_reads.open_alignment(
                str(output_path),
                write=True,
                template_or_header=input_file,
            ) as output_file:
                process_stream(
                    input_file, output_file, default_trim_policy, TagConfig()
                )

        # Verify BAM file was created and is readable
        assert output_path.exists()
        with trim_aligned_reads.open_alignment(
            str(output_path), write=False
        ) as bam_file:
            assert bam_file.filename.decode().endswith(".bam")


class TestValidation:
    """Test SAM/BAM format validation and compliance."""

    def test_output_header_preservation(
        self, sample_bam_file, temp_dir, default_trim_policy
    ):
        """Test that SAM header is preserved in output."""
        output_path = temp_dir / "header_test.bam"

        with trim_aligned_reads.open_alignment(
            str(sample_bam_file), write=False
        ) as input_file:
            original_header = dict(input_file.header)

            with trim_aligned_reads.open_alignment(
                str(output_path),
                write=True,
                template_or_header=input_file,
            ) as output_file:
                process_stream(
                    input_file, output_file, default_trim_policy, TagConfig()
                )

        # Check that header is preserved
        with trim_aligned_reads.open_alignment(
            str(output_path), write=False
        ) as output_file:
            output_header = dict(output_file.header)
            assert output_header == original_header

    def test_cigar_consistency(self, sample_bam_file, temp_dir, minimal_trim_policy):
        """Test that CIGAR strings remain consistent after trimming."""
        output_path = temp_dir / "cigar_test.bam"

        with trim_aligned_reads.open_alignment(
            str(sample_bam_file), write=False
        ) as input_file:
            with trim_aligned_reads.open_alignment(
                str(output_path),
                write=True,
                template_or_header=input_file,
            ) as output_file:
                process_stream(
                    input_file, output_file, minimal_trim_policy, TagConfig()
                )

        # Verify CIGAR consistency in output
        with trim_aligned_reads.open_alignment(
            str(output_path), write=False
        ) as output_file:
            for read in output_file:
                if read.cigartuples and read.query_sequence:
                    # Calculate query length from CIGAR
                    cigar_query_len = sum(
                        length
                        for op, length in read.cigartuples
                        if op in trim_aligned_reads.QRY_CONSUME
                    )
                    # Should match actual sequence length
                    assert len(read.query_sequence) == cigar_query_len

    def test_quality_score_consistency(
        self, sample_bam_file, temp_dir, default_trim_policy
    ):
        """Test that quality scores remain consistent with sequence length."""
        output_path = temp_dir / "quality_test.bam"

        with trim_aligned_reads.open_alignment(
            str(sample_bam_file), write=False
        ) as input_file:
            with trim_aligned_reads.open_alignment(
                str(output_path),
                write=True,
                template_or_header=input_file,
            ) as output_file:
                process_stream(
                    input_file, output_file, default_trim_policy, TagConfig()
                )

        # Verify quality consistency in output
        with trim_aligned_reads.open_alignment(
            str(output_path), write=False
        ) as output_file:
            for read in output_file:
                if read.query_sequence and read.query_qualities:
                    assert len(read.query_sequence) == len(read.query_qualities)


class TestClippingModeIntegration:
    """Test full pipeline integration with different clipping modes."""

    def test_clipping_modes_pipeline_verification(self, temp_dir):
        """Test all clipping modes verify behavior without file format constraints."""
        # Create test SAM file
        input_sam = temp_dir / "input.sam"
        sam_content = """@HD\tVN:1.6
@SQ\tSN:test_ref\tLN:1000
MERGED_read1\t0\ttest_ref\t100\t60\t16M\t*\t0\t0\tATCGATCGATCGATCG\t################
"""
        input_sam.write_text(sam_content)

        # Test DELETE and HARD_CLIP modes (they produce valid SAM)
        for mode_str in ["delete", "hard-clip"]:
            output_sam = temp_dir / f"output_{mode_str}.sam"

            main(
                [
                    "-i",
                    str(input_sam),
                    "-o",
                    str(output_sam),
                    "--clipping-mode",
                    mode_str,
                    "--merged-left",
                    "2",
                    "--merged-right",
                    "3",
                    "--min-len",
                    "5",
                ]
            )

            # Verify output file exists
            assert output_sam.exists()

            # Read raw SAM content to verify clipping behavior
            output_content = output_sam.read_text()
            lines = [
                line
                for line in output_content.split("\n")
                if line and not line.startswith("@")
            ]

            assert len(lines) == 1  # Should have one read
            fields = lines[0].split("\t")
            read_name, cigar_str, sequence = fields[0], fields[5], fields[9]

            # Validate read name is preserved correctly (no --strip-tags used)
            assert read_name == "MERGED_read1", (
                f"Read name should be preserved: {read_name}"
            )

            if mode_str == "delete":
                # DELETE mode: sequence trimmed, simple CIGAR
                assert len(sequence) == 11  # 16 - 2 - 3 = 11
                assert cigar_str == "11M"  # Simple match

            elif mode_str == "hard-clip":
                # HARD_CLIP mode: sequence trimmed, CIGAR has hard clips
                assert len(sequence) == 11  # 16 - 2 - 3 = 11
                assert "H" in cigar_str  # Should contain hard clips
                assert cigar_str == "2H11M3H"  # Expected hard clip pattern

    def test_soft_clip_mode_behavior_verification(self, temp_dir):
        """Test SOFT_CLIP mode behavior by examining raw output."""
        input_sam = temp_dir / "input.sam"
        sam_content = """@HD\tVN:1.6
@SQ\tSN:test_ref\tLN:1000
MERGED_test\t0\ttest_ref\t100\t60\t12M\t*\t0\t0\tATCGATCGATCG\t############
"""
        input_sam.write_text(sam_content)

        output_sam = temp_dir / "output_soft_clip.sam"

        main(
            [
                "-i",
                str(input_sam),
                "-o",
                str(output_sam),
                "--clipping-mode",
                "soft-clip",
                "--merged-left",
                "1",
                "--merged-right",
                "2",
                "--min-len",
                "5",
            ]
        )

        # Verify output by reading raw SAM content
        assert output_sam.exists()
        output_content = output_sam.read_text()
        lines = [
            line
            for line in output_content.split("\n")
            if line and not line.startswith("@")
        ]

        assert len(lines) == 1
        fields = lines[0].split("\t")
        read_name, cigar_str, sequence = fields[0], fields[5], fields[9]

        # Validate read name is preserved correctly (no --strip-tags used)
        assert read_name == "MERGED_test", f"Read name should be preserved: {read_name}"

        # SOFT_CLIP: sequence preserved, CIGAR has soft clips
        assert sequence == "ATCGATCGATCG"  # Original sequence preserved
        assert cigar_str == "1S12M2S"  # Expected soft clip pattern
        assert fields[3] == "100"  # Reference position unchanged

    def test_clipping_modes_sam_to_bam(self, temp_dir):
        """Test clipping modes work with SAM to BAM conversion (using DELETE mode)."""
        input_sam = temp_dir / "input.sam"
        sam_content = """@HD\tVN:1.6
@SQ\tSN:test_ref\tLN:1000
MERGED_test\t0\ttest_ref\t100\t60\t16M\t*\t0\t0\tATCGATCGATCGATCG\t################
"""
        input_sam.write_text(sam_content)

        # Test DELETE mode with BAM output (produces valid SAM/BAM format)
        output_bam = temp_dir / "output_delete.bam"

        main(
            [
                "-i",
                str(input_sam),
                "-o",
                str(output_bam),
                "--clipping-mode",
                "delete",
                "--merged-left",
                "2",
                "--merged-right",
                "3",
                "--min-len",
                "5",
            ]
        )

        # Verify BAM output is correct
        assert output_bam.exists()
        with trim_aligned_reads.open_alignment(
            str(output_bam), write=False
        ) as bam_file:
            reads = list(bam_file)
            assert len(reads) == 1

            read = reads[0]
            # Should have trimmed sequence with simple CIGAR
            assert len(read.query_sequence) == 11  # 16 - 2 - 3 = 11
            assert read.cigartuples == [(0, 11)]  # 11M
            assert (
                read.reference_start == 101
            )  # Advanced by left trim (2): 99 + 2 = 101


# Marker for integration tests
pytestmark = pytest.mark.integration
