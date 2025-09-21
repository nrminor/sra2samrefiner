#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
#     "pysam",
# ]
# ///

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from enum import Enum, auto
from typing import TYPE_CHECKING, NamedTuple

import pysam
from loguru import logger

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Sequence

# ------------------------------- CONSTANTS -------------------------------- #

# CIGAR op codes
# 0:M, 1:I, 2:D, 3:N, 4:S, 5:H, 6:P, 7:=, 8:X
REF_CONSUME = {0, 2, 3, 7, 8}
QRY_CONSUME = {0, 1, 4, 7, 8}
BOTH_CONSUME = {0, 7, 8}

# Emit a progress debug line after processing this many primary reads
DEBUG_EVERY: int = 100_000


# ------------------------------- DATA TYPES -------------------------------- #


@dataclass(frozen=True)
class TagConfig:
    """Configuration for read name tagging and classification."""

    merged_prefix: str = "MERGED_"
    unmerged_prefix: str = "UNMERGED_"
    read1_suffix: str = "/1"
    read2_suffix: str = "/2"

    def strip_tags(self, qname: str) -> str:
        """Remove prefixes from read names."""
        if qname.startswith(self.merged_prefix):
            return qname[len(self.merged_prefix) :]
        if qname.startswith(self.unmerged_prefix):
            return qname[len(self.unmerged_prefix) :]
        return qname


class ClippingMode(Enum):
    """Defines how trimmed bases should be handled."""

    DELETE = auto()  # Current behavior - remove completely
    SOFT_CLIP = auto()  # Convert to soft clips (S operations)
    HARD_CLIP = auto()  # Convert to hard clips (H operations)


@dataclass(frozen=True)
class TrimPolicy:
    """
    Trim amounts in query-bases for each read category and clipping behavior.

    Read categories:
    - merged_*: Reads successfully merged by upstream tools (e.g., bbmerge)
    - r1_*/r2_*: Unmerged paired-end reads (R1 gets left trim, R2 gets right trim)
    - single_*: Single-end reads or untagged reads

    Clipping modes:
    - DELETE: Remove trimmed bases completely (default, most compatible)
    - SOFT_CLIP: Mark trimmed regions as soft clips (S), preserve sequence
    - HARD_CLIP: Mark trimmed regions as hard clips (H), remove sequence
    """

    merged_left: int = 30  # Left trim for merged reads
    merged_right: int = 30  # Right trim for merged reads
    r1_left: int = 30  # Left trim for unmerged R1 reads
    r2_right: int = 30  # Right trim for unmerged R2 reads
    single_left: int = 30  # Left trim for single-end/untagged reads
    single_right: int = 30  # Right trim for single-end/untagged reads
    min_len: int = 20  # Drop reads shorter than this AFTER trimming
    clipping_mode: ClippingMode = ClippingMode.DELETE  # How to handle trimmed bases


class ReadCategory(Enum):
    """Closed set of categories derived from read names, with trim logic."""

    MERGED = auto()
    UNMERGED_R1 = auto()
    UNMERGED_R2 = auto()
    OTHER = auto()

    @staticmethod
    def classify(qname: str, tag_config: TagConfig) -> ReadCategory:
        """Classify by upstream naming scheme."""
        if qname.startswith(tag_config.merged_prefix):
            return ReadCategory.MERGED
        if qname.startswith(tag_config.unmerged_prefix) and qname.endswith(
            tag_config.read1_suffix
        ):
            return ReadCategory.UNMERGED_R1
        if qname.startswith(tag_config.unmerged_prefix) and qname.endswith(
            tag_config.read2_suffix
        ):
            return ReadCategory.UNMERGED_R2
        return ReadCategory.OTHER

    def trim_extents(self, policy: TrimPolicy) -> tuple[int, int]:
        """Return (left_trim, right_trim) for this category."""
        match self:
            case ReadCategory.MERGED:
                return max(0, policy.merged_left), max(0, policy.merged_right)
            case ReadCategory.UNMERGED_R1:
                return max(0, policy.r1_left), 0
            case ReadCategory.UNMERGED_R2:
                return max(0, policy.r2_right), 0
            case ReadCategory.OTHER:
                # Apply single-end trimming for untagged reads (single-end, Nanopore, etc.)
                return max(0, policy.single_left), max(0, policy.single_right)


# ----------------------------- LOGGING SETUP ------------------------------- #


def configure_logging(verbose: int, quiet: int) -> None:
    """
    Base at SUCCESS (0). Positive → louder (more verbose), negative → quieter.
    Map:
      +3.. = TRACE
      +2   = DEBUG
      +1   = INFO
       0   = SUCCESS
      -1   = WARNING
      -2   = ERROR
      <=-3 = CRITICAL
    """
    logger.remove()
    delta = verbose - quiet
    match delta:
        case d if d >= 3:  # noqa: PLR2004
            level_str = "TRACE"
        case 2:
            level_str = "DEBUG"
        case 1:
            level_str = "INFO"
        case 0:
            level_str = "SUCCESS"
        case -1:
            level_str = "WARNING"
        case -2:
            level_str = "ERROR"
        case d if d <= -3:  # noqa: PLR2004
            level_str = "CRITICAL"
    logger.add(sys.stderr, level=level_str)
    logger.debug(f"Logger configured at level: {level_str}")


# ---------------------------- CIGAR UTILITIES ------------------------------ #


class CigarOp(NamedTuple):
    """One CIGAR run: (operation code, run length)."""

    op: int
    length: int

    @staticmethod
    def from_tuple(t: tuple[int, int]) -> CigarOp:
        """Convert a raw (op, len) tuple to CigarOp."""
        op, ln = t
        return CigarOp(op, ln)

    @staticmethod
    def to_tuple(run: CigarOp) -> tuple[int, int]:
        """Convert a CigarOp back to a raw (op, len) tuple."""
        return (run.op, run.length)


class Cigar(list[CigarOp]):
    """A list of CigarOp with helpers for conversion and compaction."""

    @classmethod
    def from_pysam(cls, cig_raw: list[tuple[int, int]] | None) -> Cigar | None:
        """
        Convert pysam's list[(op, len)] to a Cigar. Returns None if input is None.
        """
        if cig_raw is None:
            return None
        return cls(CigarOp.from_tuple(t) for t in cig_raw)

    def to_pysam(self) -> list[tuple[int, int]]:
        """Convert this Cigar back to list[(op, len)] for pysam."""
        return [CigarOp.to_tuple(run) for run in self]

    def push_compact(self, op: int, ln: int) -> None:
        """
        Append (op, ln), merging with the last run if `op` matches.
        Ignores non-positive lengths.
        """
        # Positive invariant: operation code must be valid CIGAR operation (0-8)
        assert 0 <= op <= 8, (
            f"Invalid CIGAR operation code {op}: must be 0-8 (M,I,D,N,S,H,P,=,X)"
        )  # noqa: PLR2004

        # Negative invariant: length must not cause integer overflow when added
        if self and self[-1].op == op:
            last = self[-1]
            assert last.length + ln <= (1 << 28) - 1, (
                f"CIGAR length overflow: {last.length} + {ln} exceeds maximum safe integer"
            )

        if ln <= 0:
            return
        if self and self[-1].op == op:
            last = self[-1]
            self[-1] = CigarOp(op, last.length + ln)
            return
        self.append(CigarOp(op, ln))


def _consume_from_left(cig: Cigar, trim_q: int) -> tuple[Cigar, int]:
    """
    Consume `trim_q` query-bases from the LEFT edge of `cig`.

    Returns
    -------
    new_cigar : Cigar
        The CIGAR after trimming from the left (already compacted).
    ref_advance : int
        How many reference bases were consumed from {M,=,X} while trimming.
        Add this to `reference_start`.
    """
    # Positive invariant: trim amount must be non-negative
    assert trim_q >= 0, f"Left trim amount must be non-negative, got {trim_q}"

    # Negative invariant: CIGAR must contain valid operations only
    if cig:
        assert all(0 <= run.op <= 8 and run.length > 0 for run in cig), (  # noqa: PLR2004
            f"CIGAR contains invalid operations: {[(run.op, run.length) for run in cig]}"
        )

    if trim_q <= 0 or not cig:
        return cig, 0

    # Calculate total query bases available for validation
    total_query_bases = sum(run.length for run in cig if run.op in QRY_CONSUME)

    out = Cigar()
    remaining = trim_q
    ref_advance = 0
    total_consumed_query = 0
    i = 0
    '''
    3 operations:
    trim run, advance start, reduce remaining trim len
    trim : all
    advance: REF_CONSUME
    reduce: QRY_CONSUME
    
    M: trim, advance, reduce
    I: trim, reduce
    D: trim, advance
    N: trim, advance
    S: trim, reduce
    H: trim
    P: trim
    =: trim, advance, reduce
    X: trim, advance, reduce

    4 cases:
    BOTH_CONSUME: trim, advance, reduce
    REF_CONSUME and NOT QRY_CONSUME: trim, advance
    QRY_CONSUME and NOT REF_CONSUME: trim, reduce
    NOT REF_CONSUME or QRY_CONSUME: trim
    '''
    while i < len(cig) and remaining > 0:
        run = cig[i]
        i += 1
        # H/P: trimmed from CIGAR, no further actions needed
        if not (run.op in REF_CONSUME or run.op in QRY_CONSUME):
            continue

        take = min(run.length, remaining)
        # advance start position for REF_CONSUME: M/D/N/=/X
        if run.op in REF_CONSUME:
            if run.op in BOTH_CONSUME:
                ref_advance += take
            else:
                ref_advance += run.length
        # reduce remaining for QRY_CONSUME: M/I/S/=/X
        # remove run otherwise
        if not run.op in QRY_CONSUME:
            continue
        keep_len = run.length - take
        remaining -= take
        total_consumed_query += take
        if keep_len > 0: # only keep if run still exists
            out.push_compact(run.op, keep_len)

    # Append untouched tail
    for j in range(i, len(cig)):
        tail = cig[j]
        out.push_compact(tail.op, tail.length)

    # Positive invariant: reference advance must be non-negative and bounded
    assert ref_advance >= 0, f"Reference advance cannot be negative: {ref_advance}"
    assert total_consumed_query <= min(trim_q, total_query_bases), (
        f"Consumed {total_consumed_query} query bases but trim_q={trim_q}, total_available={total_query_bases}"
    )

    # Negative invariant: output CIGAR must not have zero-length operations
    if out:
        assert all(run.length > 0 for run in out), (
            f"Output CIGAR contains zero-length operations: {[(run.op, run.length) for run in out]}"
        )

    return out, ref_advance


def _consume_from_right(cig: Cigar, trim_q: int) -> Cigar:
    """
    Consume `trim_q` query-bases from the RIGHT edge of `cig`.
    Reference start is unaffected by definition.
    """
    # Positive invariant: trim amount must be non-negative
    assert trim_q >= 0, f"Right trim amount must be non-negative, got {trim_q}"

    # Negative invariant: CIGAR must contain valid operations only
    if cig:
        assert all(0 <= run.op <= 8 and run.length > 0 for run in cig), (  # noqa: PLR2004
            f"CIGAR contains invalid operations: {[(run.op, run.length) for run in cig]}"
        )

    if trim_q <= 0 or not cig:
        return cig

    # Calculate total query bases available for validation
    total_query_bases = sum(run.length for run in cig if run.op in QRY_CONSUME)

    out_rev = Cigar()
    remaining = trim_q
    total_consumed_query = 0
    '''
    2 operations:
    trim run, reduce remaining trim len
    trim : all
    reduce: QRY_CONSUME
    
    M: trim, reduce
    I: trim, reduce
    D: trim
    N: trim
    S: trim, reduce
    H: trim
    P: trim
    =: trim reduce
    X: trim, reduce

    2 cases:
    QRY_CONSUME: trim, reduce
    NOT QRY_CONSUME: trim
    '''
    # Build reversed output with compaction
    for run in reversed(cig):
        if remaining > 0:
            if run.op in QRY_CONSUME:
                take = min(run.length, remaining)
                keep_len = run.length - take
                remaining -= take
                total_consumed_query += take
                if keep_len > 0:
                    out_rev.push_compact(run.op, keep_len)
        else:
            out_rev.push_compact(run.op, run.length)

    # Reverse back (preserving compaction)
    out = Cigar()
    for run in reversed(out_rev):
        out.push_compact(run.op, run.length)

    # Positive invariant: consumed bases must not exceed available or requested
    assert total_consumed_query <= min(trim_q, total_query_bases), (
        f"Right trim consumed {total_consumed_query} query bases but trim_q={trim_q}, total_available={total_query_bases}"
    )

    # Negative invariant: output CIGAR must not have zero-length operations
    if out:
        assert all(run.length > 0 for run in out), (
            f"Output CIGAR contains zero-length operations: {[(run.op, run.length) for run in out]}"
        )

    return out


def _convert_to_soft_clips(cig: Cigar, left_trim: int, right_trim: int) -> Cigar:
    """
    Convert trim amounts to soft clips (S operations) at CIGAR ends.

    Args:
        cig: Original CIGAR operations
        left_trim: Bases to soft clip from left (query coordinates)
        right_trim: Bases to soft clip from right (query coordinates)

    Returns:
        New CIGAR with soft clips added at appropriate ends

    Note:
        Soft clips (S) consume query but not reference space.
    """
    # Input validation
    assert left_trim >= 0, f"Left trim must be non-negative: {left_trim}"
    assert right_trim >= 0, f"Right trim must be non-negative: {right_trim}"

    out = Cigar()

    # Add left soft clip if needed
    if left_trim > 0:
        out.push_compact(4, left_trim)  # S operation

    # Copy existing CIGAR operations
    for op in cig:
        out.push_compact(op.op, op.length)

    # Add right soft clip if needed
    if right_trim > 0:
        out.push_compact(4, right_trim)  # S operation

    return out


def _convert_to_hard_clips(cig: Cigar, left_trim: int, right_trim: int) -> Cigar:
    """
    Convert trim amounts to hard clips (H operations) at CIGAR ends.

    Args:
        cig: Original CIGAR operations
        left_trim: Bases to hard clip from left (query coordinates)
        right_trim: Bases to hard clip from right (query coordinates)

    Returns:
        New CIGAR with hard clips added at appropriate ends

    Note:
        Hard clips (H) consume neither query nor reference space.
    """
    # Input validation
    assert left_trim >= 0, f"Left trim must be non-negative: {left_trim}"
    assert right_trim >= 0, f"Right trim must be non-negative: {right_trim}"

    out = Cigar()

    # Add left hard clip if needed
    if left_trim > 0:
        out.push_compact(5, left_trim)  # H operation

    # Copy existing CIGAR operations
    for op in cig:
        out.push_compact(op.op, op.length)

    # Add right hard clip if needed
    if right_trim > 0:
        out.push_compact(5, right_trim)  # H operation

    return out


def trim_alignment_in_place(  # noqa: C901, PLR0912, PLR0915
    aln: pysam.AlignedSegment,
    left: int,
    right: int,
    mode: ClippingMode = ClippingMode.DELETE,
) -> None:
    """
    Trim query-bases from the read in-place and update sequence, qualities, CIGAR,
    and reference_start consistently. Strand-aware:
      - For reverse reads, the 5' end is on the RIGHT of the CIGAR.
      - Sequence slicing is always in read orientation.

    Args:
        aln: The aligned segment to trim in-place
        left: Bases to trim from left (read orientation)
        right: Bases to trim from right (read orientation)
        mode: How to handle trimmed bases (DELETE/SOFT_CLIP/HARD_CLIP)
    """
    # Positive invariant: trim amounts are normalized to non-negative
    left = max(left, 0)
    right = max(right, 0)
    assert left >= 0 and right >= 0, (  # noqa: PT018
        f"Normalized trim amounts must be non-negative: left={left}, right={right}"
    )

    seq: str | None = aln.query_sequence
    if seq is None:
        logger.debug(
            f"Skip trim: the read '{aln.query_name}' has no query_sequence (likely unmapped).",
        )
        return

    # Positive invariant: sequence and qualities must have consistent lengths
    qual = aln.query_qualities  # array('B') or None
    if qual is not None:
        assert len(qual) == len(seq), (
            f"Sequence/quality length mismatch for '{aln.query_name}': seq={len(seq)}, qual={len(qual)}"
        )

    cig_raw = aln.cigartuples
    if cig_raw is None:
        # Handle reads without CIGAR (unmapped) - use same safe boundary logic as main path
        seq_len = len(seq)
        cut_left = min(left, seq_len)
        cut_right = min(right, max(seq_len - cut_left, 0))
        keep_end = seq_len - cut_left - cut_right

        # Apply safe slicing with boundary checking
        new_seq = "" if keep_end <= 0 else seq[cut_left : cut_left + keep_end]
        new_qual = (
            None
            if qual is None
            else (None if keep_end <= 0 else qual[cut_left : cut_left + keep_end])
        )

        # Negative invariant: trimming cannot result in negative-length sequences
        expected_new_len = max(0, seq_len - left - right)
        assert len(new_seq) == expected_new_len, (
            f"Sequence trimming error for '{aln.query_name}': expected {expected_new_len}, got {len(new_seq)}"
        )

        aln.query_sequence = new_seq
        aln.query_qualities = new_qual
        logger.debug(
            f"Trimmed sequence/qualities only (no CIGAR) for '{aln.query_name}'.",
        )
        return

    cig: Cigar = Cigar.from_pysam(cig_raw)
    ref_start_before = aln.reference_start

    # Positive invariant: reference start must be non-negative if set
    if ref_start_before is not None:
        assert ref_start_before >= 0, (
            f"Invalid reference_start for '{aln.query_name}': {ref_start_before} (must be non-negative)"
        )

    # Strand-aware CIGAR trimming amounts
    cig_left = right if aln.is_reverse else left
    cig_right = left if aln.is_reverse else right

    if cig_left or cig_right:
        logger.debug(
            f"CIGAR trim for '{aln.query_name}': is_rev={aln.is_reverse}, left={cig_left}, right={cig_right}, "
            f"ref_start_before={ref_start_before}, mode={mode.name}",
        )

    # Mode-specific trimming behavior
    match mode:
        case ClippingMode.DELETE:
            # Current behavior - consume from CIGAR and trim sequence
            if cig_left > 0:
                cig, ref_advance = _consume_from_left(cig, cig_left)
                if ref_start_before is not None:
                    new_ref_start = ref_start_before + ref_advance
                    # Negative invariant: new reference start cannot be negative
                    assert new_ref_start >= 0, (
                        f"Reference start underflow for '{aln.query_name}': {ref_start_before} + {ref_advance} = {new_ref_start}"
                    )
                    aln.reference_start = new_ref_start
                # If ref_start_before is None, leave it as None (unmapped read)
            if cig_right > 0:
                cig = _consume_from_right(cig, cig_right)

            # Trim sequence/qualities in read orientation
            qlen = len(seq)
            cut_left = min(cig_left, qlen)
            cut_right = min(cig_right, max(qlen - cut_left, 0))
            keep_end = qlen - cut_left - cut_right

            # Positive invariant: trimming calculations must be arithmetically sound
            assert cut_left >= 0 and cut_right >= 0 and keep_end >= 0, (  # noqa: PT018
                f"Invalid trim calculations for '{aln.query_name}': cut_left={cut_left}, cut_right={cut_right}, keep_end={keep_end}"
            )
            assert cut_left + cut_right + keep_end == qlen, (
                f"Trim arithmetic error for '{aln.query_name}': {cut_left} + {cut_right} + {keep_end} != {qlen}"
            )

            new_seq = "" if keep_end <= 0 else seq[cut_left : cut_left + keep_end]
            new_qual = None if qual is None else qual[cut_left : cut_left + keep_end]

        case ClippingMode.SOFT_CLIP:
            # Add soft clips, keep full sequence/qualities
            cig = _convert_to_soft_clips(cig, cig_left, cig_right)
            new_seq = seq  # Keep full sequence
            new_qual = qual  # Keep full qualities
            # reference_start unchanged (soft clips don't consume reference)

        case ClippingMode.HARD_CLIP:
            # Add hard clips, but first consume CIGAR like DELETE mode to match trimmed sequence
            if cig_left > 0:
                cig, ref_advance = _consume_from_left(cig, cig_left)
                # NOTE: ref_start is NOT updated for hard clips (they don't consume reference)
            if cig_right > 0:
                cig = _consume_from_right(cig, cig_right)

            # Now add hard clips to the consumed CIGAR
            cig = _convert_to_hard_clips(cig, cig_left, cig_right)

            # Trim sequence/qualities in read orientation (same as DELETE)
            qlen = len(seq)
            cut_left = min(cig_left, qlen)
            cut_right = min(cig_right, max(qlen - cut_left, 0))
            keep_end = qlen - cut_left - cut_right

            # Same validation as DELETE mode
            assert cut_left >= 0 and cut_right >= 0 and keep_end >= 0, (  # noqa: PT018
                f"Invalid trim calculations for '{aln.query_name}': cut_left={cut_left}, cut_right={cut_right}, keep_end={keep_end}"
            )
            assert cut_left + cut_right + keep_end == qlen, (
                f"Trim arithmetic error for '{aln.query_name}': {cut_left} + {cut_right} + {keep_end} != {qlen}"
            )

            new_seq = "" if keep_end <= 0 else seq[cut_left : cut_left + keep_end]
            new_qual = None if qual is None else qual[cut_left : cut_left + keep_end]
            # reference_start unchanged (hard clips don't consume reference)

    # Verify CIGAR/sequence consistency after trimming (mode-dependent)
    new_cigar_query_len = sum(op[1] for op in cig.to_pysam() if op[0] in QRY_CONSUME)

    # Mode-specific validation of CIGAR/sequence consistency
    match mode:
        case ClippingMode.DELETE:
            # For DELETE: sequence length should match CIGAR query consumption
            assert len(new_seq) == new_cigar_query_len, (
                f"CIGAR/sequence mismatch after {mode.name} trimming '{aln.query_name}': "
                f"seq_len={len(new_seq)}, cigar_query_len={new_cigar_query_len}"
            )
        case ClippingMode.HARD_CLIP:
            # For HARD_CLIP: sequence trimmed, but CIGAR query consumption only from non-H operations
            # Hard clips (H) don't consume query, so cigar_query_len should equal trimmed sequence length
            assert len(new_seq) == new_cigar_query_len, (
                f"HARD_CLIP sequence/CIGAR mismatch for '{aln.query_name}': "
                f"seq_len={len(new_seq)}, non_hard_clip_query_len={new_cigar_query_len}"
            )
        case ClippingMode.SOFT_CLIP:
            # For SOFT_CLIP: sequence unchanged, CIGAR query consumption includes original + soft clips
            assert len(new_seq) == len(seq), (
                f"SOFT_CLIP should preserve sequence length for '{aln.query_name}': "
                f"original={len(seq)}, new={len(new_seq)}"
            )
            # For soft clips: total query consumption = original CIGAR consumption + soft clips
            original_cigar_query_len = sum(
                op[1] for op in cig_raw if op[0] in QRY_CONSUME
            )
            expected_cigar_query_len = original_cigar_query_len + cig_left + cig_right
            assert new_cigar_query_len == expected_cigar_query_len, (
                f"SOFT_CLIP CIGAR query consumption should be original_cigar + trims for '{aln.query_name}': "
                f"expected={expected_cigar_query_len} (orig_cigar={original_cigar_query_len}+left={cig_left}+right={cig_right}), "
                f"actual={new_cigar_query_len}"
            )

    # Positive invariant: quality array must match sequence length if present
    if new_qual is not None:
        assert len(new_qual) == len(new_seq), (
            f"Quality/sequence length mismatch after trimming '{aln.query_name}': qual={len(new_qual)}, seq={len(new_seq)}"
        )

    aln.cigartuples = cig.to_pysam()
    aln.query_sequence = new_seq
    aln.query_qualities = new_qual

    logger.debug(
        f"Trim complete: kept={len(new_seq)} bases; "
        f"ref_start {ref_start_before} -> {aln.reference_start}",
    )


# ----------------------------- I/O UTILITIES ------------------------------- #


def _io_mode_from_ext(path: str, write: bool) -> str:  # noqa: FBT001
    """Determine pysam open mode from filename extension."""
    lower = path.lower()
    if lower.endswith(".sam"):
        return "w" if write else "r"
    if lower.endswith(".bam"):
        return "wb" if write else "rb"
    if lower.endswith(".cram"):
        return "wc" if write else "rc"
    msg = "Output/input must end with .sam, .bam, or .cram"
    logger.error(msg)
    raise ValueError(msg)


def open_alignment(
    path: str,
    write: bool,  # noqa: FBT001
    template_or_header: pysam.AlignmentFile | dict | None = None,
    reference: str | None = None,
) -> pysam.AlignmentFile:
    """
    Open SAM/BAM/CRAM with correct mode. For CRAM, pass a reference filename.
    - If write=True and template_or_header is an AlignmentFile, we use 'template=...'
      to preserve header (lossless).
    - Otherwise, pass a header dict.
    """
    # Positive invariant: path must be a non-empty string
    assert isinstance(path, str) and len(path) > 0, (
        f"Path must be non-empty string, got: {path!r}"
    )  # noqa: PT018

    # Negative invariant: reference path must be valid if provided
    if reference is not None:
        assert isinstance(reference, str) and len(reference) > 0, (  # noqa: PT018
            f"Reference path must be non-empty string if provided, got: {reference!r}"
        )

    mode = _io_mode_from_ext(path, write)

    # Positive invariant: mode must be valid pysam mode string
    valid_modes = {"r", "rb", "rc", "w", "wb", "wc"}
    assert mode in valid_modes, f"Invalid pysam mode '{mode}' for path '{path}'"

    kwargs = {}
    if path.lower().endswith(".cram") and reference is None:
        logger.warning(
            f"Opening CRAM without explicit reference: {path}. "
            "Decoding may fail unless the reference is resolvable.",
        )
    if path.lower().endswith(".cram") and reference is not None:
        kwargs["reference_filename"] = reference

    action = "write" if write else "read"
    logger.debug(f"Opening for {action}: {path} (mode={mode})")
    if write:
        # Negative invariant: writing requires proper template or header
        assert template_or_header is not None, (
            f"Writing to '{path}' requires template_or_header but got None"
        )

        if isinstance(template_or_header, pysam.AlignmentFile):
            return pysam.AlignmentFile(
                path,
                mode,
                template=template_or_header,
                **kwargs,
            )
        if isinstance(template_or_header, dict):
            return pysam.AlignmentFile(path, mode, header=template_or_header, **kwargs)
        msg = f"Writing requires either a template AlignmentFile or a header dict, got {type(template_or_header)}"
        logger.error(msg)
        raise ValueError(msg)
    return pysam.AlignmentFile(path, mode, **kwargs)


def batched(
    iterable: Iterable[pysam.AlignedSegment],
    batch_size: int,
) -> Iterator[list[pysam.AlignedSegment]]:
    """
    Yield lists of reads up to batch_size. Keeps memory bounded and provides
    a simple place to insert batch-wise operations if ever needed.
    """
    batch: list[pysam.AlignedSegment] = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= batch_size:
            # Hand batch to caller
            yield batch

            # Caller has finished when we resume here; drop refs eagerly
            batch.clear()
    if batch:
        yield batch
        batch.clear()


# ------------------------------ CORE LOGIC --------------------------------- #


def process_stream(  # noqa: C901, PLR0912, PLR0913, PLR0915
    inp: pysam.AlignmentFile,
    outp: pysam.AlignmentFile,
    policy: TrimPolicy,
    tag_config: TagConfig,
    batch_size: int = 10000,
    drop_untagged: bool = False,  # noqa: FBT001, FBT002
    strip_tags: bool = False,  # noqa: FBT001, FBT002
) -> tuple[int, int, int]:
    """
    Stream input -> output in batches, trimming per read-name class.

    Processing behavior:
    - Drops unmapped, secondary, and supplementary alignments
    - Classifies reads by tag prefixes (MERGED_*/UNMERGED_*/OTHER)
    - Applies category-specific trimming based on policy
    - Optionally drops untagged reads (if drop_untagged=True)
    - Optionally strips tag prefixes from read names (if strip_tags=True)

    Clipping modes (from policy.clipping_mode):
    - DELETE: Remove trimmed bases completely (default, most compatible)
    - SOFT_CLIP: Mark trimmed regions as soft clips (S), preserve sequence
    - HARD_CLIP: Mark trimmed regions as hard clips (H), remove sequence

    Args:
        inp: Input alignment file
        outp: Output alignment file
        policy: Trimming policy including clipping mode
        tag_config: Tag configuration for read classification
        batch_size: Number of reads to process per batch
        drop_untagged: If True, drop reads without tag prefixes
        strip_tags: If True, remove tag prefixes from output read names

    Returns:
        Tuple of (kept_reads, dropped_unmapped_or_nonprimary, dropped_too_short)
    """
    # Positive invariant: batch size must be positive
    assert batch_size > 0, f"Batch size must be positive, got {batch_size}"

    # Extract clipping mode from policy
    clipping_mode = policy.clipping_mode

    # Positive invariant: policy must have valid trim values
    assert policy.min_len >= 0, (
        f"Policy min_len must be non-negative, got {policy.min_len}"
    )
    assert policy.merged_left >= 0 and policy.merged_right >= 0, (  # noqa: PT018
        f"Policy merged trim values must be non-negative: left={policy.merged_left}, right={policy.merged_right}"
    )
    assert policy.r1_left >= 0 and policy.r2_right >= 0, (  # noqa: PT018
        f"Policy unmerged trim values must be non-negative: r1_left={policy.r1_left}, r2_right={policy.r2_right}"
    )
    assert policy.single_left >= 0 and policy.single_right >= 0, (  # noqa: PT018
        f"Policy single-end trim values must be non-negative: single_left={policy.single_left}, single_right={policy.single_right}"
    )

    kept = 0
    dropped_flag = 0  # unmapped / secondary / supplementary
    dropped_short = 0  # too short after (or predicted) trimming
    dropped_primary_other = 0  # primary dropped due to --drop-untagged
    seen_primary = 0

    for batch in batched(inp, batch_size):
        logger.debug(f"Processing batch of size {len(batch)}")
        for aln in batch:
            # non-primary → drop & count
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                dropped_flag += 1
                continue

            # primary from here on
            seen_primary += 1
            if seen_primary % DEBUG_EVERY == 0:
                logger.debug(
                    f"Progress: primary={seen_primary}, kept={kept}, "
                    f"dropped_nonprimary={dropped_flag}, dropped_primary_other={dropped_primary_other}, "
                    f"dropped_short={dropped_short}",
                )
        
            category = ReadCategory.classify(aln.query_name, tag_config)

            if category is ReadCategory.OTHER and drop_untagged:
                dropped_primary_other += 1
                logger.debug(
                    f"Dropping untagged primary read '{aln.query_name}' due to --drop-untagged.",
                )
                continue

            # Determine trim extents
            if category is ReadCategory.OTHER:
                # Apply single-end policy to untagged reads that we keep
                left, right = max(0, policy.single_left), max(0, policy.single_right)
            else:
                left, right = category.trim_extents(policy)

            # Fast-path: no trimming → check length & write/drop
            if left == 0 and right == 0:
                qlen = aln.query_length or (
                    len(aln.query_sequence) if aln.query_sequence else 0
                )
                if qlen < policy.min_len:
                    dropped_short += 1
                    logger.debug(
                        f"Dropping short (no-trim) read '{aln.query_name}': length={qlen} < min_len={policy.min_len}",
                    )
                    continue
                outp.write(aln)
                kept += 1
                continue

            # Early-drop heuristic: predicted post-trim length (mode-dependent)
            qlen0 = aln.query_length or (
                len(aln.query_sequence) if aln.query_sequence else 0
            )
            # Mode-specific length prediction
            match clipping_mode:
                case ClippingMode.DELETE | ClippingMode.HARD_CLIP:
                    predicted_len = max(0, qlen0 - left - right)  # Length will decrease
                case ClippingMode.SOFT_CLIP:
                    predicted_len = qlen0  # Length stays the same

            if predicted_len < policy.min_len:
                dropped_short += 1
                if clipping_mode == ClippingMode.SOFT_CLIP:
                    logger.debug(
                        f"Dropping short (SOFT_CLIP) read '{aln.query_name}': "
                        f"original_len={predicted_len} < min_len={policy.min_len}",
                    )
                else:
                    logger.debug(
                        f"Dropping short (predicted post-trim) read '{aln.query_name}': "
                        f"{qlen0} - {left} - {right} = {predicted_len} < min_len={policy.min_len}",
                    )
                continue

            # Perform trimming (sequence/quals/CIGAR/ref_start updated consistently)
            if category is ReadCategory.OTHER:
                logger.debug(
                    f"Trimming untagged/single-end read '{aln.query_name}': left={left} right={right}",
                )
            else:
                logger.debug(
                    f"Trimming {category.name} read '{aln.query_name}': left={left} right={right}",
                )

            pre_trim_seq_len = len(aln.query_sequence) if aln.query_sequence else 0
            trim_alignment_in_place(aln, left, right, clipping_mode)

            # Strip tags if requested
            if strip_tags:
                aln.query_name = tag_config.strip_tags(aln.query_name)

            # Post-trim validation and final length check
            new_seq = aln.query_sequence or ""
            new_qual = aln.query_qualities
            new_len = len(new_seq)

            if new_qual is not None:
                assert len(new_qual) == new_len, (
                    f"Post-trim sequence/quality mismatch for '{aln.query_name}': seq={new_len}, qual={len(new_qual)}"
                )
            assert new_len <= pre_trim_seq_len, (
                f"Trimming increased sequence length for '{aln.query_name}': {pre_trim_seq_len} -> {new_len}"
            )

            if new_len < policy.min_len:
                dropped_short += 1
                logger.debug(
                    f"Dropping short (post-trim) read: length={new_len} < min_len={policy.min_len}",
                )
                continue

            outp.write(aln)
            kept += 1

    # Final invariants: counters must be consistent
    total_processed = kept + dropped_flag + dropped_short + dropped_primary_other
    assert (
        kept >= 0
        and dropped_flag >= 0
        and dropped_short >= 0
        and dropped_primary_other >= 0
    ), (  # noqa: PT018
        f"Counters cannot be negative: kept={kept}, dropped_nonprimary={dropped_flag}, "
        f"dropped_primary_other={dropped_primary_other}, dropped_short={dropped_short}"
    )
    assert seen_primary == kept + dropped_short + dropped_primary_other, (
        f"Primary read count inconsistency: seen={seen_primary}, kept={kept}, "
        f"dropped_short={dropped_short}, dropped_primary_other={dropped_primary_other}"
    )

    logger.info(
        f"Process totals — processed={total_processed}, kept={kept}, "
        f"dropped_nonprimary={dropped_flag}, dropped_primary_other={dropped_primary_other}, "
        f"dropped_short={dropped_short}",
    )

    # Preserve original return API: (kept, dropped_unmapped_or_nonprimary, dropped_too_short)
    return kept, dropped_flag, dropped_short


# --------------------------------- CLI ------------------------------------- #


def build_parser() -> argparse.ArgumentParser:
    """
    CLI:
      -v / -vv / -vvv : increase verbosity (INFO -> DEBUG -> TRACE)
      -q / -qq / -qqq : decrease verbosity (INFO -> WARNING -> ERROR -> CRITICAL)
    (Mutually exclusive.)
    """
    p = argparse.ArgumentParser(
        description=(
            "Trim aligned reads in SAM/BAM/CRAM based on read-name labels:\n"
            "  - merged:   <merged_prefix><name>\n"
            "  - unmerged: <unmerged_prefix><name>/1 or <unmerged_prefix><name>/2\n"
            "Updates CIGAR and reference_start consistently; removes unmapped/secondary/supplementary.\n"
            "Supports three clipping modes: delete (remove bases), soft-clip (mark in CIGAR), hard-clip (mark + remove)."
        ),
    )

    # I/O
    p.add_argument(
        "-i",
        "--in",
        dest="in_path",
        required=True,
        help="Input SAM/BAM/CRAM",
    )
    p.add_argument(
        "-o",
        "--out",
        dest="out_path",
        required=True,
        help="Output SAM/BAM/CRAM",
    )
    p.add_argument(
        "--ref",
        dest="reference",
        default=None,
        help="Reference FASTA (required/recommended for CRAM read/write)",
    )

    # Trimming policy
    p.add_argument(
        "--merged-left",
        type=int,
        default=30,
        help="Left trim for merged reads (query bases)",
    )
    p.add_argument(
        "--merged-right",
        type=int,
        default=30,
        help="Right trim for merged reads (query bases)",
    )
    p.add_argument(
        "--r1-left",
        type=int,
        default=30,
        help="Left trim for unmerged R1 reads",
    )
    p.add_argument(
        "--r2-right",
        type=int,
        default=30,
        help="Right trim for unmerged R2 reads",
    )
    p.add_argument(
        "--single-left",
        type=int,
        default=30,
        help="Left trim for single-end/untagged reads",
    )
    p.add_argument(
        "--single-right",
        type=int,
        default=30,
        help="Right trim for single-end/untagged reads",
    )
    p.add_argument(
        "--min-len",
        type=int,
        default=20,
        help="Minimum read length after trimming",
    )

    # Streaming
    p.add_argument(
        "--batch-size",
        type=int,
        default=10000,
        help="Batch size for streaming",
    )

    # Tag configuration
    tag_group = p.add_argument_group("Tag Configuration")
    tag_group.add_argument(
        "--merged-prefix",
        default="MERGED_",
        help="Prefix for merged reads (default: MERGED_)",
    )
    tag_group.add_argument(
        "--unmerged-prefix",
        default="UNMERGED_",
        help="Prefix for unmerged reads (default: UNMERGED_)",
    )
    tag_group.add_argument(
        "--strip-tags",
        action="store_true",
        help="Remove prefixes from read names in output",
    )

    # Clipping configuration
    clipping_group = p.add_argument_group("Clipping Configuration")
    clipping_group.add_argument(
        "--clipping-mode",
        choices=["delete", "soft-clip", "hard-clip"],
        default="delete",
        help=(
            "How to handle trimmed bases:\n"
            "  delete: Remove completely (default, best for most tools)\n"
            "  soft-clip: Mark in CIGAR as soft clips (S), preserve sequence\n"
            "  hard-clip: Mark in CIGAR as hard clips (H), remove sequence"
        ),
    )

    # Selection
    p.add_argument(
        "--drop-untagged",
        action="store_true",
        help="If set, drop reads that are not labeled with prefixes (otherwise single-end trimming is applied)",
    )

    # Verbosity: -v/-vv/-vvv or -q/-qq/-qqq (mutually exclusive)
    g = p.add_mutually_exclusive_group()
    g.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity (use up to -vvv).",
    )
    g.add_argument(
        "-q",
        "--quiet",
        action="count",
        default=0,
        help="Decrease verbosity (use up to -qqq).",
    )

    return p


def main(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    configure_logging(args.verbose, args.quiet)
    logger.info("Starting trimming run.")

    # Convert CLI string to ClippingMode enum
    clipping_mode_map = {
        "delete": ClippingMode.DELETE,
        "soft-clip": ClippingMode.SOFT_CLIP,
        "hard-clip": ClippingMode.HARD_CLIP,
    }
    try:
        clipping_mode = clipping_mode_map[args.clipping_mode]
    except KeyError:
        logger.error(f"Invalid clipping mode: {args.clipping_mode}")
        logger.error(f"Valid options: {list(clipping_mode_map.keys())}")
        sys.exit(1)

    policy = TrimPolicy(
        merged_left=max(0, args.merged_left),
        merged_right=max(0, args.merged_right),
        r1_left=max(0, args.r1_left),
        r2_right=max(0, args.r2_right),
        single_left=max(0, args.single_left),
        single_right=max(0, args.single_right),
        min_len=max(0, args.min_len),
        clipping_mode=clipping_mode,
    )
    logger.debug(f"TrimPolicy: {policy}")

    tag_config = TagConfig(
        merged_prefix=args.merged_prefix,
        unmerged_prefix=args.unmerged_prefix,
        # read1_suffix and read2_suffix keep defaults
    )
    logger.debug(f"TagConfig: {tag_config}")

    input_alignment = open_alignment(
        args.in_path,
        write=False,
        reference=args.reference,
    )
    try:
        output_alignment = open_alignment(
            args.out_path,
            write=True,
            template_or_header=input_alignment,
            reference=args.reference,
        )
    except (OSError, ValueError):
        input_alignment.close()
        raise

    try:
        kept, dropped_flag, dropped_short = process_stream(
            inp=input_alignment,
            outp=output_alignment,
            policy=policy,
            tag_config=tag_config,
            batch_size=max(1, args.batch_size),
            drop_untagged=bool(args.drop_untagged),
            strip_tags=bool(args.strip_tags),
        )
    finally:
        output_alignment.close()
        input_alignment.close()

    logger.success(
        f"Kept: {kept} | Dropped (unmapped/secondary/supplementary): {dropped_flag} | "
        f"Dropped (too short after trim): {dropped_short}",
    )
    logger.info("Trimming run complete.")


if __name__ == "__main__":
    main()
