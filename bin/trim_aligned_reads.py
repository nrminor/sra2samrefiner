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

MERGED_PREFIX: str = "MERGED_"
UNMERGED_PREFIX: str = "UNMERGED_"
READ1_SUFFIX: str = "/1"
READ2_SUFFIX: str = "/2"

# CIGAR op codes
# 0:M, 1:I, 2:D, 3:N, 4:S, 5:H, 6:P, 7:=, 8:X
REF_CONSUME = {0, 2, 3, 7, 8}
QRY_CONSUME = {0, 1, 4, 7, 8}
BOTH_CONSUME = {0, 7, 8}

# Emit a progress debug line after processing this many primary reads
DEBUG_EVERY: int = 100_000


# ------------------------------- DATA TYPES -------------------------------- #


@dataclass(frozen=True)
class TrimPolicy:
    """Trim amounts in query-bases for each category."""

    merged_left: int = 30
    merged_right: int = 30
    r1_left: int = 30
    r2_right: int = 30
    min_len: int = 20  # drop reads shorter than this AFTER trimming


class ReadCategory(Enum):
    """Closed set of categories derived from read names, with trim logic."""

    MERGED = auto()
    UNMERGED_R1 = auto()
    UNMERGED_R2 = auto()
    OTHER = auto()

    @staticmethod
    def classify(qname: str) -> ReadCategory:
        """Classify by upstream naming scheme."""
        if qname.startswith(MERGED_PREFIX):
            return ReadCategory.MERGED
        if qname.startswith(UNMERGED_PREFIX) and qname.endswith(READ1_SUFFIX):
            return ReadCategory.UNMERGED_R1
        if qname.startswith(UNMERGED_PREFIX) and qname.endswith(READ2_SUFFIX):
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
                return 0, max(0, policy.r2_right)
            case ReadCategory.OTHER:
                return 0, 0


# ----------------------------- LOGGING SETUP ------------------------------- #


def configure_logging(verbose: int, quiet: int) -> None:
    """
    Base at INFO (0). Positive → louder, negative → quieter.
    Map: +1=DEBUG, +2..=TRACE | -1=WARNING, -2=ERROR, -3..=CRITICAL
    """
    logger.remove()
    level_str = "INFO"
    match verbose - quiet:
        case delta if delta >= 2:  # noqa: PLR2004
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
        case delta if delta < -2:  # noqa: PLR2004
            level_str = "CRITICAL"
        case _:
            level_str = "INFO"
            logger.add(sys.stderr, level=level_str)
            logger.warning("Invalid verbosity settings provided; defaulting to INFO.")
            return
    logger.add(sys.stderr, level=level_str)
    logger.debug(f"Logger configured at level: {level_str}")


# ---------------------------- CIGAR UTILITIES ------------------------------ #


class CigarOp(NamedTuple):
    """One CIGAR run: (operation code, run length)."""

    op: int
    length: int

    @staticmethod
    def from_tuple(t: tuple[int, int]) -> "CigarOp":
        """Convert a raw (op, len) tuple to CigarOp."""
        op, ln = t
        return CigarOp(op, ln)

    @staticmethod
    def to_tuple(run: "CigarOp") -> tuple[int, int]:
        """Convert a CigarOp back to a raw (op, len) tuple."""
        return (run.op, run.length)


class Cigar(list[CigarOp]):
    """A list of CigarOp with helpers for conversion and compaction."""

    @classmethod
    def from_pysam(cls, cig_raw: list[tuple[int, int]] | None) -> "Cigar | None":
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
    if trim_q <= 0 or not cig:
        return cig, 0

    out = Cigar()
    remaining = trim_q
    ref_advance = 0
    i = 0

    while i < len(cig) and remaining > 0:
        run = cig[i]
        if run.op in QRY_CONSUME:
            take = min(run.length, remaining)
            keep_len = run.length - take
            remaining -= take
            if run.op in BOTH_CONSUME:
                ref_advance += take
            out.push_compact(run.op, keep_len)
        else:
            # D/N/H/P: keep fully; they don't consume query so `remaining` stays.
            out.push_compact(run.op, run.length)
        i += 1

    # Append untouched tail
    for j in range(i, len(cig)):
        tail = cig[j]
        out.push_compact(tail.op, tail.length)

    return out, ref_advance


def _consume_from_right(cig: Cigar, trim_q: int) -> Cigar:
    """
    Consume `trim_q` query-bases from the RIGHT edge of `cig`.
    Reference start is unaffected by definition.
    """
    if trim_q <= 0 or not cig:
        return cig

    out_rev = Cigar()
    remaining = trim_q

    # Build reversed output with compaction
    for run in reversed(cig):
        if remaining > 0 and run.op in QRY_CONSUME:
            take = min(run.length, remaining)
            keep_len = run.length - take
            remaining -= take
            out_rev.push_compact(run.op, keep_len)
        else:
            out_rev.push_compact(run.op, run.length)

    # Reverse back (preserving compaction)
    out = Cigar()
    for run in reversed(out_rev):
        out.push_compact(run.op, run.length)
    return out


def trim_alignment_in_place(aln: pysam.AlignedSegment, left: int, right: int) -> None:
    """
    Trim query-bases from the read in-place and update sequence, qualities, CIGAR,
    and reference_start consistently. Strand-aware:
      - For reverse reads, the 5' end is on the RIGHT of the CIGAR.
      - Sequence slicing is always in read orientation.
    """
    left = max(left, 0)
    right = max(right, 0)

    seq: str | None = aln.query_sequence
    if seq is None:
        logger.debug(
            f"Skip trim: the read '{aln.query_name}' has no query_sequence (likely unmapped)."
        )
        return

    qual = aln.query_qualities  # array('B') or None
    cig_raw = aln.cigartuples
    if cig_raw is None:
        new_seq = seq[left : len(seq) - right if right else None]
        new_qual = None if qual is None else qual[left : len(seq) - right if right else None]
        aln.query_sequence = new_seq
        aln.query_qualities = new_qual
        logger.debug(f"Trimmed sequence/qualities only (no CIGAR) for '{aln.query_name}'.")
        return

    cig: Cigar = Cigar.from_pysam(cig_raw)
    ref_start_before = aln.reference_start

    # Strand-aware CIGAR trimming amounts
    cig_left = right if aln.is_reverse else left
    cig_right = left if aln.is_reverse else right
    if cig_left or cig_right:
        logger.debug(
            f"CIGAR trim for '{aln.query_name}': is_rev={aln.is_reverse}, left={cig_left}, right={cig_right}, "
            f"ref_start_before={ref_start_before}",
        )

    if cig_left > 0:
        cig, ref_advance = _consume_from_left(cig, cig_left)
        aln.reference_start = ref_start_before + ref_advance
    if cig_right > 0:
        cig = _consume_from_right(cig, cig_right)

    # Now trim sequence/qualities in read orientation
    qlen = len(seq)
    cut_left = min(left, qlen)
    cut_right = min(right, max(qlen - cut_left, 0))
    keep_end = qlen - cut_left - cut_right
    new_seq = "" if keep_end <= 0 else seq[cut_left : cut_left + keep_end]
    new_qual = None if qual is None else qual[cut_left : cut_left + keep_end]

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
    mode = _io_mode_from_ext(path, write)
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
        if isinstance(template_or_header, pysam.AlignmentFile):
            return pysam.AlignmentFile(path, mode, template=template_or_header, **kwargs)
        if isinstance(template_or_header, dict):
            return pysam.AlignmentFile(path, mode, header=template_or_header, **kwargs)
        msg = "Writing requires either a template AlignmentFile or a header dict."
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


def process_stream(
    inp: pysam.AlignmentFile,
    outp: pysam.AlignmentFile,
    policy: TrimPolicy,
    batch_size: int = 10000,
    drop_untagged: bool = False,  # noqa: FBT001, FBT002
) -> tuple[int, int, int]:
    """
    Stream input -> output in batches, trimming per read-name class.
    - Drops unmapped, secondary, and supplementary alignments.
    - Optionally drops untagged reads (if you only want labeled ones).
    Returns counts: (kept, dropped_unmapped_or_nonprimary, dropped_too_short)
    """
    kept = 0
    dropped_flag = 0
    dropped_short = 0
    seen_primary = 0

    for batch in batched(inp, batch_size):
        logger.debug(f"Processing batch of size {len(batch)}")
        for aln in batch:
            # primary-only
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                dropped_flag += 1
                continue

            seen_primary += 1
            if seen_primary % DEBUG_EVERY == 0:
                logger.debug(
                    f"Progress: primary={seen_primary}, kept={kept}, "
                    f"dropped_flag={dropped_flag}, dropped_short={dropped_short}",
                )

            category = ReadCategory.classify(aln.query_name)
            if category is ReadCategory.OTHER and drop_untagged:
                dropped_flag += 1
                logger.debug(f"Dropping untagged read '{aln.query_name}' due to --drop-untagged.")
                continue

            left, right = category.trim_extents(policy)
            if left == 0 and right == 0:
                qlen = aln.query_length or (len(aln.query_sequence) if aln.query_sequence else 0)
                if qlen < policy.min_len:
                    dropped_short += 1
                    logger.debug(
                        f"Dropping short (no-trim) read '{aln.query_name}': length={qlen} < min_len={policy.min_len}",
                    )
                    continue
                outp.write(aln)
                kept += 1
                continue

            logger.debug(f"Trimming {category.name}: left={left} right={right}")
            trim_alignment_in_place(aln, left, right)

            new_len = len(aln.query_sequence or "")
            if new_len < policy.min_len:
                dropped_short += 1
                logger.debug(
                    f"Dropping short (post-trim) read: length={new_len} < min_len={policy.min_len}",
                )
                continue

            outp.write(aln)
            kept += 1

    logger.info(
        f"Process summary — kept={kept}, dropped_flag={dropped_flag}, dropped_short={dropped_short}",
    )
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
            f"  - merged:   {MERGED_PREFIX}<name>\n"
            f"  - unmerged: {UNMERGED_PREFIX}<name>{READ1_SUFFIX} or {UNMERGED_PREFIX}<name>{READ2_SUFFIX}\n"
            "Updates CIGAR and reference_start consistently; removes unmapped/secondary/supplementary."
        ),
    )

    # I/O
    p.add_argument("-i", "--in", dest="in_path", required=True, help="Input SAM/BAM/CRAM")
    p.add_argument("-o", "--out", dest="out_path", required=True, help="Output SAM/BAM/CRAM")
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
    p.add_argument("--r1-left", type=int, default=30, help="Left trim for unmerged R1 reads")
    p.add_argument("--r2-right", type=int, default=30, help="Right trim for unmerged R2 reads")
    p.add_argument("--min-len", type=int, default=20, help="Minimum read length after trimming")

    # Streaming
    p.add_argument("--batch-size", type=int, default=10000, help="Batch size for streaming")

    # Selection
    p.add_argument(
        "--drop-untagged",
        action="store_true",
        help="If set, drop reads that are not labeled MERGED_/UNMERGED_",
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

    policy = TrimPolicy(
        merged_left=max(0, args.merged_left),
        merged_right=max(0, args.merged_right),
        r1_left=max(0, args.r1_left),
        r2_right=max(0, args.r2_right),
        min_len=max(0, args.min_len),
    )
    logger.debug(f"TrimPolicy: {policy}")

    input_alignment = open_alignment(args.in_path, write=False, reference=args.reference)
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
            batch_size=max(1, args.batch_size),
            drop_untagged=bool(args.drop_untagged),
        )
    finally:
        output_alignment.close()
        input_alignment.close()

    logger.success(
        f"Kept: {kept} | Dropped (unmapped/secondary/supplementary/untagged): {dropped_flag} | "
        f"Dropped (too short after trim): {dropped_short}",
    )
    logger.info("Trimming run complete.")


if __name__ == "__main__":
    main()
