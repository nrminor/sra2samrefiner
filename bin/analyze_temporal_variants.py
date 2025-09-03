#!/usr/bin/env python3
"""
Temporal variant analysis for SRATracking functionality.

Replaces the functionality of Comp_time_windows.py with proper batch processing
that faithfully replicates the original domain logic.
"""

import argparse
import gzip
import re
import sys
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Literal, NewType, cast, overload

import polars as pl
from loguru import logger
from pydantic import Field, ValidationInfo, field_validator
from pydantic.dataclasses import dataclass

# Constants to eliminate magic numbers
ORF1AB_SPLIT_POSITION = 4401  # Position where ORF1ab splits into ORF1a/ORF1b
MIN_COVARS_COLUMNS = 3  # SAMRefiner: variant, count, abundance
MIN_NTCALLS_COLUMNS = 11  # SAMRefiner: need at least through "Counts" column
NTCALLS_COUNTS_COLUMN = 10  # Column index for "Counts" in SAMRefiner nt_calls (0-indexed)
THRESHOLD_100_READS = 100  # Threshold for samples_100_reads statistic
THRESHOLD_1K_READS = 1000  # Threshold for samples_1k_reads statistic
MIN_TSV_COLUMNS_FOR_LINEAGE = 3  # Minimum columns in lineage TSV files
MIN_MUTATION_PARTS = 2  # Minimum parts when splitting mutation on ":"
EXPECTED_WINDOW_PAIR_SIZE = 2  # Expected size for time window pairs

# Define statically analyzable SARS-CoV-2 ORF types
Orf = Literal[
    "S",
    "N",
    "ORF1a",
    "ORF1b",
    "ORF3a",
    "ORF6",
    "ORF7a",
    "ORF7b",
    "ORF8",
    "ORF10",
    "Del",
]

# Semantic domain types for static analysis
OrfName = NewType("OrfName", Orf)
AAChange = NewType("AAChange", str)
LineageName = NewType("LineageName", str)


@dataclass(frozen=True)
class MutationKey:
    """Immutable key for mutation lookups with type safety."""

    orf: OrfName
    aa_change: AAChange

    def __hash__(self) -> int:
        return hash((self.orf, self.aa_change))


@dataclass
class LineageMapping:
    """Represents lineages associated with a specific mutation."""

    mutation_key: MutationKey
    lineages: list[LineageName]
    cdc_ranking: int = Field(default=-1)  # -1 if not in CDC data


class MutationLookup:
    """
    Type-safe mutation lookup with multiple access patterns.

    Supports both dict-style and method-based access with static type checking.
    Replaces the weakly-typed nested Dict[str, Dict[str, List[str]]] structure.
    """

    def __init__(self) -> None:
        self._data: dict[MutationKey, LineageMapping] = {}

    # Dict-style access with MutationKey
    def __getitem__(self, key: MutationKey) -> LineageMapping:
        return self._data[key]

    def __contains__(self, key: MutationKey) -> bool:
        return key in self._data

    def __setitem__(self, key: MutationKey, value: LineageMapping) -> None:
        self._data[key] = value

    # Method-based access with automatic type conversion
    @overload
    def get(self, orf: OrfName, aa_change: AAChange) -> LineageMapping | None: ...

    @overload
    def get(self, orf: Orf, aa_change: str) -> LineageMapping | None: ...

    def get(self, orf, aa_change):
        """Get lineage mapping with automatic type casting and validation."""
        # Cast to proper NewTypes with static analysis support
        orf_typed = cast("OrfName", orf)
        aa_change_typed = cast("AAChange", aa_change)

        key = MutationKey(orf_typed, aa_change_typed)
        return self._data.get(key)

    def add_mapping(
        self,
        orf: Orf,
        aa_change: str,
        lineages: list[str],
        cdc_ranking: int = -1,
    ) -> None:
        """Add a new mapping with compile-time ORF validation."""
        key = MutationKey(
            cast("OrfName", orf),  # Type checker enforces valid ORF literals
            AAChange(aa_change),
        )
        lineage_names = [LineageName(lineage) for lineage in lineages]
        self._data[key] = LineageMapping(key, lineage_names, cdc_ranking)

    def get_all_orfs(self) -> list[OrfName]:
        """Get all ORFs present in the lookup."""
        return list({key.orf for key in self._data})

    def get_lineages_for_orf(self, orf: Orf) -> list[LineageMapping]:
        """Get all mutations and their lineages for a specific ORF."""
        orf_typed = cast("OrfName", orf)
        return [mapping for key, mapping in self._data.items() if key.orf == orf_typed]

    def size(self) -> int:
        """Get the total number of mutations in the lookup."""
        return len(self._data)


@dataclass
class AnalysisConfig:
    """Configuration for temporal variant analysis."""

    delta_threshold: float = Field(default=0.02, ge=0.0, le=1.0)
    count_threshold: int = Field(default=100, ge=1)
    window_weeks: int = Field(default=3, ge=1, le=52)
    min_samples_per_window: int = Field(default=20, ge=1)


@dataclass
class TimeWindow:
    """Represents a time window for variant comparison."""

    name: str = Field(min_length=1)
    start_date: date = Field()
    end_date: date = Field()

    @field_validator("end_date")
    @classmethod
    def end_after_start(cls, v: date, info: ValidationInfo) -> date:
        if info.data and "start_date" in info.data and v <= info.data["start_date"]:
            error_msg = "end_date must be after start_date"
            raise ValueError(error_msg)
        return v


@dataclass
class VariantStatistics:
    """Four-component statistics tracking for variants (matching original)."""

    read_count: int = Field(ge=0)  # Total read count
    sample_count: int = Field(ge=0)  # Number of samples
    samples_100_reads: int = Field(ge=0)  # Samples with >=100 reads
    samples_1k_reads: int = Field(ge=0)  # Samples with >=1000 reads


@dataclass
class TemporalComparisonData:
    """Bundle of data needed for temporal comparison to reduce function parameters."""

    class Config:
        arbitrary_types_allowed = True  # Allow polars DataFrames

    covars_data: pl.DataFrame
    ntcalls_data: pl.DataFrame
    sample_metadata: pl.DataFrame
    baseline_variants: list[str]
    ranked_variants: list[str]
    mutation_lookup: MutationLookup


@dataclass
class AnalysisInputs:
    """Bundle of input paths and configuration for temporal analysis."""

    input_dir: Path
    sample_metadata_path: Path
    lineage_mutations_path: Path
    variant_surveillance_path: Path
    output_paths: dict[str, Path]


def extract_position_from_variant(variant_string: str) -> int:
    """
    Extract genomic position from variant string.

    Replicates get_pos() function from original code.
    """
    match = re.search(r"\d+", variant_string)
    return int(match.group(0)) if match else -1


def normalize_orf_name(orf: str, aa_change: str) -> Orf:
    """
    Normalize ORF names following original logic with type safety.

    Handles ORF1ab splitting based on AA position (Borf function logic).
    """
    if "ORF1ab" in orf:
        # Extract position from AA change for ORF1ab splitting
        aa_position = extract_position_from_variant(aa_change)
        return cast("Orf", "ORF1b" if aa_position > ORF1AB_SPLIT_POSITION else "ORF1a")

    if "ORF" in orf:
        normalized = orf.split("_")[0]
        # Validate against known ORFs
        valid_orfs = [
            "S",
            "N",
            "ORF1a",
            "ORF1b",
            "ORF3a",
            "ORF6",
            "ORF7a",
            "ORF7b",
            "ORF8",
            "ORF10",
            "Del",
        ]
        return cast("Orf", normalized if normalized in valid_orfs else "S")  # Default fallback

    # Handle single character ORFs (e.g., "S" -> "S")
    single_char = orf[0].upper() if orf else "S"
    return cast("Orf", single_char if single_char in ["S", "N"] else "S")  # Safe fallback


def load_baseline_jn1_variants(
    lineage_mutations_path: Path,
    target_lineage: str = "JN.1",
) -> dict[str, list[str]]:
    """
    Load baseline JN.1 variants to exclude from analysis.

    Replicates the JN_changes and JN_AAchanges logic from original code.
    """
    baseline_variants = []
    baseline_aa_changes = {}

    if not lineage_mutations_path.exists():
        logger.warning("No lineage mutations file found for baseline exclusion")
        return {"variants": [], "aa_changes": {}}

    # Find JN.1 baseline mutations with lazy loading
    jn1_mutations = (
        pl.scan_csv(lineage_mutations_path, separator="\t")  # Lazy loading
        .filter(pl.col("lineage") == target_lineage)
        .select("mutations")
        .collect()  # Materialize only the filtered data
        .get_column("mutations")
        .to_list()
    )

    if not jn1_mutations:
        logger.warning(f"No baseline mutations found for {target_lineage}")
        return {"variants": [], "aa_changes": {}}

    # Parse JN.1 mutations
    for mutation_set in jn1_mutations:
        for mutation in mutation_set.split(";"):
            if "|" not in mutation:
                continue

            baseline_variants.append(mutation)
            # Parse AA changes by ORF
            for pm_change_raw in mutation.split("|")[1:]:
                pm_change = pm_change_raw.strip("()")
                if ":" not in pm_change:
                    continue

                orf = pm_change.split(":")[0]
                aa_change = pm_change.split("(")[-1].strip(")")

                if orf not in baseline_aa_changes:
                    baseline_aa_changes[orf] = []
                baseline_aa_changes[orf].append(aa_change)

    logger.info(f"Loaded {len(baseline_variants)} baseline variants for exclusion")
    return {"variants": baseline_variants, "aa_changes": baseline_aa_changes}


def load_cdc_variant_ranking(variant_surveillance_path: Path) -> list[str]:
    """
    Load and rank variants by CDC circulation data.

    Replicates the ranked_vars logic from original code.
    """
    if not variant_surveillance_path.exists():
        logger.warning("No CDC surveillance data found for variant ranking")
        return []

    # Load and process CDC data with lazy evaluation
    return (
        pl.scan_csv(variant_surveillance_path, separator="\t")  # Lazy loading
        .filter(pl.col("usa_or_hhsregion") == "USA")
        .filter(pl.col("modeltype") == "smoothed")
        .filter(pl.col("share") > 0)
        .group_by("variant")
        .agg(pl.col("share").mean().alias("avg_share"))
        .with_columns(pl.col("variant") + "*")  # Add asterisk matching original
        .sort("avg_share", descending=True)
        .collect()  # Materialize only at the end
        .get_column("variant")
        .to_list()
    )


def parse_single_covars_file(covars_file: Path) -> list[dict[str, str | int]]:
    """
    Parse a single SAMRefiner covars file and convert to 4-component statistics.

    SAMRefiner format: Co-Variants	Count	Abundance
    Converts to SRATracking 4-component format matching lines 280-285.
    """
    accession = covars_file.stem.replace(".SARS2.wg", "")  # Extract clean accession
    covars_data = []

    # Handle both gzipped and plain text files
    file_opener = gzip.open if str(covars_file).endswith(".gz") else open

    with file_opener(covars_file, "rt") as f:
        # Skip first 2 header lines (SAMRefiner format)
        f.readline()  # Sample header: SRR32418827.SARS2.wg(476428)
        f.readline()  # Column headers: Co-Variants	Count	Abundance

        for line in f:
            if "Reference" in line:
                continue

            parts = line.strip().split("\t")
            if len(parts) < MIN_COVARS_COLUMNS:  # Need variant, count, abundance
                continue

            try:
                variant = parts[0]  # Full variant string with AA changes
                raw_count = int(parts[1])  # SAMRefiner count
                # Ignore abundance (parts[2]) - we calculate from counts

                # Convert to 4-component statistics (replicating original lines 280-285)
                four_component_stats = convert_to_4component_stats(raw_count)

                covars_data.append(
                    {
                        "accession": accession,
                        "variant": variant,
                        **four_component_stats,
                    },
                )

            except (ValueError, IndexError):
                continue

    return covars_data


def convert_to_4component_stats(raw_count: int) -> dict[str, int]:
    """
    Convert SAMRefiner count to 4-component statistics.

    Replicates the logic from Comp_time_windows.py lines 280-285:
    cur_stats = [count, 1, 0, 0]  # [read_count, sample_count, 100+, 1k+]
    """
    stats = {
        "read_count": raw_count,
        "sample_count": 1,  # Each file represents 1 sample
        "samples_100_reads": 0,
        "samples_1k_reads": 0,
    }

    # Apply thresholds matching original logic
    if raw_count >= THRESHOLD_100_READS:
        stats["samples_100_reads"] = 1
        if raw_count >= THRESHOLD_1K_READS:
            stats["samples_1k_reads"] = 1

    return stats


def parse_single_ntcalls_file(ntcalls_file: Path) -> list[dict[str, str | int]]:
    """
    Parse a single SAMRefiner nt_calls file and extract coverage data.

    SAMRefiner format: Position	ref NT	AAs	A	T	C	G	-	N	Total	Primary NT	Counts	Abundance...
    Extracts position and coverage matching SRATracking line 218.
    """
    accession = ntcalls_file.stem.replace(".SARS2.wg", "")  # Extract clean accession
    ntcalls_data = []

    # Handle both gzipped and plain text files
    file_opener = gzip.open if str(ntcalls_file).endswith(".gz") else open

    with file_opener(ntcalls_file, "rt") as f:
        # Skip first 2 header lines (SAMRefiner format)
        f.readline()  # Sample header: SRR32418827.SARS2.wg(476428)
        f.readline()  # Column headers: Position	ref NT	AAs	A	T	C	G	-	N	Total	Primary NT	Counts...

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < MIN_NTCALLS_COLUMNS:  # Need at least through "Counts" column
                continue

            # Check column 2 truthy (AAs column) - matching original line 217
            if not parts[2]:
                continue

            try:
                position = int(parts[0])  # Position column
                coverage = int(parts[NTCALLS_COUNTS_COLUMN])  # "Counts" column (Primary NT Counts)
            except (ValueError, IndexError):
                continue

            ntcalls_data.append(
                {
                    "accession": accession,
                    "position": position,
                    "coverage": coverage,
                },
            )

    return ntcalls_data


def load_samrefiner_outputs(input_dir: Path) -> dict[str, pl.DataFrame]:
    """
    Load and parse all SAMRefiner output files with proper tidyverse-style polars.

    Supports both daily and retrospective modes by loading all files in directory.
    Handles mixed current + historical data for temporal analysis.
    """
    logger.info(f"Loading SAMRefiner outputs from {input_dir}")

    # Match actual SAMRefiner output file patterns (works for both modes)
    covars_files = list(input_dir.glob("*_covars.tsv*"))  # Both .tsv and .tsv.gz
    ntcalls_files = list(input_dir.glob("*_nt_calls.tsv*"))  # Both .tsv and .tsv.gz

    if not covars_files or not ntcalls_files:
        logger.warning(f"No SAMRefiner files found in {input_dir}")
        return {"covars": pl.DataFrame(), "ntcalls": pl.DataFrame()}

    logger.info(f"Found {len(covars_files)} covars files, {len(ntcalls_files)} ntcalls files")
    logger.info("Processing files from both current run and historical data (if present)")

    # Load all covars data using helper function
    all_covars_data = [
        record for covars_file in covars_files for record in parse_single_covars_file(covars_file)
    ]

    # Load all ntcalls data using helper function
    all_ntcalls_data = [
        record
        for ntcalls_file in ntcalls_files
        for record in parse_single_ntcalls_file(ntcalls_file)
    ]

    logger.info(f"Loaded {len(all_covars_data)} variant records and {len(all_ntcalls_data)} coverage records")

    return {
        "covars": pl.DataFrame(all_covars_data),
        "ntcalls": pl.DataFrame(all_ntcalls_data),
    }


def enrich_samrefiner_with_metadata(
    samrefiner_data: dict[str, pl.DataFrame],
    sample_metadata: pl.DataFrame,
) -> dict[str, pl.DataFrame]:
    """
    Enrich SAMRefiner data with collection dates from metadata.

    Critical for dual-mode operation where dates determine temporal grouping,
    not file timestamps or directory structure.
    """
    if samrefiner_data["covars"].height == 0:
        logger.warning("No SAMRefiner covars data to enrich")
        return samrefiner_data

    if sample_metadata.height == 0:
        logger.warning("No sample metadata available for date enrichment")
        return samrefiner_data

    # Enrich covars data with collection dates
    enriched_covars = (
        samrefiner_data["covars"]
        .join(
            sample_metadata.select(["Sample", "Collection Date"]),
            left_on="accession",
            right_on="Sample",
            how="left",
        )
        .filter(pl.col("Collection Date").is_not_null())  # Only keep samples with known dates
    )

    # Enrich ntcalls data with collection dates
    enriched_ntcalls = (
        samrefiner_data["ntcalls"]
        .join(
            sample_metadata.select(["Sample", "Collection Date"]),
            left_on="accession",
            right_on="Sample",
            how="left",
        )
        .filter(pl.col("Collection Date").is_not_null())  # Only keep samples with known dates
    )

    logger.info(f"Enriched data: {enriched_covars.height} covars, {enriched_ntcalls.height} ntcalls with dates")

    # Log any samples that couldn't be matched (important for debugging)
    unmatched_covars = samrefiner_data["covars"].height - enriched_covars.height
    if unmatched_covars > 0:
        logger.warning(f"{unmatched_covars} covars records lacked collection date metadata")

    return {
        "covars": enriched_covars,
        "ntcalls": enriched_ntcalls,
    }


def create_time_windows_exact(weekly_data: pl.DataFrame) -> list[list[TimeWindow]]:
    """
    Create exact sliding time window comparisons matching original triweeks logic.

    Replicates the precise date arithmetic from weekly_out() function.
    """
    weeks = weekly_data.get_column("collection_week").to_list()
    comparison_sets = []

    # Need at least 6 consecutive weeks (matching original range check)
    for j in range(len(weeks) - 5):
        # Replicate exact original triweek calculation
        anchor_week = weeks[j]
        anchor_week3 = weeks[j + 3]

        # Convert to datetime for arithmetic (matching original strptime/strftime)
        anchor_dt = (
            datetime.combine(anchor_week, datetime.min.time())
            if isinstance(anchor_week, date)
            else anchor_week
        )
        anchor3_dt = (
            datetime.combine(anchor_week3, datetime.min.time())
            if isinstance(anchor_week3, date)
            else anchor_week3
        )

        # Exact original formula replication
        triweek1_start = (anchor_dt - timedelta(days=14)).date()
        triweek1_end = (anchor_dt + timedelta(days=6)).date()
        triweek2_start = (anchor3_dt - timedelta(days=14)).date()
        triweek2_end = (anchor3_dt + timedelta(days=6)).date()

        window_pair = [
            TimeWindow(
                name=f"{triweek1_start.strftime('%Y-%m-%d')}--{triweek1_end.strftime('%Y-%m-%d')}",
                start_date=triweek1_start,
                end_date=triweek1_end,
            ),
            TimeWindow(
                name=f"{triweek2_start.strftime('%Y-%m-%d')}--{triweek2_end.strftime('%Y-%m-%d')}",
                start_date=triweek2_start,
                end_date=triweek2_end,
            ),
        ]

        comparison_sets.append(window_pair)

    logger.info(f"Created {len(comparison_sets)} exact comparison window sets")
    return comparison_sets


def calculate_window_variant_stats(
    covars_data: pl.DataFrame,
    ntcalls_data: pl.DataFrame,
    sample_metadata: pl.DataFrame,
    window: TimeWindow,
    baseline_variants: list[str],
) -> pl.DataFrame:
    """
    Calculate variant statistics for a single time window.

    Uses proper tidyverse-style polars chaining and replicates original 4-component stats.
    """
    return (
        covars_data.join(
            sample_metadata.select(["Sample", "Collection Date"]),
            left_on="accession",
            right_on="Sample",
        )
        .with_columns(
            pl.col("Collection Date")
            .str.to_date("%Y-%m-%d", strict=False)
            .alias("collection_date"),
        )
        .filter(pl.col("collection_date").is_not_null())
        .filter(
            (pl.col("collection_date") >= window.start_date)
            & (pl.col("collection_date") <= window.end_date),
        )
        .filter(~pl.col("variant").is_in(baseline_variants))  # Exclude JN.1 baseline variants
        .filter(pl.col("variant").str.contains(r"\|"))  # Only variants with "|" (matching original)
        .filter(~pl.col("variant").str.ends_with("N"))  # Exclude N mutations (matching original)
        .with_columns(
            pl.col("variant")
            .map_elements(lambda x: extract_position_from_variant(x), return_dtype=pl.Int32)
            .alias("position"),
        )
        .group_by(["variant", "position"])
        .agg(
            [
                pl.col("read_count").sum().alias("total_read_count"),
                pl.col("sample_count").sum().alias("total_sample_count"),
                pl.col("samples_100_reads").sum().alias("total_samples_100"),
                pl.col("samples_1k_reads").sum().alias("total_samples_1k"),
            ],
        )
        .join(
            ntcalls_data.join(
                sample_metadata.select(["Sample", "Collection Date"]),
                left_on="accession",
                right_on="Sample",
            )
            .with_columns(
                pl.col("Collection Date")
                .str.to_date("%Y-%m-%d", strict=False)
                .alias("collection_date"),
            )
            .filter(pl.col("collection_date").is_not_null())
            .filter(
                (pl.col("collection_date") >= window.start_date)
                & (pl.col("collection_date") <= window.end_date),
            )
            .group_by("position")
            .agg(pl.col("coverage").sum().alias("total_coverage")),
            on="position",
        )
        .with_columns(
            [
                (pl.col("total_read_count") / pl.col("total_coverage")).alias("abundance"),
                pl.lit(window.name).alias("time_window"),
            ],
        )
        .filter(pl.col("abundance") > 0)  # Only keep non-zero abundances
    )


def extract_variant_components(variant: str) -> dict[str, str]:
    """
    Extract all components from variant string with full original parsing logic.

    Replicates the complex ORF parsing from original weekly_out() function.
    """
    if "|" not in variant:
        return {
            "nt_change": variant,
            "orf": "",
            "aa_change": "",
            "position": str(extract_position_from_variant(variant)),
        }

    parts = variant.split("|")
    nt_change = parts[0]
    position = extract_position_from_variant(nt_change)

    # Extract AA changes (there can be multiple)
    orfs = []
    aa_changes = []

    for pm_change_raw in parts[1:]:
        pm_change = pm_change_raw.strip("()")
        aa_change = pm_change.split("(")[-1].strip(")") if "(" in pm_change else pm_change

        orf = pm_change.split(":")[0] if ":" in pm_change else ""

        # Apply ORF normalization (matching original logic)
        orf = normalize_orf_name(orf, aa_change)

        # Skip silent mutations (original logic: aa_change[0] == aa_change[-1])
        # Check for silent mutations (same first and last amino acid)
        min_aa_length_for_silent_check = 2
        if len(aa_change) >= min_aa_length_for_silent_check and aa_change[0] == aa_change[-1]:
            continue

        orfs.append(orf)
        aa_changes.append(aa_change)

    return {
        "nt_change": nt_change,
        "orf": ", ".join(list(set(orfs))) if orfs else "",
        "aa_change": ", ".join(aa_changes) if aa_changes else "",
        "position": str(position),
    }


def build_lineage_mutation_lookup(lineage_mutations: pl.DataFrame) -> MutationLookup:
    """
    Build type-safe lineage mutation lookup.

    Replicates JN_clade_PM_dict from original code with modern type system.
    """
    lookup = MutationLookup()

    if lineage_mutations.height == 0:
        logger.warning("No lineage mutations provided")
        return lookup

    # Process each lineage and its mutations
    for row in lineage_mutations.iter_rows(named=True):
        lineage = row["lineage"]
        if not lineage.endswith("*"):
            lineage += "*"  # Add asterisk matching original

        mutations = row["mutations"].split(";") if row["mutations"] else []

        for mutation in mutations:
            if ":" not in mutation:
                continue

            parts = mutation.split(":")
            if len(parts) < MIN_MUTATION_PARTS:
                continue

            orf_str = parts[0]
            aa_change = parts[1]

            # Validate ORF against known/supported SARS-CoV-2 ORFs
            valid_orfs = [
                "S",
                "N",
                "ORF1a",
                "ORF1b",
                "ORF3a",
                "ORF6",
                "ORF7a",
                "ORF7b",
                "ORF8",
                "ORF10",
                "Del",
            ]
            if orf_str not in valid_orfs:
                logger.debug(f"Skipping unknown ORF: {orf_str}")
                continue

            # Get existing mapping or create new one
            existing = lookup.get(orf_str, aa_change)
            if existing:
                # Add to existing lineage list if not already present
                if LineageName(lineage) not in existing.lineages:
                    existing.lineages.append(LineageName(lineage))
            else:
                # Create new mapping
                lookup.add_mapping(orf_str, aa_change, [lineage])

    logger.info(
        f"Built mutation lookup with {len(lookup.get_all_orfs())} ORFs, {lookup.size()} total mutations",
    )
    return lookup


def rank_associated_variants(
    variant_components: dict[str, str],
    ranked_variants: list[str],
    mutation_lookup: MutationLookup,
) -> str:
    """
    Rank associated variants by CDC circulation data using type-safe lookup.

    Replicates the variant ranking logic from original code with modern types.
    """
    orf = variant_components["orf"]
    aa_change = variant_components["aa_change"]
    nt_change = variant_components["nt_change"]

    # Handle deletion cases (matching original logic)
    if "del" in nt_change.lower():
        orf_key = "Del"
        aa_key = nt_change.strip("ATCGdel")
    else:
        orf_key = orf
        aa_key = aa_change[1:] if len(aa_change) > 1 else aa_change  # Strip first character

    # Look up associated lineages with type safety
    mapping = mutation_lookup.get(orf_key, aa_key)  # Type-safe access!
    if not mapping or not mapping.lineages:
        return ""

    # Rank by CDC circulation data (matching original ranked_matched_vars logic)
    matched_lineages = [
        str(lineage_name) for lineage_name in mapping.lineages
    ]  # Convert back to strings for comparison
    ranked_matches = []

    # First, add lineages that appear in CDC ranking
    ranked_matches = [
        ranked_variant for ranked_variant in ranked_variants if ranked_variant in matched_lineages
    ]

    # Then add remaining lineages not in CDC ranking
    ranked_matches.extend(
        [lineage for lineage in matched_lineages if lineage not in ranked_matches],
    )

    return ",".join(ranked_matches)


def process_temporal_comparison(
    window_pair: list[TimeWindow],
    comparison_data: TemporalComparisonData,
    config: AnalysisConfig,
) -> pl.DataFrame:
    """
    Process temporal comparison between two time windows.

    Replicates the weekly_out() function logic with proper tidyverse polars.
    """
    if len(window_pair) != EXPECTED_WINDOW_PAIR_SIZE:
        return pl.DataFrame()

    window1, window2 = window_pair

    # Calculate statistics for both windows (data already enriched with dates)
    window1_stats = calculate_window_variant_stats(
        comparison_data.covars_data,  # Already enriched with Collection Date
        comparison_data.ntcalls_data,  # Already enriched with Collection Date
        comparison_data.sample_metadata,  # For any additional joins needed
        window1,
        comparison_data.baseline_variants,
    )
    window2_stats = calculate_window_variant_stats(
        comparison_data.covars_data,  # Already enriched with Collection Date
        comparison_data.ntcalls_data,  # Already enriched with Collection Date
        comparison_data.sample_metadata,  # For any additional joins needed
        window2,
        comparison_data.baseline_variants,
    )

    if window1_stats.height == 0 or window2_stats.height == 0:
        return pl.DataFrame()

    # Find variants present in both windows for comparison
    return (
        window1_stats.join(window2_stats, on=["variant", "position"], how="inner", suffix="_w2")
        .with_columns(
            [(pl.col("abundance") - pl.col("abundance_w2")).abs().alias("abundance_diff")],
        )
        .filter(pl.col("abundance_diff") > config.delta_threshold)  # Significance threshold
        .filter(
            pl.col("total_samples_1k") + pl.col("total_samples_1k_w2") > 1,
        )  # 1k+ samples threshold
        .with_columns(
            [
                pl.col("variant")
                .map_elements(lambda x: extract_variant_components(x), return_dtype=pl.Object)
                .alias("components"),
            ],
        )
        .with_columns(
            [
                pl.col("components")
                .map_elements(lambda x: x["nt_change"], return_dtype=pl.Utf8)
                .alias("nt_change"),
                pl.col("components")
                .map_elements(lambda x: x["orf"], return_dtype=pl.Utf8)
                .alias("orf"),
                pl.col("components")
                .map_elements(lambda x: x["aa_change"], return_dtype=pl.Utf8)
                .alias("aa_change"),
            ],
        )
        .filter(pl.col("aa_change") != "")  # Only non-silent mutations
        .with_columns(
            [
                pl.struct(["nt_change", "orf", "aa_change"])
                .map_elements(
                    lambda row: rank_associated_variants(
                        {
                            "nt_change": row["nt_change"],
                            "orf": row["orf"],
                            "aa_change": row["aa_change"],
                        },
                        comparison_data.ranked_variants,
                        comparison_data.mutation_lookup,
                    ),
                    return_dtype=pl.Utf8,
                )
                .alias("associated_variants"),
            ],
        )
        .select(
            [
                "position",
                "nt_change",
                "orf",
                "aa_change",
                "associated_variants",
                pl.col("time_window").alias("time_span_1"),
                pl.col("total_read_count").alias("count_1"),
                pl.col("abundance").alias("abundance_1"),
                pl.col("total_coverage").alias("coverage_1"),
                pl.col("total_sample_count").alias("samples_1"),
                pl.col("total_samples_1k").alias("samples_1k_1"),
                pl.col("time_window_w2").alias("time_span_2"),
                pl.col("total_read_count_w2").alias("count_2"),
                pl.col("abundance_w2").alias("abundance_2"),
                pl.col("total_coverage_w2").alias("coverage_2"),
                pl.col("total_sample_count_w2").alias("samples_2"),
                pl.col("total_samples_1k_w2").alias("samples_1k_2"),
                "abundance_diff",
                pl.lit(datetime.now(tz=timezone.utc).strftime("%Y-%m-%d")).alias(
                    "process_end_date",
                ),
            ],
        )
    )


def reshape_to_long_format(comparison_results: pl.DataFrame) -> pl.DataFrame:
    """
    Reshape comparison results to long format.

    Replicates the pandas reshaping logic at the end of Comp_time_windows.py.
    """
    if comparison_results.height == 0:
        return pl.DataFrame()

    # Split into two dataframes for each time window
    window1_data = comparison_results.select(
        [
            "position",
            "nt_change",
            "orf",
            "aa_change",
            "associated_variants",
            "time_span_1",
            "count_1",
            "abundance_1",
            "coverage_1",
            "samples_1",
            "samples_1k_1",
            "process_end_date",
        ],
    ).rename(
        {
            "time_span_1": "time_span",
            "count_1": "count",
            "abundance_1": "abundance",
            "coverage_1": "total_coverage",
            "samples_1": "pos_samples",
            "samples_1k_1": "samples_1k_reads",
        },
    )

    window2_data = comparison_results.select(
        [
            "position",
            "nt_change",
            "orf",
            "aa_change",
            "associated_variants",
            "time_span_2",
            "count_2",
            "abundance_2",
            "coverage_2",
            "samples_2",
            "samples_1k_2",
            "process_end_date",
        ],
    ).rename(
        {
            "time_span_2": "time_span",
            "count_2": "count",
            "abundance_2": "abundance",
            "coverage_2": "total_coverage",
            "samples_2": "pos_samples",
            "samples_1k_2": "samples_1k_reads",
        },
    )

    # Combine and add time span start/end columns (matching original)
    return (
        pl.concat([window1_data, window2_data])
        .with_columns(
            [
                pl.col("time_span").str.split("--").list.get(0).alias("time_span_start"),
                pl.col("time_span").str.split("--").list.get(1).alias("time_span_end"),
            ],
        )
        .drop("time_span")
        .select(
            [
                "position",
                "nt_change",
                "orf",
                "aa_change",
                "associated_variants",
                "time_span_start",
                "time_span_end",
                "count",
                "abundance",
                "total_coverage",
                "pos_samples",
                "samples_1k_reads",
                "process_end_date",
            ],
        )
        .unique()  # Remove duplicates
        .sort(["position", "time_span_start"])
    )


def analyze_temporal_variants(
    config: AnalysisConfig,
    inputs: AnalysisInputs,
) -> None:
    """Execute complete temporal variant analysis pipeline."""
    logger.info("Starting temporal variant analysis")

    # Load all required data
    logger.info("Loading input data")
    raw_samrefiner_data = load_samrefiner_outputs(inputs.input_dir)

    if raw_samrefiner_data["covars"].height == 0:
        logger.error("No SAMRefiner covars data found")
        return

    # Use lazy loading for sample metadata
    sample_metadata = pl.scan_csv(inputs.sample_metadata_path, separator="\t").collect()

    # Enrich SAMRefiner data with collection dates from metadata
    logger.info("Enriching SAMRefiner data with collection dates")
    samrefiner_data = enrich_samrefiner_with_metadata(raw_samrefiner_data, sample_metadata)

    # Load lineage mutations and build lookup with lazy evaluation
    lineage_mutations = (
        pl.scan_csv(inputs.lineage_mutations_path, separator="\t").collect()
        if inputs.lineage_mutations_path.exists()
        else pl.DataFrame()
    )
    mutation_lookup = build_lineage_mutation_lookup(lineage_mutations)

    # Load baseline variants for exclusion
    baseline_data = load_baseline_jn1_variants(inputs.lineage_mutations_path)
    baseline_variants = baseline_data["variants"]

    # Load CDC variant ranking
    ranked_variants = load_cdc_variant_ranking(inputs.variant_surveillance_path)

    # Group samples by collection week
    logger.info("Grouping samples by collection week")
    weekly_data = (
        sample_metadata.with_columns(
            [
                pl.col("Collection Date")
                .str.to_date("%Y-%m-%d", strict=False)
                .alias("collection_date"),
                pl.col("Collection Date")
                .str.to_date("%Y-%m-%d", strict=False)
                .dt.truncate("1w")
                .alias("collection_week"),
            ],
        )
        .filter(pl.col("collection_date").is_not_null())
        .group_by("collection_week")
        .agg([pl.col("Sample").alias("accessions"), pl.len().alias("sample_count")])
        .filter(pl.col("sample_count") >= config.min_samples_per_window)
        .sort("collection_week")
    )

    if weekly_data.height == 0:
        logger.error("No samples found with valid collection dates")
        return

    # Create exact time window comparisons
    logger.info("Creating time window comparisons")
    comparison_window_sets = create_time_windows_exact(weekly_data)

    if not comparison_window_sets:
        logger.error("Insufficient weeks for temporal comparison")
        return

    # Bundle data for temporal comparison
    comparison_data = TemporalComparisonData(
        covars_data=samrefiner_data["covars"],
        ntcalls_data=samrefiner_data["ntcalls"],
        sample_metadata=sample_metadata,
        baseline_variants=baseline_variants,
        ranked_variants=ranked_variants,
        mutation_lookup=mutation_lookup,
    )

    # Process all window comparisons
    logger.info("Processing temporal comparisons")
    all_comparisons = []

    for window_pair in comparison_window_sets:
        comparison_result = process_temporal_comparison(
            window_pair,
            comparison_data,
            config,
        )

        # Early continue for empty results
        if comparison_result.height == 0:
            continue

        all_comparisons.append(comparison_result)

    if not all_comparisons:
        logger.warning("No significant temporal changes found")
        return

    # Combine all comparison results
    combined_results = pl.concat(all_comparisons)

    # Reshape to long format (matching original pandas logic)
    logger.info("Reshaping results to long format")
    long_format_results = reshape_to_long_format(combined_results)

    # Save results
    logger.info("Saving analysis results")
    save_temporal_variants(long_format_results, inputs.output_paths["temporal"])
    save_weekly_summary(weekly_data, inputs.output_paths["weekly"])
    save_significant_changes(combined_results, inputs.output_paths["significant"])

    logger.info("Temporal variant analysis completed successfully")


def save_temporal_variants(variants_df: pl.DataFrame, output_path: Path) -> None:
    """Save temporal variant analysis results in exact original format."""
    if variants_df.height == 0:
        logger.warning("No temporal variants to save")
        return

    variants_df.write_csv(output_path, separator="\t")
    logger.info(f"Saved {variants_df.height} temporal variant records to {output_path}")


def save_weekly_summary(weekly_data: pl.DataFrame, output_path: Path) -> None:
    """Save weekly sample summary data."""
    weekly_data.write_csv(output_path, separator="\t")
    logger.info(f"Saved weekly summary to {output_path}")


def save_significant_changes(changes_df: pl.DataFrame, output_path: Path) -> None:
    """Save significant temporal changes."""
    if changes_df.height == 0:
        logger.warning("No significant changes to save")
        return

    changes_df.write_csv(output_path, separator="\t")
    logger.info(f"Saved {changes_df.height} significant changes to {output_path}")


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze temporal variants from SAMRefiner outputs",
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing SAMRefiner outputs",
    )
    parser.add_argument("--sample-metadata", required=True, help="Sample metadata TSV file")
    parser.add_argument("--lineage-mutations", required=True, help="Lineage mutations TSV file")
    parser.add_argument(
        "--variant-surveillance",
        required=True,
        help="CDC surveillance data TSV file",
    )
    parser.add_argument(
        "--delta-threshold",
        type=float,
        default=0.02,
        help="Abundance change threshold",
    )
    parser.add_argument("--count-threshold", type=int, default=100, help="Minimum count threshold")
    parser.add_argument("--window-weeks", type=int, default=3, help="Time window size in weeks")
    parser.add_argument(
        "--output-temporal",
        required=True,
        help="Output path for temporal variants",
    )
    parser.add_argument("--output-weekly", required=True, help="Output path for weekly summary")
    parser.add_argument(
        "--output-significant",
        required=True,
        help="Output path for significant changes",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    config = AnalysisConfig(
        delta_threshold=args.delta_threshold,
        count_threshold=args.count_threshold,
        window_weeks=args.window_weeks,
    )

    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        logger.error(f"Input directory does not exist: {input_dir}")
        sys.exit(1)

    output_paths = {
        "temporal": Path(args.output_temporal),
        "weekly": Path(args.output_weekly),
        "significant": Path(args.output_significant),
    }

    # Ensure output directories exist
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)

    inputs = AnalysisInputs(
        input_dir=input_dir,
        sample_metadata_path=Path(args.sample_metadata),
        lineage_mutations_path=Path(args.lineage_mutations),
        variant_surveillance_path=Path(args.variant_surveillance),
        output_paths=output_paths,
    )

    try:
        analyze_temporal_variants(config, inputs)
    except (ValueError, OSError, pl.exceptions.PolarsError) as e:
        logger.error(f"Temporal variant analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
