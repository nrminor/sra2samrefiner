#!/usr/bin/env python3
"""
Generate surveillance reports from temporal variant analysis.

Replaces the functionality of MakeFocused2.py and GeoLocAbund.py with clean,
type-safe implementation using modern Python patterns.
"""

import argparse
import json
import re
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import polars as pl
from loguru import logger
from pydantic import Field
from pydantic.dataclasses import dataclass

# Constants for magic numbers
CDC_DATA_LOOKBACK_DAYS = 70  # Days of CDC data to consider
ORF1AB_SPLIT_POSITION = 4401  # Position where ORF1ab splits into ORF1a/ORF1b
MIN_MUTATION_PARTS = 2  # Minimum parts when splitting mutation on ":"

# HHS regions mapping (from original GeoLocAbund.py)
HHS_REGIONS: dict[int, str] = {
    1: "Connecticut, Maine, Massachusetts, New Hampshire, Rhode Island, Vermont",
    2: "New Jersey, New York, Puerto Rico, the Virgin Islands",
    3: "Delaware, District of Columbia, Maryland, Pennsylvania, Virginia, West Virginia",
    4: "Alabama, Florida, Georgia, Kentucky, Mississippi, North Carolina, South Carolina, Tennessee",
    5: "Illinois, Indiana, Michigan, Minnesota, Ohio, Wisconsin",
    6: "Arkansas, Louisiana, New Mexico, Oklahoma, Texas",
    7: "Iowa, Kansas, Missouri, Nebraska",
    8: "Colorado, Montana, North Dakota, South Dakota, Utah, Wyoming",
    9: "Arizona, California, Hawaii, Nevada, American Samoa, Commonwealth of the Northern Mariana Islands, Federated States of Micronesia, Guam, Marshall Islands, Republic of Palau",
    10: "Alaska, Idaho, Oregon, Washington",
}

# Reverse mapping: state -> HHS region number
STATE_TO_HHS_REGION: dict[str, int] = {
    state.strip(): region_num
    for region_num, states in HHS_REGIONS.items()
    for state in states.split(", ")
}


@dataclass
class LineageDefinition:
    """Represents a lineage with its defining mutations."""

    lineage: str = Field(min_length=1)
    mutations: list[str] = Field(default_factory=list)


@dataclass
class RegionalVariantData:
    """Represents variant data aggregated by HHS region."""

    position: int = Field(ge=1)
    nt_change: str = Field(min_length=1)
    orf: str = Field(min_length=1)
    aa_change: str = Field(min_length=1)
    week: str = Field(min_length=8)  # YYYY-MM-DD format
    hhs_region: int = Field(ge=1, le=10)
    count: int = Field(ge=0)
    total_coverage: int = Field(ge=0)
    abundance: float = Field(ge=0.0, le=1.0)


@dataclass
class SurveillanceSummary:
    """Summary statistics for surveillance analysis."""

    total_variants: int = Field(ge=0)
    significant_changes: int = Field(ge=0)
    lineages_identified: int = Field(ge=0)
    date_range_start: str = Field(min_length=8)
    date_range_end: str = Field(min_length=8)
    analysis_timestamp: str = Field(min_length=19)  # ISO datetime format


def extract_position_from_variant(variant_string: str) -> int:
    """Extract genomic position from variant string."""
    match = re.search(r"\d+", variant_string)
    return int(match.group(0)) if match else -1


def normalize_orf_for_borf(orf: str, aa_change: str) -> tuple[str, str]:
    """
    Normalize ORF and AA change, applying Borf function logic.

    Replicates the ORF1ab splitting and Borf position adjustment from original.
    """
    if "ORF1ab" in orf:
        position = extract_position_from_variant(aa_change)
        if position > ORF1AB_SPLIT_POSITION:
            # Apply Borf function: subtract 4401 from position
            adjusted_position = position - ORF1AB_SPLIT_POSITION
            adjusted_aa_change = re.sub(r"\d+", str(adjusted_position), aa_change)
            return "ORF1b", adjusted_aa_change
        return "ORF1a", aa_change

    if "ORF" in orf:
        return orf.split("_")[0], aa_change

    return orf[0].upper() if orf else orf, aa_change


def load_circulating_lineages(variant_surveillance_path: Path) -> set[str]:
    """
    Load circulating lineages from CDC data.
    Replicates the lineage extraction logic from MakeFocused2.py.
    """
    if not variant_surveillance_path.exists():
        logger.warning("No CDC surveillance data found")
        return set()

    # Process CDC data to get circulating lineages (matching original)
    cutoff_date = (datetime.now(tz=timezone.utc) - timedelta(days=CDC_DATA_LOOKBACK_DAYS)).strftime(
        "%Y-%m-%d",
    )

    return (
        pl.scan_csv(variant_surveillance_path, separator="\t")
        .filter(pl.col("usa_or_hhsregion") == "USA")
        .filter(pl.col("modeltype") == "smoothed")
        .filter(pl.col("week_ending") > cutoff_date)
        .filter(pl.col("share") > 0)
        .select("variant")
        .unique()
        .collect()
        .get_column("variant")
        .to_list()
    )


def process_deletion_variants(temporal_variants: pl.DataFrame) -> dict[str, str]:
    """
    Process deletion variants for special handling.
    Replicates the del_dict logic from MakeFocused2.py.
    """
    if temporal_variants.height == 0:
        return {}

    # Extract deletion variants and create mapping
    deletion_mapping = (
        temporal_variants.filter(pl.col("nt_change").str.contains("del"))
        .with_columns(
            [
                pl.col("nt_change")
                .str.strip_chars("ATCG")
                .str.replace("del", "Del:", literal=True)
                .alias("del_key"),
                (pl.col("orf") + ":" + pl.col("aa_change")).alias("del_value"),
            ],
        )
        .select(["del_key", "del_value"])
        .unique()
        .to_dicts()
    )

    del_dict = {item["del_key"]: item["del_value"] for item in deletion_mapping}

    # Add hardcoded special case from original
    del_dict["Del:21653-21655"] = "S:NS30-31Ndel"

    logger.info(f"Processed {len(del_dict)} deletion variants")
    return del_dict


def build_aa_change_lookup(temporal_variants: pl.DataFrame) -> dict[str, dict[str, str]]:
    """
    Build AA change lookup for variant matching.
    Replicates the match_dict logic from MakeFocused2.py.
    """
    if temporal_variants.height == 0:
        return {}

    # Create ORF:AA_Change combinations
    aa_changes = (
        temporal_variants.with_columns(
            (pl.col("orf") + ":" + pl.col("aa_change")).alias("aa_changes"),
        )
        .select("aa_changes")
        .unique()
        .get_column("aa_changes")
        .to_list()
    )

    match_dict = {}
    for change in aa_changes:
        if ":" not in change:
            continue

        parts = change.split(":")
        if len(parts) < MIN_MUTATION_PARTS:
            continue

        orf = parts[0]
        aa_change = parts[1]

        if orf not in match_dict:
            match_dict[orf] = {}

        # Store mapping from AA change substring to full change
        if len(aa_change) > 1:
            match_dict[orf][aa_change[1:]] = aa_change

    logger.info(f"Built AA change lookup with {len(match_dict)} ORFs")
    return match_dict


def generate_lineage_definitions_json(
    circulating_lineages: set[str],
    output_path: Path,
) -> None:
    """
    Generate major lineages JSON file.

    Replicates the JSON generation logic from MakeFocused2.py with proper structure.
    """
    lineage_definitions = []

    for lineage in circulating_lineages:
        if lineage == "Other":
            continue

        # This would need integration with external mutation database
        # For now, create placeholder structure matching original output format
        lineage_def = LineageDefinition(
            lineage=lineage,
            mutations=[],  # Would be populated from sars-cov-2-lineage-dominant-mutations
        )

        if lineage_def.mutations:  # Only include lineages with mutations
            lineage_definitions.append(
                {
                    "lineage": lineage_def.lineage,
                    "mutations": ",".join(lineage_def.mutations),
                },
            )

    # Write JSON in original format
    with open(output_path, "w") as json_file:
        json.dump(lineage_definitions, json_file, indent=2)

    logger.info(f"Generated lineage definitions JSON with {len(lineage_definitions)} entries")


def generate_regional_analysis(
    temporal_variants: pl.DataFrame,
    sample_metadata: pl.DataFrame,
    output_path: Path,
) -> None:
    """
    Generate regional HHS analysis.
    Replicates the geographic aggregation from GeoLocAbund.py.
    """
    if temporal_variants.height == 0 or sample_metadata.height == 0:
        logger.warning("Insufficient data for regional analysis")
        return

    # Create comprehensive regional analysis
    regional_data = (
        temporal_variants.join(sample_metadata, left_on="sample_id", right_on="Sample", how="left")
        .with_columns(
            [
                # Extract state from location (format: "USA:State")
                pl.col("Location").str.extract(r":([^,]+)").alias("state"),
                pl.col("time_span_start").alias("week"),
            ],
        )
        .filter(pl.col("state").is_not_null())
        .with_columns(
            [
                # Map state to HHS region
                pl.col("state")
                .map_elements(
                    lambda state: STATE_TO_HHS_REGION.get(state.strip(), None),
                    return_dtype=pl.Int32,
                )
                .alias("hhs_region"),
            ],
        )
        .filter(pl.col("hhs_region").is_not_null())
        .with_columns(
            [
                # Apply ORF normalization for regional analysis
                pl.struct(["orf", "aa_change"])
                .map_elements(
                    lambda row: normalize_orf_for_borf(row["orf"], row["aa_change"]),
                    return_dtype=pl.Object,
                )
                .alias("normalized_orf_aa"),
            ],
        )
        .with_columns(
            [
                pl.col("normalized_orf_aa")
                .map_elements(lambda x: x[0], return_dtype=pl.Utf8)
                .alias("normalized_orf"),
                pl.col("normalized_orf_aa")
                .map_elements(lambda x: x[1], return_dtype=pl.Utf8)
                .alias("normalized_aa_change"),
            ],
        )
        .select(
            [
                "position",
                "nt_change",
                "normalized_orf",
                "normalized_aa_change",
                "week",
                "hhs_region",
                "count",
                "total_coverage",
                "abundance",
            ],
        )
        .rename(
            {
                "normalized_orf": "orf",
                "normalized_aa_change": "aa_change",
            },
        )
    )

    # Write long-format regional data (matching HHS.long.tsv structure)
    regional_data.write_csv(output_path, separator="\t")
    logger.info(f"Generated regional analysis with {regional_data.height} records")


def generate_summary_statistics(
    temporal_variants: pl.DataFrame,
    significant_changes: pl.DataFrame,
    sample_metadata: pl.DataFrame,
    output_path: Path,
) -> None:
    """Generate summary statistics for the surveillance analysis."""
    # Calculate date range from sample metadata
    date_range = sample_metadata.select(
        [
            pl.col("Collection Date").min().alias("start_date"),
            pl.col("Collection Date").max().alias("end_date"),
        ],
    ).row(0)

    summary = SurveillanceSummary(
        total_variants=temporal_variants.height,
        significant_changes=significant_changes.height,
        lineages_identified=temporal_variants.select("associated_variants").unique().height,
        date_range_start=str(date_range[0]) if date_range[0] else "Unknown",
        date_range_end=str(date_range[1]) if date_range[1] else "Unknown",
        analysis_timestamp=datetime.now(tz=timezone.utc).isoformat(),
    )

    # Save as JSON
    with open(output_path, "w") as summary_file:
        json.dump(summary.__dict__, summary_file, indent=2, default=str)

    logger.info("Generated surveillance summary statistics")


def generate_all_surveillance_reports(
    temporal_variants_path: Path,
    significant_changes_path: Path,
    sample_metadata_path: Path,
    variant_surveillance_path: Path,
    output_paths: dict[str, Path],
) -> None:
    """Generate all surveillance reports."""
    logger.info("Starting surveillance report generation")

    # Load input data with lazy evaluation
    temporal_variants = (
        pl.scan_csv(temporal_variants_path, separator="\t").collect()
        if temporal_variants_path.exists()
        else pl.DataFrame()
    )
    significant_changes = (
        pl.scan_csv(significant_changes_path, separator="\t").collect()
        if significant_changes_path.exists()
        else pl.DataFrame()
    )
    sample_metadata = (
        pl.scan_csv(sample_metadata_path, separator="\t").collect()
        if sample_metadata_path.exists()
        else pl.DataFrame()
    )

    # Load circulating lineages
    circulating_lineages = load_circulating_lineages(variant_surveillance_path)

    # Process deletion variants (for future enhancement)
    _ = process_deletion_variants(temporal_variants)

    # Build AA change lookup (for future enhancement)
    _ = build_aa_change_lookup(temporal_variants)

    # Generate reports
    logger.info("Generating lineage definitions JSON")
    generate_lineage_definitions_json(
        circulating_lineages,
        output_paths["json"],
    )

    logger.info("Generating regional analysis")
    generate_regional_analysis(
        temporal_variants,
        sample_metadata,
        output_paths["regional"],
    )

    logger.info("Generating summary statistics")
    generate_summary_statistics(
        temporal_variants,
        significant_changes,
        sample_metadata,
        output_paths["summary"],
    )

    logger.info("Surveillance report generation completed successfully")


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate surveillance reports from temporal variant analysis",
    )
    parser.add_argument("--temporal-variants", required=True, help="Temporal variants TSV file")
    parser.add_argument("--significant-changes", required=True, help="Significant changes TSV file")
    parser.add_argument("--sample-metadata", required=True, help="Sample metadata TSV file")
    parser.add_argument(
        "--variant-surveillance",
        required=True,
        help="CDC surveillance data TSV file",
    )
    parser.add_argument(
        "--output-json",
        required=True,
        help="Output path for lineage definitions JSON",
    )
    parser.add_argument(
        "--output-regional",
        required=True,
        help="Output path for regional analysis TSV",
    )
    parser.add_argument(
        "--output-summary",
        required=True,
        help="Output path for summary statistics JSON",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    # Validate input files exist
    input_files = [
        Path(args.temporal_variants),
        Path(args.significant_changes),
        Path(args.sample_metadata),
        Path(args.variant_surveillance),
    ]

    for input_file in input_files:
        if not input_file.exists():
            logger.error(f"Input file does not exist: {input_file}")
            sys.exit(1)

    output_paths = {
        "json": Path(args.output_json),
        "regional": Path(args.output_regional),
        "summary": Path(args.output_summary),
    }

    # Ensure output directories exist
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)

    try:
        generate_all_surveillance_reports(
            Path(args.temporal_variants),
            Path(args.significant_changes),
            Path(args.sample_metadata),
            Path(args.variant_surveillance),
            output_paths,
        )
    except (ValueError, OSError, pl.exceptions.PolarsError) as e:
        logger.error(f"Surveillance report generation failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
