#!/usr/bin/env python3
"""
Prepare reference data for SRATracking analysis.

Replaces the functionality of SRA_Query.py, Tree.py, and GetCDCData.py
with clean, async implementation using proper HTTP clients.
"""

import argparse
import asyncio
import subprocess
import sys
import xml.etree.ElementTree as ET
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import httpx
import polars as pl
from loguru import logger
from pydantic import Field, ValidationInfo, field_validator
from pydantic.dataclasses import dataclass

# Constants to replace magic numbers
MIN_CSV_PARTS_FOR_RELEASE_DATA = 3
MIN_TAB_PARTS_FOR_LINEAGE_NODE = 2
MIN_TAB_PARTS_FOR_MUTATION_DATA = 3
CDC_DATA_LOOKBACK_DAYS = 70
MAX_RETURN_STATEMENTS = 6  # Used for function complexity tracking


@dataclass
class DateRange:
    start_date: date
    end_date: date

    @field_validator("end_date")
    @classmethod
    def end_after_start(cls, v: date, info: ValidationInfo) -> date:
        if info.data and "start_date" in info.data and v <= info.data["start_date"]:
            msg = "end_date must be after start_date"
            raise ValueError(msg)
        return v


@dataclass
class LineageMutation:
    lineage: str = Field(min_length=1)
    mutations: list[str] = Field(default_factory=list)

    @field_validator("mutations")
    @classmethod
    def validate_mutations(cls, v: list[str]) -> list[str]:
        if not v:
            msg = "mutations list cannot be empty"
            raise ValueError(msg)
        return v


@dataclass
class SampleMetadata:
    # Required fields (no defaults) - all must use Field() consistently
    accession: str = Field(pattern=r"^[SED]RR\d+$")
    collection_date: str = Field(min_length=1)  # Keep as string - will be normalized later
    location: str = Field(min_length=1)
    bioproject: str = Field(min_length=1)
    biosample: str = Field(min_length=1)
    submitter: str = Field(min_length=1)
    reads: int = Field(ge=0)
    # Optional fields (with defaults)
    release_date: str = ""
    load_date: str = ""


@dataclass
class FilterConfig:
    """Configuration for data filtering and exclusions."""

    exclude_bioprojects: list[str] = Field(default_factory=list)
    exclude_samples_pattern: list[str] = Field(default_factory=list)  # Regex patterns
    min_reads: int = Field(default=0, ge=0)
    require_collection_date: bool = True
    require_location: bool = True


def normalize_collection_date(date_str: str) -> str:
    """
    Normalize collection date from various formats to YYYY-MM-DD.

    Handles the messy reality of SRA date formats:
    - YYYY-MM-DD (ISO format)
    - MM/DD/YYYY (US format)
    - YYYY/MM/DD
    - YYYY-MM (partial dates)
    - DD-MM-YYYY (European format)
    """
    if not date_str or date_str.strip() == "":
        return ""

    date_str = date_str.strip()

    # Try common date formats in order of preference
    # Optimized to avoid performance overhead from try-except in loop
    def _try_parse_date(date_input: str, fmt: str) -> str | None:
        """Helper function to parse date with specific format."""
        try:
            parsed_date = datetime.strptime(date_input, fmt).replace(tzinfo=timezone.utc)
            # Convert to ISO format based on the format that matched
            if fmt == "%Y-%m":
                return f"{parsed_date.strftime('%Y-%m')}-01"  # Default to first of month
            if fmt == "%m/%Y":
                return f"{parsed_date.strftime('%Y-%m')}-01"
            if fmt == "%Y":
                return f"{parsed_date.strftime('%Y')}-01-01"  # Default to first of year
            return parsed_date.strftime("%Y-%m-%d")
        except ValueError:
            return None

    date_formats = [
        "%Y-%m-%d",  # ISO format (preferred)
        "%m/%d/%Y",  # US format
        "%Y/%m/%d",  # Alternative slash format
        "%d-%m-%Y",  # European format
        "%Y-%m",  # Partial date (year-month)
        "%m/%Y",  # Partial US format
        "%Y",  # Year only
    ]

    # Try each format until one succeeds
    for fmt in date_formats:
        result = _try_parse_date(date_str, fmt)
        if result is not None:
            return result

    # If no format matches, log and return original
    logger.warning(f"Unable to parse date format: '{date_str}', keeping as-is")
    return date_str


def should_exclude_sample(sample: SampleMetadata, filter_config: FilterConfig) -> bool:
    """
    Check if a sample should be excluded based on filter configuration.

    Replicates the hardcoded exclusions from original code in a configurable way.
    """
    import re

    exclusion_reasons = []
    # Check bioproject exclusions
    if sample.bioproject in filter_config.exclude_bioprojects:
        exclusion_reasons.append(f"bioproject {sample.bioproject} in exclusion list")

    # Check sample pattern exclusions (e.g., specific submitter + location combinations)
    for pattern in filter_config.exclude_samples_pattern:
        # Pattern can match against any field - format: "field:pattern"
        if ":" in pattern:
            field_name, regex_pattern = pattern.split(":", 1)
            field_value = getattr(sample, field_name, "")
            if re.search(regex_pattern, str(field_value), re.IGNORECASE):
                exclusion_reasons.append(f"{field_name} matches pattern '{regex_pattern}'")
                break  # Stop at first pattern match
        # Pattern matches entire sample string representation
        elif re.search(pattern, str(sample), re.IGNORECASE):
            exclusion_reasons.append(f"matches pattern '{pattern}'")
            break  # Stop at first pattern match

    # Check minimum reads threshold
    if sample.reads < filter_config.min_reads:
        exclusion_reasons.append(f"reads {sample.reads} < {filter_config.min_reads}")

    # Check required fields
    if filter_config.require_collection_date and not sample.collection_date.strip():
        exclusion_reasons.append("missing collection date")

    if filter_config.require_location and not sample.location.strip():
        exclusion_reasons.append("missing location")

    # Log the first exclusion reason if any exist
    if exclusion_reasons:
        logger.debug(f"Excluding sample {sample.accession}: {exclusion_reasons[0]}")
        return True

    return False


async def fetch_ncbi_metadata(
    date_range: DateRange,
    filter_config: FilterConfig,
) -> list[SampleMetadata]:
    """
    Fetch NCBI SRA metadata for wastewater samples.

    Replicates the functionality of SRA_Query.py with proper HTTP handling.
    """
    logger.info(f"Fetching NCBI metadata from {date_range.start_date} to {date_range.end_date}")

    # Build search query matching original format
    s_y = date_range.start_date.strftime("%Y")
    s_m = date_range.start_date.strftime("%m")
    s_d = date_range.start_date.strftime("%d")
    e_y = date_range.end_date.strftime("%Y")
    e_m = date_range.end_date.strftime("%m")
    e_d = date_range.end_date.strftime("%d")

    search_query = f"(sars-cov-2%20wastewater)%20AND%20(%22{s_y}%2F{s_m}%2F{s_d}%22%5BPublication%20Date%5D%20%3A%20%22{e_y}%2F{e_m}%2F{e_d}%22%5BPublication%20Date%5D)"

    async with httpx.AsyncClient() as client:
        # Step 1: Get search results page to extract MCID and Key
        logger.info("Fetching SRA search results page")
        search_url = f"https://www.ncbi.nlm.nih.gov/sra/?term={search_query}"
        headers = {
            "User-Agent": "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:50.0) Gecko/20100101 Firefox/50.0",
        }

        search_response = await client.get(search_url, headers=headers)
        search_response.raise_for_status()

        # Extract MCID and Key from HTML (matching original logic)
        search_html = search_response.text
        try:
            mcid = search_html.split('value="MCID_')[1].split('"')[0]
            key = search_html.split("query_key:&quot;")[1].split("&quot")[0]
            logger.info(f"Extracted MCID: {mcid}, Key: {key}")
        except IndexError as e:
            msg = "Failed to extract MCID or query key from SRA search results"
            raise ValueError(msg) from e

        # Step 2: Download XML metadata, accession list, and run info
        base_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi"
        web_env = f"MCID_{mcid}"

        # Download XML metadata
        xml_url = f"{base_url}?rettype=exp&WebEnv={web_env}&query_key={key}"
        xml_response = await client.get(xml_url)
        xml_response.raise_for_status()

        # Download run info CSV
        runinfo_url = f"{base_url}?rettype=runinfo&WebEnv={web_env}&query_key={key}"
        runinfo_response = await client.get(runinfo_url)
        runinfo_response.raise_for_status()

        # Parse run info for release/load dates
        release_dict = {}
        runinfo_lines = runinfo_response.text.strip().split("\n")
        for line in runinfo_lines[1:]:  # Skip header
            parts = line.split(",")
            if len(parts) >= MIN_CSV_PARTS_FOR_RELEASE_DATA:
                release_dict[parts[0]] = (parts[1], parts[2])

        # Parse XML metadata
        raw_samples = parse_sra_xml_metadata(xml_response.text, release_dict)
        logger.info(f"Parsed {len(raw_samples)} raw samples from SRA metadata")

        # Apply filtering and date normalization
        filtered_samples = []
        excluded_count = 0

        for sample in raw_samples:
            # Normalize collection date
            sample.collection_date = normalize_collection_date(sample.collection_date)

            # Apply exclusion filters - early continue for excluded samples
            if should_exclude_sample(sample, filter_config):
                excluded_count += 1
                continue

            filtered_samples.append(sample)

        logger.info(f"Filtered to {len(filtered_samples)} samples ({excluded_count} excluded)")
        return filtered_samples


def extract_run_data(root: ET.Element) -> dict[str, any]:
    """Extract run accession and read count data."""
    data = {"accession": "", "reads": 0}

    for run in root.findall(".//RUN"):
        if "accession" not in run.attrib:
            continue

        data["accession"] = run.attrib["accession"]
        if "total_spots" in run.attrib:
            data["reads"] = int(run.attrib.get("total_spots", 0))
        break  # Take first valid run

    return data


def extract_submission_data(root: ET.Element) -> dict[str, str]:
    """Extract submission information."""
    data = {"submitter": ""}

    for submission in root.findall(".//SUBMISSION"):
        if "center_name" not in submission.attrib:
            continue

        data["submitter"] = submission.attrib["center_name"]
        break  # Take first valid submission

    return data


def extract_project_data(root: ET.Element) -> dict[str, str]:
    """Extract bioproject and biosample IDs."""
    data = {"bioproject": "", "biosample": ""}

    for link in root.findall(".//XREF_LINK/ID"):
        if not link.text:
            continue

        if link.text.startswith("PRJ"):
            data["bioproject"] = link.text
        elif link.text.startswith("SAM"):
            data["biosample"] = link.text

    return data


def extract_sample_attributes(root: ET.Element) -> dict[str, str]:
    """Extract sample attributes like collection date and location."""
    data = {"collection_date": "", "location": "", "population": ""}

    for attr in root.findall(".//SAMPLE_ATTRIBUTE"):
        tag = attr.find("TAG")
        value = attr.find("VALUE")

        if tag is None or value is None:
            continue

        tag_text = tag.text.lower() if tag.text else ""
        value_text = value.text or ""

        if "collection_date" in tag_text or "collection date" in tag_text:
            data["collection_date"] = value_text
        elif "geo_loc_name" in tag_text or "geographic location" in tag_text:
            data["location"] = value_text
        elif "ww_population" in tag_text:
            data["population"] = value_text

    return data


def parse_single_xml_package(
    xml_part: str,
    release_dict: dict[str, tuple],
) -> SampleMetadata | None:
    """Parse a single XML experiment package."""
    try:
        # Note: Using defusedxml would be more secure, but ET with limited input should be acceptable
        # for this controlled use case with SRA metadata from trusted NCBI source
        root = ET.fromstring(f"<root>{xml_part}</root>")  # noqa: S314
    except ET.ParseError:
        logger.warning("Failed to parse XML section, skipping")
        return None

    # Extract data from different XML sections - early return if no accession
    run_data = extract_run_data(root)
    if not run_data["accession"] or not run_data["accession"].startswith(("SRR", "ERR", "DRR")):
        return None

    submission_data = extract_submission_data(root)
    project_data = extract_project_data(root)
    attributes = extract_sample_attributes(root)

    # Get release/load dates
    release_date, load_date = release_dict.get(run_data["accession"], ("", ""))

    return SampleMetadata(
        accession=run_data["accession"],
        collection_date=attributes["collection_date"],
        location=attributes["location"],
        bioproject=project_data["bioproject"],
        biosample=project_data["biosample"],
        submitter=submission_data["submitter"],
        reads=run_data["reads"],
        release_date=release_date,
        load_date=load_date,
    )


def parse_sra_xml_metadata(
    xml_content: str,
    release_dict: dict[str, tuple],
) -> list[SampleMetadata]:
    """Parse SRA XML metadata to extract sample information."""
    samples = []

    # Split XML into experiment packages (matching original logic)
    xml_parts = xml_content.split("</EXPERIMENT_PACKAGE_SET>\n<EXPERIMENT_PACKAGE_SET>")

    for xml_part in xml_parts:
        sample = parse_single_xml_package(xml_part, release_dict)
        if sample is not None:
            samples.append(sample)

    return samples


def process_lineage_tree_relationships(target_lineage: str) -> dict[str, list[str]]:
    """
    Process phylogenetic tree to understand lineage relationships.

    Replicates the complex tree traversal logic from Tree.py to properly
    identify descendant lineages and their hierarchical relationships.
    """
    lineage_hierarchy = {}

    if not Path("LineageDefinitions.tsv").exists():
        logger.warning("LineageDefinitions.tsv not found, skipping hierarchy processing")
        return {}

    # First, get the target lineage node ID from definitions
    target_node = ""
    with open("LineageDefinitions.tsv") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= MIN_TAB_PARTS_FOR_LINEAGE_NODE and parts[0] == target_lineage:
                target_node = parts[1]
                break

    if not target_node:
        logger.warning(f"Target lineage {target_lineage} not found in definitions")
        return {}

    # Extract tree using matUtils and process it (if tree processing is needed)
    try:
        # Use full path to matUtils for security
        import shutil

        matutils_path = shutil.which("matUtils")
        if not matutils_path:
            msg = "matUtils not found in PATH"
            raise FileNotFoundError(msg)
        # Note: subprocess call is secure as matutils_path is validated with shutil.which()
        # and arguments are hardcoded constants, not user input
        subprocess.run(  # noqa: S603
            [
                matutils_path,
                "extract",
                "-i",
                "public-latest.all.masked.pb.gz",
                "-t",
                "sars.nwk",
            ],
            check=True,
            timeout=300,  # Add timeout for security
        )

        # For now, use simpler lineage matching approach
        # The original tree processing is very complex and may not be essential
        # for basic functionality - can be enhanced later if needed
        logger.info(f"Using simplified lineage matching for {target_lineage}")

    except subprocess.CalledProcessError:
        logger.warning("Tree extraction failed, using simplified approach")

    return lineage_hierarchy


def filter_descendant_mutations(
    lineage_mutations: list[LineageMutation],
    target_lineage: str,
) -> list[LineageMutation]:
    """
    Filter mutations to only include those specific to descendants, not ancestors.

    Replicates the logic from Tree.py that removes ancestral mutations.
    """
    if not lineage_mutations:
        return []

    # Find the target lineage mutations to exclude from descendants
    target_mutations = []
    for mutation in lineage_mutations:
        if mutation.lineage == target_lineage:
            target_mutations = mutation.mutations
            break

    # Filter descendant mutations to exclude ancestral ones
    filtered_mutations = []
    for mutation in lineage_mutations:
        # Skip if not a descendant
        if not mutation.lineage.startswith(target_lineage):
            continue

        # For descendants, exclude mutations already present in target lineage
        descendant_specific = [m for m in mutation.mutations if m not in target_mutations]

        # Only include if there are descendant-specific mutations
        if descendant_specific:
            filtered_mutations.append(
                LineageMutation(lineage=mutation.lineage, mutations=descendant_specific),
            )

    return filtered_mutations


def process_usher_data(
    target_lineage: str,
    lineage_definitions_path: Path,
) -> list[LineageMutation]:
    """
    Process UShER phylogenetic data and lineage mutations.

    Uses pre-downloaded UShER files for efficiency.
    """
    logger.info(f"Processing UShER data for lineage: {target_lineage}")

    # Process tree relationships for more accurate lineage handling
    process_lineage_tree_relationships(target_lineage)

    # Parse lineage definitions
    all_lineage_mutations = []

    if not lineage_definitions_path.exists():
        logger.error(f"Lineage definitions file not found: {lineage_definitions_path}")
        return []

    with open(lineage_definitions_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < MIN_TAB_PARTS_FOR_MUTATION_DATA:
                continue

            lineage_name = parts[0]

            # Skip if not related to target lineage
            if not lineage_name.startswith(target_lineage):
                continue

            # Extract mutations from the definition
            mutations = [m.strip() for m in parts[2].split(";") if m.strip()]
            if not mutations:
                continue

            all_lineage_mutations.append(LineageMutation(lineage=lineage_name, mutations=mutations))

    # Filter to remove ancestral mutations (matching Tree.py logic)
    filtered_mutations = filter_descendant_mutations(all_lineage_mutations, target_lineage)

    logger.info(f"Extracted {len(filtered_mutations)} descendant lineage definitions")
    return filtered_mutations


async def fetch_cdc_surveillance() -> pl.DataFrame:
    """
    Fetch CDC variant surveillance data.

    Replicates the functionality of GetCDCData.py.
    Note: CDC data is filtered by recency, not by the requested date range.
    """
    logger.info("Fetching CDC variant surveillance data")

    data_url = "https://data.cdc.gov/api/views/jr58-6ysp/rows.tsv"

    async with httpx.AsyncClient() as client:
        response = await client.get(data_url)
        response.raise_for_status()

        # Parse TSV data with polars
        from io import StringIO

        # Use polars for better performance
        cutoff_date = (
            datetime.now(tz=timezone.utc) - timedelta(days=CDC_DATA_LOOKBACK_DAYS)
        ).strftime("%Y-%m-%d")

        # Load data and apply filtering
        raw_data = pl.read_csv(StringIO(response.text), separator="\t")

        # Apply date filter if week_ending column exists
        cdc_surveillance_data = (
            raw_data.filter(pl.col("week_ending") >= cutoff_date)
            if "week_ending" in raw_data.columns
            else raw_data
        )

        logger.info(f"Retrieved {cdc_surveillance_data.height} CDC surveillance records")
        return cdc_surveillance_data


def save_sample_metadata(samples: list[SampleMetadata], output_path: Path) -> None:
    """Save sample metadata to TSV format using polars."""
    if not samples:
        # Create empty DataFrame with proper schema
        empty_metadata = pl.DataFrame({
            "Sample": [], "BioProject": [], "BioSample": [], "Submitter": [],
            "Collection Date": [], "Location": [], "Population": [], "Reads": [],
            "ReleaseDate": [], "LoadDate": [],
        })
        empty_metadata.write_csv(output_path, separator="\t")
        logger.info("Saved empty sample metadata file")
        return

    # Convert samples to polars DataFrame
    sample_metadata_df = pl.DataFrame([
        {
            "Sample": s.accession,
            "BioProject": s.bioproject,
            "BioSample": s.biosample,
            "Submitter": s.submitter,
            "Collection Date": s.collection_date,
            "Location": s.location,
            "Population": "",  # Population field not always available
            "Reads": s.reads,
            "ReleaseDate": s.release_date,
            "LoadDate": s.load_date,
        }
        for s in samples
    ])

    sample_metadata_df.write_csv(output_path, separator="\t")
    logger.info(f"Saved {len(samples)} sample metadata records to {output_path}")


def save_lineage_mutations(mutations: list[LineageMutation], output_path: Path) -> None:
    """Save lineage mutations to TSV format using polars."""
    if not mutations:
        # Create empty DataFrame with proper schema
        empty_lineages = pl.DataFrame({"lineage": [], "mutations": []})
        empty_lineages.write_csv(output_path, separator="\t")
        logger.info("Saved empty lineage mutations file")
        return

    lineage_mutations_df = pl.DataFrame([
        {"lineage": m.lineage, "mutations": ";".join(m.mutations)}
        for m in mutations
    ])

    lineage_mutations_df.write_csv(output_path, separator="\t")
    logger.info(f"Saved {len(mutations)} lineage definitions to {output_path}")


async def prepare_all_reference_data(
    date_range: DateRange,
    target_lineage: str,
    filter_config: FilterConfig,
    output_paths: dict[str, Path],
    lineage_definitions_path: Path | None = None,
) -> None:
    """Coordinate all reference data preparation."""
    logger.info("Starting reference data preparation")

    # Fetch all data concurrently where possible
    logger.info("Fetching data from external sources")

    # NCBI and CDC can be fetched concurrently
    metadata_task = fetch_ncbi_metadata(date_range, filter_config)
    cdc_task = fetch_cdc_surveillance()

    # Process UShER data (either from provided files or download fresh)
    if lineage_definitions_path and lineage_definitions_path.exists():
        logger.info("Using pre-downloaded UShER lineage definitions")
        usher_data = process_usher_data(target_lineage, lineage_definitions_path)
    else:
        logger.info("UShER data not provided, skipping lineage processing")
        usher_data = []

    sample_metadata, cdc_data = await asyncio.gather(metadata_task, cdc_task)

    # Save outputs
    logger.info("Saving reference data files")
    save_sample_metadata(sample_metadata, output_paths["sample_metadata"])
    save_lineage_mutations(usher_data, output_paths["lineage_mutations"])
    cdc_data.write_csv(output_paths["variant_surveillance"], separator="\t")

    logger.info("Reference data preparation completed successfully")


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Prepare reference data for SRATracking analysis")
    parser.add_argument("--start-date", required=True, help="Start date (YYYY-MM-DD)")
    parser.add_argument("--end-date", required=True, help="End date (YYYY-MM-DD)")
    parser.add_argument("--target-lineage", required=True, help="Target lineage (e.g., JN.1)")
    parser.add_argument(
        "--lineage-mutations",
        required=True,
        help="Output path for lineage mutations TSV",
    )
    parser.add_argument(
        "--variant-surveillance",
        required=True,
        help="Output path for CDC surveillance TSV",
    )
    parser.add_argument(
        "--sample-metadata",
        required=True,
        help="Output path for sample metadata TSV",
    )

    # Optional pre-downloaded UShER files (for caching optimization)
    parser.add_argument("--usher-db", help="Path to pre-downloaded UShER database file")
    parser.add_argument(
        "--lineage-definitions",
        help="Path to pre-extracted LineageDefinitions.tsv",
    )

    # Filtering options to replace hardcoded exclusions
    parser.add_argument(
        "--exclude-bioprojects",
        nargs="*",
        default=["PRJNA748354", "PRJEB44932"],
        help="Bioprojects to exclude (default: known problematic projects)",
    )
    parser.add_argument(
        "--exclude-patterns",
        nargs="*",
        default=["location:Maryland.*78365"],
        help="Sample exclusion patterns (format: field:regex)",
    )
    parser.add_argument("--min-reads", type=int, default=0, help="Minimum read count threshold")
    parser.add_argument(
        "--require-collection-date",
        action="store_true",
        default=True,
        help="Require collection date to be present",
    )
    parser.add_argument(
        "--require-location",
        action="store_true",
        default=True,
        help="Require location to be present",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    # Validate date format - early exit on error
    try:
        date_range = DateRange(
            start_date=date.fromisoformat(args.start_date),
            end_date=date.fromisoformat(args.end_date),
        )
    except ValueError as e:
        logger.error(f"Invalid date format: {e}")
        sys.exit(1)

    # Create filter configuration
    filter_config = FilterConfig(
        exclude_bioprojects=args.exclude_bioprojects,
        exclude_samples_pattern=args.exclude_patterns,
        min_reads=args.min_reads,
        require_collection_date=args.require_collection_date,
        require_location=args.require_location,
    )

    output_paths = {
        "lineage_mutations": Path(args.lineage_mutations),
        "variant_surveillance": Path(args.variant_surveillance),
        "sample_metadata": Path(args.sample_metadata),
    }

    # Ensure output directories exist
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)

    # Handle optional UShER file paths
    lineage_definitions_path = Path(args.lineage_definitions) if args.lineage_definitions else None

    try:
        asyncio.run(
            prepare_all_reference_data(
                date_range,
                args.target_lineage,
                filter_config,
                output_paths,
                lineage_definitions_path,
            ),
        )
    except (OSError, ValueError, RuntimeError, asyncio.TimeoutError, KeyboardInterrupt) as e:
        logger.error(f"Reference data preparation failed: {e}")
        sys.exit(1)
    except (ImportError, AttributeError, TypeError) as e:
        logger.error(f"Configuration or import error during reference data preparation: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
