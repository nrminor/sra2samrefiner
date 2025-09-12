#!/usr/bin/env python3
"""
Test the updated single-end trimming logic
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'bin'))

from trim_aligned_reads import ReadCategory, TrimPolicy

def test_single_end_logic():
    """Test that ReadCategory.OTHER now applies single-end trimming."""
    
    print("=== TESTING SINGLE-END TRIMMING LOGIC ===\n")
    
    # Create a test policy
    policy = TrimPolicy(
        merged_left=10, merged_right=15,
        r1_left=20, r2_right=25,  
        single_left=30, single_right=35,  # NEW parameters
        min_len=20
    )
    
    print(f"Test policy: {policy}\n")
    
    # Test different read categories
    test_cases = [
        ("MERGED_read123", ReadCategory.MERGED, (10, 15)),
        ("UNMERGED_read456/1", ReadCategory.UNMERGED_R1, (20, 0)), 
        ("UNMERGED_read789/2", ReadCategory.UNMERGED_R2, (0, 25)),
        ("SRR12345678.1", ReadCategory.OTHER, (30, 35)),  # Single-end read
        ("nanopore_read_001", ReadCategory.OTHER, (30, 35)),  # Nanopore read
        ("pacbio_m54321", ReadCategory.OTHER, (30, 35)),  # PacBio read
    ]
    
    print("Testing read classification and trimming:")
    for read_name, expected_category, expected_trim in test_cases:
        # Test classification
        actual_category = ReadCategory.classify(read_name)
        assert actual_category == expected_category, \
            f"Classification failed for '{read_name}': expected {expected_category}, got {actual_category}"
        
        # Test trim extents
        actual_trim = actual_category.trim_extents(policy)
        assert actual_trim == expected_trim, \
            f"Trim extents failed for '{read_name}': expected {expected_trim}, got {actual_trim}"
        
        status = "✓ PASS" if actual_category != ReadCategory.OTHER else "✓ PASS (NOW TRIMMED!)"
        print(f"  {read_name:<20} → {actual_category.name:<12} → trim {actual_trim} {status}")
    
    print("\n=== KEY IMPROVEMENT ===")
    print("Before: ReadCategory.OTHER reads got (0, 0) - NO TRIMMING")  
    print("After:  ReadCategory.OTHER reads get (30, 35) - PROPER TRIMMING")
    print("\nThis now supports:")
    print("- Illumina single-end reads")
    print("- Oxford Nanopore reads")  
    print("- PacBio reads")
    print("- Ion Torrent reads")
    print("- Any untagged sequencing data")
    
    print("\n✅ All tests passed! Single-end trimming is now working.")

if __name__ == "__main__":
    test_single_end_logic()