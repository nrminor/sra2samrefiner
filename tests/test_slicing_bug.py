#!/usr/bin/env python3
"""
Test the sequence slicing logic for unmapped reads to find the bug.
"""

def test_unmapped_slicing_logic():
    """Test the unmapped read slicing logic for correctness."""
    
    print("=== TESTING UNMAPPED READ SLICING LOGIC ===\n")
    
    # The current logic from the code:
    def current_slicing(seq, left, right):
        """Replicate the current slicing logic exactly."""
        new_seq = seq[left : len(seq) - right if right else None]
        return new_seq
    
    # What the logic SHOULD be:
    def correct_slicing(seq, left, right):
        """Correct slicing logic."""
        seq_len = len(seq)
        # Use the same logic as the main trimming path
        cut_left = min(left, seq_len)
        cut_right = min(right, max(seq_len - cut_left, 0))
        keep_end = seq_len - cut_left - cut_right
        new_seq = "" if keep_end <= 0 else seq[cut_left : cut_left + keep_end]
        return new_seq
    
    test_cases = [
        ("ATCGATCG", 0, 0, "No trimming"),
        ("ATCGATCG", 2, 0, "Left only"),
        ("ATCGATCG", 0, 3, "Right only"), 
        ("ATCGATCG", 2, 3, "Both ends"),
        ("ATCGATCG", 10, 0, "Left exceeds length"),
        ("ATCGATCG", 0, 10, "Right exceeds length"),
        ("ATCGATCG", 5, 5, "Total exceeds length"),
        ("A", 0, 1, "Single base, right trim"),
        ("A", 1, 0, "Single base, left trim"),
        ("", 1, 1, "Empty sequence"),
    ]
    
    print("Testing sequence slicing logic:")
    for seq, left, right, desc in test_cases:
        current_result = current_slicing(seq, left, right)
        correct_result = correct_slicing(seq, left, right)
        
        match = "✅" if current_result == correct_result else "❌ BUG!"
        
        print(f"\n{desc}: '{seq}', left={left}, right={right}")
        print(f"  Current: '{current_result}' (len={len(current_result)})")
        print(f"  Correct: '{correct_result}' (len={len(correct_result)})")
        print(f"  {match}")
        
        if current_result != correct_result:
            print(f"  ⚠️  DIFFERENCE DETECTED!")

def test_edge_case_deep_dive():
    """Deep dive into the specific edge cases."""
    
    print("\n=== DEEP DIVE INTO EDGE CASES ===\n")
    
    # Test the problematic expression: len(seq) - right if right else None
    test_cases = [
        (8, 0, "right=0 case"),
        (8, 3, "right=3 case"), 
        (8, 8, "right=length case"),
        (8, 10, "right>length case"),
        (0, 0, "empty sequence"),
        (1, 1, "right=length on single base"),
    ]
    
    for seq_len, right, desc in test_cases:
        # The problematic expression
        slice_end = seq_len - right if right else None
        print(f"{desc}: len={seq_len}, right={right}")
        print(f"  slice_end = len(seq) - right if right else None = {slice_end}")
        
        # What happens with different left values
        for left in [0, 2, seq_len, seq_len + 1]:
            if left <= seq_len:
                simulated_seq = "A" * seq_len
                try:
                    result = simulated_seq[left:slice_end]
                    print(f"    seq[{left}:{slice_end}] = '{result}' (len={len(result)})")
                except Exception as e:
                    print(f"    seq[{left}:{slice_end}] = ERROR: {e}")
        print()

if __name__ == "__main__":
    test_unmapped_slicing_logic()
    test_edge_case_deep_dive()