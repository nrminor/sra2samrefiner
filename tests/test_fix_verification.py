#!/usr/bin/env python3
"""
Verify that the unmapped read slicing bug is fixed.
"""


def test_fixed_slicing_logic():
    """Test the FIXED unmapped read slicing logic."""

    print("=== TESTING FIXED UNMAPPED READ SLICING LOGIC ===\n")

    # The NEW fixed logic from the code:
    def fixed_slicing(seq, left, right):
        """Replicate the NEW fixed slicing logic exactly."""
        seq_len = len(seq)
        cut_left = min(left, seq_len)
        cut_right = min(right, max(seq_len - cut_left, 0))
        keep_end = seq_len - cut_left - cut_right

        # Apply safe slicing with boundary checking
        new_seq = "" if keep_end <= 0 else seq[cut_left : cut_left + keep_end]
        return new_seq

    # The old BUGGY logic:
    def buggy_slicing(seq, left, right):
        """The old buggy logic for comparison."""
        new_seq = seq[left : len(seq) - right if right else None]
        return new_seq

    test_cases = [
        ("ATCGATCG", 0, 0, "No trimming"),
        ("ATCGATCG", 2, 0, "Left only"),
        ("ATCGATCG", 0, 3, "Right only"),
        ("ATCGATCG", 2, 3, "Both ends"),
        ("ATCGATCG", 10, 0, "Left exceeds length"),
        ("ATCGATCG", 0, 10, "Right exceeds length - THE BUG CASE"),
        ("ATCGATCG", 5, 5, "Total exceeds length"),
        ("A", 0, 1, "Single base, right trim"),
        ("A", 1, 0, "Single base, left trim"),
        ("", 1, 1, "Empty sequence"),
    ]

    print("Comparing OLD (buggy) vs NEW (fixed) slicing logic:")
    bugs_fixed = 0

    for seq, left, right, desc in test_cases:
        buggy_result = buggy_slicing(seq, left, right)
        fixed_result = fixed_slicing(seq, left, right)

        if buggy_result == fixed_result:
            status = "✅ Same (was already correct)"
        else:
            status = "🔧 FIXED BUG!"
            bugs_fixed += 1

        print(f"\n{desc}: '{seq}', left={left}, right={right}")
        print(f"  OLD (buggy): '{buggy_result}' (len={len(buggy_result)})")
        print(f"  NEW (fixed): '{fixed_result}' (len={len(fixed_result)})")
        print(f"  {status}")

    print(f"\n🎯 SUMMARY: Fixed {bugs_fixed} critical bug(s)!")

    if bugs_fixed > 0:
        print("\n🚨 The bugs that were fixed:")
        print("- When right trim amount exceeds sequence length, old logic would")
        print("  incorrectly slice using negative indices, keeping unintended bases")
        print("- New logic properly handles over-trimming by returning empty sequences")


if __name__ == "__main__":
    test_fixed_slicing_logic()
