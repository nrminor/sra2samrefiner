#!/usr/bin/env python3
"""
Comprehensive edge case analysis for trim_aligned_reads.py
Tests quality score handling, boundary conditions, and potential arithmetic errors.
"""


def analyze_quality_score_handling():
    """Analyze the quality score trimming logic for correctness."""

    print("=== QUALITY SCORE HANDLING ANALYSIS ===\n")

    # From trim_aligned_reads.py:288-293
    print("Current quality trimming logic:")
    print("  qlen = len(seq)")
    print("  cut_left = min(left, qlen)")
    print("  cut_right = min(right, max(qlen - cut_left, 0))")
    print("  keep_end = qlen - cut_left - cut_right")
    print("  new_qual = qual[cut_left : cut_left + keep_end]")
    print()

    # Test cases
    test_cases = [
        (10, 2, 1, "Normal case"),
        (10, 5, 8, "Over-trimming both ends"),
        (10, 15, 0, "Left trim > sequence length"),
        (10, 0, 12, "Right trim > sequence length"),
        (5, 2, 4, "Total trim > sequence length"),
        (0, 1, 1, "Empty sequence"),
        (1, 0, 0, "No trimming"),
    ]

    print("Testing boundary conditions:")
    for qlen, left, right, desc in test_cases:
        cut_left = min(left, qlen)
        cut_right = min(right, max(qlen - cut_left, 0))
        keep_end = qlen - cut_left - cut_right

        # Simulate slicing
        slice_start = cut_left
        slice_end = cut_left + keep_end

        print(f"{desc}: qlen={qlen}, left={left}, right={right}")
        print(f"  cut_left={cut_left}, cut_right={cut_right}, keep_end={keep_end}")
        print(
            f"  slice=[{slice_start}:{slice_end}], length={max(0, slice_end - slice_start)}"
        )

        # Check for issues
        if keep_end < 0:
            print(f"  ⚠️  NEGATIVE LENGTH: keep_end={keep_end}")
        elif slice_end < slice_start:
            print(f"  ⚠️  INVALID SLICE: end < start")
        elif slice_start < 0 or slice_end < 0:
            print(f"  ⚠️  NEGATIVE INDICES")
        else:
            print(f"  ✅ OK")
        print()


def analyze_cigar_edge_cases():
    """Test CIGAR manipulation for edge cases that could cause errors."""

    print("=== CIGAR EDGE CASE ANALYSIS ===\n")

    # Edge cases that could break CIGAR logic
    edge_cases = [
        ("Empty CIGAR", []),
        ("Only hard clips", [(5, 10), (5, 5)]),
        ("Only insertions", [(1, 20)]),
        ("Only deletions", [(2, 15)]),
        (
            "Mixed with zero-length ops",
            [(0, 10), (1, 0), (0, 5)],
        ),  # Invalid but possible
        ("All soft clips", [(4, 30)]),
        ("Complex pattern", [(4, 2), (0, 5), (1, 3), (0, 10), (2, 2), (0, 8), (4, 1)]),
    ]

    REF_CONSUME = {0, 2, 3, 7, 8}
    QRY_CONSUME = {0, 1, 4, 7, 8}

    print("Testing CIGAR patterns:")
    for desc, cigar in edge_cases:
        ref_consumed = sum(length for op, length in cigar if op in REF_CONSUME)
        qry_consumed = sum(length for op, length in cigar if op in QRY_CONSUME)

        print(f"{desc}: {cigar}")
        print(f"  Reference consumed: {ref_consumed}")
        print(f"  Query consumed: {qry_consumed}")

        # Check for issues
        if any(length <= 0 for op, length in cigar):
            print(f"  ⚠️  ZERO OR NEGATIVE LENGTH OPERATIONS")
        if not cigar:
            print(f"  ⚠️  EMPTY CIGAR (unmapped read?)")
        if qry_consumed == 0 and ref_consumed == 0:
            print(f"  ⚠️  NO CONSUMPTION (unmapped?)")
        print()


def analyze_compaction_logic():
    """Test CIGAR compaction logic for correctness."""

    print("=== CIGAR COMPACTION ANALYSIS ===\n")

    # Simulate push_compact logic from lines 158-169
    def push_compact(cigar_list, op, ln):
        """Simulate the push_compact method"""
        if ln <= 0:
            return cigar_list
        if cigar_list and cigar_list[-1][0] == op:
            last_op, last_len = cigar_list[-1]
            cigar_list[-1] = (op, last_len + ln)
        else:
            cigar_list.append((op, ln))
        return cigar_list

    # Test compaction scenarios
    test_scenarios = [
        ("Adjacent same ops", [(0, 5), (0, 3)], None),
        ("Different ops", [(0, 5), (1, 3)], None),
        ("Zero length add", [(0, 5)], (0, 0)),
        ("Negative length add", [(0, 5)], (0, -2)),
        ("Empty list add", [], (0, 5)),
    ]

    print("Testing CIGAR compaction:")
    for desc, initial, add_op in test_scenarios:
        cigar = initial.copy()
        if add_op:
            op, ln = add_op
            original = cigar.copy()
            result = push_compact(cigar, op, ln)
            print(f"{desc}: {original} + ({op},{ln}) = {result}")
        else:
            print(f"{desc}: {cigar}")
        print()


def check_reference_start_arithmetic():
    """Verify reference_start updates are arithmetically sound."""

    print("=== REFERENCE START ARITHMETIC ===\n")

    print("Key insight: reference_start must be updated when consuming")
    print("reference bases from the LEFT side of CIGAR only.")
    print()

    # Test the ref_advance calculation from _consume_from_left
    test_cigars = [
        [(0, 10), (1, 5), (0, 10)],  # M-I-M: only M ops advance reference
        [(4, 3), (0, 15), (4, 2)],  # S-M-S: only M advances reference
        [(1, 8), (2, 4), (0, 12)],  # I-D-M: D and M advance reference
        [(5, 10), (0, 20)],  # H-M: only M advances reference
    ]

    BOTH_CONSUME = {0, 7, 8}  # M, =, X
    QRY_CONSUME = {0, 1, 4, 7, 8}

    print("Testing reference advancement calculation:")
    for i, cigar in enumerate(test_cigars, 1):
        print(f"Test {i}: {cigar}")

        # Simulate consuming from left with different trim amounts
        for trim_q in [0, 3, 8, 15]:
            remaining = trim_q
            ref_advance = 0
            consumed_ops = []

            for op, length in cigar:
                if remaining <= 0:
                    break
                if op in QRY_CONSUME:
                    take = min(length, remaining)
                    remaining -= take
                    consumed_ops.append(f"{take}/{length} of op {op}")
                    if op in BOTH_CONSUME:
                        ref_advance += take

            print(
                f"  trim_q={trim_q}: ref_advance={ref_advance}, consumed={consumed_ops}"
            )
        print()


if __name__ == "__main__":
    analyze_quality_score_handling()
    analyze_cigar_edge_cases()
    analyze_compaction_logic()
    check_reference_start_arithmetic()
