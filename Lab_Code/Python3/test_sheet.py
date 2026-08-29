#!/usr/bin/env python3

# FYI: Run this with either "python -m pytest" or just "make test" if you prefer.

# Only covers the helper functions. I need to add something to actually test
# the terminal output if this is going to be legitimately useful.

import pytest
import sheet as test


# =================================================================
# Test truncateLongCell() and padStrToLength()
# -----------------------------------------------------------------
def test_truncateLongCell_works():
    assert test.truncateLongCell("1234", 4, "...") == "1234"
    assert test.truncateLongCell("12345", 4, "...") == "1..."
    assert test.truncateLongCell("12345678", 8, "...")  == "12345678"
    assert test.truncateLongCell("123456789", 8, "...") == "12345..."

def test_truncateLongCell_paradoxically_makes_longer_strings_in_pathological_cases():
    assert test.truncateLongCell("1234", 1, "super_long_suffix") == "super_long_suffix"  # Maybe we should rework this
    assert test.truncateLongCell("1234", -1, "super_long_suffix") == "super_long_suffix"  # Maybe we should rework this

def test_padStrToLength_works():
    assert test.padStrToLength("ab", 5, "~") == "ab~~~"
    assert test.padStrToLength("abcde", 5, "~") == "abcde"
    assert test.padStrToLength("abcdefg", 5, "~") == "abcdefg"  # (Doesn't trim)
