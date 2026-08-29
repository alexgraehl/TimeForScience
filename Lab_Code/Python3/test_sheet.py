#!/usr/bin/env python3

# FYI: Run this with either "python -m pytest" or just "make test" if you prefer.

# Only covers the helper functions. I need to add something to actually test
# the terminal output if this is going to be legitimately useful.

import sheet as test

# =================================================================
# Test padStrToLength()
# -----------------------------------------------------------------
def test_padStrToLength_works():
    assert test.padStrToLength("ab", 5, "~") == "ab~~~"
    assert test.padStrToLength("abcde", 5, "~") == "abcde"
    assert test.padStrToLength("abcdefg", 5, "~") == "abcdefg"  # (Doesn't trim)


# =================================================================
# Test AGW_File_Data Regex handling
# -----------------------------------------------------------------
def test_AGW_File_Data_initial_regex_is_off():
    f = test.AGW_File_Data("abc.txt")
    assert not f.regexIsActive()
    assert f.getRegexString() == ""

def test_AGW_File_Data_search_enabled_then_matches():
    f = test.AGW_File_Data("file.txt")
    f.changeCurrentSearchTerm("abc")
    assert bool(f.regexIsActive())
    assert f.getRegexString() == "abc"
    assert f.stringDoesMatchRegex("qwerty_abcdefg")
    assert f.stringDoesMatchRegex("qwerty_ABCdefg")  # Doesn't match (case sensitive matches only)
    assert not f.stringDoesMatchRegex("-qwerty_zzzdefg")

def test_AGW_File_Data_default_search_is_case_insensitive():
    f = test.AGW_File_Data("file.txt")
    f.changeCurrentSearchTerm("qwerty")  # Default is case INSENSITIVE
    assert f.stringDoesMatchRegex("QWErty")
    f.changeCurrentSearchTerm("qwerty", None)  # 'None' doesn't change anything
    assert f.stringDoesMatchRegex("QWErty")

def test_AGW_File_Data_case_sensitivity_works():
    f = test.AGW_File_Data("file.txt")
    f.changeCurrentSearchTerm("zzz", None)  # Default is case INSENSITIVE
    assert f.stringDoesMatchRegex("zzz")
    assert f.stringDoesMatchRegex("ZZZ")  # Matches the differing-case letters too
    f.changeCurrentSearchTerm("zzz", True)  # Now it's case SENSITIVE
    assert f.stringDoesMatchRegex("zzz") == True # Now it only matches 'zzz'
    assert f.stringDoesMatchRegex("ZZZ") == False
    f.changeCurrentSearchTerm("zzz", None)  # No-op since "None" doesn't change the case sensitivity (still case-insensitive)
    assert f.stringDoesMatchRegex("ZZZ") == False  # Still ONLY matches the lower-case letters
    f.changeCurrentSearchTerm("zzz", False)  # Switch back to case-insensitive
    assert f.stringDoesMatchRegex("zzz")
    assert f.stringDoesMatchRegex("ZZZ")

def test_AGW_File_Data_no_match_when_no_search_term_set():
    f = test.AGW_File_Data("a.txt")
    assert not f.stringDoesMatchRegex("")

def test_AGW_File_Data_appendToCurrentSearchTerm_works():
    f = test.AGW_File_Data("f.txt")
    f.appendToCurrentSearchTerm("1")
    f.appendToCurrentSearchTerm("2")
    f.appendToCurrentSearchTerm("3")
    assert f.getRegexString() == "123"

def test_AGW_File_Data_trimRegex_removes_trailing_chars():
    f = test.AGW_File_Data("c.txt")
    f.changeCurrentSearchTerm("abc")
    f.trimRegex(1)
    assert f.getRegexString() == "ab"

def test_AGW_File_Data_search_becomes_inactive_if_trimmed_to_zero_length():
    f = test.AGW_File_Data("somefile.txt")
    f.changeCurrentSearchTerm("a")
    f.trimRegex(1)
    assert f.regexIsActive() == False
    f.trimRegex(99)  # Not really valid, but should still be tolerated
    assert f.regexIsActive() == False
    assert not f.getRegexString()

def test_AGW_File_Data_clearCurrentSearchTerm_works():
    f = test.AGW_File_Data("somefile.txt")
    f.changeCurrentSearchTerm("searchtime")
    assert f._compiledRegex is not None
    f.clearCurrentSearchTerm()
    assert not f.regexIsActive()
    assert not f.getRegexString()
    assert f._compiledRegex is None  # Should have been cleared out as well

def test_AGW_File_Data_search_becomes_inactive_if_empty_string_is_set():
    f = test.AGW_File_Data("somefile.txt")
    f.changeCurrentSearchTerm("")
    assert not f.regexIsActive()
    assert not f.getRegexString()
