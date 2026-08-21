#!/usr/bin/env python3
import argparse
import pytest
import fastq_remove_unpaired_reads as test

# =================================================================
# Test validate_args()
# -----------------------------------------------------------------
def test_validation_ok():
    # Normal acceptable input
    good_args = argparse.Namespace(
            out1="out1.fq.gz",
            out2="out2.fq.gz",
            remainder=["input_fastq_1.fastq.gz", "input_fastq_2.fastq.gz"],
            never_symlink=True, verbose=True)
    assert test.validate_args(good_args)

def test_validation_fails_if_input_fastqs_have_same_filename_part():
    # Even if the output paths are different, the FILENAMES themselves must be distinct.
    fake_args = argparse.Namespace(
            out1="a_different_dir/same_basename.fq.gz",
            out2="other_directory/same_basename.fq.gz",
            remainder=["input_fastq_1.fastq.gz", "input_fastq_2.fastq.gz"])
    with pytest.raises(SystemExit):
        test.validate_args(fake_args)  # Blows up since both files are named same.fastq.gz, even though the FULL PATHS are different.


# =================================================================
# Test scrub_name()
# -----------------------------------------------------------------
def test_scrub_name_removes_trailing_comments():
    assert test.scrub_name("@READ_ALPHA/1 junk comment here") == "@READ_ALPHA/_ANY_"

def test_scrub_name_resolves_matching_records_to_same_name():
    assert test.scrub_name("@READ_BETA/1 some comment") == test.scrub_name("@READ_BETA/2 unrelated comment")

def test_scrub_name_unpaired_is_unchanged_but_comment_is_still_removed():
    # If there's no '/1' or '/2', we don't add '/_ANY_'
    assert test.scrub_name("@READ_GAMMA cool comment") == "@READ_GAMMA"  # The comment is removed, but note the lack of "ANY" here.
# =================================================================

