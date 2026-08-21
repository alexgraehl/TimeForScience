#!/usr/bin/env python3
import fastq_remove_unpaired_reads as test

# ---------------------------------------------------------------------------
# scrub_name
def test_scrub_name_removes_trailing_comments():
    assert test.scrub_name("@READ_ALPHA/1 junk comment here") == "@READ_ALPHA/_ANY_"

def test_scrub_name_resolves_matching_records_to_same_name():
    assert test.scrub_name("@READ_BETA/1 some comment") == test.scrub_name("@READ_BETA/2 unrelated comment")

def test_scrub_name_unpaired_is_unchanged_but_comment_is_still_removed():
    # If there's no '/1' or '/2', we don't add '/_ANY_'
    assert test.scrub_name("@READ_GAMMA cool comment") == "@READ_GAMMA"  # The comment is removed, but note the lack of "ANY" here.
# ---------------------------------------------------------------------------

