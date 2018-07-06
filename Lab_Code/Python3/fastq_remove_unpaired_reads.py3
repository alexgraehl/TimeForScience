#!/usr/bin/env python3
'''
This script requires python3. 

WARNING: it reads all the IDs frmo a fastq file into memory, so you need to have at least enough memory to store the read names!
If the files are huge, then you might need a lot of space! 12 GB was insufficient for two 33 GB (compressed) fastq files.

You can test it like so:

# should have all the same reads
   python3 `which fastq_remove_unpaired_reads.py3` --out1=a.fq.gz --out2=b.fq.gz --verbose ~/workspace/DATA_RESCOMP/fastq_test_data/A_human_pair*.5000.fq.gz

# should have NONE of the same reads
   python3 `which fastq_remove_unpaired_reads.py3` --out1=a.fq.gz --out2=b.fq.gz --verbose ~/workspace/DATA_RESCOMP/fastq_test_data/A_human_pair1.5000.fq.gz ~/workspace/DATA_RESCOMP/fastq_test_data/B_human_pair1.4000.fq.gz
'''
import os.path
import sys
import re
import bz2
import gzip
from collections import OrderedDict
from subprocess import call
# import ipdb; ipdb.set_trace()
import argparse

def open_compressed_agw(filename, mode, verbose=False):
    if verbose: sys.stderr.write("[fastq_remove_unpaired_reads.py3]: Parsing maybe-compressed file named: " + str(filename) + "\n")
    if re.search(r"[.](gz|gzip|Z)$", filename, flags=re.IGNORECASE):
        try:
            return gzip.open(filename=filename, mode=mode)
        except OSError as e:
            sys.stderr.write("Hey! Note from Alex here: this MAY BE A PROBLEM with python having issues with certain gzip format files. try 'gzip --test' on that file! See here: https://stackoverflow.com/questions/4928560/how-can-i-work-with-gzip-files-which-contain-extra-data")
            # you can consider running 'gzip --test' on that file to see if it has 'trailing garbage' data
            raise # re-raise
    if re.search(r"[.](bz2|bzip2)$", filename, flags=re.IGNORECASE):
        return bz2.BZ2File(filename=filename, mode=mode)
    else:
        return open(filename=filename, mode=mode)

def argErrorAndExit(msg="(No additional information given)"):
    raise SystemExit("[ERROR] in arguments to this script: " + msg)


def scrub_name(s):
    s2 = re.sub(r"\s.*", "", s.strip())
    s3 = re.sub(r"[/][12]", "/_ANY_", s2)
    return s3

def populate_OrderedDict_with_names(filename, verbose=False):
    # https://docs.python.org/3/library/collections.html#collections.OrderedDict
    sss = OrderedDict()
    if verbose: sys.stderr.write("[fastq_remove_unpaired_reads.py3]: Opening file named --> " + str(filename) + "\n")
    with open_compressed_agw(filename, 'rt', verbose=verbose) as fff:
        for linenum,line in enumerate(fff):
            if (linenum % 4 == 0):
                xid = scrub_name(line)
                if not xid.startswith("@"): raise("Invalid FASTQ line!")
                if xid in sss:
                    sys.stder.write("[Warning]: DUPLICATED RECORD NAME: the record name <" + str(xid) + "> appeared TWICE, with the second appearance on line number <" + str(linenum) + "> in file <" + filename + ">. Continuing anyway...\n")
                    pass
                sss[xid] = True
                pass
            pass
        pass
    return sss

def write_fastqs_in_set(infile, bothset, destname, verbose):
    assert isinstance(infile, str)
    assert isinstance(destname, str)
    assert isinstance(bothset, set)
    n_printed, n_omitted = 0, 0
    if (infile == destname): raise Exception("Infile cannot be the same as destination name.")
    with open_compressed_agw(infile, 'rt') as fff, open_compressed_agw(destname, 'wt') as dest:
        linenum = 0 # do NOT try to use enumerate here!
        for line in fff:
            if (linenum % 4 == 0):
                xid = scrub_name(line)
                if xid in bothset:
                    dest.write(line) # name
                    dest.write(fff.readline()) # sequence
                    dest.write(fff.readline()) # quality header
                    dest.write(fff.readline()) # quality score
                    linenum   += 3 # <-- remember that we just wrote THREE more lines!
                    n_printed += 1
                    #if verbose: sys.stderr.write("[OK]   Writing  ID " + xid + "\n")
                else:
                    if verbose: sys.stderr.write("[OMIT] OMITTING ID " + xid + "\n")
                    n_omitted += 1
                pass
            linenum += 1 # do NOT move this to an enumerate!
            pass
        pass
    sys.stderr.write("[fastq_remove_unpaired_reads.py3]: Wrote " + str(n_printed) + " records (and omitted " + str(n_omitted) + ") to --> " + destname + "\n")
    return (n_printed, n_omitted)


def main():
    PAIR_PLACEHOLDER_SYMBOL='@@@'
    parser = argparse.ArgumentParser(description="%(prog)s: fastq filterer in python 3. Removes any missing partial pairs, e.g. ReadA in forward pair file only, but ReadB in reverse pair file only---those would be removed.",
                                     epilog='''Example usage: python %(prog)s (examples go here)''',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--out1",  dest="out1", type=str, default=None, metavar="FILENAME", required=True, help="Output filename for fq1. Mutually exclusive with 'outpattern' above.")
    parser.add_argument("--out2",  dest="out2", type=str, default=None, metavar="FILENAME", required=True, help="Output filename for fq2. Mutually exclusive with 'outpattern' above.")
    parser.add_argument("--never-symlink", "--never_symlink", "--no-symlinks", dest="never_symlink", action="store_true", required=False, help="If specified, we do NOT just create a symlink when the files have identical reads. Default is to just symlink the files when all the reads at the same.")
    parser.add_argument("-v", "--verbose", dest="verbose",  action="store_true", help="Print verbose status messages to stderr")
    parser.add_argument("remainder", nargs=argparse.REMAINDER) # get the REMAINING un-parsed arguments (for example, a bunch of filenames)
    args = parser.parse_args()
    (out1,out2) = (args.out1, args.out2)
    if not re.search("[.](fq|fastq)[.](gz|gzip)$", out1, flags=re.I): argErrorAndExit("ARGUMENT ERROR: Your output filenames must both end in '.fastq.gz' or '.fq.gz'. They MUST be gzipped, we do not support uncompressed files! Offending not-gzipped-fastq name was: " + str(out2))
    if not re.search("[.](fq|fastq)[.](gz|gzip)$", out2, flags=re.I): argErrorAndExit("ARGUMENT ERROR: Your output filenames must both end in '.fastq.gz' or '.fq.gz'. They MUST be gzipped, we do not support uncompressed files! Offending not-gzipped-fastq name was: " + str(out2))
    if not 2 == len(args.remainder): argErrorAndExit("ARGUMENT ERROR: You must specify TWO input FASTQ files to filter out reads from.")
    (fq1, fq2) = (args.remainder[0], args.remainder[1])
    for f in [fq1, fq2]:
        if not os.path.isfile(f): raise Exception("Failed to find input file <"+f+">")
        pass
    dict1 = populate_OrderedDict_with_names(fq1, verbose=args.verbose)
    dict2 = populate_OrderedDict_with_names(fq2, verbose=args.verbose)
    if args.verbose: sys.stderr.write("Num records in file 1 (" + str(fq1) + "): " + str(len(dict1)) + "\n")
    if args.verbose: sys.stderr.write("Num records in file 2 (" + str(fq2) + "): " + str(len(dict2)) + "\n")
    files_have_all_matching_reads = (list(dict1.keys()) == list(dict2.keys()))
    if     files_have_all_matching_reads and args.verbose: sys.stderr.write("Files have matching read names, in the same order! Perfect.\n")
    if not files_have_all_matching_reads and args.verbose: sys.stderr.write("Files did NOT have all matching reads. Only printing matches.\n")
    if args.never_symlink:                                 sys.stderr.write("Since we are NOT supposed to symlink, we will be writing the files out again no matter what.\n")
    should_write_new_output_files = args.never_symlink or (not files_have_all_matching_reads)
    if should_write_new_output_files:
        keys_in_both = set(dict1.keys()).intersection(set(dict2.keys()))    
        (ok1, omit1) = write_fastqs_in_set(infile=fq1, bothset=keys_in_both, destname=out1, verbose=args.verbose)
        (ok2, omit2) = write_fastqs_in_set(infile=fq2, bothset=keys_in_both, destname=out2, verbose=args.verbose)
        assert ok1 == ok2
    else:
        call(["ln", "-s", os.path.realpath(fq1), out1]) # just symlink... don't re-write the file
        call(["ln", "-s", os.path.realpath(fq2), out2]) # just symlink... don't re-write the file
        sys.stderr.write("fastq_remove_unapired_reads has Generated two ABSOLUTE PATH symlinks: <" + out1 + "> and <" + out2 + ">\n")
        pass
    return # end of 'main'

# Must come at the VERY END!
if __name__ == "__main__":
    #assert sys.version_info >= (3,5), "You need to run this with python3.5 or later."
    main()
    pass

