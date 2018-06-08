#!/usr/bin/env python3
'''
This script requires python3. 
You can test it like so:  fastq_remove_unpaired_reads.py3 -o 'out_pair@@@.fq.gz' a_1.fastq.gz a_2.fastq.gz
# Resulting files from that invocation above would be 'out_pair1.fq.gz' and 'out_pair2.fq.gz'
Note that the '@@@' is special and is replaced by the pair number

'''
import os.path
import sys
import re
import gzip
import bz2
#sys.path.insert(0, os.path.join(os.environ['TIME_FOR_SCIENCE_DIR'], "Lab_Code","Python3"))

# import ipdb; ipdb.set_trace()
#from .open_compressed_agw import Open_compressed_agw
import argparse
#import sys
#import pdb   #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques
#import time # we usually just want the "sleep" function

def open_compressed_agw(filename, mode):
    if re.search(r"[.](gz|gzip|Z)$", filename, flags=re.IGNORECASE):
        return gzip.open(filename=filename, mode=mode)
    if re.search(r"[.](bz2|bzip2)$", filename, flags=re.IGNORECASE):
        return bz2.BZ2File(filename=filename, mode=mode)
    else:
        return open(filename=filename, mode=mode)

def argErrorAndExit(msg="(No additional information given)"):
    raise SystemExit("[ERROR] in arguments to this script: " + msg)

def main():
    parser = argparse.ArgumentParser(description="%(prog)s: fastq filterer in python 3.",
                                     epilog='''Example usage: python %(prog)s (examples go here)''',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--outpattern", "-o",  dest="outpattern", type=str, default=None, metavar="PATTERN", required=True, help="Output pattern. Should have a '#' somewhere, which will become a 1 or 2 in the output.")
    #parser.add_argument("-N", "--name",  dest="username", type=str,             default="Dave",               required=False, help="specify a username to run as")
    #parser.add_argument("-p", "--port",  dest="portnum",  type=int,             default=80,                                   help="port number to run on")
    parser.add_argument("-v", "--verbose", dest="verbose",  action="store_true", help="Print verbose status messages to stderr")
    parser.add_argument("remainder", nargs=argparse.REMAINDER) # get the REMAINING un-parsed arguments (for example, a bunch of filenames)
    args = parser.parse_args()

    PAIR_PLACEHOLDER_SYMBOL='@@@'
    if args.outpattern is None or not re.search(PAIR_PLACEHOLDER_SYMBOL, args.outpattern):
        print("ARGUMENT ERROR: You must specify an OUTPUT PATTERN with a " + PAIR_PLACEHOLDER_SYMBOL + " in it, like this:   -o MYFILE_Pair_Number_" + PAIR_PLACEHOLDER_SYMBOL + ".fastq.gz")
        print("                The " + PAIR_PLACEHOLDER_SYMBOL + " will be replaced with a '1' or '2'.")
        print("                The pattern you specified was: " + str(args.outpattern))
        sys.exit(1)

    if not re.search("[.](fq|fastq)[.](gz|gzip)$", args.outpattern, flags=re.IGNORECASE):
        print("ARGUMENT ERROR: Your output pattern must end in '.fastq.gz' or '.fq.gz'.")
        sys.exit(1)

    if not 2 == len(args.remainder):
        print("ARGUMENT ERROR: You must specify TWO input FASTQ files to filter out reads from.")
        sys.exit(1)

    (fq1, fq2) = (args.remainder[0], args.remainder[1])
    if not os.path.isfile(fq1): raise Exception("Failed to find input file <"+fq1+">")
    if not os.path.isfile(fq2): raise Exception("Failed to find input file <"+fq2+">")

    out1 = args.outpattern.replace(PAIR_PLACEHOLDER_SYMBOL, "1") # e.g. "_pair1"
    out2 = args.outpattern.replace(PAIR_PLACEHOLDER_SYMBOL, "2") # e.g. "_pair2"
    
    def scrub_name(s):
        s2 = re.sub(r"\s.*", "", s.strip())
        s3 = re.sub(r"[/][12]", "/_ANY_", s2)
        return s3

    def populate_set_with_names(filename):
        sss = set()
        with open_compressed_agw(filename, 'rt') as fff:
            for linenum,line in enumerate(fff):
                if (linenum % 4 == 0):
                    xid = scrub_name(line)
                    if not xid.startswith("@"): raise("Invalid FASTQ line!")
                    sss.add(xid)
                    #print("\t".join([xid, yid, str(xid == yid)]))
                    pass
                pass
            pass
        return sss

    in1, in2 = populate_set_with_names(fq1), populate_set_with_names(fq2)
    inboth = in1.intersection(in2)

    if args.verbose: sys.stderr.write("Num records in file 1: " + str(len(in1)) + "\n")
    if args.verbose: sys.stderr.write("Num records in file 2: " + str(len(in2)) + "\n")
    if args.verbose: sys.stderr.write("Num records in the INTERSECTION (in both sets): " + str(len(inboth)) + "\n")
    def write_fastqs_in_set(filename, sss, destname):
        n_printed, n_omitted = 0, 0
        if (filename == destname): raise Exception("Filename cannot be the same as destination name.")
        with open_compressed_agw(filename, 'rt') as fff, open_compressed_agw(destname, 'wt') as dest:
            linenum = 0 # do NOT try to use enumerate here!
            for line in fff:
                if (linenum % 4 == 0):
                    xid = scrub_name(line)
                    if xid in sss:
                        dest.write(line) # name
                        dest.write(fff.readline()) # sequence
                        dest.write(fff.readline()) # quality header
                        dest.write(fff.readline()) # quality score
                        linenum   += 3 # <-- remember that we just wrote THREE more lines!
                        n_printed += 1
                        if args.verbose: sys.stderr.write("[OK]   Writing  ID " + xid + "\n")
                    else:
                        if args.verbose: sys.stderr.write("[OMIT] OMITTING ID " + xid + "\n")
                        n_omitted += 1
                    pass
                linenum += 1 # do NOT move this to an enumerate!
                pass
            pass
        sys.stderr.write("[fastq_remove_unpaired_reads.py3]: Wrote " + str(n_printed) + " records (and omitted " + str(n_omitted) + ") to --> " + destname + "\n")
        return (n_printed, n_omitted)

    (ok1, omit1) = write_fastqs_in_set(fq1, inboth, destname=out1)
    (ok2, omit2) = write_fastqs_in_set(fq2, inboth, destname=out2)

    if ok1 != ok2: raise Exception("Should never occur by definition")
    #if ok1 != ok2:
    #    sys.stderr.print("WARNING: Differing numbers of output lines in inputs. FQ file 1 ({fq1}) had {ok1} printed lines and {omit1} omitted lines, while
    
    return # end of 'main'

# Must come at the VERY END!
if __name__ == "__main__":
    main()
    pass





