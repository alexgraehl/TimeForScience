#!/usr/bin/python

USAGE = '''%prog [options] <infile.txt>   >   table.tab

%prog takes a possibly-ragged-edged input file (differing number of items per line)
and outputs a totally rectangular tabular file (same number of items per line).

It's particularly useful if you're going to use <paste> to combine files, as each
input line must have the same number of lines for this to work reliably.

Input can be a single file, or blank to read from STDIN. Output is to STDOUT.

The special filename '-' (hyphen) will also force reading from STDIN.

Example 1:  tabular_from_ragged.py --delim="\t" --pad="NO_VALUE" myfile.txt > table.txt

Example 2:  zcat myfile.gz | tabular_from_ragged.py > table.txt
   * Reading from STDIN generates (and deletes) a temporary file in /tmp.
   * (This could be slow for huge files, or if /tmp is on another filesystem.)
'''

# Note: we use the 'optparse' module to parse command line arguments, becuase the superior 'argparse' requires Python 2.7 (which is not present in older distributions, like Ubuntu 10.x) http://docs.python.org/release/2.5.2/lib/module-optparse.html

import textwrap
import sys
import optparse
import os
import pdb; #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques

#import time # we just want "sleep"
#import os.path


globalOptions = None
globalArgs    = None

def handleCommandLineOptions():
    global globalArgs
    global globalOptions

    parser = optparse.OptionParser(USAGE, version='%prog version 1.0, generates a tabular-style file with the same number of elements in each line, when given a (possibly) ragged-edged file. You can change the input/output delimiter with --delim=<DELIM>. You can change what text gets used for padding out empty cells with --pad=<PADDING>.')
    parser.add_option("-d", "--delim", dest="delim", default="\t", type="string", help="Specify the delimiter character. Default is a <tab>.")
    parser.add_option("-p", "--pad", dest="pad", default="", type="string", help="Specify the text to print in a padded-out cell. Default is nothing. 'NA' is a popular option.")
    (globalOptions, globalArgs) = parser.parse_args()

    #pdb.set_trace()
    for attr, value in globalOptions.__dict__.iteritems():
        #print attr, value
        pass

    if len(globalArgs) < 1:
        sys.stderr.write("table-no-ragged.py STDERR message: Since no filenames were specified, we are reading from STDIN.\n")
        #parser.error("We need one un-parsed argument (a filename to read). This program is NOT currently able to handle reading from STDIN (standard in)!")
        pass

    if len(globalArgs) > 1:
        parser.error("You can only specify ONE filename on the command line! You specified more than one.")
        pass

    return

# Must come at the VERY END!
if __name__ == "__main__":
    handleCommandLineOptions()

    readingFromSTDIN = ((len(globalArgs) == 0) or (globalArgs[0] == '-')) ## "is a filename present on the command line? If not, then we're reading from STDIN"
    
    if readingFromSTDIN:
        TEMPFILENAME = "/tmp/table-no-ragged-temp-file-agw-table-random-file-delete-me-please.tmp"
        # Ok, NO filenames were specified on the command line, so we'll
        # write STDIN to a temporary file and then read that file.
        # This is because table-no-ragged requires TWO passes of the data, so we need a file (we can only read STDIN once).
        theFile = TEMPFILENAME
        ftemp = open(theFile, 'w')
        for line in sys.stdin:
            ftemp.write(line)
            pass
        ftemp.close()
    else:
        theFile = globalArgs[0] ## <-- the filename on the command line, if present
        pass

    maxNumElementsInLine = -1

    # Note to programmers: don't print anything here, or you'll break the fact that this outputs to stdout! You could print to stderr, though.
    
    lineNum = 0

    fff = open(theFile, 'r')
    for lineStr in fff:
        lineNum += 1
        l = lineStr.split(globalOptions.delim)
        maxNumElementsInLine = max(maxNumElementsInLine, len(l))
        pass
    fff.close()

    fff = open(theFile, 'r')
    for lineStr in fff:
        lineNum += 1
        l = lineStr.split(globalOptions.delim)
        numElementsToAddToThisLine = (maxNumElementsInLine - len(l))
        sys.stdout.write(lineStr.rstrip('\n') + ((globalOptions.delim + globalOptions.pad)*numElementsToAddToThisLine) + "\n")
        pass
    fff.close()

    # Try to remove the temp file, if we generated one.
    if readingFromSTDIN:
        if os.access(TEMPFILENAME, os.W_OK):
            os.remove(TEMPFILENAME)
        else:
            sys.stderr.write("WARNING: Could NOT remove the temp file we generated named " + TEMPFILENAME + " !\n")
            pass
        pass



    pass

