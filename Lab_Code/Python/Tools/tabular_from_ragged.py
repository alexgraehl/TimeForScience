#!/usr/bin/python

'''
This is a program that takes a possibly-ragged-edged input file and outputs a totally rectangular tabular file.
It's useful when programs output differing numbers of items for different lines, but you want to operate as if all the lines had the same number of elements.

Note: input MUST be a file, it cannot be STDIN (i.e., you cannot directly send any UNIX pipes to this program). Output is STDOUT.

Usage:   tabular_from_ragged.py  --delim="\t"  --pad="NO_VALUE"   YOUR_FILE.txt

Note: we use the 'optparse' module to parse command line arguments, becuase the superior 'argparse' requires Python 2.7 (which is not present in older distributions, like Ubuntu 10.x) http://docs.python.org/release/2.5.2/lib/module-optparse.html

A python program.
'''

import sys
import optparse

import pdb; #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques

import textwrap
#import time # we just want "sleep"
#import os.path


globalOptions = None
globalArgs    = None

def handleCommandLineOptions():
    global globalArgs
    global globalOptions

    parser = optparse.OptionParser("usage: %prog [options]", version='%prog version 1.0, generates a tabular-style file with the same number of elements in each line, when given a (possibly) ragged-edged file. You can change the input/output delimiter with --delim=<DELIM>. You can change what text gets used for padding out empty cells with --pad=<PADDING>.')
    parser.add_option("-d", "--delim", dest="delim", default="\t", type="string", help="Specify the delimiter character. Default is a <tab>.")
    parser.add_option("-p", "--pad", dest="pad", default="", type="string", help="Specify the text to print in a padded-out cell. Default is nothing. 'NA' is a popular option.")
    (globalOptions, globalArgs) = parser.parse_args()

    #pdb.set_trace()
    for attr, value in globalOptions.__dict__.iteritems():
        #print attr, value
        pass

    if len(globalArgs) != 1:
        parser.error("We need one un-parsed argument (a filename to read). This program is NOT currently able to handle reading from STDIN (standard in)!")
        pass
    return

# Must come at the VERY END!
if __name__ == "__main__":
    handleCommandLineOptions()

    theFile = globalArgs[0]

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

    pass





