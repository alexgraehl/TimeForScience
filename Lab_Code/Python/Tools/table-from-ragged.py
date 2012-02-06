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





# #!/usr/bin/python

# '''
# A python program that takes an input tab-delimited file and makes it so it doesn't have ragged ends. Everything is padded out with BLANK entries so that all lines have the same number of tabs.

# Options:
# --delim (-d): The input and output delimiter. Default is a tab
# --fill (-x): The thing to fill ragged columns with. Default is nothing (''). Can be NA or whatever you want.

# Example:

# Sample Input:
#   this
#   is_a    file   with
#   ragged  edges

# Output of --fill="NA":
#   this    NA     NA
#   is_a    file   with
#   ragged  edges  NA

# Note: we use the 'optparse' module to parse command line arguments, becuase the superior 'argparse' requires Python 2.7 (which is not present in older distributions, like Ubuntu 10.x) http://docs.python.org/release/2.5.2/lib/module-optparse.html

# '''

# import sys
# import optparse
# import string

# import pdb #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques

# import textwrap
# #import time # we just want "sleep"
# #import os.path


# globalOptions = None
# globalArgs    = None

# def handleCommandLineOptions():
#     global globalArgs    ## must have this here in order to ASSIGN globally!
#     global globalOptions ## must have this here in order to ASSIGN globally!
#     parser = optparse.OptionParser("usage: %prog [-d <delim>] [-x <filler>] FILENAME", version='%prog version 1.0')
#     parser.add_option("-d", "--delim", dest="delim", default="\t", type="string", help="A delimiter between elements. Default is [tab]. Output and input delimiter are always the same.")
#     parser.add_option("-x", "--fill", dest="fill", default="", type="string", help="The text to fill blank cells with. Default is [nothing at all]. One common value for this is NA.")
#     (globalOptions, globalArgs) = parser.parse_args()

#     if len(globalArgs) != 1:
#         parser.error("You need to specify exactly one filename for this program to operate on.")
#         pass
#     return

# # Must come at the VERY END!
# if __name__ == "__main__":
#     handleCommandLineOptions()
#     fn = globalArgs[0] # filename
#     lineNum = 0

#     maxItemsPerLine = 0

#     ## read the file once to JUST count items per line
#     with open(fn, 'r') as fff:
#         for line in fff:
#             lineNum += 1
#             splitup = line.split( globalOptions.delim )
#             if (len(splitup) > maxItemsPerLine):
#                 maxItemsPerLine = len(splitup)
#                 pass
#             pass
#         pass

#     ## read the file a second time to actually print the output
#     with open(fn, 'r') as fff:
#         for line in fff:
#             splitup = line.split( globalOptions.delim )
#             numItemsToAdd = maxItemsPerLine - len(splitup) ## number of extra columns to pad here
            
#             sys.stdout.write( string.join( splitup, globalOptions.delim).rstrip("\n\r") ) ## print the original string
#             sys.stdout.write(numItemsToAdd * (globalOptions.delim + globalOptions.fill)) ## <-- print the right number of delimiters!
#             sys.stdout.write("\n")
#             pass
#         pass

#     pass





