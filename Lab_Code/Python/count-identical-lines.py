#!/usr/bin/python


ALEX_PROGRAM_USAGE_TEXT='''
Duplicate line counter. For fasta / fastq / csfasta files.

Usage from a file: count-identical-lines.py YOURFILE.txt > out.txt

Usage from STDIN:  cat YOURFILE.txt | count-identical-lines.py > out.txt

You may have to use "cut" beforehand if you only want ONE column out of a file.

Note that this is ONLY useful for incredibly large files. Otherwise, just use UNIX sort, like so:
 sort YOURFILE | uniq -c

(It only takes about 5 times longer to use UNIX sort on a 2.2 GB file using "sort THEFILE | uniq -c")

Remember to redirect STDOUT to a file (that is, don't forget the '> out.txt' part of the command).
'''

import getopt
import sys
import os
import re
import random
import textwrap

#import pdb; #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques

TERMINAL_WIDTH = 80

def usageAndQuit(exitCode, message=None):
    message = textwrap.fill(message, TERMINAL_WIDTH)
    if (message is not None):
        print("")
        print("simulator: ")
        print(message)
        print("simulator: Printing usage information below.")
        print("*"*TERMINAL_WIDTH) ; print("*"*TERMINAL_WIDTH) ;
        pass
    print(ALEX_PROGRAM_USAGE_TEXT) # at the very bottom of this file
    if (message is not None):
        print("(End of usage information)")
        print("*"*TERMINAL_WIDTH) ; print("*"*TERMINAL_WIDTH) ;
        print("simulator: " + message)
        print("*"*TERMINAL_WIDTH) ; print("*"*TERMINAL_WIDTH) ;
        print("[Program Terminated]")
        pass
    sys.exit(exitCode)
    return


if __name__ == "__main__":
    sys.stderr.write("Note that you can get the same results, only sorted, with the following UNIX commands: sort YOURFILE | uniq -c > OUTPUT_FILE\nThat is about 5 times slower on a 2.2 GB file (15 minutes vs 3 minutes), but that isn't typically a huge deal.\n")
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "h", ["help"])
        # Docs for getopt: http://docs.python.org/library/getopt.html
    except (getopt.GetoptError):
        usageAndQuit(1, "Encountered an unknown command line option!\n")
        raise
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usageAndQuit(0, "Printed the HELP information, since the --help option was supplied.")
            pass
        print("Unprocessed arguments:" , args)
        print("Unprocessed options:" , opts)
        
        pass


    if (len(args) == 0 or args[0] == '-'): # hyphen means 'read from stdin'
        theFile = sys.stdin
    else:
        if (len(args) != 1):
            print(args)
            print(len(args))
            usageAndQuit(1, "ARGUMENTS ERROR: We need exactly ONE filename passed in on the command line.")
            raise

        inputFilename = args[0]
        try:
            theFile = open(inputFilename, 'r') ## The annotated bed file MUST have the "GTF" gene annotation so we know which reads are associated with which genes.
        except:
            sys.stderr.write("ERROR: Could not open specified file " + inputFilename)
            raise
        pass

    ddd = dict()
    lineNum = 0
    for line in theFile:
        lineNum += 1
        if (line in ddd):
            ddd[line] += 1
        else:
            ddd[line] = 1
            pass
        pass
    theFile.close()
    
    for key, value in ddd.iteritems():
        sys.stdout.write(str(value) + "\t" + key.rstrip() + "\n") # rstrip removes whitespace from right side of the key. This is important!
        pass

    sys.stderr.write("[Done -- Read a total of " + str(lineNum) + " lines.]")
    pass



##
