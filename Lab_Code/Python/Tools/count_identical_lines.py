#!/usr/bin/python


ALEX_PROGRAM_USAGE_TEXT='''
Duplicate counter. For fasta / fastq / csfasta files.

Remember to output to STDOUT!
'''

import getopt
import sys
import os
import re
import random
import textwrap

import pdb; #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques

log = None
LOG_FILE_NAME = "sim.log.txt"
TERMINAL_WIDTH = 80

def logWrite(message):
    global log
    if (log is None):
        try:
            log = open(LOG_FILE_NAME, 'w')
        except:
            print("ERROR (Possibly due to lack of permissions in this directory?): Could not open the following output file required for logging status: " + LOG_FILE_NAME)
            raise
        pass
    #message = textwrap.fill(message, TERMINAL_WIDTH)
    log.write(message + "\n") ## actually write to the log FILE
    sys.stdout.write(message + "\n") ## ...and print to the screen
    return

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
    random.seed()

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "hwi:d"
                                       , ["help", "warn"
                                          ]
                                       )
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

    if (len(args) == 0):
        theFile = sys.stdin
    else:
        if (len(args) != 1):
            print(args)
            print(len(args))
            usageAndQuit(1, "ARGUMENTS ERROR: We need exactly ONE filename passed in on the command line.")
            raise

        inputFilename = args[0] #"b5_test_rand.bed"
        try:
            theFile = open(inputFilename, 'r') ## The annotated bed file MUST have the "GTF" gene annotation so we know which reads are associated with which genes.
        except:
            logWrite("ERROR (Possibly due to lack of permissions in this directory?): Could not open the following output file required for logging status: " + LOG_FILE_NAME)
            raise
        pass

    theHash = dict()

    lineNum = 0
    for line in theFile:
        '''Read through each gene in the annotated BED file.'''
        lineNum += 1
        
        if (line in theHash):
            theHash[line] += 1
        else:
            theHash[line] = 1
            pass

        pass
    theFile.close()
    
    for key, value in theHash.iteritems():
        sys.stdout.write("" + key.rstrip() + "\t" + str(value) + "\n")
        pass

    logWrite("[Done -- Successful exit from simulator. Read a total of " + str(lineNum) + " lines.]")
    log.close() ## Finally, close the diagnostic log file that we've been writing messages to this whole time
    pass



##
