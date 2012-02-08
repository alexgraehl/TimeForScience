#!/usr/bin/python

'''
Note: we use the 'optparse' module to parse command line arguments, becuase the superior 'argparse' requires Python 2.7 (which is not present in older distributions, like Ubuntu 10.x) http://docs.python.org/release/2.5.2/lib/module-optparse.html

A python program.
'''

from subprocess import call
import sys
import optparse
import os
import re

import pdb #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques

import textwrap
#import time # we just want "sleep"
#import os.path


globalOptions = None
globalArgs    = None

def handleCommandLineOptions():
    global globalArgs    ## must have this here in order to ASSIGN globally!
    global globalOptions ## must have this here in order to ASSIGN globally!
    parser = optparse.OptionParser("usage: %prog [options]", version='%prog version 1.0')
#    parser.add_option("-f", "--file", dest="filename", help="write report to FILE", metavar="FILE")
#    parser.add_option("-N", "--name", dest="username", default="John Doe", type="string", help="specify a username to run as")
#    parser.add_option("-p", "--port", dest="portnum", default=80, type="int", help="port number to run on")
#    parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")
    (globalOptions, globalArgs) = parser.parse_args()

    if len(globalArgs) < 1:
        print("You need to add at least one directory or filename to check for proper backing up!")
        pass

    for i in range(len(globalArgs)):
        if (not re.search("^/", globalArgs[i])):
            parser.error("All paths to check must be COMPLETE full paths, starting from the root level (i.e., starting with a / ). You specified at least one argument that was not a full path from the root level, namely: " + globalArgs[i])
            pass

        if (not re.search("/$", globalArgs[i])): ## If there's no trailing slash...
            globalArgs[i] += "/" ## ...then add a trailing slash no matter what
            pass

        pass


    return


def writeSizesToFile(inWhere, toFile):
    tempName = "check-backup-times-1.tmp"
    acmd = "find " + inWhere + " -type d -name '[a-zA-Z0-9]*' -print0 | xargs -0 du > " + tempName

    call(acmd, shell=True) ## find all the files and calculate their sizes

    regexpReadyName = re.sub("/", "[/]", inWhere) ## escape slashes so that SED can use them below
    regexpReadyName = re.sub(" ", "[ ]", regexpReadyName)    ## escape literal spaces so that SED can use them below
    replacePathsCmd = " sed 's/" + regexpReadyName + "//' " + tempName + " > " + toFile

    print(replacePathsCmd)
    call(replacePathsCmd, shell=True)
    return


# Must come at the VERY END!
if __name__ == "__main__":
    handleCommandLineOptions()
    #print("Note that that the global variable containing that attribute CANNOT BE MODIFIED unless you set 'global globalOptions' or 'global globalArgs' in the code below.")


    for i in range(len(globalArgs) - 1):
        if (i % 2 != 0):
            continue  ## <-- skip the odd input arguments---note that we compare these TWO at a time!
        

        fileA = "check-backup-times-z1.tmp"
        fileB = "check-backup-times-z2.tmp"

        aaa = globalArgs[i]
        bbb = globalArgs[i+1]

        print(aaa)
        writeSizesToFile(aaa, fileA)

        print(bbb)
        writeSizesToFile(bbb, fileB)


        diffCmd = "diff " + fileA + " " + fileB + " > " + " diff.out"
        call(diffCmd , shell=True)

        pass


    pass





