#!/usr/bin/python

'''
Note: we use the 'optparse' module to parse command line arguments, becuase the 'argparse' requires Python 2.7 (which is not present in older distributions, like Ubuntu 10.x) http://docs.python.org/release/2.5.2/lib/module-optparse.html
A python program.
'''

from __future__ import print_function
from __future__ import division

import sys
import optparse
import pdb; #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques
import textwrap
#import time # we just want "sleep"
#import os.path

globalOptions = None
globalArgs    = None

def argErrorAndExit(msg="(No additional information given)"):
    raise SystemExit("[ERROR] in arguments to this script: " + msg)

def handleCommandLineOptions():
    global globalArgs    ## must have this here in order to ASSIGN globally!
    global globalOptions ## must have this here in order to ASSIGN globally!
    parser = optparse.OptionParser("usage: %prog [options]", version='%prog version 1.0')
    parser.add_option("-f", "--file", dest="filename", help="write report to FILE", metavar="FILE")
    parser.add_option("-N", "--name", dest="username", default="John Doe", type="string", help="specify a username to run as")
    parser.add_option("-p", "--port", dest="portnum", default=80, type="int", help="port number to run on")
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")
    (globalOptions, globalArgs) = parser.parse_args()

    #pdb.set_trace()
    print("There were this many named command line arguments (including DEFAULT values that were not specified by the user): " + str(len(globalOptions.__dict__)))
    for attr, value in globalOptions.__dict__.iteritems():
        print("   Argument <" + attr + "> = " + str(value))
        pass

    print("There were also this many un-parsed command line arguments: " + str(len(globalArgs)))
    for item in globalArgs:
        print("   Un-parsed command line argument: " + item)
        pass


    if len(globalArgs) == 1:
        print("Probaly you want the first un-parsed input on the command line to be a filename that the user specified. So this is just globalArgs[0]. You can use this in the main function below.")
        pass

    if len(globalArgs) != 1:
        parser.error("Args needs to be 1 for some reason: incorrect number of arguments. In other words, you didn't specify a filename (or some additional parameter) at the end of the command line invocation.")
        pass

    return

# Must come at the VERY END!
if __name__ == "__main__":
    print("Getting ready to handle command line arguments...")
    handleCommandLineOptions()
    print("Handled the command line arguments!")
    print("Here is the value for globalOptions.filename: " + str(globalOptions.filename))
    print("Note that that the global variable containing that attribute CANNOT BE MODIFIED unless you set 'global globalOptions' or 'global globalArgs' in the code below.")
    lineNum = 0
    try:
        with open('some_file', 'r') as fff:
            for line in fff:
                lineNum += 1
                ldelim = line.split("\t")
                if (lineNum % 100 == 0):
                    print("Writing every 100th line...")
                    pass
                pass # end 'for'
            pass # end 'with'
    except:
        print("Failed to open the example test file. Normally we should re-raise this exception, probably.")
        pass  #raise

    pass





