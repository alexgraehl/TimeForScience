#!/usr/bin/python

'''
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
    global globalArgs    ## must have this here in order to ASSIGN globally!
    global globalOptions ## must have this here in order to ASSIGN globally!
    parser = optparse.OptionParser("usage: %prog [options]", version='%prog version 1.0')
    parser.add_option("-f", "--file", dest="filename", help="write report to FILE", metavar="FILE")
    parser.add_option("-N", "--name", dest="username", default="John Doe", type="string", help="specify a username to run as")
    parser.add_option("-p", "--port", dest="portnum", default=80, type="int", help="port number to run on")
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")
    (globalOptions, globalArgs) = parser.parse_args()

    #pdb.set_trace()
    for attr, value in globalOptions.__dict__.iteritems():
        print attr, value
        pass

    if len(globalArgs) != 2:
        parser.error("Args needs to be 2 for some reason: incorrect number of arguments")
        pass
    return

# Must come at the VERY END!
if __name__ == "__main__":
    print "Getting ready to handle command line arguments..."
    handleCommandLineOptions()
    print "Handled the command line arguments!"
    print("Here is the value for globalOptions.filename: " + globalOptions.filename)
    print("Note that that the global variable containing that attribute cannot be modified unless you set 'global globalOptions' in the code")
    lineNum = 0
    try:
        fff = open('some_file', 'r')
    except:
        raise

    for line in fff:
        lineNum += 1
        ldelim = line.split("\t")
        if (lineNum % 1000 == 0):
            print("Writing every 1000th line...")
            pass
        pass
    


    fff.close()



    pass





