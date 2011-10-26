#!/usr/bin/python

'''
Note: argparse requires Python 2.7. Hence, we use optparse here so that our scripts work with older Python installations (like Ubuntu 10.x)

http://docs.python.org/release/2.5.2/lib/module-optparse.html

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
    pass





