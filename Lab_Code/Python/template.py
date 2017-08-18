#!/usr/bin/env python
'''
This script should work on both python2
                            AND python3 directly,
    as long as you have the 'future' module installed! (e.g.: "pip install future")
See here for details on making a python 2/3 script: http://python-future.org/compatible_idioms.html

Requires at least python 2.7+ due to 'argparse' (which is not in python 2.6).

As a programmer, the main thing to watch out for is that you must use
            "future.utils.iteritems( mydictionary )" to iterate over a dictionary
            instead of "mydictionary.iteritems()"
            (You can make this shorter with 'from future.utils import iteritems as YOUR_NAME_HERE')

Almost everything else "just works."

Try it out:
        python2 ./template.py some_filename x y z
        python3 ./template.py some_filename x y z

'''
# =================== These are part of the modules that handle python 2/3 compatibility ==================
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
import future          # pip install future
import future.utils    # (installed with above)    #from future.utils import iteritems
import future.builtins # pip install future
from   future.builtins import range # example: mylist = list(range(5))     assert mylist == [0, 1, 2, 3, 4]
import past     # pip install future
#import six           # pip install six  Useful for: if isinstance(value, six.string_types):
try:   basestring
except NameError:   basestring = str  # useful for 'isinstance(value, basestring)'
# ================ Below are the NORMAL modules that you'd want for regular programming ====================

# import ipdb; ipdb.set_trace()
import argparse
import os
import sys
import pdb   #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques
import textwrap
#import time # we usually just want the "sleep" function

def argErrorAndExit(msg="(No additional information given)"):
    raise SystemExit("[ERROR] in arguments to this script: " + msg)

def main():
    print("Getting ready to handle command line arguments...")
    parser = argparse.ArgumentParser(description="%(prog)s: template for python 2 and python 3.",
                                     epilog='''Example usage: python %(prog)s (no example yet)\n(some examples go here)\n(More examples go here)''',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--file",  dest="filename", type=str,             default=None, metavar="FILE", required=False, help="Minimum x (length) value.")
    parser.add_argument("-N", "--name",  dest="username", type=str,             default="Dave",               required=False, help="specify a username to run as")
    parser.add_argument("-p", "--port",  dest="portnum",  type=int,             default=80,                                   help="port number to run on")
    parser.add_argument("-q", "--quiet", dest="verbose",  action="store_false",                                               help="don't print status messages to stdout")
    parser.add_argument("remainder", nargs=argparse.REMAINDER) # get the REMAINING un-parsed arguments (for example, a bunch of filenames)
    args = parser.parse_args()

    if not len(args.remainder) >= 1:
        print("You should really specify at least one file to this script, but oh well!")
        pass

    print("There were also this many un-parsed command line arguments: " + str(len(args.remainder)))
    for i,extra_argument in enumerate(args.remainder):
        print("The extra argument number " + str(i) + " was: " + extra_argument)
        filename = extra_argument
        if os.path.isfile(filename):
            try:
                with open(filename, str("r")) as fff:
                    for linenum,line in enumerate(fff):
                        ldelim = line.split("\t")
                        if (linenum % 10 == 0):
                            print("Writing every 10th line...")
                            pass
                        if linenum > 100:
                            print("After 100 lines, we stop reading this file.")
                            break

                        pass  # end 'for'
                    pass  # end 'with'
            except Exception as e:
                print("Failed to open the example test file. Normally we should re-raise this exception, probably. Exception: " + str(e))
                #raise #future.utils.raise_with_traceback(e)
                pass  # raise
            pass
        pass

    for x in range(10):
        sys.stdout.write(str(x)+"\t")
        pass
    print("")

    for a,b in zip([1,2,3],[10,20,30,40]): # note that the 40 is ***OMITTED***!!!
        print("{a}, {b}".format(a=str(a), b=str(b)))
        pass

    #pdb.set_trace()
    d = dict()
    for k,v in future.utils.iteritems(d): # <-- python 2/3-compatible version of 'iteritems'
        print(str(k) + str(v))
        pass

    byte_str = b'This is a BYTE string, not a unicode one! We rarely want to use this feature.'

    assert isinstance("mystring", basestring), "Uh oh, something went wrong if this gets triggered"
    
    #
    print("Handled the command line arguments!")
    return # end of 'main'

# Must come at the VERY END!
if __name__ == "__main__":
    main()
    pass





