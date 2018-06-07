#!/usr/bin/env python3
'''
This script requires python3
'''

# import ipdb; ipdb.set_trace()
import argparse
import os
import sys
#import pdb   #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques
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
    for k,v in d.items(): # <-- python 2/3-compatible version of 'iteritems'
        print(str(k) + str(v))
        pass

    byte_str = b'This is a BYTE string, not a unicode one! We rarely want to use this feature.'

    assert isinstance("mystring", str), "Uh oh, something went wrong if this gets triggered"
    
    #
    print("Handled the command line arguments!")
    return # end of 'main'

# Must come at the VERY END!
if __name__ == "__main__":
    main()
    pass





