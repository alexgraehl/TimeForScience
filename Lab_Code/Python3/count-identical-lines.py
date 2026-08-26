#!/usr/bin/env python3

import getopt
import sys
import textwrap
from collections import defaultdict

# import pdb; #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques

ALEX_PROGRAM_USAGE_TEXT = """
Duplicate line counter by Alex Williams. Intended For FASTA/FASTQ/CSFASTA files.

Provides the number of times a specific line is seen in a file, for example, a file that looks like this:

```text
ZZZZ
A
A
A
ZZZZ
BB
BB
CCC
A
```

Would result in this output, where an additional (leftmost) column now shows the number of times a specific line was seen.

```text
2  ZZZZ
4  A
2  BB
1  CCC
```

The original line order is preserved, except that duplicates are removed.

This is a niche utility: in most cases, you'd be better off just running 'sort YOURFILE | uniq -c', unless you
actually care about preserving the order of the input lines.

For certain large FASTQ files (benchmarked on a 2.2 GB file), this is about 5x faster.

Usage with a file input:
   count-identical-lines.py YOURFILE.txt > out.txt

Usage when the data is piped from STDIN:
   cat YOURFILE.txt | count-identical-lines.py > out.txt
   (Remember to redirect STDOUT to a file (that is, don't forget the '> out.txt' part of the command).)

You may have to use "cut" beforehand if you only want ONE column out of a file.
"""

TERMINAL_WIDTH = 80

def usageAndQuit(exitCode: int, message: str="") -> None:
    message = textwrap.fill(message, TERMINAL_WIDTH)
    if message is not None:
        print("")
        print("simulator: ")
        print(message)
        print("simulator: Printing usage information below.")
        print("*" * TERMINAL_WIDTH)
        print("*" * TERMINAL_WIDTH)
        pass
    print(ALEX_PROGRAM_USAGE_TEXT)  # at the very bottom of this file
    if message is not None:
        print("(End of usage information)")
        print("*" * TERMINAL_WIDTH)
        print("*" * TERMINAL_WIDTH)
        print("simulator: " + message)
        print("*" * TERMINAL_WIDTH)
        print("*" * TERMINAL_WIDTH)
        print("[Program Terminated]")
        pass
    sys.exit(exitCode)  # No return value

if __name__ == "__main__":
    sys.stderr.write(
        "## FYI: Note that you can get similar but SORTED counts with the UNIX command --> sort YOURFILE | uniq -c > OUTPUT_FILE\n"
        "## This script is about 5x faster, and may be more useful if you want the counts in ORIGINAL line order (of first line occurrence).\n"
    )
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
        print("Unprocessed arguments:", args)
        print("Unprocessed options:", opts)
        pass

    if len(args) == 0 or args[0] == "-":  # hyphen means 'read from stdin'
        theFile = sys.stdin
    else:
        if len(args) != 1:
            print(args)
            print(len(args))
            usageAndQuit(1, "ARGUMENTS ERROR: We need exactly ONE filename passed in on the command line.")
            raise

        inputFilename = args[0]
        try:
            theFile = open(inputFilename, "r")
        except:
            sys.stderr.write("ERROR: Could not open the specified input file: " + inputFilename)
            raise
        pass

    ddd = defaultdict(int)  # (FYI: 0 is the default for int)
    for lineNum, line in enumerate(theFile):
        ddd[line] += 1  # Increment the number of "duplicates" (possibly the first occurrence)
        pass
    theFile.close()

    # Note that we write to 'sys.stdout' here to make it easier to redirect the (non-stderr) output.
    for key, value in ddd.items():
        # rstrip removes whitespace from right side of the key. This is important!
        sys.stdout.write(str(value) + "\t" + key.rstrip() + "\n")  
        pass

    # These status messages are intentionally written to STDERR.
    sys.stderr.write(f"[Done -- Read a total of {lineNum} lines.]\n")
    pass
