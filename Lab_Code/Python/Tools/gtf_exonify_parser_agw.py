#!/usr/bin/python

USAGE = '''%prog [options] <infile.txt>   >   table.tab

%prog ????????

Example ??????????
   * ????????
'''

# Note: we use the 'optparse' module to parse command line arguments, becuase the superior 'argparse' requires Python 2.7 (which is not present in older distributions, like Ubuntu 10.x) http://docs.python.org/release/2.5.2/lib/module-optparse.html
import textwrap
import sys
import optparse
import os
import pdb; #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques
import re
#import os.path

globalOptions = None
globalArgs    = None

def handleCommandLineOptions():
    global globalArgs
    global globalOptions

    parser = optparse.OptionParser(USAGE, version='%prog version 1.0, makes exon views for transcripts from a GTF.')
    #parser.add_option("-d", "--delim", dest="delim", default="\t", type="string", help="Specify the delimiter character. Default is a <tab>.")
    #parser.add_option("-p", "--pad", dest="pad", default="", type="string", help="Specify the text to print in a padded-out cell. Default is nothing. 'NA' is a popular option.")
    (globalOptions, globalArgs) = parser.parse_args()

    #pdb.set_trace()
    for attr, value in globalOptions.__dict__.iteritems():
        #print attr, value
        pass

    if len(globalArgs) < 1:
        sys.stderr.write("Since no filenames were specified, we are reading from STDIN...\n")
        #parser.error("We need one un-parsed argument (a filename to read). This program is NOT currently able to handle reading from STDIN (standard in)!")
        pass

    if len(globalArgs) > 1:
        parser.error("You can only specify ONE filename on the command line! You specified more than one.")
        pass

    return

# Must come at the VERY END!
if __name__ == "__main__":
    handleCommandLineOptions()
    readingFromSTDIN = ((len(globalArgs) == 0) or (globalArgs[0] == '-')) ## "is a filename present on the command line? If not, then we're reading from STDIN"
    gtfDelim = "\t"
    gtfFilename = "<STDIN>" if readingFromSTDIN else globalArgs[0]
    fff = sys.stdin if readingFromSTDIN else open(gtfFilename, 'r')

    # Note to programmers: don't print any diagnostic messages to STDOUT, use STDERR instead!

    # Columns in a GTF file look like this:
    # 0000    11111111111111111111    2222    3333333 4444444 5       6       7       88888888888888888888888888888888888
    # chr1    processed_transcript    exon    3195982 3197398 .       -       .       exon_number "2"; gene_biotype "protein_coding"; gene_id "ENSMUSG00000051951"; gene_name "Xkr4"; transcript_id "ENSMUST00000162897"; transcript_name "Xkr4-003"; tss_id "TSS49891";
    FEATURE_TYPE_COLUMN = 2 # i.e., "exon" or whatever
    EXON_LEFT_COLUMN = 3  # Numbered from ZERO, not 1!
    EXON_RIGHT_COLUMN = 4
    STRAND_COLUMN     = 6
    FREE_COLUMN       = 8

    exonNumPat = re.compile("exon_number \"(\d+)\"", re.IGNORECASE)
    geneIdPat = re.compile("gene_id \"([^\"]+)\"", re.IGNORECASE)
    tNamePat = re.compile("transcript_name \"([^\"]+)\"", re.IGNORECASE)


    gdict = {};

    lineNum = 0;
    for lineStr in fff:
        lineNum += 1
        s = lineStr.split(gtfDelim)

        if (len(s) < (FREE_COLUMN+1)):
            sys.stderr.write("WARNING: line " + str(lineNum) + " in the input GTF file \"" + gtfFilename + "\" does not have " + str(FREE_COLUMN+1) + " elements as we expected it to. Skipping it. This may indicate a malformed GTF input file!\n")
            continue # Don't do anything with this weird malformed line
            pass

        if (s[FEATURE_TYPE_COLUMN] != 'exon'):
            sys.stderr.write("Skipping non-exon on line " + str(lineNum) + "...\n")
            continue

        freeText = s[FREE_COLUMN]
        strand = s[STRAND_COLUMN]
        posLeft = int(s[EXON_LEFT_COLUMN])
        posRight = int(s[EXON_RIGHT_COLUMN])

        exSearch = exonNumPat.search(freeText)
        gSearch = geneIdPat.search(freeText)
        tSearch = tNamePat.search(freeText)

        ex  = exSearch.group(1) if exSearch else '<somehow, no exon number!>' # ".group(1)" (not zero!) is how you get the first captured text in parens
        gid = gSearch.group(1) if gSearch else '<somehow, no gene id!>'
        tn  = tSearch.group(1) if tSearch else '<somehow, no transcript id!>'

        if (not exSearch or not gSearch or not tSearch):
            sys.stderr.write("Skipping the odd/invalid line " + str(lineNum) + "...")
            continue
            pass

        if (gid not in gdict): # gid = gene ID
            gdict[gid] = {} # new hash...
            pass

        if (tn not in gdict[gid]): # tn = transcript ID
            gdict[gid][tn] = {} # new hash...
            pass

        gdict[gid][tn][ex] = {"left":posLeft, "right":posRight, "size":abs(posLeft-posRight), "num":int(ex)}

        dbug = False

        if dbug:
            sys.stdout.write(s[0] + "\t")
            sys.stdout.write(ex + "\t")
            sys.stdout.write(gid + "\t")
            sys.stdout.write(tn + "\t")
            sys.stdout.write("\n")
            pass
        

        #sys.stdout.write(lineStr.rstrip('\n') + ((globalOptions.delim + globalOptions.pad)*numElementsToAddToThisLine) + "\n")
        pass

    
    fff.close()


    for gkey in gdict:
        sys.stdout.write(":" + gkey + ":"  + "\n")
        for tkey in gdict[gkey]:
            sys.stdout.write(" * " + tkey + ":"  + "\n")
            for ekey,edict in gdict[gkey][tkey].iteritems():
                sys.stdout.write("       E " + ekey + ":" + str(edict['size']) + " bases \n")


                pass # end "for exon key"
            pass # end "for transcript key"
        pass # end "for gene key"


    pass

