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
import operator
#import os.path

globalOptions = None
globalArgs    = None


MINIMUM_EXON_SIZE_TO_PRINT = 4
INTRON_SPACING_TO_PRINT    = 1

def warnToConsole(msg):
    sys.stderr.write(msg + "\n")
    return

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
        warnToConsole("[Note] Since no filenames were specified, we are reading from STDIN...")
        #parser.error("We need one un-parsed argument (a filename to read). This program is NOT currently able to handle reading from STDIN (standard in)!")
        pass

    if len(globalArgs) > 1:
        parser.error("[ERROR in arguments] You can only specify ONE filename on the command line! You specified more than one.")
        pass

    return

def initNonexistent(theHash, theKeyToCheck, initThing):     # Puts "initThing" into theHash if theHash doesn't ALREADY contain something at theKeyToCheck
    if (theKeyToCheck not in theHash):
        theHash[theKeyToCheck] = initThing
    return

def getAllExonKeysForGene(theGeneHash):
    assert(isinstance(theGeneHash, dict))
    assert(isinstance(theGeneHash['totalExons'], int))
    xkeys = []
    for i in range(theGeneHash['totalExons']):
        xkeys.append( str(i+1) )
        pass
    return xkeys # returns something like "1","2","3","4"...etc up to the total number of exon keys for this gene in ANY transcript

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
    GENERAL_TYPE_COLUMN = 1
    SUBFEATURE_TYPE_COLUMN = 2 # i.e., "exon" or whatever
    EXON_LEFT_COLUMN = 3  # Numbered from ZERO, not 1!
    EXON_RIGHT_COLUMN = 4
    STRAND_COLUMN     = 6
    FREE_COLUMN       = 8

    exonNumPat = re.compile("exon_number \"(\d+)\"", re.IGNORECASE)
    geneIdPat = re.compile("gene_id \"([^\"]+)\"", re.IGNORECASE)
    tIdPat = re.compile("transcript_id \"([^\"]+)\"", re.IGNORECASE)


    gdict = {};

    lineNum = 0;
    for lineStr in fff:
        lineNum += 1
        s = lineStr.split(gtfDelim)

        if (len(s) < (FREE_COLUMN+1)):
            warnToConsole("WARNING: line " + str(lineNum) + " in the input GTF file \"" + gtfFilename + "\" does not have " + str(FREE_COLUMN+1) + " elements as we expected it to. Skipping it. This may indicate a malformed GTF input file!")
            continue # Don't do anything with this weird malformed line

        if (s[SUBFEATURE_TYPE_COLUMN] != 'exon'):
            warnToConsole("Skipping non-exon on line " + str(lineNum) + "...")
            continue

        freeText = s[FREE_COLUMN]
        strand = s[STRAND_COLUMN]
        posLeft = int(s[EXON_LEFT_COLUMN])
        posRight = int(s[EXON_RIGHT_COLUMN])
        generalType = s[GENERAL_TYPE_COLUMN]
        
        exSearch = exonNumPat.search(freeText)
        gSearch = geneIdPat.search(freeText)
        tSearch = tIdPat.search(freeText)

        ex  = exSearch.group(1) if exSearch else '<somehow, no exon number!>' # ".group(1)" (not zero!) is how you get the first captured text in parens
        gid = gSearch.group(1) if gSearch else '<somehow, no gene id!>'
        tid  = tSearch.group(1) if tSearch else '<somehow, no transcript id!>'

        if (not exSearch or not gSearch or not tSearch):
            warnToConsole("Skipping the odd/invalid line " + str(lineNum) + "...")
            continue
            pass

        initNonexistent(gdict                              , gid, {'totalExons':-999, 'sizesPerExonHash':{}, 'trHash':{} }) # new hash...
        initNonexistent(gdict[gid]['trHash']               , tid, {'strand':strand, 'general_type':generalType, 'trSize':-999, 'exHash':{} })
        initNonexistent(gdict[gid]['trHash'][tid]['exHash'], ex , {"left":posLeft, "right":posRight, "exSize":abs(posLeft-posRight), "exRelativeSize":-999, "num":int(ex)} )

        dbug = False
        if dbug:
            sys.stdout.write(s[0] + "\t")
            sys.stdout.write("Ex#" + ex + "\t")
            sys.stdout.write(gid + "\t")
            sys.stdout.write(tid + "\t")
            sys.stdout.write("\n")
            pass
        

        #sys.stdout.write(lineStr.rstrip('\n') + ((globalOptions.delim + globalOptions.pad)*numElementsToAddToThisLine) + "\n")
        pass
    
    fff.close()

    # Calculate the TRANSCRIPT TOTAL SIZES and MAX NUM EXONS FOR THIS GENE by adding up the exons
    for gKey,gHash in gdict.iteritems():
        #sys.stdout.write(":" + gKey + ":"  + "\n")
        for tKey,trHash in gdict[gKey]['trHash'].iteritems():
            tStrand = trHash['strand']
            tType   = trHash['general_type']
            for eKey,exHash in trHash['exHash'].iteritems():
                initNonexistent(gdict[gKey]['sizesPerExonHash'], eKey, []) # new list
                gHash['sizesPerExonHash'][eKey].append( { 'size':exHash['exSize'], 'assocTr':tKey } ) # save this exon's size and which transcript it belongs to
                pass
                # end "for eKey,exHash..."
            pass
        # end "for tKey,trHash..."
        gHash['totalExons'] = len(gHash['sizesPerExonHash']) # number of unique numbered/named exons that are theoretically possible for a transcript, even if no single transcript actually contains all the named/numbered exons
        for eKey in getAllExonKeysForGene(gHash):
            sizesList        = gHash['sizesPerExonHash'][eKey]
            orderedSizesList = sorted(sizesList, key=operator.itemgetter('size'))  # sorted by size
            sortedUniqSizesList = sorted(set(map(lambda k: k['size'], sizesList)))
            for xHash in gHash['sizesPerExonHash'][eKey]:
                thisSize = xHash['size']
                associatedTranscript = xHash['assocTr']
                assert(associatedTranscript in gHash['trHash'])
                relativeSizeOrderForThisExon = sortedUniqSizesList.index(thisSize) # get the RELATIVE ranking of sizes for this exon. 0 = shortest exon, 1 = slightly longer, etc etc. If all exons are the same size, they will all have "0" as their value.
                assert(relativeSizeOrderForThisExon >= 0)
                gHash['trHash'][associatedTranscript]['exHash'][eKey]['exRelativeSize'] = relativeSizeOrderForThisExon
                pass
            pass
        # end "for eKey..."
        #sys.stdout.write("Gene " + gKey + " had " + str(gdict[gKey]['totalExons']) + " exons in the most-exony transcript.\n")
        pass
     # end "for gene key"


    for gKey,gHash in gdict.iteritems():
        for tKey,trHash in gHash['trHash'].iteritems():
            tStrand = trHash['strand']
            tType   = trHash['general_type']
            maxExonsThisGene = gHash['totalExons']

            sys.stdout.write(gKey + "\t" + tKey + "\t")
            for eKey in getAllExonKeysForGene(gHash):
                hasThisExon = eKey in trHash['exHash']
                listOfAllSizesForThisExon = map(lambda k: k['size'], gHash['sizesPerExonHash'][eKey])
                numDifferentLengthsForThisExon = len(set(listOfAllSizesForThisExon))
                biggestExonRank = (numDifferentLengthsForThisExon-1)
                totalExonSpace = (biggestExonRank + MINIMUM_EXON_SIZE_TO_PRINT)
                xToPrint       = 0 if not hasThisExon else trHash['exHash'][eKey]['exRelativeSize'] + MINIMUM_EXON_SIZE_TO_PRINT
                blankToPrint   = totalExonSpace - xToPrint
                leadingIntronToPrint   = 0 if eKey == "1" else INTRON_SPACING_TO_PRINT
                sys.stdout.write(("_" * leadingIntronToPrint) + "X"*xToPrint + "-"*blankToPrint)
                pass # end "for exon key"
            sys.stdout.write("\n") # This part is after this transcript's exons are ALL processed
            pass # end "for transcript key"
        pass # end "for gene key"


    pass # end of the <<if __name__ == "__main__">>  part

