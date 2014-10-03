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



def initNonexistent(theHash, theKeyToCheck, initThing):
    # Puts "initThing" into theHash if theHash doesn't ALREADY contain something at theKeyToCheck
    if (theKeyToCheck not in theHash):
        theHash[theKeyToCheck] = initThing
        pass
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
            sys.stderr.write("WARNING: line " + str(lineNum) + " in the input GTF file \"" + gtfFilename + "\" does not have " + str(FREE_COLUMN+1) + " elements as we expected it to. Skipping it. This may indicate a malformed GTF input file!\n")
            continue # Don't do anything with this weird malformed line

        if (s[SUBFEATURE_TYPE_COLUMN] != 'exon'):
            sys.stderr.write("Skipping non-exon on line " + str(lineNum) + "...\n")
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
            sys.stderr.write("Skipping the odd/invalid line " + str(lineNum) + "...\n")
            continue
            pass

        initNonexistent(gdict                              , gid, {'totalExons':-999, 'sizesPerExon':{}, 'trHash':{} }) # new hash...
        initNonexistent(gdict[gid]['trHash']               , tid, {'strand':strand, 'general_type':generalType, 'trSize':-999, 'exHash':{} })
        initNonexistent(gdict[gid]['trHash'][tid]['exHash'], ex , {"left":posLeft, "right":posRight, "exSize":abs(posLeft-posRight), "num":int(ex)} )

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
        sys.stdout.write(":" + gKey + ":"  + "\n")
        for tKey,trHash in gdict[gKey]['trHash'].iteritems():
            tStrand = trHash['strand']
            tType   = trHash['general_type']
            for eKey,exHash in trHash['exHash'].iteritems():
                initNonexistent(gdict[gKey]['sizesPerExon'], eKey, []) # new list
                gHash['sizesPerExon'][eKey].append( { 'size':exHash['exSize'], 'assocTr':tKey, 'sizeOrder':-999 } ) # save this exon's size and which transcript it belongs to
                # 'sizeOrder' is the exon's relative size in THIS transcript compared to the shortest transcript: 0 = shortest, 1 = slightly longer, etc... if this exon is the same length in all transcripts, then all values will be zero!
                pass
                # end "for eKey,exHash..."
            pass
        # end "for tKey,trHash..."
        gHash['totalExons'] = len(gHash['sizesPerExon']) # number of unique numbered/named exons that are theoretically possible for a transcript, even if no single transcript actually contains all the named/numbered exons

        for eKey in getAllExonKeysForGene(gHash):
            sizesList        = gHash['sizesPerExon'][eKey]
            orderedSizesList = sorted(sizesList, key=operator.itemgetter('size'))  # sorted by size
            print sizesList
            print orderedSizesList
            sortedUniqSizesList = sorted(set(map(lambda k: k['size'], sizesList)))
            print sortedUniqSizesList
            relativeExonSizeOrdering = []

            for xHash in gHash['sizesPerExon'][eKey]:
                thisSize = xHash['size']
                xHash['sizeOrder'] = relativeExonSizeOrdering.append( sortedUniqSizesList.index(thisSize)  ) # get the RELATIVE ranking of sizes for this exon. 0 = shortest exon, 1 = slightly longer, etc etc. If all exons are the same size, they will all have "0" as their value.
                pass
            #print sizesForThisExon
            #print sortedUniqueSizes
            #print relativeExonSizeOrdering
            print "^^^^^^^^^^^^^"
            #assert(len(relativeExonSizeOrdering) == len(sizesForThisExon))
            pass
        # end "for eKey..."
        sys.stdout.write("Gene " + gKey + " had " + str(gdict[gKey]['totalExons']) + " exons in the most-exony transcript.\n")
        pass
     # end "for gene key"


    for gkey in gdict:
        maxExonsSeenSoFar = 0
        sys.stdout.write(":" + gkey + ":"  + "\n")
        for tKey,trHash in gdict[gkey]['trHash'].iteritems():
            sys.stdout.write("   * " + tKey + ":"  + "\n")
            tStrand = trHash['strand']
            tType   = trHash['general_type']
            maxExonsThisGene = gdict[gkey]['totalExons']
            for i in range(maxExonsThisGene):
                exonIdentifier = str(i+1)
                hasThisExon = exonIdentifier in trHash['exHash']
                if (hasThisExon):
                    sys.stdout.write("--XXXX--")
                else:
                    sys.stdout.write("--------")
                    pass
                pass # end "for exon key"

            # This part is after this transcript's exons are ALL processed
            sys.stdout.write("\n")

            pass # end "for transcript key"
        pass # end "for gene key"


                #for eKey,exHash in trHash['exHash'].iteritems():
                #sys.stdout.write("       Exon #" + eKey + ":" + str(exHash['exSize']) + " bases (strand " + tStrand + ") -- general type is " + tType +" \n")



    pass # end of the <<if __name__ == "__main__">>  part
