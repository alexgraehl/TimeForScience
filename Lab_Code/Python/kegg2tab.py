#!/usr/bin/python
"""KEGG XML relation parser

    Converts a KGML file into a tab-delimited file with one relation
    per line.  Output is written to standard output.  If no arguments
    are given, input is read from standard input.

Usage:

    kegg2tab.py [options] kegg_input.xml

Options:

    -h, --help    print this usage statement
"""
from xml.dom import minidom
from xml.parsers.expat import ExpatError
import getopt
import sys

class KgmlLink:
    def __init__(self, domelement):
        self.entry1 = domelement.getAttribute("entry1")
        self.entry2 = domelement.getAttribute("entry2")
        self.type   = domelement.getAttribute("type")
        self.subtype = ''
        self.subvalue = ''
        if len(domelement.getElementsByTagName("subtype")) == 1:
            s = domelement.getElementsByTagName("subtype")[0]
            self.subtype = s.getAttribute("name")
            self.subvalue = s.getAttribute("value")            

class KgmlRelationships:
    def __init__(self, source=None):
        self.entities = {}
        self.links = []
        self.accession = ''
        self.title = ''
        try:
            xmldoc = minidom.parse(getInputStream(source)).documentElement
            self.loadXML(xmldoc)
        except ExpatError:
            print "Couldn't parse an XML document from:", source
            return None
    def loadXML(self, xmldoc):
        self.accession = xmldoc.getAttribute("name")
        self.title = xmldoc.getAttribute("title")
        self.parseEntities(xmldoc)
        self.parseRelationships(xmldoc)
    def parseEntities(self, xmldoc):
        for n in xmldoc.getElementsByTagName("entry"):
            if n.hasAttribute("id") and n.hasAttribute("name"):
                names = n.getAttribute("name").split()
                self.entities[n.getAttribute("id")] = names
    def parseRelationships(self, xmldoc):
        for n in xmldoc.getElementsByTagName("relation"):
            self.links.append(KgmlLink(n))
    def linksByAllNames(self):
        """Return a link element for the cross of all names named in a
        relationship.  For example if there are two names in the entry
        referenced by entry1 and 3 names in the entry referenced by
        entry2 for a relationship, then this function will yield six
        separate links.  Also, this function will iterate over names
        for compounds referred to by the \"compound\" and \"hidden
        compound\" subtypes."""
        for l in self.links:
            if l.entry1 not in self.entities:
                print "WARNING: couldn't find referenced entry", l.entry2
            elif l.entry2 not in self.entities:
                print "WARNING: couldn't find referenced entry", l.entry2
            else:
                for n1 in self.entities[l.entry1]:
                    for n2 in self.entities[l.entry2]:
                        if (l.subtype in ("compound", "hidden compound")
                            and l.subvalue in self.entities):
                            for c in self.entities[l.subvalue]:
                                yield [n1, n2, l.type, l.subtype, c]
                        else:
                            yield [n1, n2, l.type, l.subtype, l.subvalue]
    def allEntryNames(self):
        for entry in self.entities.values():
            for name in entry:
                yield name


def getInputStream(source):
    """Try to open a file, first using stdin if the source is the
    string -, second trying the file as a URL, third try as a regular
    file on a file system, and finally treat it just as a raw string.
    Based on toolbox.openAnything() at
    http://diveintopython.org/xml_processing/index.html"""
    if hasattr(source, "read"):
        return source
    import urllib
    if source == "-":
        return sys.stdin
    try:
        return urllib.urlopen(source)
    except (IOError, OSError):
        pass
    try:                                  
        return open(source)               
    except (IOError, OSError):            
        pass                              
    import StringIO                       
    return StringIO.StringIO(str(source)) 
        
def usage():
    print __doc__

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hn", ["help"])
    except getopt.GetoptError:
        print "Error parsing arguments\n"
        usage()
        sys.exit(2)
    for opt, arg, in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
    if len(argv) == 0:
        argv = ["-"]
    for a in args:
        kgml = KgmlRelationships(a)
        if ('-n', '') in opts:
            for n in kgml.allEntryNames():
                print "\t".join([n, kgml.accession, kgml.title])
            print opts
            print args
        else: 
            for l in kgml.linksByAllNames():
                print "\t".join(l + [kgml.accession, kgml.title])

if __name__ == "__main__":
    main(sys.argv[1:])
