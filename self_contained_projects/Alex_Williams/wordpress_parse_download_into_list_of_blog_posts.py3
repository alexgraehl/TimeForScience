#!/usr/bin/env python3
'''
Requires python3.

'''


''' Parses a wordpress archive and makes an index page for it. '''
''' Currently specific to worst ideas, but should basically work for anything. '''

# how to make a border=0 image on wordpress:   put this into the IMG in the HTML editor:  style="border:0"

import argparse
import os
import sys
#import ipdb  # note: ipdb also prints weird stuff to STDOUT
import html
#import textwrap
from collections import defaultdict

def argErrorAndExit(msg="(No additional information given)"):
    raise SystemExit("[ERROR] in arguments to this script: " + msg)

def main():
    postdict   = dict()
    catdict = defaultdict(list) # key = category, value = list with TITLES of the thing in pdict with this key
    parser = argparse.ArgumentParser(description="%(prog)s: template for python 2 and python 3.",
                                     epilog='''Example usage: python %(prog)s (no example yet)\n(some examples go here)\n(More examples go here)''',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--file",  dest="filename", type=str,             default=None, metavar="FILE", required=False, help="Minimum x (length) value.")
    parser.add_argument("-N", "--name",  dest="username", type=str,             default="Dave",               required=False, help="specify a username to run as")
    parser.add_argument("-p", "--port",  dest="portnum",  type=int,             default=80,                                   help="port number to run on")
    parser.add_argument("-q", "--quiet", dest="verbose",  action="store_false",                                               help="don't print status messages to stdout")
    parser.add_argument("remainder", nargs=argparse.REMAINDER) # get the REMAINING un-parsed arguments (for example, a bunch of filenames)
    args = parser.parse_args()

    if not len(args.remainder) == 1:
        raise Exception("You must specify exactly one XML file to this script!")

    filename = args.remainder[0]
    if not os.path.isfile(filename): raise Exception(f"Failed to open file <{filename}>")

    import re
    import xml.etree.ElementTree as ET
    orig_xstr            = open(filename, encoding="UTF-8").read()
    invalid_xml          = re.compile(u'[\x00-\x08\x0B-\x0C\x0E-\x1F\x7F]')
    (xstr, n_invalid) = invalid_xml.subn('', orig_xstr)
    if n_invalid > 0:
        print(f"Removed {n_invalid} invalid characters from input XML (usually \\x00).", file=sys.stderr)
        pass
    root = ET.fromstring(xstr) #, parser=parser)
    #print(root)
    #print("Parse time...")

    # topmost thing is the CHANNEL
    # then we have the ITEMS

    def print_elem(elem, indent=0):
        print("    "*indent + "«"+str(elem.tag) + "»" + ": " + str(elem.text))  # , d.attrib)

    #things_we_want = {"title":"human-readable title"
    #                  , "link":"URL for this thing"
    #                  , "pubDate":"publication date in a weird format, like Mon, 20 Apr 2015 12:00:00 +0000"
    #                  }

    ignore_us = {"{http://wordpress.org/export/1.2/}tag"
        , "{http://wordpress.org/export/1.2/}category"
                 , "{http://wordpress.org/export/1.2/}author"
                 }

    #def find_or_this_text(el, findme, not_found_text)
    #    x = el.find(findme)
    #    if x is None: return not_found_text
    #    else return

    def htmlify(s):
        s = html.escape(s, quote=True)
        s = re.sub(r"[“]", "&ldquo;", s) # html escape DOES NOT handle these for some reason
        s = re.sub(r"[”]", "&rdquo;", s)
        s = re.sub(r"[‘]", "&lsquo;", s)
        s = re.sub(r"[’]", "&rsquo;", s)
        s = re.sub(r"[™]", "&trade;", s)
        s = re.sub(r"[–]", "&ndash;", s)
        s = re.sub(r"[—]", "&mdash;", s)
        s = re.sub(r"((\S)|^)&quot;", "\\1&rdquo;", s) # >>"like this" <-- the first one
        s = re.sub(r"&quot;((\S)|$)", "&ldquo;\\1", s) # "like this"<< <-- the second one
        s = re.sub(r"&#x27;", "&rsquo;", s) # the apostrophe in "can't". &apos; is rendered as a STRAIGHT QUOTE for some reason in Chrome.
        s = re.sub(r"[š]", "&scaron;", s)
        s = re.sub(r"[æ]", "&aelig;", s)
        s = re.sub(r"[Æ]", "&AElig;", s)
        return s


    POST_TYPE_TAG   = "{http://wordpress.org/export/1.2/}post_type"
    POST_STATUS_TAG = "{http://wordpress.org/export/1.2/}status"  # throw it out if it's not "publish" (e.g. if it's draft)
    def delve(elem, indent=0): # <-- recursively delves into the XML tree
        if elem.tag in ignore_us:
            return # we don't care about this tag (or its children!), so continue

        #print_elem(elem, indent)
        #ipdb.set_trace()
        # we want to include <wp:post_type>post
        # but ignore <wp:post_type>attachment
        # and probably ignore <wp:post_type>page
        pt = elem.find(POST_TYPE_TAG)
        if pt is not None and pt.text not in ["post"]:
            #print(f"#COMMENT#  Skipping a non-post item, and all of its children. The post type was: {pt.text}")
            return

        pstat = elem.find(POST_STATUS_TAG)
        if pstat is not None and pstat.text not in ["publish"]:
            #print(f"#COMMENT# Skipping a non-'publish' item, and all of its children. The status was: {pstat.text}")
            return

        if str(elem.tag) == "item":
            # ok, we should fish out the things I want:
            title = elem.find("title").text
            title = htmlify(title)
            link  = elem.find("link").text
            #date  = elem.find("pubDate").text
            date = elem.find("{http://wordpress.org/export/1.2/}post_date").text
            category_elements_list = elem.findall("category")
            categories = [c.text for c in category_elements_list ]

            from datetime import datetime
            try:
                #print("Got this date: " + date)
                dobj = datetime.strptime(date, "%Y-%m-%d %H:%M:%S")
                #2017-12-04 05:25:44
                #dobj = datetime.strptime(date,               "%a, %d %b %Y %H:%M:%S %z")
                # Expects the input to look like this --> Fri, 04 Jul 2014 08:24:48 +0000
            except ValueError as e:
                #print(f"Failed to parse this date: {date}       <-- Exception message: {str(e)}")
                dobj = None
                pass
            #print(str(dobj))
            #print(f"""<LI><A HREF="{link}">{title}</A> ({date})""")
            postdict[title] = {"link":link, "title":title, "date":dobj, "categories":categories}
            for one_category in categories:
                catdict[one_category].append(title)
                pass

            #print("----------")

            #if g is not None:
            #    print(elem.find(t))
            pass

        #if str(elem.tag) == "title":
        #    # ok guys, this is it! this is the title
        #    print("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ")

        for child in elem:
            delve(child, indent+1) # recursively delve into the tree
            pass

        return # end of 'delve'

    delve(root) # <-- recursively delves into the XML tree

    # Old iterative way of debugging this
    # for b in root:
    #     pelem(b, 1)
    #     for c in b:
    #         pelem(c, 2)
    #         for d in c:
    #             pelem(d, 3)
    #             pass
    #         pass
    #     pass
    #import ipdb
    #ipdb.set_trace()

    #sorted_dates = [datetime.datetime.strptime(ts, "%Y-%m-%d") for ts in timestamps]

    def print_item_to_html(elem):
        #datestr = elem['date'].strftime("%B %-d, %Y") # e.g. "May 4, 2017" (%-d means day of month WITHOUT a leading zero!)
        datestr = elem['date'].strftime("%b %d, %Y") # e.g. "Apr 04, 2017"
        #ipdb.set_trace()
        print(f"""<LI><CODE>{datestr}</CODE>: <A HREF="{elem['link']}">{elem['title']}</A>""")
        pass

    for title,elem in postdict.items():
        pass

    clist = [x for _,x in sorted(zip([len(v) for v in catdict.values()], catdict.keys()))]
    clist.remove("Uncategorized") # remove this category

    top_categories = ['Film / TV', 'Finance', 'Education', 'Architecture', 'Beasts', 'Art', 'Diet', 'Fitness', 'Language', 'Games', 'Politics', 'Computers', 'Law', 'Health', 'Household', 'Automotive', 'Small Business', 'Culture', 'Technology', 'UI / UX', 'Design', 'Communication', 'Fashion', 'Government']
    #print(clist)

    titles_we_saw_set = set() # have we seen this title already?

    num_titles_total = len(postdict)

    print(f"<H1>Posts Listed by Category</H1>")
    print(f"<UL>")
    print(f"""<LI>Total number of unique posts: {num_titles_total}""")
    print(f"</UL>")

    print(f"<P>")
    def print_entire_category(category, post_title_list):
        titleized = category.title()
        if titleized == "Ui / Ux": titleized = "User Interface (UI / UX)"
        if titleized == "Film / Tv": titleized = "Film and Television"
        if titleized == "Automotive": titleized = "Automotive / Cars / Traffic"
        if titleized == "Diet": titleized = "Food &amp; Diet"
        if titleized == "Games": titleized = "Games (Board Games &amp; Video Games)"
        if titleized == "Beasts": titleized = "Beasts (Pets &amp; Other Animals)"
        #titleized = re.sub(" ", "&nbsp;", titleized)

        print(f"<H1>{titleized}</H1>")
        print(f"<UL>")
        for t in post_title_list:
            print_item_to_html(elem=postdict[t])
            titles_we_saw_set.add(t)
            pass
        print(f"</UL>")
        return

    for category_str in top_categories:
        print_entire_category(category=category_str, post_title_list=catdict[category_str])
        pass

    # print all the UNSEEN categories that weren't included in any of the categories above
    orphaned_titles = set(postdict.keys()).difference(titles_we_saw_set)
    print_entire_category(category="Miscellaneous", post_title_list=orphaned_titles)

    for category,_ in catdict.items():
        #print(f" * CATEGORY: {category}")
        pass

    print("Handled the command line arguments!", file=sys.stderr)
    return # end of 'main'

# Must come at the VERY END!
if __name__ == "__main__":
    main()
    pass

# ipython3 --pdb parse.py3 a.xml





