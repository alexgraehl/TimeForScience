#!/usr/bin/env python3

import sys
import curses  # <-- docs at http://docs.python.org/library/curses.html
import curses.ascii  # <-- docs at http://docs.python.org/library/curses.ascii.html
import curses.textpad
import gzip
import bz2
import re

# import pdb  # pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques

# @COMMENT@ sheet.py is an interactive terminal-based spreadsheet viewer for tab-delimited files. Frequency-of-use rating: 9/10. See also "sheet.pl" if you want a non-interactive version that behaves similarly to the UNIX "column" tool.

# from __future__ import print_function # lets you do print("...") like python 3
# from __future__ import division # no more "1/2 == 0"

# Instructions for pylint below:
# pylint: disable=line-too-long, superfluous-parens, bad-whitespace, unused-wildcard-import, trailing-whitespace, unnecessary-pass, missing-docstring, invalid-name, global-statement, multiple-statements, too-many-locals, too-many-statements, too-many-branches, too-few-public-methods, too-many-lines, too-many-instance-attributes, too-many-arguments, wildcard-import

"""
%prog%: a python3 program for viewing tab-delimited files in spreadsheet-like format. Only a viewer--you cannot edit files with it!
Uses the "CURSES" terminal interaction library to talk to the terminal. Updated for Python3 in 2022.
by Alex Williams
Example 1:   sheet.py  yourFile.tab   or else:    cat someFile | sheet.py
Example 2: "sheet.py -i ~/T" (tab-delimted file)
"""

try:
    import numpy as np  # The matrices use this. Probably this dependency should not exist.
except ImportError:
    print("[sheet.py]: Unable to load the module 'numpy'. Probably it is not installed")
    print(
        "[sheet.py]: [TERMINATING EARLY]: sheet.py cannot run, because the module 'numpy' was not installed."
    )
    print("[sheet.py]:                      On Linux, you should be able to install it with the command:")
    print("[sheet.py]:                      sudo apt-get install python-numpy")
    print("[sheet.py]: Try installing numpy, and then running sheet.py again...")
    sys.exit(1)
    pass

GLOB = 1

kNUM_BYTES_TO_READ_AT_A_TIME = (
    100000  # read 100kb at a time. Hopefully that's enough to populate at least one screen full of data!
)
# For debugging, you can set this number very low (say to 100) to see how the loading only happens every time a user moves the cursor.

kWANT_TO_ADJUST_CURSOR = 0
kWANT_TO_DIRECTLY_MOVE_CURSOR = 1
kWANT_TO_MOVE_TO_SEARCH_RESULT = 2
kWANT_TO_CHANGE_FILE = 3

ROW_HEADER_MAXIMUM_COLUMN_FRACTION_OF_SCREEN = 0.25  # The row header column (i.e., the leftmost column) cannot be any wider than this fraction of the total screen width. 1.0 means "do not change--it can be the entire screen," 0.25 means "one quarter of the screen is the max, etc. 0.5 was the default before.

WIDTH_WHEN_TRANSPOSED = (
    25  # This is a hack! Currently colwidth isn't properly computed for transposed matrices
)
MAX_COL_WIDTH = 30  # If a column is wider than this, then clip it to this width

# ===============
class SearchBoard:
    def __init__(self, nrow, ncol):
        self.hit = np.empty((nrow, ncol), dtype=bool)
        self.dirty = np.ones((nrow, ncol), dtype=bool)
        pass

    def clear(self):
        self.dirty = np.ones((self.numRows(), self.numCols()), dtype=bool)
        pass

    def set(self, row, col, value):
        self.hit[row][col] = value
        self.cleanify(row, col)  # This cell has been legitimately set...
        pass

    def matches(self, row, col):
        if self.dirty[row][col]:
            raise  # It's not a "clean" value yet
        else:
            return self.hit[row][col]
        pass

    def isDirty(self, row, col):
        return self.dirty[row][col]

    def dirtify(self, row, col):
        self.dirty[row][col] = True
        pass

    def cleanify(self, row, col):
        self.dirty[row][col] = False
        pass

    def numRows(self):
        return self.hit.shape[0]

    def numCols(self):
        return self.hit.shape[1]

    pass  # End of SearchBoard class


# ===============

import textwrap

# import time # we just want "sleep"
import os.path
import getopt

import re  # Regexp: http://docs.python.org/library/re.html

kSTDIN_HYPHEN = "-"

KEY_MODE_NORMAL_INPUT = 0
KEY_MODE_SEARCH_INPUT = 1

HEADER_NUM_DELIMITER_STRING = ": "  # the separator between "ROW_1" and "here is the header". Example: "ROW_982-->Genomes". In that case, the string would have been "-->"

KEYS_TOGGLE_HIGHLIGHT_NUMBERS_MODE = (ord("!"),)

SHIFT_UP_KEY_ID = 65  # ASCII id for "Shift key + arrow key"
SHIFT_DOWN_KEY_ID = 66
SHIFT_RIGHT_KEY_ID = 67
SHIFT_LEFT_KEY_ID = 68

FAST_UP_KEY_IN_LESS = ord("b")
FAST_DOWN_KEY_IN_LESS = curses.ascii.SP  # space

KEYS_MOVE_TO_TOP_IN_LESS = ord("g")
KEYS_MOVE_TO_BOTTOM_IN_LESS = ord("G")

FAST_DOWN_KEY_IN_EMACS = curses.ascii.ctrl(ord("v"))
FAST_UP_KEY_IN_EMACS = curses.ascii.alt(ord("v"))

GOTO_LINE_START_KEY_IN_EMACS = curses.ascii.ctrl(ord("a"))
GOTO_LINE_END_KEY_IN_EMACS = curses.ascii.ctrl(ord("e"))

KEYS_MOVE_TO_TOP = (curses.KEY_HOME, KEYS_MOVE_TO_TOP_IN_LESS)
KEYS_MOVE_TO_BOTTOM = (curses.KEY_END, KEYS_MOVE_TO_BOTTOM_IN_LESS)

KEYS_MOVE_LEFT = (curses.KEY_LEFT, ord("j"), curses.ascii.ctrl(ord("b")))
KEYS_MOVE_UP = (curses.KEY_UP, ord("i"), curses.ascii.ctrl(ord("p")))

KEYS_MOVE_DOWN = (curses.KEY_DOWN, ord("k"), curses.ascii.ctrl(ord("n")))
KEYS_MOVE_RIGHT = (curses.KEY_RIGHT, ord("l"), curses.ascii.ctrl(ord("f")))

KEYS_MOVE_RIGHT_FAST = (ord("L"), SHIFT_RIGHT_KEY_ID)
KEYS_MOVE_LEFT_FAST = (ord("J"), SHIFT_LEFT_KEY_ID)
KEYS_MOVE_UP_FAST = (
    ord("I"),
    ord("v"),
    SHIFT_UP_KEY_ID,
    curses.KEY_PPAGE,
    FAST_UP_KEY_IN_LESS,
    FAST_UP_KEY_IN_EMACS,
)
KEYS_MOVE_DOWN_FAST = (
    ord("K"),
    SHIFT_DOWN_KEY_ID,
    curses.KEY_NPAGE,
    FAST_DOWN_KEY_IN_LESS,
    FAST_DOWN_KEY_IN_EMACS,
)

KEYS_GOTO_LINE_START = (ord("a"), GOTO_LINE_START_KEY_IN_EMACS)
KEYS_GOTO_LINE_END = (ord("e"), GOTO_LINE_END_KEY_IN_EMACS)

KEYS_NEXT_FILE = (ord("."), ord(">"))
KEYS_PREVIOUS_FILE = (ord(","), ord("<"))
KEYS_TRANSPOSE = (ord("t"), ord("T"))

KEYS_QUIT = (
    ord("q"),
    ord("Q"),
    curses.ascii.ESC,
)  # , curses.ascii.ctrl(ord('q')), curses.ascii.ctrl(ord('d')), curses.ascii.ctrl(ord('c')), curses.ascii.ESC)


KEYS_GOTO_NEXT_MATCH = (ord("n"),)
KEYS_GOTO_PREVIOUS_MATCH = (ord("N"),)

KEYS_WANT_TO_ENTER_SEARCH_MODE = (ord("/"), curses.ascii.ctrl(ord("m")))

# ==== HERE ARE KEYS THAT ARE SPECIFIC TO SEARCH MODE ====
KEYS_SEARCH_MODE_FINISHED = (curses.KEY_ENTER, curses.ascii.LF, curses.ascii.CR)

KEYS_SEARCH_MODE_CANCEL = (curses.ascii.ESC, curses.ascii.ctrl(ord("g")), curses.KEY_CANCEL)

KEYS_SEARCH_MODE_BACKSPACE = (curses.KEY_BACKSPACE, curses.ascii.BS)
KEYS_SEARCH_MODE_DELETE_FORWARD = (curses.KEY_DC, curses.ascii.DEL)

# ==== END OF KEYS THAT ARE SPECIFIC TO SEARCH MODE ====

STANDARD_BG_COLOR = curses.COLOR_BLACK

# If a table ends with a "ragged end," and some rows aren't even the proper
# length, then the straggling cells get this color.
# You will see it a lot in ragged-end files, like list files.
RAGGED_END_ID = 1
RAGGED_END_TEXT_COLOR = curses.COLOR_GREEN  # BLUE #WHITE
RAGGED_END_BG_COLOR = STANDARD_BG_COLOR  # BLUE

SELECTED_CELL_ID = 2
SELECTED_CELL_TEXT_COLOR = curses.COLOR_CYAN
SELECTED_CELL_BG_COLOR = curses.COLOR_MAGENTA

COL_HEADER_ID = 3
COL_HEADER_TEXT_COLOR = curses.COLOR_YELLOW  # BLACK
COL_HEADER_BG_COLOR = STANDARD_BG_COLOR  # curses.COLOR_BLUE #BLACK #YELLOW

ROW_HEADER_ID = 4
ROW_HEADER_TEXT_COLOR = curses.COLOR_GREEN  # BLACK
ROW_HEADER_BG_COLOR = STANDARD_BG_COLOR  # curses.COLOR_BLUE #BLACK #GREEN

BOX_COLOR_ID = 5  # The borders of the cells
BOX_COLOR_TEXT_COLOR = curses.COLOR_YELLOW
BOX_COLOR_BG_COLOR = STANDARD_BG_COLOR

BLANK_COLOR_ID = 6
BLANK_COLOR_TEXT_COLOR = curses.COLOR_CYAN
BLANK_COLOR_BG_COLOR = STANDARD_BG_COLOR

SEARCH_MATCH_COLOR_ID = 7  # Highlighted search results
SEARCH_MATCH_COLOR_TEXT_COLOR = curses.COLOR_YELLOW
SEARCH_MATCH_COLOR_BG_COLOR = curses.COLOR_RED

WARNING_COLOR_ID = 8  # "Error message" color
WARNING_COLOR_TEXT_COLOR = curses.COLOR_YELLOW
WARNING_COLOR_BG_COLOR = curses.COLOR_RED

NUMERIC_NEGATIVE_COLOR_ID = 9  # Negative numbers are this style
NUMERIC_NEGATIVE_COLOR_TEXT_COLOR = curses.COLOR_RED
NUMERIC_NEGATIVE_COLOR_BG_COLOR = STANDARD_BG_COLOR

NUMERIC_POSITIVE_COLOR_ID = 10  # Positive numbers are this style
NUMERIC_POSITIVE_COLOR_TEXT_COLOR = curses.COLOR_CYAN
NUMERIC_POSITIVE_COLOR_BG_COLOR = STANDARD_BG_COLOR

ACTIVE_FILENAME_COLOR_ID = 11
ACTIVE_FILENAME_COLOR_TEXT_COLOR = curses.COLOR_YELLOW
ACTIVE_FILENAME_COLOR_BG_COLOR = STANDARD_BG_COLOR

HELP_AREA_ID = 12
HELP_AREA_TEXT_COLOR = curses.COLOR_RED
HELP_AREA_BG_COLOR = curses.COLOR_BLUE

LAST_USED_PROPERTY_ID = HELP_AREA_ID

# dictionary
COLOR_PROPERTIES = {"NA": {"fg": curses.COLOR_BLUE, "bg": STANDARD_BG_COLOR}}

# NONE_COLOR = curses.COLOR_YELLOW
# NONE_BKG = curses.COLOR_YELLOW

# RAGGED_END_ID = 1
# RAGGED_END_TEXT_COLOR   = curses.COLOR_WHITE
# RAGGED_END_BG_COLOR     = curses.COLOR_BLUE

# SELECTED_CELL_ID = 2
# SELECTED_CELL_TEXT_COLOR = curses.COLOR_BLACK
# SELECTED_CELL_BG_COLOR = curses.COLOR_MAGENTA

# COL_HEADER_ID = 3
# COL_HEADER_TEXT_COLOR = curses.COLOR_BLACK
# COL_HEADER_BG_COLOR = curses.COLOR_YELLOW

# ROW_HEADER_ID = 4
# ROW_HEADER_TEXT_COLOR = curses.COLOR_BLACK
# ROW_HEADER_BG_COLOR = curses.COLOR_GREEN


GlobalCurrentFile = None
GlobalCurrentFilename = None
GlobalCurrentNumLinesLoaded = 0

DEFAULT_HIGHLIGHT_NUMBERS_SETTING = True  # By default, color numbers + as green and - as red. (Based on the settings in NUMERIC_POSITIVE_COLOR_ID and NUMERIC_NEGATIVE_COLOR_ID.)

_DEBUG = False  # Run the program with "-w" for debugging to be enabled

_debugFilename = "DEBUG.sheet.py.out.tmp"


def DebugPrint(argStr="", nl="\n"):
    f = open(_debugFilename, "a")
    f.write(str(argStr) + nl)
    f.close()
    print(str(argStr) + nl)
    return


def DebugPrintTable(table):
    for entireRow in table:
        for cell in entireRow:
            DebugPrint(stringFromAny(cell) + ", ", nl="")
            pass
        DebugPrint()
        pass
    return


class Point:
    def __init__(self, argX, argY):
        self.x = argX
        self.y = argY


# class RowCol: # Like a Point, but ROW (Y) comes first, then COL (X)
#     def __init__(self, argRow, argCol):
#         self.row = argRow
#         self.col = argCol

#     def __init__(self, rowColList): # Used to init size from "getmaxyx()" apparently python doesn't like multiple constructors!
#         (self.row, self.col) = rowColList


class Size:
    def __init__(self, argWidth, argHeight):
        self.width = argWidth
        self.height = argHeight

    def resizeYX(self, rowColList):
        (self.height, self.width) = rowColList


class AGW_File_Data_Collection:
    """
    This object stores all the state about a file. It's a lot like a "window" on a normal GUI. Confusingly, it is NOT the same as an AGW_Win, which is a "CURSES" terminal sub-window. Sorry for the confusion.
    """

    def __init__(self):
        self.collection = []  # <-- contants: a bunch of AGW_File_Data objects
        self.currentFileIdx = None  # which one are we currently looking at? probably should be not part of the object, come to think of it.
        pass

    def size(self):
        """Basically, how many files did we load."""
        return len(self.collection)

    def getCurrent(self):  # AGW_File_Data_Collection
        """Get the current data object that backs the spreadsheet we are looking at right this second."""
        return self.getInfoAtIndex(self.currentFileIdx)

    def getInfoAtIndex(self, i):
        return self.collection[i]

    def addFileInfo(self, fileInfoObj):
        if not isinstance(fileInfoObj, AGW_File_Data):
            print("Uh oh, someone tried to add some random object into the file info collection.")
            raise
        else:
            self.collection.append(fileInfoObj)
            pass


class AGW_File_Data:
    def __init__(self, argFilename):
        self.filename = argFilename
        self.table = AGW_Table()
        self.defaultCellProperty = curses.A_NORMAL  # REVERSE
        self.hasColHeader = False
        self.hasRowHeader = False
        self.cursorPos = Point(0, 0)  # what is the selected cell
        self.__regex = None
        self.__compiledRegex = None
        self.regexIsCaseSensitive = False
        self.boolHighlightNumbers = DEFAULT_HIGHLIGHT_NUMBERS_SETTING
        pass

    def getNumCols(self):
        return self.table.getNumCols()

    def getNumRows(self):
        return self.table.getNumRows()

    def getActiveCellX(self):
        return self.cursorPos.x

    def getActiveCellY(self):
        return self.cursorPos.y

    def toggleNumericHighlighting(self):
        self.boolHighlightNumbers = not self.boolHighlightNumbers
        setCommandStr("Toggled highlighting of numeric values.")
        return

    def getRegexString(self):
        return stringFromAny(self.__regex)

    def appendToCurrentSearchTerm(self, newThing):
        self.changeCurrentSearchTerm(stringFromAny(self.__regex) + stringFromAny(newThing))
        return

    # In class AGW_File_Data
    def changeCurrentSearchTerm(self, argSearchString, argIsCaseSens=None):
        self.__regex = argSearchString
        # setWarning("just set regex to: " + str(self.__regex))

        if argIsCaseSens is not None:
            self.regexIsCaseSensitive = argIsCaseSens
            pass

        if lenFromAny(self.__regex) <= 0:
            self.__regex = None
            self.__compiledRegex = None
            pass
        else:
            if self.regexIsCaseSensitive:
                self.__compiledRegex = re.compile(self.__regex)
                pass
            else:
                self.__compiledRegex = re.compile(self.__regex, re.IGNORECASE)
                pass
            pass
        return

    def clearCurrentSearchTerm(self):
        self.changeCurrentSearchTerm(None, None)
        pass

    def trimRegex(self, numChars=1):
        """Remove the last <numChars> characters from the end of the search term (i.e., it is like pressing backspace). Clears the search term (which sets it to None) if it is going to be zero-length."""
        currentRegexLength = lenFromAny(self.__regex)
        if currentRegexLength <= 1:
            self.clearCurrentSearchTerm()
        else:
            self.changeCurrentSearchTerm(self.__regex[: (currentRegexLength - numChars)])
        return

    def regexIsActive(self):
        """Tells us whether we should be highlighting the search terms or not"""
        return self.__regex is not None

    # In class AGW_File_Data
    def stringDoesMatchRegex(self, stringToCheck):
        if self.__regex is None or stringToCheck is None:
            return False

        if self.__compiledRegex.search(
            stringToCheck
        ):  # Note the difference between *search* and *match*(match does not do partial results)
            # DebugPrint("Got a match with: " + stringToCheck + " from the compiled regex: " + self.__regex)
            return True
        else:
            # DebugPrint("NO match for: " + stringToCheck + " from the compiled regex: " + self.__regex)
            return False


# =======================================
# End of class AGW_File_Data
# =======================================

# =======================================
# Start of class AGW_Win
# These are "curses" (terminal) "windows"--basically just regions in the terminal
# that can be written to separately with their own coordinate systems.
# =======================================
class AGW_Win:
    def __init__(self):
        self.win = None
        self.pos = Point(0, 0)  # where is the top-left of this window?
        self.windowWidth = 0
        self.windowHeight = 0
        pass

    # Make a new window with the height, width, and at the location specified
    def initWindow(self, argHeight, argWidth, atY, atX):
        self.pos = Point(atX, atY)
        self.windowWidth = argWidth
        self.windowHeight = argHeight
        try:
            self.win = curses.newwin(argHeight, argWidth, atY, atX)
            pass
        except:
            raise "Cannot allocate a curses window with height " + str(argHeight) + " and width " + str(
                argWidth
            ) + " at Y=" + str(atY) + " and X=" + str(
                atX
            )  # + ". Error was: " + err.message
        pass

    # safeAddCh: Safely adds a single character to an AGW_Win object
    # It is "safe" because it does not throw an error if it overruns
    # the pad (instead, it just doesn't draw anything at all)
    def safeAddCh(self, y, x, argChar, attr=0):
        if y >= self.windowHeight or x >= self.windowWidth - 1:  # out of bounds!
            return

        try:
            self.win.addch(y, x, argChar, attr)
            pass
        except (curses.error):
            # 1, "safeAddCh is messed up: " + err.message)
            print(
                "ERROR: Unable to print a character at",
                y,
                x,
                "with window dimensions (in chars): ",
                self.windowHeight,
                self.windowWidth,
            )
            raise
        except:
            raise

    # safeAddStr: Safely adds a string to a CURSES "win" object.
    # It is "safe" because it does not throw an error if it overruns
    # Additionally, it does NOT WRAP TEXT. This is different from default addstr.
    def safeAddStr(self, y, x, string, attr=0):
        try:
            if string is None:
                return
            if x < 0 or y < 0:
                return

            if y >= self.windowHeight:
                return  # off screen!

            if x + len(string) >= self.windowWidth:  # Runs off to the right.
                newLength = max(0, (self.windowWidth - x - 1))
                string = string[:newLength]
                pass

            if len(string) > 0:
                self.win.addstr(y, x, string, attr)

        except (curses.error):
            print("safeAddStr is messed up!")
            raise
        except:
            raise


## Data window (the main window with the cells in it)
class AGW_DataWin(AGW_Win):
    def __init__(self):
        AGW_Win.__init__(self)  # parent constructor
        self.info = None
        self.defaultCellProperty = curses.A_NORMAL  # REVERSE
        pass

    def getTable(self):
        return self.info.table

    def setInfo(self, whichInfo):
        if not isinstance(whichInfo, AGW_File_Data):
            print("### Someone passed in a not-an-AGW_File_Data object to AGW_DataWin--->setInfo()\n")
            raise
        self.info = whichInfo
        return

    def getInfo(self):
        return self.info

    def drawTable(
        self,
        whichInfo,
        topCell,
        leftCell,
        nRowsToDraw=None,
        nColsToDraw=None,
        boolPrependRowCoordinate=False,
        boolPrependColCoordinate=False,
    ):
        self.win.erase()  # or clear()

        theTable = self.getTable()

        if nColsToDraw is not None:
            numToDrawX = nColsToDraw
        else:
            numToDrawX = theTable.getNumCols()

        if nRowsToDraw is not None:
            numToDrawY = nRowsToDraw
        else:
            numToDrawY = theTable.getNumRows()

        start = Point(leftCell, topCell)  # "which cells to draw"
        end = Point(
            min(leftCell + numToDrawX, theTable.getNumCols()),
            min(topCell + numToDrawY, theTable.getNumRows()),
        )

        if end.y >= theTable.getNumRows():  # See if we need to read more lines from the file
            theTable.readFromCurrentFile()  # Ok, read some more lines!

        cellTextPos = Point(None, None)
        for r in range(start.y, end.y):  # Go through each row, one at a time, top to bottom
            cellTextPos.y = gCellBorders.height + (gCellBorders.height + cellHeight) * (
                r - start.y
            )  # <-- all cells are the same HEIGHT, so this can be computed in one equation

            if cellTextPos.y >= self.windowHeight:
                break

            cellTextPos.x = gCellBorders.width  # initialization! (update is way down below)

            for c in range(start.x, end.x):  # Go through ALL the cells on this row
                if cellTextPos.x >= self.windowWidth:
                    break

                cell = theTable.cellValue(r, c)  # <-- this is the actual text that is in the cell
                maxLenForThisCell = theTable.getColWidth(c)  # This is how long this cell text can be.

                if boolPrependColCoordinate:
                    cell = str(c + 1) + HEADER_NUM_DELIMITER_STRING + stringFromAny(cell)
                    pass

                if boolPrependRowCoordinate:
                    cell = str(r + 1) + HEADER_NUM_DELIMITER_STRING + stringFromAny(cell)
                    pass

                selectedPt = whichInfo.cursorPos  # Where is the cursor...

                cellIsSelected = c == selectedPt.x and r == selectedPt.y
                shouldHighlightCell = cellIsSelected  # or (self.highlightEntireCol and c == selectedPt.x) or (self.highlightEntireRow and r == selectedPt.y)

                cellAttr = self.defaultCellProperty

                if whichInfo.boolHighlightNumbers:
                    cellAttr = attributeForNumeric(cell, self.defaultCellProperty)

                if cell in ("NA", "na", "N/A", "n/a", "NaN"):
                    cellAttr = curses.color_pair(COLOR_PROPERTIES["NA"]["id"])
                    pass

                if cell is None:
                    # This doesn't work for some reason, maybe it's getting overwritten?
                    # Note: this is NOT exactly the same as a cell with no data in it! It's a 'ragged end' situation
                    cell = padStrToLength("", maxLenForThisCell, "~")  #  distinct from an *empty* cell)
                    cellAttr = curses.color_pair(
                        WARNING_COLOR_ID
                    )  # (indicate that there isn't a cell here at all
                    pass

                drawCheckerboard = False
                if r == 0 and c == 0 and (self.info.hasRowHeader and self.info.hasColHeader):
                    cellAttr = curses.A_NORMAL  # there is an "odd man out" in the top left, for files
                    drawCheckerboard = True  # with both a column AND a row header
                    pass
                elif r == 0 and self.info.hasColHeader:
                    cellAttr = curses.color_pair(COL_HEADER_ID)
                    drawCheckerboard = True
                    pass
                elif c == 0 and self.info.hasRowHeader:
                    cellAttr = curses.color_pair(ROW_HEADER_ID)
                    drawCheckerboard = True
                    pass

                if shouldHighlightCell:
                    cellAttr = curses.color_pair(SELECTED_CELL_ID) + curses.A_NORMAL
                    pass

                # setWarning(str(self.getColWidth(c)) + " is the length for col " + str(c))

                cell = padStrToLength(
                    cell, maxLenForThisCell, " "
                )  # <-- needed so that the cells don't randomly get garbage everywhere when they don't re-draw all the way. Very annoying. If you remove this, then long cells will randomly overwrite text of shorter cells when you move left and right.
                # cell = truncateLongCell(cell, maxLenForThisCell, truncationSuffix)

                if drawCheckerboard:
                    bgAttr = cellAttr
                    self.win.attron(bgAttr)  # curses.color_pair(BOX_COLOR_ID))
                    hLineLength = max(0, min(self.windowWidth - cellTextPos.x, maxLenForThisCell))
                    self.win.hline(cellTextPos.y, cellTextPos.x, " ", hLineLength)  # checkerboard character
                    self.win.attroff(bgAttr)  # curses.color_pair(BOX_COLOR_ID))
                    pass
                # curses.ACS_HLINE
                # curses.ACS_DIAMOND
                # curses.ACS_CKBOARD # checkerboard

                # whichInfo.changeCurrentSearchTerm("F")

                if (
                    whichInfo.regexIsActive() and self is sheetWin
                ):  # <<<<<<< HORRIBLE HACK!!!! FIX LATER!!! should work on all windows!! not just the main one
                    # "self is sheetWin" is a horrible hack! Fix it eventually
                    theTable.initRegexTable()
                    if theTable.regTab.isDirty(r, c):
                        # calcluate the regex...
                        # DebugPrint("calc: " + str(r) + ", " + str(c))
                        boolRegexMatched = whichInfo.stringDoesMatchRegex(cell)
                        # DebugPrint("  match status is: " + str(boolRegexMatched))
                        theTable.regTab.set(row=r, col=c, value=boolRegexMatched)
                        pass

                    if theTable.regTab.matches(row=r, col=c):
                        # DebugPrint("This was the result at " + str(r) + ", " + str(c) + ": " + whichInfo.regex + " <-- " + cell)
                        cellAttr = curses.color_pair(
                            SEARCH_MATCH_COLOR_ID
                        )  # + curses.A_REVERSE #curses.A_BLINK
                        pass

                    pass

                self.safeAddStr(cellTextPos.y, cellTextPos.x, cell, cellAttr)  # Draw text
                self.win.attron(curses.color_pair(BOX_COLOR_ID))

                if gCellBorders.height > 0:  # horizontal lines
                    hLineLength = max(0, min(self.windowWidth - cellTextPos.x - 1, maxLenForThisCell))

                    try:
                        self.win.hline(cellTextPos.y - 1, cellTextPos.x, curses.ACS_HLINE, hLineLength)
                        pass
                    except:
                        print(
                            "problem when cellTextPos is y="
                            + str(cellTextPos.y)
                            + ", x="
                            + str(cellTextPos.x)
                            + " and also the rows and cols are "
                            + str(r)
                            + " and "
                            + str(c)
                            + " and terminal width is "
                            + str(gTermSize.width)
                            + " and height is "
                            + str(gTermSize.height)
                        )
                        raise

                    pass

                if gCellBorders.width > 0:  # vertical lines
                    vLineLength = max(0, min(self.windowHeight - cellTextPos.y, gCellBorders.width))
                    try:
                        self.win.vline(cellTextPos.y, cellTextPos.x - 1, curses.ACS_VLINE, vLineLength)
                        pass
                    except:
                        print(
                            "problem when cellTextPos is y="
                            + str(cellTextPos.y)
                            + ", x="
                            + str(cellTextPos.x)
                            + " and also the rows and cols are "
                            + str(r)
                            + " and "
                            + str(c)
                            + " and terminal width is "
                            + str(gTermSize.width)
                            + " and height is "
                            + str(gTermSize.height)
                        )
                        raise
                    pass

                if gCellBorders.width > 0 and gCellBorders.height > 0:
                    ch = calculateBorderChar(
                        r=r,
                        c=c,
                        topRow=0,
                        leftCol=0,
                        bottomRow=theTable.getNumRows(),
                        rightCol=theTable.getNumCols(),
                    )
                    self.safeAddCh(cellTextPos.y - 1, cellTextPos.x - 1, ch)
                    pass

                self.win.attroff(curses.color_pair(BOX_COLOR_ID))

                cellTextPos.x += gCellBorders.width + maxLenForThisCell
                pass
            pass

        self.win.refresh()
        pass  # end of "drawTable"
        # theWin.attroff(curses.A_REVERSE)

    #             Attributes: A_BLINK, A_BOLD, A_DIM, A_REVERSE, A_STANDOUT, A_UNDERLINEUnderlined


class AGW_Table:
    def __init__(self):
        self.clearTable()
        pass

    def clearTable(self):
        self.__nCells = Size(0, 0)  # size in cells
        self.__cells = []
        self.colWidth = []  # maximum char length in a col
        self.isRagged = False  # Does the table have "ragged" ends? (differing col counts). Ragged means "some rows have more cols than others"
        self.regTab = None  # this will be a table of "None" (not yet checked) "True" or "False", depending on whether the cell currently matches the regex!
        pass

    def getColWidth(self, colIdx):
        # Report how wide each column is
        if gIsTransposed:
            return WIDTH_WHEN_TRANSPOSED  # Transposed cells are always the same width. This is kind of a hack, as we just don't compute the cell widths
        try:
            if self.colWidth[colIdx] > MAX_COL_WIDTH:
                return MAX_COL_WIDTH
            else:
                return self.colWidth[colIdx]
        except:
            print(
                "### Someone passed in an invalid column index, "
                + str(colIdx)
                + ". Max was "
                + str(self.getNumCols())
                + ".\n"
            )
            raise  # return 0 #raise #return 0

    def getNumCols(self):
        if gIsTransposed:
            return self.__nCells.height  # <-- if displaying transposed!
        else:
            return self.__nCells.width  # table width in number of cells

    def getNumRows(self):
        if gIsTransposed:
            return self.__nCells.width  # <-- if displaying transposed!
        else:
            return self.__nCells.height  # table height in number of cells

    def getHeaderCellForCol(self, colIdx):
        return self.cellValue(0, colIdx)

    def getHeaderCellForRow(self, rowIdx):
        return self.cellValue(rowIdx, 0)

    def cellValue(self, row, col):
        try:
            if gIsTransposed:
                return self.__cells[col][row]  # If it's transposed, then flip the row and column!
            else:
                return self.__cells[row][col]
        except:
            return "~~~~~"  # ("R" + str(row) + ", C" + str(col) + " out of bounds")

    def initRegexTable(self):
        self.regTab = SearchBoard(nrow=self.getNumRows(), ncol=self.getNumCols())
        return

    def appendRowOfCellContents(self, contents):
        # Turns a line from a file into a properly-formatted internal
        # representation of row of cells.

        self.__cells.append(contents)

        if self.__nCells.width > 0 and self.__nCells.width != len(contents):
            self.isRagged = True  # Ragged table! (not a "true" table)
            pass

        numItemsInThisNewRow = len(contents)
        self.__nCells.width = max(
            self.__nCells.width, numItemsInThisNewRow
        )  # Set the table width to that of the MAXIMUM number of cols that a row has
        self.__nCells.height += 1  # Read another row...

        for c in range(0, self.__nCells.width):
            # Ok, find the **longest** item in each column. This could be slow for a tall column!
            cellWidth = None
            if c >= numItemsInThisNewRow or (contents[c] is None):
                cellWidth = 0
            else:
                cellWidth = len(contents[c])

            numColsWithWidth = len(self.colWidth)
            if c == 0 and c < numColsWithWidth:  # and self.hasRowHeader):
                # It's the ROW header (the leftmost column)! Gotta account for the ": " in the row header
                cellWidth += len(str(self.getNumCols())) + len(
                    HEADER_NUM_DELIMITER_STRING
                )  # + gCellBorders.width*2
                pass

            if c >= numColsWithWidth:
                cellWidth += len(str(c + 1)) + len(
                    HEADER_NUM_DELIMITER_STRING
                )  # Accounting for the col header!
                self.colWidth.append(cellWidth)
                pass
            else:
                self.colWidth[c] = max(self.colWidth[c], cellWidth)
                pass
            pass
        pass

    def loadNewFile(self, theFilename):
        self.closeCurrentFile()  # Close the old table first, if there is one
        global GlobalCurrentNumLinesLoaded
        global GlobalCurrentFile
        global GlobalCurrentFilename
        self.clearTable()
        # We want to see which file we should read from.
        # 1. We either got a HYPHEN, which is the UNIX shorthand for "read from STDIN"
        # 2. or we got a real filename, in which case we try to read the file from the filesystem
        if theFilename is kSTDIN_HYPHEN:
            GlobalCurrentFile = sys.stdin
        elif re.search(r"[.](Z|gz|gzip)$", theFilename, re.IGNORECASE):
            GlobalCurrentFile = gzip.open(theFilename)  # read GZIP
        elif re.search(r"[.](bz2|bzip2)$", theFilename, re.IGNORECASE):
            GlobalCurrentFile = bz2.BZ2File(theFilename)  # read BZIP2
        elif re.search(r"[.](zip)$", theFilename, re.IGNORECASE):
            raise "Sadly we cannot look at a .ZIP file directly!"
        else:
            GlobalCurrentFile = open(theFilename, "r")  # read NORMAL TEXT
        print("Loaded the file named " + theFilename)
        GlobalCurrentFilename = theFilename
        GlobalCurrentNumLinesLoaded = 0
        pass

    def closeCurrentFile(self):
        global GlobalCurrentNumLinesLoaded
        global GlobalCurrentFile
        global GlobalCurrentFilename
        if (
            GlobalCurrentFilename is kSTDIN_HYPHEN
        ):  # If we were reading from STDIN, then we have some cleanup to do.
            os.close(0)  # <-- closes stdin (required!)
            sys.stdin = open("/dev/tty", "r")  # <-- more STDIN cleanup...
            pass
        else:
            if GlobalCurrentFile is not None:
                GlobalCurrentFile.close()
            pass
        GlobalCurrentFile = None
        GlobalCurrentFilename = None
        GlobalCurrentNumLinesLoaded = 0
        pass

    def readFromCurrentFile(self):
        global GlobalCurrentNumLinesLoaded
        global GlobalCurrentFile
        global GlobalCurrentFilename

        local_delim = (
            global_delimiter if global_delimiter else guess_delimiter_from_filename(GlobalCurrentFilename)
        )

        try:
            # Read lines from a file into our data structure
            # probably should use xreadlines!
            lines = GlobalCurrentFile.readlines(
                kNUM_BYTES_TO_READ_AT_A_TIME
            )  # <-- note that we only read a small number of lines every time the user interacts with the file!
            if lines:
                for line in lines:
                    line = line.rstrip("\n")  # <-- like Perl "chomp"--remove the trailing newlines
                    self.appendRowOfCellContents(line.split(local_delim))
                    GlobalCurrentNumLinesLoaded += 1
                    pass
                pass
        except:
            print(
                "[sheet.py]: Error: Problem attempting to read from the file named " + GlobalCurrentFilename
            )
            print(
                "[sheet.py]:        We were able to read "
                + str(GlobalCurrentNumLinesLoaded)
                + " lines from it."
            )
            print("[sheet.py]:        Verify that this file exists and is readable.")
            raise

        return  # end of function


def lenFromAny(argThing):  # returns 0 for None's length. Otherwise you get an exception
    if argThing is None:
        return 0
    else:
        return len(argThing)


def stringFromAny(argString):  # Returns a string, even if given None as an input
    if argString is None:
        return ""
    else:
        return str(argString)


def cursesClearLine(window, lineYPos):
    window.move(lineYPos, 0)
    window.clrtoeol()
    pass


def agwEnglishPlural(string, numOf, suffix="s"):
    """You pass in a string like "squid" and a number
    indicating how many squid there are. If the number
    is one, then "squid" is returned, otherwise "squids"
    is returned."""
    if 1 == numOf:
        return string
    else:
        return string + suffix
    pass


def attributeForNumeric(theProspectiveNumber, defaultAttr):
    # Give it something that might be a number.
    # if it *is* a numnber, it returns the format for that type of
    # number. If it isn't, it returns "defaultAttr"
    try:
        if float(theProspectiveNumber) < 0:
            return curses.color_pair(NUMERIC_NEGATIVE_COLOR_ID)
        elif float(theProspectiveNumber) > 0:
            return curses.color_pair(NUMERIC_POSITIVE_COLOR_ID)
        else:  # for zero, maybe NaN too
            return defaultAttr
    except (ValueError, TypeError):
        return defaultAttr  # not a number
    pass


gTermSize = Size(None, None)  # Size of the actual xterm terminal window that sheet.py lives in

gStandardScreen = None  # The main screen

gIsTransposed = False  # Whether to show the data "as-is" or transposed

mainInfo = AGW_File_Data_Collection()  # An array of AGW_File_Datas

sheetWin = AGW_DataWin()
infoWin = AGW_Win()
helpWin = AGW_Win()

colHeaderWin = AGW_DataWin()  # "Window" for the column headers at the top of the screen
rowHeaderWin = AGW_DataWin()  # "Window" for the row headers at the left side of the screen

cellHeight = 1  # How many lines of text each cell takes up, vertically. 1 is normal. Integer.

gCellBorders = Size(1, 0)  # Size(1,1) # Width, Height
# Size(1, 0) means there is a border to the left/right of each cell, but no border above/below

gCurrentMode = KEY_MODE_NORMAL_INPUT  # Start in normal input (not search) mode

gWantToQuit = False  # False, since, by default, the user does not want to immediately quit!

# windowPos = Point(0,0) # Initial position in the table of data: Point(0,0) is the top left

truncationSuffix = "..."  # Suffix for "this cell is too long to fit on screen"

global_delimiter = None  # automatically guess with each file "\t"

fastMoveSpeed = Point(10, 10)  # How many cells to scroll when the user is moving with shift-arrows

gWarningMessage = None  # A message that we can display in the "warning" region in the info pane
gCommandStr = None  # The current command, which we sometimes display in the info pane


def setWarning(string):
    global gWarningMessage
    gWarningMessage = string
    return


def setCommandStr(string):
    global gCommandStr
    gCommandStr = string


def clearCommandStr():
    global gCommandStr
    gCommandStr = None


def usageAndQuit(exitCode, message=None):
    fillWidth = 80
    message = textwrap.fill(message, fillWidth)
    if message is not None:
        print("")
        print("sheet.py: ")
        print(message)
        print("sheet.py: Printing usage information below.")
        print("*" * fillWidth)
        print("*" * fillWidth)
        pass
    print(ALEX_PROGRAM_USAGE_TEXT)  # at the very bottom of this file
    if message is not None:
        print("(End of usage information for sheet.py.)")
        print("*" * fillWidth)
        print("*" * fillWidth)
        print("sheet.py: ")
        print(message)
        pass
    sys.exit(exitCode)
    return


def initializeWindowSettings(aScreen, fileInfoToReadFrom):
    # if (fileInfoToReadFrom is None):
    #    usageAndQuit(1, "Missing a command-line argument: We did not have at least one file passed in as an argument on the command line! Pass in at least one valid file on the command line. Maybe you passed in a file that could not be read for some reason, or passed in a directory name.\n")
    #    raise
    if not isinstance(fileInfoToReadFrom, AGW_File_Data):
        print(
            "### Init window: Someone passed in a not-an-AGW_File_Data object to initializeWindowSettings\n"
        )
        raise

    global sheetWin, infoWin, colHeaderWin, rowHeaderWin, helpWin

    # =========== Specify the dimensions of the windows here! ==================
    INFO_PANEL_HEIGHT = 6  # How many vertical lines of terminal is the info panel in total?
    COL_HEADER_HEIGHT = 2  # How many vertical lines of terminal is the column header? Answer, it should be 2: one line for the header itself, and one below for the horizontal line that separates it from the main window (sheetWin).
    HELP_WIN_HEIGHT = 2  # How tall is the "help" window at the bottom?
    # =========== Specify the dimensions of the windows here! ==================

    gTermSize.resizeYX(aScreen.getmaxyx())  # Resize the curses window to use the whole available screen
    sheetWin.setInfo(fileInfoToReadFrom)
    fileInfoToReadFrom.table.loadNewFile(fileInfoToReadFrom.filename)
    fileInfoToReadFrom.table.readFromCurrentFile()

    if fileInfoToReadFrom.hasColHeader and fileInfoToReadFrom.getNumRows() >= 2:
        fileInfoToReadFrom.cursorPos.y = 1  # start off with the column to the RIGHT of the col header selected, instead of having the column header show up twice
        pass

    if fileInfoToReadFrom.hasRowHeader and fileInfoToReadFrom.getNumCols() >= 2:
        fileInfoToReadFrom.cursorPos.x = 1  # start off with the first not-a-header row highlighted
        pass

    # ===============================================
    # Figure out how wide the row header should be.
    rowHeaderMaxWidth = 0
    if sheetWin.getTable().getNumCols() > 0:
        # If there is at least one column, then the row header is as wide as the first column
        # (but note that the maximum size is adjusted below)
        rowHeaderMaxWidth = sheetWin.getTable().getColWidth(0) + (gCellBorders.width * 2)
        pass

    # However: the maximum allowed rowHeaderWidth is some fraction of the total screen width
    if rowHeaderMaxWidth > (gTermSize.width * ROW_HEADER_MAXIMUM_COLUMN_FRACTION_OF_SCREEN):
        rowHeaderMaxWidth = int(gTermSize.width * ROW_HEADER_MAXIMUM_COLUMN_FRACTION_OF_SCREEN)
        pass
    # Done figuring out how wide the row header should be
    # ===============================================

    SHEET_WIN_HEIGHT = (
        gTermSize.height - INFO_PANEL_HEIGHT - COL_HEADER_HEIGHT - HELP_WIN_HEIGHT
    )  # How tall is the main data window?
    sheetWin.initWindow(
        SHEET_WIN_HEIGHT,  # height
        (gTermSize.width - rowHeaderMaxWidth),  # width
        (INFO_PANEL_HEIGHT + COL_HEADER_HEIGHT),  # location (y)
        rowHeaderMaxWidth,
    )  # location (x)

    infoWin.initWindow(
        INFO_PANEL_HEIGHT, gTermSize.width, 0, 0  # height  # width  # location--y
    )  # location--x

    helpWin.initWindow(
        HELP_WIN_HEIGHT,  # height
        gTermSize.width,  # width
        (INFO_PANEL_HEIGHT + COL_HEADER_HEIGHT + SHEET_WIN_HEIGHT),  # location (y)
        0,
    )  # location (x)

    # Goes along the left side!
    ROW_HEADER_HEIGHT = SHEET_WIN_HEIGHT
    ROW_HEADER_LOC_X = 0
    ROW_HEADER_LOC_Y = INFO_PANEL_HEIGHT + COL_HEADER_HEIGHT
    rowHeaderWin.initWindow(ROW_HEADER_HEIGHT, rowHeaderMaxWidth, ROW_HEADER_LOC_Y, ROW_HEADER_LOC_X)
    # rowHeaderWin.defaultCellProperty = curses.color_pair(ROW_HEADER_ID)

    # Goes along the TOP
    COL_HEADER_WIDTH = gTermSize.width - rowHeaderMaxWidth
    COL_HEADER_LOC_X = rowHeaderMaxWidth
    COL_HEADER_LOC_Y = INFO_PANEL_HEIGHT
    colHeaderWin.initWindow(COL_HEADER_HEIGHT, COL_HEADER_WIDTH, COL_HEADER_LOC_Y, COL_HEADER_LOC_X)
    # colHeaderWin.defaultCellProperty = curses.color_pair(COL_HEADER_ID)

    colHeaderWin.setInfo(fileInfoToReadFrom)
    rowHeaderWin.setInfo(fileInfoToReadFrom)

    fastMoveSpeed.x = 5  # gTermSize.width // 2 # <-- measured in CELLS, not characters!
    fastMoveSpeed.y = sheetWin.windowHeight // 2  # measured in CELLS, not characters!
    pass


def guess_delimiter_from_filename(filename):
    if re.search("[.]csv([.](bz2|bzip|gz|gzip|xz|Z))?", filename, flags=re.IGNORECASE):
        return ","  # comma delimiter
    if re.search("[.](tsv|tab|txt)([.](bz2|bzip|gz|gzip|xz|Z))?", filename, flags=re.IGNORECASE):
        return "\t"  # tab delimited
    return "\t"  # default delimiter is tab


def processCommandLineArgs(argv):
    try:
        opts, args = getopt.gnu_getopt(
            argv, "hwd:i:", ["help", "warn", "delim=", "input="]
        )  # Docs for getopt: http://docs.python.org/library/getopt.html
    except (getopt.GetoptError):
        # usageAndQuit(1, "Encountered an unknown command line option!\n" + "sheet.py: Please remember that long names must have double-dashes!\n" + "sheet.py: (i.e. -warn generates an error, but --warn is correct)\n")
        raise

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usageAndQuit(0, "Printed the HELP information, since the --help option was supplied.")
            pass
        elif opt in ("-w", "--warn"):
            global _DEBUG
            _DEBUG = True
            pass
        elif opt in ("-d", "--delim"):
            global global_delimiter
            if len(arg) != 1:
                usageAndQuit(
                    1,
                    "Your delimiter, which is between the parentheses here --> ("
                    + arg
                    + ") must be EXACTLY ONE CHARACTER. Yours was not, which is not allowed.",
                )
                pass
            global_delimiter = arg
            pass
        else:
            pass

    if _DEBUG:
        print("Unprocessed arguments:", args)
        print("Unprocessed options:", opts)
        pass

    arrFilenamesToRead = []

    if sys.stdin.isatty() is False:
        # If the user has piped in a file, add that to the "things to read"
        # Note that STDIN must be the VERY FIRST thing in the list of filenames
        # to read.
        if len(arrFilenamesToRead) > 0:
            setWarning(
                "Since input was piped in through STDIN, we ignored the files that were passed on the command line."
            )
            pass
        arrFilenamesToRead.extend(kSTDIN_HYPHEN)  # Add to the list
        pass
    else:
        arrFilenamesToRead.extend(args)  # Add the unread files... AFTER handling stdin.
        pass

    # print("Files to read = <" , arrFilenamesToRead , ">") ; sys.exit(1)
    if not arrFilenamesToRead:
        usageAndQuit(
            0, "ARGUMENT ERROR: You must specify at least one filename (to read from) on the command line!"
        )
        pass

    for filename in arrFilenamesToRead:
        # Check to make sure the specified filenames are valid.
        # We also allow a hyphen, which is UNIX shorthand for "read from STDIN."
        if filename is kSTDIN_HYPHEN or os.path.isfile(filename):
            # If we got here, then it's a valid file or STDIN to read from...
            singleFileInfo = AGW_File_Data(filename)
            singleFileInfo.hasColHeader = True
            singleFileInfo.hasRowHeader = True
            mainInfo.addFileInfo(singleFileInfo)
            pass
        else:
            print(
                "[sheet.py]: Error in reading a file: Tried to read <"
                + filename
                + ">, but that file either did not exist, or could not be read!"
            )
            pass
        pass

    mainInfo.currentFileIdx = 0
    if mainInfo.size() == 0:
        usageAndQuit(
            0,
            "[sheet.py]: No files that were specified on the command line could be read. Maybe you specified a directory (instead of a list of files). If you want to list all the files in a directory, try:\tsheet.py your_directory/*\n",
        )
    return  # End of command-line-reading function


def drawHelpWin(theScreen):
    helpWin.win.erase()
    lineAttr = curses.color_pair(HELP_AREA_ID)  # <-- set the border color
    helpWin.safeAddStr(
        0,
        0,
        (
            "Help: [Q]uit   [ijkl]: Move cursor (faster with [Shift])"
            + "   [/]: Search             [T]ranspose table                  "
        ),
        curses.color_pair(HELP_AREA_ID),
    )
    helpWin.safeAddStr(
        1,
        0,
        (
            "      [a/e/g/G]: Go to left/right/top/bottom of table   "
            + "   [</>]: Prev/next file   [!]: Toggle number highlighting    "
        ),
        curses.color_pair(HELP_AREA_ID),
    )
    helpWin.win.refresh()
    return


# Example of what the Info Window shows:
# File " " is 2 rows X 3 cols. (Ragged ends)
# File list:
# Row 1: <FDSFSDF>
# Col 8: <LKJOIJWE>
# Value: <Value>
# > Command
def drawInfoWin(theScreen, theInfo, inTab):
    """inTab: the actual table of data that is going to be drawn"""

    activeCellPos = theInfo.cursorPos

    infoWin.win.erase()

    FILE_INFO_ROW = 0
    FILE_LIST_ROW = 1
    ROW_HEADER_ROW = 2
    COL_HEADER_ROW = 3
    VALUE_ROW = 4
    COMMAND_ROW = 5

    if theInfo.filename is kSTDIN_HYPHEN:
        filename = "<STDIN>"
    else:
        filename = '"' + theInfo.filename + '"'

    if gIsTransposed:
        filename = "Transposed file " + filename
    else:
        filename = "File " + filename

    if inTab.isRagged:
        raggedText = " The table is ragged--some rows have differing numbers of columns."
    else:
        raggedText = ""

    if mainInfo.size() >= 1:
        fileNumberStr = " (#" + str(mainInfo.currentFileIdx + 1) + " of " + str(mainInfo.size()) + ")"
    else:
        fileNumberStr = ""

    fileStatusStr = (
        filename
        + fileNumberStr
        + " has "
        + str(inTab.getNumRows())
        + agwEnglishPlural(" row", inTab.getNumRows())
        + " and "
        + str(inTab.getNumCols())
        + agwEnglishPlural(" column", inTab.getNumCols())
        + "."
        + raggedText
    )

    infoWin.safeAddStr(FILE_INFO_ROW, 0, fileStatusStr)  # infoWin is the topmost "window" pane

    rowStr1 = (
        "Row #"
        + str(activeCellPos.y + 1)
        + ": "
        + stringFromAny(inTab.getHeaderCellForRow(activeCellPos.y))
    )
    cursesClearLine(infoWin.win, ROW_HEADER_ROW)
    infoWin.safeAddStr(ROW_HEADER_ROW, 0, rowStr1, curses.color_pair(ROW_HEADER_ID))

    colStr1 = (
        "Col #"
        + str(activeCellPos.x + 1)
        + ": "
        + stringFromAny(inTab.getHeaderCellForCol(activeCellPos.x))
    )
    cursesClearLine(infoWin.win, COL_HEADER_ROW)
    infoWin.safeAddStr(COL_HEADER_ROW, 0, colStr1, curses.color_pair(COL_HEADER_ID))

    cellText = stringFromAny(inTab.cellValue(activeCellPos.y, activeCellPos.x))
    valueStr1 = "(R_" + str(activeCellPos.y + 1) + ", C_" + str(activeCellPos.x + 1) + ") = " + cellText
    infoWin.safeAddStr(VALUE_ROW, 0, valueStr1)

    if mainInfo.size() < 2:
        # If there's only ONE file that the user specified, then don't both showing a file list.
        fileListStr1 = ""  # only one file to read... no need for a list
        pass
    else:
        # Draw a list of the currently files that we can switch between
        xOffset = 0
        fileListStr1 = "File List: "
        infoWin.safeAddStr(FILE_LIST_ROW, 0, fileListStr1)
        xOffset += len(fileListStr1)
        for i in range(0, mainInfo.size()):
            fileListStr1 = mainInfo.getInfoAtIndex(i).filename
            if i == mainInfo.currentFileIdx:
                fileListAttr = curses.color_pair(ACTIVE_FILENAME_COLOR_ID)
                infoWin.safeAddStr(FILE_LIST_ROW, xOffset, fileListStr1, fileListAttr)
                xOffset += len(fileListStr1)
                fileListStr1 = ""
                pass
            fileListAttr = curses.A_NORMAL
            if i < mainInfo.size() - 1:
                fileListStr1 += ", "
                pass
            infoWin.safeAddStr(FILE_LIST_ROW, xOffset, fileListStr1, fileListAttr)
            xOffset += len(fileListStr1)
            pass

        pass

    commandStr1 = "" + stringFromAny(gCommandStr)
    infoWin.safeAddStr(COMMAND_ROW, 0, commandStr1, curses.A_NORMAL)

    #     #setWarning(str(curses.COLORS)) # <-- number of colors the current terminal can support

    if gWarningMessage is not None:
        for i in range(0, 4):
            cursesClearLine(infoWin.win, i)
            attr = None
            if i % 2 == 0:
                attr = curses.color_pair(WARNING_COLOR_ID)
            else:
                attr = curses.color_pair(WARNING_COLOR_ID)
            infoWin.safeAddStr(i, 0, gWarningMessage, attr)
            pass

        setWarning(None)  # And then clear the warning

        pass

    infoWin.win.refresh()

    return


def drawEverything(theScreen):
    activeFI = mainInfo.getCurrent()  # get the current meta-info storage thing for the main window
    activeCellPos = activeFI.cursorPos  # Ask where the cursor is...
    theScreen.refresh()  # Refresh it first...
    drawInfoWin(theScreen, activeFI, sheetWin.getTable())  # The "info" pane at the top...
    drawHelpWin(theScreen)  # The "help" window at the bottom...
    sheetWin.drawTable(activeFI, activeCellPos.y, activeCellPos.x)  # The main table...
    rowHeaderWin.drawTable(
        activeFI, activeCellPos.y, 0, nColsToDraw=1, boolPrependRowCoordinate=True
    )  # Row header... (the leftmost, vertically-oriented pane along the left edge of the table)
    colHeaderWin.drawTable(
        activeFI, 0, activeCellPos.x, nRowsToDraw=1, boolPrependColCoordinate=True
    )  # Column header... (the column header along the top of the table)
    lineAttr = curses.color_pair(BOX_COLOR_ID)  # <-- set the border color
    colHeaderWin.win.attron(lineAttr)  # "start drawing lines"
    colHeaderWin.win.hline(
        colHeaderWin.windowHeight - 1, 0, curses.ACS_CKBOARD, colHeaderWin.windowWidth
    )  # Draw a horizontal line below the column header...
    colHeaderWin.win.attroff(lineAttr)  # "stop drawing lines"
    colHeaderWin.win.refresh()
    return


def mainScreenHandlingLoop(theScreen):
    setUpCurses()
    initializeWindowSettings(theScreen, mainInfo.getCurrent())  # load the first file...
    # theScreen.nodelay(True) # <-- makes "getch" non-blocking
    ch = None
    while not gWantToQuit:
        try:
            drawEverything(theScreen)  # Draw the screen...
            ch = theScreen.getch()  # Now get input from the user...
            # theScreen.nodelay(False)  # <-- makes "getch" *blocking*
            # DebugPrint(str(ch))
            # global GLOB ; GLOB = (GLOB + 1)
            # infoWin.win.addstr(0, 0, "Read char <" + str(ch) + ">" + str(GLOB), curses.color_pair(COL_HEADER_ID))
            if gCurrentMode == KEY_MODE_SEARCH_INPUT:
                handleKeysForSearchMode(ch, sheetWin.getTable(), theScreen)
                pass
            elif gCurrentMode == KEY_MODE_NORMAL_INPUT:
                handleKeysForNormalMode(ch, sheetWin.getTable(), theScreen)
                pass

        except KeyboardInterrupt:
            break  # Exit the program on a Ctrl-C as well. Regular terminal printing is automatically restored by "curses.wrapper"
        except:
            raise  # Something unexpected has happened. Better report it!

        pass

    return  # end of mainScreenHandlingLoop


def setUpCurses():  # initialize the curses environment
    if not curses.has_colors():
        print(
            "UH OH, this terminal does not support color! We might crash. Quitting now anyway until I figure out what to do. Sorry. This might not actually be a problem, but I will need to test it to see what happens in a non-color terminal!"
        )
        sys.exit(1)
        pass
    curses.start_color()
    curses.init_pair(RAGGED_END_ID, RAGGED_END_TEXT_COLOR, RAGGED_END_BG_COLOR)
    curses.init_pair(SELECTED_CELL_ID, SELECTED_CELL_TEXT_COLOR, SELECTED_CELL_BG_COLOR)
    curses.init_pair(COL_HEADER_ID, COL_HEADER_TEXT_COLOR, COL_HEADER_BG_COLOR)
    curses.init_pair(ROW_HEADER_ID, ROW_HEADER_TEXT_COLOR, ROW_HEADER_BG_COLOR)
    curses.init_pair(BOX_COLOR_ID, BOX_COLOR_TEXT_COLOR, BOX_COLOR_BG_COLOR)
    curses.init_pair(BLANK_COLOR_ID, BLANK_COLOR_TEXT_COLOR, BLANK_COLOR_BG_COLOR)
    curses.init_pair(SEARCH_MATCH_COLOR_ID, SEARCH_MATCH_COLOR_TEXT_COLOR, SEARCH_MATCH_COLOR_BG_COLOR)
    curses.init_pair(WARNING_COLOR_ID, WARNING_COLOR_TEXT_COLOR, WARNING_COLOR_BG_COLOR)
    curses.init_pair(
        NUMERIC_NEGATIVE_COLOR_ID, NUMERIC_NEGATIVE_COLOR_TEXT_COLOR, NUMERIC_NEGATIVE_COLOR_BG_COLOR
    )
    curses.init_pair(
        NUMERIC_POSITIVE_COLOR_ID, NUMERIC_POSITIVE_COLOR_TEXT_COLOR, NUMERIC_POSITIVE_COLOR_BG_COLOR
    )
    curses.init_pair(
        ACTIVE_FILENAME_COLOR_ID, ACTIVE_FILENAME_COLOR_TEXT_COLOR, ACTIVE_FILENAME_COLOR_BG_COLOR
    )
    curses.init_pair(HELP_AREA_ID, HELP_AREA_TEXT_COLOR, HELP_AREA_BG_COLOR)

    propid = LAST_USED_PROPERTY_ID + 1
    for propname, props in COLOR_PROPERTIES.items():
        curses.init_pair(propid, props["fg"], props["bg"])
        COLOR_PROPERTIES[propname]["id"] = propid  # save the identifying number
        propid = propid + 1
        pass

    CURSES_INVISIBLE_CURSOR = 0
    CURSES_VISIBLE_CURSOR = 1
    CURSES_HIGHLIGHTED_CURSOR = 2
    try:
        curses.curs_set(CURSES_INVISIBLE_CURSOR)  # Don't show a blinking cursor
    except (curses.error):
        print('Unable to set cursor state to "invisible"')
        pass

    curses.meta(1)  # Allow 8-bit chars
    return


def truncateLongCell(argString, argMaxlen, argTruncString):
    if len(argString) > argMaxlen:
        truncateToThisLen = argMaxlen - len(argTruncString)
        return argString[:truncateToThisLen] + argTruncString
    else:
        return argString
    pass


def padStrToLength(argString, argMaxlen, padChar):
    # pad out the length with blanks so that the ENTIRE
    # length is taken up. However! curses does not like to
    # draw just plain things unfortunately, so we have to
    # trick it with a -. This is very hackish. Maybe I
    # should draw colored boxes instead?
    numBlankSpacesToAdd = argMaxlen - len(argString)
    if numBlankSpacesToAdd > 0:
        return argString + str(padChar * numBlankSpacesToAdd)
    else:
        return argString
    return


def calculateBorderChar(r, c, topRow, leftCol, bottomRow, rightCol):
    if r == 0:  # TOP R
        if c == 0:
            ch = curses.ACS_ULCORNER
        elif c == rightCol:
            ch = curses.ACS_URCORNER
        else:
            ch = curses.ACS_TTEE
        pass
    elif r == bottomRow:
        if c == 0:
            ch = curses.ACS_LLCORNER
        elif c == rightCol:
            ch = curses.ACS_LRCORNER
        else:
            ch = curses.ACS_BTEE
        pass
    else:
        if c == 0:
            ch = curses.ACS_LTEE
        elif c == rightCol:
            ch = curses.ACS_RTEE
        else:
            ch = curses.ACS_PLUS
        pass

    return ch


def handleKeysForSearchMode(argCh, currentTable, theScreen):
    # If we are in search mode, then when the user types, that text is added to the query.
    finishedSearch = False
    cancelSearch = False

    if argCh in KEYS_SEARCH_MODE_FINISHED:  # User wants to STOP entering search text
        # ----------------------------------
        if len(mainInfo.getCurrent().getRegexString()) > 0:
            finishedSearch = True  # If there *is* a search string
            pass
        else:
            cancelSearch = True  # Search string is blank--user cancelled the search
            pass
        # ----------------------------------
        pass
    elif (
        argCh in KEYS_SEARCH_MODE_CANCEL
    ):  # User wants to CANCEL searching, deleting any query that was there
        cancelSearch = True
        # ----------------------------------
        pass
    elif argCh in KEYS_SEARCH_MODE_BACKSPACE or argCh in KEYS_SEARCH_MODE_DELETE_FORWARD:
        # we don't support forward-delete yet. sorry.
        # Pretend it's regular delete for now...
        if mainInfo.getCurrent().regexIsActive() is False:
            # The user deleted PAST the beginning of the string
            cancelSearch = True
            pass
        else:
            mainInfo.getCurrent().trimRegex(numChars=1)
            currentTable.regTab.clear()
            pass
        # ----------------------------------
        pass
    else:
        # Append whatever the user typed to the search string...
        try:
            charToAdd = chr(argCh)
            mainInfo.getCurrent().appendToCurrentSearchTerm(charToAdd)
            currentTable.regTab.clear()
            pass
        except (ValueError):
            # we don't care if "charToAdd" is not in the range of add-able chars
            pass
        # ----------------------------------
        pass

    if cancelSearch:
        setCommandStr("Search Cancelled")
        mainInfo.getCurrent().clearCurrentSearchTerm()
        currentTable.regTab.clear()
        GLOBAL_exitSearchMode()
        pass
    elif finishedSearch:
        setCommandStr('Searching for "' + mainInfo.getCurrent().getRegexString() + '"')
        GLOBAL_exitSearchMode()
        pass
    else:
        setCommandStr("Search (press enter when done): " + mainInfo.getCurrent().getRegexString())
        pass

    return


def handleKeysForNormalMode(argCh, currentTable, theScreen):
    # Handle user keyboard input when we are *not* in search mode.
    # This will handle the majority of user interactions.
    theAction = None  # What did the user want to do?
    wantToMove = Point(0, 0)  # See if we need to move the cursor / scroll up / down / left /right
    wantToChangeFileIdx = None  # Do we want to go to the next/previous file?
    activeCellPos = mainInfo.getCurrent().cursorPos  # Where are we currently, in the table?
    # if (argCh is not None and argCh is not ''): raise "ArgCh: " + str(argCh)
    if argCh in KEYS_QUIT:
        global gWantToQuit
        gWantToQuit = True
    elif argCh in KEYS_TOGGLE_HIGHLIGHT_NUMBERS_MODE:
        mainInfo.getCurrent().toggleNumericHighlighting()
    elif argCh in KEYS_MOVE_TO_TOP:
        activeCellPos.y = 0
    elif argCh in KEYS_MOVE_TO_BOTTOM:
        activeCellPos.y = currentTable.getNumRows() - 1
    elif argCh in KEYS_MOVE_RIGHT:
        wantToMove.x = 1
    elif argCh in KEYS_MOVE_LEFT:
        wantToMove.x = -1
    elif argCh in KEYS_MOVE_UP:
        wantToMove.y = -1
    elif argCh in KEYS_MOVE_DOWN:
        wantToMove.y = 1
    elif argCh in KEYS_MOVE_RIGHT_FAST:
        wantToMove.x = fastMoveSpeed.x
    elif argCh in KEYS_MOVE_LEFT_FAST:
        wantToMove.x = -fastMoveSpeed.x
    elif argCh in KEYS_MOVE_UP_FAST:
        wantToMove.y = -fastMoveSpeed.y
    elif argCh in KEYS_MOVE_DOWN_FAST:
        wantToMove.y = fastMoveSpeed.y
    elif argCh in KEYS_PREVIOUS_FILE:
        wantToChangeFileIdx = -1
    elif argCh in KEYS_NEXT_FILE:
        wantToChangeFileIdx = 1
    elif argCh in KEYS_GOTO_NEXT_MATCH:
        theAction = kWANT_TO_MOVE_TO_SEARCH_RESULT
        theActionParam = -1
    elif argCh in KEYS_GOTO_PREVIOUS_MATCH:
        theAction = kWANT_TO_MOVE_TO_SEARCH_RESULT
        theActionParam = 1
    elif argCh in KEYS_TRANSPOSE:  # Display the file in transposed format
        global gIsTransposed
        gIsTransposed = ~gIsTransposed
        temp = activeCellPos.x
        activeCellPos.x = activeCellPos.y
        activeCellPos.y = temp
        if gIsTransposed:
            setCommandStr("Now displaying the file in TRANSPOSED format.")
            pass
        else:
            setCommandStr("Now displaying the file in NON-TRANSPOSED format again.")
            pass
        drawEverything(theScreen)
    elif argCh in KEYS_GOTO_LINE_END:
        activeCellPos.x = currentTable.getNumCols() - 1
    elif argCh in KEYS_GOTO_LINE_START:
        activeCellPos.x = 0
    # elif argCh in (ord('w'),): sheetWin.win.mvwin(20,0)
    # elif argCh in (ord('R'),): gCellBorders.width = (1 - gCellBorders.width)  # Add/remove column borders (not useful)
    # elif argCh in (ord('r'),): gCellBorders.height = (1 - gCellBorders.height) ## Add row borders (not useful)
    elif argCh in KEYS_WANT_TO_ENTER_SEARCH_MODE:
        GLOBAL_setUserInteractionMode(KEY_MODE_SEARCH_INPUT)
        mainInfo.getCurrent().clearCurrentSearchTerm()
        currentTable.initRegexTable()
        setCommandStr("Search (press enter when done): ")
        pass
    else:
        pass  # unrecognized key

    if theAction is None:
        pass
    elif theAction == kWANT_TO_MOVE_TO_SEARCH_RESULT:
        setCommandStr("Sorry! Not implemented yet.")
        pass

    if wantToMove.x != 0 or wantToMove.y != 0:
        activeCellPos.y = max(0, min(currentTable.getNumRows() - 1, (activeCellPos.y + wantToMove.y)))
        activeCellPos.x = max(0, min(currentTable.getNumCols() - 1, (activeCellPos.x + wantToMove.x)))
        clearCommandStr()  # whatever the previous command was, it no longer applies now that we have moved
        pass

    if wantToChangeFileIdx is not None:
        possibleNewIndex = mainInfo.currentFileIdx + wantToChangeFileIdx
        if mainInfo.getCurrent().filename is kSTDIN_HYPHEN:
            setWarning(
                ">>> Cannot change files, because we read input from an STDIN pipe. Multiple files cannot be read when reading from a pipe."
            )
            pass
        elif possibleNewIndex < 0:
            setWarning(
                ">>> Cannot go to the previous file, because we are already at the beginning of the file list."
            )
            pass
        elif possibleNewIndex >= mainInfo.size():
            setWarning(
                ">>> Cannot go to the next file, because we are already at the end of the file list."
            )
            pass
        else:
            sheetWin.win.clear()
            infoWin.win.clear()
            colHeaderWin.win.clear()
            colHeaderWin.win.refresh()
            rowHeaderWin.win.clear()
            mainInfo.currentFileIdx = possibleNewIndex  # update which file we are currently examining
            initializeWindowSettings(theScreen, mainInfo.getCurrent())
            setCommandStr('Changed to the file named "' + mainInfo.getCurrent().filename + '".')
            currentTable.loadNewFile(mainInfo.getCurrent().filename)
        pass
    return


def GLOBAL_exitSearchMode():
    if gCurrentMode != KEY_MODE_SEARCH_INPUT:
        raise Exception("Uh oh, tried to exit search mode... but we were not even IN search mode!!")
    else:
        GLOBAL_setUserInteractionMode(KEY_MODE_NORMAL_INPUT)
        pass
    return


def GLOBAL_setUserInteractionMode(argNewMode):
    global gCurrentMode
    gCurrentMode = argNewMode
    return


ALEX_PROGRAM_USAGE_TEXT = """sheet.py
A program for displaying tab-delimited files in spreadsheet format.

This is a file-viewing program only, not an editor.

Written for Python 2.5.1+ by Alex Williams, 2010-2015.

Apr 2015: Added transparent support for gzip and bzip2 file reading.

Dec 2018: Added support for automatic delimiter guessing.

Requires the "curses" terminal module, which is built-in
with most python distributions.


Usage:
   sheet.py INPUT_FILENAMES

or to read from STDIN:

   cat INPUT_FILE | sheet.py

Example usage:

sheet.pl myfile.tab
	This reads in the tab-delimited file "myfile.tab"

sheet.pl file1.tab anotherFile.tab
	This reads in two files. "file1.tab" is the first
	one you will see. You would have to press '>' in
	order to go to "anotherFile.tab".

Options and arguments:
-h or --help:
	Print this usage/help message

-w or --warn:
	Enable debugging warnings

-d or --delim="DELIMITER"
        Manually specify a delimiter (e.g. -d '|' for the pipe, or -d ';')
        Not necessary if your filename already ends in '.tsv/.tab/.txt' or 'csv'.

Controls:

Navigation keys:
[  i  ]      [ Arrow  ]
[ jkl ]  or  [  Keys  ]:   (Hold shift to move faster)

	* Move the cursor around. Note that "IJKL" makes an
	  inverted-T shape, just like the arrow keys on most
	  keyboards.

	* Hold "shift" along with IJKL or the arrow keys in
          order to move faster.

	* Space Bar (down) and "b" (up) work to move
          quickly, just like in LESS.

	* Also supported is the standard terminal/emacs
          Ctrl-P/N/B/F for moving the cursor around.

Fast navigation keys:
[   g  ]
[ a   e]
[   G  ]:
	* Move to the top (g) or bottom (G) of a file (same
          keybindings as in LESS).

	* Move to the leftmost column (a) or rightmost
          column (e), with keybindings similar to emacs's.

[ < ] and [ > ]:

	* Switch between files, when more than one file was
          specified. Not valid when the input is STDIN.

[ / ]:
	* Set the search term. Text that matches will be
          highlighted. Hit "enter" to finish the search.
          Press [ / ] again to cancel.

[ q ]:
	* Quit the program.


BUGS / TO DO:
  * When you load more than one file, the header appears again about 30 rows down. Weird.
  * Should add a "--vi" switch to give vi-style hjkl / jkl; navigation.
"""


# class Dispatcher:
#     def __init__(self):
#         self._keys = []

#     def add_listener(self, key_list, method):
#         for key in key_list:
#             self._keys[key] = method

#     def move_up(self, key):
#         pass

#     def dispatch(self, key):
#         try:
#             method = self._keys[key]
#             return self.__getattr__(method)(self, key)
#         except:
#             return


# d = Dispatcher()
# d.add_listener([33, 34 35], 'move_up')
# d.add_listener([33, 34 35], 'move_down')


# def main_loop():
#     while True:
#         key = get_key()
#         d.dispatch(key)

# Must come at the VERY END!
if __name__ == "__main__":
    processCommandLineArgs(sys.argv[1:])
    curses.wrapper(mainScreenHandlingLoop)
    pass
