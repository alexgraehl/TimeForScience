#!/usr/bin/python

'''
Gene Test

A python program.
'''

import sys

import curses # <-- docs at: http://docs.python.org/library/curses.html
import curses.ascii   # <-- docs at http://docs.python.org/library/curses.ascii.html
import curses.textpad
import curses.wrapper

import pdb; #pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques

import textwrap
#import time # we just want "sleep"
#import os.path
import getopt

try: 
    import numpy #from numpy import * # The matrices use this
except ImportError:
    print("[TERMINATING EARLY] because the Python module 'numpy' was not installed.")
    sys.exit(1)
    pass

class GameGrid:
    def __init__(self, xsize, ysize):
        self.xsize = xsize
        self.ysize = ysize

        self.grid = numpy.zeros(shape=(xsize,ysize), dtype=int)
        pass

gg = GameGrid(xsize=20,ysize=15)

gridWin = None

def agwInitCurses(): # initialize the curses environment                                                                           
    if (not curses.has_colors()):
        print("UH OH, this terminal does not support color! Quitting.")
        sys.exit(1)
        pass
    curses.start_color()
    ## Initialize some colors below. The first number is an ID and must be UNIQUE.
    curses.init_pair(1, curses.COLOR_GREEN, curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_RED, curses.COLOR_BLACK)
    ## Done initializing some colors
    CURSES_INVISIBLE_CURSOR = 0
    CURSES_VISIBLE_CURSOR = 1
    CURSES_HIGHLIGHTED_CURSOR = 2
    try:
        curses.curs_set(CURSES_INVISIBLE_CURSOR) # Don't show a blinking cursor                                                  
    except(curses.error):
        print("Unable to set cursor state to \"invisible\"")
        pass
    curses.meta(1)  # Allow 8-bit chars
    return




def mainScreenHandlingLoop(theScreen):
    agwInitCurses()

    sizeyx = theScreen.getmaxyx()
    
    gridWin = curses.newwin(sizeyx[0], sizeyx[1], 0, 0)
    
    ch = None
    while (1):
        try:
            theScreen.refresh() # Refresh it first..

            for xxx in range(0, gg.grid.shape[0]):
                for yyy in range(0, gg.grid.shape[1]):
                    attr = (curses.color_pair(((xxx+yyy)%2) + 1) + curses.A_NORMAL)
                    gridWin.attron(attr)
                    gridWin.addch(yyy, xxx, 'a', attr)
                    gridWin.attroff(attr)
                    pass
                pass
            gridWin.refresh()
            ch = theScreen.getch()    # Now get input from the user...
            
        except KeyboardInterrupt:
            break # Exit the program on a Ctrl-C as well. Regular terminal printing is automatically restored by "curses.wrapper"
        except:
            raise # Something unexpected has happened. Better report it!

        pass
    
    return # end of mainScreenHandlingLoop


# class KeyDispatcher:
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


# d = KeyDispatcher()
# d.add_listener([33, 34 35], 'move_up')
# d.add_listener([33, 34 35], 'move_down')


# def main_loop():
#     while True:
#         key = get_key()
#         d.dispatch(key)



# Must come at the VERY END!
if __name__ == "__main__":
    curses.wrapper(mainScreenHandlingLoop)
    pass





