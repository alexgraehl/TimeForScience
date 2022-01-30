#!/usr/bin/env python3

"""
A python3 micro-program to show how CURSES (terminal "drawing" output) works. Just shows a grid of colored 'a's.
"""

import sys

import curses  # <-- docs at: http://docs.python.org/library/curses.html
import curses.ascii  # <-- docs at http://docs.python.org/library/curses.ascii.html
import curses.textpad

# import textwrap
# import getopt

# import pdb
# pdb.set_trace() ## Python Debugger! See: http://aymanh.com/python-debugging-techniques


# import time # we just want "sleep"
# import os.path

grid_win = None


def agwInitCurses():  # initialize the curses environment
    if not curses.has_colors():
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
        curses.curs_set(CURSES_INVISIBLE_CURSOR)  # Don't show a blinking cursor
    except (curses.error):
        print('Unable to set cursor state to "invisible"')
        pass
    curses.meta(1)  # Allow 8-bit chars
    return


def mainScreenHandlingLoop(the_screen):
    agwInitCurses()
    sizeyx = the_screen.getmaxyx()
    grid_win = curses.newwin(sizeyx[0], sizeyx[1], 0, 0)
    ch = None

    x_size = 80
    y_size = 20
    while True:  # <-- infinite loop!
        try:
            the_screen.refresh()  # Refresh it first..

            for xxx in range(x_size):
                for yyy in range(y_size):
                    attr = (
                        curses.color_pair(((xxx + yyy) % 2) + 1) + curses.A_NORMAL
                    )  # Set text color (note the '% 2' modular division to alternate between red and green)
                    grid_win.attron(attr)
                    grid_win.addch(yyy, xxx, "a", attr)  # Just add the letter 'a' at this position
                    grid_win.attroff(attr)
                    pass
                pass
            grid_win.addstr("\n")
            grid_win.addstr("🔥🔥🔥 (Press Ctrl-C to exit) 🔥🔥🔥")
            grid_win.refresh()
            ch = the_screen.getch()  # Now get input from the user... (unused in this example code)

        except KeyboardInterrupt:
            break  # Exit the program on a Ctrl-C as well. Regular terminal printing is automatically restored.
        except:
            raise  # Something unexpected has happened. Better report it!

        pass

    return  # end of mainScreenHandlingLoop


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
