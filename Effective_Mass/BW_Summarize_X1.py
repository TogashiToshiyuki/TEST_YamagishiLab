#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import datetime
import functools
import glob
import math
import os
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from matplotlib.animation import FuncAnimation
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Pt
from scipy.interpolate import griddata
from scipy.optimize import differential_evolution
from tabulate import tabulate

print = functools.partial(print, flush=True)


def main():
    args, before = get_args()


def get_args():
    before = time.time()
    help_text = StandardPhrases.ProgramAbst

    parser = argparse.ArgumentParser(description=help_text, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('MaterName',
                        help='Material name')
    parser.add_argument('--debug', '-d', '-D',
                        help="Start in debug mode",
                        action='store_true')
    parser.add_argument('--gif', '-g', '-G',
                        help="Create GJF File from All File",
                        action='store_true')
    args = parser.parse_args()

    if args.debug:
        print(f"{Color.RED}*************** Caution!!! Debug Started!!! ***************{Color.RESET}\n")

    return args, before


class StandardPhrases:
    def __init__(self):
        self._StandardPhrases = "StandardPhrases"

    ProgramAbst = ("\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n"
                   "*   This program collects and processes calculation results of                *\n"
                   "*       optimized aggregate structures of organic semiconductors              *\n"
                   "*       with HerringBone (HB) structures.                                     *\n"
                   "*   It performs the following tasks:                                          *\n"
                   "*       - Calculates effective masses                                         *\n"
                   "*       - Plots band diagrams                                                 *\n"
                   "*       - Generates summary slides of the results                             *\n"
                   "*                                                                             *\n"
                   "*   To run this program, please ensure the following files are available:     *\n"
                   "*       - “MoleculeName”_3molp1_t”N”d_results                                 *\n"
                   "*       - “MoleculeName”_3molp2_t”N”d_results                                 *\n"
                   "*       - “MoleculeName”_3molp3_t”N”d_results                                 *\n"
                   "*                                                                             *\n"
                   "*   Replace “MoleculeName” with the actual name of your molecule,             *\n"
                   "*       and “N” with the appropriate value used in your calculations.         *\n"
                   "*                                                                             *\n"
                   "*   When you start the program, you will be prompted to enter your name,      *\n"
                   "*       - which will be used as the slide creator’s name.                     *\n"
                   "*                                                                             *\n"
                   "*                                 created by Toshiyuki Togashi 2024/11/27     *\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n")
    AbnormalEnd = "\n************* Programme DID NOT terminate successfully. *************\n"
    HelpText = ("\n"
                "The required file may not exist.\n"
                "Please check.")


class Color:
    """
    ANSI escape code for color
    """

    def __init__(self):
        self._Color = "Color"

    BLACK = '\033[30m'  # black
    RED = '\033[31m'  # red
    GREEN = '\033[32m'  # green
    YELLOW = '\033[33m'  # yellow
    BLUE = '\033[34m'  # blue
    MAGENTA = '\033[35m'  # magenta
    CYAN = '\033[36m'  # cyan
    LIME = '\033[92m'  # lime
    WHITE = '\033[37m'  # white
    COLOR_DEFAULT = '\033[39m'  # default
    BOLD = '\033[1m'  # bold
    UNDERLINE = '\033[4m'  # underline
    INVISIBLE = '\033[08m'  # invisible
    REVERSE = '\033[07m'  # reverse
    BG_BLACK = '\033[40m'  # (background)black
    BG_RED = '\033[41m'  # (background)red
    BG_GREEN = '\033[42m'  # (background)green
    BG_YELLOW = '\033[43m'  # (background)yellow
    BG_BLUE = '\033[44m'  # (background)blue
    BG_MAGENTA = '\033[45m'  # (background)magenta
    BG_CYAN = '\033[46m'  # (background)cyan
    BG_WHITE = '\033[47m'  # (background)white
    BG_DEFAULT = '\033[49m'  # (background)default
    RESET = '\033[0m'  # reset
