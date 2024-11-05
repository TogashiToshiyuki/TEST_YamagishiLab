#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import functools
import os
import numpy as np
import math
import datetime
import subprocess
import sys
import time
import glob

print = functools.partial(print, flush=True)

class Stereotyped:
    def __init__(self):
        self._stereotype = "Stereotyped"

    ProgramAbst = ("\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n"
                   "*   To Search the energetically stable aggregation                      *\n"
                   "*            for Herringbone (HB) aggregation                           *\n"
                   "*                                                                       *\n"
                   "*   The programme requires the following arguments.                     *\n"
                   "*     - xyz file of the molecule from which the calculations are made.  *\n"
                   "*                                                                       *\n"
                   "*   The following files are also required.                              *\n"
                   "*     - CalcSetting_HB.txt                                              *\n"
                   "*                                                                       *\n"
                   "*    Option                                                             *\n"
                   "*      - Tilt angle can be added to the molecule                        *\n"
                   "*                          created by Toshiyuki Togashi 2024/10/31      *\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n")
    AbnormalEnd = "\n************* Programme DID NOT terminate successfully. *************\n"
    HelpText = ("\n"
                "The required file may not exist.\n"
                "Please check.")
    # memo: Headers and text for NITT WS, need to be adjusted if used elsewhere
    Sh_txt = ("#!/bin/sh\n"
              "\n"
              "#$ -S /bin/sh\n"
              "#$ -cwd\n"
              "#$ -V\n"
              "#$ -pe gau 12\n"
              "#$ -q all.q\n"
              "\n"
              "module load gaussian/g16\n"
              "export GAUSS_SCRDIR=/scr/$JOB_ID\n"
              "mkdir /scr/$JOB_ID\n"
              "\n"
              "g16 xxx.gjf\n"
              "rm -rf /scr/$JOB_ID\n"
              "\n")
    tcal_sh_txt = ("#!/bin/sh\n"
                   "\n"
                   "#$ -S /bin/sh\n"
                   "#$ -cwd\n"
                   "#$ -V\n"
                   "#$ -pe gau 12\n"
                   "#$ -q all.q\n"
                   "\n"
                   "tcal *.xyz\n"
                   "\n"
                   "g16 xxx.gjf\n"
                   "rm -rf /scr/$JOB_ID\n"
                   "\n")
    Header_2mol = ("\n%mem=24GB"
                   "\n%nprocshared=12"
                   "\n%chk=test_2mol.chk"
                   "\n#p pbepbe/6-31g(d) empiricaldispersion=gd3 counterpoise=2"
                   "\n"
                   "\nTitle Card Required"
                   "\n"
                   "\n0 1"
                   "\n")
    Header_3mol = ("\n"
                   "%mem=24GB\n"
                   "%nprocshared=12\n"
                   "%chk=test_3mol.chk\n"
                   "#p pbepbe/6-31g(d) empiricaldispersion=gd3 counterpoise=3\n"
                   "\n"
                   "Title Card Required\n"
                   "\n"
                   "0 1\n")
    CalcSet_BW_template = ("Edge Axis: [x, y, or z]\n"
                           "Faceon Axis: [x, y, or z]\n"
                           "Mol3 Other_Transition: A.AA\n"
                           "\n"
                           "[Comment]\n")
    InitialCondition_2mol_Temp = "8.0\n"
    InitialCondition_3mol_Temp = ("-3.0 3.9\n"
                                      "-2.5 3.9\n"
                                      "-2.0 3.9\n"
                                      "-1.5 3.9\n"
                                      "-1.0 3.9\n"
                                      "-0.5 3.9\n"
                                      "0 3.9\n"
                                      "0.5 3.9\n"
                                      "1.0 3.9\n"
                                      "1.5 3.9\n"
                                      "2.0 3.9\n"
                                      "2.5 3.9\n"
                                      "3.0 3.9")


class Color:
    """
    ANSI escape code for color
    """
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


class Constant:
    Cn = 1
    d_another = 0
    CycleCondition_n_02 = 1
    CycleCondition_n_01 = 1
    CycleCondition_n_005 = 1


class CheckRequired(argparse.Action):
    """
    Check if the option is required.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        if not getattr(namespace, 'chk', False):
            parser.error(f"{option_string} requires --chk")
        setattr(namespace, self.dest, values)