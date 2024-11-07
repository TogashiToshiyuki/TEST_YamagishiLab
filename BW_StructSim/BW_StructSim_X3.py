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


def main():
    # program start
    # View Program Overview
    print(Stereotyped.ProgramAbst)

    # argument parser
    args, before = arg_parser()

    # create an instance of the BrickWork class
    bw = BrickWork(args)

    # Check if the required files exist
    bw.File_Set_Check()


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
                   "*     - CalcSetting_BW.txt                                              *\n"
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
    InitialCondition_3mol_Temp = ("-3.0\n"
                                  "-2.5\n"
                                  "-2.0\n"
                                  "-1.5\n"
                                  "-1.0\n"
                                  "-0.5\n"
                                  "0\n"
                                  "0.5\n"
                                  "1.0\n"
                                  "1.5\n"
                                  "2.0\n"
                                  "2.5\n"
                                  "3.0")


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


def arg_parser():
    """
    Argument parser
    :return: args, before
    :rtype: argparse.Namespace, float
    """
    before = time.time()
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('MaterName',
                        help="Molecular Name")
    parser.add_argument('--debug', '-d', '-D',
                        help="Start the programme in Debug mode.",
                        action="store_true")
    parser.add_argument('--tcal', '-t',
                        help="Argument for not calculating tcal.",
                        action="store_false")
    parser.add_argument('--chk', '--Check', '-c',
                        help="Check the structure.",
                        action="store_true")
    parser.add_argument('--xyz', '--XYZ',
                        action=CheckRequired,
                        nargs='?',
                        const=True, default=False,
                        help='Create .xyz files')
    parser.add_argument('--manual', '-m',
                        action=CheckRequired,
                        nargs='?',
                        const=True, default=False,
                        help='Create files of any condition')

    args = parser.parse_args()

    if args.debug:
        print(f"{Color.RED}*************** Caution!!! Debug Started!!! ***************{Color.RESET}")
        print(f"In debug mode, the following files are not deleted\n"
              f"\t - sh files\n"
              f"\t - chk files\n"
              f"\t - log files\n"
              f"\t - gjf files\n")
    else:
        pass

    return args, before


class BrickWork:
    def __init__(self, args):
        self.MaterName = args.MaterName
        self.Flag_xyz = args.xyz
        while True:
            mol_pos_number = input(f"\n"
                                   f"Calculation from 3 molecules have been selected.\n"  # 3molでの計算が選ばれました
                                   f"{Color.GREEN}Please select a structure from 1~3..{Color.RESET}\n"
                                   f"\t 1: H-H\n"
                                   f"\t 2: T-T\n"
                                   f"\t 3: H-T\n"
                                   f"\t >>> ")
            try:
                mol_pos_number = int(mol_pos_number)
            except (IndexError, ValueError):
                print(f"{Color.RED}\t Error: Enter a number.{Color.RESET}")
                continue
            if mol_pos_number in {1, 2, 3}:
                break
            else:
                print(f"{Color.RED}\t Error: Incorrect input.{Color.RESET}")
                continue
        self.mol_pos = f"p{mol_pos_number}"
        self.dirpath = f"./{self.MaterName}_3mol{self.mol_pos}"
        Operator = input(f"\n{Color.GREEN}Retrieve the operator name.{Color.RESET}\n"
                         "\t>>> The name entered here will be used to identify the operator.\n"
                         f"\tPlease enter the operator.\n"
                         f"\t{Color.GREEN}>>> {Color.RESET}")
        if Operator == "":
            Operator = "ONE"
        else:
            pass
        print("")
        self.Operator = Operator
        self.messages, self.HelpList = [], []

    def help_check_exit(self):
        """
        Check if the required files exist.
        :return:
        """
        self.message_show()
        if True in self.HelpList:
            print(f"{Color.RED}{Stereotyped.AbnormalEnd}{Color.RESET}")
            exit()
        else:
            pass
        return None

    def message_show(self):
        """
        Show messages.
        :return:
        """
        for message in self.messages:
            print(message)
        self.messages.clear()
        return None

    def File_Set_Check(self):
        """
        Check if the required files exist.
        :return:
        """
        print(f"{Color.GREEN}Checking the required files...{Color.RESET}")

        if os.path.exists(f"./{self.MaterName}.xyz"):
            self.messages.append(f"\t>>> {self.MaterName}.xyz: {Color.GREEN}Exist.{Color.RESET}")
            self.HelpList.append(False)
        else:
            self.messages.append(f"\t>>> {Color.RED}{self.MaterName}.xyz: Does not Exist.{Color.RESET}")
            self.HelpList.append(True)

        if os.path.exists(f"./CalcSetting_BW.txt"):
            self.messages.append(f"\t>>> CalcSetting_BW.txt: {Color.GREEN}Exist.{Color.RESET}")
            self.HelpList.append(False)
        else:
            self.messages.append(f"\t>>> {Color.RED}CalcSetting_BW.txt: Does not Exist.{Color.RESET}")
            self.HelpList.append(True)
            with open(f"./CalcSetting_BW.txt", "w") as f:
                f.write(Stereotyped.CalcSet_BW_template)
            self.messages.append(f"\t>>> CalcSetting_BW.txt: {Color.GREEN}Created.{Color.RESET}\n"
                                 f"\t>>> Please set the calculation conditions in the file.")
        self.help_check_exit()

        if os.path.exists(f"./ConditionList_3mol{self.mol_pos}.txt"):
            self.messages.append(f"\t>>> ConditionList_3mol{self.mol_pos}.txt: {Color.GREEN}Exist.{Color.RESET}")
            self.HelpList.append(False)
        elif os.path.exists(f"./InitialCondition_3mol{self.mol_pos}"):
            self.messages.append(f"\t>>> InitialCondition_3mol{self.mol_pos}.txt: {Color.GREEN}Exist.{Color.RESET}")
            self.HelpList.append(False)
        else:
            self.messages.append(
                f"\t>>> {Color.RED}ConditionList_3mol{self.mol_pos}.txt: Does not Exist.{Color.RESET}\n"
                f"\t>>> {Color.RED}InitialCondition_3mol{self.mol_pos}.txt: Does not Exist.{Color.RESET}\n")
            self.HelpList.append(True)
            with open(f"./InitialCondition_3mol{self.mol_pos}.txt", "w") as f:
                f.write(Stereotyped.InitialCondition_3mol_Temp)
            self.messages.append(f"\t>>> InitialCondition_3mol{self.mol_pos}.txt: {Color.GREEN}Created.{Color.RESET}\n"
                                 f"\t>>> Please set the calculation conditions in the file.")
        self.help_check_exit()
        return None


if __name__ == "__main__":
    main()
