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


def main():
    # program start
    # View Program Overview
    print(Stereotyped.ProgramAbst)

    # Define variables
    messages, HelpList, MaterName = [], [], ""

    # Argument parsing
    args, MaterName, Debug, Nmol, calculation_tcal_Flag = arg_parser(messages, HelpList)

    # Acquisition of structure
    mol_pos = structure_acquisition(Nmol)

    # Get the operator
    Operator = getOperator()


# Check if help has occurred, if it has occurred, exit, if it has not occurred, clear messages
def help_check_exit(messages, HelpList):
    """
    Check if help has occurred, if it has occurred, exit, if it has not occurred, clear messages
    :param messages: list of messages to be displayed
    :type messages: list
    :param HelpList: list of help flags
    :type HelpList: list
    """
    message_show(messages)
    if True in HelpList:
        print(f"{Color.RED}{Stereotyped.AbnormalEnd}{Color.RESET}")
        exit()
    else:
        pass
    return None


# Display message
def message_show(messages):
    """
    Display message
    :param messages: list of messages to be displayed
    :type messages: list
    """
    for message in messages:
        print(message)
    messages.clear()
    return None


# Argument parsing
def arg_parser(messages, HelpList):
    """
    Argument parsing
    :param messages: list of messages to be displayed
    :param HelpList: list of help flags
    :return:
    """
    Help_Text = "現在作成中です。"

    parser = argparse.ArgumentParser(description=Help_Text, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('MaterNameXYZ',
                        help="MaterNameXYZ.xyz")
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

    # Create a mutually exclusive group that requires one argument
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--two_mol', '-2', '--2mol',
                       help="Selection of 2mol calculations",
                       action="store_true")
    group.add_argument('--three_mol', '-3', '--3mol',
                       help="Selection of 3mol calculations",
                       action="store_true")

    args = parser.parse_args()
    if args.debug:
        Debug = True
        messages.append(f"{Color.RED}*************** Caution!!! Debug Started!!! ***************{Color.RESET}\n")
        messages.append(f"In debug mode, the following files are not deleted\n"
                        f"\t - sh files\n"
                        f"\t - chk files\n"
                        f"\t - log files\n"
                        f"\t - gjf files\n")
    else:
        Debug = False

    if os.path.exists(args.MaterNameXYZ) and args.MaterNameXYZ.endswith(".xyz"):
        MaterName = args.MaterNameXYZ[:-4]
        HelpList.append(False)
    else:
        MaterName = ""
        messages.append(f"{Color.RED} file selected. [molecule name].xyz is required.{Color.RESET}")
        HelpList.append(True)

    if args.chk:
        mkCheckFile(args, Debug, messages, HelpList, MaterName)

    if args.two_mol and args.three_mol:
        messages.append(f"\t>>>{Color.RED}2mol and 3mol cannot be selected at the same time.{Color.RESET}\n")
        HelpList.append(True)

    if not args.two_mol and not args.three_mol:
        messages.append(f"\t>>>{Color.RED}2mol or 3mol must be selected.{Color.RESET}\n")
        HelpList.append(True)

    calculation_tcal_Flag = args.tcal
    if args.tcal:
        pass
    else:
        messages.append(f"{Color.RED}Runs the program without transfer integral calculations.{Color.RESET}")

    Nmol = ""
    if args.two_mol:
        Nmol = "2mol"
    if args.three_mol:
        Nmol = "3mol"

    help_check_exit(messages, HelpList)

    return args, MaterName, Debug, Nmol, calculation_tcal_Flag


# 実装は待機待ち中（優先度低）
# Create a file for structural verification
def mkCheckFile(args, Debug, messages, HelpList, MaterName):
    print(f"{Color.RED}"
          f"********** Starts generating files for structural verification **********\n"
          f"{Color.RESET}")


# Acquisition of structure
def structure_acquisition(Nmol):
    """
    function to acquire the structure
    :param Nmol:
    :type Nmol: str
    :rtype: str
    :return:
    """
    if "3mol" in Nmol:
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
        mol_pos = f"p{mol_pos_number}"
    else:
        mol_pos = ""

    return mol_pos


# Check if the required files exist
def File_Set_Check(messages, HelpList, MaterName, Nmol, mol_pos):
    print(f"\n{Color.GREEN}File set check...{Color.RESET}")
    if os.path.exists(f"./{MaterName}.xyz"):
        messages.append(f"\t>>> {MaterName}.xyz: Found")
        HelpList.append(False)
    else:
        messages.append(f"\t{Color.RED}>>> {MaterName}.xyz: NOT Found{Color.RESET}")
        HelpList.append(True)

    if os.path.exists(f"./CalcSetting_BW.txt"):
        messages.append(f"\t>>> CalcSetting_BW.txt: Found")
        HelpList.append(False)
    else:
        messages.append(f"\t{Color.RED}>>> CalcSetting_BW.txt: NOT Found{Color.RESET}\n"
                        f"\tMake the correct CalcSetting_BW.txt and then restart the program.")
        HelpList.append(True)
        with open("CalcSetting_BW.txt", "w") as f:
            f.write(Stereotyped.CalcSet_BW_template)
        messages.append(f"\t>>> CalcSetting_BW.txt: Created")
    if True in HelpList:
        messages.append(f"{Color.RED}Required file set DOES NOT exist in the correct directory.{Color.RESET}")
    else:
        pass
    help_check_exit(messages, HelpList)

    if os.path.exists(f"./ConditionList_{Nmol}{mol_pos}.txt"):
        messages.append(f"\t>>> ./ConditionList_{Nmol}{mol_pos}.txt: Found")
        HelpList.append(False)
    elif os.path.exists(f"./InitialCondition_{Nmol}{mol_pos}d.txt"):
        messages.append(f"\t>>> ./InitialCondition_{Nmol}{mol_pos}.txt: Found")
        HelpList.append(False)
    else:
        messages.append(f"\t>>> {Color.RED}./ConditionList_{Nmol}{mol_pos}.txt and "
                        f"./InitialCondition_{Nmol}{mol_pos}.txt: Not Found{Color.RESET}\n"
                        f"\t    Make the correct InitialCondition_{Nmol}{mol_pos}.txt "
                        f"and then restart the program.")
        with open(f"InitialCondition_{Nmol}{mol_pos}.txt", "w") as file:
            if "2mol" in Nmol:
                f.write(Stereotyped.InitialCondition_2mol_Temp)
            else:
                f.write(Stereotyped.InitialCondition_3mol_Temp)
        HelpList.append(True)

    # MaterName_2mol_min
    if "3mol" in Nmol:
        if os.path.exists(f"./{MaterName}_2mol_min.txt"):
            messages.append(f"\t>>> {MaterName}_2mol_min.txt: Found")
            HelpList.append(False)
        else:
            messages.append(f"\t>>> {Color.RED}{MaterName}_2mol_min.txt: NOT Found{Color.RESET}")
            HelpList.append(True)
    elif "2mol" in Nmol:
        pass

    if True in HelpList:
        messages.append(f"{Color.RED}Required file set DOES NOT exist in the correct directory.{Color.RESET}")
    else:
        messages.append(f"Required file set {Color.GREEN}EXIST{Color.RESET} in the correct directory.")

    help_check_exit(messages, HelpList)

    return None


# Get the operator
def getOperator():
    """
    Get the operator
    :rtype: str
    :return: operator
    """
    print(f"\n{Color.GREEN}Retrieve the operator name.{Color.RESET}\n"
          "\t>>> The name entered here will be used to identify the operator.")
    Operator = input(f"\tPlease enter the operator.\n"
                     f"\t{Color.GREEN}>>> {Color.RESET}")
    if Operator == "":
        Operator = "ONE"
    else:
        pass
    return Operator


if __name__ == "__main__":
    main()