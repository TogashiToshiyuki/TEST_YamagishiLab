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

    args, MaterName, Debug, Nmol, calculation_tcal_Flag, before = arg_parser(messages, HelpList)

    bw = BrickWall(args.MaterNameXYZ, args)

    # Define variables
    messages, HelpList, MaterName = [], [], ""

    bw.File_Set_Check()

    RefLines, dev, which = bw.mkFirstCondition()


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


def arg_parser(messages, HelpList):
    """
    Argument parser
    """
    before = time.time()
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('MaterNameXYZ',
                        help="Molecular_Name.xyz")
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
        mkCheckFile(args, Debug, messages, HelpList, MaterName, before)

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

    return args, MaterName, Debug, Nmol, calculation_tcal_Flag, before


def mkCheckFile(args, Debug, messages, HelpList, MaterName, before):
    print(f"{Color.RED}"
          f"********** Starts generating files for structural verification **********\n"
          f"{Color.RESET}")


class BrickWall:
    def __init__(self, file, args):
        self._base_path = os.path.splitext(file)[0]
        self.MaterName = file[:-4]
        if args.two_mol:
            self.Nmol = "2mol"
        if args.three_mol:
            self.Nmol = "3mol"
        if self.Nmol == "2mol":
            self.mol_pos = ""
        elif self.Nmol == "3mol":
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
        else:
            self.mol_pos = ""
        self.dirpath = f"./{self.MaterName}_{self.Nmol}{self.mol_pos}"
        self.tcalpath = f"./{self.MaterName}_{self.Nmol}{self.mol_pos}_tcal"
        self.messages = []
        self.HelpList = []
        print(f"\n{Color.GREEN}Retrieve the operator name.{Color.RESET}\n"
              "\t>>> The name entered here will be used to identify the operator.")
        Operator = input(f"\tPlease enter the operator.\n"
                         f"\t{Color.GREEN}>>> {Color.RESET}")
        if Operator == "":
            Operator = "ONE"
        else:
            pass
        print("")
        self.Operator = Operator
        self.Flag_XYZ = args.xyz

    def help_check_exit(self):
        self.message_show()
        if True in self.HelpList:
            print(f"{Color.RED}{Stereotyped.AbnormalEnd}{Color.RESET}")
            exit()
        else:
            pass
        return None

    def message_show(self):
        for message in self.messages:
            print(message)
        self.messages.clear()
        return None

    def File_Set_Check(self):
        print(f"{Color.GREEN}File set Check...{Color.RESET}")

        if os.path.exists(f"./{self.MaterName}.xyz"):
            self.messages.append(f"\t>>> {self.MaterName}.xyz: Exist.")
            self.HelpList.append(False)
        else:
            self.messages.append(f"\t{Color.RED}>>> {self.MaterName}.xyz: Not exist.{Color.RESET}")
            self.HelpList.append(True)

        if os.path.exists(f"./CalcSetting_BW.txt"):
            self.messages.append(f"\t>>> CalcSetting_BW.txt: Exist.")
            self.HelpList.append(False)
        else:
            self.messages.append(f"\t{Color.RED}>>> CalcSetting_BW.txt: Not exist.{Color.RESET}")
            self.HelpList.append(True)
            with open(f"./CalcSetting_BW.txt", "w") as f:
                f.write(Stereotyped.CalcSet_BW_template)
            self.messages.append(f"\t>>> CalcSetting_BW.txt: {Color.GREEN}Created.{Color.RESET}\n"
                                 f"\t>>> Please set the calculation conditions in the file.\n")
        self.help_check_exit()

        if os.path.exists(f"./ConditionList_{self.Nmol}{self.mol_pos}.txt"):
            self.messages.append(f"\t>>> ConditionList_{self.Nmol}{self.mol_pos}.txt: Exist.")
            self.HelpList.append(False)
        elif os.path.exists(f"./InitialCondition_{self.Nmol}{self.mol_pos}.txt"):
            self.messages.append(f"\t>>> InitialCondition_{self.Nmol}{self.mol_pos}.txt: Exist.")
            self.HelpList.append(False)
        else:
            self.messages.append(
                f"\t{Color.RED}>>> ConditionList_{self.Nmol}{self.mol_pos}.txt: Not exist.{Color.RESET}")
            self.messages.append(
                f"\t{Color.RED}>>> InitialCondition_{self.Nmol}{self.mol_pos}.txt: Not exist.{Color.RESET}")
            self.HelpList.append(True)
            with open(f"./InitialCondition_{self.Nmol}{self.mol_pos}.txt", "w") as f:
                if "2mol" in self.Nmol:
                    f.write(Stereotyped.InitialCondition_2mol_Temp)
                elif "3mol" in self.Nmol:
                    f.write(Stereotyped.InitialCondition_3mol_Temp)
            self.messages.append(
                f"\t>>> InitialCondition_{self.Nmol}{self.mol_pos}.txt: {Color.GREEN}Created.{Color.RESET}\n"
                f"\t>>> Please set the initial conditions in the file.\n")

        if "3mol" in self.Nmol:
            if os.path.exists(f"./{self.MaterName}_2mol_min.txt"):
                self.messages.append(f"\t>>> {self.MaterName}_2mol_min.txt: Exist.")
                self.HelpList.append(False)
            else:
                self.messages.append(f"\t{Color.RED}>>> {self.MaterName}_2mol_min.txt: Not exist.{Color.RESET}")
                self.HelpList.append(True)
                self.messages.append(f"\t>>> Please set the minimum energy structure of 2mol in the file.\n")
        else:
            pass

        self.help_check_exit()
        return None

    def mkFirstCondition(self):
        which, RefLines = "", []
        if "2mol" in self.Nmol:
            which = "D_edge"
            dev = 0.1
        elif "3mol" in self.Nmol:
            which = "D_faceon"
            RefLines = self.getRefLines(f"./{self.MaterName}_2mol_min.txt")
            dev = 0.2
            os.makedirs(self.tcalpath, exist_ok=True)
        else:
            which = "D_edge"
            dev = 0.1

        os.makedirs(self.dirpath, exist_ok=True)

        self.mkConditionFile(which, RefLines, dev)
        return RefLines, dev, which

    @staticmethod
    def getRefLines(FileName):
        with open(FileName, "r") as f:
            lines = f.readlines()
        RefLines = []
        for line in lines:
            contents = line.strip().split()
            try:
                float(contents[0])
                RefLines.append(line)
            except ValueError:
                pass
        return RefLines

    def mkConditionFile(self, which, RefLines, dev):
        print(f"\n{Color.GREEN}Create a Condition file...{Color.RESET}")
        if os.path.exists(f"./ConditionList_{self.Nmol}{self.mol_pos}.txt"):
            self.messages.append(f"\t>>> ./ConditionList_Tilt_{self.Nmol}{self.mol_pos}.txt: Exist.\n"
                                 f"\t>>> Performs calculations based on conditions that exist.")
        else:
            self.messages.append(f"\t>>> ./ConditionList_Tilt_{self.Nmol}{self.mol_pos}.txt: Not exist.\n"
                                 f"\t>>> {Color.GREEN}Create a new ConditionList.txt{Color.RESET} "
                                 f"based on the InitialCondition.")
            self.help_check_exit()
            NewConditions = []
            with open(f"./InitialCondition_{self.Nmol}{self.mol_pos}.txt", "r") as f:
                lines = f.readlines()
            for line in lines:
                if "2mol" in self.Nmol:
                    Other = 0
                    Val = float(line.strip())
                    RefValues = ["na", "na"]
                else:
                    Other = float(line.strip().split()[0])
                    Val = float(line.strip().split()[1])
                    RefValues = [RefLines[0].strip().split()[0], RefLines[0].strip().split()[1]]

                NewCondition = self.mkNewCondition(Val - (2 * dev), Other, RefValues, which)
                NewConditions.append(NewCondition)
                NewCondition = self.mkNewCondition(Val - (1 * dev), Other, RefValues, which)
                NewConditions.append(NewCondition)
                NewCondition = self.mkNewCondition(Val, Other, RefValues, which)
                NewConditions.append(NewCondition)
                NewCondition = self.mkNewCondition(Val + (1 * dev), Other, RefValues, which)
                NewConditions.append(NewCondition)
                NewCondition = self.mkNewCondition(Val + (2 * dev), Other, RefValues, which)
                NewConditions.append(NewCondition)
            NewConditions = list(set(NewConditions))
            NewConditions.sort()
            with open(f"./ConditionList_{self.Nmol}{self.mol_pos}.txt", "w") as f:
                f.writelines(NewConditions)
            self.messages.append(f"\t>>> ./ConditionList_Tilt_{self.Nmol}{self.mol_pos}.txt: "
                                 f"{Color.GREEN}Created.{Color.RESET}\n")
        self.help_check_exit()
        return None

    def mkNewCondition(self, Val, Other, RefValues, which):
        Val, Other = int(Val * 100), int(Other * 100)
        Other = int(Other * 100)
        if "2mol" in self.Nmol:
            D_Edge = Val
            NewCondition = f"{D_Edge}\n"
        elif "3mol" in self.Nmol:
            if which == "D_edge":
                D_Edge = Val
                D_Faceon = int(RefValues[1] * 100)
            elif which == "D_faceon":
                D_Edge = int(RefValues[0] * 100)
                D_Faceon = Val
            else:
                D_Edge = 0
                D_Faceon = 0
            NewCondition = f"{Other}_{D_Edge}_{D_Faceon}\n"
        else:
            NewCondition = ""
        print(f"\t\t{NewCondition.strip()}")
        return NewCondition

    def Most_Stable_Structure_Search(self):
        self.getTemporaryStructure()

    def getTemporaryStructure(self):
        judge = False
        while not judge:
            qsubList = []
            Conditions = self.getCondition(f"./ConditionList_{self.Nmol}{self.mol_pos}.txt")

            with open(f"{self.dirpath}/G.sh", "w") as f:
                f.write(Stereotyped.Sh_txt)

            for Condition in Conditions:
                File_Name = f"{self.MaterName}_{self.Nmol}{self.mol_pos}_{Condition}"
                if os.path.exists(f"{self.dirpath}/{File_Name}.log"):
                    pass
                else:
                    self.mkFiles(Condition, File_Name)

    @staticmethod
    def getCondition(FileName):
        with open(FileName, "r") as f:
            lines = f.readlines()
        Conditions = []
        for line in lines:
            Conditions.append(line.strip())
        return Conditions

    def mkFiles(self, Condition, File_Name):
        CHK_Name = f"{File_Name}.chk"
        GJF_Name = f"{File_Name}.gjf"
        SH_Name = f"G-{self.Operator}_{Condition}.sh"

        ConditionList = Condition.strip().split("_")
        Edge_Axis, Face_Axis, Other_Axis, Matrix_Mol3, rotate = self.Axis_Setting_BW(self.messages, self.HelpList)

    def Axis_Setting_BW(self, messages, HelpList):
        try:
            with open(f"./CalcSetting_BW.txt", "r") as f:
                lines = f.readlines()
            HelpList.append(False)
        except FileNotFoundError:
            print(f"{Color.RED}CalcSetting_BW.txt: Not exist.{Color.RESET}")
            with open(f"./CalcSetting_BW.txt", "w") as f:
                f.write(Stereotyped.CalcSet_BW_template)
            print(f"{Color.GREEN}CalcSetting_BW.txt: Created.{Color.RESET}")
            HelpList.append(True)
        self.help_check_exit()

        params = {"edge axis": "", "faceon axis": "", "mol3 other_transition": "", "other axis": ""}
        for line in lines[:4]:
            key = line.split(":")[0].strip().lower()
            if key in params:
                params[key] = line.split(":")[1].strip()

        if params["edge axis"] == "x":
            if params["faceon axis"] == "z":
                params["other axis"] = "y"
            elif params["faceon axis"] == "y":
                params["other axis"] = "z"
        if params["edge axis"] == "y":
            if params["faceon axis"] == "x":
                params["other axis"] = "z"
            elif params["faceon axis"] == "z":
                params["other axis"] = "x"
        if params["edge axis"] == "z":
            if params["faceon axis"] == "x":
                params["other axis"] = "y"
            elif params["faceon axis"] == "y":
                params["other axis"] = "x"

        axis1 = params["edge axis"]
        axis2 = params["faceon axis"]
        axis3 = params["other axis"]
        rotate = f"{axis1}{axis2}{axis3}"

        Mol3_Other = self.mkDirection(params["mol3 other_transition"], params["other axis"],
                                      "Mol3 Other Transition", messages, HelpList)

        Matrix_Mol3 = np.array(Mol3_Other)

        return params["edge axis"], params["faceon axis"], params["other axis"], Matrix_Mol3, rotate

    def mkDirection(self, axis_direction, input_axis, axis_name, messages, HelpList):
        """
        Make the direction
        :param axis_direction:
        :param input_axis:
        :param axis_name:
        :param messages:
        :param HelpList:
        :return:
        """
        try:
            axis_direction = float(axis_direction)
        except ValueError:
            messages.append(f"{Color.RED}Error: Invalid {axis_name} direction.{Color.RESET}")
            HelpList.append(True)
        help_check_exit(messages, HelpList)
        directions = {"x": [axis_direction, 0.0, 0.0], "y": [0.0, axis_direction, 0.0], "z": [0.0, 0.0, axis_direction]}
        if input_axis not in directions:
            messages.append(f"{Color.RED}Error: Invalid {axis_name} axis.{Color.RESET}")
            HelpList.append(True)
        self.help_check_exit()
        return np.array(directions[input_axis])


if __name__ == "__main__":
    main()
