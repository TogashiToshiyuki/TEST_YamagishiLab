#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import datetime
import functools
import glob
import math
import os
import random
import subprocess
import sys
import time

import numpy as np

print = functools.partial(print, flush=True)


def main():
    # program start
    # View Program Overview
    print(Stereotyped.ProgramAbst)

    # argument parser
    args, before = arg_parser()

    # create an instance of the BrickWork class
    bw = BrickWork(args, before)

    # Create the Initial condition
    bw.mkInitialCondition()

    # Search for the most stable structure
    MinConditions = bw.Most_Stable_Search()

    # Calculate TIs
    bw.Calculate_TI(MinConditions)

    # Save the results
    bw.Result_Data_set(before)

    # program end
    print(f"{Color.GREEN}"
          f"************************* ALL PROCESSES END *************************"
          f"{Color.RESET}\n")
    return


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
                           "Initial offset Edge: 2.6\n"
                           "Initial offset Faceon: 3.8\n"
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
    """
    Constants
    """
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


# arg_parser
def arg_parser():
    """
    Argument parser
    :return: args, before
    :rtype: argparse.Namespace, float
    """
    before = time.time()
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('MaterNameXYZ',
                        help="Molecular_Name.xyz")
    parser.add_argument('--debug', '-d', '-D',
                        help="Start the programme in Debug mode.",
                        action="store_false")
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
    def __init__(self, args, before):
        self.MaterName = args.MaterNameXYZ[:-4]
        self.Flag_xyz = args.xyz
        self.Debug = args.debug
        self.chk = args.chk

        with open(f"{self.MaterName}.xyz", "r") as f:
            self.NinMol = f.readline()
            f.readline()
            self.AtomList = f.readlines()

        # Select the structure of the molecule
        if not self.chk:
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
            self.mol_pos = "p1"

        self.dirpath = f"./{self.MaterName}_3mol{self.mol_pos}"
        self.tcalpath = f"./{self.MaterName}_3mol{self.mol_pos}_tcal"
        if not self.chk:
            os.makedirs(self.dirpath, exist_ok=True)
        else:
            pass
        self.calculation_tcal_flag = args.tcal

        # Retrieve the operator name
        if not self.chk:
            Operator = input(f"\n{Color.GREEN}Retrieve the operator name.{Color.RESET}\n"
                             "\t>>> The name entered here will be used to identify the operator.\n"
                             f"\tPlease enter the operator.\n"
                             f"\t{Color.GREEN}>>> {Color.RESET}")
            if Operator == "":
                Operator = "ONE"
            else:
                pass
            print("")
        else:
            Operator = "ONE"
        self.Operator = Operator

        self.messages, self.HelpList = [], []

        self.File_Set_Check()

        # Check the axis setting
        try:
            with open(f"./CalcSetting_BW.txt", "r") as f:
                lines = f.readlines()
            self.HelpList.append(False)
        except FileNotFoundError:
            print(f"{Color.RED}CalcSetting_BW.txt: Not exist.{Color.RESET}")
            with open(f"./CalcSetting_BW.txt", "w") as f:
                f.write(Stereotyped.CalcSet_BW_template)
            print(f"{Color.GREEN}CalcSetting_BW.txt: Created.{Color.RESET}")
            self.HelpList.append(True)
        self.help_check_exit()
        params = {"edge axis": "",
                  "faceon axis": "",
                  "mol3 other_transition": "",
                  "other axis": "",
                  "initial offset edge": "",
                  "initial offset faceon": ""}
        for line in lines[:6]:
            key = line.split(":")[0].strip().lower()
            if key in params:
                params[key] = line.split(":")[1].strip()
        params["other axis"] = {
            ("x", "z"): "y",
            ("x", "y"): "z",
            ("y", "x"): "z",
            ("y", "z"): "x",
            ("z", "x"): "y",
            ("z", "y"): "x"
        }.get((params["edge axis"], params["faceon axis"]))
        axis1 = params["edge axis"]
        axis2 = params["faceon axis"]
        axis3 = params["other axis"]
        rotate = f"{axis1}{axis2}{axis3}"
        axis_pairs = [(axis1, axis2, "Edge", "Faceon"),
                      (axis1, axis3, "Edge", "Other"),
                      (axis2, axis3, "Faceon", "Other")]
        for a, b, name1, name2 in axis_pairs:
            if a == b:
                self.messages.append(f"{Color.RED}Error: Invalid axis setting.{Color.RESET}\n"
                                     f"\t>>> The '{Color.RED}{name1}{Color.RESET}' axis and "
                                     f"the '{Color.RED}{name2}{Color.RESET}' axis cannot be the same.\n"
                                     f"\t>>> The current setting Edge:   '{Color.RED}{a}{Color.RESET}'.\n"
                                     f"\t>>> The current setting Faceon: '{Color.RED}{b}{Color.RESET}'.")
                self.HelpList.append(True)
        if rotate not in {"xyz", "xzy", "yxz", "yzx", "zxy", "zyx"}:
            self.messages.append(f"{Color.RED}Error: Invalid axis setting.{Color.RESET}\n"
                                 f"\t>>> rotate must be one of 'xyz', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx'.\n"
                                 f"\t>>> The current setting is '{Color.RED}{rotate}{Color.RESET}'.")
            self.HelpList.append(True)
        self.help_check_exit()
        self.Edge_Axis = params["edge axis"]
        self.Faceon_Axis = params["faceon axis"]
        self.Other_Axis = params["other axis"]
        self.initial_offset_Edge = float(params["initial offset edge"])
        self.initial_offset_Faceon = float(params["initial offset faceon"])
        self.rotate = rotate
        self.Mol3_Other = self.mkDirection(
            params["mol3 other_transition"], params["other axis"], "Mol3 Other Transition")
        Debug_message_List = [f"Edge Axis: {self.Edge_Axis}",
                              f"Faceon Axis: {self.Faceon_Axis}",
                              f"Other Axis: {self.Other_Axis}",
                              f"Rotate: {rotate}",
                              f"Mol3 Other Transition: {self.Mol3_Other}",
                              f"Initial offset Edge: {self.initial_offset_Edge}",
                              f"Initial offset Faceon: {self.initial_offset_Faceon}"]
        self.debug_message(Debug_message_List)

        self.mkCheckFile(args, before)

    # Displays messages and exits if unexpected behavior is detected.
    def help_check_exit(self):
        """
        Checks for unexpected behaviors and displays accumulated messages.

        This function performs two main tasks:

        1. It displays accumulated messages stored in `self.message`, allowing
           for effective message management within the program.
        2. It checks the `HelpList` attribute for any instances of `True`.
           If `True` is found, indicating an unexpected behavior, the program
           terminates after displaying an error message.

        :returns: None
        :raises SystemExit: Terminates the program if `True` is found in
                            `HelpList`, signaling an abnormal end due to an
                            unexpected behavior.
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

    def debug_message(self, Debug_message_List):
        """
        Show messages in debug mode.
        :param Debug_message_List:
        :return:
        """
        if self.Debug:
            print(f"\n{Color.RED}Debug Message{Color.RESET}")
            for message in Debug_message_List:
                print(message)
            print(f"{Color.RED}End Debug Message{Color.RESET}\n")
        else:
            pass
        return None

    # checks for required files, halts if missing, and creates templates as needed.
    def File_Set_Check(self):
        """
        Verify the existence of required files and create templates if missing.

        This function checks for the existence of specific required files for
        further program execution. If any required file is missing, the program
        will halt, display the missing files, and optionally create template
        files as needed. The created templates help guide users in setting up
        necessary conditions for calculations.

        :return: None
        :raises FileNotFoundError: If essential files are missing and the user
            needs to set conditions in the created templates before continuing.
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
                                 f"\t>>> {Color.GREEN}"
                                 f"Please set the calculation conditions in the "
                                 f"{Color.UNDERLINE}'CalcSetting_BW.txt.'{Color.RESET}")
        self.help_check_exit()

        if not self.chk:
            if os.path.exists(f"./ConditionList_3mol{self.mol_pos}.txt"):
                self.messages.append(f"\t>>> ConditionList_3mol{self.mol_pos}.txt: {Color.GREEN}Exist.{Color.RESET}")
                self.HelpList.append(False)
            elif os.path.exists(f"./InitialCondition_3mol{self.mol_pos}.txt"):
                self.messages.append(f"\t>>> InitialCondition_3mol{self.mol_pos}.txt: {Color.GREEN}Exist.{Color.RESET}")
                self.HelpList.append(False)
            else:
                self.messages.append(
                    f"\t>>> {Color.RED}ConditionList_3mol{self.mol_pos}.txt: Does not Exist.{Color.RESET}\n"
                    f"\t>>> {Color.RED}InitialCondition_3mol{self.mol_pos}.txt: Does not Exist.{Color.RESET}\n")
                self.HelpList.append(True)
                with open(f"./InitialCondition_3mol{self.mol_pos}.txt", "w") as f:
                    f.write(Stereotyped.InitialCondition_3mol_Temp)
                self.messages.append(f"\t>>> InitialCondition_3mol{self.mol_pos}.txt: "
                                     f"{Color.GREEN}Created.{Color.RESET}\n"
                                     f"\t>>> {Color.GREEN}"
                                     f"Please set the calculation conditions in the "
                                     f"{Color.UNDERLINE}'InitialCondition_3mol{self.mol_pos}.txt'.{Color.RESET}")
            self.help_check_exit()
        else:
            pass
        return None

    def mkCheckFile(self, args, before):
        if self.chk:
            pass
        else:
            return None
        print(f"{Color.RED}"
              f"********** Starts generating files for structural verification **********\n"
              f"{Color.RESET}")

        filepath = f"./StructCheck-{self.MaterName}"
        os.makedirs(filepath, exist_ok=True)
        Temp_SHs = []

        with open(f"{filepath}/G.sh", "w") as f:
            f.write(Stereotyped.Sh_txt)
        Temp_SHs.append(f"G.sh")
        if not args.manual:
            print(f"Automatically generate files for structural verification...")
            self.mol_pos = "p1"
            Temp_SHs.append(self.mkFiles("0_1800_400", filepath).replace("qsub ", ""))
            self.messages.append(f"\t>>> {self.MaterName}_3mol{self.mol_pos}_0_1800_400.gjf: Created.")

            self.mol_pos = "p2"
            self.mkFiles("0_1800_400", filepath)
            self.messages.append(f"\t>>> {self.MaterName}_3mol{self.mol_pos}_0_1800_400.gjf: Created.")

            self.mol_pos = "p3"
            self.mkFiles("0_1800_400", filepath)
            self.messages.append(f"\t>>> {self.MaterName}_3mol{self.mol_pos}_0_1800_400.gjf: Created.")
            self.help_check_exit()

        print("\nDelete unnecessary files...")
        Temp_SHs = list(set(Temp_SHs))
        for SH in Temp_SHs:
            try:
                subprocess.run(["rm", SH], timeout=10, cwd=filepath, check=True, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError:
                self.messages.append(f"\t>>> {Color.RED}{SH}: Failed to delete.{Color.RESET}")
                self.HelpList.append(True)
            else:
                self.messages.append(f"\t>>> {SH}: {Color.GREEN}Deleted.{Color.RESET}")
                self.HelpList.append(False)
        self.help_check_exit()
        if self.Debug:
            print("*************** Debug Finished!!! ***************")

        after = time.time()

        # End of the program
        print(f"\n"
              f"Elapsed Time: {(after - before):.0f} s\n{Color.GREEN}"
              f"************************* ALL PROCESSES END *************************"
              f"{Color.RESET}\n")
        exit()

    def mkInitialCondition(self):
        # Create the first condition
        self.messages.append(f"\n{Color.GREEN}Creating the condition for first cycle...{Color.RESET}")
        if os.path.exists(f"./ConditionList_3mol{self.mol_pos}.txt"):
            self.messages.append(f"\t>>> ConditionList_3mol{self.mol_pos}.txt: {Color.GREEN}Exist.{Color.RESET}\n"
                                 f"\t>>> Performs calculations based on conditions that exist.")
        else:
            Coordinate_X, Coordinate_Y, Coordinate_Z = [], [], []
            for Atom in self.AtomList:
                List = Atom.strip().split()
                Coordinate_X.append(float(List[1]))
                Coordinate_Y.append(float(List[2]))
                Coordinate_Z.append(float(List[3]))
            List = {
                "x": [max(Coordinate_X), min(Coordinate_X)],
                "y": [max(Coordinate_Y), min(Coordinate_Y)],
                "z": [max(Coordinate_Z), min(Coordinate_Z)]
            }

            Edge_Max = self.transform_number(List[self.Edge_Axis][0])
            Edge_Min = self.transform_number(List[self.Edge_Axis][1])
            Faceon_Max = self.transform_number(List[self.Faceon_Axis][0])
            Faceon_Min = self.transform_number(List[self.Faceon_Axis][1])
            Debug_Message_List = [
                f"Edge Max: {Edge_Max}",
                f"Edge Min: {Edge_Min}",
                f"Faceon Max: {Faceon_Max}",
                f"Faceon Min: {Faceon_Min}"
            ]
            self.debug_message(Debug_Message_List)
            First_Edge = (self.transform_number((Edge_Max - Edge_Min) / 2) +
                          Edge_Max + self.initial_offset_Edge)
            First_Faceon = (self.transform_number((Faceon_Max - Faceon_Min) / 2)
                            + Faceon_Max + self.initial_offset_Faceon)
            Debug_Message_List = [
                f"First Edge: {First_Edge}",
                f"First Faceon: {First_Faceon}"
            ]
            self.debug_message(Debug_Message_List)
            NewConditions, dev = [], 0.2
            self.messages.append(f"\t>>> ConditionList_3mol{self.mol_pos}.txt: Does not Exist.\n"
                                 f"\t>>> {Color.GREEN}Create a new ConditionList.txt{Color.RESET} "
                                 f"based on the InitialCondition.")
            self.help_check_exit()
            with open(f"./InitialCondition_3mol{self.mol_pos}.txt", "r") as f:
                lines = f.readlines()
            for line in lines:
                Other = float(line.strip().split()[0])
                Debug_Message_List = [f"Edge: {First_Edge}, Faceon: {First_Faceon}, Other: {Other}"]
                self.debug_message(Debug_Message_List)
                NewConditions.append(
                    self.mkNewCondition(Other, First_Edge - (2 * dev), "Edge", [First_Edge, First_Faceon]))
                NewConditions.append(
                    self.mkNewCondition(Other, First_Edge - dev, "Edge", [First_Edge, First_Faceon]))
                NewConditions.append(
                    self.mkNewCondition(Other, First_Edge, "Edge", [First_Edge, First_Faceon]))
                NewConditions.append(
                    self.mkNewCondition(Other, First_Edge + dev, "Edge", [First_Edge, First_Faceon]))
                NewConditions.append(
                    self.mkNewCondition(Other, First_Edge + (2 * dev), "Edge", [First_Edge, First_Faceon]))
            NewConditions = sorted(set(NewConditions))
            with open(f"./ConditionList_3mol{self.mol_pos}.txt", "w") as f:
                for Condition in NewConditions:
                    f.write(f"{Condition}")
            self.messages.append(f"\t>>> ConditionList_3mol{self.mol_pos}.txt: {Color.GREEN}Created.{Color.RESET}")
        self.help_check_exit()
        return None

    @staticmethod
    def transform_number(number):
        if number % 1 >= 0.5:
            return math.ceil(number * 2) / 2
        else:
            return round(number * 2) / 2

    @staticmethod
    def mkNewCondition(Other, Val, which, RefValues):
        if which == "Edge":
            Edge = round(float(Val), 2)
            Faceon = float(RefValues[1])
        elif which == "Faceon":
            Edge = float(RefValues[0])
            Faceon = round(float(Val), 2)
        else:
            Edge = RefValues[0]
            Faceon = RefValues[1]
        NewCondition = f"{int(Other * 100)}_{int(round(Edge * 100))}_{int(round(Faceon * 100))}\n"
        print(f"\t\t{NewCondition.strip()}")
        return NewCondition

    def Most_Stable_Search(self):
        self.getTemporaryStructure("Edge", 0.2)

        temp_structure = []
        # 0.2
        MostStable = False
        dev = 0.2
        while not MostStable:
            RefLines = self.getRefLines(f"./{self.MaterName}_3mol{self.mol_pos}_min.txt")
            temp_structure.append(RefLines)
            self.mkCycleConditions("Faceon", dev)
            self.getTemporaryStructure("Faceon", dev)
            self.mkCycleConditions("Edge", dev)
            self.getTemporaryStructure("Edge", dev)
            RefLines = self.getRefLines(f"./{self.MaterName}_3mol{self.mol_pos}_min.txt")
            MostStable = self.CompareStructures(RefLines, temp_structure[-1])

        # 0.1
        print(f"\n{Color.GREEN}**********\nTransition in 0.1-Å increments.\n{Color.RESET}")
        MostStable = False
        dev = 0.1
        while not MostStable:
            RefLines = self.getRefLines(f"./{self.MaterName}_3mol{self.mol_pos}_min.txt")
            temp_structure.append(RefLines)
            self.mkCycleConditions("Faceon", dev)
            self.getTemporaryStructure("Faceon", dev)
            self.mkCycleConditions("Edge", dev)
            self.getTemporaryStructure("Edge", dev)
            RefLines = self.getRefLines(f"./{self.MaterName}_3mol{self.mol_pos}_min.txt")
            MostStable = self.CompareStructures(RefLines, temp_structure[-1])

        # 0.05
        print(f"\n{Color.GREEN}**********\nTransition in 0.05-Å increments.\n{Color.RESET}")
        MostStable = False
        dev = 0.05
        while not MostStable:
            RefLines = self.getRefLines(f"./{self.MaterName}_3mol{self.mol_pos}_min.txt")
            temp_structure.append(RefLines)
            self.mkCycleConditions("Faceon", dev)
            self.getTemporaryStructure("Faceon", dev)
            self.mkCycleConditions("Edge", dev)
            self.getTemporaryStructure("Edge", dev)
            RefLines = self.getRefLines(f"./{self.MaterName}_3mol{self.mol_pos}_min.txt")
            MostStable = self.CompareStructures(RefLines, temp_structure[-1])
        MinConditions = self.getMinCondition()

        print(f"\n{Color.GREEN}The most stable structure has been found at '{len(MinConditions)}' steps.{Color.RESET}")
        print("\tMost stable Conditions:")
        for i in range(len(MinConditions)):
            print(f"\t\t{MinConditions[i]}")
        if len(temp_structure) < 500:
            with open(f"{self.MaterName}_3mol{self.mol_pos}_mins.hist", "w") as f:
                for i in range(len(temp_structure)):
                    f.write(f"***** Structures after the '{i + 1}'th cycle *****\n")
                    Lines = temp_structure[i]
                    f.writelines(Lines)
        print(f"\n\t>>> {self.MaterName}_3mol{self.mol_pos}_mins.hist: Created.")
        return MinConditions

    def getTemporaryStructure(self, which, dev):
        judge = False
        while not judge:
            qsubList = []
            Conditions = self.getConditions(f"./ConditionList_3mol{self.mol_pos}.txt")
            with open(f"{self.dirpath}/G.sh", "w") as f:
                f.write(Stereotyped.Sh_txt)

            for Condition in Conditions:
                if os.path.exists(f"{self.dirpath}/{self.MaterName}_3mol{self.mol_pos}_{Condition}.log"):
                    pass
                else:
                    qsubList.append(self.mkFiles(Condition, self.dirpath))
            self.job_submission(qsubList, which, self.dirpath)
            if self.Debug:
                pass
            else:
                self.rmWildCards(f"{self.dirpath}/*.sh*")
                self.rmWildCards(f"{self.dirpath}/*.chk")
            print("\n**********\nReading Data...\n")
            self.readEnergy()
            judge = self.mkNewConditionLists(which, dev)
        return None

    @staticmethod
    def getConditions(FileName):
        Condition = []
        with open(FileName, "r") as f:
            lines = f.readlines()
        for line in lines:
            Condition.append(line.strip())
        return Condition

    def mkFiles(self, Condition, dirpath):
        FileName = f"{self.MaterName}_3mol{self.mol_pos}_{Condition}"
        CHK_FileName = f"{FileName}.chk"
        GJF_FileName = f"{FileName}.gjf"
        SH_FileName = f"G-{self.Operator}_{Condition}.sh"

        Condition = Condition.strip().split("_")
        Matrix_Mol3 = {
                          "x": [int(Condition[0]) / 100, 0, 0],
                          "y": [0, int(Condition[0]) / 100, 0],
                          "z": [0, 0, int(Condition[0]) / 100]
                      }[self.Edge_Axis] + self.Mol3_Other
        Transitions = self.mkTransition(self.Edge_Axis, self.Faceon_Axis, Condition[1], Condition[2], Matrix_Mol3)

        Angles = {
            "p1": {
                "x": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                "y": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                "z": [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            },
            "p2": {
                "x": [[180, 0, 0], [180, 0, 0], [180, 0, 0]],
                "y": [[0, 180, 0], [0, 180, 0], [0, 180, 0]],
                "z": [[0, 0, 180], [0, 0, 180], [0, 0, 180]]
            },
            "p3": {
                "x": [[180, 0, 0], [180, 0, 0], [180, 0, 0]],
                "y": [[0, 180, 0], [0, 180, 0], [0, 180, 0]],
                "z": [[0, 0, 180], [0, 0, 180], [0, 0, 180]]
            }
        }.get(self.mol_pos).get(self.Other_Axis)

        Element, Mol1_pos, Mol2_pos, Mol3_pos = self.mkAtomList(Angles[0], Angles[1], Angles[2], Transitions)

        with open(f"Temp_Header_3mol.txt", "w") as temp_header:
            temp_header.write(Stereotyped.Header_3mol)
        with open(f"Temp_Header_3mol.txt", "r") as temp_header:
            Headers = temp_header.readlines()
        os.remove(f"Temp_Header_3mol.txt")
        Headers[3] = f"%chk={CHK_FileName}\n"

        self.write_gjf_file(f"{dirpath}/{FileName}.gjf",
                            Headers, Element, Mol1_pos, Mol2_pos, Mol3_pos)
        if self.Flag_xyz:
            self.write_xyz_file(f"{dirpath}/{FileName}.xyz",
                                Element, Mol1_pos, Mol2_pos, Mol3_pos)
            self.messages.append(f"\t>>> {FileName}.xyz: Created.")

        with open(f"{dirpath}/G.sh", "r") as f:
            lines = f.readlines()
        lines[12] = f"g16 {GJF_FileName}\n"
        with open(f"{dirpath}/{SH_FileName}", "w") as f:
            f.writelines(lines)
        qsub_temp = f"qsub {SH_FileName}"
        return qsub_temp

    def mkDirection(self, axis_direction, input_axis, axis_name):
        """
        Make the direction
        :param axis_direction:
        :param input_axis:
        :param axis_name:
        :return:
        """
        try:
            axis_direction = float(axis_direction)
        except ValueError:
            self.messages.append(f"{Color.RED}Error: Invalid {axis_name} direction.{Color.RESET}")
            self.HelpList.append(True)
        self.help_check_exit()
        directions = {"x": [axis_direction, 0.0, 0.0], "y": [0.0, axis_direction, 0.0], "z": [0.0, 0.0, axis_direction]}
        if input_axis not in directions:
            self.messages.append(f"{Color.RED}Error: Invalid {axis_name} axis.{Color.RESET}")
            self.HelpList.append(True)
        self.help_check_exit()
        return np.array(directions[input_axis])

    def mkTransition(self, Axis_Edge, Axis_Faceon, Direction_Edge, Direction_Faceon, Matrix_Mol3):
        Edge_transl = self.mkDirection(Direction_Edge, Axis_Edge, "Edge Direction")
        Faceon_transl = self.mkDirection(Direction_Faceon, Axis_Faceon, "Faceon Direction")

        return [Edge_transl, Faceon_transl, Matrix_Mol3]

    def mkAtomList(self, Mol1_Angles, Mol2_Angles, Mol3_Angles, Transitions):
        """
        Make the Atom List
        :param Mol1_Angles:
        :param Mol2_Angles:
        :param Mol3_Angles:
        :param Transitions:
        :return:
        """
        Mol1, Mol2, Mol3, Elist = [], [], [], []
        for Atom in self.AtomList:
            Contents = Atom.split()
            Elist.append(Contents.pop(0))
            Position = np.array(list(map(float, Contents)))
            atm_m1 = self.Rotate(Position, Mol1_Angles[0], Mol1_Angles[1], Mol1_Angles[2])
            Mol1.append(atm_m1)
            atm_m2 = self.Rotate(Position, Mol1_Angles[0], Mol2_Angles[1], Mol2_Angles[2]) + (Transitions[0] / 100)
            Mol2.append(atm_m2)
            atm_m3 = (self.Rotate(Position, Mol3_Angles[0], Mol3_Angles[1], Mol3_Angles[2]) + Transitions[0] / 200
                      + Transitions[1] / 100 + Transitions[2])
            Mol3.append(atm_m3)
        return Elist, Mol1, Mol2, Mol3

    def Rotate(self, Current, Tx, Ty, Tz):
        Tx, Ty, Tz = map(math.radians, [Tx, Ty, Tz])
        Rx = np.array([[1, 0, 0], [0, math.cos(Tx), -math.sin(Tx)], [0, math.sin(Tx), math.cos(Tx)]])
        Ry = np.array([[math.cos(Ty), 0, math.sin(Ty)], [0, 1, 0], [-math.sin(Ty), 0, math.cos(Ty)]])
        Rz = np.array([[math.cos(Tz), -math.sin(Tz), 0], [math.sin(Tz), math.cos(Tz), 0], [0, 0, 1]])
        rotation_matrices = {'x': Rx, 'y': Ry, 'z': Rz}
        R = np.eye(3)
        for axis in self.rotate:
            R = np.dot(rotation_matrices[axis], R)
        return np.dot(R, Current)

    def write_gjf_file(self, filename, headers, elements, positions1, positions2=None, positions3=None):
        with open(filename, "w") as file:
            for header in headers:
                file.write(header)
            for elem, pos in zip(elements, positions1):
                file.write(f" {elem:<2}  {self.format_coordinate(pos)}  1\n")
            if positions2 is not None:
                for elem, pos in zip(elements, positions2):
                    file.write(f" {elem:<2}  {self.format_coordinate(pos)}  2\n")
            if positions3 is not None:
                for elem, pos in zip(elements, positions3):
                    file.write(f" {elem:<2}  {self.format_coordinate(pos)}  3\n")
            file.write("\n")
        return

    def write_xyz_file(self, filename, elements, positions1, positions2=None, positions3=None):
        with open(filename, "w") as file:
            if positions3 is None:
                file.write(f"{len(elements) * 2}\n")
            else:
                file.write(f"{len(elements) * 3}\n")
            file.write("00000001\n")
            for elem, pos in zip(elements, positions1):
                file.write(f" {elem:<2}  {self.format_coordinate(pos)}\n")
            if positions2 is not None:
                for elem, pos in zip(elements, positions2):
                    file.write(f" {elem:<2}  {self.format_coordinate(pos)}\n")
            if positions3 is not None:
                for elem, pos in zip(elements, positions3):
                    file.write(f" {elem:<2}  {self.format_coordinate(pos)}\n")
            file.write("\n")
        return

    @staticmethod
    def format_coordinate(coord):
        """
        Format the coordinate
        :param coord:
        :return:
        """
        x, y, z = map(float, [coord[0], coord[1], coord[2]])
        x = f'{x:.10f}'.rjust(15, ' ')
        y = f'{y:.10f}'.rjust(15, ' ')
        z = f'{z:.10f}'.rjust(15, ' ')
        return f'{x} {y} {z}'

    def job_submission(self, qsubList, which, dirpath):
        """
        Submit the job
        :param qsubList:
        :param which:
        :param dirpath:
        :return:
        """
        print("\n**********\nJobs are submitting...")
        if len(qsubList) == 0:
            self.messages.append(f"\t>>> Job was not submitted.\n"
                                 f"\t>>> {Color.GREEN}Calculations with the conditions might be finished.{Color.RESET}")
            self.help_check_exit()
        else:
            for qsub in qsubList:
                qsub = qsub.split()
                subprocess.run(qsub, cwd=dirpath)

            Wait_minutes = 2
            MyJobIDList = self.My_JobIDList(qsubList)
            MyJobIDList.sort()
            RunningJobIDList = self.Running_JobIDList()

            Flag, wait_job_count = self.check_jobs(RunningJobIDList, MyJobIDList)
            start_time = datetime.datetime.now()

            if Flag:
                formated_ST = start_time.strftime("%m/%d %H:%M:%S")
                if which == "Edge":
                    term = "the stable distance in Edge direction"
                elif which == "Faceon":
                    term = "the stable distance in Faceon direction"
                elif which == "tcal":
                    term = "tcal"
                else:
                    term = ""
                print(f"{Color.GREEN}\t>>> '{int(len(MyJobIDList))}' calculations for '{term}' was submitted!!"
                      f" {Color.RESET}at {formated_ST}")
                print(f"\n"
                      f"{Color.GREEN}Wait until jobID {MyJobIDList[-1]}!!{Color.RESET}\n"
                      f"start ID: {RunningJobIDList[0]}\n"
                      f"\n")
            else:
                pass

            while Flag:
                Flag, wait_job_count = self.check_jobs(RunningJobIDList, MyJobIDList)
                formated_NOW, elapsed_time = self.getElapsedTime(start_time)
                Now_Time = datetime.datetime.now()
                formated_end_time = (
                    (Now_Time + datetime.timedelta(minutes=(wait_job_count * Wait_minutes))).strftime("%m/%d %H:%M:%S"))
                sys.stdout.write(
                    "\033[1F\033[G%s" %
                    f"\t{formated_NOW} ({elapsed_time} min. passed): Job {RunningJobIDList[0]} is in progress.       \n"
                    f"\tNext Check >>> {Wait_minutes * wait_job_count} minute later! forecast: {formated_end_time}    ")
                sys.stdout.flush()
                time.sleep(Wait_minutes * wait_job_count * 60)
                RunningJobIDList = self.Running_JobIDList()
                Flag, wait_job_count = self.check_jobs(RunningJobIDList, MyJobIDList)
            else:
                print(f"{Color.GREEN}\n\n"
                      f"Calculation cycles for {which} until JobID {MyJobIDList[-1]} were finished.{Color.RESET}")
        return

    @staticmethod
    def My_JobIDList(qsubList):
        """
        Get the list of job IDs
        :param qsubList:
        :return:
        """
        jobCount = len(qsubList)
        result = subprocess.run(["qstat"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
        lines = result.stdout.strip().splitlines()
        running_job_count = len(lines) - 2
        if not len(lines) == 0:
            runningJobID = []
            lastJobID = int(lines[-1].strip().split()[0])
            for i in range(running_job_count):
                JobID = int(lines[i + 2].strip().split()[0])
                runningJobID.append(JobID)
            myStartID = int(lastJobID - jobCount + 1)
            myJobID = list(range(int(myStartID), lastJobID + 1))
            myJobIDs = []
            for A in myJobID:
                myJobIDs.append(A)
            myJobIDs.sort()
            return myJobIDs
        else:
            return [0]

    @staticmethod
    def Running_JobIDList():
        """
        Get the list of running job IDs
        :return:
        """
        result = subprocess.run(["qstat"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
        lines = result.stdout.strip().splitlines()
        running_job_count = len(lines) - 2
        if len(lines) == 0:
            runningJobID = []
        else:
            runningJobID = []
            for i in range(running_job_count):
                JobID = int(lines[i + 2].strip().split()[0])
                runningJobID.append(JobID)
            runningJobID.sort()
        return runningJobID

    @staticmethod
    def check_jobs(current_jobs, my_jobs):
        """
        Check the jobs
        :param current_jobs:
        :param my_jobs:
        :return:
        """
        contains_my_jobs = any(job in current_jobs for job in my_jobs)
        my_job_positions = [i for i, job in enumerate(current_jobs) if job in my_jobs]
        if my_job_positions:
            last_position = max(my_job_positions)
            remaining_time = last_position + 1
        else:
            remaining_time = 0
        if remaining_time == 2 or remaining_time == 1:
            remaining_time = 0.5
        else:
            pass
        return contains_my_jobs, remaining_time

    @staticmethod
    def getElapsedTime(start_time):
        """
        Get the elapsed time
        :param start_time:
        :return:
        """
        now = datetime.datetime.now()
        elapsed_time = now - start_time
        elapsed_time = round(elapsed_time.total_seconds() / 60, 1)
        formated_NOW = now.strftime("%m/%d %H:%M")
        return formated_NOW, elapsed_time

    @staticmethod
    def rmWildCards(wildcard):
        """
        Remove the wildcard
        :param wildcard:
        :return:
        """
        lines = glob.glob(wildcard)
        for line in lines:
            subprocess.run(["rm", line], timeout=10)
        return

    def readEnergy(self):
        OtherList = []
        LogList = glob.glob(f"{self.dirpath}/{self.MaterName}_3mol{self.mol_pos}_*.log")
        LogList.sort()
        for LogName in LogList:
            FileName, Edge, Faceon, Other = self.getVAL_fromLogName(LogName)
            OtherList.append(Other)
        OtherList = sorted(set(OtherList))

        with open(f"{self.MaterName}_3mol{self.mol_pos}_all.txt", "w") as AllData:
            with open(f"{self.MaterName}_3mol{self.mol_pos}_min.txt", "w") as MinData:
                header = ("Distance in Other direction (Å)\tDistance in Edge direction (Å)\t"
                          "Distance in Faceon direction (Å)\tCounterpoise corrected energy (A.U)"
                          "\tBSSE energy (A.U)")
                AllData.write(f"*****  {self.MaterName}_3mol{self.mol_pos} All Results *****\n")
                MinData.write(
                    f"***** {self.MaterName}_3mol{self.mol_pos} Minimum Energy at each Angle *****\n")
                AllData.write(f" \t{header}\n")
                MinData.write(f"{header}\n")
                MinimumL = []

                print(f"\t{header}")

                for Other in OtherList:
                    CPE_Dict = {}
                    VAL_Dict = {}
                    Keys = []
                    AllData.write("******\t******\t******\t******\t******\t******\n")
                    for LogName in LogList:
                        FileName, VEdge, VFaceon, VOther = self.getVAL_fromLogName(LogName)
                        if Other == VOther:
                            with open(LogName, "r") as f:
                                data = f.read()
                            if "Normal termination" in data:
                                CPE, BSE = self.getEnergy(data)
                                CPE_Dict[FileName] = CPE
                                VAL_Dict[FileName] = f"{VOther}\t{VEdge}\t{VFaceon}\t{CPE}\t{BSE}\n"
                                Keys.append(FileName)
                            else:
                                print(f"{Color.RED}Error: {LogName} did not finish normally.{Color.RESET}")
                                with open(f"./{self.MaterName}_3mol{self.mol_pos}_Error.log", "a") as f:
                                    f.write(f"{LogName}\n")
                                with open(f"./{self.MaterName}_3mol{self.mol_pos}_Error.log", "r") as f:
                                    error_lines = f.readlines()
                                if len(error_lines) != len(set(error_lines)):
                                    print(f"{Color.RED}Error: {LogName} is duplicated.{Color.RESET}")
                                    self.messages.append(f"{Color.RED}Error: {LogName} is duplicated.{Color.RESET}")
                                    self.HelpList.append(True)
                                self.help_check_exit()
                                self.rmWildCards(f"{LogName}")

                    minkey = min(CPE_Dict, key=CPE_Dict.get)
                    MinimumL.append(minkey)
                    MinData.write(f"{VAL_Dict[minkey]}")
                    Keys.sort()
                    for Key in Keys:
                        if Key == minkey:
                            sentence = f"*\t{VAL_Dict.get(Key)}"
                            AllData.write(sentence)
                            print(sentence.strip())

                        else:
                            sentence = f"-\t{VAL_Dict.get(Key)}"
                            AllData.write(sentence)
                            print(sentence.strip())
        return None

    @staticmethod
    def getVAL_fromLogName(LogName):
        FileName = LogName.split("/")[-1][:-4]
        Condition = LogName[:-4].split("/")[-1].split("_")
        Edge = int(Condition[3]) / 100
        Faceon = int(Condition[4]) / 100
        Other = int(Condition[2]) / 100
        return FileName, Edge, Faceon, Other

    @staticmethod
    def getEnergy(data):
        """
        Get the energy
        :param data:
        :return:
        """
        CPE = data.split("Counterpoise corrected energy =")[1]
        CPE = CPE.split("BSSE energy")[0]
        CPE = CPE.strip()
        CPE = float(CPE)
        BSE = data.split("BSSE energy =")[1]
        BSE = BSE.split("sum of fragments")[0]
        BSE = BSE.strip()
        BSE = float(BSE)
        return CPE, BSE

    @staticmethod
    def getRefLines(FileName):
        RefLines = []
        with open(FileName, "r") as f:
            lines = f.readlines()

        for line in lines:
            contents = line.strip().split()
            try:
                float(contents[0])
                RefLines.append(line.strip())
            except ValueError:
                pass
        return RefLines

    @staticmethod
    def getRefValues(Other, RefLines):
        RefValues = ''
        for RefLine in RefLines:
            Contents = RefLine.strip().split("\t")
            try:
                RefOther = round(float(Contents[0]), 2)
                if Other == RefOther:
                    RefValues = [round(float(Contents[1]), 2),
                                 round(float(Contents[2]), 2)]
                else:
                    pass
            except (ValueError, IndexError):
                pass
        return RefValues

    def mkNewConditionLists(self, which, dev):
        with open(f"./{self.MaterName}_3mol{self.mol_pos}_all.txt", "r") as f:
            AllLines = f.readlines()

        OtherList, Judges, Compel_Other, Compel_Val, NewConditions = [], [], [], [], []
        for AllLine in AllLines:
            Contents = AllLine.strip().split("\t")
            try:
                OtherList.append(round(float(Contents[1]), 2))
            except (ValueError, IndexError):
                pass
        OtherList = sorted(set(OtherList))
        RefLines = self.getRefLines(f"./{self.MaterName}_3mol{self.mol_pos}_min.txt")

        for Other in OtherList:
            SV, ValueList = 0, []
            RefValues = self.getRefValues(Other, RefLines)
            if which == "Edge":
                SV = float(RefValues[1])
            elif which == "Faceon":
                SV = float(RefValues[0])

            for AllLine in AllLines:
                Contents = AllLine.strip().split("\t")
                try:
                    DataOther = round(float(Contents[1]), 2)
                    if DataOther == Other:
                        DataDEdge, DataDFaceon = float(Contents[2]), float(Contents[3])
                        if which == "Edge":
                            RefDFaceon = round(float(RefValues[1]), 2)
                            if DataDFaceon == RefDFaceon:
                                ValueList.append(DataDEdge)
                            else:
                                pass
                            if Contents[0] == "*":
                                SV = DataDEdge
                        elif which == "Faceon":
                            RefDEdge = round(float(RefValues[0]), 2)
                            if DataDEdge == RefDEdge:
                                ValueList.append(DataDFaceon)
                            else:
                                pass
                            if Contents[0] == "*":
                                SV = DataDFaceon
                except (ValueError, IndexError):
                    pass
            if round(SV + 0.05, 2) in ValueList and round(SV - 0.05, 2) in ValueList:
                Judges.append("complete")
                print(f"\nThe local minimum by 0.05Å step for {Other} was Found;\t{SV}.")
                print(f"{Color.GREEN}\tCOMPLETE!!{Color.RESET}")
                Compel_Other.append(Other)
                Compel_Val.append(SV)
            elif round(SV + dev, 2) in ValueList and round(SV - dev, 2) in ValueList:
                Judges.append("complete")
                print(f"\nThe local minimum by {dev}Å step for {Other} was Found;\t{SV}.")
                print(f"{Color.GREEN}\tCOMPLETE!!{Color.RESET}")
                Compel_Other.append(Other)
                Compel_Val.append(SV)
            elif SV == min(ValueList):
                Judges.append("not complete")
                print(f"\nLocal minimum for {Other} was NOT Found in the cycle")
                print(f"{Color.RED}\tNOT COMPLETE(1)!!{Color.RESET}"
                      f"{Constant.Cn} new conditions bellow were appended.")
                for i in range(Constant.Cn):
                    NewConditions.append(self.mkNewCondition(Other, SV - (dev * (i + 1)), which, RefValues))
            elif SV == max(ValueList):
                Judges.append("not complete")
                print(f"\nLocal minimum for {Other} was NOT Found in the cycle")
                print(f"{Color.RED}\tNOT COMPLETE(2)!!{Color.RESET}"
                      f"{Constant.Cn} new conditions bellow were appended.")
                for i in range(Constant.Cn):
                    NewConditions.append(self.mkNewCondition(Other, SV + (dev * (i + 1)), which, RefValues))
        with open(f"./ConditionList_3mol{self.mol_pos}.txt", "r") as f:
            orgCondition = f.readlines()
        NewList = sorted(set(orgCondition + NewConditions))

        with open(f"./ConditionList_3mol{self.mol_pos}.txt", "w") as f:
            for NewCondition in NewList:
                f.write(f"{NewCondition}")
        if len(Compel_Other) == 0:
            print(f"{Color.RED}\nAll the local minimums were NOT found in the cycle.{Color.RESET}")
        else:
            print(f"\nLocal minimum was FOUND in {len(Compel_Other)} conditions.")
            for i in range(len(Compel_Other)):
                print(f"\t{Compel_Other[i]}: {Compel_Val[i]}")

        Judges = list(set(Judges))
        if "not complete" in Judges:
            judge = False
        elif "complete" in Judges and len(Judges) == 1:
            judge = True
        else:
            judge = False

        return judge

    def mkCycleConditions(self, which, dev):
        Others, NewConditions = [], []
        if dev == 0.2:
            n = Constant.CycleCondition_n_02
        elif dev == 0.1:
            n = Constant.CycleCondition_n_01
        elif dev == 0.05:
            n = Constant.CycleCondition_n_005
        else:
            n = 1

        RefLines = self.getRefLines(f"./{self.MaterName}_3mol{self.mol_pos}_min.txt")
        for RefLine in RefLines:
            Contents = RefLine.strip().split()
            try:
                Others.append(round(float(Contents[0]), 2))
            except (ValueError, IndexError):
                pass
        print(f"\n{Color.GREEN}Creating the new conditions...{Color.RESET}")
        print("\tNew conditions for the next cycle:")
        for Other in Others:
            ValueList = []
            RefValues = self.getRefValues(Other, RefLines)
            if which == "Edge":
                SV = float(RefValues[1])
            elif which == "Faceon":
                SV = float(RefValues[0])
            else:
                SV = int()
            with open(f"./{self.MaterName}_3mol{self.mol_pos}_all.txt", "r") as All:
                AllLines = All.readlines()
            for AllLine in AllLines:
                Contents = AllLine.strip().split("\t")
                try:
                    DataOther = round(float(Contents[1]), 2)
                    if DataOther == Other:
                        DataDEdge, DataDFaceon = float(Contents[2]), float(Contents[3])
                        if which == "Edge":
                            RefDFaceon = round(float(RefValues[1]), 2)
                            if DataDFaceon == RefDFaceon:
                                ValueList.append(DataDEdge)
                            else:
                                pass
                            if Contents[0] == "*":
                                SV = DataDEdge
                        elif which == "Faceon":
                            RefDEdge = round(float(RefValues[0]), 2)
                            if DataDEdge == RefDEdge:
                                ValueList.append(DataDFaceon)
                            else:
                                pass
                            if Contents[0] == "*":
                                SV = DataDFaceon
                except (ValueError, IndexError):
                    pass
            if round(SV + 0.05, 2) in ValueList and round(SV - 0.05, 2) in ValueList:
                pass
            else:
                Debug_message = [f"Other, SV: {Other}, {SV}"]
                self.debug_message(Debug_message)
                for i in range(n):
                    NewConditions.append(self.mkNewCondition(Other, SV - (dev * (i + 1)), which, RefValues))
                    NewConditions.append(self.mkNewCondition(Other, SV + (dev * (i + 1)), which, RefValues))
        with open(f"./ConditionList_3mol{self.mol_pos}.txt", "r") as f:
            orgCondition = f.readlines()
        NewList = sorted(set(orgCondition + NewConditions))
        with open(f"./ConditionList_3mol{self.mol_pos}.txt", "w") as f:
            for NewCondition in NewList:
                f.write(f"{NewCondition}")
        return None

    @staticmethod
    def CompareStructures(RefLines, LinesBefore):
        for RefLine in RefLines:
            RefContents = RefLine.strip().split()
            RefOther, RefDEdge, RefDFaceon = RefContents[0], RefContents[1], RefContents[2]
            for LineBefore in LinesBefore:
                Contents = LineBefore.strip().split()
                Other = Contents[0]
                if Other == RefOther:
                    if Contents[1] != RefDEdge or Contents[2] != RefDFaceon:
                        return False
        return True

    def getMinCondition(self):
        with open(f"./{self.MaterName}_3mol{self.mol_pos}_min.txt", "r") as f:
            lines = f.readlines()
        del lines[0:2]
        MinConditions = []
        for line in lines:
            contents = line.strip().split()
            Other = int(round(float(contents[0]) * 100))
            Edge = int(round(float(contents[1]) * 100))
            Faceon = int(round(float(contents[2]) * 100))
            MinConditions.append(f"{Other}_{Edge}_{Faceon}")
        return MinConditions

    def Calculate_TI(self, MinConditions):
        if self.calculation_tcal_flag:
            os.makedirs(self.tcalpath, exist_ok=True)
            qsubList = []
            print(f"\n"
                  f"Copying the com files for minimum energies...")
            for Condition in MinConditions:
                command = ["rename", "gjf", "com",
                           f"{self.dirpath}/{self.MaterName}_3mol{self.mol_pos}_{Condition}.gjf"]
                self.execute(command, False)
                command = ["cp", f"{self.dirpath}/{self.MaterName}_3mol{self.mol_pos}_{Condition}.com",
                           f"{self.tcalpath}/{self.MaterName}_3mol{self.mol_pos}_{Condition}.com"]
                self.execute(command, False)
                command = ["rename", "com", "gjf",
                           f"{self.dirpath}/{self.MaterName}_3mol{self.mol_pos}_{Condition}.com"]
                self.execute(command, False)
            print(f"\t>>> {Color.GREEN}Completed!!{Color.RESET}")
            self.mkXYZFiles()
            print(f"\n**********\n"
                  f"{Color.GREEN}Calculating transfer integrals...\n{Color.RESET}")
            if not os.path.exists(f"./{self.tcalpath}/{self.MaterName}_3mol{self.mol_pos}_tcal.log"):
                XYZs = glob.glob(f"{self.tcalpath}/*.xyz")
                for XYZ in XYZs:
                    if "_m1.xyz" in XYZ or "_m2.xyz" in XYZ or "-12.xyz" in XYZ or "-23.xyz" in XYZ or "-31.xyz" in XYZ:
                        XYZ = XYZ.replace("./", "")
                        os.remove(f"{XYZ}")
                    else:
                        pass
                self.XYZ_3mol_to_XYZ_2mol()
                with open(f"{self.tcalpath}/tcal.sh", "w") as f:
                    f.write(Stereotyped.tcal_sh_txt)
                qsubList.append("qsub tcal.sh")

                self.job_submission(qsubList, "tcal", self.tcalpath)

                if self.Debug:
                    pass
                else:
                    self.rmWildCards(f"{self.tcalpath}/*.sh*")
                subprocess.run(["rename", "tcal", f"{self.MaterName}_3mol{self.mol_pos}_tcal", "tcal.log"],
                               cwd=self.tcalpath)
                self.readlog()
            else:
                print(f"\t>>> tcal.log: {Color.GREEN}Already exists!!{Color.RESET}")
                print(f"\t>>> {Color.GREEN}Calculation of transfer integrals was skipped.{Color.RESET}")
        else:
            print(f"\n"
                  f"**********\n"
                  f"{Color.GREEN}Calculation of transfer integrals was skipped.{Color.RESET}")
        self.PhaseCheck()
        return

    @staticmethod
    def execute(command_list, TEXT, directory=None):
        command = ' '.join(command_list)
        if TEXT:
            print(f'> {command}')

        res = subprocess.run(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                             cwd=directory)
        # python3.6

        # res = subprocess.run(command_list, capture_output=True, text=True)
        # python3.7 or later

        if res.returncode:
            print(res.stderr.strip())
        return

    def mkXYZFiles(self):
        F_paths = glob.glob(f"{self.tcalpath}/*.com")
        print(f"\n"
              f"**********\n"
              f"Making XYZ files...")
        if not F_paths:
            print(f"\t>>> There is NO .com file in the {self.tcalpath} folder.")
            return

        for F_path in F_paths:
            self.process_file(F_path)

        ComFileNames = glob.glob(f"{self.tcalpath}/*.com")
        XYZFileNameList = []
        for ComFileName in ComFileNames:
            FormatedComfn = (ComFileName[:-4]).strip().split("/")[-1]
            XYZFileName = f"{FormatedComfn}.xyz"
            if self.Debug:
                print(f"\t{FormatedComfn}.com -> {XYZFileName}")
            XYZFileNameList.append(XYZFileName)
            subprocess.run(["newzmat", "-icart", "-oxyz", FormatedComfn, XYZFileName],
                           timeout=2, cwd=self.tcalpath)
        Ls = []
        for XYZFileName in XYZFileNameList:
            with open(f"{self.tcalpath}/{XYZFileName}", "r") as xyzf:
                lines = xyzf.readlines()
            Ls.append(len(lines))
            with open(f"{self.tcalpath}/{XYZFileName}", "w") as xyzf:
                xyzf.write(f"{len(lines)}\n")
                xyzf.write("00000001\n")
                for line in lines:
                    xyzf.write(line)
        print(f"\t>>> {self.tcalpath}/*.xyz: {Color.GREEN}Created!!{Color.RESET} ({len(Ls)} files)")
        return

    @staticmethod
    def process_file(F_path):
        """
        Process the file
        :param F_path:
        :return:
        """
        with open(F_path, "r") as f:
            lines = f.readlines()
        new_lines = [
            line.replace("X", "C") if line.strip().split() and line.strip().split()[0] == "X"
            else line
            for line in lines
        ]
        with open(F_path[:-4] + ".com", "w") as file:
            file.writelines(new_lines)
        return

    def XYZ_3mol_to_XYZ_2mol(self):
        filepaths = glob.glob(f"{self.tcalpath}/*_3mol*.xyz")

        Faults = []
        print("\nConverting 3mol to 2mol XYZ files...")
        for filepath in filepaths:
            Mol1, Mol2, Mol3 = [], [], []
            file_core = filepath[:-4]
            with open(filepath, "r") as f:
                # Read the first line (the number of atoms)
                number_atoms = f.readline()
                Comment = f.readline()
                coordinates = f.readlines()
            self.messages.append("\t" + filepath.strip())

            if int(float(number_atoms)) % 3 != 0:
                Faults.append(filepath)
                self.messages.append(f"\t>>> {Color.RED}Error: {filepath} did not divide by 3.{Color.RESET}")
            else:
                Atoms_inMol = int(round(float(number_atoms) / 3, 5))
                for i in range(Atoms_inMol):
                    Mol1.append(coordinates[i])
                    Mol2.append(coordinates[i + Atoms_inMol])
                    Mol3.append(coordinates[i + 2 * Atoms_inMol])
                Atoms_in_TcalXYZFile = int(round(Atoms_inMol * 2, 3))

                with open(f"{file_core}-12.xyz", "w") as newF12, open(f"{file_core}-23.xyz", "w") as newF23, open(
                        f"{file_core}-31.xyz", "w") as newF31:
                    newF12.write(f"{Atoms_in_TcalXYZFile}\n"
                                 f"{Comment}")
                    newF23.write(f"{Atoms_in_TcalXYZFile}\n"
                                 f"{Comment}")
                    newF31.write(f"{Atoms_in_TcalXYZFile}\n"
                                 f"{Comment}")
                    for i in range(Atoms_inMol):
                        newF12.write(Mol1[i])
                        newF23.write(Mol2[i])
                        newF31.write(Mol3[i])
                    for i in range(Atoms_inMol):
                        newF12.write(Mol2[i])
                        newF23.write(Mol3[i])
                        newF31.write(Mol1[i])
                self.messages.append(f"\t>>> {file_core}-12.xyz: Created!!")
                self.messages.append(f"\t>>> {file_core}-23.xyz: Created!!")
                self.messages.append(f"\t>>> {file_core}-31.xyz: Created!!")
                os.rename(filepath, f"{file_core}.all")
        if len(Faults) != 0:
            for Fault in Faults:
                self.messages.append(f"\t>>> {Color.RED}Error: {Fault} did not divide by 3.{Color.RESET}")
            self.HelpList.append(True)
        else:
            self.messages.append(f"\n\t>>> {Color.GREEN}XYZ 3mol to 2mol: Succeeded!!{Color.RESET}")
        self.help_check_exit()
        return

    def readlog(self):
        print(f"\nReading the Tcal log file...")
        filepath = f"{self.tcalpath}/{self.MaterName}_3mol{self.mol_pos}_tcal.log"
        output_filepath = f"{self.tcalpath}/{self.MaterName}_3mol{self.mol_pos}_TIs.txt"
        with open(filepath, "r") as f:
            lines = f.readlines()

        DataList = []
        keywords = ["Input File Name:", "NLUMO", "LUMO", "HOMO", "NHOMO"]
        for line in lines:
            for keyword in keywords:
                if keyword in line:
                    DataList.append(line)
                    break
        if len(DataList) % 5 != 0:
            self.messages.append(f"\t>>> {Color.RED}Error: UNEXPECTED ERROR in read Tcal log process.{Color.RESET}")
            self.HelpList.append(True)
            self.help_check_exit()

        with open(output_filepath, "w") as file:
            header = (f"***** Transfer Integrals in "
                      f"{self.tcalpath}/{self.MaterName}_3mol{self.mol_pos}_tcal.log *****\n")
            columns = "input file\tNLUMO (meV)\tLUMO (meV)\tHOMO (meV)\tNHOMO (meV)\n"
            file.write(header + columns)
            print(header.strip())
            print(columns.strip())
            for i in range(0, len(DataList), 5):
                input_file = DataList[i].split(": ")[1].split(".xyz")[0]
                NLUMO = self.extract_value(DataList[i + 1], "NLUMO")
                LUMO = self.extract_value(DataList[i + 2], "LUMO")
                HOMO = self.extract_value(DataList[i + 3], "HOMO")
                NHOMO = self.extract_value(DataList[i + 4], "NHOMO")
                file.write(f"{input_file}\t{NLUMO}\t{LUMO}\t{HOMO}\t{NHOMO}\n")
                print(f"{input_file}\t{NLUMO}\t{LUMO}\t{HOMO}\t{NHOMO}")
        return

    @staticmethod
    def extract_value(line, keyword):
        return float(line.split(keyword)[-1].split()[0])

    def PhaseCheck(self):
        print("\n**********\nPhase Checking...")
        chkKEY = "-12"
        CalCores = self.mkCalCoreList(chkKEY)
        if len(CalCores) == 0:
            print(f"\t>>> Any chk file is NOT found in the specified folder.\n"
                  f"\t>>> The process for phase check was {Color.GREEN}skipped.{Color.RESET}")
        else:
            with open(f"{self.tcalpath}/PhaseCheck.txt", "w") as file:
                file.write("Name\tfor LUMO\tfor HOMO\n")
                for CalCore in CalCores:
                    self.mkCubeFile(f"{CalCore}_m1.chk")
                    self.mkCubeFile(f"{CalCore}_m2.chk")

                    HomoChk = self.ComparePhase(CalCore, "HOMO")
                    LumoChk = self.ComparePhase(CalCore, "LUMO")

                    file.write(f"{CalCore}\t{LumoChk}\t{HomoChk}\n")
            if self.Debug:
                pass
            else:
                self.rmWildCards(f"{self.tcalpath}/*.chk")
                self.rmWildCards(f"{self.tcalpath}/*.log")
                self.rmWildCards(f"{self.tcalpath}/*.gjf")
        return

    def mkCalCoreList(self, chkKEY):
        FileList = glob.glob(f"{self.tcalpath}/*{chkKEY}*")
        CalCores = []
        for file in FileList:
            if "_m1.chk" in file or "_m2.chk" in file:
                CalCores.append(file[file.rfind("/") + 1:-7])
            else:
                pass
        CalCores = list(set(CalCores))
        return CalCores

    def mkCubeFile(self, chkfile):
        def run_command(command, description):
            subprocess.run(command, cwd=self.tcalpath, timeout=10)
            print(f"\t{description}: {Color.GREEN}Completed!!{Color.RESET}")

        print(f"For {chkfile} ...")
        fchk = f"{chkfile[:-4]}.fch"
        HomoCub = f"{chkfile[:-4]}_HOMO.cub"
        LumoCub = f"{chkfile[:-4]}_LUMO.cub"

        run_command(["formchk", chkfile, fchk], f"File Conversion ({chkfile} -> {fchk})")
        run_command(["cubegen", "0", "MO=Homo", fchk, HomoCub, "-2", "h"], f"{HomoCub} has been created")
        run_command(["cubegen", "0", "MO=Lumo", fchk, LumoCub, "-2", "h"], f"{LumoCub} has been created")
        run_command(["rm", fchk], f"Removing {fchk}")
        return

    def ComparePhase(self, CalCore, HomoLumo):
        """
        :param CalCore:
        :param HomoLumo:
        :return:
        """
        Val1 = self.ReadCube(f"{CalCore}_m1_{HomoLumo}.cub")
        Val2 = self.ReadCube(f"{CalCore}_m2_{HomoLumo}.cub")

        if len(Val1) != len(Val2):
            return "not available(1) (Different number of grids in the file)"

        PClist = []
        print(f"{CalCore} {HomoLumo} Data Excerpt")
        for i in random.sample(range(len(Val1)), min(20, len(Val1))):
            if abs(Val1[i]) > 1e-7:
                Val = Val1[i] * Val2[i]
                print(f"m1 value = {Val1[i]}\tm2 value = {Val2[i]}")
                PClist.append("Opposit" if Val < 0 else "Same")

        PC = list(set(PClist))
        if len(PC) == 1:
            return PC[0]
        elif len(PC) > 1:
            return "not available(2)"
        else:
            return "not available(3)"

    def ReadCube(self, Cubefile):
        """
        Read the cube file
        :param Cubefile:
        :return:
        """
        with open(f"{self.tcalpath}/{Cubefile}", "r") as file:
            for _ in range(2):  # skip the first two lines (comment and number of atoms)
                next(file)

            NumE, Xgrid, Ygrid, Zgrid = (int(next(file).split()[0]) for _ in range(4))
            for _ in range(NumE):  # skip the number of electrons
                next(file)

            values = [float(x) for _ in range(Xgrid * Ygrid)
                      for line in file for x in line.split()]

        subprocess.run(["rm", Cubefile], cwd=self.tcalpath, timeout=10)
        return values

    def Result_Data_set(self, before):
        result_name = f"{self.MaterName}_3mol{self.mol_pos}"
        result_path = f"{result_name}_result"
        MinFileName = f"{result_name}_min.txt"
        TIFileName = f"{result_name}_TIs.txt"
        AllFileName = f"{result_name}_all.txt"
        PCFileName = "PhaseCheck.txt"
        Lacks = []
        print(f"\n**********\n{Color.GREEN}Resulting Data Set for {result_name}...{Color.RESET}")
        os.makedirs(result_path, exist_ok=True)

        if os.path.isfile(MinFileName) and os.path.isfile(f"{self.tcalpath}/{TIFileName}"):
            self.combineData(TIFileName, MinFileName, PCFileName, result_path)
        else:
            self.copy_file(MinFileName, f"{result_path}/{MinFileName}", Lacks)
            print("\n")

        print("\tCopying Files required for recalculation...")
        self.copy_file(AllFileName, f"{result_path}/{AllFileName}", Lacks)

        # MaterName.xyz
        self.copy_file(f"{self.MaterName}.xyz", f"{result_path}/{self.MaterName}.xyz", Lacks)

        # ConditionList_3mol{self.mol_pos}.txt
        self.copy_file(f"ConditionList_3mol{self.mol_pos}.txt",
                       f"{result_path}/ConditionList_3mol{self.mol_pos}.txt", Lacks)

        # CalcSetting_BW.txt
        self.copy_file("CalcSetting_BW.txt", f"{result_path}/CalcSetting_BW.txt", Lacks)

        # min file
        self.copy_file(f"{result_name}_min.txt",
                       f"{result_path}/{result_name}_min.txt", Lacks)
        self.copy_file(f"{result_name}_mins.hist",
                       f"{result_path}/{result_name}_mins.hist", Lacks)

        # structure files
        structure_files = glob.glob(f"{self.tcalpath}/*.all")
        if not structure_files:
            print(f"\t>>> There is NO .all file in the {self.tcalpath} folder.")
            Lacks.append("Structure files")
        else:
            print(f"\n\tCopying the structure files...")
            for file in structure_files:
                name = os.path.basename(file).replace(".all", "")
                self.copy_file(file, f"{result_path}/{name}.xyz", Lacks)
            print(f"\t>>> {Color.GREEN}Structural Files were copied into the {result_path} folder.{Color.RESET}")

        if Lacks:
            self.messages.append(f"\nSeveral result files were not found in the specific folder:")
            for lack in Lacks:
                self.messages.append(f"\t>>> {lack}")
            self.HelpList.append(True)
        else:
            self.messages.append(f"\n\t{Color.GREEN}All files were copied successfully!!{Color.RESET}")
        self.help_check_exit()

        after = time.time()
        print(f"\n"
              f"Elapsed Time: {(after - before):.0f} s\n{Color.GREEN}")
        return None

    def copy_file(self, src, dest, lacks_list):
        if os.path.isfile(src):
            self.execute(["cp", src, dest], False)
            print(f"\t>>> {src} -> {dest}: {Color.GREEN}Complete{Color.RESET}")
        else:
            lacks_list.append(src)
        return

    def combineData(self, TIFileName, MinFileName, PCFileName, result_path):
        print(f"\tCombining Data...\n")
        with open(f"{self.tcalpath}/{TIFileName}", "r") as f:
            TI_lines = f.readlines()
        del TI_lines[0:2]

        with open(f"{MinFileName}", "r") as f:
            Min_lines = f.readlines()
        del Min_lines[0:2]

        if os.path.exists(f"{self.tcalpath}/{PCFileName}"):
            with open(f"{self.tcalpath}/{PCFileName}", "r") as f:
                PC_lines = f.readlines()
            del PC_lines[0:2]
        else:
            PC_lines = []

        if len(TI_lines) == len(Min_lines) * 3:
            pass
        else:
            self.messages.append(f"\t>>>{Color.RED}Error: The numbers of data lines in {MinFileName} "
                                 f"and {TIFileName} DO NOT match.{Color.RESET}\n"
                                 f"A file was NOT be changed.")
            self.HelpList.append(True)
        self.help_check_exit()

        CombLines = []
        CombLines12 = []
        CombLines23 = []
        CombLines31 = []
        print("Entry\tAngle\tDcol\tDtrv\tCpCE\tBSE\t**"
              "\tTI-NLUMO\tTI-LUMO\tTI-HOMO\tTI-NHOMO\t**\tPC (LUMO)\tTI-LUMO\tPC (HOMO)\tTI-HOMO")

        for Min_line in Min_lines:
            MinData = Min_line.split()

            Entry = (f"{self.MaterName}_3mol{self.mol_pos}"
                     f"_{int(round(float(MinData[0]) * 100, 5))}"
                     f"_{int(round(float(MinData[1]) * 100, 5))}"
                     f"_{int(round(float(MinData[2]) * 100, 5))}")

            for TI_line in TI_lines:
                TIData = TI_line.split()
                if (TIData[0] == Entry
                        or TIData[0] == f"{Entry}-12"
                        or TIData[0] == f"{Entry}-23"
                        or TIData[0] == f"{Entry}-31"):
                    CombLine_temp = (
                        f"{TIData[0]}\t{MinData[0]}\t{MinData[1]}\t{MinData[2]}\t{MinData[3]}\t{MinData[4]}"
                        f"\t**\t{TIData[1]}\t{TIData[2]}\t{TIData[3]}\t{TIData[4]}")

                    HomoChk = "yet"
                    LumoChk = "yet"
                    if len(PC_lines) != 0:
                        for PCLine in PC_lines:
                            PCData = PCLine.split()
                            if PCData[0].strip() == TIData[0].strip():
                                LumoChk = PCData[1]
                                HomoChk = PCData[2]
                    else:
                        pass

                    correctTI_LUMO = self.correctTI(LumoChk, TIData[2])
                    correctTI_HOMO = self.correctTI(HomoChk, TIData[3])
                    CombLine = f"{CombLine_temp}\t**\t{LumoChk}\t{correctTI_LUMO}\t{HomoChk}\t{correctTI_HOMO}\n"
                    if TIData[0] == Entry:
                        CombLines.append(CombLine)
                    elif TIData[0] == f"{Entry}-12":
                        CombLines12.append(CombLine)
                    elif TIData[0] == f"{Entry}-23":
                        CombLines23.append(CombLine)
                    elif TIData[0] == f"{Entry}-31":
                        CombLines31.append(CombLine)
                    else:
                        pass
                else:
                    pass
        self.help_check_exit()
        CombLines.sort()
        CombLines12.sort()
        CombLines23.sort()
        CombLines31.sort()
        self.saveCombData(f"{self.MaterName}_3mol{self.mol_pos}",
                          f"{result_path}/{self.MaterName}_3mol{self.mol_pos}_min-TIs.txt", CombLines)
        self.saveCombData(f"{self.MaterName}_3mol{self.mol_pos}-12",
                          f"{result_path}/{self.MaterName}_3mol{self.mol_pos}_min-TIs-12.txt", CombLines12)
        self.saveCombData(f"{self.MaterName}_3mol{self.mol_pos}-23",
                          f"{result_path}/{self.MaterName}_3mol{self.mol_pos}_min-TIs-23.txt", CombLines23)
        self.saveCombData(f"{self.MaterName}_3mol{self.mol_pos}-31",
                          f"{result_path}/{self.MaterName}_3mol{self.mol_pos}_min-TIs-31.txt", CombLines31)
        print(f"\n"
              f"\t>>> {Color.GREEN}Combining Data: Succeeded!!{Color.RESET}\n"
              f"\t>>> {MinFileName} and {TIFileName} were combined into "
              f"{result_path}/{self.MaterName}_3mol{self.mol_pos}_min-TIs*.txt.\n")
        return

    @staticmethod
    def correctTI(PhaseChk, TI):
        TI = float(TI)
        if PhaseChk == "Same":
            CorrectTI = TI
        elif PhaseChk == "Opposit":
            CorrectTI = -1 * TI
        else:
            CorrectTI = TI
        return CorrectTI

    @staticmethod
    def saveCombData(Name, path, List):
        header = (
            f"*************** {Name} Minimum Energy and Transfer Integrals ***************\n"
            "Entry\tAngle (deg.)\tD in col. (Angstrom)\tD in transv. (Angstrom)\t"
            "Counterpoise corrected energy (AU)\tBSSE energy (AU)\t**\t"
            "TI-NLUMO (meV)\tTI-LUMO (meV)\tTI-HOMO (meV)\tTI-NHOMO (meV)\t**\t"
            "PhaseChk (LUMO)\tTI (LUMO)\tPhaseChk (HOMO)\tTI (HOMO)\n"
        )
        if len(List) == 0:
            pass
        else:
            with open(path, "w") as File:
                File.write(header)
                print(f"\n{header}")
                for line in List:
                    File.write(line)
                    print(line.strip())
        return


if __name__ == "__main__":
    main()
