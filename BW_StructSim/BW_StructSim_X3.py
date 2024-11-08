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

    # Create the first condition
    bw.mkFirstCondition()

    # Search for the most stable structure
    bw.Most_Stable_Search()


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
    initial_offset = 3.0


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
        self.MaterName = args.MaterName[:-4]
        self.Flag_xyz = True
        self.Debug = args.debug

        with open(f"{self.MaterName}.xyz", "r") as f:
            self.NinMol = f.readline()
            f.readline()
            self.AtomList = f.readlines()

        # Select the structure of the molecule
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
        self.tcalpath = f"./{self.MaterName}_3mol{self.mol_pos}_tcal"
        os.makedirs(self.dirpath, exist_ok=True)

        # Retrieve the operator name
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
        params = {"edge axis": "", "faceon axis": "", "mol3 other_transition": "", "other axis": ""}
        for line in lines[:4]:
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
        self.rotate = rotate
        self.Mol3_Other = self.mkDirection(
            params["mol3 other_transition"], params["other axis"], "Mol3 Other Transition")
        Debug_message_List = [f"Edge Axis: {self.Edge_Axis}",
                              f"Faceon Axis: {self.Faceon_Axis}",
                              f"Other Axis: {self.Other_Axis}",
                              f"Rotate: {rotate}",
                              f"Mol3 Other Transition: {self.Mol3_Other}"]
        self.debug_message(Debug_message_List)

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
            self.messages.append(f"\t>>> InitialCondition_3mol{self.mol_pos}.txt: {Color.GREEN}Created.{Color.RESET}\n"
                                 f"\t>>> Please set the calculation conditions in the file.")
        self.help_check_exit()
        return None

    def mkFirstCondition(self):
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
        First_Edge = self.transform_number((Edge_Max - Edge_Min) / 2) + Edge_Max + Constant.initial_offset
        First_Faceon = self.transform_number((Faceon_Max - Faceon_Min) / 2) + Faceon_Max + Constant.initial_offset
        Debug_Message_List = [
            f"First Edge: {First_Edge}",
            f"First Faceon: {First_Faceon}"
        ]
        self.debug_message(Debug_Message_List)

        # Create the first condition
        self.messages.append(f"\n{Color.GREEN}Creating the first condition...{Color.RESET}")
        if os.path.exists(f"./ConditionList_3mol{self.mol_pos}.txt"):
            self.messages.append(f"\t>>> ConditionList_3mol{self.mol_pos}.txt: {Color.GREEN}Exist.{Color.RESET}\n"
                                 f"\t>>> Performs calculations based on conditions that exist.")
        else:
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
                    f.write(f"{Condition}\n")
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
            Edge = round(Val, 2)
            Faceon = RefValues[1]
        elif which == "Faceon":
            Edge = RefValues[0]
            Faceon = round(Val, 2)
        else:
            Edge = RefValues[0]
            Faceon = RefValues[1]
        NewCondition = f"{int(Other * 100)}_{int(round(Edge * 100))}_{int(round(Faceon * 100))}"
        print(f"\t\t{NewCondition.strip()}")
        return NewCondition

    def Most_Stable_Search(self):
        self.getTemporaryStructure()

    def getTemporaryStructure(self):
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
                    qsubList.append(self.mkFiles(Condition))
            self.job_submission(qsubList, "DEdge")
            if self.Debug:
                pass
            else:
                self.rmWildCards(f"{self.dirpath}/*.sh*")
                self.rmWildCards(f"{self.dirpath}/*.chk")
            print("\n**********\nReading Data...\n")
            self.readEnergy()
            exit()
        pass

    @staticmethod
    def getConditions(FileName):
        Condition = []
        with open(FileName, "r") as f:
            lines = f.readlines()
        for line in lines:
            Condition.append(line.strip())
        return Condition

    def mkFiles(self, Condition):
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

        self.write_gjf_file(f"{self.dirpath}/{FileName}.gjf",
                            Headers, Element, Mol1_pos, Mol2_pos, Mol3_pos)
        print(f"{self.dirpath}/{FileName}.gjf have been created.")
        if self.Flag_xyz:
            self.write_xyz_file(f"{self.dirpath}/{FileName}.xyz",
                                Element, Mol1_pos, Mol2_pos, Mol3_pos)
            print(f"{self.dirpath}/{FileName}.xyz have been created.")

        with open(f"{self.dirpath}/G.sh", "r") as f:
            lines = f.readlines()
        lines[12] = f"g16 {GJF_FileName}\n"
        with open(f"{self.dirpath}/{SH_FileName}", "w") as f:
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

    # job submission
    def job_submission(self, qsubList, which):
        """
        Submit the job
        :param qsubList:
        :param which:
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
                subprocess.run(qsub, cwd=self.dirpath)

            Wait_minutes = 2
            MyJobIDList = self.My_JobIDList(qsubList)
            MyJobIDList.sort()
            RunningJobIDList = self.Running_JobIDList()

            Flag, wait_job_count = self.check_jobs(RunningJobIDList, MyJobIDList)
            start_time = datetime.datetime.now()

            if Flag:
                formated_ST = start_time.strftime("%m/%d %H:%M:%S")
                if which == "DEdge":
                    term = "the stable distance in Edge direction"
                elif which == "DFaceon":
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

    def mkNewConditionLists(self):
        with open(f"./{self.MaterName}_3mol{self.mol_pos}_all.txt", "r") as f:
            Alllines = f.readlines()


if __name__ == "__main__":
    main()
