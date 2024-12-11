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
    CalcSetting_HB = ("Column Direction: [x, y, or z]\n"
                      "Transverse Direction: [x, y, or z]\n"
                      "Tilt Axis: [x, y, or z]\n"
                      "Tilt Angle: [Number]\n"
                      "\n"
                      "[Comment]\n")
    InitialCondition_2mol_Temp = ("10 7.0\n"
                                  "20 6.3\n"
                                  "30 5.7\n"
                                  "40 5.0\n"
                                  "50 4.5\n"
                                  "60 4.2\n"
                                  "70 3.9\n"
                                  "80 3.9\n"
                                  "90 3.9")
    InitialCondition_3mol_Temp = ("10 3.9\n"
                                  "20 3.9\n"
                                  "30 4.2\n"
                                  "40 4.5\n"
                                  "50 5.0\n"
                                  "60 5.7\n"
                                  "70 6.3\n"
                                  "80 7.0\n"
                                  "90 7.5")


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
    args, MaterName, Debug, Nmol, calculation_tcal_Flag, before = arg_parser(messages, HelpList)

    # Get the tilt angle
    Tilt, Formated_Tilt = get_TiltAngle(messages, HelpList, MaterName)

    # Acquisition of structure
    mol_pos = structure_acquisition(Nmol)

    # Check if the required files exist
    File_Set_Check(messages, HelpList, MaterName, Nmol, mol_pos, Formated_Tilt)

    # Get the operator
    Operator = getOperator()

    # Make the condition file
    which, dev, dirpath, tcal_path, RefLines = mkConditionFile(MaterName, Nmol, mol_pos,
                                                               Formated_Tilt, messages, HelpList)

    # Search for the most stable structure
    Most_Stable_Structure_Search(MaterName, Nmol, mol_pos, Tilt, Formated_Tilt,
                                 dirpath, which, dev, RefLines, Operator, messages, HelpList, Debug, tcal_path, args)

    # Calculate the transfer integral
    Calculate_TI(calculation_tcal_Flag, tcal_path, MaterName, Nmol, mol_pos,
                 Formated_Tilt, Debug, messages, HelpList)

    # Save the results
    Result_Data_set(MaterName, Nmol, Formated_Tilt, mol_pos, tcal_path, messages, HelpList)

    after = time.time()
    # End of the program
    print(f"\n"
          f"Elapsed Time: {(after - before):.0f} s\n{Color.GREEN}"
          f"************************* ALL PROCESSES END *************************"
          f"{Color.RESET}\n")
    return


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
    before = time.time()

    Help_Text = ("Structural optimisation and transfer integration of HerringBone structures.\n"
                 "Pass [molecule name].xyz as the first argument.\n"
                 "Pass the number of mol as the first argument.\n"
                 "If required, you will be asked to enter the Tilt angle, a name that uniquely identifies you, "
                 "and a structural pattern.")

    parser = argparse.ArgumentParser(description=Help_Text, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('MaterNameXYZ',
                        help="MaterNameXYZ.xyz")
    parser.add_argument('--debug', '-d', '-D',
                        help="Start the programme in Debug mode.",
                        action="store_true")

    parser.add_argument('--notcal', '--nt',
                        help="Argument for not calculating tcal.",
                        action="store_true")
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
    parser.add_argument('--energy', '-e',
                        action="store_true")

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

    calculation_tcal_Flag = args.notcal
    if args.notcal:
        messages.append(f"{Color.RED}Runs the program without transfer integral calculations.{Color.RESET}")
    else:
        pass

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

    # get the tilt angle
    Tilt, Formated_Tilt = get_TiltAngle(messages, HelpList, MaterName)
    file_path = f"./StructCheck-{MaterName}-t{Formated_Tilt}d"
    os.makedirs(file_path, exist_ok=True)
    Flag_XYZ = args.xyz
    Temp_SHs = []

    with open(f"{file_path}/G.sh", "w") as f:
        f.write(Stereotyped.Sh_txt)
    Temp_SHs.append(f"G.sh")

    if not args.manual:
        print(f"Automatically generate files for structural verification...")

        File_Name = f"{MaterName}_2mol_t{Formated_Tilt}d_90d-400"
        Temp_SH = mkFiles(MaterName, "2mol", "", "90d-400", "Temp",
                          file_path, Tilt, Formated_Tilt, Flag_XYZ, messages, HelpList)
        Temp_SHs.append(Temp_SH.replace("qsub", "").strip())
        messages.append(f"\t>>> {File_Name}.gjf: {Color.GREEN}Complete.{Color.RESET}")

        File_Name = f"{MaterName}_2mol_t{Formated_Tilt}d_60d-420"
        Temp_SH = mkFiles(MaterName, "2mol", "", "60d-420", "Temp",
                          file_path, Tilt, Formated_Tilt, Flag_XYZ, messages, HelpList)
        Temp_SHs.append(Temp_SH.replace("qsub", "").strip())
        messages.append(f"\t>>> {File_Name}.gjf: {Color.GREEN}Complete.{Color.RESET}")

        File_Name = f"{MaterName}_2mol_t{Formated_Tilt}d_30d-600"
        Temp_SH = mkFiles(MaterName, "2mol", "", "30d-600", "Temp",
                          file_path, Tilt, Formated_Tilt, Flag_XYZ, messages, HelpList)
        Temp_SHs.append(Temp_SH.replace("qsub", "").strip())
        messages.append(f"\t>>> {File_Name}.gjf: {Color.GREEN}Complete.{Color.RESET}")

        File_Name = f"{MaterName}_3molp1_t{Formated_Tilt}d_60d-420-600"
        Temp_SH = mkFiles(MaterName, "3mol", "p1", "60d-420-600", "Temp",
                          file_path, Tilt, Formated_Tilt, Flag_XYZ, messages, HelpList)
        Temp_SHs.append(Temp_SH.replace("qsub", "").strip())
        messages.append(f"\t>>> {File_Name}.gjf: {Color.GREEN}Complete.{Color.RESET}")

        File_Name = f"{MaterName}_3molp2_t{Formated_Tilt}d_60d-420-600"
        Temp_SH = mkFiles(MaterName, "3mol", "p2", "60d-420-600", "Temp",
                          file_path, Tilt, Formated_Tilt, Flag_XYZ, messages, HelpList)
        Temp_SHs.append(Temp_SH.replace("qsub", "").strip())
        messages.append(f"\t>>> {File_Name}.gjf: {Color.GREEN}Complete.{Color.RESET}")

        File_Name = f"{MaterName}_3molp3_t{Formated_Tilt}d_60d-420-600"
        Temp_SH = mkFiles(MaterName, "3mol", "p3", "60d-420-600", "Temp",
                          file_path, Tilt, Formated_Tilt, Flag_XYZ, messages, HelpList)
        Temp_SHs.append(Temp_SH.replace("qsub", "").strip())
        messages.append(f"\t>>> {File_Name}.gjf: {Color.GREEN}Complete.{Color.RESET}")

        help_check_exit(messages, HelpList)
    else:
        print(f"{Color.RED}Manual mode selected...{Color.RESET}")
        while True:
            structure = input(f"{Color.GREEN}Please enter the structure you want to create.{Color.RESET}\n"
                              f"\t1: 2mol\n"
                              f"\t2: 3molp1 (H-H)\n"
                              f"\t3: 3molp2 (T-T)\n"
                              f"\t4: 3molp3(H-T)\n"
                              f"\t5: exit\n"
                              f"\t{Color.GREEN}>>> {Color.RESET}")

            if structure == "5":
                EXIT = input("\n"
                             "Do you want to exit the programme?\n"
                             f"\ty/n {Color.GREEN}>>> {Color.RESET}")
                if EXIT == "y":
                    print(f"\t>>> {Color.RED}Exit the programme.{Color.RESET}")
                    break
                else:
                    continue

            while True:
                Rotate_Angle = input(f"\n"
                                     f"Please Enter the Rotate Angle.\n"
                                     f"\t{Color.GREEN}>>> {Color.RESET}")
                try:
                    Rotate_Angle = int(Rotate_Angle)
                    break
                except (IndexError, ValueError):
                    print(f"{Color.RED}\tError: Enter a number.{Color.RESET}")
                    continue

            while True:
                D_Col = input(f"\n"
                              f"Please Enter the Column Direction.{Color.GREEN}(Å){Color.RESET}\n"
                              f"\t{Color.GREEN}>>> {Color.RESET}")
                try:
                    D_Col = float(D_Col)
                    break
                except (IndexError, ValueError):
                    print(f"{Color.RED}\tError: Enter a number.{Color.RESET}")
                    continue

            if structure == "1":
                messages.append(f"\n{Color.GREEN}Create a file...{Color.RESET}")
                Temp_Condition = f"{Rotate_Angle}d-{int(D_Col * 100)}"
                Temp_SH = mkFiles(MaterName, "2mol", "", Temp_Condition, "Temp",
                                  file_path, Tilt, Formated_Tilt, Flag_XYZ, messages, HelpList)
                File_Name = f"{MaterName}_2mol_t{Formated_Tilt}d_{Temp_Condition}"
                Temp_SHs.append(Temp_SH.replace("qsub", "").strip())
                messages.append(f"\t>>> {File_Name}.gjf: {Color.GREEN}Complete.{Color.RESET}")
                messages.append(f"\n"
                                f"\t>>> {Color.GREEN}Complete{Color.RESET}\n")

                help_check_exit(messages, HelpList)
                continue
            else:
                while True:
                    D_Transv = input(f"\n"
                                     f"Please Enter the Transverse Direction.{Color.GREEN}(Å){Color.RESET}\n"
                                     f"\t{Color.GREEN}>>> {Color.RESET}")
                    try:
                        D_Transv = float(D_Transv)
                        break
                    except (IndexError, ValueError):
                        print(f"{Color.RED}\tError: Enter a number.{Color.RESET}")
                        continue
                mol_pos = ""
                if structure == "2":
                    mol_pos = "p1"
                elif structure == "3":
                    mol_pos = "p2"
                elif structure == "4":
                    mol_pos = "p3"
                Temp_Condition = f"{Rotate_Angle}d-{int(D_Col * 100)}-{int(D_Transv * 100)}"
                messages.append(f"\n{Color.GREEN}Create a file...{Color.RESET}")
                Temp_SH = mkFiles(MaterName, "3mol", mol_pos, Temp_Condition, "Temp",
                                  file_path, Tilt, Formated_Tilt, Flag_XYZ, messages, HelpList)
                File_Name = f"{MaterName}_2mol_t{Formated_Tilt}d_{Temp_Condition}"
                Temp_SHs.append(Temp_SH.replace("qsub", "").strip())
                messages.append(f"\t>>> {File_Name}.gjf: {Color.GREEN}Complete.{Color.RESET}")
                messages.append(f"\n"
                                f"\t>>> {Color.GREEN}Complete{Color.RESET}\n")

                help_check_exit(messages, HelpList)
                continue

    print("\nDelete unnecessary files...")
    Temp_SHs = list(set(Temp_SHs))
    for SH in Temp_SHs:
        try:
            subprocess.run(["rm", SH], timeout=10, cwd=file_path, check=True, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            messages.append(f"\t>>> {Color.RED}{SH}: Failed to delete.{Color.RESET}")
            HelpList.append(True)
        else:
            messages.append(f"\t>>> {SH}: {Color.GREEN}Deleted.{Color.RESET}")
            HelpList.append(False)
    help_check_exit(messages, HelpList)
    if Debug:
        print("*************** Debug Finished!!! ***************")

    after = time.time()

    # End of the program
    print(f"\n"
          f"Elapsed Time: {(after - before):.0f} s\n{Color.GREEN}"
          f"************************* ALL PROCESSES END *************************"
          f"{Color.RESET}\n")
    exit()


# Get the tilt angle
def get_TiltAngle(messages, HelpList, MoleculeName):
    """
    function to get the tilt angle
    :param messages: list of messages to be displayed
    :type messages: list
    :param HelpList: list of help flags
    :type HelpList: list
    :param MoleculeName: name of the molecule
    :type MoleculeName: str
    :rtype: float
    :rtype: int
    :return:
    """
    print("Checking if 'CalcSetting_HB.txt' exists.")
    if os.path.exists("./CalcSetting_HB.txt"):
        messages.append("\t>>> CalcSetting_HB.txt: Found")
        HelpList.append(False)
    else:
        messages.append(f"\t>>> {Color.RED}CalcSetting_HB.txt: NOT Found{Color.RESET}\n"
                        "\t>>> Make the correct CalcSetting_HB.txt and then restart the program.")
        HelpList.append(True)
        with open("CalcSetting_HB.txt", "w") as file:
            file.write(Stereotyped.CalcSetting_HB)

    help_check_exit(messages, HelpList)

    with open("CalcSetting_HB.txt", "r") as f:
        lines = f.readlines()
    params = {"tilt angle": ""}
    for line in lines[:5]:
        key = line.split(":")[0].strip().lower()
        if key in params:
            params[key] = line.split(":")[1].strip()
    Tilt = params["tilt angle"]
    while True:
        try:
            Tilt = float(Tilt)
            break
        except (IndexError, ValueError):
            Tilt = input("\nPlease enter the Tilt angle at degree."
                         "\n>>> ")
            continue
    print(f"\t>>> {Color.GREEN}{Tilt}-degree{Color.RESET} has been entered as Tilt")
    Formated_Tilt = int(Tilt * 10)  # For File Name
    print(f"\n{Color.RED}"
          f"*************** Caution!!! The Tilt angle multiplied by 10 is used in the file name. ***************"
          f"{Color.RESET}"
          f"\n*************** Like: {MoleculeName}_3molp1_{Color.RED}t{Formated_Tilt}d{Color.RESET}_10d-500-800.gjf\n")

    return Tilt, Formated_Tilt


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
def File_Set_Check(messages, HelpList, MaterName, Nmol, mol_pos, Formated_Tilt):
    """
    Check if the required files exist
    :param messages: list of messages to be displayed
    :type messages: list
    :param HelpList: list of help flags
    :type HelpList: list
    :param MaterName: name of the molecule
    :type MaterName: str
    :param Nmol: number of molecules
    :type Nmol: str
    :param mol_pos: position of the molecule
    :type mol_pos: str
    :param Formated_Tilt: formatted tilt angle
    :type Formated_Tilt: int
    """
    print(f"\n{Color.GREEN}File set check...{Color.RESET}")
    if os.path.exists(f"./{MaterName}.xyz"):
        messages.append(f"\t>>> {MaterName}.xyz: Found")
        HelpList.append(False)
    else:
        messages.append(f"\t{Color.RED}>>> {MaterName}.xyz: NOT Found{Color.RESET}")
        HelpList.append(True)

    if os.path.exists(f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt"):
        messages.append(f"\t>>> ./ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt: Found")
        HelpList.append(False)
    elif os.path.exists(f"./InitialCondition_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt"):
        messages.append(f"\t>>> ./InitialCondition_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt: Found")
        HelpList.append(False)
    else:
        messages.append(f"\t>>> {Color.RED}./ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt and "
                        f"./InitialCondition_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt: Not Found{Color.RESET}\n"
                        f"\t    Make the correct InitialCondition_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt "
                        f"and then restart the program.")
        with open(f"InitialCondition_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt", "w") as file:
            if "2mol" in Nmol:
                file.write(Stereotyped.InitialCondition_2mol_Temp)
            elif "3mol" in Nmol:
                file.write(Stereotyped.InitialCondition_3mol_Temp)
            else:
                file.write(Stereotyped.InitialCondition_2mol_Temp)
        HelpList.append(True)

    if "3mol" in Nmol:
        if os.path.exists(f"./{MaterName}_2mol_t{Formated_Tilt}d_min.txt"):
            messages.append(f"\t>>> {MaterName}_2mol_t{Formated_Tilt}d_min.txt: Found")
            HelpList.append(False)
        else:
            messages.append(f"\t{Color.RED}>>> {MaterName}_2mol_t{Formated_Tilt}d_min.txt: Not Found{Color.RESET}")
            HelpList.append(True)
    elif "2mol" in Nmol:
        HelpList.append(False)

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


# make a condition file
def mkConditionFile(MaterName, Nmol, mol_pos, Formated_Tilt, messages, HelpList):
    """
    Make the condition file
    :param MaterName:
    :param Nmol:
    :param mol_pos:
    :param Formated_Tilt:
    :param messages:
    :param HelpList:
    :return:
    """
    RefLines, which = [], ""
    if Nmol == "2mol":
        which = "Dcol"
    elif Nmol == "3mol":
        which = "Dtrv"
        RefLines = getRefLines(f"{MaterName}_2mol_t{Formated_Tilt}d_min.txt")

    if "3mol" in Nmol:
        dev = 0.2
    else:
        dev = 0.1

    dirpath = f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d"
    tcalpath = f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_tcal"
    os.makedirs(dirpath, exist_ok=True)
    if "3mol" in Nmol:
        os.makedirs(tcalpath, exist_ok=True)

    messages.append(f"\n{Color.GREEN}Create a GJF file...{Color.RESET}")
    if os.path.exists(f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt"):
        messages.append(f"\t>>> ./ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt: Exist.\n"
                        f"\t>>> Performs calculations based on conditions that exist.")
    else:
        NewConditions = []
        messages.append(f"\t>>> ./ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt: Not Found.\n"
                        f"\t>>> {Color.GREEN}Create a new ConditionList.txt{Color.RESET} "
                        f"based on the InitialCondition.")
        help_check_exit(messages, HelpList)
        with open(f"./InitialCondition_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt", "r") as f:
            lines = f.readlines()
        for line in lines:
            Deg = float(line.strip().split()[0])
            Val = float(line.strip().split()[1])
            RefValues = getRefValues(Deg, RefLines)
            NewCondition = mkNewCondition(Nmol, Deg, Val - (2 * dev), which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val - (1 * dev), which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val + (1 * dev), which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val + (2 * dev), which, RefValues)
            NewConditions.append(NewCondition)
        NewConditions = list(set(NewConditions))
        NewConditions.sort()
        with open(f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt", "w") as f:
            for NewCondition in NewConditions:
                f.write(NewCondition)
        messages.append(f"\t>>> ./ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt: "
                        f"{Color.GREEN}Created.{Color.RESET}")

    help_check_exit(messages, HelpList)

    return which, dev, dirpath, tcalpath, RefLines


def getRefLines(FileName):
    """
    Get the reference lines from the file
    :param FileName: name of the file
    :type FileName: str
    :return: list of reference lines
    """
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


def getRefValues(Deg, RefLines):
    """
    Get the reference values
    :param Deg: degree
    :type Deg: float
    :param RefLines: list of reference lines
    :type RefLines: list
    :return: list of reference values
    """
    RefValues = ''
    if len(RefLines) == 0:
        RefValues = ["na", "na"]
    else:
        for RefLine in RefLines:
            Contents = RefLine.strip().split('\t')
            try:
                RefDeg = round(float(Contents[0]), 1)
                if RefDeg == Deg:
                    RefValues = [Contents[1], Contents[2]]
                else:
                    pass
            except (ValueError, IndexError):
                pass

    return RefValues


def mkNewCondition(Nmol, Deg, Val, which, RefValues):
    """
    Make new condition
    :param Nmol:
    :param Deg:
    :param Val:
    :param which:
    :param RefValues:
    :return:
    """
    if "2mol" in Nmol:
        Dcol = int(round(Val * 100, 2))
        Deg = int(Deg)
        NewCondition = f"{Deg}d-{Dcol}\n"
    elif "3mol" in Nmol:
        Deg = int(Deg)
        if which == "Dcol":
            Dcol = int(round(Val * 100, 2))
            Dtrv = float(RefValues[1]) * 100
            Dtrv = int(round(Dtrv, 2))
        elif which == "Dtrv":
            Dtrv = int(round(Val * 100, 2))
            Dcol = float(RefValues[0]) * 100
            Dcol = int(round(Dcol, 2))
        else:
            Dtrv = Dcol = 0
        NewCondition = f"{Deg}d-{Dcol}-{Dtrv}\n"
    else:
        NewCondition = f""
    print(f"\t\t{NewCondition.strip()}")
    return NewCondition


# Search for the most stable structure
def Most_Stable_Structure_Search(MaterName, Nmol, mol_pos, Tilt, Formated_Tilt,
                                 dirpath, which, dev, RefLines, Operator, messages, HelpList, Debug, tcalpath, args):
    getTemporaryStructure(MaterName, Nmol, mol_pos, Formated_Tilt, dirpath, Operator, Tilt, which,
                          dev, RefLines, messages, HelpList, Debug, args)

    if "2mol" in Nmol:
        print(f"\n\t>>> {Color.GREEN}Calculations for 2mol were successfully finished.{Color.RESET}")
    elif "3mol" in Nmol:
        temp_Structures = []
        MostStable, dev = False, 0.2
        while not MostStable:
            RefLines = getRefLines(f"{MaterName}_3mol{mol_pos}_t{Formated_Tilt}d_min.txt")
            temp_Structures.append(RefLines)
            mkCycleConditions(RefLines, "Dcol", dev, Nmol, Formated_Tilt, mol_pos, MaterName)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Formated_Tilt, dirpath, Operator, Tilt, "Dcol",
                                  dev, RefLines, messages, HelpList, Debug, args)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            mkCycleConditions(RefLines, "Dtrv", dev, Nmol, Formated_Tilt, mol_pos, MaterName)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Formated_Tilt, dirpath, Operator, Tilt, "Dtrv",
                                  dev, RefLines, messages, HelpList, Debug, args)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            MostStable = CompareStructures(RefLines, temp_Structures[-1])

        print(f"\n{Color.GREEN}**********\nTransition in 0.1-Å increments.\n{Color.RESET}")
        MostStable, dev = False, 0.1
        while not MostStable:
            RefLines = getRefLines(f"{MaterName}_3mol{mol_pos}_t{Formated_Tilt}d_min.txt")
            temp_Structures.append(RefLines)
            mkCycleConditions(RefLines, "Dcol", dev, Nmol, Formated_Tilt, mol_pos, MaterName)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Formated_Tilt, dirpath, Operator, Tilt, "Dcol",
                                  dev, RefLines, messages, HelpList, Debug, args)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            mkCycleConditions(RefLines, "Dtrv", dev, Nmol, Formated_Tilt, mol_pos, MaterName)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Formated_Tilt, dirpath, Operator, Tilt, "Dtrv",
                                  dev, RefLines, messages, HelpList, Debug, args)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            MostStable = CompareStructures(RefLines, temp_Structures[-1])

        print(f"\n{Color.GREEN}**********\nTransition in 0.05-Å increments.\n{Color.RESET}")
        MostStable, dev = False, 0.05
        while not MostStable:
            RefLines = getRefLines(f"{MaterName}_3mol{mol_pos}_t{Formated_Tilt}d_min.txt")
            temp_Structures.append(RefLines)
            mkCycleConditions(RefLines, "Dcol", dev, Nmol, Formated_Tilt, mol_pos, MaterName)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Formated_Tilt, dirpath, Operator, Tilt, "Dcol",
                                  dev, RefLines, messages, HelpList, Debug, args)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            mkCycleConditions(RefLines, "Dtrv", dev, Nmol, Formated_Tilt, mol_pos, MaterName)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Formated_Tilt, dirpath, Operator, Tilt, "Dtrv",
                                  dev, RefLines, messages, HelpList, Debug, args)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            MostStable = CompareStructures(RefLines, temp_Structures[-1])
        MinConditions = getMinConditions(MaterName, Nmol, Formated_Tilt, mol_pos)

        for Condition in MinConditions:
            command = ["rename", "gjf", "com",
                       f"{dirpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_{Condition}.gjf"]
            execute(command, False)

            command = ["cp", f"{dirpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_{Condition}.com",
                       f"{tcalpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_{Condition}.com"]
            execute(command, False)

            command = ["rename", "com", "gjf",
                       f"{dirpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_{Condition}.com"]
            execute(command, False)
        if len(temp_Structures) < 100:
            with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_mins.hist", "w") as file:
                for i in range(len(temp_Structures)):
                    file.write(f"***** Structures after the '{i + 1}'th cycle *****\n")
                    Lines = temp_Structures[i]
                    for Line in Lines:
                        file.write(Line)
        else:
            pass

        print(f"\n"
              f"\t>>> {Color.GREEN}Local minimum values were successfully found "
              f"at '{len(MinConditions)}' angles.{Color.RESET}\n"
              f"\t>>> Com files for minimum energies were copied into {tcalpath} folder")
        mkXYZfile(tcalpath, Debug)
    return


def getTemporaryStructure(MaterName, Nmol, mol_pos, Formated_Tilt, dirpath, Operator, Tilt,
                          which, dev, RefLines, messages, HelpList, Debug, args):
    judge = False
    while not judge:
        qsubList = []
        Conditions = getConditions(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt")
        with open(f"{dirpath}/G.sh", "w") as f:
            f.write(Stereotyped.Sh_txt)

        for Condition in Conditions:
            if os.path.exists(f"{dirpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_{Condition}.log"):
                pass
            else:
                qsub_temp = mkFiles(MaterName, Nmol, mol_pos, Condition, Operator, dirpath,
                                    Tilt, Formated_Tilt, False, messages, HelpList)
                qsubList.append(qsub_temp)
        job_submission(messages, HelpList, qsubList, dirpath, Nmol, which)
        if Debug:
            pass
        else:
            rmWildCards(f"{dirpath}/*.sh*")
            rmWildCards(f"{dirpath}/*.chk")

        print("\n**********\nReading Data...\n")
        readEnergy(dirpath, MaterName, Nmol, Formated_Tilt, mol_pos, messages, HelpList)
        if args.energy:
            exit()
        judge = mkNewConditionLists(MaterName, Nmol, Formated_Tilt, which, dev, RefLines, mol_pos)
    return


def getConditions(FileName):
    """
    Get the conditions from the file
    :param FileName: name of the file
    :type FileName: str
    :return: list of conditions
    """
    Condition = []
    with open(FileName, "r") as f:
        lines = f.readlines()
    for line in lines:
        Condition.append(line.strip())
    return Condition


def mkFiles(MaterName, Nmol, mol_pos, Condition, Operator, dirpath,
            Tilt_Angle, Formated_Tilt, Flag_XYZ, messages, HelpList):
    """
    Generate input files for molecular simulations based on provided parameters.

    This function creates various input files required for molecular simulations,
    including Gaussian job files (.gjf), checkpoint files (.chk), and shell
    script files (.sh). The generated files are tailored to the specific
    conditions and molecular configurations provided by the user.

    :param str MaterName: Name of the material or molecule.
    :param str Nmol: The number of molecules involved in the simulation
                     (e.g., '2mol', '3mol').
    :param str mol_pos: The position type of the molecules.
    :param str Condition: Condition string encoding the rotation angle,
                          column direction displacement, and transverse
                          direction displacement.
    :param str Operator: Operator name or identifier for the simulation.
    :param str dirpath: Directory path where the generated files will be saved.
    :param float Tilt_Angle: Angle of tilt applied to the molecules.
    :param int Formated_Tilt: Formatted string of the tilt angle for file naming.
    :param bool Flag_XYZ: Flag to determine whether XYZ coordinate files should
                          be generated.
    :param list messages: List of messages to be displayed.
    :param list HelpList: List of help flags.

    :returns: Command string for submitting the generated shell script
              to a job scheduler.
    :rtype: str
    """
    File_Name = f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_{Condition}"
    CHK_Name = File_Name + ".chk"
    GJF_Name = File_Name + ".gjf"
    SH_Name = f"G-{Operator}_t{Formated_Tilt}d_{Condition}.sh"

    ConditionList = Condition.strip().split("-")
    Rotate_Angle = float(ConditionList[0].replace("d", ""))
    D_Col = float(ConditionList[1]) / 100
    if Nmol == "3mol":
        D_Transv = float(ConditionList[2]) / 100
    else:
        D_Transv = 0

    Direction_Col, Direction_Transv, Rotate_Axis, Tilt_Axis, rotate = Axis_Setting_HB(messages, HelpList)
    Angles = {
        "": {
            "x": {"y": [[Rotate_Angle, Tilt_Angle, 0], [Rotate_Angle, Tilt_Angle, 0], [-Rotate_Angle, Tilt_Angle, 0]],
                  "z": [[Rotate_Angle, 0, Tilt_Angle], [Rotate_Angle, 0, Tilt_Angle], [-Rotate_Angle, 0, Tilt_Angle]]},
            "y": {"x": [[Tilt_Angle, Rotate_Angle, 0], [Tilt_Angle, Rotate_Angle, 0], [Tilt_Angle, -Rotate_Angle, 0]],
                  "z": [[0, Rotate_Angle, Tilt_Angle], [0, Rotate_Angle, Tilt_Angle], [0, -Rotate_Angle, Tilt_Angle]]},
            "z": {"x": [[Tilt_Angle, 0, Rotate_Angle], [Tilt_Angle, 0, Rotate_Angle], [Tilt_Angle, 0, -Rotate_Angle]],
                  "y": [[0, Tilt_Angle, Rotate_Angle], [0, Tilt_Angle, Rotate_Angle], [0, Tilt_Angle, -Rotate_Angle]]}
        },
        "p1": {
            "x": {"y": [[Rotate_Angle, Tilt_Angle, 0], [Rotate_Angle, Tilt_Angle, 0], [-Rotate_Angle, Tilt_Angle, 0]],
                  "z": [[Rotate_Angle, 0, Tilt_Angle], [Rotate_Angle, 0, Tilt_Angle], [-Rotate_Angle, 0, Tilt_Angle]]},
            "y": {"x": [[Tilt_Angle, Rotate_Angle, 0], [Tilt_Angle, Rotate_Angle, 0], [Tilt_Angle, -Rotate_Angle, 0]],
                  "z": [[0, Rotate_Angle, Tilt_Angle], [0, Rotate_Angle, Tilt_Angle], [0, -Rotate_Angle, Tilt_Angle]]},
            "z": {"x": [[Tilt_Angle, 0, Rotate_Angle], [Tilt_Angle, 0, Rotate_Angle], [Tilt_Angle, 0, -Rotate_Angle]],
                  "y": [[0, Tilt_Angle, Rotate_Angle], [0, Tilt_Angle, Rotate_Angle], [0, Tilt_Angle, -Rotate_Angle]]}
        },
        "p2": {
            "x": {"y": [[-Rotate_Angle, Tilt_Angle, 0], [-Rotate_Angle, Tilt_Angle, 0], [Rotate_Angle, Tilt_Angle, 0]],
                  "z": [[-Rotate_Angle, 0, Tilt_Angle], [-Rotate_Angle, 0, Tilt_Angle], [Rotate_Angle, 0, Tilt_Angle]]},
            "y": {"x": [[Tilt_Angle, -Rotate_Angle, 0], [Tilt_Angle, -Rotate_Angle, 0], [Tilt_Angle, Rotate_Angle, 0]],
                  "z": [[0, -Rotate_Angle, Tilt_Angle], [0, -Rotate_Angle, Tilt_Angle], [0, Rotate_Angle, Tilt_Angle]]},
            "z": {"x": [[Tilt_Angle, 0, -Rotate_Angle], [Tilt_Angle, 0, -Rotate_Angle], [Tilt_Angle, 0, Rotate_Angle]],
                  "y": [[0, Tilt_Angle, -Rotate_Angle], [0, Tilt_Angle, -Rotate_Angle], [0, Tilt_Angle, Rotate_Angle]]}
        },
        "p3": {
            "x": {"y": [[Rotate_Angle, Tilt_Angle, 0], [Rotate_Angle, Tilt_Angle, 0],
                        [-180 - Rotate_Angle, Tilt_Angle, 0]],
                  "z": [[Rotate_Angle, 0, Tilt_Angle], [Rotate_Angle, 0, Tilt_Angle],
                        [-180 - Rotate_Angle, 0, Tilt_Angle]]},
            "y": {"x": [[Tilt_Angle, Rotate_Angle, 0], [Tilt_Angle, Rotate_Angle, 0],
                        [Tilt_Angle, -180 - Rotate_Angle, 0]],
                  "z": [[0, Rotate_Angle, Tilt_Angle], [0, Rotate_Angle, Tilt_Angle],
                        [0, -180 - Rotate_Angle, Tilt_Angle]]},
            "z": {"x": [[Tilt_Angle, 0, Rotate_Angle], [Tilt_Angle, 0, Rotate_Angle],
                        [Tilt_Angle, 0, -180 - Rotate_Angle]],
                  "y": [[0, Tilt_Angle, Rotate_Angle], [0, Tilt_Angle, Rotate_Angle],
                        [0, Tilt_Angle, -180 - Rotate_Angle]]}
        }
    }.get(mol_pos, {}).get(Rotate_Axis, {}).get(Direction_Col, [])

    col_transl = mkDirection(D_Col, Direction_Col, "column direction", messages, HelpList)
    transv_transl = mkDirection(D_Transv, Direction_Transv, "transverse direction", messages, HelpList)
    Direction_Other = "xyz".replace(Direction_Col, "").replace(Direction_Transv, "")
    a_transl = mkDirection(Constant.d_another, Direction_Other, "other direction", messages, HelpList)
    Transitions = [col_transl, transv_transl, a_transl]

    Element, Mol1_pos, Mol2_pos, Mol3_pos, NinMol = mkAtomList(MaterName, Angles[0], Angles[1],
                                                               Angles[2], rotate, Transitions)

    if "2mol" in Nmol:
        with open(f"Temp_Header_2mol.txt", "w") as temp_header:
            temp_header.write(Stereotyped.Header_2mol)
        with open(f"Temp_Header_2mol.txt", "r") as temp_header:
            Headers = temp_header.readlines()
        os.remove(f"Temp_Header_2mol.txt")
        Headers[3] = f"%chk={CHK_Name}\n"
    elif "3mol" in Nmol:
        with open(f"Temp_Header_3mol.txt", "w") as temp_header:
            temp_header.write(Stereotyped.Header_3mol)
        with open(f"Temp_Header_3mol.txt", "r") as temp_header:
            Headers = temp_header.readlines()
        os.remove(f"Temp_Header_3mol.txt")
        Headers[3] = f"%chk={CHK_Name}\n"

    if "2mol" in Nmol:
        write_gjf_file(f"{dirpath}/{File_Name}.gjf",
                       Headers, Element, Mol1_pos, Mol2_pos)
        if Flag_XYZ:
            write_xyz_file(f"{dirpath}/{File_Name}.xyz",
                           Element, Mol1_pos, Mol2_pos)
            messages.append(f"\t>>> {File_Name}.xyz: {Color.GREEN}Complete.{Color.RESET}")
    elif "3mol" in Nmol:
        write_gjf_file(f"{dirpath}/{File_Name}.gjf",
                       Headers, Element, Mol1_pos, Mol2_pos, Mol3_pos)
        if Flag_XYZ:
            write_xyz_file(f"{dirpath}/{File_Name}.xyz",
                           Element, Mol1_pos, Mol2_pos, Mol3_pos)
            messages.append(f"\t>>> {File_Name}.xyz: {Color.GREEN}Complete.{Color.RESET}")
    with open(f"{dirpath}/G.sh", "r") as orgSH:
        lines = orgSH.readlines()
    lines[12] = f"g16 {GJF_Name}\n"
    with open(f"{dirpath}/{SH_Name}", "w") as newSH:
        for line in lines:
            newSH.write(line)
    qsub_temp = f"qsub {SH_Name}"
    return qsub_temp


def Axis_Setting_HB(messages, HelpList):
    """
    Set the axis for the HerringBone structure
    :param messages: list of messages to be displayed
    :param HelpList: list of help flags
    :return: column direction, transverse direction, other axis, tilt axis, rotate
    """
    try:
        with open("CalcSetting_HB.txt", "r") as File:
            lines = File.readlines()
    except FileNotFoundError:
        messages.append(f"{Color.RED}CalcSetting_HB.txt: NOT Found{Color.RESET}")
        HelpList.append(True)
        with open("CalcSetting_HB.txt", "w") as file:
            file.write(Stereotyped.CalcSetting_HB)
        messages.append(f"{Color.RED}CalcSetting_HB.txt: Created{Color.RESET}")
    help_check_exit(messages, HelpList)

    params = {"column direction": "", "transverse direction": "", "Tilt Angle": "",
              "tilt axis": "z", "other axis": "y"}
    for line in lines[:5]:
        key = line.split(":")[0].strip().lower()
        if key in params:
            params[key] = line.split(":")[1].strip()
    if params["column direction"] == "x":
        if params["transverse direction"] == "z":
            params["other axis"] = "y"
        elif params["transverse direction"] == "y":
            params["other axis"] = "z"
    if params["column direction"] == "y":
        if params["transverse direction"] == "x":
            params["other axis"] = "z"
        elif params["transverse direction"] == "z":
            params["other axis"] = "x"
    if params["column direction"] == "z":
        if params["transverse direction"] == "x":
            params["other axis"] = "y"
        elif params["transverse direction"] == "y":
            params["other axis"] = "x"

    axis1 = params["other axis"]
    axis2 = params["transverse direction"]
    axis3 = params["column direction"]
    rotate = f"{axis1}{axis2}{axis3}"

    if params["column direction"] == params["transverse direction"]:
        messages.append(f"{Color.RED}Error: Column direction and transverse direction are the same.{Color.RESET}")
        HelpList.append(True)
    if params["column direction"] == params["other axis"]:
        messages.append(f"{Color.RED}Error: Column direction and other axis are the same.{Color.RESET}")
        HelpList.append(True)
    if params["transverse direction"] == params["other axis"]:
        messages.append(f"{Color.RED}Error: Transverse direction and other axis are the same.{Color.RESET}")
        HelpList.append(True)
    if params["tilt axis"] not in {"x", "y", "z"}:
        messages.append(f"{Color.RED}Error: Invalid tilt axis.{Color.RESET}")
        HelpList.append(True)
    if rotate not in {"xyz", "xzy", "yxz", "yzx", "zxy", "zyx"}:
        messages.append(f"{Color.RED}Error: Invalid rotation axis.{Color.RESET}")
        HelpList.append(True)
    help_check_exit(messages, HelpList)

    return params["column direction"], params["transverse direction"], params["other axis"], params["tilt axis"], rotate


def mkDirection(axis_direction, input_axis, axis_name, messages, HelpList):
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
    help_check_exit(messages, HelpList)
    return np.array(directions[input_axis])


def mkAtomList(MaterName, Mol1_Angles, Mol2_Angles, Mol3_Angles, rotate, Translations):
    """
    Make the atom list
    :param MaterName:
    :param Mol1_Angles:
    :param Mol2_Angles:
    :param Mol3_Angles:
    :param rotate:
    :param Translations:
    :return:
    """
    with open(f"./{MaterName}.xyz", "r") as f:
        NinMol = int(f.readline())
        f.readline()
        AtomList = f.readlines()
    Mol1, Mol2, Mol3, Elist = [], [], [], []
    for Atom in AtomList:
        Contents = Atom.split()
        Elist.append(Contents.pop(0))
        Position = np.array(list(map(float, Contents)))
        atm_m1 = Rotate(Position, *Mol1_Angles, rotate)
        Mol1.append(atm_m1)
        atm_m2 = Rotate(Position, *Mol2_Angles, rotate) + Translations[0]
        Mol2.append(atm_m2)
        atm_m3 = Rotate(Position, *Mol3_Angles, rotate) + Translations[0] / 2 + Translations[1]
        Mol3.append(atm_m3)
    return Elist, Mol1, Mol2, Mol3, NinMol


def Rotate(Current, Tx, Ty, Tz, rotation):
    """
    Rotate the current position
    :param Current:
    :param Tx:
    :param Ty:
    :param Tz:
    :param rotation:
    :return:
    """
    Tx, Ty, Tz = map(math.radians, [Tx, Ty, Tz])
    Rx = np.array([[1, 0, 0], [0, math.cos(Tx), -math.sin(Tx)], [0, math.sin(Tx), math.cos(Tx)]])
    Ry = np.array([[math.cos(Ty), 0, math.sin(Ty)], [0, 1, 0], [-math.sin(Ty), 0, math.cos(Ty)]])
    Rz = np.array([[math.cos(Tz), -math.sin(Tz), 0], [math.sin(Tz), math.cos(Tz), 0], [0, 0, 1]])
    rotation_matrices = {'x': Rx, 'y': Ry, 'z': Rz}
    R = np.eye(3)
    for axis in rotation:
        R = np.dot(rotation_matrices[axis], R)
    return np.dot(R, Current)


def write_gjf_file(filename, headers, elements, positions1, positions2=None, positions3=None):
    """
    Write the Gaussian job file
    :param filename:
    :param headers:
    :param elements:
    :param positions1:
    :param positions2:
    :param positions3:
    :return:
    """
    with open(filename, "w") as file:
        for header in headers:
            file.write(header)
        for elem, pos in zip(elements, positions1):
            file.write(f" {elem:<2}  {format_coordinate(pos)}  1\n")
        if positions2 is not None:
            for elem, pos in zip(elements, positions2):
                file.write(f" {elem:<2}  {format_coordinate(pos)}  2\n")
        if positions3 is not None:
            for elem, pos in zip(elements, positions3):
                file.write(f" {elem:<2}  {format_coordinate(pos)}  3\n")
        file.write("\n")
    return


def write_xyz_file(filename, elements, positions1, positions2=None, positions3=None):
    """
    Write the XYZ file
    :param filename:
    :param elements:
    :param positions1:
    :param positions2:
    :param positions3:
    :return:
    """
    with open(filename, "w") as file:
        if positions3 is None:
            file.write(f"{len(elements) * 2}\n")
        else:
            file.write(f"{len(elements) * 3}\n")
        file.write("00000001\n")
        for elem, pos in zip(elements, positions1):
            file.write(f" {elem:<2}  {format_coordinate(pos)}\n")
        if positions2 is not None:
            for elem, pos in zip(elements, positions2):
                file.write(f" {elem:<2}  {format_coordinate(pos)}\n")
        if positions3 is not None:
            for elem, pos in zip(elements, positions3):
                file.write(f" {elem:<2}  {format_coordinate(pos)}\n")
        file.write("\n")
    return


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
def job_submission(messages, HelpList, qsubList, dirpath, Nmol, which):
    """
    Submit the job
    :param messages:
    :param HelpList:
    :param qsubList:
    :param dirpath:
    :param Nmol:
    :param which:
    :return:
    """
    print("\n**********\nJobs are submitting...")
    if len(qsubList) == 0:
        messages.append(f"\t>>> Job was not submitted.\n"
                        f"\t>>> {Color.GREEN}Calculations with the conditions might be finished.{Color.RESET}")
        help_check_exit(messages, HelpList)
    else:
        for qsub in qsubList:
            qsub = qsub.split()
            subprocess.run(qsub, cwd=f"./{dirpath}")

        if "2mol" in Nmol:
            Wait_minutes = 1
        elif "3mol" in Nmol:
            Wait_minutes = 2
        else:
            Wait_minutes = 1
        MyJobIDList = My_JobIDList(qsubList)
        MyJobIDList.sort()
        RunningJobIDList = Running_JobIDList()

        Flag, wait_job_count = check_jobs(RunningJobIDList, MyJobIDList)
        start_time = datetime.datetime.now()

        if Flag:
            formated_ST = start_time.strftime("%m/%d %H:%M:%S")
            if which == "Dcol":
                term = "the stable distance in column direction"
            elif which == "Dtrv":
                term = "the stable distance in transverse direction"
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
            Flag, wait_job_count = check_jobs(RunningJobIDList, MyJobIDList)
            formated_NOW, elapsed_time = getElapsedTime(start_time)
            Now_Time = datetime.datetime.now()
            formated_end_time = (
                (Now_Time + datetime.timedelta(minutes=(wait_job_count * Wait_minutes))).strftime("%m/%d %H:%M:%S"))
            sys.stdout.write(
                "\033[1F\033[G%s" %
                f"\t{formated_NOW} ({elapsed_time} min. passed): Job {RunningJobIDList[0]} is in progress.           \n"
                f"\tNext Check >>> {Wait_minutes * wait_job_count} minute later! forecast: {formated_end_time}    ")
            sys.stdout.flush()
            time.sleep(Wait_minutes * wait_job_count * 60)
            RunningJobIDList = Running_JobIDList()
            Flag, wait_job_count = check_jobs(RunningJobIDList, MyJobIDList)
        else:
            print(f"{Color.GREEN}\n\n"
                  f"Calculation cycles for {which} until JobID {MyJobIDList[-1]} were finished.{Color.RESET}")
    return


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


def readEnergy(dirpath, MaterName, Nmol, Formated_Tilt, mol_pos, messages, HelpList):
    """
    Read the energy
    :param dirpath:
    :param MaterName:
    :param Nmol:
    :param Formated_Tilt:
    :param mol_pos:
    :param messages:
    :param HelpList:
    :return:
    """
    DegList = []
    LogList = glob.glob(f"{dirpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_*.log")
    LogList.sort()
    for Log in LogList:
        FileName, deg, Vdcol, Vdtrv = getVAL_fromLogName(Nmol, Log, messages, HelpList)
        DegList.append(deg)
    DegList = list(set(DegList))
    DegList.sort()

    with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_all.txt", "w") as AllData:
        with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt", "w") as MinData:
            header = ("Angle \tDistance in column direction (Å)\tDistance in transverse direction (Å)"
                      "\tCounterpoise corrected energy (A.U)\tBSSE energy (A.U)")
            AllData.write(f"*****  {MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d All Results *****\n")
            MinData.write(f"***** {MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d Minimum Energy at each Angle *****\n")
            AllData.write(f" \t{header}\n")
            MinData.write(f"{header}\n")
            MinimumL = []

            print(f"\t{header}")
            for Deg in DegList:
                CPE_Dict = {}
                VAL_Dict = {}
                Keys = []
                AllData.write("******\t******\t******\t******\t******\t******\n")
                for Log in LogList:
                    Condition = getCondition_fromName(Log)
                    deg = Condition[:Condition.find("d-")]
                    Vdeg = round(float(deg), 5)
                    if Vdeg == Deg:
                        with open(Log, "r") as file:
                            data = file.read()
                        if "Normal termination" in data:
                            FileName, Vdeg, Vdcol, Vdtrv = getVAL_fromLogName(Nmol, Log, messages, HelpList)
                            CPE, BSE = getEnergy(data)
                            CPE_Dict[FileName] = CPE
                            VAL_Dict[FileName] = f"{Vdeg}\t{Vdcol}\t{Vdtrv}\t{CPE}\t{BSE}\n"
                            Keys.append(FileName)
                        else:
                            print(f"{Color.RED}Error: {Log} did not finish normally.{Color.RESET}")
                            LogFileName = Log
                            with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_error.log", "a") as error_log:
                                Log = Log.split("/")[-1]
                                Log = Log.split(".")[0]
                                Log = Log.split("_")[-1]
                                error_log.write(f"{Log}\n")
                            with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_error.log", "r") as error_log:
                                error_lines = error_log.readlines()
                            if len(error_lines) != len(set(error_lines)):
                                print(f"{Color.RED}Error: {Log} is duplicated.{Color.RESET}")
                                messages.append(f"{Color.RED}Error: {Log} is duplicated.{Color.RESET}")
                                HelpList.append(True)
                            help_check_exit(messages, HelpList)
                            rmWildCards(LogFileName)

                minkey = min(CPE_Dict, key=CPE_Dict.get)
                MinimumL.append(minkey)
                MinData.write(VAL_Dict.get(minkey))
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
    return


def getVAL_fromLogName(Nmol, Log, messages, HelpList):
    deg, dcol, dtrv = 0, 0, 0
    Condition = getCondition_fromName(Log)
    FileName = Log[Log.rfind("/") + 1:-4]

    Conditions = Condition.split('-')

    if "2mol" in Nmol and Condition.count("-") == 1:
        deg = Conditions[0].replace('d', '')
        dcol = Conditions[1]
        dtrv = "na"
    elif "3mol" in Nmol and Condition.count("-") == 2:
        deg = Conditions[0].replace('d', '')
        dcol = Conditions[1]
        dtrv = Conditions[2]
    else:
        messages.append(f"{Color.RED}Error: Invalid Condition.{Color.RESET}")
        HelpList.append(True)

    help_check_exit(messages, HelpList)

    Vdeg = round(float(deg), 5)
    Vdcol = round(float(dcol) / 100, 5)
    try:
        Vdtrv = round(float(dtrv) / 100, 5)
    except ValueError:
        Vdtrv = dtrv
    return FileName, Vdeg, Vdcol, Vdtrv


def getCondition_fromName(Name):
    Condition = Name[Name.rfind("_") + 1:Name.rfind(".")]
    return Condition


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


def mkNewConditionLists(MaterName, Nmol, Formated_Tilt, which, dev, RefLines, mol_pos):
    with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_all.txt", "r") as All:
        AllLines = All.readlines()

    DegList = []
    for AllLine in AllLines:
        Contents = AllLine.strip().split('\t')
        try:
            Deg = round(float(Contents[1]), 1)
            DegList.append(Deg)
        except IndexError:
            pass
        except ValueError:
            pass
    DegList = list(set(DegList))
    DegList.sort()

    Judges, ComplDeg, ComplVal, NewConditions = [], [], [], []

    for Deg in DegList:
        SV = 0
        RefValues = getRefValues(Deg, RefLines)
        ValueList = []
        if "3mol" in Nmol and which == "Dcol":
            SV = round(float(RefValues[1]), 2)
        elif "3mol" in Nmol and which == "Dtrv":
            SV = round(float(RefValues[0]), 2)

        for AllLine in AllLines:
            Contents = AllLine.strip().split("\t")

            try:
                DataDeg = round(float(Contents[1]), 1)
                if DataDeg == Deg and "2mol" in Nmol:
                    DataDcol = round(float(Contents[2]), 2)
                    ValueList.append(DataDcol)
                    if Contents[0] == "*":
                        SV = DataDcol
                elif DataDeg == Deg and "3mol" in Nmol:
                    DataDcol = round(float(Contents[2]), 2)
                    DataDtrv = round(float(Contents[3]), 2)

                    if which == "Dcol":
                        RefDtrv = round(float(RefValues[1]), 2)
                        if DataDtrv == RefDtrv:
                            ValueList.append(DataDcol)
                        if Contents[0] == "*":
                            SV = DataDcol
                    if which == "Dtrv":
                        RefDcol = round(float(RefValues[0]), 2)
                        if DataDcol == RefDcol:
                            ValueList.append(DataDtrv)
                        if Contents[0] == "*":
                            SV = DataDtrv
            except ValueError:
                pass
            except IndexError:
                pass

        if round(SV + 0.05, 2) in ValueList and round(SV - 0.05, 2) in ValueList:
            Judges.append("complete")
            print(f"\nThe local minimum by 0.05 step for {Deg} degree was Found;\t{SV}.")
            print(f"{Color.GREEN}\tCOMPLETE!!{Color.RESET}")
            ComplDeg.append(Deg)
            ComplVal.append(SV)
        elif round(SV + dev, 2) in ValueList and round(SV - dev, 2) in ValueList:
            Judges.append("complete")
            print(f"\nThe local minimum by {dev} step for {Deg} degree was Found;\t{SV}.")
            print(f"{Color.GREEN}\tCOMPLETE!!{Color.RESET}")
            ComplDeg.append(Deg)
            ComplVal.append(SV)
        elif SV == min(ValueList):
            Judges.append("not complete")
            print(f"\nLocal minimum for {Deg} degree was NOT FOUND in the cycle.")
            print(f"\tNOT COMPLETE(1)!! {Constant.Cn} new conditions bellow were appended to the ConditionList.txt.")
            for i in range(Constant.Cn):
                NewCondition = mkNewCondition(Nmol, Deg, SV - dev * (i + 1), which, RefValues)
                NewConditions.append(NewCondition)
        elif SV == max(ValueList):
            Judges.append("not complete")
            print(f"\nLocal minimum for {Deg} degree was NOT FOUND in the cycle.")
            print(f"\tNOT COMPLETE(2)!! {Constant.Cn} new conditions bellow were appended to the ConditionList.txt.")
            for i in range(Constant.Cn):
                NewCondition = mkNewCondition(Nmol, Deg, SV + dev * (i + 1), which, RefValues)
                NewConditions.append(NewCondition)
    with open(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt", "r") as f:
        orgCondition = f.readlines()
    NewList = orgCondition + NewConditions
    NewList = list(set(NewList))
    NewList.sort()

    with open(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt", "w") as f:
        for content in NewList:
            f.write(content)
    if len(ComplDeg) == 0:
        print("\nLocal minimum have NOT been FOUND in any angle.\n")
    else:
        print(f"\nLocal minimum was FOUND in {len(ComplDeg)} angles.\n")
        if "2mol" in Nmol or "3mol" in Nmol:
            print(f"Angle\t{which}")
        else:
            pass
        for i in range(len(ComplDeg)):
            print(f"{ComplDeg[i]}\t{ComplVal[i]}")

    Judges = list(set(Judges))
    if "not complete" in Judges:
        judge = False
    elif "complete" in Judges and len(Judges) == 1:
        judge = True
    else:
        judge = False

    return judge


def mkCycleConditions(RefLines, which, dev, Nmol, Formated_Tilt, mol_pos, MaterName):
    Degs = []
    NewConditions = []
    if dev == 0.2:
        n = Constant.CycleCondition_n_02
    elif dev == 0.1:
        n = Constant.CycleCondition_n_01
    elif dev == 0.05:
        n = Constant.CycleCondition_n_005
    else:
        n = 1

    for RefLine in RefLines:
        Contents = RefLine.strip().split()
        Degs.append(round(float(Contents[0]), 1))
    print("\nCreating new conditions for the next cycle...")
    print("\tNew conditions for the next cycle:")
    for Deg in Degs:
        ValueList = []
        RefValues = getRefValues(Deg, RefLines)
        if which == "Dcol":
            SV = round(float(RefValues[0]), 2)
        elif which == "Dtrv":
            SV = round(float(RefValues[1]), 2)
        else:
            SV = int()

        with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_all.txt", "r") as All:
            AllLines = All.readlines()
            for AllLine in AllLines:
                Contents = AllLine.strip().split("\t")
                try:
                    DataDeg = round(float(Contents[1]), 1)
                    if DataDeg == Deg:
                        DataDcol = round(float(Contents[2]), 2)
                        DataDtrv = round(float(Contents[3]), 2)

                        if which == "Dcol":
                            RefDtrv = round(float(RefValues[1]), 2)
                            if DataDtrv == RefDtrv:
                                ValueList.append(DataDcol)
                            if Contents[0] == "*":
                                SV = DataDcol
                        if which == "Dtrv":
                            RefDcol = round(float(RefValues[0]), 2)
                            if DataDcol == RefDcol:
                                ValueList.append(DataDtrv)
                            if Contents[0] == "*":
                                SV = DataDtrv
                except ValueError:
                    pass
                except IndexError:
                    pass

        if round(SV + 0.05, 2) in ValueList and round(SV - 0.05, 2) in ValueList:
            pass
        else:
            for i in range(n):
                NewCondition = mkNewCondition(Nmol, float(int(Deg)), SV - dev * (i + 1), which, RefValues)
                NewConditions.append(NewCondition)
                NewCondition = mkNewCondition(Nmol, float(int(Deg)), SV + dev * (i + 1), which, RefValues)
                NewConditions.append(NewCondition)

    with open(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt", "r") as file:
        orgCondition = file.readlines()
    NewList = orgCondition + NewConditions
    NewList = list(set(NewList))
    NewList.sort()
    with open(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt", "w") as file:
        for content in NewList:
            file.write(content)
    return


def CompareStructures(RefLines, LinesBefore):
    for RefLine in RefLines:
        RefContents = RefLine.strip().split()
        RefDeg, RefDcol, RefDtrv = RefContents[0], RefContents[1], RefContents[2]
        for LineBefore in LinesBefore:
            Contents = LineBefore.strip().split()
            Deg = Contents[0]
            if Deg == RefDeg:
                if Contents[1] != RefDcol or Contents[2] != RefDtrv:
                    return False
    return True


def getMinConditions(MaterName, Nmol, Formated_Tilt, mol_pos):
    with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt") as file:
        lines = file.readlines()
    del lines[0:2]

    MinConditions = []
    for line in lines:
        contents = line.strip().split()
        Deg = int(float(contents[0]))
        Dcol = int(round(float(contents[1]) * 100, 5))
        if contents[2] == "na":
            Condition = f"{Deg}d-{Dcol}"
        else:
            Dtrsv = int(round(float(contents[2]) * 100, 5))
            Condition = f"{Deg}d-{Dcol}-{Dtrsv}"
        MinConditions.append(Condition)
    return MinConditions


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


def mkXYZfile(tcalpath, Debug):
    """
    Make the XYZ file
    :param tcalpath:
    :param Debug:
    :return:
    """
    F_paths = glob.glob(f"{tcalpath}/*.com")
    print(f"\n"
          f"**********\n"
          f"Making XYZ files...")
    if not F_paths:
        print(f"\t>>> There is NO .com file in the {tcalpath} folder.")
        return

    for F_path in F_paths:
        process_file(F_path)

    ComFileNames = glob.glob(f"{tcalpath}/*.com")
    XYZFileNameList = []
    for ComFileName in ComFileNames:
        FormatedComfn = (ComFileName[:-4]).strip().split("/")[-1]
        XYZFileName = f"{FormatedComfn}.xyz"
        if Debug:
            print(f"{ComFileName} -> {XYZFileName}")
        XYZFileNameList.append(XYZFileName)
        subprocess.run(["newzmat", "-icart", "-oxyz", FormatedComfn, XYZFileName], timeout=2, cwd=tcalpath)
    Ls = []
    for XYZFileName in XYZFileNameList:
        with open(f"{tcalpath}/{XYZFileName}", "r") as xyzf:
            lines = xyzf.readlines()
        Ls.append(len(lines))
        with open(f"{tcalpath}/{XYZFileName}", "w") as xyzf:
            xyzf.write(f"{len(lines)}\n")
            xyzf.write("00000001\n")
            for line in lines:
                xyzf.write(line)
    print(f"\t>>> {tcalpath}/*.xyz: {Color.GREEN}Created!!{Color.RESET} ({len(Ls)} files)")
    return


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


def Calculate_TI(calculation_tcal_flag, tcal_path, MaterName, Nmol, mol_pos,
                 Formated_Tilt, Debug, messages, HelpList):
    if calculation_tcal_flag or "2mol" in Nmol:
        pass
    else:
        print(f"\n**********\n"
              f"{Color.GREEN}Calculating transfer integrals...\n{Color.RESET}")
        if not os.path.exists(f"./{tcal_path}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_tcal.log"):
            qsubList = []
            XYZs = glob.glob(f"{tcal_path}/*.xyz")
            for XYZ in XYZs:
                if "_m1.xyz" in XYZ or "_m2.xyz" in XYZ or "-12.xyz" in XYZ or "-23.xyz" in XYZ or "-31.xyz" in XYZ:
                    XYZ = XYZ.replace("./", "")
                    os.remove(f"{XYZ}")
                else:
                    pass
            XYZ_3mol_to_XYZ_2mol(tcal_path, Debug, messages, HelpList)

            with open(f"{tcal_path}/tcal.sh", "w") as f:
                f.write(Stereotyped.tcal_sh_txt)
            qsubList.append("qsub tcal.sh")
            job_submission(messages, HelpList, qsubList, tcal_path, Nmol, "tcal")
            if Debug:
                pass
            else:
                rmWildCards(f"{tcal_path}/*.sh*")
            subprocess.run(["rename", "tcal", f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_tcal", "tcal.log"],
                           cwd=tcal_path)
            readlog(tcal_path, MaterName, Nmol, mol_pos, Formated_Tilt, messages, HelpList)
        else:
            print(f"\t>>> tcal.log: {Color.GREEN}Already exists!!{Color.RESET}")
            print(f"\t>>> {Color.GREEN}Calculation of transfer integrals was skipped.{Color.RESET}")
    # Phase check
    PhaseCheck(tcal_path, Debug)
    return


def XYZ_3mol_to_XYZ_2mol(tcal_path, Debug, messages, HelpList):
    filepaths = glob.glob(f"{tcal_path}/*_3mol*.xyz")

    Faults = []
    print("\nConverting 3mol to 2mol XYZ files...")
    for filepath in filepaths:
        Mol1, Mol2, Mol3 = [], [], []
        filecore = filepath[:-4]
        with open(filepath, "r") as f:
            # Read the first line ( the number of atoms )
            number_atoms = f.readline()
            Comment = f.readline()
            coordinates = f.readlines()
        if Debug:
            messages.append("\t>>> " + filepath.strip())
            messages.append(f"\t>>> Number of atoms: {number_atoms.strip()}")

        if int(float(number_atoms)) % 3 != 0:
            Faults.append(filepath)
            messages.append(f"\t>>> {Color.RED}Error: {filepath} did not divide by 3.{Color.RESET}")
        else:
            Atoms_inMol = int(round(float(number_atoms) / 3, 5))
            for i in range(Atoms_inMol):
                Mol1.append(coordinates[i])
                Mol2.append(coordinates[i + Atoms_inMol])
                Mol3.append(coordinates[i + 2 * Atoms_inMol])
            Atoms_in_TcalXYZFile = int(round(Atoms_inMol * 2, 3))

            with open(f"{filecore}-12.xyz", "w") as newF12, open(f"{filecore}-23.xyz", "w") as newF23, open(
                    f"{filecore}-31.xyz", "w") as newF31:
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
            messages.append(f"\t>>> {filecore}-12.xyz: Created!!")
            messages.append(f"\t>>> {filecore}-23.xyz: Created!!")
            messages.append(f"\t>>> {filecore}-31.xyz: Created!!")
            os.rename(filepath, f"{filecore}.all")
    if len(Faults) != 0:
        for Fault in Faults:
            messages.append(f"\t>>> {Color.RED}Error: {Fault} did not divide by 3.{Color.RESET}")
        HelpList.append(True)
    else:
        messages.append(f"\n\t>>> {Color.GREEN}XYZ 3mol to 2mol: Succeeded!!{Color.RESET}")
    help_check_exit(messages, HelpList)
    return


def readlog(tcal_path, MaterName, Nmol, mol_pos, Formated_Tilt, messages, HelpList):
    print(f"\nReading the Tcal log file...")
    filepath = f"{tcal_path}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_tcal.log"
    with open(filepath, "r") as file:
        lines = file.readlines()

    DataList = []
    keywords = ["Input File Name:", "NLUMO", "LUMO", "HOMO", "NHOMO"]
    for line in lines:
        for keyword in keywords:
            if keyword in line:
                DataList.append(line)
                break

    if len(DataList) % 5 != 0:
        print(f"\t>>> {Color.RED}Error: UNEXPECTED ERROR in read Tcal log process.{Color.RESET}")
        HelpList.append(True)
        help_check_exit(messages, HelpList)

    output_filepath = f"{tcal_path}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_TIs.txt"
    with open(output_filepath, "w") as file:
        header = (f"***** Transfer Integrals in "
                  f"{tcal_path}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_tcal.log *****\n")
        columns = "input file\tNLUMO (meV)\tLUMO (meV)\tHOMO (meV)\tNHOMO (meV)\n"
        file.write(header + columns)
        print(header.strip())
        print(columns.strip())

        for i in range(0, len(DataList), 5):
            input_file = DataList[i].split(": ")[1].split(".xyz")[0]
            NLUMO = extract_value(DataList[i + 1], "NLUMO")
            LUMO = extract_value(DataList[i + 2], "LUMO")
            HOMO = extract_value(DataList[i + 3], "HOMO")
            NHOMO = extract_value(DataList[i + 4], "NHOMO")
            file.write(f"{input_file}\t{NLUMO}\t{LUMO}\t{HOMO}\t{NHOMO}\n")
            print(f"{input_file}\t{NLUMO}\t{LUMO}\t{HOMO}\t{NHOMO}")
    return


def extract_value(line, keyword):
    return float(line.split(keyword)[-1].split()[0])


def PhaseCheck(tcal_path, Debug):
    print("\n**********\nPhase Checking...")
    chkKEY = "-12"
    CalCores = mkCalCoreList(tcal_path, chkKEY)
    if len(CalCores) == 0:
        print(f"\t>>> Any chk file is NOT found in the specified folder.\n"
              f"\t>>> The process for phase check was {Color.GREEN}skipped.{Color.RESET}")
    else:
        with open(f"{tcal_path}/PhaseCheck.txt", "w") as file:
            file.write("Name\tfor LUMO\tfor HOMO\n")
            for CalCore in CalCores:
                mkCubeFile(tcal_path, f"{CalCore}_m1.chk")
                mkCubeFile(tcal_path, f"{CalCore}_m2.chk")

                HomoChk = ComparePhase(tcal_path, CalCore, "HOMO")
                LumoChk = ComparePhase(tcal_path, CalCore, "LUMO")

                file.write(f"{CalCore}\t{LumoChk}\t{HomoChk}\n")
        if Debug:
            pass
        else:
            rmWildCards(f"{tcal_path}/*.chk")
            rmWildCards(f"{tcal_path}/*d-*.log")
            rmWildCards(f"{tcal_path}/*.gjf")
    return


def mkCalCoreList(tcal_path, chkKEY):
    FileList = glob.glob(f"{tcal_path}/*{chkKEY}*")
    CalCores = []
    for file in FileList:
        if "_m1.chk" in file or "_m2.chk" in file:
            CalCores.append(file[file.rfind("/") + 1:-7])
        else:
            pass
    CalCores = list(set(CalCores))
    return CalCores


def mkCubeFile(tcal_path, chkfile):
    def run_command(command, description):
        subprocess.run(command, cwd=tcal_path, timeout=1000)
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


def ComparePhase(tcal_path, CalCore, HomoLumo):
    """

    :param tcal_path:
    :param CalCore:
    :param HomoLumo:
    :return:
    """
    Val1 = ReadCube(tcal_path, f"{CalCore}_m1_{HomoLumo}.cub")
    Val2 = ReadCube(tcal_path, f"{CalCore}_m2_{HomoLumo}.cub")

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


def ReadCube(tcal_path, Cubefile):
    """

    :param tcal_path:
    :param Cubefile:
    :return:
    """
    with open(f"{tcal_path}/{Cubefile}", "r") as file:
        for _ in range(2):  # skip the first two lines (comment and number of atoms)
            next(file)

        NumE, Xgrid, Ygrid, Zgrid = (int(next(file).split()[0]) for _ in range(4))
        for _ in range(NumE):  # skip the number of electrons
            next(file)

        values = [float(x) for _ in range(Xgrid * Ygrid)
                  for line in file for x in line.split()]

    subprocess.run(["rm", Cubefile], cwd=tcal_path, timeout=10)
    return values


def Result_Data_set(MaterName, Nmol, Formated_Tilt, mol_pos, tcal_path, messages, HelpList):
    result_name = f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d"
    print(f"\n**********\n{Color.GREEN} Resulting Data Set for {result_name}...{Color.RESET}")

    result_path = f"{result_name}_results"
    os.makedirs(result_path, exist_ok=True)

    MinFileName = f"{result_name}_min.txt"
    TIFileName = f"{result_name}_TIs.txt"
    AllFileName = f"{result_name}_all.txt"
    PCFileName = "PhaseCheck.txt"

    Lacks = []

    if os.path.isfile(MinFileName) and os.path.isfile(f"{tcal_path}/{TIFileName}"):
        combineData(MaterName, Nmol, tcal_path, result_path, MinFileName, TIFileName,
                    PCFileName, Formated_Tilt, mol_pos, messages, HelpList)
    else:
        copy_file(MinFileName, f"{result_path}/{MinFileName}", Lacks, HelpList)
        print("\n")

    print("\tCopying Files required for recalculation...")
    copy_file(AllFileName, f"{result_path}/{AllFileName}", Lacks, HelpList)

    copy_file(f"{MaterName}.xyz", f"{result_path}/{MaterName}.xyz", Lacks, HelpList)

    condition_list = f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt"
    copy_file(condition_list, f"{result_path}/{condition_list}", Lacks, HelpList)

    copy_file("CalcSetting_HB.txt", f"{result_path}/CalcSetting_HB.txt", Lacks, HelpList)

    print("\n\tCopying min file...")
    copy_file(f"{result_name}_min.txt", f"{result_path}/{result_name}_min.txt", Lacks, HelpList)

    if "3mol" in Nmol:
        copy_file(f"{result_name}_mins.hist", f"{result_path}/{result_name}_mins.hist", Lacks, HelpList)

    if "2mol" in Nmol:
        struct_files = glob.glob(f"{tcal_path}/*.xyz")
    elif "3mol" in Nmol:
        struct_files = glob.glob(f"{tcal_path}/*.all")
    else:
        struct_files = []

    if not struct_files:
        print("\t>>> Any file for aggregation structure is not found in the specific folder.")
        Lacks.append("Structural xyz files")
    else:
        print("\n\tCopying structural xyz files...")
        for file in struct_files:
            name = os.path.basename(file).replace(".all", "")
            copy_file(file, f"{result_path}/{name}.xyz", Lacks, HelpList)
        print(f"\t>>> {Color.GREEN}Structural Files were copied into the {result_path} folder.{Color.RESET}")

    if Lacks:
        messages.append(f"\nSeveral result files were not found in the specific folder:")
        for lack in Lacks:
            messages.append(f"\t>>> {lack}")
    else:
        messages.append(f"\n\t{Color.GREEN}All files were copied successfully!!{Color.RESET}")

    help_check_exit(messages, HelpList)
    return


def combineData(MaterName, Nmol, tcal_path, result_path, MinFileName, TIFileName, PCFileName,
                Formated_Tilt, mol_pos, messages, HelpList):
    print(f"\tCombining Data...\n")
    with open(f"{tcal_path}/{TIFileName}", "r") as TIFile:
        TILines = TIFile.readlines()
        del TILines[0:2]

    with open(MinFileName, "r") as MinFile:
        MinLines = MinFile.readlines()
        del MinLines[0:2]

    if os.path.exists(f"{tcal_path}/{PCFileName}"):
        with open(f"{tcal_path}/{PCFileName}", "r") as PCFile:
            PCLines = PCFile.readlines()
            del PCLines[0]
    else:
        PCLines = []

    if len(TILines) == len(MinLines) and "2mol" in Nmol:
        pass
    elif len(TILines) == len(MinLines) * 3 and "3mol" in Nmol:
        pass
    else:
        messages.append(f"\t>>>{Color.RED}Error: The numbers of data lines in {MinFileName} "
                        f"and {TIFileName} DO NOT match.{Color.RESET}\n"
                        f"A file was NOT be changed.")
        HelpList.append(True)
    help_check_exit(messages, HelpList)

    CombLines = []
    CombLines12 = []
    CombLines23 = []
    CombLines31 = []
    print("Entry\tAngle\tDcol\tDtrv\tCpCE\tBSE\t**"
          "\tTI-NLUMO\tTI-LUMO\tTI-HOMO\tTI-NHOMO\t**\tPC (LUMO)\tTI-LUMO\tPC (HOMO)\tTI-HOMO")
    for MinLine in MinLines:
        MinData = MinLine.split()

        if "2mol" in Nmol:
            Entry = (f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d-{int(float(MinData[0]))}d"
                     f"-{int(round(float(MinData[1]) * 100, 5))}")
        elif "3mol" in Nmol:
            Entry = (f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_{int(float(MinData[0]))}d"
                     f"-{int(round(float(MinData[1]) * 100, 5))}"
                     f"-{int(round(float(MinData[2]) * 100, 5))}")
        else:
            messages.append(f"\t>>>{Color.RED}Error: UNEXPECTED (Nmol not specified){Color.RESET}")
            HelpList.append(True)
            help_check_exit(messages, HelpList)
            exit()

        for TILine in TILines:
            TIData = TILine.split()
            if (TIData[0] == Entry
                    or TIData[0] == f"{Entry}-12"
                    or TIData[0] == f"{Entry}-23"
                    or TIData[0] == f"{Entry}-31"):
                CombLine_temp = (f"{TIData[0]}\t{MinData[0]}\t{MinData[1]}\t{MinData[2]}\t{MinData[3]}\t{MinData[4]}"
                                 f"\t**\t{TIData[1]}\t{TIData[2]}\t{TIData[3]}\t{TIData[4]}")

                HomoChk = "yet"
                LumoChk = "yet"
                if len(PCLines) != 0:
                    for PCLine in PCLines:
                        PCData = PCLine.split()
                        if PCData[0].strip() == TIData[0].strip():
                            LumoChk = PCData[1]
                            HomoChk = PCData[2]
                else:
                    pass

                correctTI_LUMO = correctTI(LumoChk, TIData[2])
                correctTI_HOMO = correctTI(HomoChk, TIData[3])
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
    help_check_exit(messages, HelpList)
    CombLines.sort()
    CombLines12.sort()
    CombLines23.sort()
    CombLines31.sort()
    saveCombData(f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d",
                 f"{result_path}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min-TIs.txt", CombLines)
    saveCombData(f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d-12",
                 f"{result_path}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min-TIs-12.txt", CombLines12)
    saveCombData(f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d-23",
                 f"{result_path}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min-TIs-23.txt", CombLines23)
    saveCombData(f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d-31",
                 f"{result_path}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min-TIs-31.txt", CombLines31)
    print(f"\n"
          f"\t>>> {Color.GREEN}Combining Data: Succeeded!!{Color.RESET}\n"
          f"\t>>> {MinFileName} and {TIFileName} were combined into "
          f"{result_path}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min-TIs*.txt.\n")
    return


def correctTI(PhaseChk, TI):
    TI = float(TI)
    if PhaseChk == "Same":
        CorrectTI = TI
    elif PhaseChk == "Opposit":
        CorrectTI = -1 * TI
    else:
        CorrectTI = TI
    return CorrectTI


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


def copy_file(src, dest, lacks_list, HelpList):
    if os.path.isfile(src):
        execute(["cp", src, dest], False)
        print(f"\t>>> {src} copy to {dest}: {Color.GREEN}Complete{Color.RESET}")
    else:
        lacks_list.append(src)
        HelpList.append(True)
    return


if __name__ == "__main__":
    main()
