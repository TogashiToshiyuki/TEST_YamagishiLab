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
    AbnStop = ("\n"
               "*************** PROGRAM ABNORMALLY STOPPED ***************")
    HelpText = ("\n"
                "The required file may not exist.\n"
                "Please check.")
    CalcSet_BW_template = ("Edge Axis: [x, y, or z]\n"
                           "Faceon Axis: [x, y, or z]\n"
                           "Mol3 Other_Transition: A.AA\n"
                           "\n"
                           "[Comment]\n")
    InitialCondition_2mol_template = "8.0\n"
    InitialCondition_3mol_template = ("-3.0 3.9\n"
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


class CheckRequired(argparse.Action):
    """
    Custom argparse action to enforce the presence of a dependent argument.

    This class defines a custom action for the argparse module, which ensures
    that a specific argument (`--chk`) must be provided when another argument
    using this action is specified. If `--chk` is not present, the parser will
    raise an error and terminate the program.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        if not getattr(namespace, 'chk', False):
            parser.error(f"{option_string} requires --chk")
        setattr(namespace, self.dest, values)


def main():
    messages, HelpList, Debug = [], [], False

    # argument processing
    args, Debug, MaterName, calculation_tcal_Flag, Nmol, Flag_XYZ = getArgument(HelpList, messages)

    # File set check
    mol_pos = CheckFiles(HelpList, messages, Nmol, MaterName)

    # get Operator Name
    Operator = getOperator()

    # Create conditions for the first calculation
    RefLines, dev, dirpath, tcalpath, which = MakeFirstCondition(Nmol, MaterName, mol_pos)

    # Search for stable structures
    getTemporaryStructure(MaterName, Nmol, mol_pos, which, RefLines, dev, dirpath, Debug, Operator)


def printfList(List):
    """リストの要素を順に出力する。

    Args:
    List (list): 出力する要素を持つリスト。

    Returns:
    None
    """
    for content in List:
        print(content)
    return


def CheckPoint(HelpList):
    """
    Check if any element in the HelpList is True and, if so, print help
    text and terminate the program.

    This function iterates through the provided list, `HelpList`, and
    checks for the presence of `True`. If a `True` value is found, the
    function prints predefined help text and an abnormal stop message
    before terminating the program using `exit()`. If no `True` value
    is present, the function simply passes without further action.

    :param HelpList: A list of boolean values to check for `True`.
    :type HelpList: list[bool]
    :returns: None
    :raises SystemExit: If a `True` value is found in `HelpList`.
    """
    if True in HelpList:
        print(Stereotyped.HelpText)
        print(Stereotyped.AbnStop)
        exit()
    else:
        pass
    return


def getArgument(HelpList, messages):
    """
    Parse and validate command-line arguments for molecular calculations.

    This function handles the parsing of command-line arguments for a
    molecular calculation program. It ensures that required arguments
    are provided and validates the format of the molecular name file.
    The function also manages debug mode and checks for mutually
    exclusive options related to molecular calculations.

    :param HelpList: A list that stores error flags related to
                     argument validation.
    :type HelpList: list
    :param messages: A list that stores error or warning messages
                     encountered during argument validation.
    :type messages: list
    :returns: A tuple containing the parsed arguments (`args`),
              a debug flag (`Debug`), the molecular name (`MaterName`),
              a flag for tcal calculation (`calculation_tcal_Flag`),
              and the selected molecular calculation type (`Nmol`).
    :rtype: tuple
    :raises ArgumentError: If the required molecular name file is not
                           provided or if invalid combinations of
                           arguments are selected.
    """
    Debug, MaterName = False, ""
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

    # Check if at least one of --chk, --two_mol, or --three_mol is provided
    if not (args.chk or args.two_mol or args.three_mol):
        parser.error("One of --chk, --two_mol, or --three_mol must be provided.")

    if os.path.exists(args.MaterNameXYZ) and args.MaterNameXYZ.endswith(".xyz"):
        MaterName = args.MaterNameXYZ[:-4]
    else:
        parser.error("Incorrect file selected. [molecule name].xyz is required.")

    if args.debug:
        Debug = True
        print("*************** Caution!!! Debug Started!!! ***************\n")

    if args.two_mol and args.three_mol:
        messages.append("2mol and 3mol cannot be selected at the same time!!!")
        HelpList.append(True)
    if not args.two_mol and not args.three_mol:
        messages.append("2mol and 3mol cannot be selected at the same time.")
        HelpList.append(True)

    calculation_tcal_Flag = args.tcal

    Nmol = ""
    if args.two_mol:
        Nmol = "2mol"
    if args.three_mol:
        Nmol = "3mol"

    return args, Debug, MaterName, calculation_tcal_Flag, Nmol, args.chk


def CheckFiles(HelpList, messages, Nmol, MaterName):
    """
    Check the existence of required files and guide the user through
    necessary steps if any are missing.

    This function verifies the presence of critical files such as
    'CalcSetting_BW.txt', molecular structure files, and condition lists.
    It also assists the user in selecting a molecular structure if the
    calculation involves three molecules. If any required files are missing,
    default templates are used to create them, and the user is prompted
    to restart the program.

    :param HelpList: A list that tracks whether help is needed based on
                     the existence of required files.
    :type HelpList: list
    :param messages: A list of messages to be displayed to the user,
                     documenting the status of each file check.
    :type messages: list
    :param Nmol: A string representing the number of molecules involved
                 in the calculation, e.g., "2mol" or "3mol".
    :type Nmol: str
    :param MaterName: The base name of the material files to check,
                      such as 'MaterName.xyz'.
    :type MaterName: str
    :return: The molecular position identifier if '3mol' is selected,
             otherwise an empty string.
    :rtype: str
    :raises ValueError: If an invalid molecular position number is
                        entered by the user.
    """
    # Check CalcSetting_BW.txt
    print("\nChecking if 'CalcSetting_BW.txt' exists...")
    if os.path.exists("CalcSetting_BW.txt"):
        messages.append("\tCalcSetting_BW.txt: Found")
        HelpList.append(False)
    else:
        messages.append("\tCalcSetting_BW.txt: NOT Found\n"
                        "\tMake the correct CalcSetting_BW.txt and then restart the program.")
        HelpList.append(True)
        with open("CalcSetting_BW.txt", "w") as file:
            file.write(Stereotyped.CalcSet_BW_template)

    printfList(messages)
    CheckPoint(HelpList)
    messages.clear()
    HelpList.clear()

    if "3mol" in Nmol:
        while True:
            mol_pos_number = input("\n"
                                   "Calculation from 3 molecules have been selected.\n"
                                   "Please select a structure from 1~3..\n"
                                   "1: H-H\n"
                                   "2: T-T\n"
                                   "3: H-T\n"
                                   ">>> ")
            try:
                mol_pos_number = int(mol_pos_number)
            except (IndexError, ValueError):
                print("Please enter a number.")
                continue
            if mol_pos_number in {1, 2, 3}:
                break
            else:
                print("Error: Incorrect input.")
                continue
        mol_pos = f"p{mol_pos_number}"
    else:
        mol_pos = ""

    print("\nFile set check...")
    # MaterName.xyz
    if os.path.exists(f"./{MaterName}.xyz"):
        messages.append(f"\t{MaterName}.xyz: Found")
        HelpList.append(False)
    else:
        messages.append(f"\t{MaterName}.xyz: NOT Found")
        HelpList.append(True)

    # ConditionList or InitialCondition
    if os.path.exists(f"./ConditionList_{Nmol}{mol_pos}.txt"):
        messages.append(f"\tConditionList_{Nmol}{mol_pos}.txt: Found")
        HelpList.append(False)
    elif os.path.exists(f"InitialCondition_{Nmol}{mol_pos}.txt"):
        messages.append(f"\tInitialCondition_{Nmol}{mol_pos}.txt: Found")
    else:
        messages.append(f"\tBoth ConditionList_{Nmol}{mol_pos}.txt and InitialCondition_{Nmol}{mol_pos}.txt "
                        f"could not be Found.\n"
                        f"\tMake the correct InitialCondition_{Nmol}{mol_pos}.txt and then restart the program.")
        HelpList.append(True)
        with open(f"InitialCondition_{Nmol}{mol_pos}.txt", "w") as f:
            if "2mol" in Nmol:
                f.write(Stereotyped.InitialCondition_2mol_template)
            else:
                f.write(Stereotyped.InitialCondition_3mol_template)

    # MaterName_2mol_min
    if "3mol" in Nmol:
        if os.path.exists(f"./{MaterName}_2mol_min.txt"):
            messages.append(f"\t{MaterName}_2mol_min.txt: Found")
            HelpList.append(False)
        else:
            messages.append(f"\n"
                            f"\t{MaterName}_2mol_min.txt: NOT Found\n"
                            f"\tIt is possible that the calculation of 2mol has not been carried out.\n")
            HelpList.append(True)
    elif "2mol" in Nmol:
        HelpList.append(False)

    if True in HelpList:
        messages.append("Required file set DOES NOT exist in the correct directory.")
    else:
        messages.append("Required file set EXISTs in the correct directory.")

    printfList(messages)
    CheckPoint(HelpList)

    return mol_pos


def getOperator():
    Operator = input("\n"
                     "Enter your name.\n"
                     ">>> ")
    if Operator == "":
        Operator = "ONE"
    print("\n")
    return Operator


def MakeFirstCondition(Nmol, MaterName, mol_pos):
    which = ""
    RefLines = []
    if "2mol" in Nmol:
        which = "D_edge"
    elif "3mol" in Nmol:
        which = "D_faceon"
        RefLines = getRefLines(f"./{MaterName}_2mol_min.txt")
    if "3mol" in Nmol:
        dev = 0.2
    else:
        dev = 0.1

    mkConditionFile(Nmol, mol_pos, which, RefLines, dev)

    dirpath = f"./{MaterName}_{Nmol}"
    tcalpath = f"./{MaterName}_{Nmol}_tcal"
    os.makedirs(dirpath, exist_ok=True)
    print("")
    return RefLines, dev, dirpath, tcalpath, which


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


def mkConditionFile(Nmol, mol_pos, which, RefLines, dev):
    if os.path.exists(f"./ConditionList_{Nmol}{mol_pos}.txt"):
        print(f"ConditionList_{Nmol}{mol_pos}.txt: EXIST!!!")
        return
    else:
        print(f"ConditionList_{Nmol}{mol_pos}.txt: NOT EXIST!!!")
        with open(f"InitialCondition_{Nmol}{mol_pos}.txt", "r") as f:
            lines = f.readlines()
        NewConditions = []
        for line in lines:
            RefValues, Val, Other = [], 0, 0

            if "2mol" in Nmol:
                Val = float(line.strip())
                RefValues = ["na", "na"]
            elif "3mol" in Nmol:
                Other = float(line.strip().split()[0])
                Val = float(line.strip().split()[1])
                print(Val)
                print(RefLines[0])
                RefValues = [RefLines[0], RefLines[1]]

            NewCondition = mkNewCondition(Nmol, Other, Val - (2 * dev), which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Other, Val - (1 * dev), which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Other, Val, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Other, Val + (1 * dev), which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Other, Val + (2 * dev), which, RefValues)
            NewConditions.append(NewCondition)
        NewConditions = list(set(NewConditions))
        NewConditions.sort()
        with open(f"./ConditionList_{Nmol}{mol_pos}.txt", "w") as f:
            for Condition in NewConditions:
                f.write(Condition)
        print(f"New ConditionList_{Nmol}{mol_pos}.txt was written from InitialCondition_{Nmol}{mol_pos}.txt!!!")
        return


def getRefValues(Other, RefLines):
    RefValues = ''
    if len(RefLines) == 0:
        RefValues = ["na", "na"]
    else:
        for RefLine in RefLines:
            Contents = RefLine.strip().split('\t')
            try:
                RefOther = round(float(Contents[0]), 1)
                if RefOther == Other:
                    RefValues = [Contents[1], Contents[2]]
                else:
                    pass
            except (ValueError, IndexError):
                pass

    return RefValues


def mkNewCondition(Nmol, Other, Val, which, RefValues):
    Other = int(Other * 100)
    if "2mol" in Nmol:
        D_edge = int(round(Val * 100, 2))
        NewCondition = f"{D_edge}\n"
    elif "3mol" in Nmol:
        if which == "D_edge":
            D_edge = int(round(Val * 100, 2))
            D_faceon = int(float(RefValues[1]) * 100)
        elif which == "D_faceon":
            D_faceon = int(round(Val * 100, 2))
            D_edge = int(float(RefValues[0]) * 100)
        else:
            D_edge = D_faceon = 0
        print(f"D_edge: {D_edge}")
        print(f"D_faceon: {D_faceon}")
        NewCondition = f"{Other}_{D_edge}_{D_faceon}\n"
    else:
        NewCondition = ""
    print(f"\t\t{NewCondition.strip()}")
    return NewCondition


def getTemporaryStructure(MaterName, Nmol, mol_pos, which, RefLines, dev, dirpath, Debug, Operator):
    judge = False
    while not judge:
        with open(f"{dirpath}/G.sh", "w") as f:
            f.write(Stereotyped.Sh_txt)
        Conditions = getConditions(f"ConditionList_{Nmol}{mol_pos}.txt")
        qsubList = []
        for Condition in Conditions:
            if os.path.exists(f"./{dirpath}/{MaterName}_{Nmol}{mol_pos}_{Condition}.log"):
                pass
            else:
                qsub_temp = mkFiles(MaterName, Nmol, mol_pos, Condition, Operator, dirpath, False)
                qsubList.append(qsub_temp)
        print("\n**********\nJobs are submitting...")

        if len(qsubList) == 0:
            print("Any job was not submitted. Calculations with the conditions might be finished.")
        else:
            if "2mol" in Nmol:
                Wait_minutes = 1
            elif "3mol" in Nmol:
                Wait_minutes = 2
            else:
                Wait_minutes = 1

            for qsub in qsubList:
                qsub = qsub.split()
                subprocess.run(qsub, cwd=f"./{dirpath}")

            MyJobIDList = My_JobIDList(qsubList)
            MyJobIDList.sort()
            RunningList = Running_JobIDList()

            Flag, wait_job_count = check_jobs(RunningList, MyJobIDList)
            start_time = datetime.datetime.now()
            if Flag:
                formated_ST = start_time.strftime("%m/%d %H:%M:%S")
                if which == "Dcol":
                    term = "the stable distance in column direction"
                elif which == "Dtrv":
                    term = "the stable distance in transverse direction"
                else:
                    term = ""
                print(f"\n'{int(len(MyJobIDList))}' calculations for '{term}' was submitted!! at {formated_ST}")
                print(f"\n"
                      f"Wait until jobID {MyJobIDList[-1]}!!\n"
                      f"start ID: {RunningList[0]}\n"
                      f"\n")
            else:
                pass

            while Flag:
                Flag, wait_job_count = check_jobs(RunningList, MyJobIDList)
                formated_NOW, elapsed_time = getElapsedTime(start_time)
                Now_time = datetime.datetime.now()
                formated_end_time = (
                    (Now_time + datetime.timedelta(minutes=(wait_job_count * Wait_minutes))).strftime("%m/%d %H:%M:%S"))
                sys.stdout.write(
                    "\033[1F\033[G%s" %
                    f"\t{formated_NOW} ({elapsed_time} min. passed): Job {RunningList[0]} is in progress.            \n"
                    f"\tNext Check >>> {Wait_minutes * wait_job_count} minute later! forecast: {formated_end_time}    ")
                sys.stdout.flush()
                time.sleep(Wait_minutes * wait_job_count * 60)
                RunningList = Running_JobIDList()
                Flag, wait_job_count = check_jobs(RunningList, MyJobIDList)
            else:
                print(f"\n\nCalculation cycles for {which} until JobID {MyJobIDList[-1]} were finished.")
                rmWildCards(f"{dirpath}/*.sh*")
                if Debug:
                    pass
                else:
                    rmWildCards(f"{dirpath}/*.chk")
        print("\n**********\nReading Data...\n")
        X_readEnergies(dirpath, MaterName, Nmol, mol_pos)
        judge = True

    if "2mol" in Nmol:
        print("Calculations for 2mol were successfully finished.")
    else:
        pass
    return


def getConditions(FileName):
    Conditions = []
    with open(FileName, "r") as f:
        lines = f.readlines()
    for line in lines:
        Conditions.append(line.strip())
    return Conditions


def mkFiles(MaterName, Nmol, mol_pos, Condition, Operator, dirpath, Flag_XYZ):
    File_Name = f"{MaterName}_{Nmol}{mol_pos}_{Condition}"
    CHK_Name = File_Name + ".chk"
    GJF_Name = File_Name + ".gjf"
    SH_Name = f"G-{Operator}_{Condition}.sh"
    Flag_XYZ = True

    ConditionList = Condition.strip().split("_")
    Edge_Axis, Face_Axis, Other_Axis, Matrix_Mol3, rotate = Axis_Setting_BW()
    if "3mol" in Nmol:
        Matrix_Mol3 = {
                          "x": [int(ConditionList[0]), 0, 0],
                          "y": [0, int(ConditionList[0]), 0],
                          "z": [0, 0, int(ConditionList[0])]
                      }[Face_Axis] + Matrix_Mol3
        Transitions = mkTransition(Edge_Axis, Face_Axis, ConditionList[1], ConditionList[2], Matrix_Mol3)
    else:
        Transitions = mkTransition(Edge_Axis, Face_Axis, ConditionList[0], 0, Matrix_Mol3)

    Angles = {
        "": {
            "x": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            "y": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            "z": [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        },
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
    }.get(mol_pos).get(Other_Axis)

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
            print(f"{dirpath}/{File_Name}.xyz have been created.")
    elif "3mol" in Nmol:
        write_gjf_file(f"{dirpath}/{File_Name}.gjf",
                       Headers, Element, Mol1_pos, Mol2_pos, Mol3_pos)
        if Flag_XYZ:
            write_xyz_file(f"{dirpath}/{File_Name}.xyz",
                           Element, Mol1_pos, Mol2_pos, Mol3_pos)
            print(f"{dirpath}/{File_Name}.xyz have been created.")

    with open(f"{dirpath}/G.sh", "r") as orgSH:
        lines = orgSH.readlines()
    lines[12] = f"g16 {GJF_Name}\n"
    with open(f"{dirpath}/{SH_Name}", "w") as newSH:
        for line in lines:
            newSH.write(line)
    qsub_temp = f"qsub {SH_Name}"
    return qsub_temp


def Axis_Setting_BW():
    """
    CalcSetting_BW.txtファイルを読み取り、軸設定を行います。

    AxSetting.txtファイルを読み取り、各軸の設定を行います。
    ファイルが存在しない場合は新しいファイルを作成し、初期設定を記述します。
    Edge Axis, Faceon Axis, Mol3 Edge Transition, Mol3 Other Transitionの設定を
    解析し、対応する軸を計算して返します。

    パラメータ:
        なし

    戻り値:
        edge axis (str): Edge Axisの設定値
        faceon axis (str): Faceon Axisの設定値
        Matrix_Mol3 (str): Mol3 Edge TransitionとMol3 Other Transitionを合成した値

    発生する例外:
        FileNotFoundError: AxSetting.txtファイルが存在しない場合
    """
    try:
        with open("CalcSetting_BW.txt", "r") as File:
            lines = File.readlines()
    except FileNotFoundError:
        print("\nThe AxSetting.txt is not found in the current directory.\n"
              "Make the correct AxSetting.txt and then restart the program.\n")
        with open("CalcSetting_BW.txt", "w") as file:
            file.write("Edge Axis: [x, y, or z]\n"
                       "Faceon Axis: [x, y, or z]\n"
                       "Mol3 Other_Transition: A.AA\n"
                       "\n"
                       "[Comment]")
        print(Stereotyped.AbnStop)
        exit()

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

    Mol3_Other = mkDirection(params["mol3 other_transition"], params["other axis"], "Mol3 Other Transition")

    Matrix_Mol3 = np.array(Mol3_Other)

    return params["edge axis"], params["faceon axis"], params["other axis"], Matrix_Mol3, rotate


def mkDirection(axis_direction, input_axis, axis_name):
    try:
        axis_direction = float(axis_direction)
    except ValueError:
        print(f"\nThe axis direction for {axis_name} is INVALID: {axis_direction}")
        print(Stereotyped.AbnStop)
        exit()

    directions = {"x": [axis_direction, 0.0, 0.0], "y": [0.0, axis_direction, 0.0], "z": [0.0, 0.0, axis_direction]}
    if input_axis not in directions:
        print(f"\nThe axis of {axis_name} is INVALID.")
        print(Stereotyped.AbnStop)
        exit()
    return np.array(directions[input_axis])


def mkTransition(Axis_Edge, Axis_Faceon, Direction_Edge, Direction_Faceon, Matrix_Mol3):
    if Direction_Edge == Direction_Faceon:
        print("\nThe axis of Edge direction is same with the axis of Faceon direction.")
        print(Stereotyped.AbnStop)
        exit()

    Edge_transl = mkDirection(Direction_Edge, Axis_Edge, "Edge Direction")
    Faceon_transl = mkDirection(Direction_Faceon, Axis_Faceon, "Faceon Direction")

    return [Edge_transl, Faceon_transl, Matrix_Mol3]


def mkAtomList(MaterName, Mol1_Angles, Mol2_Angles, Mol3_Angles, rotate, Translations):
    """
    原子リストを生成する関数。

    指定されたファイルから原子の座標を読み込み、3つの分子構造に対して
    回転および平行移動を行い、原子リストを生成します。

    Parameters:
    MaterName (str): 分子の名前（拡張子なし）。
    Mol1_Angles (tuple): 分子1に対する回転角度（ラジアン）。
    Mol2_Angles (tuple): 分子2に対する回転角度（ラジアン）。
    Mol3_Angles (tuple): 分子3に対する回転角度（ラジアン）。
    rotate (str): 回転順を表す(xyz,zxyなど)
    Translations (list): 分子2および分子3に対する平行移動ベクトル。

    Returns:
    tuple: 原子の元素記号リスト、分子1、分子2、分子3の各原子位置リスト、
           読み込んだ原子の数。

    Raises:
    FileNotFoundError: 指定されたファイルが存在しない場合に発生します。
    ValueError: ファイルのフォーマットが正しくない場合に発生します。
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
        atm_m2 = Rotate(Position, *Mol2_Angles, rotate) + (Translations[0]/100)
        Mol2.append(atm_m2)
        atm_m3 = Rotate(Position, *Mol3_Angles, rotate) + Translations[0] / 200 + Translations[1]/100 + Translations[2]
        Mol3.append(atm_m3)
    return Elist, Mol1, Mol2, Mol3, NinMol


def Rotate(Current, Tx, Ty, Tz, rotation):
    """3次元座標を指定された軸と角度で回転させる。

    Args:
        Current (np.array): 回転させたい3次元座標 (形状: (3,)).
        Tx (float): x軸周りの回転角度 (単位: 度).
        Ty (float): y軸周りの回転角度 (単位: 度).
        Tz (float): z軸周りの回転角度 (単位: 度).
        rotation (str): 回転軸の順序を表す文字列 ('x', 'y', 'z' の組み合わせ).

    Returns:
        np.array: 回転後の3次元座標 (形状: (3,)).

    Raises:
        ValueError: rotation に無効な軸が含まれている場合.
        TypeError: Current が np.array でない場合.
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
    """Gaussian input file (gjf)を作成する。

    複数の座標セット（最大3つ）とダミー原子座標を扱うことができる。

    Args:
    filename (str): 出力ファイル名 (.gjf 拡張子)
    headers (list of str): ヘッダー行のリスト
    elements (list of str): 元素記号のリスト
    positions1 (list of float): 1つ目の座標セット ([[x1, y1, z1], [x2, y2, z2], ...])
    positions2 (list of float, optional): 2つ目の座標セット
    positions3 (list of float, optional): 3つ目の座標セット
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
    return f"{coord[0]: 15.10f}     {coord[1]: 15.10f}     {coord[2]: 15.10f}"


def My_JobIDList(qsubList):
    """qsubした自分のジョブのIDリストを作成する。

    qstatコマンドを実行し、ジョブIDリストを作成する。ジョブが実行されていない場合は、[0]を返す。

    Args:
        qsubList: qsubジョブのリスト。

    Returns:
        qsubジョブのIDリスト。
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
    現在実行中のジョブIDのリストを取得する関数。

    この関数はqstatコマンドを実行し、実行中のジョブIDを取得してリストに
    格納し、ソートしたリストを返します。

    戻り値:
        list: 実行中のジョブIDを含むソートされたリスト。

    例外:
        なし。ただし、qstatコマンドが失敗した場合や、出力フォーマットが
        期待通りでない場合、正しい結果が得られない可能性があります。
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
    """自分のジョブが現在のジョブリストに含まれているか、終了時間を推定する。

    Args:
        current_jobs (list): 現在実行中のジョブのIDリスト。
        my_jobs (list): 自分のジョブのIDリスト。

    Returns:
        tuple:
            - bool: 自分のジョブが含まれているかどうか。
            - int: 自分のジョブが終了するまでの推定時間（ジョブ数）。含まれていない場合は0。
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
    """開始時刻からの経過時間（分）と現在時刻を計算する。

    Args:
        start_time (datetime.datetime): 経過時間を計測する開始時刻。

    Returns:
        tuple: 現在時刻（"MM/DD HH:MM"形式）と開始時刻からの経過時間（分、小数点以下1桁）のタプル。
    """
    now = datetime.datetime.now()
    elapsed_time = now - start_time
    elapsed_time = round(elapsed_time.total_seconds() / 60, 1)
    formated_NOW = now.strftime("%m/%d %H:%M")
    return formated_NOW, elapsed_time


def rmWildCards(wildcard):
    """
    指定されたワイルドカードに一致するファイルを削除する。

    指定されたワイルドカードに一致するファイルをすべて検索し、
    見つかったファイルを削除する。この関数は、システムの
    コマンド `rm` を使用してファイルを削除する。

    パラメータ:
    wildcard (str): ワイルドカード文字列。

    戻り値:
    なし

    発生する例外:
    subprocess.TimeoutExpired: 削除処理がタイムアウトした場合。
    """
    lines = glob.glob(wildcard)
    for line in lines:
        subprocess.run(["rm", line], timeout=10)
    return


def X_readEnergies(dir_path, MaterName, Nmol, mol_pos):
    LogList = glob.glob(f"{dir_path}/{MaterName}_{Nmol}{mol_pos}_*.log")
    LogList.sort()
    List = []
    for Log in LogList:
        FileName, VEdge, VFaceon, VOther = getVALfromLogName(Nmol, Log)
        if "2mol" in Nmol:
            List.append(VEdge)
        elif "3mol" in Nmol:
            List.append(VOther)
    List = list(set(List))
    List.sort()

    with open(f"./{MaterName}_{Nmol}{mol_pos}_all.txt", "w") as AllData:
        with open(f"./{MaterName}_{Nmol}{mol_pos}_min.txt", "w") as MinData:
            header = ("Distance in Edge direction (Å)\tDistance in Faceone direction (Å)"
                      "\tDistance in Other direction (Å)\tCounterpoise corrected energy (A.U)\tBSSE energy (A.U)")
            AllData.write(f"*****  {MaterName}_{Nmol}{mol_pos} All Results *****\n")
            MinData.write(f"***** {MaterName}_{Nmol}{mol_pos} Minimum Energy at each Angle *****\n")
            AllData.write(f" \t{header}\n")
            MinData.write(f"{header}\n")
            MinimunL = []

            print(f"\t{header}")
            for Lists in List:
                CPE_Dict = {}
                VAL_Dict = {}
                Keys = []
                AllData.write("******\t******\t******\t******\t******\t******\n")
                for Log in LogList:
                    Condition = getCondition_fromName(Log)
                    if "2mol" in Nmol:
                        print(f"Log: {Log}")
                        print(f"Lists: {str(int(Lists*100))}")
                        if str(int(Lists*100)) in Log:
                            print("aaaaaaaaaaaaaaaaaaa")
                            with open(Log, "r") as file:
                                data = file.read()
                            if "Normal termination" in data:
                                FileName, VEdge, VFaceon, VOther = getVALfromLogName(Nmol, Log)
                                CPE, BSE = getEnergy(data)
                                CPE_Dict[FileName] = CPE
                                VAL_Dict[FileName] = VEdge
                                print(f"CPE: {CPE}")
                                print(f"V_Edge: {VEdge}")
                                Keys.append(FileName)
                            else:
                                print(f"\n{Log} was not terminated normally.\n")
                                os.remove(Log)
                    elif "3mol" in Nmol:
                        Condition = getCondition_fromName(Log)
                        Other = round(float(Condition.split("_")[2]), 5)
                        if Other == Lists:
                            with open(Log, "r") as file:
                                data = file.read()
                            if "Normal termination" in data:
                                FileName, VEdge, VFaceon, VOther = getVALfromLogName(Nmol, Log)
                                CPE, BSE = getEnergy(data)
                                CPE_Dict[FileName] = CPE
                                VAL_Dict[FileName] = VEdge
                                Keys.append(FileName)
                minkey = min(CPE_Dict, key=CPE_Dict.get)
                MinimunL.append(minkey)
                print(f"Val_Dict: {VAL_Dict}")
                print(f"minkey: {minkey}")
                MinData.write(str(VAL_Dict.get(minkey)))
                Keys.sort()
                print(Keys)
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


def X_getValfromLogName(Nmol, Log):
    """ログファイル名から角度、傾き、分子間距離などのパラメータを取得する。

    Args:
        Nmol (str): "2mol" または "3mol"。分子の数を表す。
        Log (str): ログファイルのパス。

    Returns:
        tuple: (FileName, Vdeg, Vdcol, Vdtrv, VTilt)
            - FileName (str): ファイル名。
            - Vdeg (float): 角度。
            - Vdcol (float): 傾き。
            - Vdtrv (float): 分子間距離。
            - VTilt (str): 傾きの方向。
    """
    try:
        FileName = Log.split("/")[-1].split(".")[0]
        Condition = FileName.split("_")[-1]

        if "2mol" in Nmol:
            VEdge = float(Condition.split("_")[2])
            VFaceon = "na"
            VOther = "na"

        elif "3mol" in Nmol:
            VEdge = float(Condition.split("_")[3])
            VFaceon = float(Condition.split("_")[4])
            VOther = float(Condition.split("_")[2])
        else:
            print("\nUNEXPECTED ERROR happens in the function of getVALfromLogName!!")
            exit()
        VEdge = round(float(VEdge), 5)
        try:
            VFaceon = round(float(VFaceon), 5)
        except ValueError:
            VFaceon = "na"
        try:
            VOther = round(float(VOther), 5)
        except ValueError:
            VOther = "na"
    except ValueError:
        print(f"\nUNEXPECTED ERROR happens in the function of getVALfromLogName!!\n")
        exit()
    return FileName, VEdge, VFaceon, VOther


def readEnergies(dir_path, MaterName, Nmol, Formated_Tilt, mol_pos):
    """
    Reads log files to extract and summarize molecular energy data.

    This function searches for log files within the specified directory,
    filters them based on naming conventions, and extracts energy data for
    each molecular orientation. The results are written to two output files:
    one containing all energy data and another containing only the minimum
    energy values for each angle.

    :param dir_path: Path to the directory containing log files.
    :type dir_path: The
    :param MaterName: Name of the material to process.
    :type MaterName: Str
    :param Nmol: Number of molecules in the simulation.
    :type Nmol: Int
    :param Formated_Tilt: Formatted tilt angle.
    :type Formated_Tilt: Str
    :param mol_pos: Molecular position identifier.
    :type mol_pos: Str
    :return: None :rtype:
    :raises FileNotFoundError: If any log file cannot be found.
    :raises IOError: If there is an issue reading a log file or writing output.
    :raises ValueError: If a log file does not contain expected data.
    """
    LogList = glob.glob(f"{dir_path}/{MaterName}_{Nmol}{mol_pos}_*.log")
    LogList.sort()
    DegList = []
    for Log in LogList:
        FileName, deg, Vdcol, Vdtrv = getVALfromLogName(Nmol, Log)
        DegList.append(deg)
    DegList = list(set(DegList))
    DegList.sort()

    with open(f"./{MaterName}_{Nmol}{mol_pos}_all.txt", "w") as AllData:
        with open(f"./{MaterName}_{Nmol}{mol_pos}_min.txt", "w") as MinData:
            header = ("Angle \tDistance in column direction (Å)\tDistance in transverse direction (Å)"
                      "\tCounterpoise corrected energy (A.U)\tBSSE energy (A.U)")
            AllData.write(f"*****  {MaterName}_{Nmol}{mol_pos} All Results *****\n")
            MinData.write(f"***** {MaterName}_{Nmol}{mol_pos} Minimum Energy at each Angle *****\n")
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
                            FileName, Vdeg, Vdcol, Vdtrv = getVALfromLogName(Nmol, Log)
                            CPE, BSE = getEnergy(data)
                            CPE_Dict[FileName] = CPE
                            VAL_Dict[FileName] = f"{Vdeg}\t{Vdcol}\t{Vdtrv}\t{CPE}\t{BSE}\n"
                            Keys.append(FileName)
                        else:
                            print(f"\t{Log} was NOT normally terminated. It was removed.")
                            List = []
                            try:
                                with open(f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_TerminatedLog.txt", "r") as f:
                                    Lines = f.readlines()
                                    for Line in Lines:
                                        List.append(Line)
                            except FileNotFoundError:
                                pass
                            List.append(Log)
                            List.sort()

                            with open(f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_TerminatedLog.txt", "w") as f:
                                for line in List:
                                    f.write(line)
                            os.remove(Log)

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


def getVALfromLogName(Nmol, Log):
    """ログファイル名と分子数から、角度、傾き、分子間距離などのパラメータを取得する。

    Args:
        Nmol (str): "2mol" または "3mol"。分子の数を表す。
        Log (str): ログファイルのパス。

    Returns:
        tuple: (FileName, Vdeg, Vdcol, Vdtrv, VTilt)
            - FileName (str): ログファイル名（拡張子なし）。
            - Vdeg (float): 角度（単位: 度）。
            - Vdcol (float): 分子間距離 dcol（単位: オングストローム）。
            - Vdtrv (float or str): 分子間距離 dtrv（単位: オングストローム）。"2mol" の場合は "na"。
            - VTilt (int): 傾き（単位: 度）。

    Raises:
        SystemExit: 予期しないエラーが発生した場合（`Nmol` が "2mol" でも "3mol" でもない場合）。

    Note:
        この関数は、ログファイル名に含まれる条件情報からパラメータを抽出します。
        ログファイル名の形式が特定の形式であることを前提としています。
        予期しない形式のログファイル名が入力された場合、エラーが発生します。
    """
    Condition = getCondition_fromName(Log)
    FileName = Log[Log.rfind("/") + 1:-4]

    Conditions = Condition.split('-')
    print(Condition)

    if "2mol" in Nmol and Condition.count("-") == 0:
        dEdge = Conditions[0].replace('d', '')
        dFaceon = "na"
        dOther = "na"
    elif "3mol" in Nmol and Condition.count("-") == 2:
        dEdge = Conditions[0].replace('d', '')
        dFaceon = Conditions[1]
        dOther = Conditions[2]
    else:
        print("\nUNEXPECTED ERROR happens in the function of getVALfromLogName!!")
        exit()

    VEdge = round(float(dEdge) / 100, 5)
    try:
        VFaceon = round(float(dFaceon) / 100, 5)
    except ValueError:
        VFaceon = dFaceon
    try:
        VOther = round(float(dOther) / 100, 5)
    except ValueError:
        VOther = dOther
    return FileName, VEdge, VFaceon, VOther


def getCondition_fromName(Name):
    Condition = Name[Name.rfind("_") + 1:Name.rfind(".")]
    return Condition


def getEnergy(data):
    """文字列からCounterpoise補正エネルギーとBSSEエネルギーを抽出する。

    Gaussianの出力ファイルなどから、Counterpoise補正エネルギーとBSSEエネルギーを
    抽出します。抽出された値は、それぞれ浮動小数点数に変換されます。

    Args:
        data (str): Counterpoise補正エネルギーとBSSEエネルギーを含む文字列。

    Returns:
        tuple: Counterpoise補正エネルギー(float)とBSSEエネルギー(float)のタプル。

    Raises:
        ValueError: dataが適切な形式でない場合。
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


def mkNewConditionList(MaterName, Nmol, which, dev, RefLines, mol_pos):
    with open(f"./{MaterName}_{Nmol}{mol_pos}_all.txt", "r") as All:
        AllLines = All.readlines()

    SetList = []
    for AllLine in AllLines:
        Contents = AllLine.strip().split("\t")


if __name__ == "__main__":
    main()
