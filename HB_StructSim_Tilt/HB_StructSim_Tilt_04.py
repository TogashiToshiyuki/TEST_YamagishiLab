#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import math
import datetime
import subprocess
import time
import glob
import random
import functools
import argparse

"""
HB_StructSim_Tilt
Created by: Toshiyuki Togashi (2024/7/11)

This program assumes that Gaussian 16 and AGE are installed on the workstation.
It performs the following five processes as required.
1. Create files: ConditionList, gjf, shell scripts
2. Perform energy integral and transfer integral calculations
3. Convert the gjf file with the optimised structure into xyz file
4. Check the phase of the molecular orbitals in the transfer integral calculation
5. Collect the files resulting from the calculations in the result folder.

Program overview.
The columnar and transverse distances are optimised alternately on the basis of single point energies.
If the obtained structure is the same as the structure obtained in the previous optimization cycle,
the obtained structure is considered the optimised structure.
The transfer integral is calculated on the basis of the optimised structure after the structural simulation.
When all calculations have been completed, the result and condition files are collected in the results' folder.
If condition files are used, the entire calculation process can be reproduced.

<Update information>
01: Development version X5 is terminated and moved to official version 01.
    Official version 01 is released as a direct copy of X5.
02: Integration of HB_StructCheck_XYZtoGJF_08.py
03: Add a programme to run with A built in.
    Algorithm for stable structure search changed from 0.2 to 0.1 to 0.05 Å.
    Other changes to argparse behaviour.

<Confirmation of operation>
Confirmed up to 2mol and 3mol structure optimisation and calculation in tcal.
NITT: Work Sation@10.24.27.56 (2024/8/7)

"""

printf = functools.partial(print, flush=True)


class Stereotyped:
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
                   "*                          created by Toshiyuki Togashi 2024/07/11      *\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n")
    AbnStop = ("\n"
               "*************** PROGRAM ABNORMALLY STOPPED ***************")
    AbnormalEnd = "\n************* Programme DID NOT terminate successfully. *************\n"
    HelpText = ("\n"
                "The required file may not exist.\n"
                "Please check.")
    # memo: NITT WS用のヘッダー、テキストのため、他で使用する場合は調整が必要
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


class Constant:
    # def mkNewConditionListの中でNew Conditionを一回でいくつ作るか
    # 増やすと状況によっては計算が早くなる
    Cn = 1
    d_another = 0


class Color:
    """
    文字色や背景色を定義するクラス
    """
    BLACK = '\033[30m'  # (文字)黒
    RED = '\033[31m'  # (文字)赤
    GREEN = '\033[32m'  # (文字)緑
    YELLOW = '\033[33m'  # (文字)黄
    BLUE = '\033[34m'  # (文字)青
    MAGENTA = '\033[35m'  # (文字)マゼンタ
    CYAN = '\033[36m'  # (文字)シアン
    LIME = '\033[92m'  # (文字)ライム
    WHITE = '\033[37m'  # (文字)白
    COLOR_DEFAULT = '\033[39m'  # 文字色をデフォルトに戻す
    BOLD = '\033[1m'  # 太字
    UNDERLINE = '\033[4m'  # 下線
    INVISIBLE = '\033[08m'  # 不可視
    REVERSE = '\033[07m'  # 文字色と背景色を反転
    BG_BLACK = '\033[40m'  # (背景)黒
    BG_RED = '\033[41m'  # (背景)赤
    BG_GREEN = '\033[42m'  # (背景)緑
    BG_YELLOW = '\033[43m'  # (背景)黄
    BG_BLUE = '\033[44m'  # (背景)青
    BG_MAGENTA = '\033[45m'  # (背景)マゼンタ
    BG_CYAN = '\033[46m'  # (背景)シアン
    BG_WHITE = '\033[47m'  # (背景)白
    BG_DEFAULT = '\033[49m'  # 背景色をデフォルトに戻す
    RESET = '\033[0m'  # 全てリセット


class CheckRequired(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if not getattr(namespace, 'chk', False):
            parser.error(f"{option_string} requires --chk")
        setattr(namespace, self.dest, values)


def main():
    """
    Main function for executing the molecular calculation workflow.

    This function manages the entire molecular calculation process, from
    initializing the structure based on a given .xyz file to performing
    structural optimization, electronic state calculations, and transfer
    integral calculations. The results are saved in a specified directory
    for further analysis.

    Detailed Description:
    The function first validates the presence and format of the molecular
    structure file. If the file is valid, it proceeds to structural checks
    and may generate additional files depending on the user's input. The
    user can choose between automatic or manual modes for generating the
    necessary calculation files. The function also handles the validation
    of input arguments and manages errors related to file handling and
    external program execution.

    Args:
        N/A: Parameters are obtained from command-line arguments.

    Returns:
        None: Results are output to standard output and written to files.

    Raises:
        IndexError: If insufficient command-line arguments are provided.
        ValueError: If an invalid file format (non-.xyz) is specified.
        FileNotFoundError: If required files are missing.
        Subprocess.CalledProcessError: If an external program fails to run.
        Exception: Other potential exceptions during the calculation process.

    Dependencies:
        This function depends on the following external modules and custom
        functions:
        - sys
        - os
        - argparse
        - Stereotyped (custom module)
        - printfList (custom function)
        - CheckPoint (custom function)
        - getRefLines (custom function)
        - mkConditionFile (custom function)
        - getTemporaryStructure (custom function)
        - CompareStructures (custom function)
        - getMinConditions (custom function)
        - mkXYZfile (custom function)
        - glob
        - subprocess
        - datetime
        - Running_JobIDList (custom function)
        - check_jobs (custom function)
        - getElapsedTime (custom function)
        - time
        - rmWildCards (custom function)
        - rename (Unix command)
        - readlog (custom function)
        - mkCalCoreList (custom function)
        - mkCubeFile (custom function)
        - ComparePhase (custom function)
        - combineData (custom function)
        - copy_file (custom function)
    """
    printf(Stereotyped.ProgramAbst)
    messages, HelpList, Debug, MaterName = [], [], False, ""

    args, MaterName = arg_parser(messages, HelpList)

    if args.debug:
        Debug = True
        messages.append(f"{Color.RED}*************** Caution!!! Debug Started!!! ***************{Color.RESET}\n")
        messages.append(f"In debug mode, the following files are not deleted\n"
                        f"\t - sh files\n"
                        f"\t - chk files\n"
                        f"\t - log files\n"
                        f"\t - gjf files\n")

    if args.chk:
        mkCheckFile(args, Debug, messages, HelpList, MaterName)

    if args.two_mol and args.three_mol:
        messages.append(f"{Color.RED}2mol and 3mol cannot be selected at the same time.{Color.RESET}")
        HelpList.append(True)
    if not args.two_mol and not args.three_mol:
        messages.append(f"{Color.RED}2mol and 3mol cannot be selected at the same time.{Color.RESET}")
        HelpList.append(True)

    calculation_tcal_Flag = args.tcal

    Nmol = ""
    if args.two_mol:
        Nmol = "2mol"
    if args.three_mol:
        Nmol = "3mol"

    printf("Checking if 'CalcSetting_HB.txt' exists.")
    if os.path.exists("./CalcSetting_HB.txt"):
        messages.append("\t>>> CalcSetting_HB.txt: Found")
        HelpList.append(False)
    else:
        messages.append(f"\t>>> {Color.RED}.txt: NOT Found{Color.RESET}\n"
                        "\t>>> Make the correct CalcSetting_HB.txt and then restart the program.")
        HelpList.append(True)
        with open("CalcSetting_HB.txt", "w") as file:
            file.write("Column Direction: [x, y, or z]\n"
                       "Transverse Direction: [x, y, or z]\n"
                       "Tilt Axis: [x, y, or z]\n"
                       "Tilt Angle: [Number]\n"
                       "\n"
                       "[Comment]\n")

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
            printf("Please enter a number.")
            Tilt = input("\nPlease enter the Tilt angle at degree."
                         "\n>>> ")
            continue
    printf(f"\n{Tilt}-degree has been entered as Tilt")
    Formated_Tilt = int(Tilt * 10)  # For File Name
    printf(f"\n{Color.RED}"
           f"*************** Caution!!! The Tilt angle multiplied by 10 is used in the file name. ***************"
           f"{Color.RESET}"
           f"\n*************** Like: {MaterName}_{Nmol}p1_{Color.RED}t{Formated_Tilt}d{Color.RESET}_10d-500-800.gjf\n")

    help_check_exit(messages, HelpList)

    if "3mol" in Nmol:
        while True:
            mol_pos_number = input(f"\n"
                                   f"Calculation from 3 molecules have been selected.\n"  # 3molでの計算が選ばれました
                                   f"{Color.GREEN}Please select a structure from 1~3..{Color.RESET}\n"
                                   f"\t1: H-H\n"
                                   f"\t2: T-T\n"
                                   f"\t3: H-T\n"
                                   f"\t>>> ")
            try:
                mol_pos_number = int(mol_pos_number)
            except (IndexError, ValueError):
                printf(f"{Color.RED} enter a number.{Color.RESET}")
                continue
            if mol_pos_number in {1, 2, 3}:
                break
            else:
                printf(f"{Color.RED}Error: Incorrect input.{Color.RESET}")
                continue
        mol_pos = f"p{mol_pos_number}"
    else:
        mol_pos = ""
        calculation_tcal_Flag = False

    printf(f"\n\n{Color.GREEN}File set check...{Color.RESET}")
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
            file.write("10 3.9\n"
                       "20 4.2\n"
                       "30 5.0\n")
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

    Operator = input(f"\n{Color.GREEN}Enter your name. >>> {Color.RESET}")
    if Operator == "":
        Operator = "ONE"
    else:
        pass
    printf("\n")

    which = ""
    RefLines = []
    if Nmol == "2mol":
        which = "Dcol"
    elif Nmol == "3mol":
        which = "Dtrv"
        RefLines = getRefLines(f"./{MaterName}_2mol_t{Formated_Tilt}d_min.txt")

    if "3mol" in Nmol:
        dev = 0.2
    else:
        dev = 0.1
    mkConditionFile(Nmol, mol_pos, which, RefLines, Formated_Tilt, dev)
    dirpath = f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d"
    tcalpath = f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_tcal"
    os.makedirs(dirpath, exist_ok=True)
    printf("\n")

    getTemporaryStructure(MaterName, Nmol, mol_pos, Tilt, which, RefLines, dev, dirpath, Debug, Operator, Formated_Tilt)
    if "2mol" in Nmol:
        printf(f"{Color.GREEN}Calculations for 2mol were successfully finished.{Color.RESET}")
    else:
        pass

    temp_Structures = []
    if "3mol" in Nmol:
        os.makedirs(tcalpath, exist_ok=True)

    if "3mol" in Nmol and not os.path.exists(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_mins.hist"):
        MostStable = False
        while not MostStable:
            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            temp_Structures.append(RefLines)
            mkCycleConditions(RefLines, "Dcol", 0.2, Nmol, Formated_Tilt, mol_pos)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Tilt, "Dcol", RefLines, 0.2, dirpath,
                                  Debug, Operator, Formated_Tilt)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            mkCycleConditions(RefLines, "Dtrv", 0.2, Nmol, Formated_Tilt, mol_pos)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Tilt, "Dtrv", RefLines, 0.2, dirpath,
                                  Debug, Operator, Formated_Tilt)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            MostStable = CompareStructures(RefLines, temp_Structures[-1])

        printf(f"\n{Color.GREEN}**********\nTransition in 0.1-Å increments.\n{Color.RESET}")
        MostStable = False
        while not MostStable:
            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            temp_Structures.append(RefLines)
            mkCycleConditions(RefLines, "Dcol", 0.1, Nmol, Formated_Tilt, mol_pos)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Tilt, "Dcol", RefLines, 0.1, dirpath,
                                  Debug, Operator, Formated_Tilt)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            mkCycleConditions(RefLines, "Dtrv", 0.1, Nmol, Formated_Tilt, mol_pos)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Tilt, "Dtrv", RefLines, 0.1, dirpath,
                                  Debug, Operator, Formated_Tilt)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            MostStable = CompareStructures(RefLines, temp_Structures[-1])

        MostStable = False
        printf(f"{Color.GREEN}\n**********\nTransition in 0.05-Å increments.\n{Color.RESET}")
        while not MostStable:
            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            temp_Structures.append(RefLines)
            mkCycleConditions(RefLines, "Dcol", 0.05, Nmol, Formated_Tilt, mol_pos)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Tilt, "Dcol", RefLines, 0.05, dirpath,
                                  Debug, Operator, Formated_Tilt)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min.txt")
            mkCycleConditions(RefLines, "Dtrv", 0.05, Nmol, Formated_Tilt, mol_pos)
            getTemporaryStructure(MaterName, Nmol, mol_pos, Tilt, "Dtrv", RefLines, 0.05, dirpath,
                                  Debug, Operator, Formated_Tilt)

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

        printf(f"\n"
               f"{Color.GREEN}Local minimum values were successfully found at '{len(MinConditions)}' com files "
               f"for minimum energies were copied into {tcalpath} folder{Color.RESET}\n"
               f"\n"
               f"**********\n"
               f"Making XYZ files...\n")
        mkXYZfile(tcalpath, Debug)

    if calculation_tcal_Flag:
        if not os.path.exists(f"{tcalpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_tcal.log"):
            printf("**********\n"
                   f"{Color.GREEN}Calculating transfer integrals...\n{Color.RESET}")
            XYZs = glob.glob("{tcalpath}/*.xyz")
            for XYZ in XYZs:
                if "_m1.xyz" in XYZ or "_m2.xyz" in XYZ or "-12.xyz" in XYZ or "-23.xyz" in XYZ or "-31.xyz" in XYZ:
                    execute(["rm", "-r", XYZ], True)
                else:
                    pass
            XYZ_3mol_to_XYZ_2mol(tcalpath, Debug)

            with open(f"{tcalpath}/tcal.sh", "w") as f:
                f.write(Stereotyped.tcal_sh_txt)
            subprocess.run(["qsub", "tcal.sh"], cwd=tcalpath)
            MyJobIDList = [Running_JobIDList()[-1]]
            RunningList = Running_JobIDList()
            Wait_minutes = 2

            Flag, wait_job_count = check_jobs(RunningList, MyJobIDList)
            start_time = datetime.datetime.now()
            if Flag:
                formated_ST = start_time.strftime("%m/%d %H:%M:%S")
                printf(f"Calculations for transfer integrals in {MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d "
                       f"were submitted at {formated_ST}.")
                printf(f"\n"
                       f"{Color.GREEN}Wait until jobID {MyJobIDList[-1]}!!{Color.RESET}\n"
                       f"start ID: {RunningList[0]}\n"
                       f"\n")
            else:
                pass
            while Flag:
                Flag, wait_job_count = check_jobs(RunningList, MyJobIDList)
                formated_NOW, elapsed_time = getElapsedTime(start_time)
                Now_time = datetime.datetime.now()
                formated_end_time = (
                    (Now_time +
                     datetime.timedelta(minutes=(wait_job_count * Wait_minutes))).strftime("%m/%d %H:%M:%S"))
                sys.stdout.write(
                    "\033[1F\033[G%s" %
                    f"\t{formated_NOW} ({elapsed_time} min. passed): Job {RunningList[0]} is in progress.            \n"
                    f"\tNext Check >>> {Wait_minutes * wait_job_count} minute later! forecast: "
                    f"{formated_end_time}         ")
                sys.stdout.flush()
                time.sleep(Wait_minutes * wait_job_count * 60)
                RunningList = Running_JobIDList()
                Flag, wait_job_count = check_jobs(RunningList, MyJobIDList)
            else:
                printf(f"\n{Color.GREEN}Calculations until JobID {MyJobIDList[-1]} were finished.\n\n\n{Color.RESET}")
                if Debug :
                    pass
                else:
                    rmWildCards(f"{tcalpath}/*.sh*")
            subprocess.run(["rename", "tcal", f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_tcal", "tcal.log"],
                           cwd=tcalpath)
            readlog(tcalpath, MaterName, Nmol, mol_pos, Formated_Tilt)
        else:
            printf()
    else:
        pass

    printf("\n**********\nPhase Checking...")
    chkKEY = "-12"
    CalCores = mkCalCoreList(tcalpath, chkKEY)
    if len(CalCores) == 0:
        printf("Any chk file is NOT found in the specified folder.\n"
               "The process for phase check was skipped.")
    else:
        with open(f"{tcalpath}/PhaseCheck.txt", "w") as file:
            file.write("Name\tfor LUMO\tfor HOMO\n")
            for CalCore in CalCores:
                mkCubeFile(tcalpath, f"{CalCore}_m1.chk")
                mkCubeFile(tcalpath, f"{CalCore}_m2.chk")

                HomoChk = ComparePhase(tcalpath, CalCore, "HOMO")
                LumoChk = ComparePhase(tcalpath, CalCore, "LUMO")

                file.write(f"{CalCore}\t{LumoChk}\t{HomoChk}\n")
        if Debug:
            pass
        else:
            rmWildCards(f"{tcalpath}/*.chk")
            rmWildCards(f"{tcalpath}/*d-*.log")
            rmWildCards(f"{tcalpath}/*.gjf")

    printf(f"\n**********\n{Color.GREEN} Resulting Data Set...{Color.RESET}\n")
    result_name = f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d"
    printf(f"For {result_name}")
    resultpath = f"./{result_name}_results"
    os.makedirs(resultpath, exist_ok=True)

    MinFileName = f"{result_name}_min.txt"
    TIFileName = f"{result_name}_TIs.txt"
    AllFileName = f"{result_name}_all.txt"
    PCFileName = "PhaseCheck.txt"

    Lacks = []

    if os.path.isfile(MinFileName) and os.path.isfile(f"{tcalpath}/{TIFileName}"):
        combineData(MaterName, Nmol, tcalpath, resultpath, MinFileName, TIFileName, PCFileName, Formated_Tilt, mol_pos)
    else:
        copy_file(MinFileName, f"{resultpath}/{MinFileName}", Lacks)
        printf(f"\n{MinFileName} was copied into the {resultpath} directory. ")

    copy_file(AllFileName, f"{resultpath}/{AllFileName}", Lacks)

    copy_file(f"{MaterName}.xyz", f"{resultpath}/{MaterName}.xyz", Lacks)

    condition_list = f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt"
    copy_file(condition_list, f"{resultpath}/{condition_list}", Lacks)

    copy_file("CalcSetting_HB.txt", f"{resultpath}/CalcSetting_HB.txt", Lacks)

    copy_file(f"{result_name}_min.txt", f"{resultpath}/{result_name}_min.txt", Lacks)

    if "3mol" in Nmol:
        copy_file(f"{result_name}_mins.hist", f"{resultpath}/{result_name}_mins.hist", Lacks)

    if "2mol" in Nmol:
        struct_files = glob.glob(f"{tcalpath}/*.xyz")
    elif "3mol" in Nmol:
        struct_files = glob.glob(f"{tcalpath}/*.all")
    else:
        struct_files = []

    if not struct_files:
        printf("Any file for aggregation structure is not found in the specific folder.")
        Lacks.append("Structural xyz files")
    else:
        printf(f"Files were copied into the {resultpath} folder.")
        for file in struct_files:
            name = os.path.basename(file).replace(".all", "")
            copy_file(file, f"{resultpath}/{name}.xyz", Lacks)

    if Lacks:
        printf("\n\n***********\nSeveral result files were not found in the specific folder.")
        for content in Lacks:
            printf(f"\t{content}")
    else:
        printf(f"\n{Color.GREEN}The {MaterName}_{Nmol} data set were collected in {resultpath}.{Color.RESET}")

    printf("\n************************* ALL PROCESSES END *************************\n")


def printfList(List):
    """リストの要素を順に出力する。

    Args:
    List (list): 出力する要素を持つリスト。

    Returns:
    None
    """
    for content in List:
        printf(content)
    return


def execute(command_list, TEXT):
    """
    Executes a command in the shell using the subprocess module.

    This function takes a list of command arguments and executes it as a
    shell command. If the `TEXT` flag is set to `True`, the command string
    is printed before execution. The function captures the standard output
    and standard error streams and returns them. If the command returns a
    non-zero exit code, the error output is printed.

    :param command_list: List of command and arguments to be executed.
    :type command_list: list of str
    :param TEXT: A flag to print the command before execution.
    :type TEXT: bool
    :return: The result of the executed command, including output and error.
    :rtype: subprocess.CompletedProcess
    :raises subprocess.CalledProcessError: If the command returns a non-zero
                                           exit code.
    """
    command = ' '.join(command_list)
    if TEXT:
        printf(f'> {command}')

    res = subprocess.run(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    # python3.6

    # res = subprocess.run(command_list, capture_output=True, text=True)
    # python3.7 or later

    if res.returncode:
        printf(res.stderr.strip())


# ヘルプが発生しているかを調べ、発生している場合は終了し、もしヘルプが発生していなければ、リストをクリアする関数
def help_check_exit(messages, HelpList):
    """
    ヘルプが発生しているかを調べ、発生している場合は終了し、もしヘルプが発生していなければ、リストをクリアする関数
    :param messages: 表示するメッセージ
    :type messages: list
    :param HelpList: ヘルプリスト
    :type HelpList: list
    """
    message_show(messages)
    if True in HelpList:
        print(f"{Color.RED}{Stereotyped.AbnormalEnd}{Color.RESET}")
        exit()
    else:
        pass
    return None


# messageを表示する関数
def message_show(messages):
    """
    messageを表示する関数
    :param messages: 表示するメッセージ
    :type messages: list
    """
    for message in messages:
        print(message)
    messages.clear()
    return None


def arg_parser(messages, HelpList):
    """
    Parse the command line arguments.
    :param messages:
    :param HelpList:
    :return:
    """
    Help_Text = ("Structural optimisation and transfer integration of HerringBone structures.\n"
                 "Pass [molecule name].xyz as the first argument.\n"
                 "Pass the number of mol as the first argument.\n"
                 "If required, you will be asked to enter the Tilt angle, a name that uniquely identifies you, "
                 "and a structural pattern.")

    parser = argparse.ArgumentParser(description=Help_Text, formatter_class=argparse.RawTextHelpFormatter)
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
        HelpList.append(False)
    else:
        MaterName = ""
        messages.append(f"{Color.RED} file selected. [molecule name].xyz is required.{Color.RESET}")
        HelpList.append(True)

    help_check_exit(messages, HelpList)
    return args, MaterName


def mkCheckFile(args, Debug, messages, HelpList, MaterName):
    printf("*************** Structural checks are in place. ***************\n")
    printf("Checking if 'CalcSetting_HB.txt' exists.")
    if os.path.exists("./CalcSetting_HB.txt"):
        messages.append("CalcSetting_HB.txt: Found")
        HelpList.append(False)
    else:
        messages.append("CalcSetting_HB.txt: NOT Found\n"
                        "Make the correct CalcSetting_HB.txt and then restart the program.")
        HelpList.append(True)
        with open("CalcSetting_HB.txt", "w") as file:
            file.write("Column Direction: [x, y, or z]\n"
                       "Transverse Direction: [x, y, or z]\n"
                       "Tilt Axis: [x, y, or z]\n"
                       "Tilt Angle: [Number]\n"
                       "\n"
                       "[Comment]\n")

    Flag_XYZ = args.xyz

    help_check_exit(messages, HelpList)

    Tilt_Angle = 0

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
            Tilt_Angle = float(Tilt)
            break
        except (IndexError, ValueError):
            printf("Please enter a number.")
            Tilt = input("\nPlease enter the Tilt angle at degree."
                         "\n>>> ")
            continue
    printf(f"\n{Tilt_Angle}-degree has been entered as Tilt")
    Formated_Tilt = int(Tilt_Angle * 10)  # For File Name
    printf(f"\n"
           f"*************** Caution!!! The Tilt angle multiplied by 10 is used in the file name. "
           f"***************\n"
           f"*************** Like: {MaterName}_3molp1_t{Formated_Tilt}d_10d-500-800.gjf\n")

    help_check_exit(messages, HelpList)
    file_path = f"./StructCheck-{MaterName}-t{Formated_Tilt}d"
    os.makedirs(file_path, exist_ok=True)

    if not args.manual:
        file_path = f"./StructCheck-{MaterName}-t{Formated_Tilt}d"
        os.makedirs(file_path, exist_ok=True)
        with open(f"{file_path}/G.sh", "w") as f:
            f.write(Stereotyped.Sh_txt)
        Temp_SHs = []
        Temp_SH = mkFiles(MaterName, "2mol", "", "90d-400", "",
                          file_path, Tilt_Angle, Formated_Tilt, Flag_XYZ)
        File_Name = f"{MaterName}_2mol_t{Formated_Tilt}d_90d-400"
        printf(f"{file_path}/{File_Name}.gjf have been created.\n")
        Temp_SHs.append(Temp_SH.replace("qsub", "").strip())

        Temp_SH = mkFiles(MaterName, "2mol", "", "60d-420", "",
                          file_path, Tilt_Angle, Formated_Tilt, Flag_XYZ)
        File_Name = f"{MaterName}_2mol_t{Formated_Tilt}d_60d-420"
        printf(f"{file_path}/{File_Name}.gjf have been created.\n")
        Temp_SHs.append(Temp_SH.replace("qsub", "").strip())

        Temp_SH = mkFiles(MaterName, "2mol", "", "30d-600", "",
                          file_path, Tilt_Angle, Formated_Tilt, Flag_XYZ)
        File_Name = f"{MaterName}_2mol_t{Formated_Tilt}d_30d-600"
        printf(f"{file_path}/{File_Name}.gjf have been created.\n")
        Temp_SHs.append(Temp_SH.replace("qsub", "").strip())

        Temp_SH = mkFiles(MaterName, "3mol", "p1", "60d-420-600", "",
                          file_path, Tilt_Angle, Formated_Tilt, Flag_XYZ)
        File_Name = f"{MaterName}_3molp1_t{Formated_Tilt}d_60d-420-600"
        printf(f"{file_path}/{File_Name}.gjf have been created.\n")
        Temp_SHs.append(Temp_SH.replace("qsub", "").strip())

        mkFiles(MaterName, "3mol", "p2", "60d-420-600", "",
                file_path, Tilt_Angle, Formated_Tilt, Flag_XYZ)

        File_Name = f"{MaterName}_3molp2_t{Formated_Tilt}d_60d-420-600"
        printf(f"{file_path}/{File_Name}.gjf have been created.\n")

        mkFiles(MaterName, "3mol", "p3", "60d-420-600", "",
                file_path, Tilt_Angle, Formated_Tilt, Flag_XYZ)
        File_Name = f"{MaterName}_3molp3_t{Formated_Tilt}d_60d-420-600"
        printf(f"{file_path}/{File_Name}.gjf have been created.\n")

        Temp_SHs.append("G.sh")
        for SH in Temp_SHs:
            subprocess.run(["rm", SH], timeout=10, cwd=file_path)
    else:
        printf("Start manual mode.")
        with open(f"{file_path}/G.sh", "w") as f:
            f.write(Stereotyped.Sh_txt)
        while True:
            structure = input("\n"
                              "Select the structure to be created.\n"
                              "1: 2mol\n"
                              "2: 3molp1 (H-H)\n"
                              "3: 3molp2 (T-T)\n"
                              "4: 3molp3(H-T)\n"
                              "5: exit\n"
                              ">>> ")
            if structure == "5":
                EXIT = input("\n"
                             "Do you want to exit the programme?\n"
                             "y/n >>> ")
                if EXIT == "y":
                    subprocess.run(["rm", "G.sh"], timeout=10, cwd=file_path)
                    if Debug:
                        printf("*************** Debug Finished!!! ***************")
                    printf("\n************************* ALL PROCESSES END *************************\n")
                    exit()
                else:
                    continue
            file_path = f"./StructCheck-{MaterName}-t{Formated_Tilt}d"
            os.makedirs(file_path, exist_ok=True)

            while True:
                Rotate_Angle = input("\n"
                                     "Please Enter the Rotate Angle.\n"
                                     ">>> ")
                try:
                    Rotate_Angle = int(Rotate_Angle)
                    break
                except (IndexError, ValueError):
                    printf("Please enter a number.")
                    continue
            while True:
                D_Col = input("\n"
                              "Please Enter the Column Direction.\n"
                              ">>> ")
                try:
                    D_Col = float(D_Col)
                    break
                except (IndexError, ValueError):
                    printf("Please enter a number.")
                    continue
            if structure == "1":
                Temp_Condition = f"{Rotate_Angle}d-{int(D_Col * 100)}"
                Temp_SH = mkFiles(MaterName, "2mol", "", Temp_Condition, "",
                                  file_path, Tilt_Angle, Formated_Tilt, Flag_XYZ)
                File_Name = f"{MaterName}_2mol_t{Formated_Tilt}d_{Temp_Condition}"
                printf(f"{file_path}/{File_Name}.gjf have been created.\n")
                subprocess.run(["rm", Temp_SH.replace("qsub", "").strip()], timeout=10, cwd=file_path)
                continue
            else:
                while True:
                    D_Transv = input("\n"
                                     "Please Enter the Transverse Direction.\n"
                                     ">>> ")
                    try:
                        D_Transv = float(D_Transv)
                        break
                    except (IndexError, ValueError):
                        printf("Please enter a number.")
                        continue
                mol_pos = ""
                if structure == "2":
                    mol_pos = "p1"
                elif structure == "3":
                    mol_pos = "p2"
                elif structure == "4":
                    mol_pos = "p3"
                Temp_Condition = f"{Rotate_Angle}d-{int(D_Col * 100)}-{int(D_Transv * 100)}"
                Temp_SH = mkFiles(MaterName, "3mol", mol_pos, Temp_Condition, "",
                                  file_path, Tilt_Angle, Formated_Tilt, Flag_XYZ)
                File_Name = f"{MaterName}_3mol{mol_pos}_t{Formated_Tilt}d_{Temp_Condition}"
                printf(f"{file_path}/{File_Name}.gjf have been created.\n")
                subprocess.run(["rm", Temp_SH.replace("qsub", "").strip()], timeout=10, cwd=file_path)
    printf("\n************************* ALL PROCESSES END *************************\n")
    exit()


def getRefLines(FileName):
    """
    指定されたファイルから参照行を取得する関数

    ファイルを読み込み、各行をリストに追加します。ValueErrorが発生した場合は無視します。

    Args:
        FileName (str): 参照行を取得するファイル名

    Returns:
        list: 参照行のリスト

    Raises:
        FileNotFoundError: 指定されたファイルが存在しない場合

    Examples:
        >>> getRefLines("ref_data.txt")
        ['line1\n', 'line2\n', 'line3\n']
    """
    with open(FileName, "r") as file:
        lines = file.readlines()

    RefLines = []
    for line in lines:
        contents = line.strip().split()
        try:
            float(contents[0])  # フロートに変換できない行(コメントとか)を飛ばすため
            RefLines.append(line)
        except ValueError:
            pass
    return RefLines


def mkConditionFile(Nmol, mol_pos, which, RefLines, Formated_Tilt, dev):
    """条件ファイルを作成する。

    InitialConditionファイルから、指定した傾斜角(Tilt)における新しい条件ファイル(ConditionList)を作成する。
    既存の条件ファイルが存在する場合は、何もしない。

    Args:
        Nmol (str): 分子数
        mol_pos (str): 分子の位置
        which (str): 条件の種類
        RefLines (list): 参照線のリスト
        Formated_Tilt (float): ファイル名に使われている傾斜角
        dev(float): 何度刻みで計算するか

    Returns:
        None
    """
    if os.path.exists(f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt"):
        printf(f"\nConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt is EXIST!!!")
        return
    else:
        printf(f"\nConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt is NOT EXIST!!!")
        with open(f"./InitialCondition_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt", "r") as f:
            lines = f.readlines()
        NewConditions = []
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
            for Condition in NewConditions:
                f.write(Condition)
        printf(f"New ConditionList_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt was written from "
               f"InitialCondition_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt!!")
        return


def getRefValues(Deg, RefLines):
    """タブ区切りファイルから指定角度に一致する参照値を取得する。

    指定された角度 (Deg) と一致する行をタブ区切りファイル (RefLines) から探し、
    対応する参照値をリスト形式で返します。一致する行がない場合は ["na", "na"] を返します。
    ファイル形式のエラーや指定角度に一致する行がない場合は空のリストを返します。

    Args:
        Deg (float): 取得したい参照値に対応する角度。
        RefLines (list): タブ区切りファイルの各行を要素とするリスト。

    Returns:
        list: 参照値 [Value1, Value2] または ["na", "na"]。

    Raises:
        ValueError: RefLines の要素が不正な形式の場合。
        IndexError: RefLines の要素がタブ区切りでない場合。
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
    """新しい実験条件の文字列表現を生成する。

    この関数は、分子の種類 (Nmol)、二面角 (Deg)、指定された値 (Val)、チルト角 (Tilt)、
    指定された値の種類 (which)、参照値 (RefValues) をもとに、新しい実験条件の文字列表現を生成します。
    生成された文字列は標準出力に表示され、関数からの戻り値としても返されます。

    Args:
    Nmol (str): 分子の種類 ("2mol" または "3mol")。
    Deg (float): 二面角。
    Val (float): 指定された値 (Dcol または Dtrv)。
    Tilt (float): チルト角。
    which (str): 指定された値の種類 ("Dcol" または "Dtrv")。
    RefValues (list): 参照値のリスト。

    Returns:
    str: 新しい実験条件の文字列表現。

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
    printf(f"\t\t{NewCondition.strip()}")
    return NewCondition


def getTemporaryStructure(MaterName, Nmol, mol_pos, Tilt, which, RefLines, dev, dirpath, Debug, Operator,
                          Formated_Tilt):
    """
    Generates and manages temporary structures and job submissions for
    material simulations.

    This function checks if the required simulation logs exist for a given
    material configuration. If not, it generates the necessary files, submits jobs,
    and monitors their progress until completion. The function also reads
    simulation results and generates new conditions for further calculations
    based on the completed jobs.

    :param MaterName: Name of the material to be simulated.
    :type MaterName: the
    :param Nmol: Number of molecules in the simulation, used in file naming.
    :type Nmol: str
    :param mol_pos: Molecular position information used in file naming.
    :type mol_pos: str
    :param Tilt: Tilt angle for the simulation.
    :type Tilt: float
    :param which: Direction indicator ('Dcol' or 'Dtrv') for the calculation.
    :type which: str
    :param RefLines: Reference lines for condition generation.
    :type RefLines: list
    :param dev: Deviation value used in generating new conditions.
    :type dev: float
    :param dirpath: a Directory path where the files and logs are stored.
    :type dirpath: str
    :param Debug: Flag to indicate whether to keep or remove debug files.
    :type Debug: bool
    :param Operator: Operator information used in file generation.
    :type Operator: str
    :param Formated_Tilt: Formatted string of the tilt angle.
    :type Formated_Tilt: int
    :return: None
    :rtype: NoneType

    :raises FileNotFoundError: If the necessary condition files are missing.
    :raises subprocess.CalledProcessError: If the job submission fails.
    """
    judge = False
    while not judge:
        Conditions = getConditions(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Formated_Tilt}d.txt")
        qsubList = []
        with open(f"{dirpath}/G.sh", "w") as f:
            f.write(Stereotyped.Sh_txt)
        for Condition in Conditions:
            if os.path.exists(f"{dirpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_{Condition}.log"):
                pass
            else:
                qsub_temp = mkFiles(MaterName, Nmol, mol_pos, Condition, Operator, dirpath, Tilt, Formated_Tilt, False)
                qsubList.append(qsub_temp)
        printf("\n**********\nJobs are submitting...")
        if len(qsubList) == 0:
            printf(f"{Color.GREEN}Job was not submitted. Calculations with the conditions might be finished.{Color.RESET}")
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
                printf(f"\n'{int(len(MyJobIDList))}' calculations for '{term}' was submitted!! at {formated_ST}")
                printf(f"\n"
                       f"{Color.GREEN}Wait until jobID {MyJobIDList[-1]}!!{Color.RESET}\n"
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
                printf(f"{Color.GREEN}\n\n"
                       f"Calculation cycles for {which} until JobID {MyJobIDList[-1]} were finished.{Color.RESET}")
                if Debug:
                    pass
                else:
                    rmWildCards(f"{dirpath}/*.sh*")
                    rmWildCards(f"{dirpath}/*.chk")
        printf("\n**********\nReading Data...\n")
        readEnergies(dirpath, MaterName, Nmol, Formated_Tilt, mol_pos)
        judge = mkNewConditionList(MaterName, Nmol, Formated_Tilt, which, dev, RefLines, mol_pos)
    return


def getConditions(FileName):
    """ファイルから条件を読み込み、リストで返す。

    指定されたファイルを読み込み、各行を条件としてリストに格納します。
    ファイルはテキスト形式であることを想定しており、各条件は1行ずつ記述されている必要があります。

    Args:
        FileName (str): 条件が記述されたファイルのパス。

    Returns:
        list: ファイルから読み込まれた条件のリスト。各要素は文字列型です。

    Raises:
        FileNotFoundError: 指定されたファイルが存在しない場合に発生します。
        PermissionError: 指定されたファイルへの読み込み権限がない場合に発生します。
        IOError: ファイルの読み込み中にエラーが発生した場合に発生します。
    """
    Conditions = []
    with open(FileName, "r") as f:
        lines = f.readlines()
    for line in lines:
        Conditions.append(line.strip())
    return Conditions


def mkFiles(MaterName, Nmol, mol_pos, Condition, Operator, dirpath, Tilt_Angle, Formated_Tilt, Flag_XYZ):
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

    Direction_Col, Direction_Transv, Rotate_Axis, Tilt_Axis, rotate = Axis_Setting_HB()
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

    col_transl = mkDirection(D_Col, Direction_Col, "column direction")
    transv_transl = mkDirection(D_Transv, Direction_Transv, "transverse direction")
    Direction_Other = "xyz".replace(Direction_Col, "").replace(Direction_Transv, "")
    a_transl = mkDirection(Constant.d_another, Direction_Other, "other direction")
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
            printf(f"{dirpath}/{File_Name}.xyz have been created.")
    elif "3mol" in Nmol:
        write_gjf_file(f"{dirpath}/{File_Name}.gjf",
                       Headers, Element, Mol1_pos, Mol2_pos, Mol3_pos)
        if Flag_XYZ:
            write_xyz_file(f"{dirpath}/{File_Name}.xyz",
                           Element, Mol1_pos, Mol2_pos, Mol3_pos)
            printf(f"{dirpath}/{File_Name}.xyz have been created.")
    with open(f"{dirpath}/G.sh", "r") as orgSH:
        lines = orgSH.readlines()
    lines[12] = f"g16 {GJF_Name}\n"
    with open(f"{dirpath}/{SH_Name}", "w") as newSH:
        for line in lines:
            newSH.write(line)
    qsub_temp = f"qsub {SH_Name}"
    return qsub_temp


def Axis_Setting_HB():
    """
    CalcSetting_HB.txtから軸設定を読み込む。

    この関数は、CalcSetting_HB.txtファイルから列方向、横方向、ダミーの座標、
    回転軸、および傾斜軸の情報を読み込み、適切な軸の組み合わせを設定します。
    CalcSetting_HB.txtが存在しない場合、新しいファイルを作成し、適切なフォーマット
    の指示を書き込みます。

    Returns:
        tuple: 列方向、横方向、ダミーの座標、回転軸、傾斜軸、回転の組み合わせ。

    Raises:
        FileNotFoundError: CalcSetting_HB.txtが見つからない場合、新しいファイルを作成し
                          プログラムを終了します。
    """
    try:
        with open("CalcSetting_HB.txt", "r") as File:
            lines = File.readlines()
    except FileNotFoundError:
        printf("\nThe CalcSetting_HB.txt is not found in the current directory.\n"
               "Make the correct CalcSetting_HB.txt and then restart the program.\n")
        with open("CalcSetting_HB.txt", "w") as file:
            file.write("Column Direction: [x, y, or z]\n"
                       "Transverse Direction: [x, y, or z]\n"
                       "Tilt Axis: [x, y, or z]\n"
                       "Tilt Angle: [Number]\n"
                       "\n"
                       "[Comment]\n")
        printf(Stereotyped.AbnStop)
        exit()

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
        printf("INVALID PARAMETERS!! Column direction and Transverse direction should be different!!")
        printf(Stereotyped.AbnStop)
        exit()

    return params["column direction"], params["transverse direction"], params["other axis"], params["tilt axis"], rotate


def mkDirection(axis_direction, input_axis, axis_name):
    """指定された軸方向に基づいて、3次元方向ベクトルを生成する。

    Args:
        axis_direction (float): 軸方向の大きさ。
        input_axis (str): 軸方向 ("x", "y", "z")。
        axis_name (str): 軸の名前 (エラーメッセージ用)。

    Returns:
        np.ndarray: 3次元方向ベクトル (例: [1.0, 0.0, 0.0])。

    Raises:
        SystemExit: 指定された軸方向が不正な場合。
    """
    try:
        axis_direction = float(axis_direction)
    except ValueError:
        printf(f"\nThe axis direction for {axis_name} is INVALID: {axis_direction}")
        printf(Stereotyped.AbnStop)
        exit()
    directions = {"x": [axis_direction, 0.0, 0.0], "y": [0.0, axis_direction, 0.0], "z": [0.0, 0.0, axis_direction]}
    if input_axis not in directions:
        printf(f"\nThe axis of {axis_name} is INVALID.")
        printf(Stereotyped.AbnStop)
        exit()
    return np.array(directions[input_axis])


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
        atm_m2 = Rotate(Position, *Mol2_Angles, rotate) + Translations[0]
        Mol2.append(atm_m2)
        atm_m3 = Rotate(Position, *Mol3_Angles, rotate) + Translations[0] / 2 + Translations[1]
        Mol3.append(atm_m3)
    return Elist, Mol1, Mol2, Mol3, NinMol


def Rotate(Current, Tx, Ty, Tz, rotation):
    """3次元座標を指定された軸と角度で回転させる。

    Args:
        Current (np.ndarray): 回転させたい3次元座標 (形状: (3,)).
        Tx (float): x軸周りの回転角度 (単位: 度).
        Ty (float): y軸周りの回転角度 (単位: 度).
        Tz (float): z軸周りの回転角度 (単位: 度).
        rotation (str): 回転軸の順序を表す文字列 ('x', 'y', 'z' の組み合わせ).

    Returns:
        np.ndarray: 回転後の3次元座標 (形状: (3,)).

    Raises:
        ValueError: rotation に無効な軸が含まれている場合.
        TypeError: Current が np.ndarray でない場合.
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


def readEnergies(dir_path, MaterName, Nmol, Formated_Tilt, mol_pos):
    """指定ディレクトリ内のログファイルからエネルギー情報を読み込み、ファイルに出力する。

    指定されたディレクトリ内のログファイルを解析し、角度ごとのエネルギー情報を抽出し、
    全結果ファイルと最小エネルギーファイルに出力します。
    正常終了していないログファイルは削除されます。

    Args:
        dir_path (str): ログファイルが保存されているディレクトリパス。
        MaterName (str): 解析対象の物質名。
        Nmol (str): 解析対象の分子数。
        Formated_Tilt (int); チルト角
        mol_pos (str): 3molの場合の構造

    Returns:
        None

    Raises:
        FileNotFoundError: 指定されたディレクトリパスが存在しない場合。
        ValueError: ログファイル名が不正な形式の場合。
    """
    LogList = glob.glob(f"{dir_path}/{MaterName}_{Nmol}{mol_pos}_*.log")
    LogList.sort()
    DegList = []
    for Log in LogList:
        FileName, deg, Vdcol, Vdtrv = getVALfromLogName(Nmol, Log)
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

            printf(f"\t{header}")
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
                            printf(f"{Color.RED}\t{Log} was NOT normally terminated. It was removed.{Color.RESET}")
                            # 正常終了しなかったものをログとして残す
                            with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_error.log", "a") as error_log:
                                error_log.write(f"{Log.split("/")[-1].split(".")[0].split("_")[-1]}\n")
                            # 正常終了しなかったログを確認し、同じログが２個以上あったらプログラムを終了する
                            with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_error.log", "r") as error_log:
                                error_lines = error_log.readlines()
                            error_lines = list(set(error_lines))
                            if len(error_lines) != len(set(error_lines)):
                                printf(f"{Color.RED}\nThere are the same error logs in the error log file. The program will be terminated.{Color.RESET}")
                                exit()

                            os.remove(Log)

                minkey = min(CPE_Dict, key=CPE_Dict.get)
                MinimumL.append(minkey)
                MinData.write(VAL_Dict.get(minkey))
                Keys.sort()
                for Key in Keys:
                    if Key == minkey:
                        sentence = f"*\t{VAL_Dict.get(Key)}"
                        AllData.write(sentence)
                        printf(sentence.strip())

                    else:
                        sentence = f"-\t{VAL_Dict.get(Key)}"
                        AllData.write(sentence)
                        printf(sentence.strip())
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

    if "2mol" in Nmol and Condition.count("-") == 1:
        deg = Conditions[0].replace('d', '')
        dcol = Conditions[1]
        dtrv = "na"
    elif "3mol" in Nmol and Condition.count("-") == 2:
        deg = Conditions[0].replace('d', '')
        dcol = Conditions[1]
        dtrv = Conditions[2]
    else:
        printf("\nUNEXPECTED ERROR happens in the function of getVALfromLogName!!")
        exit()

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


def mkNewConditionList(MaterName, Nmol, Formated_Tilt, which, dev, RefLines, mol_pos):
    """探索範囲の各角度に対して、指定された基準値との比較を行い、条件リストを更新する。

    Args:
        MaterName (str): 解析対象の物質名。
        Nmol (str): 解析対象の分子数 ("2mol" or "3mol")。
        Formated_Tilt (int): ファイル名に使われるTilt角
        which (str): 解析対象の物理量 ("Dcol" or "Dtrv")。
        dev (float): 基準値からの許容範囲。
        RefLines (list): 基準値が記載されたファイルの行リスト。
        mol_pos (str): 3モルの場合の回転角

    Returns:
        bool: 全ての角度で極小値が見つかった場合はTrue、そうでない場合はFalse。
    """
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

    Judges = []
    ComplDeg = []
    ComplVal = []
    NewConditions = []

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

        if round(SV + dev, 2) in ValueList and round(SV - dev, 2) in ValueList:
            Judges.append("complete")
            printf(f"\nThe local minimum by {dev} step for {Deg} degree was Found;\t{SV}.")
            printf(f"{Color.GREEN}\tCOMPLETE!!{Color.RESET}")
            ComplDeg.append(Deg)
            ComplVal.append(SV)
        elif SV == min(ValueList):
            Judges.append("not complete")
            printf(f"\nLocal minimum for {Deg} degree was NOT FOUND in the cycle.")
            printf(f"\tNOT COMPLETE(1)!! {Constant.Cn} new conditions bellow were appended to the ConditionList.txt.")
            for i in range(Constant.Cn):
                NewCondition = mkNewCondition(Nmol, Deg, SV - dev * (i + 1), which, RefValues)
                NewConditions.append(NewCondition)
        elif SV == max(ValueList):
            Judges.append("not complete")
            printf(f"\nLocal minimum for {Deg} degree was NOT FOUND in the cycle.")
            printf(f"\tNOT COMPLETE(2)!! {Constant.Cn} new conditions bellow were appended to the ConditionList.txt.")
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
        printf("\nLocal minimum have NOT been FOUND in any angle.\n")
    else:
        printf(f"\nLocal minimum was FOUND in {len(ComplDeg)} angles.\n")
        if "2mol" in Nmol or "3mol" in Nmol:
            printf(f"Angle\t{which}")
        else:
            pass
        for i in range(len(ComplDeg)):
            printf(f"{ComplDeg[i]}\t{ComplVal[i]}")

    Judges = list(set(Judges))
    if "not complete" in Judges:
        judge = False
    elif "complete" in Judges and len(Judges) == 1:
        judge = True
    else:
        judge = False

    return judge


def mkCycleConditions(RefLines, which, dev, Nmol, Formated_Tilt, mol_pos):
    """次のサイクルの実験条件を生成し、既存条件と合わせてファイルに保存する。

    既存の実験条件ファイルを読み込み、指定されたパラメータに基づいて新たな実験条件を生成します。
    生成された条件は既存の条件と重複しないように調整され、ソートされた上でファイルに書き戻されます。

    Args:
        RefLines (list): 参照となる実験条件のリスト（各要素は "角度 実験値1 実験値2..." の形式の文字列）。
        which (str): 生成する実験条件の種類 ("Dcol" or "Dtrv")。
        dev (float): 実験値からの偏差幅。
        Nmol (str): サンプルのモル数。
        Formated_Tilt (int): サンプルの傾斜角度。
        mol_pos (str): 3molの場合の配置

    Returns:
        None

    Raises:
        ValueError: which が "Dcol" または "Dtrv" 以外の値の場合。
    """
    Degs = []
    NewConditions = []
    if dev == 0.1:
        n = 1
    elif dev == 0.05:
        n = 1
    else:
        n = 1

    for RefLine in RefLines:
        Contents = RefLine.strip().split()
        Degs.append(round(float(Contents[0]), 1))
    printf("\nNew conditions for the next cycle:")
    for Deg in Degs:
        RefValues = getRefValues(Deg, RefLines)
        if which == "Dcol":
            SV = round(float(RefValues[0]), 2)
        elif which == "Dtrv":
            SV = round(float(RefValues[1]), 2)
        else:
            SV = int()

        for i in range(n):
            NewCondition = mkNewCondition(Nmol, float(int(Deg)), SV - dev * (n - i), which, RefValues)
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
    """
    比較対象の構造リストを参照構造リストと比較し、最も安定しているかを判定する。

    Args:
        RefLines (list of str): 参照構造リスト。各要素はスペースで区切られた文字列。
        LinesBefore (list of str): 比較対象の構造リスト。各要素はスペースで区切られた文字列。

    Returns:
        bool: すべての比較対象が参照構造と一致している場合は `True`、そうでない場合は `False`。

    Example:
        RefLines = ["1 A B", "2 C D"]
        LinesBefore = ["1 A B", "2 C D"]
        CompareStructures(RefLines, LinesBefore)
        True

        RefLines = ["1 A B", "2 C D"]
        LinesBefore = ["1 A B", "2 X Y"]
        CompareStructures(RefLines, LinesBefore)
        False
    """
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
    """ファイルから最小条件のリストを取得する。

    指定されたファイル (MaterName_Nmol_Tilt_min.txt) から、
    最小条件を表す文字列のリストを取得します。
    ファイルの最初の2行はヘッダーとして無視されます。

    Args:
        MaterName (str): 材料名。
        Nmol (str): 分子数。
        Formated_Tilt (int): Tiltの値。
        mol_pos (str): 3molの場合の構造

    Returns:
        list: 最小条件のリスト。各要素は "Deg-tTilt-Dcol" または
              "Deg-tTilt-Dcol-Dtrsv" の形式の文字列。

    Raises:
        FileNotFoundError: 指定されたファイルが存在しない場合。
    """
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


def mkXYZfile(tcalpath, Debug):
    """
    指定されたパス内のすべての.gjfファイルを.xyzファイルに変換する。

    指定されたフォルダ内の.gjfファイルをすべて検索し、それらを.xyz
    ファイルに変換します。まず.gjfファイルを読み込み、一部の行の内容を
    置換して.gjfファイルとして保存し、その後、外部コマンドを使用して
    .xyzファイルに変換します。

    パラメータ:
        tcalpath (str): .comファイルが格納されているフォルダのパス。

    例外:
        ファイルが存在しない場合、または変換に失敗した場合にエラーが発生
        する可能性があります。
    """
    F_paths = glob.glob(f"{tcalpath}/*.com")
    if not F_paths:
        printf(f"There is NO .com file in the {tcalpath} folder.")
        return

    printf(f"'{len(F_paths)}' files were transformed to XYZ files!!")
    for F_path in F_paths:
        process_file(F_path)
    """
    # cartcomtoxyzが存在する場合はこちらの方が安定する
    printf(f"Converting .com files to .xyz\n")
    subprocess.run(["cartcomtoxyz"], cwd=tcalpath, timeout=3)
    """

    # cartcomtoxyzが存在しない場合or組み込むことにメリットがある場合はこっち
    # 動作が不安定（原因究明中）のため利用を停止して従来の方法で実行中
    printf(f"Converting .com files to .xyz\n")

    ComFileNames = glob.glob(f"{tcalpath}/*.com")
    XYZFileNameList = []
    for ComFileName in ComFileNames:
        FormatedComfn = (ComFileName[:-4]).strip().split("/")[-1]
        XYZFileName = f"{FormatedComfn}.xyz"
        if Debug:
            printf(f"\t{tcalpath}/{XYZFileName}")
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


def process_file(F_path):
    """
    ダミー原子を「C」に置換する

    指定されたファイルを読み込み、行の先頭が'X'で始まる行の'X'を'C'に置換し、
    新しいファイルに保存する。

    パラメータ:
        F_path (str): 処理対象のファイルのパス。

    戻り値:
        なし

    発生する例外:
        FileNotFoundError: 指定されたファイルが見つからない場合に発生。
        IOError: 入出力操作に失敗した場合に発生。
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


def XYZ_3mol_to_XYZ_2mol(tcalpath, Debug):
    """3分子系のXYZファイルから2分子系のXYZファイルを3つ作成する。

    指定されたディレクトリ内の"_3mol"を含むXYZファイルを処理し、各分子ペアのXYZファイルを生成します。
    処理されたファイルは".all"にリネームされ、原子数が3の倍数でないファイルはエラーメッセージと共に報告されます。

    Args:
    tcalpath (str): 対象のXYZファイルが格納されているディレクトリへのパス。

    Returns:
    None
    """
    filepathes = glob.glob(f"{tcalpath}/*_3mol*.xyz")

    Faults = []
    for filepath in filepathes:
        Mol1 = []
        Mol2 = []
        Mol3 = []
        filecore = filepath[:-4]
        with open(filepath, "r") as f:
            number_atoms = f.readline()
            if Debug:
                printf(f"filename: {filepath}")
                printf(f"number_atoms: {number_atoms}")
            Comment = f.readline()
            coordinates = f.readlines()
        if int(float(number_atoms)) % 3 == 0:
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
                if Debug:
                    printf(f"\t{filecore}-12.xyz")
                    printf(f"\t{filecore}-23.xyz")
                    printf(f"\t{filecore}-31.xyz")
            os.rename(filepath, f"{filecore}.all")
            if Debug:
                printf(f"\t{filecore}.xyz has been renamed {filecore}.all.")
        else:
            Faults.append(filepath)
    if len(Faults) == 0:
        printf("\nAll xyz files for 3mol were devided into 3 xyz files for 2mol.")
    else:
        for Fault in Faults:
            printf(f"{Fault} may NOT be an xyz file for 3mol. Check it.")
    return


def readlog(tcalpath, MaterName, Nmol, mol_pos, Formated_Tilt):
    """
    指定されたログファイルを読み取り、必要なデータを抽出する。

    指定されたパスにあるログファイルを読み込み、特定のキーワード
    ("Input File Name:", "NLUMO", "LUMO", "HOMO", "NHOMO")が含まれる行を
    抽出してリストに格納する。データが正しく抽出された後、結果を
    新しいテキストファイルに出力する。

    Parameters:
        tcalpath (str): ログファイルのディレクトリパス。
        MaterName (str): 材料名。
        Nmol (str): 分子の数。
        mol_pos (str): 分子の位置。
        Formated_Tilt (int): ファイル名に使われるTilt角

    Returns:
        None

    Raises:
        SystemExit: データの抽出過程でエラーが発生した場合に発生。
    """
    filepath = f"{tcalpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_tcal.log"
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
        printf("There is an UNEXPECTED ERROR in read log process.")
        exit()

    output_filepath = f"{tcalpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_TIs.txt"
    with open(output_filepath, "w") as file:
        header = (f"***** Transfer Integrals in "
                  f"{tcalpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_tcal.log *****\n")
        columns = "input file\tNLUMO (meV)\tLUMO (meV)\tHOMO (meV)\tNHOMO (meV)\n"
        file.write(header + columns)
        printf(header.strip())
        printf(columns.strip())

        for i in range(0, len(DataList), 5):
            input_file = DataList[i].split(": ")[1].split(".xyz")[0]
            NLUMO = extract_value(DataList[i + 1], "NLUMO")
            LUMO = extract_value(DataList[i + 2], "LUMO")
            HOMO = extract_value(DataList[i + 3], "HOMO")
            NHOMO = extract_value(DataList[i + 4], "NHOMO")
            file.write(f"{input_file}\t{NLUMO}\t{LUMO}\t{HOMO}\t{NHOMO}\n")
            printf(f"{input_file}\t{NLUMO}\t{LUMO}\t{HOMO}\t{NHOMO}")


def extract_value(line, keyword):
    return float(line.split(keyword)[-1].split()[0])


def mkCalCoreList(tcalpath, chkKEY):
    """
    指定されたパス内のファイルリストから計算コアリストを作成する。

    指定されたディレクトリ内のファイル名に基づいて、特定のキーを含む
    ファイルを検索し、"_m1.chk"または"_m2.chk"を含むファイル名の
    部分文字列を抽出して、重複を排除したリストを返す。

    パラメータ:
    tcalpath (str): ファイル検索を行うディレクトリのパス。
    chkKEY (str): 検索対象のファイル名に含まれるキー。

    戻り値:
    list: "_m1.chk"または"_m2.chk"を含むファイル名から抽出された部分文字列
          の重複を排除したリスト。

    例外:
    発生しない。
    """
    FileList = glob.glob(f"{tcalpath}/*{chkKEY}*")
    CalCores = []
    for file in FileList:
        if "_m1.chk" in file or "_m2.chk" in file:
            CalCores.append(file[file.rfind("/") + 1:-7])
        else:
            pass
    CalCores = list(set(CalCores))
    return CalCores


def mkCubeFile(tcalpath, chkfile):
    """
    指定されたchkファイルからHOMOおよびLUMOのcubファイルを生成する。

    この関数は、与えられたchkファイルをformchkを用いてfchファイルに変換し、
    そのfchファイルからcubegenを用いてHOMOおよびLUMOのcubファイルを生成します。
    最後に生成されたfchファイルを削除します。

    パラメータ:
    tcalpath (str): コマンドの実行ディレクトリのパス。
    chkfile (str): 入力となるchkファイルのパス。

    戻り値:
    なし

    発生する例外:
    subprocess.TimeoutExpired: コマンドの実行がタイムアウトした場合。
    """

    def run_command(command, description):
        subprocess.run(command, cwd=tcalpath, timeout=10)
        printf(f"\t{description}: {Color.GREEN}Completed!!{Color.RESET}")

    printf(f"For {chkfile} ...")
    fchk = f"{chkfile[:-4]}.fch"
    HomoCub = f"{chkfile[:-4]}_HOMO.cub"
    LumoCub = f"{chkfile[:-4]}_LUMO.cub"

    run_command(["formchk", chkfile, fchk], f"File Conversion ({chkfile} -> {fchk})")
    run_command(["cubegen", "0", "MO=Homo", fchk, HomoCub, "-2", "h"], f"{HomoCub} has been created")
    run_command(["cubegen", "0", "MO=Lumo", fchk, LumoCub, "-2", "h"], f"{LumoCub} has been created")
    run_command(["rm", fchk], f"Removing {fchk}")
    return


def ComparePhase(tcalpath, CalCore, HomoLumo):
    """2つのCUBEファイルの位相を比較し、結果を返す。

    tcalpathディレクトリ内の `CalCore_m1_{HomoLumo}.cub` と
    `CalCore_m2_{HomoLumo}.cub` の位相を比較し、結果を文字列で返す。
    比較はランダムに選んだ20個以下の要素に対して行い、
    "Same" (同位相), "Opposit" (逆位相),
    または "not available" (判定不可) のいずれかを返す。

    Args:
        tcalpath (str): CUBEファイルが格納されているディレクトリのパス。
        CalCore (str): ファイル名の共通部分。
        HomoLumo (str): "HOMO" または "LUMO" を指定。

    Returns:
        str: 位相比較の結果。以下のいずれかの値を返す。
            - "Same": 同位相
            - "Opposit": 逆位相
            - "not available(1)": ファイルのグリッド数が異なる
            - "not available(2)": 位相が混在している
            - "not available(3)": 予期せぬエラーが発生
    """

    Val1 = ReadCube(tcalpath, f"{CalCore}_m1_{HomoLumo}.cub")
    Val2 = ReadCube(tcalpath, f"{CalCore}_m2_{HomoLumo}.cub")

    if len(Val1) != len(Val2):
        return "not available(1)"

    PClist = []
    printf(f"{CalCore} {HomoLumo} Data Excerpt")
    for i in random.sample(range(len(Val1)), min(20, len(Val1))):
        if abs(Val1[i]) > 1e-7:
            Val = Val1[i] * Val2[i]
            printf(f"m1 value = {Val1[i]}\tm2 value = {Val2[i]}")
            PClist.append("Opposit" if Val < 0 else "Same")

    PC = list(set(PClist))
    if len(PC) == 1:
        return PC[0]
    elif len(PC) > 1:
        return "not available(2)"
    else:
        return "not available(3)"


def ReadCube(tcalpath, Cubefile):
    """Cubeファイルを読み込み、値のリストを返す。

    指定されたパスにあるCubeファイルを読み込み、不要なコメント行とヘッダー情報をスキップした後、
    数値データ部分を抽出してリスト形式に変換します。読み込み後、Cubeファイルを削除します。

    Args:
        tcalpath (str): Cubeファイルが格納されているディレクトリのパス。
        Cubefile (str): 読み込むCubeファイルの名前。

    Returns:
        list: Cubeファイルから抽出された数値データのリスト。
    """
    with open(f"{tcalpath}/{Cubefile}", "r") as file:
        for _ in range(2):  # 最初の2行（コメント）をスキップ
            next(file)

        NumE, Xgrid, Ygrid, Zgrid = (int(next(file).split()[0]) for _ in range(4))
        for _ in range(NumE):  # NumEの数だけ行をスキップ
            next(file)

        values = [float(x) for _ in range(Xgrid * Ygrid)
                  for line in file for x in line.split()]

    subprocess.run(["rm", Cubefile], cwd=tcalpath, timeout=10)
    return values


def combineData(MaterName, Nmol, tcalpath, resultpath, MinFileName, TIFileName, PCFileName, Formated_Tilt, mol_pos):
    """複数のファイルからデータを読み込み、結合して保存する。

    MinFileName、TIFileName、PCFileNameで指定されたファイルからデータを読み込み、
    条件に応じて結合し、resultpathに保存する。

    Args:
        MaterName (str): マテリアル名。
        Nmol (str): 分子数 ("2mol" or "3mol")。
        tcalpath (str): 入力ファイルのディレクトリパス。
        resultpath (str): 出力ファイルのディレクトリパス。
        MinFileName (str): Minファイル名。
        TIFileName (str): TIファイル名。
        PCFileName (str): PCファイル名 (存在する場合)。
        Formated_Tilt (int):
        mol_pos (str):

    Raises:
        SystemExit:
            - MinFileNameとTIFileNameのデータ行数が一致しない場合。
            - Nmolが"2mol"でも"3mol"でもない場合。
    """
    with open(f"{tcalpath}/{TIFileName}", "r") as TIFile:
        TILines = TIFile.readlines()
        del TILines[0:2]

    with open(MinFileName, "r") as MinFile:
        MinLines = MinFile.readlines()
        del MinLines[0:2]

    if os.path.exists(f"{tcalpath}/{PCFileName}"):
        with open(f"{tcalpath}/{PCFileName}", "r") as PCFile:
            PCLines = PCFile.readlines()
            del PCLines[0]
    else:
        PCLines = []

    if len(TILines) == len(MinLines) and "2mol" in Nmol:
        pass
    elif len(TILines) == len(MinLines) * 3 and "3mol" in Nmol:
        pass
    else:
        printf(f"The numbers of data lines in {MinFileName} and {TIFileName} DO NOT match.")
        printf("Any file was NOT be changed.")
        exit()

    CombLines = []
    CombLines12 = []
    CombLines23 = []
    CombLines31 = []
    printf("Entry\tAngle\tDcol\tDtrv\tCpCE\tBSE\t**"
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
            printf("UNEXPECTED ERROR occurs in the combineData function.")
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

    CombLines.sort()
    CombLines12.sort()
    CombLines23.sort()
    CombLines31.sort()
    saveCombData(f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d",
                 f"{resultpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min-TIs.txt", CombLines)
    saveCombData(f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d-12",
                 f"{resultpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min-TIs-12.txt", CombLines12)
    saveCombData(f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d-23",
                 f"{resultpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min-TIs-23.txt", CombLines23)
    saveCombData(f"{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d-31",
                 f"{resultpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min-TIs-31.txt", CombLines31)
    printf(f"\n\n{MinFileName} and {TIFileName} were combined into "
           f"{resultpath}/{MaterName}_{Nmol}{mol_pos}_t{Formated_Tilt}d_min-TIs*.txt.")
    return


def correctTI(PhaseChk, TI):
    """位相情報に基づいて時間積分値を補正する。

    Args:
        PhaseChk (str): 位相情報 ("Same" or "Opposit")。
        TI (str): 補正前の時間積分値。

    Returns:
        float: 補正後の時間積分値。

    Raises:
        ValueError: PhaseChkが"Same"または"Opposit"でない場合。
    """
    TI = float(TI)
    if PhaseChk == "Same":
        CorrectTI = TI
    elif PhaseChk == "Opposit":
        CorrectTI = -1 * TI
    else:
        CorrectTI = TI
    return CorrectTI


def saveCombData(Name, path, List):
    """計算結果をファイルと標準出力に出力する。

    指定されたパスにファイルを作成し、リストの内容を書き込みます。
    リストが空の場合は何も出力しません。
    出力内容は、ヘッダー行とリストの各要素で構成されます。
    ヘッダー行は、計算結果の種類と各列の項目を示します。
    リストの各要素は、計算結果の各行に対応します。
    各行はタブ区切りで、各列は対応する項目の値を示します。

    Args:
        Name (str): 出力ファイルに書き込まれる計算結果の名前。
        path (str): 出力ファイルのパス。
        List (list): 書き込むデータのリスト。リストの各要素は文字列で、
            タブ区切りで各項目の値を表す。

    Returns:
        None
    """
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
            printf(f"\n{header}")
            for line in List:
                File.write(line)
                printf(line.strip())
    return


def copy_file(src, dest, lacks_list):
    """
    指定されたファイルをコピーする。

    ファイルが存在する場合は指定された宛先にコピーし、コピーが成功したことを
    コンソールに表示する。ファイルが存在しない場合は、ファイルパスを
    lacks_listに追加する。

    Parameters:
    src (str): コピー元のファイルパス。
    dest (str): コピー先のディレクトリパス。
    lacks_list (list): 存在しないファイルのパスを格納するリスト。

    Returns:
    None
    """
    if os.path.isfile(src):
        execute(["cp", src, dest], False)
        printf(f"\n{src} was copied into the {dest}")
    else:
        lacks_list.append(src)
    return


if __name__ == "__main__":
    main()
