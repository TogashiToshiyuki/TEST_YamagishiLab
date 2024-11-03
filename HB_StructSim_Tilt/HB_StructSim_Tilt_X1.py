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
from itertools import islice

"""
First writhing on Sat june 1 00:05:41
This program is writhing.

@author: Toshiyuki

dwotnisuru

"""


class Stereotyped:
    """定型文のまとめ"""
    ProgramAbst = """
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*   The aim of this programme is to check the structure for calculations at 2mol/3mol.    *
*   The programme requires the following arguments.                                       *
*   - xyz file of the molecule from which the calculations are made.                      *
*                                                                                         *
*   The following files are also required.                                                *
*   - AxSetting.txt                                                                       *
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"""
    AbnStop = """
*************** PROGRAM ABNORMALLY STOPPED ***************
    """
    sh_txt = ("#!/bin/sh\n"
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
    tcal_sh_txt = "\
    #!/bin/sh\n\
    \n\
    #$ -S /bin/sh\n\
    #$ -cwd\n\
    #$ -V\n\
    #$ -pe gau 12\n\
    #$ -q all.q\n\
    \n\
    tcal *.xyz\n\
    \n\
    g16 xxx.gjf\n\
    rm -rf /scr/$JOB_ID\n\
    \n"
    HelpText = """"""
    Header_2mol = ("%mem=24GB\n"
                   "%nprocshared=12\n"
                   "%chk=test_2mol.chk\n"
                   "#p pbepbe/6-31g(d) empiricaldispersion=gd3 counterpoise=2\n"
                   "\n"
                   "Title Card Required\n"
                   "\n"
                   "0 1\n")
    Header_2mol = ("%mem=24GB\n"
                   "%nprocshared=12\n"
                   "%chk=test_3mol.chk\n"
                   "#p pbepbe/6-31g(d) empiricaldispersion=gd3 counterpoise=3\n"
                   "\n"
                   "Title Card Required\n"
                   "\n"
                   "0 1\n")
    AxSetting_Header =("\n"
                       "The AxSetting.txt is not finded in the current directory.\n"
                       "The txt file should be described as follow.\n"
                       "*****\n"
                       "Column Direction: [x, y, or z]\n"
                       "Transverse Direction: [x, y, or z]\n"
                       "Coordinate of Dummy: x.xx, y.yy, z.zz\n"
                       "Rotate Axis: [x, y, or z]\n"
                       "Tilt Axis: [x, y, or z]\n")
    Energy_Header = ("Angle (degree)"
                     "\tDistance in column direction (Angstrom)"
                     "\tDistance in transverse direction (Angstrom)"
                     "\tCounterpoise corrected energy (A.U)"
                     "\tBSSE energy (A.U)\n")
    Combine_Header = ("Entry"
                      "\tAngle "
                      "\tD in col. (Å)"
                      "\tD in transv. (Å)"
                      "\tCounterpoise corrected energy (AU)"
                      "\tBSSE energy (AU)"
                      "\t**\tTI-NLUMO (meV)"
                      "\tTI-LUMO (meV)"
                      "\tTI-HOMO (meV)"
                      "\tTI-NHOMO (meV)\t**"
                      "\tPhaseChk (LUMO)"
                      "\tTI (LUMO)"
                      "\tPhaseChk (HOMO)"
                      "\tTI (HOMO)\n")


def main():
    Help = False
    Debug = False
    print(Stereotyped.ProgramAbst)
    messages = []
    HelpList = []
    MaterName = ""

    try:
        argument1 = sys.argv[1]
        if os.path.exists(argument1) and ".xyz" in argument1:
            MaterName = argument1[0:argument1.rfind(".xyz")]
        else:
            MaterName = ''
            Help = True
            print("""The file name is invalid in the program.
                Please check files.""")
    except IndexError:
        Help = True
    else:
        pass

    try:
        argument2 = sys.argv[2]
        if "Debug" in argument2 or "-d" in argument2 or "debug" in argument2:
            Debug = True
            print("""*************** Caution!!! Debug Started!!! ***************\n""")
        elif "Help" in argument2 or "-h" in argument2 or "help" in argument2:
            Help = True
        else:
            pass
    except IndexError:
        Debug = False
    else:
        pass

    CheckPoint(Help)

    Operator = input("\nEnter your name. >>> ")
    if Operator == "":
        Operator = "ONE"
    else:
        pass

    N = int(input("""\nPlease select the number of moles to be calculated from 2 and 3.
    >>> """))
    while N != 2 and N != 3:
        print("\nError: Incorrect input.")
        N = input("""\nPlease select the number of moles to be calculated from 2 and 3.
    Enter the number of moles to calculated.
    >>> """)
    else:
        pass
    Nmol = f"{N}mol"

    Tilt = int(input("""\nPlease enter the angle at which the molecule is tilted.
    >>> """))
    print(f"\n{Tilt} has been entered as Tilt.")

    print("\n\nFile set check...")
    if os.path.exists(f"./{MaterName}_{Nmol}_t{Tilt}d.gjf"):
        messages.append(f"{MaterName}_{Nmol}_t{Tilt}d.gjf: Found")
        HelpList.append(False)
    else:
        messages.append(f"{MaterName}_{Nmol}_t{Tilt}d.gjf: NOT Found")
        HelpList.append(True)

    if os.path.exists(f"./{MaterName}.xyz"):
        messages.append(f"{MaterName}.xyz: Found")
        HelpList.append(False)
    else:
        messages.append(f"{MaterName}.xyz: NOT Found")
        HelpList.append(True)

    if os.path.exists(f"./ConditionList_Tilt_{Nmol}_t{Tilt}d.txt"):
        messages.append(f"./ConditionList_Tilt_{Nmol}_t{Tilt}d.txt: Found")
        HelpList.append(False)
    elif os.path.exists(f"./InitialCondition_Tilt_{Nmol}_t{Tilt}d.txt"):
        messages.append(f"./InitialCondition_Tilt_{Nmol}_t{Tilt}d.txt: Found")
        HelpList.append(False)
    else:
        messages.append(
            f"./ConditionList_Tilt_{Nmol}_t{Tilt}d.txt and ./InitialCondition_Tilt_{Nmol}_t{Tilt}d.txt: Not Found")
        HelpList.append(True)

    if "3mol" in Nmol:
        if os.path.exists(f"./{MaterName}_2mol_t{Tilt}d_min.txt"):
            messages.append(f"{MaterName}_2mol_t{Tilt}d_min.txt: Found")
            HelpList.append(False)
        else:
            messages.append(f"{MaterName}_2mol_t{Tilt}d_min.txt: Not Found")
            HelpList.append(True)
    elif "2mol" in Nmol:
        HelpList.append(False)
    else:
        HelpList.append(True)
        messages.append("File name of the gjf is not correct.")

    if True in HelpList:
        Help = True
        print("Required file set DOES NOT exist in the correct directory.")
    else:
        print("Required file set EXISTs in the correct directory.")

    printList(messages)
    CheckPoint(Help)

    which = ""
    RefLines = []
    if N == 2:
        which = "Dcol"
        RefLines = []
    elif N == 3:
        which = "Dtrv"
        RefLines = getRefLines(f"./{MaterName}_2mol_t{Tilt}d_min.txt")
    else:
        Help = True

    DebugCheck([N, which], Debug, "none", False)
    CheckPoint(Help)

    if os.path.exists(f"./ConditionList_Tilt_{Nmol}_t{Tilt}d.txt"):
        with open(f"./ConditionList_Tilt_{Nmol}_t{Tilt}d.txt") as file:
            lines = file.readline()
        split_lines = lines.split('-')
        Current_Tilt = float(split_lines[1].replace('d', '').replace('t', ''))
        if not Current_Tilt == Tilt:
            print(f"""*************** Caution !!! ***************
    The current angle and the entered angle are different !!!

    Input Tilt Angle   : {Tilt}
    Current Tilt Angle : {Current_Tilt}

    The process proceeds based on the angle *Entered* !!!
    *************** Caution !!! ***************

    Processing will resume after 3 seconds.
            """)
            time.sleep(3)
        else:
            pass
    else:
        Current_Tilt = 10000

    if os.path.exists(f"./ConditionList_Tilt_{Nmol}_t{Tilt}d.txt") and Current_Tilt == Tilt:
        pass
    else:
        with open(f"./InitialCondition_Tilt_{Nmol}_t{Tilt}d.txt", "r") as file:
            lines = file.readlines()

        NewConditions = []
        for line in lines:
            Deg = float(line.strip().split()[0])
            Val = float(line.strip().split()[1])
            # 3軸目方向に傾ける角度は定義済み、Tiltで呼び出せる
            RefValues = getRefValues(Deg, RefLines)

            NewCondition = mkNewCondition(Nmol, Deg, Val, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val + 0.1, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val - 0.1, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            """
            NewCondition = mkNewCondition(Nmol, Deg, Val + 0.2, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val - 0.2, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            """
        NewConditions = list(set(NewConditions))
        NewConditions.sort()
        file = open(f"./ConditionList_Tilt_{Nmol}_t{Tilt}d.txt", "w")
        for Condition in NewConditions:
            file.write(Condition)
        file.close()
        print(f"\nNew ConditionList_Tilt_{Nmol}_t{Tilt}d.txt was written "
              f"from InitialCondition_Tilt_{Nmol}_t{Tilt}d.txt!!\n")
    dirpath = f"./{MaterName}_{Nmol}_t{Tilt}d"
    tcalpath = f"./{MaterName}_{Nmol}_t{Tilt}d_tcal"

    os.makedirs(dirpath, exist_ok=True)
    if N == 2:
        getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dcol", [], 0.1, Tilt, Operator)
        print("\nCalculations for 2mol were successfully finished.")
    elif N == 3:
        RefLines = getRefLines(f"./{MaterName}_2mol_t{Tilt}d_min.txt")
        getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dtrv", RefLines, 0.1, Tilt, Operator)
    else:
        pass

    temp_Structures = []
    if "3mol" in Nmol and not os.path.exists(f"./{MaterName}_{Nmol}_t{Tilt}d_mins.hist"):
        MostStable = False
        while not MostStable:
            RefLines = getRefLines(f"./{MaterName}_{Nmol}_t{Tilt}d_min.txt")
            temp_Structures.append(RefLines)
            mkCycleConditions(RefLines, "Dcol", 0.1, Nmol, Tilt)
            getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dcol", RefLines, 0.1, Tilt, Operator)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}_t{Tilt}d_min.txt")
            mkCycleConditions(RefLines, "Dtrv", 0.1, Nmol, Tilt)
            getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dtrv", RefLines, 0.1, Tilt, Operator)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}_t{Tilt}d_min.txt")
            MostStable = CompareStructures(RefLines, temp_Structures[-1])

        MostStable = False
        while not MostStable:
            RefLines = getRefLines(f"./{MaterName}_{Nmol}_t{Tilt}d_min.txt")
            temp_Structures.append(RefLines)
            mkCycleConditions(RefLines, "Dcol", 0.05, Nmol, Tilt)
            getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dcol", RefLines, 0.05, Tilt, Operator)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}_t{Tilt}d_min.txt")
            mkCycleConditions(RefLines, "Dtrv", 0.05, Nmol, Tilt)
            getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dtrv", RefLines, 0.05, Tilt, Operator)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}_t{Tilt}d_min.txt")
            MostStable = CompareStructures(RefLines, temp_Structures[-1])

        MinConditions = getMinConditions(MaterName, Nmol, Tilt)
        os.makedirs(tcalpath, exist_ok=True)

        for Condition in MinConditions:
            command = ["cp", f"{dirpath}/{MaterName}_{Nmol}_{Condition}.gjf",
                       f"{tcalpath}/{MaterName}_{Nmol}_{Condition}.gjf"]
            subprocess.run(command)

        if len(temp_Structures) < 100:
            file = open(f"./{MaterName}_{Nmol}_mins.hist", "w")
            for i in range(len(temp_Structures)):
                file.write(f"***** Structures after the '{i + 1}'th cycle *****")
                Lines = temp_Structures[i]
                for Line in Lines:
                    file.write(Line)
            file.close()
        else:
            pass

    print("\n**********\nMaking Resulting Data Set...\n")
    resultname = f"{MaterName}_{Nmol}_t{Tilt}d"
    print(f"For {resultname}")
    resultpath = f"./{resultname}_results"
    os.makedirs(resultpath, exist_ok=True)

    MinFileName = f"{resultname}_min.txt"
    TIFileName = f"{resultname}_TIs.txt"
    AllFileName = f"{resultname}_all.txt"
    PCFileName = "PhaseCheck.txt"

    Lacks = []

    # transfer of MinFile
    if os.path.isfile(MinFileName) and os.path.isfile(f"{tcalpath}/{TIFileName}"):
        combineData(MaterName, Nmol, tcalpath, resultpath, MinFileName, TIFileName, PCFileName)
    elif os.path.isfile(MinFileName):
        subprocess.run(["cp", MinFileName, f"{resultpath}/{MinFileName}"])
        print(f"{MinFileName} was copied into the {resultpath} forder. ")
    else:
        Lacks.append(MinFileName)

    # transfer of AllFile
    if os.path.isfile(AllFileName):
        subprocess.run(["cp", AllFileName, f"{resultpath}/{AllFileName}"])
        print(f"{AllFileName} was copied into the {resultpath} forder. ")
    else:
        Lacks.append(AllFileName)

    # transfer of XYZ file
    if os.path.isfile(f"{MaterName}.xyz"):
        subprocess.run(["cp", f"{MaterName}.xyz", f"{resultpath}/{MaterName}.xyz"])
        print(f"{MaterName}.xyz was copied into the {resultpath} forder. ")
    else:
        Lacks.append(f"{MaterName}.xyz")

    # transfer of ConditionList
    if os.path.isfile(f"ConditionList_Tilt_{Nmol}_t{Tilt}d.txt"):
        subprocess.run(
            ["cp", f"ConditionList_Tilt_{Nmol}_t{Tilt}d.txt", f"{resultpath}/ConditionList_Tilt_{Nmol}_t{Tilt}d.txt"])
        print(f"ConditionList_Tilt_{Nmol}_t{Tilt}d.txt was copied into the {resultpath} forder. ")
    else:
        Lacks.append(f"ConditionList_Tilt_{Nmol}_t{Tilt}d.txt")

    # transfer of gjf
    if os.path.isfile(f"{resultname}.gjf"):
        subprocess.run(["cp", f"{resultname}.gjf", f"{resultpath}/{resultname}.gjf"])
        print(f"{resultname}.gjf was copied into the {resultpath} forder. ")
    else:
        Lacks.append(f"{resultname}.gjf")

    # transfer of AxSetting
    if os.path.isfile("../Writhing/HB_StructSim_Tilt/AxSetting.txt"):
        subprocess.run(["cp", "AxSetting.txt", f"{resultpath}/AxSetting.txt"])
        print(f"AxSetting.txt was copied into the {resultpath} forder. ")
    else:
        Lacks.append("AxSetting.txt")

    # transfer of mins.hist
    if os.path.isfile(f"{resultname}_mins.hist") and "3mol" in Nmol:
        subprocess.run(
            ["cp", f"{resultname}_mins.hist", f"{resultpath}/{resultname}_mins.hist"])
        print(f"{resultname}_mins.hist was copied into the {resultpath} forder. ")
    else:
        if "3mol" in Nmol:
            Lacks.append(f"{resultname}_mins.hist")

    # transfer of structure files
    StructFiles = []
    if "2mol" in Nmol:
        StructFiles = glob.glob(f"{tcalpath}/*.xyz")
    elif "3mol" in Nmol:
        StructFiles = glob.glob(f"{tcalpath}/*.all")

    if len(StructFiles) == 0:
        print("Any file for aggregattion structure is not found in the specific forder.")
        Lacks.append("Stractural xyz files")
    else:
        print(f"xyz files were copied into the {resultpath} folder.")
        for File in StructFiles:
            name = File[File.rfind("/") + 1:-4]
            subprocess.run(["cp", File, f"{resultpath}/{name}.xyz"])
            print(f"\t{name}.xyz")

    if len(Lacks) != 0:
        print("\n\n***********\nSeveral result files were not found in the specific folder.")
        for content in Lacks:
            print(f"\t{content}")
    else:
        print(f"The {MaterName}_{Nmol} data set were collected in {resultpath}.")

    print("\n************************* ALL PROCESSES END *************************\n")


# 関数
def CheckPoint(HelpCheck):
    """
        ヘルプフラグをチェックし、ヘルプが必要な場合はヘルプメッセージを表示してプログラムを終了する関数

        Args:
            HelpCheck (bool): ヘルプが必要かどうかを示すフラグ

        Returns:
            None

        Raises:
            SystemExit: ヘルプフラグがTrueの場合、プログラムを終了する

        Note:
            - この関数は、プログラム開始時にヘルプオプションが指定された場合などに使用することを想定しています。
            - ヘルプメッセージの内容は、Stereotyped.HelpTextに定義されています。
            - プログラム終了時のメッセージは、Stereotyped.AbnStopに定義されています。
            - ヘルプフラグがFalseの場合は、何もしません。
        """
    if HelpCheck:
        print(Stereotyped.HelpText)
        print(Stereotyped.AbnStop)
        exit()
    else:
        pass


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
        try:
            RefLines.append(line)
        except ValueError:
            pass
    return RefLines


def getRefValues(Deg, RefLines):
    RefValues = ''
    if len(RefLines) == 0:
        RefValues = ["na", "na"]
    else:
        for RefLine in RefLines:
            Contents = RefLine.strip().split()
            try:
                RefDeg = round(float(Contents[0]), 1)
                if RefDeg == Deg:
                    RefValues = [Contents[1], Contents[2]]
                else:
                    pass
            except (ValueError, IndexError):
                pass

    return RefValues


def mkNewCondition(Nmol, Deg, Val, Tilt, which, RefValues):
    """新しく生成された条件を文字列として返す

    Args:
        Nmol (str): 分子数 ("2mol" or "3mol")
        Deg (float): 回転角度 (度)
        Val (float): DcolまたはDtrvの値
        Tilt (float): 傾斜角度 (度)
        which (str): ValがDcolかDtrvかを示す ("Dcol" or "Dtrv")
        RefValues (list): [Dcol, Dtrv] の参照値

    Returns:
        str: 生成された条件を表す文字列 (例: "30d-t15d-85-20")

    Raises:
        ValueError: 不正な引数が渡された場合

    Note:
        * Nmolが"2mol"の場合、ValはDcolとして扱われ、Dtrvは無視されます。
        * Nmolが"3mol"の場合、whichで指定された値がValとして使用され、もう一方の値はRefValuesから取得されます。
        * 生成された条件は標準出力にも表示されます。
    """
    NewCondition = ""
    Dcol = 0
    Dtrv = 0
    if "2mol" in Nmol:
        Dcol = int(round(Val * 100, 2))
        Deg = int(Deg)
        Tilt = int(Tilt)
        NewCondition = f"{Deg}d-t{Tilt}d-{Dcol}\n"
    elif "3mol" in Nmol:
        Deg = int(Deg)
        Tilt = int(Tilt)
        if which == "Dcol":
            Dcol = int(round(Val * 100, 2))
            Dtrv = int(round(float(RefValues[1]) * 100, 2))
        elif which == "Dtrv":
            Dtrv = int(round(Val * 100, 2))
            Dcol = int(round(float(RefValues[0]) * 100, 2))
        NewCondition = f"{Deg}d-t{Tilt}d-{Dcol}-{Dtrv}\n"
    print(f"\t\t{NewCondition.strip()}")
    return NewCondition


def printList(List):
    for content in List:
        print(content)
    print("\n")
    return


def DebugCheck(DebugList, Flag, memo, Exit):
    if Flag:
        bugs = int(len(DebugList))
        print("\n*************** Check flag passed ***************")
        if memo == "none":
            pass
        else:
            print(f"Memo : {memo}")
        for bug in range(bugs):
            print(DebugList[bug])
        if Exit:
            exit()
        else:
            pass
    else:
        pass
    return


def execute_command(command):
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    if result.returncode == 0:
        return result.stdout.strip()
    else:
        command = " ".join(command)
        print(f"Error executing command: {command}")
        print(f"Error message: {result.stderr.strip()}")
        return None


def getJobID():
    output = execute_command(["qstat"])
    lines = output.splitlines()

    if len(lines) == 0:
        firstJobID = "None"
        lastJobID = "None"
    else:
        firstJobID = lines[2].strip().split()[0]
        lastJobID = lines[-1].strip().split()[0]

    JobID = [firstJobID, lastJobID]
    return JobID


def getElapsedTime(start_time):
    now = datetime.datetime.now()
    elapsed_time = now - start_time
    elapsed_time = round(elapsed_time.total_seconds() / 60, 1)
    formated_NOW = now.strftime("%m/%d %H:%M")
    return formated_NOW, elapsed_time


def rmWildCards(wildcard):
    lines = glob.glob(wildcard)
    for line in lines:
        subprocess.run(["rm", line], timeout=10)
    return


def getMinConditions(MaterName, Nmol, Tilt):
    with open(f"./{MaterName}_{Nmol}_{Tilt}_min.txt") as file:
        lines = file.readlines()
    del lines[0:2]

    MinConditions = []
    for line in lines:
        contents = line.strip().split()
        Deg = int(float(contents[0]))
        Tilt = int(float(contents[1]))
        Dcol = int(round(float(contents[2]) * 100, 5))
        if contents[3] == "na":
            Condition = f"{Deg}d-t{Tilt}d-{Dcol}"
        else:
            Dtrsv = int(round(float(contents[3]) * 100, 5))
            Condition = f"{Deg}d-t{Tilt}d-{Dcol}-{Dtrsv}"
        MinConditions.append(Condition)
    return MinConditions


def getCondition_fromName(Name):
    Condition = Name[Name.rfind("_") + 1:Name.rfind(".")]
    return Condition


def getVALfromLogName(Nmol, Log):
    Condition = getCondition_fromName(Log)
    FileName = Log[Log.rfind("/") + 1:-4]

    Conditions = Condition.split('-')

    if "2mol" in Nmol and Condition.count("-") == 2:
        deg = Conditions[0].replace('d', '')
        tilt = Conditions[1].replace('d', '')
        dcol = Conditions[2]
        dtrv = "na"
    elif "3mol" in Nmol and Condition.count("-") == 3:
        deg = Conditions[0].replace('d', '')
        tilt = Conditions[1].replace('d', '')
        dcol = Conditions[2]
        dtrv = Conditions[3]
    else:
        print("\nUNEXPECTED ERROR happens in the function of getVALfromLogName!!")
        exit()

    Vdeg = round(float(deg), 5)
    Vdcol = round(float(dcol) / 100, 5)
    VTilt = round(int(tilt), 5)
    try:  # Nmol=="3mol"
        Vdtrv = round(float(dtrv) / 100, 5)
    except ValueError:  # Nmol=="2mol"
        Vdtrv = dtrv
    return FileName, Vdeg, Vdcol, Vdtrv, VTilt


def getEnergy(data):
    CPE = data.split("Counterpoise corrected energy =")[1]
    CPE = CPE.split("BSSE energy")[0]
    CPE = CPE.strip()
    CPE = float(CPE)
    BSE = data.split("BSSE energy =")[1]
    BSE = BSE.split("sum of fragments")[0]
    BSE = BSE.strip()
    BSE = float(BSE)
    return CPE, BSE


def readEnergies(dir_path, MaterName, Nmol, Tilt):
    LogList = glob.glob(f"./{dir_path}/{MaterName}_{Nmol}_*.log")
    LogList.sort()
    DegList = []
    for Log in LogList:
        FileName, deg, Vdcol, Vdtrv, Tilt = getVALfromLogName(Nmol, Log)
        DegList.append(deg)
    DegList = list(set(DegList))
    DegList.sort()

    with open(f"./{MaterName}_{Nmol}_t{Tilt}d_all.txt", "w") as AllData, open(f"./{MaterName}_{Nmol}_t{Tilt}d_min.txt",
                                                                              "w") as MinData:
        AllData.write(f"*****  {MaterName}_{Nmol} All Results *****\n")
        MinData.write(f"***** {MaterName}_{Nmol} Minimum Energy at each Angle *****\n")
        AllData.write(f" \t{Stereotyped.Energy_Header}")
        MinData.write(Stereotyped.Energy_Header)
        MinimumL = []

        print(f" \t{Stereotyped.Energy_Header}")
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
                        FileName, Vdeg, Vdcol, Vdtrv, VTilt = getVALfromLogName(Nmol, Log)
                        CPE, BSE = getEnergy(data)
                        CPE_Dict[FileName] = CPE
                        VAL_Dict[FileName] = f"{Vdeg}\t{Vdcol}\t{Vdtrv}\t{CPE}\t{BSE}\n"
                        Keys.append(FileName)

                    else:
                        print(f"\t{Log} was NOT normally terminated. It was removed.")
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


def mkNewConditionList(MaterName, Nmol, VTilt, which, dev, RefLines):
    """
    特定条件下での局所的最小値探索と、新たな探索条件リストの作成

    与えられたファイル (MaterName_Nmol_tVTiltd_all.txt) から、指定された条件 (Nmol, VTilt, which) に基づいて、
    各角度における局所的最小値を探索します。探索範囲は指定された基準値 (RefLines) からの偏差 (dev) で決定されます。
    局所的最小値が見つかった場合は "complete" と判定し、見つからなかった場合は "not complete" と判定します。
    また、"not complete" と判定された角度に対しては、新たな探索条件を生成し、既存の条件リスト (ConditionList_Tilt_Nmol_tVTiltd.txt) に追加します。
    最後に、探索結果の概要 (局所的最小値が見つかった角度の数と値) を出力します。

    Args:
        MaterName (str): 解析対象の物質名
        Nmol (str): 分子数 ("2mol" or "3mol")
        VTilt (float): 傾斜角
        which (str): 解析対象の値 ("Dcol" or "Dtrv")
        dev (float): 基準値からの偏差
        RefLines (list): 各角度における基準値のリスト

    Returns:
        bool: 全ての角度で局所的最小値が見つかった場合は True、そうでない場合は False
    """
    n = 1
    file_path = f"./{MaterName}_{Nmol}_t{VTilt}d_all.txt"
    DegList = set()

    with open(file_path, "r") as All:
        for line in All:
            try:
                Deg = round(float(line.strip().split("\t")[1]), 1)
                DegList.add(Deg)
            except (IndexError, ValueError):
                continue

    DegList = sorted(DegList)

    Judges = []
    ComplDeg = []
    ComplVal = []
    NewConditions = []

    for Deg in DegList:
        RefValues = getRefValues(Deg, RefLines)
        ValueList = []
        SV = None

        if "3mol" in Nmol:
            SV = round(float(RefValues[1 if which == "Dcol" else 0]), 2)

        with open(file_path, "r") as All:
            for line in All:
                try:
                    Contents = line.strip().split("\t")
                    DataDeg = round(float(Contents[1]), 1)

                    if DataDeg == Deg:
                        if "2mol" in Nmol:
                            DataDcol = round(float(Contents[2]), 2)
                            ValueList.append(DataDcol)
                            if Contents[0] == "*":
                                SV = DataDcol
                        elif "3mol" in Nmol:
                            DataDcol = round(float(Contents[2]), 2)
                            DataDtrv = round(float(Contents[3]), 2)
                            if which == "Dcol":
                                RefDtrv = round(float(RefValues[1]), 2)
                                if DataDtrv == RefDtrv:
                                    ValueList.append(DataDcol)
                                if Contents[0] == "*":
                                    SV = DataDcol
                            elif which == "Dtrv":
                                RefDcol = round(float(RefValues[0]), 2)
                                if DataDcol == RefDcol:
                                    ValueList.append(DataDtrv)
                                if Contents[0] == "*":
                                    SV = DataDtrv
                except (ValueError, IndexError):
                    continue

        if SV is None:
            continue

        if round(SV + dev, 2) in ValueList and round(SV - dev, 2) in ValueList:
            Judges.append("complete")
            print(f"\nThe local minimum by {dev} step for {Deg} degree was Found;\t{SV}.")
            print("\tCOMPLETE!!")
            ComplDeg.append(Deg)
            ComplVal.append(SV)
        elif round(SV + 0.1, 2) in ValueList and round(SV - 0.1, 2) in ValueList:
            Judges.append("not complete")
            print(f"\nThe local minimum by 0.1 step for {Deg} degree was FOUND;\t{SV}.")
            print("\tNOT COMPLETE(1)!! 2 new conditions bellow were appended to the ConditionList.txt.")
            ComplDeg.append(Deg)
            ComplVal.append(SV)
            NewConditions.append(mkNewCondition(Nmol, Deg, SV - 0.05, VTilt, which, RefValues))
            NewConditions.append(mkNewCondition(Nmol, Deg, SV + 0.05, VTilt, which, RefValues))
        elif SV == min(ValueList):
            Judges.append("not complete")
            print(f"\nLocal minimum for {Deg} degree was NOT FOUND in the cycle.")
            print(f"\tNOT COMPLETE(2)!! {n} new conditions bellow were appended to the ConditionList.txt.")
            for i in range(n):
                NewConditions.append(mkNewCondition(Nmol, Deg, SV - dev * (i + 1), VTilt, which, RefValues))
        elif SV == max(ValueList):
            Judges.append("not complete")
            print(f"\nLocal minimum for {Deg} degree was NOT FOUND in the cycle.")
            print(f"\tNOT COMPLETE(3)!! {n} new conditions bellow were appended to the ConditionList.txt.")
            for i in range(n):
                NewConditions.append(mkNewCondition(Nmol, Deg, SV + dev * (i + 1), VTilt, which, RefValues))
        else:
            Judges.append("not complete")
            print(f"\nLocal minimum for {Deg} degree was NOT FOUND in the cycle.")
            print("\tNOT COMPLETE(4)!! 2 new conditions bellow were appended to the ConditionList.txt.")
            ComplDeg.append(Deg)
            ComplVal.append(SV)
            NewConditions.append(mkNewCondition(Nmol, Deg, SV - 0.1, VTilt, which, RefValues))
            NewConditions.append(mkNewCondition(Nmol, Deg, SV + 0.1, VTilt, which, RefValues))

    with open(f"ConditionList_Tilt_{Nmol}_t{VTilt}d.txt", "r") as file:
        orgCondition = file.readlines()

    NewList = sorted(set(orgCondition + NewConditions))

    with open(f"ConditionList_Tilt_{Nmol}_t{VTilt}d.txt", "w") as file:
        file.writelines(NewList)

    if len(ComplDeg) == 0:
        print("\nLocal minimum have NOT been FOUND in any angle.")
    else:
        print(f"\nLocal minimum was FOUND in {len(ComplDeg)} angels.")
        if "2mol" in Nmol:
            print("Angel\tD_col")
        elif "3mol" in Nmol:
            print("Angel\tD_trv")
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


def getTemporaryStructure(Nmol, MaterName, dirpath, Debug, which, RefLines, dev, Tilt, Operator):
    judge = False

    # テスト用の生成用ファイル
    Conditions = []
    with open(f"ConditionList_Tilt_{Nmol}_t{Tilt}d.txt") as ConditionFile:
        lines = ConditionFile.readlines()
    for line in lines:
        Conditions.append(line.strip())
    qsub_List = []
    with open(f"./{dirpath}/G.sh", "w") as original_sh:
        original_sh.write(Stereotyped.sh_txt)

    for Condition in Conditions:
        if os.path.exists(f"./{dirpath}/{MaterName}_{Nmol}_{Condition}.log"):
            pass
        else:
            qsub_temp = mkRotate_Structure_File(MaterName, Nmol, Debug, dirpath, Condition, Operator)
            qsub_List.append(qsub_temp)

    '''DebugCheck(qsub_List,Debug,"qsublist", False)'''
    """
    while not judge:
        Conditions = []
        with open(f"ConditionList_Tilt_{Nmol}_t{Tilt}d.txt") as ConditionFile:
            lines = ConditionFile.readlines()
        for line in lines:
            Conditions.append(line.strip())
        qsub_List = []
        with open(f"./{dirpath}/G.sh", "w") as original_sh:
            original_sh.write(Stereotyped.sh_txt)

        for Condition in Conditions:
            if os.path.exists(f"./{dirpath}/{MaterName}_{Nmol}_{Condition}.log"):
                pass
            else:
                qsub_temp = mkRotate_Structure_File(MaterName, Nmol, Debug, dirpath, Condition, Operator)
                qsub_List.append(qsub_temp)
        JobID = ""
        print("\n**********\nJobs are submitting...")
        if len(qsub_List) == 0:
            untilID = "None"
        else:
            for qsub in qsub_List:
                qsub = qsub.split()
                subprocess.run(qsub, cwd=f"./{dirpath}")
            JobID = getJobID()
            untilID = int(JobID[1])
        os.remove(f"./{dirpath}/G.sh")

        if len(qsub_List) == 0:
            print("Any job was not submitted. Calculations with the conditions might be finished.")
        else:
            CalcEnd = False
            start_time = datetime.datetime.now()
            formated_ST = start_time.strftime("%m/%d %H:%M")
            term = ""
            if which == "Dcol":
                term = "the stable distance in column direction"
            elif which == "Dtrv":
                term = "the stable distance in transverse direction"
            print(f"\n'{len(qsub_List)}' calculations for '{term}' was submitted!! at {formated_ST}")
            if "2mol" in Nmol:
                Wait_minutes = 2
            elif "3mol" in Nmol:
                Wait_minutes = 5
            else:
                Wait_minutes = 2
            print(f"\nWait until jobID {untilID}!!")
            print(f"wait for {(int(untilID) - int(JobID[0]) + 1) * Wait_minutes} minute!")
            print(f"start ID: {JobID[0]}")

            while not CalcEnd:
                JobID = getJobID()

                try:
                    Remaining_calculations = (int(untilID) - int(JobID[0])) + 1
                    if Remaining_calculations < 1:
                        Remaining_calculations = 1
                    else:
                        pass
                except ValueError:
                    print("終了したっぽい")
                    Remaining_calculations = 0

                if "2mol" in Nmol:
                    Wait_minutes = 1
                elif "3mol" in Nmol:
                    Wait_minutes = 5
                else:
                    Wait_minutes = 1
                time.sleep(Wait_minutes * Remaining_calculations * 60)
                formated_NOW, elapsed_time = getElapsedTime(start_time)
                sys.stdout.write(
                    "\033[2k\033[G%s" % f"\t{formated_NOW} ({elapsed_time} min. passed): Job {JobID[0]} is in progress."
                                        f"\n\tThe remaining jobs are {Remaining_calculations}")
                sys.stdout.flush()

                if JobID[0] == "None" and JobID[1] == "None":
                    CalcEnd = True
                elif int(untilID) < int(JobID[0]):
                    CalcEnd = True
                elif int(untilID) >= int(JobID[0]):
                    CalcEnd = False
                else:
                    print("UNEXPECTED ERROR occured in main.")
                    print(Stereotyped.AbnStop)
                    exit()

                if not CalcEnd:
                    pass
                else:
                    print(f"Calculation cycles for {which} until JobID {untilID} were finished.")
                    rmWildCards(f"{dirpath}/*.sh*")
                    if Debug:
                        pass
                    else:
                        rmWildCards(f"{dirpath}/*.chk")
        print("\n**********\nReading Data...\n")
        readEnergies(dirpath, MaterName, Nmol, Tilt)
        judge = mkNewConditionList(MaterName, Nmol, int(Tilt), which, dev, RefLines)"""
    return


def Axis_Setting_HB():
    try:
        with open("../Writhing/HB_StructSim_Tilt/AxSetting.txt", "r") as File:
            lines = File.readlines()
    except FileNotFoundError:
        print(f"{Stereotyped.AxSetting_Header}")
        with open("../Writhing/HB_StructSim_Tilt/AxSetting.txt", "w") as file:
            file.write("Column Direction:  \n")
            file.write("Transverse Direction: \n")
            file.write("Coordinate of Dummy: x.xx,y.yy,z.zz\n")
            file.write("Rotate Axis: [x, y, or z]\n")
            file.write("Tilt Axis: [x, y, or z]\n")
        print(Stereotyped.AbnStop)
        exit()
    Direction_Col, Direction_Transv, Rotate_Axis, Tilt_Axis = 0, 0, 0, 0
    Position_Dummy = []
    for i in range(len(lines)):
        line = lines[i]
        content = line[:line.find(":")]
        if "column direction" in content.lower() and i < 5:
            Direction_Col = line[line.rfind(":") + 1:].strip()
        elif "transverse direction" in content.lower() and i < 5:
            Direction_Transv = line[line.rfind(":") + 1:].strip()
        elif "coordinate of dummy" in content.lower() and i < 5:
            XYZ = line[line.rfind(":") + 1:].strip()
            XYZ = list(np.float_(XYZ.split(",")))
            Position_Dummy = np.array(XYZ)
        elif "Rotate Axis" in content.lower() and i < 5:
            Rotate_Axis = line[line.rfind(":") + 1:].strip()
        elif "molecular2 angles" in content.lower() and i < 5:
            Tilt_Axis = line[line.rfind(":") + 1:].strip()
        else:
            Rotate_Axis = "x"
            Tilt_Axis = "z"
            pass

    if Direction_Col == Direction_Transv:
        print("INVALID PARAMETERS!! Column direction and Transverse direction should be different!!")
        print(Stereotyped.AbnStop)
        exit()

    return Direction_Col, Direction_Transv, Position_Dummy, Rotate_Axis, Tilt_Axis


def mkDirection(axis_direction, input_axis, axis_name):  # ダミー原子用の繰り返し処理を関数化
    if input_axis == "x":
        matrix = np.array([axis_direction, 0.0, 0.0])
    elif input_axis == "y":
        matrix = np.array([0.0, axis_direction, 0.0])
    elif input_axis == "z":
        matrix = np.array([0.0, 0.0, axis_direction])
    else:
        print(f"\nThe axis of {axis_name} is INVALID.")
        print(Stereotyped.AbnStop)
        exit()
    return matrix


def mkDummysPos_HB(Direction_Col, Direction_Transv, Position_Dummy, Dt1, Dc1):  # ダミー原子の座標と変換行列を返す
    if Direction_Col == Direction_Transv:  # 同じ向きに移動しようとしていないか確認
        print("""
        The axis of transverse direction is same with the axis of column direction.""")
        print(Stereotyped.AbnStop)
        exit()

    # Creation of Transition Matrix
    col_transl = mkDirection(Dc1, Direction_Col, "column direction")
    transv_transl = mkDirection(Dt1, Direction_Transv, "transverse direction")
    Direction_Other = "xyz"
    Direction_Other = Direction_Other.replace(Direction_Col, "")
    Direction_Other = Direction_Other.replace(Direction_Transv, "")
    a_transl = mkDirection(10, Direction_Other, "other direction")

    # Creation dummy atoms
    X1_pos = Position_Dummy + col_transl
    X2_pos = X1_pos + transv_transl
    X3_posc = Position_Dummy + col_transl / 2
    X3_pos = Position_Dummy + col_transl
    X4_pos = X3_pos + a_transl
    X5_pos = X3_posc + transv_transl
    X6_pos = Position_Dummy

    Dummy_pos = [X1_pos, X2_pos, X3_pos, X4_pos, X5_pos, X6_pos]
    Translations = [col_transl, transv_transl, a_transl]
    return Dummy_pos, Translations


def Angle(Current, Tx, Ty, Tz, rotation):
    Tx = math.radians(Tx)
    Ty = math.radians(Ty)
    Tz = math.radians(Tz)
    Rx = np.array([[1, 0, 0],
                   [0, math.cos(Tx), -math.sin(Tx)],
                   [0, math.sin(Tx), math.cos(Tx)]])

    Ry = np.array([[math.cos(Ty), 0, math.sin(Ty)],
                   [0, 1, 0],
                   [-math.sin(Ty), 0, math.cos(Ty)]])

    Rz = np.array([[math.cos(Tz), -math.sin(Tz), 0],
                   [math.sin(Tz), math.cos(Tz), 0],
                   [0, 0, 1]])
    rotation_matrices = {'x': Rx, 'y': Ry, 'z': Rz}
    R = np.eye(3)
    for axis in rotation:
        R = np.dot(rotation_matrices[axis], R)
    rotated_coords = np.dot(R, Current)
    return rotated_coords


def mkAtomList(Translations, MaterName, Angles, mol2_Angles, mol3_Angles, rotate):
    with open(f"./{MaterName}.xyz", "r") as f:
        NinMol = int(f.readline())
        AtomList = list(islice(f, 1, None))

    Mol1, Mol2, Mol3, Elist = [], [], [], []
    for Atom in AtomList:
        Atom = Atom.strip()
        Contents = Atom.split(None)
        Elist.append(Contents.pop(0))
        Position = list(np.float_(Contents))
        primitive_molecule = np.array(Position)

        atm_m1 = Angle(primitive_molecule, Angles[0], Angles[1], Angles[2], rotate)
        Mol1.append(atm_m1)

        atm_m2 = Angle(primitive_molecule, mol2_Angles[0], mol2_Angles[1], mol2_Angles[2], rotate)
        atm_m2 = atm_m2 + Translations[0]
        Mol2.append(atm_m2)

        atm_m3 = Angle(primitive_molecule, mol3_Angles[0], mol3_Angles[1], mol3_Angles[2], rotate)
        atm_m3 = atm_m3 + Translations[0] / 2 + Translations[1]
        Mol3.append(atm_m3)

    return Elist, Mol1, Mol2, Mol3, NinMol


def format_coordinate(coord, Number):
    formatted_coord = f"{coord[0]: 15.10f}   {coord[1]: 15.10f}   {coord[2]: 15.10f}  {Number}"
    return formatted_coord


def mkRotate_Structure_File(mk_Mater_Name, Nmol, Debug, dir_path, Condition, Operator):
    File_Name = f"{mk_Mater_Name}_{Nmol}_{Condition}"
    CHK_Name = File_Name + ".chk"
    GJF_Name = File_Name + ".gjf"
    SH_Name = f"G-{Operator}_{Condition}.sh"

    Condition_Elements = Condition.split('-')
    angle = float(Condition_Elements[0].replace('d', ''))
    Tilt = float(Condition_Elements[1].replace('d', '').replace('t', ''))
    dc1 = float(Condition_Elements[2]) / 100

    if Nmol == "3mol":
        dt1 = float(Condition_Elements[3]) / 100
    else:
        dt1 = 0

    Direction_Col, Direction_Transv, Position_Dummy, Rotate_Axis, Tilt_Axis = Axis_Setting_HB()
    Angles = []
    if Rotate_Axis == "x":
        if Tilt_Axis == "y":
            Angles = [
                [angle, Tilt, 0], [angle, Tilt, 0], [-angle, Tilt, 0]
            ]
        elif Tilt_Axis == "z":
            Angles = [
                [angle, 0, Tilt], [angle, 0, Tilt], [-angle, 0, Tilt]
            ]
        else:
            pass
    elif Rotate_Axis == "y":
        if Tilt_Axis == "x":
            Angles = [
                [Tilt, angle, 0], [Tilt, angle, 0], [Tilt, -angle, 0]
            ]
        elif Tilt_Axis == "z":
            Angles = [
                [0, angle, Tilt], [0, angle, Tilt], [0, -angle, Tilt]
            ]
        else:
            pass
    elif Rotate_Axis == "z":
        if Tilt_Axis == "x":
            Angles = [
                [Tilt, 0, angle], [Tilt, 0, angle], [Tilt, 0, -angle]
            ]
        elif Tilt_Axis == "y":
            Angles = [
                [0, Tilt, angle], [0, Tilt, angle], [0, Tilt, -angle]
            ]
        else:
            pass
    else:
        Angles = [[angle, 0, Tilt], [angle, 0, Tilt], [-angle, 0, Tilt]]
    DebugCheck(Angles, Debug, "", False)

    DummyCoords, Transitions = mkDummysPos_HB(Direction_Col, Direction_Transv, Position_Dummy, dt1, dc1)
    Element, Mol1_pos, Mol2_pos, Mol3_pos, NinMol = mkAtomList(Transitions, mk_Mater_Name, Angles[0], Angles[1],
                                                               Angles[2], "xyz")

    with open(f"Temp_Header_{Nmol}.txt", "w") as temp_header:
        if Nmol == "2mol":
            temp_header.write(Stereotyped.Header_2mol)
        elif Nmol == "3mol":
            temp_header.write(Stereotyped.Header_3mol)
    with open(f"Temp_Header_{Nmol}.txt", "r") as temp_header:
        Headers = temp_header.readlines()
    os.remove(f"Temp_Header_{Nmol}.txt")
    Headers[2] = f"%chk={CHK_Name}\n"

    if "2mol" in Nmol:
        with open(f"./{dir_path}/{mk_Mater_Name}_{Nmol}_{Condition}.gjf", "w") as zf2m:
            for Header in Headers:
                zf2m.write(Header)
            for i in range(NinMol):
                element = "{: <2}".format(Element[i])
                mol_string = format_coordinate(Mol1_pos[i], 1)
                zf2m.write(f" {element}  {mol_string}\n")
            for i in range(NinMol):
                element = "{: <2}".format(Element[i])
                mol_string = format_coordinate(Mol2_pos[i], 2)
                zf2m.write(f" {element}  {mol_string}\n")
            zf2m.write(f"\n")
    elif "3mol" in Nmol:
        with open(f"./{dir_path}/{mk_Mater_Name}_{Nmol}_{Condition}.gjf", "w") as zf3m:
            for Header in Headers:
                zf3m.write(Header)
            for i in range(NinMol):
                element = "{: <2}".format(Element[i])
                mol_string = format_coordinate(Mol1_pos[i], 1)
                zf3m.write(f" {element}  {mol_string}\n")
            for i in range(NinMol):
                element = "{: <2}".format(Element[i])
                mol_string = format_coordinate(Mol2_pos[i], 2)
                zf3m.write(f" {element}  {mol_string}\n")
            for i in range(NinMol):
                element = "{: <2}".format(Element[i])
                mol_string = format_coordinate(Mol3_pos[i], 3)
                zf3m.write(f" {element}  {mol_string}\n")
            '''
            # ダミー原子
            for i in range(6):
                mol_string = format_coordinate(DummyCoords[i], 4)
                zf3m.write(f" Ar   {mol_string}\n")
            '''
            zf3m.write(f"\n")

    orgSH = open(f"{dir_path}/G.sh", "r")
    lines = orgSH.readlines()
    lines[12] = f"g16 {GJF_Name}\n"
    orgSH.close()
    with open(f"{dir_path}/{SH_Name}", "w") as newSH:
        for line in lines:
            newSH.write(line)
    qsub_temp = f"qsub {SH_Name}"
    return qsub_temp


def mkCycleConditions(RefLines, which, dev, Nmol, Tilt):
    Degs = []
    NewConditions = []
    n = 1
    for RefLine in RefLines:
        Contents = RefLine.strip().split()
        Degs.append(round(float(Contents[0]), 1))
    print("New conditions for the next cycle:")
    for Deg in Degs:
        RefValues = getRefValues(Deg, RefLines)
        if which == "Dcol":
            SV = round(float(RefValues[0]), 2)
        elif which == "Dtrv":
            SV = round(float(RefValues[1]), 2)
        else:
            SV = int()
            CheckPoint(True)

        for i in range(n):
            NewCondition = mkNewCondition(Nmol, float(int(Deg)), SV - dev * (n - 1), Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, float(int(Deg)), SV - dev * (n - 1), Tilt, which, RefValues)
            NewConditions.append(NewCondition)

        with open(f"ConditionList_{Nmol}_t{Tilt}d.txt", "r") as file:
            orgCondition = file.readlines()
        NewList = orgCondition + NewConditions
        NewList = list(set(NewList))
        NewList.sort()
        with open(f"ConditionList_{Nmol}_t{Tilt}d.txt", "w") as file:
            for content in NewList:
                file.write(content)
    return


def CompareStructures(RefLines, LinesBefore):
    """
    比較対象の構造リストを参照構造リストと比較し、最も安定しているかを判定します。

    Args:
        RefLines (list of str): 参照構造リスト。各要素はスペースで区切られた文字列です。
        LinesBefore (list of str): 比較対象の構造リスト。各要素はスペースで区切られた文字列です。

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
        RefDeg, RefDcol, RefDtrv = RefLine.strip().split()
        for LineBefore in LinesBefore:
            Deg, Dcol, Dtrv = LineBefore.strip().split()
            if Deg == RefDeg:
                if Dcol != RefDcol or Dtrv != RefDtrv:
                    return False
    return True


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
    if len(List) == 0:
        pass
    else:
        with open(path, "w") as File:
            File.write(f"*************** {Name} Minimum Energy and Transfer Integrals ***************\n")
            File.write(Stereotyped.Combine_Header)
            print(f"*************** {Name} Minimum Energy and Transfer Integrals ***************\n")
            print(Stereotyped.Combine_Header)

            for line in List:
                File.write(line)
                print(line.strip())
    return


def combineData(MaterName, Nmol, tcalpath, resultpath, MinFileName, TIFileName, PCFileName):
    with open(f"{tcalpath}/{TIFileName}", "r") as TIFile:
        TILines = TIFile.readlines()
        del TILines[0:2]

    with (MinFileName, "r") as MinFile:
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
        print(f"The numbers of data lines in {MinFileName} and {TIFileName} DO NOT match.")
        print("Any file was NOT be changed.")
        exit()

    CombLines = []
    CombLines12 = []
    CombLines23 = []
    CombLines31 = []
    print("Entry\tAngle\tDcol\tDtrv\tCpCE\tBSE\t**\tTI-NLUMO\tTI-LUMO\tTI-HOMO\tTI-NHOMO\t**"
          "\tPC (LUMO)\tTI-LUMO\tPC (HOMO)\tTI-HOMO")
    for MinLine in MinLines:
        MinData = MinLine.split()

        if "2mol" in Nmol:
            Entry = f"{MaterName}_{Nmol}_{int(float(MinData[0]))}d-{int(round(float(MinData[1]) * 100, 5))}"
        elif "3mol" in Nmol:
            Entry = (f"{MaterName}_{Nmol}_{int(float(MinData[0]))}d-{int(round(float(MinData[1]) * 100, 5))}"
                     f"-{int(round(float(MinData[2]) * 100, 5))}")
        else:
            print("UNEXPECTED ERROR occurs in the combineData function.")
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
    saveCombData(f"{MaterName}_{Nmol}", f"{resultpath}/{MaterName}_{Nmol}_min-TIs.txt", CombLines)
    saveCombData(f"{MaterName}_{Nmol}-12", f"{resultpath}/{MaterName}_{Nmol}_min-TIs-12.txt",
                 CombLines12)
    saveCombData(f"{MaterName}_{Nmol}-23", f"{resultpath}/{MaterName}_{Nmol}_min-TIs-23.txt",
                 CombLines23)
    saveCombData(f"{MaterName}_{Nmol}-31", f"{resultpath}/{MaterName}_{Nmol}_min-TIs-31.txt",
                 CombLines31)
    print(f"{MinFileName} and {TIFileName} were combined into {resultpath}/{MaterName}_{Nmol}_min-TIs*.txt.")
    return


if __name__ == '__main__':
    main()
