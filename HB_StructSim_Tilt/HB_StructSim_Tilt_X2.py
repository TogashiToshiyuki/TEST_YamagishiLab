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

# 定数、定型文
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
Sh_txt = "\
#!/bin/sh\n\
\n\
#$ -S /bin/sh\n\
#$ -cwd\n\
#$ -V\n\
#$ -pe gau 12\n\
#$ -q all.q\n\
\n\
module load gaussian/g16\n\
export GAUSS_SCRDIR=/scr/$JOB_ID\n\
mkdir /scr/$JOB_ID\n\
\n\
g16 xxx.gjf\n\
rm -rf /scr/$JOB_ID\n\
\n"
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
Header_2mol = ("\n"
               "\n%mem=24GB"
               "\n%nprocshared=12"
               "\n%chk=test_2mol.chk"
               "\n#p pbepbe/6-31g(d) empiricaldispersion=gd3 counterpoise=2"
               "\n"
               "\nTitle Card Required"
               "\n"
               "\n0 1"
               "\n")
Header_3mol = ("\n"
               "\n%mem=24GB"
               "\n%nprocshared=12"
               "\n%chk=test_3mol.chk"
               "\n#p pbepbe/6-31g(d) empiricaldispersion=gd3 counterpoise=3"
               "\n"
               "\nTitle Card Required"
               "\n"
               "\n0 1"
               "\n")


def main():
    print(ProgramAbst)
    messages = []
    HelpList = []
    Help = False
    Debug = False

    try:
        argument1 = sys.argv[1]
        if os.path.exists(argument1) and ".xyz" in argument1:
            MaterName = argument1[0:argument1.rfind(".xyz")]
        else:
            MaterName = ''
            Help = True
            print("The file name is invalid in the program."
                  "\nPlease check files.")
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

    N = int(input("\nPlease select the number of moles to be calculated from 2 and 3."
                  "\n>>> "))
    while N != 2 and N != 3:
        print("\nError: Incorrect input.")
        N = input("\nPlease select the number of moles to be calculated from 2 and 3."
                  "\nEnter the number of moles to calculated."
                  "\n>>> ")
    else:
        pass
    Nmol = f"{N}mol"

    Tilt = int(input("\nPlease enter the angle at which the molecule is tilted."
                     "\n>>> "))
    print(f"\n{Tilt} has been entered as Tilt.")

    if "3mol" in Nmol:
        mol_pos_number = int(input("\n3 moles have been selected."
                                   "\nPlease select a structure from 1~3.."
                                   "\n>>> "))
        mol_pos = f"p{mol_pos_number}"
    else:
        mol_pos = ""

    print("\n\nFile set check...")
    if os.path.exists(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d.gjf"):
        messages.append(f"{MaterName}_{Nmol}{mol_pos}_t{Tilt}d.gjf: Found")
        HelpList.append(False)
    else:
        messages.append(f"{MaterName}_{Nmol}{mol_pos}_t{Tilt}d.gjf: NOT Found")
        HelpList.append(True)

    if os.path.exists(f"./{MaterName}.xyz"):
        messages.append(f"{MaterName}.xyz: Found")
        HelpList.append(False)
    else:
        messages.append(f"{MaterName}.xyz: NOT Found")
        HelpList.append(True)

    if os.path.exists(f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt"):
        messages.append(f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt: Found")
        HelpList.append(False)
    elif os.path.exists(f"./InitialCondition_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt"):
        messages.append(f"./InitialCondition_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt: Found")
        HelpList.append(False)
    else:
        messages.append(
            f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt and ./InitialCondition_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt: Not Found")
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

    CheckPoint(Help)

    if os.path.exists(f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt"):
        with open(f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt") as file:
            lines = file.readline()
        split_lines = lines.split('-')
        Current_Tilt = float(split_lines[1].replace('d', '').replace('t', ''))
        if not Current_Tilt == Tilt:
            print(f"*************** Caution !!! ***************"
                  f"\nThe current angle and the entered angle are different !!!"
                  f"\n"
                  f"\nInput Tilt Angle   : {Tilt}"
                  f"\nCurrent Tilt Angle : {Current_Tilt}"
                  f"\n"
                  f"\nThe process proceeds based on the angle *Entered* !!!"
                  f"\n*************** Caution !!! ***************"
                  f"\n"
                  f"\nProcessing will resume after 3 seconds.")
            time.sleep(3)
        else:
            pass
    else:
        Current_Tilt = 10000

    if os.path.exists(f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt") and Current_Tilt == Tilt:
        pass
    else:
        with open(f"./InitialCondition_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt", "r") as file:
            lines = file.readlines()

        NewConditions = []
        for line in lines:
            Deg = float((line.strip().split()[0]).replace('d', '').replace('t', ''))
            Val = float((line.strip().split()[1]).replace('d', '').replace('t', ''))
            # 3軸目方向に傾ける角度は定義済み、Tiltで呼び出せる
            RefValues = getRefValues(Deg, RefLines)

            NewCondition = mkNewCondition(Nmol, Deg, Val, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val + 0.1, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val - 0.1, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val + 0.2, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, Val - 0.2, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
        NewConditions = list(set(NewConditions))
        NewConditions.sort()
        with open(f"./ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt", "w") as file:
            for Condition in NewConditions:
                file.write(Condition)
        print(f"\nNew ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt was written "
              f"from InitialCondition_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt!!\n")
    dirpath = f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d"
    tcalpath = f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_tcal"
    os.makedirs(dirpath, exist_ok=True)
    if N == 2:
        getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dcol", [], 0.1, Tilt, Operator, mol_pos)
        print("\nCalculations for 2mol were successfully finished.")
    elif N == 3:
        RefLines = getRefLines(f"./{MaterName}_2mol_t{Tilt}d_min.txt")
        for i in range(3):
            print(f"Deg = {(i + 1) * 10}")
            RefValues = getRefValues((i + 1) * 10, RefLines)
            printList(RefValues)
            print(f"RefValues[0](SV)={RefValues[0]}")
        DebugCheck(RefLines, Debug, "RefLines before getTemp", True)
        getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dtrv", RefLines, 0.1, Tilt, Operator, mol_pos)

    temp_Structures = []
    if "3mol" in Nmol and not os.path.exists(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_mins.hist"):
        MostStable = False
        while not MostStable:
            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_min.txt")
            temp_Structures.append(RefLines)
            mkCycleConditions(RefLines, "Dcol", 0.1, Nmol, Tilt)
            getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dcol", RefLines, 0.1, Tilt, Operator, mol_pos)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_min.txt")
            mkCycleConditions(RefLines, "Dtrv", 0.1, Nmol, Tilt)
            getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dtrv", RefLines, 0.1, Tilt, Operator, mol_pos)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_min.txt")
            MostStable = CompareStructures(RefLines, temp_Structures[-1])

        MostStable = False
        while not MostStable:
            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_min.txt")
            temp_Structures.append(RefLines)
            mkCycleConditions(RefLines, "Dcol", 0.05, Nmol, Tilt)
            getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dcol", RefLines, 0.05, Tilt, Operator, mol_pos)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_min.txt")
            mkCycleConditions(RefLines, "Dtrv", 0.05, Nmol, Tilt)
            getTemporaryStructure(Nmol, MaterName, dirpath, Debug, "Dtrv", RefLines, 0.05, Tilt, Operator, mol_pos)

            RefLines = getRefLines(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_min.txt")
            MostStable = CompareStructures(RefLines, temp_Structures[-1])

        MinConditions = getMinConditions(MaterName, Nmol, Tilt)
        os.makedirs(tcalpath, exist_ok=True)

        for Condition in MinConditions:
            command = ["cp", f"{dirpath}/{MaterName}_{Nmol}{mol_pos}_{Condition}.gjf",
                       f"{tcalpath}/{MaterName}_{Nmol}{mol_pos}_{Condition}.gjf"]
            subprocess.run(command)

        if len(temp_Structures) < 100:
            file = open(f"./{MaterName}_{Nmol}{mol_pos}_mins.hist", "w")
            for i in range(len(temp_Structures)):
                file.write(f"***** Structures after the '{i + 1}'th cycle *****")
                Lines = temp_Structures[i]
                for Line in Lines:
                    file.write(Line)
            file.close()
        else:
            pass

    print("\n**********\nMaking Resulting Data Set...\n")
    resultname = f"{MaterName}_{Nmol}{mol_pos}_t{Tilt}d"
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
    if os.path.isfile(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt"):
        subprocess.run(
            ["cp", f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt", f"{resultpath}/ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt"])
        print(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt was copied into the {resultpath} forder. ")
    else:
        Lacks.append(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt")

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
        print(HelpText)
        print(AbnStop)
        exit()
    else:
        pass
    return


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
            RefLines.append(line)
        except ValueError:
            pass
    return RefLines


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
        NewCondition = f"{Deg}d-{Tilt}d-{Dcol}\n"
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
    """シェルコマンドを実行し、正常終了時は標準出力を、異常終了時はエラーメッセージを返す。

    Args:
        command (str or list): 実行するコマンド。文字列または文字列のリスト。

    Returns:
        str or None: コマンドが正常終了した場合、標準出力を返す。
                     異常終了した場合はNoneを返す。

    Raises:
        subprocess.CalledProcessError: コマンド実行中にエラーが発生した場合に発生する可能性がある。
    """
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    if result.returncode == 0:
        return result.stdout.strip()
    else:
        command = " ".join(command)
        print(f"Error executing command: {command}")
        print(f"Error message: {result.stderr.strip()}")
        return None


def has_common_element(A, B):
    return any(element in B for element in A)


def getJobID():
    """現在実行中のジョブの最初と最後のIDを取得する。

    qstatコマンドを実行し、その出力結果からジョブIDを抽出する。
    ジョブが存在しない場合は、"None"を返す。

    Returns:
        list: 最初と最後のジョブIDを含むリスト。ジョブが存在しない場合は、["None", "None"]を返す。
    """
    output = execute_command(["qstat"])
    lines = output.splitlines()
    running_job_count = len(lines) - 2
    if len(lines) == 0:
        firstJobID = "None"
        lastJobID = "None"
        JobIDList = []
    else:
        JobIDList = []
        firstJobID = lines[2].strip().split()[0]
        lastJobID = lines[-1].strip().split()[0]
        for i in range(running_job_count):
            JobID = int(lines[i + 2].strip().split()[0])
            JobIDList.append(JobID)

    JobID = [firstJobID, lastJobID]
    return JobID, JobIDList


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
    lines = glob.glob(wildcard)
    for line in lines:
        subprocess.run(["rm", line], timeout=10)
    return


def getCondition_fromName(Name):
    Condition = Name[Name.rfind("_") + 1:Name.rfind(".")]
    return Condition


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

    if "2mol" in Nmol and Condition.count("-") == 2:
        deg = Conditions[0].replace('d', '')
        tilt = (Conditions[1].replace('d', '').replace('t', ''))
        dcol = Conditions[2]
        dtrv = "na"
    elif "3mol" in Nmol and Condition.count("-") == 3:
        deg = Conditions[0].replace('d', '')
        tilt = Conditions[1].replace('d', '').replace('t', '')
        dcol = Conditions[2]
        dtrv = Conditions[3]
    else:
        print("\nUNEXPECTED ERROR happens in the function of getVALfromLogName!!")
        exit()

    Vdeg = round(float(deg), 5)
    Vdcol = round(float(dcol) / 100, 5)
    VTilt = round(int(tilt), 5)
    try:
        Vdtrv = round(float(dtrv) / 100, 5)
    except ValueError:
        Vdtrv = dtrv
    return FileName, Vdeg, Vdcol, Vdtrv, VTilt


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


def readEnergies(dir_path, MaterName, Nmol, Tilt, mol_pos):
    """指定ディレクトリ内のログファイルからエネルギー情報を読み込み、ファイルに出力する。

    指定されたディレクトリ内のログファイルを解析し、角度ごとのエネルギー情報を抽出し、
    全結果ファイルと最小エネルギーファイルに出力します。
    正常終了していないログファイルは削除されます。

    Args:
        dir_path (str): ログファイルが保存されているディレクトリパス。
        MaterName (str): 解析対象の物質名。
        Nmol (str): 解析対象の分子数。
        mol_pos (str): 3molの場合の構造

    Returns:
        None

    Raises:
        FileNotFoundError: 指定されたディレクトリパスが存在しない場合。
        ValueError: ログファイル名が不正な形式の場合。
    """
    LogList = glob.glob(f"./{dir_path}/{MaterName}_{Nmol}{mol_pos}_*.log")
    LogList.sort()
    DegList = []
    for Log in LogList:
        FileName, deg, Vdcol, Vdtrv, Tilt = getVALfromLogName(Nmol, Log)
        DegList.append(deg)
    DegList = list(set(DegList))
    DegList.sort()

    with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_all.txt", "w") as AllData:
        with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_min.txt", "w") as MinData:
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


def check_jobs(myJobID, runningJobID):
    # チェック1: runningJobIDの要素の中にmyJobIDの要素が一部でも含まれているか
    contains_my_jobs = any(job in runningJobID for job in myJobID)

    # チェック2: runningJobIDの中で、最初からmyJobIDの最後の要素までの間にいくつ要素があるか
    if myJobID[-1] in runningJobID:
        end_index = runningJobID.index(myJobID[-1]) + 1
        job_count = end_index  # runningJobIDの最初からmyJobIDの最後の要素までの間の要素数
    else:
        job_count = len([job for job in runningJobID if job <= myJobID[-1]])

    return int(contains_my_jobs), job_count


def mkNewConditionList(MaterName, Nmol, Tilt, which, dev, RefLines, mol_pos):
    """探索範囲の各角度に対して、指定された基準値との比較を行い、条件リストを更新する。

    Args:
        MaterName (str): 解析対象の物質名。
        Nmol (str): 解析対象の分子数 ("2mol" or "3mol")。
        Tilt (float): 解析対象の傾斜角度。
        which (str): 解析対象の物理量 ("Dcol" or "Dtrv")。
        dev (float): 基準値からの許容範囲。
        RefLines (list): 基準値が記載されたファイルの行リスト。
        mol_pos (str): 3モルの場合の回転角

    Returns:
        bool: 全ての角度で極小値が見つかった場合はTrue、そうでない場合はFalse。
    """
    n = 2
    with open(f"./{MaterName}_{Nmol}{mol_pos}_t{Tilt}d_all.txt", "r") as All:
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
        print(RefValues)
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
                        # if DataDtrv == RefDtrv:
                        ValueList.append(DataDcol)
                        if Contents[0] == "*":
                            SV = DataDcol
                    if which == "Dtrv":
                        RefDcol = round(float(RefValues[0]), 2)
                        # if DataDcol == RefDcol:
                        ValueList.append(DataDtrv)
                        if Contents[0] == "*":
                            SV = DataDtrv
                            print(SV)
            except ValueError:
                pass
            except IndexError:
                pass

        if round(SV + dev, 2) in ValueList and round(SV - dev, 2) in ValueList:
            Judges.append("complete")
            print(f"\nThe local minimum by {dev} step for {Deg} degree was Found;\t{SV}.")
            print("\tCOMPLETE!!")
            ComplDeg.append(Deg)
            ComplVal.append(SV)
        elif round(SV + 0.1, 2) in ValueList and round(SV - 0.2, 2) in ValueList:
            Judges.append("not complete")
            print(f"\nThe local minimum by 0.1 step for {Deg} degree was FOUND;\t{SV}.")
            print("\tNOT COMPLETE(1)!! 2 new conditions bellow were appended to the ConditionList.txt.")
            ComplDeg.append(Deg)
            ComplVal.append(SV)
            NewCondition = mkNewCondition(Nmol, Deg, SV - 0.05, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
            NewCondition = mkNewCondition(Nmol, Deg, SV + 0.05, Tilt, which, RefValues)
            NewConditions.append(NewCondition)
        elif SV == min(ValueList):
            Judges.append("not complete")
            print(f"\nLocal minimum for {Deg} degree was NOT FOUND in the cycle.")
            print(f"\tNOT COMPLETE(2)!! {n} new conditions bellow were appended to the ConditionList.txt.")
            for i in range(n):
                NewCondition = mkNewCondition(Nmol, Deg, SV - dev * (i + 1), Tilt, which, RefValues)
                NewConditions.append(NewCondition)
        elif SV == max(ValueList):
            Judges.append("not complete")
            print(f"\nLocal minimum for {Deg} degree was NOT FOUND in the cycle.")
            print(f"\tNOT COMPLETE(3)!! {n} new conditions bellow were appended to the ConditionList.txt.")
            for i in range(n):
                NewCondition = mkNewCondition(Nmol, Deg, SV + dev * (i + 1), Tilt, which, RefValues)
                NewConditions.append(NewCondition)
    print("New Condition:")
    printList(NewConditions)
    with open(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt", "r") as f:
        orgCondition = f.readlines()
    NewList = orgCondition + NewConditions
    NewList = list(set(NewList))
    NewList.sort()
    print("New List")
    printList(NewList)
    with open(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt", "w") as f:
        for content in NewList:
            f.write(content)
    if len(ComplDeg) == 0:
        print("\nLocal minimum have NOT been FOUND in any angle.")
    else:
        print(f"\nLocal minimum was FOUND in {len(ComplDeg)} angels.")
        if "2mol" in Nmol:
            print("Angel\tD_col")
        elif "3mol" in Nmol:
            print("Angel\tD_trv")
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


def getTemporaryStructure(Nmol, MaterName, dirpath, debug, which, RefLines, dev, Tilt, Operator, mol_pos=None):
    """一時的な分子構造を取得し、計算を実行する。

    指定された条件に基づいて一時的な分子構造を作成し、計算ジョブを投入・監視する。
    計算が完了すると、結果を読み込み、新しい条件リストを作成する。

    Args:
        Nmol (str): 分子の種類を表す文字列。
        MaterName (str): 材料の名前。
        dirpath (str): 作業ディレクトリのパス。
        debug (bool): デバッグモードの有効/無効。
        which (str): 計算の種類 ("Dcol" or "Dtrv")。
        RefLines (list): 基準となる行のリスト。
        dev (float): 偏差。
        Tilt (int): 傾斜角度。
        Operator (str): 演算子。
        mol_pos (int): 3モルの際の構造。(1~3)

    Returns:
        bool: 新しい条件リストが作成された場合はTrue、それ以外はFalse。
    """

    # テスト用のファイル
    """
    Conditions = []
    with open(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt") as ConditionFile:
        lines = ConditionFile.readlines()
    for line in lines:
        Conditions.append(line.strip())
    qsub_List = []
    with open(f"./{dirpath}/G.sh", "w") as original_sh:
        original_sh.write(Sh_txt)

    for Condition in Conditions:
        if os.path.exists(f"./{dirpath}/{MaterName}_{Nmol}{mol_pos}_{Condition}.log"):
            pass
        else:
            qsub_temp = mkRotate_Structure_File(MaterName, Nmol, debug, dirpath, Condition, Operator, mol_pos)
            qsub_List.append(qsub_temp)
    """
    judge = False
    while not judge:
        Conditions = []
        with open(f"ConditionList_Tilt_{Nmol}{mol_pos}_t{Tilt}d.txt") as ConditionFile:
            lines = ConditionFile.readlines()
        for line in lines:
            Conditions.append(line.strip())
        qsub_List = []
        with open(f"./{dirpath}/G.sh", "w") as original_sh:
            original_sh.write(Sh_txt)

        for Condition in Conditions:
            if os.path.exists(f"./{dirpath}/{MaterName}_{Nmol}{mol_pos}_{Condition}.log"):
                pass
            else:
                qsub_temp = mkRotate_Structure_File(MaterName, Nmol, debug, dirpath, Condition, Operator, mol_pos)
                qsub_List.append(qsub_temp)

        print("\n**********\nJobs are submitting...")

        if len(qsub_List) == 0:
            print("Any job was not submitted. Calculations with the conditions might be finished.")
        else:
            CalcEnd = False

            if "2mol" in Nmol:
                Wait_minutes = 1
            elif "3mol" in Nmol:
                Wait_minutes = 2
            else:
                Wait_minutes = 1

            for qsub in qsub_List:
                qsub = qsub.split()
                subprocess.run(qsub, cwd=f"./{dirpath}")

            JobID, JobIDList = getJobID()
            untilID = int(JobID[1])
            jobCount = len(qsub_List)
            myStartID = int(untilID - jobCount + 1)
            myJobList = list(range(myStartID, untilID + 1))
            Flag, wait_job_count = check_jobs(myJobList, JobIDList)

            start_time = datetime.datetime.now()
            formated_ST = start_time.strftime("%m/%d %H:%M:%S")
            if which == "Dcol":
                term = "the stable distance in column direction"
            elif which == "Dtrv":
                term = "the stable distance in transverse direction"
            end_time = start_time + datetime.timedelta(minutes=(wait_job_count * Wait_minutes))
            formated_end_time = end_time.strftime("%m/%d %H:%M:%S")
            print(f"\n'{jobCount}' calculations for '{term}' was submitted!! at {formated_ST}")
            print(f"\nWait until jobID {untilID}!!")
            print(f"wait for {wait_job_count * Wait_minutes} minute!"
                  f"(forecast: {formated_end_time})")
            print(f"start ID: {JobID[0]}")

            while not CalcEnd:
                JobID, JobIDList = getJobID()
                untilID = int(JobID[1])
                jobCount = len(qsub_List)
                myStartID = int(untilID - jobCount + 1)
                myJobList = list(range(myStartID, untilID + 1))
                Flag, wait_job_count = check_jobs(myJobList, JobIDList)

                formated_NOW, elapsed_time = getElapsedTime(start_time)
                Now_time = datetime.datetime.now()
                formated_end_time = (
                    (Now_time + datetime.timedelta(minutes=(wait_job_count * Wait_minutes))).strftime("%m/%d %H:%M:%S"))
                sys.stdout.write(
                    "\033[2k\033[G%s" %
                    f"\t{formated_NOW} ({elapsed_time} min. passed): Job {JobID[0]} is in progress.\n"
                    f"\tNext Check >>> {Wait_minutes * wait_job_count} minute later! forecast: {formated_end_time}")
                sys.stdout.flush()
                time.sleep(Wait_minutes * wait_job_count * 60)

                JobID, JobIDList = getJobID()
                untilID = int(JobID[1])
                jobCount = len(qsub_List)
                myStartID = int(untilID - jobCount + 1)
                myJobList = list(range(myStartID, untilID + 1))
                Flag, wait_job_count = check_jobs(myJobList, JobIDList)

                if JobID[0] == "None" and JobID[1] == "None":
                    CalcEnd = True
                elif Flag:
                    CalcEnd = True
                elif not Flag:
                    CalcEnd = False
                else:
                    print("UNEXPECTED ERROR occured in main.")
                    print(AbnStop)
                    exit()

                if not CalcEnd:
                    pass
                else:
                    print(f"\nCalculation cycles for {which} until JobID {untilID} were finished.")
                    rmWildCards(f"{dirpath}/*.sh*")
                    if debug:
                        pass
                    else:
                        rmWildCards(f"{dirpath}/*.chk")
        print("\n**********\nReading Data...\n")
        readEnergies(dirpath, MaterName, Nmol, Tilt, mol_pos)
        if "3mol" in Nmol:
            RefLines = getRefLines(f"./{MaterName}_3mol{mol_pos}_t{Tilt}d_min.txt")
        judge = mkNewConditionList(MaterName, Nmol, int(Tilt), which, dev, RefLines, mol_pos)

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


def getMinConditions(MaterName, Nmol, Tilt):
    """ファイルから最小条件のリストを取得する。

    指定されたファイル (MaterName_Nmol_Tilt_min.txt) から、
    最小条件を表す文字列のリストを取得します。
    ファイルの最初の2行はヘッダーとして無視されます。

    Args:
        MaterName (str): 材料名。
        Nmol (str): 分子数。
        Tilt (int): Tiltの値。

    Returns:
        list: 最小条件のリスト。各要素は "Deg-tTilt-Dcol" または
              "Deg-tTilt-Dcol-Dtrsv" の形式の文字列。

    Raises:
        FileNotFoundError: 指定されたファイルが存在しない場合。
    """
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


def Axis_Setting_HB():
    """軸設定を読み込む関数

    AxSetting.txtファイルから、カラム方向、横断方向、ダミー原子の座標、回転軸、チルト軸を読み込みます。
    AxSetting.txtが存在しない場合は、テンプレートファイルを作成し、プログラムを終了します。
    カラム方向と横断方向が同じ場合は、エラーメッセージを表示してプログラムを終了します。

    Returns:
        Direction_Col (str): カラム方向 ("x", "y", "z")
        Direction_Transv (str): 横断方向 ("x", "y", "z")
        Position_Dummy (np.ndarray): ダミー原子の座標 (x, y, z)
        Rotate_Axis (str): 回転軸 ("x", "y", "z")
        Tilt_Axis (str): チルト軸 ("x", "y", "z")

    Raises:
        FileNotFoundError: AxSetting.txt が見つからない場合
        SystemExit: AxSetting.txt の形式が不正な場合
    """
    try:
        with open("../Writhing/HB_StructSim_Tilt/AxSetting.txt", "r") as File:
            lines = File.readlines()
    except FileNotFoundError:
        print("\nThe AxSetting.txt is not finded in the current directory."
              "\nThe txt file should be described as follow."
              "\n*****"
              "\nColumn Direction: [x, y, or z]"
              "\nTransverse Direction: [x, y, or z]"
              "\nCoordinate of Dummy: x.xx, y.yy, z.zz"
              "\nRotate Axis: [x, y, or z]"
              "\nTilt Axis: [x, y, or z]"
              "\n*****"
              "\n"
              "\nMake the correct AxSetting.txt and then restart the program."
              "\n")
        with open("../Writhing/HB_StructSim_Tilt/AxSetting.txt", "w") as file:
            file.write("Column Direction:  \n")
            file.write("Transverse Direction: \n")
            file.write("Coordinate of Dummy: x.xx,y.yy,z.zz\n")
            file.write("Rotate Axis: [x, y, or z]\n")
            file.write("Tilt Axis: [x, y, or z]\n")
        print(AbnStop)
        exit()

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
        print(AbnStop)
        exit()

    return Direction_Col, Direction_Transv, Position_Dummy, Rotate_Axis, Tilt_Axis


def mkCycleConditions(RefLines, which, dev, Nmol, Tilt):
    """次のサイクルの実験条件を生成し、既存条件と合わせてファイルに保存する。

    既存の実験条件ファイルを読み込み、指定されたパラメータに基づいて新たな実験条件を生成します。
    生成された条件は既存の条件と重複しないように調整され、ソートされた上でファイルに書き戻されます。

    Args:
        RefLines (list): 参照となる実験条件のリスト（各要素は "角度 実験値1 実験値2..." の形式の文字列）。
        which (str): 生成する実験条件の種類 ("Dcol" or "Dtrv")。
        dev (float): 実験値からの偏差幅。
        Nmol (str): サンプルのモル数。
        Tilt (float): サンプルの傾斜角度。

    Returns:
        None

    Raises:
        ValueError: which が "Dcol" または "Dtrv" 以外の値の場合。
    """
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


def mkDirection(axis_direction, input_axis, axis_name):
    if input_axis == "x":
        matrix = np.array([axis_direction, 0.0, 0.0])
    elif input_axis == "y":
        matrix = np.array([0.0, axis_direction, 0.0])
    elif input_axis == "z":
        matrix = np.array([0.0, 0.0, axis_direction])
    else:
        print(f"\nThe axis of {axis_name} is INVALID.")
        print(AbnStop)
        exit()
    return matrix


def mkDummysPos_HB(Direction_Col, Direction_Transv, Position_Dummy, Dt1, Dc1):
    """ダミー原子の位置座標と並進ベクトルを生成する。

    カラム方向と横方向の軸に基づいて、6つのダミー原子を生成し、
    それぞれの位置座標をリストとして返します。
    また、カラム方向、横方向、その他の方向の並進ベクトルも返します。

    Args:
        Direction_Col (str): カラム方向の軸 ("x", "y", or "z")。
        Direction_Transv (str): 横方向の軸 ("x", "y", or "z")。
        Position_Dummy (list): ダミー原子の基準位置座標 [x, y, z]。
        Dt1 (float): 横方向の並進距離。
        Dc1 (float): カラム方向の並進距離。

    Returns:
        list: 6つのダミー原子の位置座標 [[x1, y1, z1], [x2, y2, z2], ...]。
        list: 3つの並進ベクトル [[tx1, ty1, tz1], [tx2, ty2, tz2], [tx3, ty3, tz3]]。

    Raises:
        SystemExit: カラム方向と横方向の軸が同じ場合に発生。
    """
    if Direction_Col == Direction_Transv:
        print("""
        The axis of transverse direction is same with the axis of column direction.""")
        print(AbnStop)
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
    """座標を指定された軸周りに回転させた後の座標を計算する。

    この関数は、3次元座標を指定された回転軸(x, y, z)周りに指定された角度で回転させ、
    回転後の新しい座標を計算します。回転はオイラー角の順序で行われます。

    Args:
        Current (numpy.ndarray): 回転させたい3次元座標 (形状: (3,))
        Tx (float): x軸周りの回転角度 (単位: 度)
        Ty (float): y軸周りの回転角度 (単位: 度)
        Tz (float): z軸周りの回転角度 (単位: 度)
        rotation (str): 回転軸の順序を表す文字列 ('x', 'y', 'z'の組み合わせ)

    Returns:
        numpy.ndarray: 回転後の新しい3次元座標 (形状: (3,))
    """
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
    """xyzファイルから原子リストを生成する。

    指定されたxyzファイルを読み込み、各原子座標を回転・並進操作で変換し、3つの分子リストを作成する。

    Args:
        Translations (list): 並進ベクトル。
        MaterName (str): xyzファイル名。
        Angles (list): 回転角度。
        mol2_Angles (list): 2番目の分子の回転角度。
        mol3_Angles (list): 3番目の分子の回転角度。
        rotate (str): 回転順。

    Returns:
        tuple: 原子種リスト、1番目の分子リスト、2番目の分子リスト、3番目の分子リスト、分子内の原子数。
    """
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


def mkRotate_Structure_File(mk_Mater_Name, Nmol, Debug, dir_path, Condition, Operator, mol_pos=None):
    """Gaussian計算用のジョブファイル（.gjf）とシェルスクリプト（.sh）を作成する関数

    指定された条件に基づいて、Gaussian計算用の入力ファイル（.gjf）と
    それを実行するためのシェルスクリプト（.sh）を作成します。
    ファイル名は、材料名、分子数、計算条件などから自動生成されます。

    Args:
        mk_Mater_Name (str): 材料名
        Nmol (str): 分子数 ("2mol" or "3mol")
        Debug (bool): デバッグモードの有効/無効
        dir_path (str): ファイルの出力先ディレクトリパス
        Condition (str): 計算条件（例: "15d0t100-100"）
        Operator (str): Gaussianの演算子

    Returns:
        str: 作成されたシェルスクリプトを実行するためのqsubコマンド

    Raises:
        SystemExit: 回転軸の設定が不正な場合
    """
    File_Name = f"{mk_Mater_Name}_{Nmol}{mol_pos}_{Condition}"
    CHK_Name = File_Name + ".chk"
    GJF_Name = File_Name + ".gjf"
    SH_Name = f"G-{Operator}_{Condition}.sh"

    Condition_Elements = Condition.split('-')
    Angle = float(Condition_Elements[0].replace('d', ''))
    Tilt = float(Condition_Elements[1].replace('d', '').replace('t', ''))
    dc1 = float(Condition_Elements[2]) / 100

    if Nmol == "3mol":
        dt1 = float(Condition_Elements[3]) / 100
    else:
        dt1 = 0

    Direction_Col, Direction_Transv, Position_Dummy, Rotate_Axis, Tilt_Axis = Axis_Setting_HB()

    if "2mol" in Nmol:
        if Rotate_Axis == "x":
            if Tilt_Axis == "y":
                Angles = [
                    [Angle, Tilt, 0], [Angle, Tilt, 0], [-Angle, Tilt, 0]
                ]
            elif Tilt_Axis == "z":
                Angles = [
                    [Angle, 0, Tilt], [Angle, 0, Tilt], [-Angle, 0, Tilt]
                ]
            else:
                pass
        elif Rotate_Axis == "y":
            if Tilt_Axis == "x":
                Angles = [
                    [Tilt, Angle, 0], [Tilt, Angle, 0], [Tilt, -Angle, 0]
                ]
            elif Tilt_Axis == "z":
                Angles = [
                    [0, Angle, Tilt], [0, Angle, Tilt], [0, -Angle, Tilt]
                ]
            else:
                pass
        elif Rotate_Axis == "z":
            if Tilt_Axis == "x":
                Angles = [
                    [Tilt, 0, Angle], [Tilt, 0, Angle], [Tilt, 0, -Angle]
                ]
            elif Tilt_Axis == "y":
                Angles = [
                    [0, Tilt, Angle], [0, Tilt, Angle], [0, Tilt, -Angle]
                ]
            else:
                pass
        else:
            Angles = [[Angle, 0, Tilt], [Angle, 0, Tilt], [-Angle, 0, Tilt]]
    elif "3mol" in Nmol:
        if "1" in mol_pos:
            if Rotate_Axis == "x":
                if Tilt_Axis == "y":
                    Angles = [
                        [Angle, Tilt, 0], [Angle, Tilt, 0], [-Angle, Tilt, 0]
                    ]
                elif Tilt_Axis == "z":
                    Angles = [
                        [Angle, 0, Tilt], [Angle, 0, Tilt], [-Angle, 0, Tilt]
                    ]
                else:
                    pass
            elif Rotate_Axis == "y":
                if Tilt_Axis == "x":
                    Angles = [
                        [Tilt, Angle, 0], [Tilt, Angle, 0], [Tilt, -Angle, 0]
                    ]
                elif Tilt_Axis == "z":
                    Angles = [
                        [0, Angle, Tilt], [0, Angle, Tilt], [0, -Angle, Tilt]
                    ]
                else:
                    pass
            elif Rotate_Axis == "z":
                if Tilt_Axis == "x":
                    Angles = [
                        [Tilt, 0, Angle], [Tilt, 0, Angle], [Tilt, 0, -Angle]
                    ]
                elif Tilt_Axis == "y":
                    Angles = [
                        [0, Tilt, Angle], [0, Tilt, Angle], [0, Tilt, -Angle]
                    ]
                else:
                    pass
            else:
                Angles = [[Angle, 0, Tilt], [Angle, 0, Tilt], [-Angle, 0, Tilt]]
        elif "2" in mol_pos:
            if Rotate_Axis == "x":
                if Tilt_Axis == "y":
                    Angles = [
                        [-Angle, Tilt, 0], [-Angle, Tilt, 0], [Angle, Tilt, 0]
                    ]
                elif Tilt_Axis == "z":
                    Angles = [
                        [-Angle, 0, Tilt], [-Angle, 0, Tilt], [Angle, 0, Tilt]
                    ]
                else:
                    pass
            elif Rotate_Axis == "y":
                if Tilt_Axis == "x":
                    Angles = [
                        [Tilt, -Angle, 0], [Tilt, -Angle, 0], [Tilt, Angle, 0]
                    ]
                elif Tilt_Axis == "z":
                    Angles = [
                        [0, -Angle, Tilt], [0, -Angle, Tilt], [0, Angle, Tilt]
                    ]
                else:
                    pass
            elif Rotate_Axis == "z":
                if Tilt_Axis == "x":
                    Angles = [
                        [Tilt, 0, -Angle], [Tilt, 0, -Angle], [Tilt, 0, Angle]
                    ]
                elif Tilt_Axis == "y":
                    Angles = [
                        [0, Tilt, -Angle], [0, Tilt, -Angle], [0, Tilt, Angle]
                    ]
                else:
                    pass
            else:
                Angles = [[Angle, 0, Tilt], [Angle, 0, Tilt], [-Angle, 0, Tilt]]
        elif "3" in mol_pos:
            if Rotate_Axis == "x":
                if Tilt_Axis == "y":
                    Angles = [
                        [Angle, Tilt, 0], [Angle, Tilt, 0], [-Angle-180, Tilt, 0]
                    ]
                elif Tilt_Axis == "z":
                    Angles = [
                        [Angle, 0, Tilt], [Angle, 0, Tilt], [-Angle-180, 0, Tilt]
                    ]
                else:
                    pass
            elif Rotate_Axis == "y":
                if Tilt_Axis == "x":
                    Angles = [
                        [Tilt, Angle, 0], [Tilt, Angle, 0], [Tilt, -Angle-180, 0]
                    ]
                elif Tilt_Axis == "z":
                    Angles = [
                        [0, Angle, Tilt], [0, Angle, Tilt], [0, -Angle-180, Tilt]
                    ]
                else:
                    pass
            elif Rotate_Axis == "z":
                if Tilt_Axis == "x":
                    Angles = [
                        [Tilt, 0, Angle], [Tilt, 0, Angle], [Tilt, 0, -Angle-180]
                    ]
                elif Tilt_Axis == "y":
                    Angles = [
                        [0, Tilt, Angle], [0, Tilt, Angle], [0, Tilt, -Angle-180]
                    ]
                else:
                    pass
            else:
                Angles = [[Angle, 0, Tilt], [Angle, 0, Tilt], [-Angle, 0, Tilt]]

    DebugCheck(Angles, Debug, "rotate axis", False)

    DummyCoords, Transitions = mkDummysPos_HB(Direction_Col, Direction_Transv, Position_Dummy, dt1, dc1)
    Element, Mol1_pos, Mol2_pos, Mol3_pos, NinMol = mkAtomList(Transitions, mk_Mater_Name, Angles[0], Angles[1],
                                                               Angles[2], "xyz")

    with open(f"Temp_Header_{Nmol}.txt", "w") as temp_header:
        temp_header.write(globals()[f"Header_{Nmol}"])
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
        print(f"\n{mk_Mater_Name}_{Nmol}_{Condition}.gjf have been created.")
    elif "3mol" in Nmol:
        with open(f"./{dir_path}/{mk_Mater_Name}_{Nmol}{mol_pos}_{Condition}.gjf", "w") as zf3m:
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
        print(f"\n{mk_Mater_Name}_{Nmol}{mol_pos}_{Condition}.gjf have been created.")

    with open(f"{dir_path}/G.sh", "r") as orgSH:
        lines = orgSH.readlines()
        lines[12] = f"g16 {GJF_Name}\n"
    with open(f"{dir_path}/{SH_Name}", "w") as newSH:
        for line in lines:
            newSH.write(line)
    print(f"{SH_Name} have been created.")
    qsub_temp = f"qsub {SH_Name}"
    return qsub_temp


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
            print(header)
            for line in List:
                File.write(line)
                print(line.strip())
    return


def combineData(MaterName, Nmol, tcalpath, resultpath, MinFileName, TIFileName, PCFileName):
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

    Raises:
        SystemExit:
            - MinFileNameとTIFileNameのデータ行数が一致しない場合。
            - Nmolが"2mol"でも"3mol"でもない場合。
    """
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
    print("Entry\tAngle\tDcol\tDtrv\tCpCE\tBSE\t**"
          "\tTI-NLUMO\tTI-LUMO\tTI-HOMO\tTI-NHOMO\t**\tPC (LUMO)\tTI-LUMO\tPC (HOMO)\tTI-HOMO")
    for MinLine in MinLines:
        MinData = MinLine.split()

        if "2mol" in Nmol:
            Entry = f"{MaterName}_{Nmol}_{int(float(MinData[0]))}d-{int(round(float(MinData[1]) * 100, 5))}"
        elif "3mol" in Nmol:
            Entry = f"{MaterName}_{Nmol}_{int(float(MinData[0]))}d-{int(round(float(MinData[1]) * 100, 5))}-{int(round(float(MinData[2]) * 100, 5))}"
        else:
            print("UNEXPECTED ERROR occurs in the combineData function.")
            exit()

        for TILine in TILines:
            TIData = TILine.split()
            if (TIData[0] == Entry
                    or TIData[0] == f"{Entry}-12"
                    or TIData[0] == f"{Entry}-23"
                    or TIData[0] == f"{Entry}-31"):
                CombLine_temp = f"{TIData[0]}\t{MinData[0]}\t{MinData[1]}\t{MinData[2]}\t{MinData[3]}\t{MinData[4]}\t**\t{TIData[1]}\t{TIData[2]}\t{TIData[3]}\t{TIData[4]}"

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


if __name__ == "__main__":
    main()
