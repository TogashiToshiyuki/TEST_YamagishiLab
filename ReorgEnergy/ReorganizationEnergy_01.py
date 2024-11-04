#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 2024
Program for calculating relocation energy.
Created by Toshiyuki Togashi

Derived from development version X2
"""

import argparse
import os
import subprocess
import datetime
import sys
import time
import glob
import functools


print = functools.partial(print, flush=True)


class Stereotyped:
    # プログラムの概要
    ProgramAbst = ("\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n"
                   "*   Programme for calculating relocation energy.            *\n"
                   "*                                                           *\n"
                   "*   This programme requires the following files.            *\n"
                   "*     - A gjf file of the molecule you want to calculate.   *\n"
                   "*                                                           *\n"
                   "*   This programme requires the following arguments.        *\n"
                   "*     - Molecule name.                                      *\n"
                   "*                created by Toshiyuki Togashi 2024/10/10    *\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n")
    # 不正終了
    AbnormalEnd = "\n************* Programme DID NOT terminate successfully. *************\n"
    Header_Template = ("%nprocshared=12\n"
                       "%mem=32GB\n"
                       "%chk=Template.chk\n"
                       "# opt freq=noraman b3lyp/6-31g(d) geom=connectivity\n"
                       "\n"
                       "Title Card Required\n"
                       "\n"
                       "0 1\n")
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
              "rm -rf /scr/$JOB_ID\n")


class Constants:
    # ハートリーからeVに変換するための定数
    Hartree_to_eV = 27.21138602
    # 気体定数[R] (J/(mol・K))
    R = 8.314462618
    # 真空中の光速[m/s]
    C = 299792458
    # 電気素量[C]
    E = 1.602176634e-19


class PeriodicTable:
    # 原子番号に対応する元素記号と原子量の辞書
    # memo: 他の原子番号も追加可能
    atomic_symbols = {
        1: ['H', 1.00794],  # 水素
        2: ['He', 4.002602],  # ヘリウム
        3: ['Li', 6.941],  # リチウム
        4: ['Be', 9.012182],  # ベリリウム
        5: ['B', 10.811],  # ホウ素
        6: ['C', 12.0107],  # 炭素
        7: ['N', 14.00674],  # 窒素
        8: ['O', 15.9994],  # 酸素
        9: ['F', 18.9984032],  # フッ素
        10: ['Ne', 20.1797],  # ネオン
        11: ['Na', 22.989770],  # ナトリウム
        12: ['Mg', 24.3050],  # マグネシウム
        13: ['Al', 26.981538],  # アルミニウム
        14: ['Si', 28.0855],  # ケイ素
        15: ['P', 30.973761],  # リン
        16: ['S', 32.066],  # 硫黄
        17: ['Cl', 35.4527],  # 塩素
        18: ['Ar', 39.948],  # アルゴン
        19: ['K', 39.0983],  # カリウム
        20: ['Ca', 40.078],  # カルシウム
        21: ['Sc', 44.955910],  # スカンジウム
        22: ['Ti', 47.867],  # チタン
        23: ['V', 50.9415],  # バナジウム
        24: ['Cr', 51.9961],  # クロム
        25: ['Mn', 54.938049],  # マンガン
        26: ['Fe', 55.845],  # 鉄
        27: ['Co', 58.933200],  # コバルト
        28: ['Ni', 58.6934],  # ニッケル
        29: ['Cu', 63.546],  # 銅
        30: ['Zn', 65.409],  # 亜鉛
        31: ['Ga', 69.723],  # ガリウム
        32: ['Ge', 72.61],  # ゲルマニウム
        33: ['As', 74.92160],  # ヒ素
        34: ['Se', 78.96],  # セレン
        35: ['Br', 79.904],  # 臭素
        36: ['Kr', 83.80],  # クリプトン
        # 必要に応じて他の原子番号も追加
    }


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


class Basefunctions:
    function = {
        1: ["b3lyp/6-31g(d)", "b3lyp_6-31Gd"],
        2: ["b3lyp/6-311+g(d,p)", "b3lyp_6-311+Gdp"],
        3: ["exit", "none"],  # 終了
    }


# メイン関数
def main():
    # プログラムスタート
    # 変数の初期化
    messages, HelpList, Debug, MaterName, EnergyList, Addition_Information = [], [], False, "", [], []

    # 引数の処理
    args, Debug = arg_parser()
    MolecularName = args.MolecularName

    # プログラムの概要を表示
    print(Stereotyped.ProgramAbst)

    # ファイルの存在確認
    check_files(args.MolecularName, messages, HelpList)

    # オペレーター名を取得
    operator_name = get_operator_name(Debug, Addition_Information)

    # 結合情報を取得
    bond_info = get_bond_info(args.MolecularName, messages, HelpList, Debug, Addition_Information)

    # 基底関数を設定
    basis_function = set_basis_function(args)

    # 0価gjfファイルを作成し、これを計算
    calculate_from_gjf(MolecularName, messages, HelpList, Debug,
                       "EG", "0", operator_name, bond_info, Addition_Information, basis_function)

    # 0価の構造最適化済みファイルから1価gjfファイルを作成し、これを計算
    calculate_from_log(MolecularName, f"{MolecularName}_0_EG_b3lyp_6-31Gd.log", messages, HelpList, Debug,
                       "EG", "+1", operator_name, bond_info, Addition_Information, basis_function)

    # 1価の構造最適化済みファイルから0価gjfファイルを作成し、これを計算
    calculate_from_log(MolecularName, f"{MolecularName}_+1_EG_b3lyp_6-31Gd.log", messages, HelpList, Debug,
                       "SP", "+0", operator_name, bond_info, Addition_Information, basis_function)

    # 0価の構造最適化済みファイルから1価gjfファイルを作成し、これを計算
    calculate_from_log(MolecularName, f"{MolecularName}_0_EG_b3lyp_6-31Gd.log", messages, HelpList, Debug,
                       "SP", "+1", operator_name, bond_info, Addition_Information, basis_function)

    # エネルギーを取得し、再配置エネルギーを計算
    EnergyList, result_meV = Calculate_ReorgEnergy(MolecularName, EnergyList, Debug, basis_function)

    # 計算結果をファイルに保存
    save_results(MolecularName, EnergyList, result_meV, Addition_Information, operator_name, basis_function, args)

    # プログラム終了
    print(f"\n{Color.GREEN}************* Programme terminated successfully. *************{Color.RESET}\n")


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
    if HelpList:
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


# 引数の処理を行う関数
def arg_parser():
    """
    引数の処理を行う関数
    :return: args : 引数の情報
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description=Stereotyped.ProgramAbst, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('MolecularName', help="Molecular name")
    parser.add_argument('--debug', '-d', '-D',
                        help="Debug mode", action='store_false')
    parser.add_argument('--function', '-f', '-F',
                        help="Change the basis function", action='store_true')
    args = parser.parse_args()
    if args.debug:
        Debug = True
        print(f"\n{Color.RED}*************** Caution!!! Debug Started!!! ***************{Color.RESET}")
        print(f"\t>>> mol_name: {args.MolecularName}\n")
    else:
        Debug = False
    return args, Debug


# 必要なファイルが揃っているかを確認する関数
def check_files(MolecularName, messages, HelpList):
    """
    Check if the required files are available.
    :return:
    """
    print("Checking if the required files are available...")
    if os.path.exists(MolecularName + ".gjf"):
        messages.append(f"\t>>> The {MolecularName}.gjf file is available.")
    else:
        messages.append(f"\t{Color.RED}>>> The {MolecularName}.gjf file is not available.{Color.RESET}")
        HelpList.append(True)
    message_show(messages)
    messages.clear()
    if HelpList:
        print(f"{Color.RED}\t>>> Some required files are not available.\n"
              f"\t>>> Please upload the required files and try again.{Color.RESET}")
    help_check_exit(messages, HelpList)
    print(f"\n\t{Color.GREEN}>>> All required files are available.{Color.RESET}")
    return None


# オペレーター名を取得する関数
def get_operator_name(Debug, Addition_Information):
    """
    Get the name of the operator.
    :return:
    """
    print("\n"
          "Retrieve the operator name.\n"
          "\t>>> The name entered here will be used to identify the operator.")
    operator_name = input("\tEnter your name. >>> ")

    if not operator_name:
        operator_name = "ONE"
        print("\n"
              "\t>>> No name entered. Operator name set to 'ONE'.")

    if Debug:
        print(f"\n\t>>> Operator name: {operator_name}")
        Addition_Information.append(("\n********************************** Additional Information "
                                     "**********************************\n"))

    return operator_name


# 結合情報を取得する関数
def get_bond_info(MolecularName, messages, HelpList, Debug, Addition_Information):
    """
    Get bond information.
    :param Addition_Information:
    :param Debug:
    :param MolecularName: 分子名
    :type MolecularName: str
    :param messages: 表示するメッセージ
    :type messages: list
    :param HelpList: ヘルプリスト
    :type HelpList: list
    :return:
    """
    print("\nCoupling information is obtained...")
    # 結合情報を格納するリスト
    bond_info = []
    try:
        with open(f"{MolecularName}.gjf", 'r') as file:
            lines = file.readlines()
            is_coordinates_section = False
            is_additional_input_section = False

            for line in lines:
                # Check for the start of coordinates section (系の電荷とスピン多重度の次の行)
                if not is_coordinates_section and line.strip() and all(
                        c.isdigit() or c.isspace() for c in line.strip().split()):
                    is_coordinates_section = True
                    continue

                # Look for the blank line after a coordinates section to start capturing additional input
                if is_coordinates_section and line.strip() == "":
                    is_additional_input_section = True
                    continue

                # Append lines to additional input if we're in the additional input section
                if is_additional_input_section and line.strip():
                    bond_info.append(line.strip())
        if Debug:
            print("****** The bond information ******")
            print("\n".join(bond_info))
            print("**********************************")
            Addition_Information.append("\n****** The bond information ******\n")
            Addition_Information.append("\n".join(bond_info))
            Addition_Information.append("\n**********************************\n")
        messages.append(f"\t{Color.GREEN}>>> The bond information was successfully obtained.{Color.RESET}")
    except FileNotFoundError:
        messages.append(f"\t{Color.RED}>>> The {MolecularName}.gjf file is not available.{Color.RESET}")
        HelpList.append(True)
    help_check_exit(messages, HelpList)
    return "\n".join(bond_info)


# 基底関数を設定する関数
def set_basis_function(args):
    """
    Set the basis function.
    :return:
    """
    print("\nSetting the basis function...")
    if args.function:
        while True:
            print(f"\n"
                  f"{Color.BLUE}Selected basis function{Color.RESET}")
            for key, value in Basefunctions.function.items():
                print(f"\t{key}: {value[0]}")
            function_set = input("\t>>> ")
            try:
                basis_function = Basefunctions.function.get(int(function_set))
                if basis_function == "exit":
                    print(f"\t{Color.RED}>>> Exit the program.{Color.RESET}")
                    exit()
                else:
                    pass
                print(f"\t>>> The basis function is set to {basis_function[1]}.")
                return basis_function
            except (ValueError, TypeError):
                print(f"{Color.RED}Error: Invalid input.{Color.RESET}")
                break
    basis_function = ["b3lyp/6-31g(d)", "b3lyp_6-31Gd"]
    print(f"\t>>> The basis function is set to {basis_function[1]}.")
    return basis_function


# gjfファイルから計算用gjfファイルを作成し、これを計算する関数
def calculate_from_gjf(MolecularName, messages, HelpList, Debug,
                       EG_or_SP, Charge, Operator, bond_info, Addition_Information, basis_function):
    print(f"{Color.BLUE}\n"
          f"Calculate '{MolecularName} {Charge} {EG_or_SP}'{Color.RESET}")
    # 計算したコンディションかどうかを確認し、計算していない場合は計算を実行する
    messages.append(
        f"Checking whether {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log is present..."
    )
    if os.path.exists(f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log"):
        messages.append(f"\t>>> {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log is available.")
        messages.append(
            f"\t>>> {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.gjf has already been calculated.")
        help_check_exit(messages, HelpList)
        mono_molecular_log_open(f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log", messages, HelpList)
    else:
        messages.append(
            f"\t>>> {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log is not available.\n"
            f"\t>>> Probably not calculated."
        )
        help_check_exit(messages, HelpList)

        # gjfファイルからの座標の抜き出し
        with open(f"{MolecularName}.gjf", "r") as file:
            data = file.read()
        # 最初の８個の要素は飛ばす
        List = data.splitlines()[6:]
        # 元素記号の取得
        data = []
        for element in List:
            try:
                coordinate = [float(element.split()[1]), float(element.split()[2]), float(element.split()[3])]
                formatted_coordinate = format_coordinate(coordinate)
                formatted_line = f" {element.split()[0]}                 {formatted_coordinate}"
                data.append(formatted_line)
            except IndexError:
                break

        # GJFファイルへの書き込み
        write_gjf_file(MolecularName, messages, HelpList, Debug, EG_or_SP, Charge, data, bond_info, basis_function)

        # shファイルの作成
        qsub_temp = write_sh_file(messages, HelpList, EG_or_SP, Charge, Operator,
                                  f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.gjf", basis_function)

        # 計算の実行
        job_execution(qsub_temp, Debug)

    # 計算が正常に終了したかを確認
    normal_termination_check(f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log", messages, HelpList, True)

    # 周波数にマイナスがないかを確認
    check_negative_frequency(f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log",
                             messages, HelpList, Debug, Addition_Information)
    return None


# 構造最適化済みlogファイルからgjfファイルを作成し、これを計算する関数
def calculate_from_log(MolecularName, LogFileName, messages, HelpList, Debug,
                       EG_or_SP, Charge, Operator, bond_info, Addition_Information, basis_function):
    """
    0価の構造最適化済みファイルからgjfファイルを作成し、これを計算する関数
    :param basis_function:
    :param Addition_Information:
    :param LogFileName:
    :param bond_info:
    :param Operator:
    :param Charge:
    :param EG_or_SP:
    :param Debug:
    :param HelpList:
    :param messages:
    :param MolecularName: 分子名
    :type MolecularName: str
    """
    print(f"{Color.BLUE}\n"
          f"Calculate '{MolecularName} {Charge} {EG_or_SP}'{Color.RESET}")
    # 計算したコンディションかどうかを確認し、計算していない場合は計算を実行する
    messages.append(
        f"Checking whether {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log is present..."
    )
    if os.path.exists(f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log"):
        messages.append(
            f"\t>>> {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log is available."
        )
        messages.append(
            f"\t>>> {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.gjf has already been calculated.")
        help_check_exit(messages, HelpList)
        mono_molecular_log_open(f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log", messages, HelpList)
    else:
        messages.append(
            f"\t>>> {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log is not available.\n"
            f"\t>>> Probably not calculated."
        )
        help_check_exit(messages, HelpList)

        # 構造最適化済みlogファイルを開く
        data = mono_molecular_log_open(f"{LogFileName}", messages, HelpList)

        # 座標データの抽出
        formatted_data = format_atomic_data(data, messages, HelpList)
        if Debug:
            print(f"****** Coordinate data generated from {MolecularName}.log ******\n"
                  f"chemical symbol\tCoordinate")
            for formatted_line in formatted_data:
                print(formatted_line)
            print("********************************************************************\n")

        # GJFファイルへの書き込み
        write_gjf_file(MolecularName, messages, HelpList, Debug,
                       EG_or_SP, Charge, formatted_data, bond_info, basis_function)

        # shファイルの作成
        qsub_temp = write_sh_file(messages, HelpList, EG_or_SP, Charge, Operator,
                                  f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.gjf", basis_function)

        # 計算の実行
        job_execution(qsub_temp, Debug)

    # 計算が正常に終了したかを確認
    normal_termination_check(f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log", messages, HelpList, True)

    # EGの場合、周波数にマイナスがないかを確認
    if EG_or_SP == "EG":
        check_negative_frequency(f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log",
                                 messages, HelpList, Debug, Addition_Information)
    return None


# 構造最適化済みlogファイルを開く関数
def mono_molecular_log_open(MaterNameLog, messages, HelpList):
    """
    単分子の構造最適化済みlogファイルを開く関数
    :param HelpList:
    :param messages:
    :param MaterNameLog: 構造最適化済みファイル
    :type MaterNameLog: str
    :return: data : 単分子の構造最適化済みlogファイルのデータ
    :rtype: str
    """
    # MaterNameLogの存在確認
    print(f"\nChecking whether {MaterNameLog} is present...")
    try:
        with open(MaterNameLog, "r"):
            messages.append(f"\t>>> {MaterNameLog} is found.")
    except FileNotFoundError:
        messages.append(f"{Color.RED}\t>>> {MaterNameLog} is NOT FOUND.\n"
                        f"\t>>> Check that the correct file is specified.{Color.RESET}")
        HelpList.append(True)
    help_check_exit(messages, HelpList)

    # 単分子の構造最適化済みlogファイルを開く
    with open(MaterNameLog, "r") as file:
        data = file.read()

    # 正常に計算終了したデータかどうかを確認
    normal_termination_check(MaterNameLog, messages, HelpList)
    return data


# 正常に計算終了したデータかどうかを確認する関数
def normal_termination_check(MaterNameLog, messages, HelpList, Debug=False):
    """
    正常に計算終了したデータかどうかを確認する関数
    :param Debug:
    :param MaterNameLog:
    :param HelpList:
    :param messages:
    """
    with open(MaterNameLog, "r") as file:
        data = file.read()
    if Debug:
        print("\nChecking that the calculation has been completed successfully...")
    if "Normal termination" in data:
        rmWildCards(f"./*.sh*")
        if Debug:
            messages.append(f"\t{Color.GREEN}>>> Successfully completed.{Color.RESET}")
    else:
        HelpList.append(True)
        if HelpList:
            messages.append(f"\t{Color.RED}>>> DID NOT COMPLETE SUCCESSFULLY.\n"
                            f"\t>>> {MaterNameLog} was not normally termination.\n"
                            f"\t>>> Please check the {MaterNameLog} file.{Color.RESET}")
    help_check_exit(messages, HelpList)
    return None


# データを加工する関数
def format_atomic_data(data, messages, HelpList):
    """
    データを加工する関数
    :param HelpList:
    :param messages:
    :param data:
    :return:
    """
    # 出力用のリスト
    formatted_data = []
    data = data.split("Standard orientation")[-1]
    data = data.split("---------------------------------------------------------------------")[2]
    for line in data.strip().splitlines():
        parts = line.split()

        # 各要素を取り出す
        atomic_number = int(parts[1])
        coordinate = [float(parts[3]), float(parts[4]), float(parts[5])]
        formatted_coordinate = format_coordinate(coordinate)

        # 元素記号を取得
        symbol = None
        try:
            symbol = PeriodicTable.atomic_symbols.get(atomic_number)[0]
        except KeyError:
            messages.append(f"{Color.RED}Error: Unknown atomic number {atomic_number}{Color.RESET}")
            HelpList.append(True)
        help_check_exit(messages, HelpList)

        # フォーマットされた行を作成
        formatted_line = f" {symbol}                 {formatted_coordinate}"
        formatted_data.append(formatted_line)

    return formatted_data


# 座標をフォーマットする関数
def format_coordinate(coord):
    """
    座標をフォーマットする関数
    :param coord:
    :return:
    """
    return f"{coord[0]: 11.8f}   {coord[1]: 11.8f}   {coord[2]: 11.8f}"


# GJFファイルへの書き込み関数
def write_gjf_file(MolecularName, messages, HelpList, Debug,
                   EG_or_SP, Charge, formatted_data, footer_data, basis_function):
    """
    GJFファイルへの書き込み関数
    :param basis_function:
    :param footer_data:
    :param formatted_data:
    :param Debug:
    :param Charge:
    :param EG_or_SP:
    :param HelpList:
    :param messages:
    :param MolecularName:
    """
    print(f"Creating {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.gjf...")
    # 仮のヘッダーを作成
    header_data = Stereotyped.Header_Template.splitlines()
    header_data[0] = header_data[0] + "\n"
    header_data[1] = header_data[1] + "\n"
    header_data[2] = f"%chk={MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.chk" + "\n"
    # EGかSPかで分岐
    if EG_or_SP == "EG":
        header_data[3] = f"# opt freq=noraman {basis_function[0]} geom=connectivity" + "\n"
    elif EG_or_SP == "SP":
        header_data[3] = f"# {basis_function[0]} geom=connectivity" + "\n"
    else:
        messages.append(f"{Color.RED}Error: Unknown option {EG_or_SP}{Color.RESET}")
        HelpList.append(True)
    help_check_exit(messages, HelpList)
    header_data[4] = header_data[4] + "\n"
    header_data[5] = header_data[5] + "\n"
    header_data[6] = header_data[6] + "\n"
    # Chargeによって分岐
    if Charge == "0" or Charge == "+0":
        header_data[7] = "0 1\n"
    elif Charge == "+1":
        header_data[7] = "1 2\n"
    if Debug:
        print(f"\n************* {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.gjfのヘッダー *************")
        for formatted_line in header_data:
            print(formatted_line, end="")
        print("********************************************************************\n")
    # GJFファイルへの書き込み
    with open(f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.gjf", "w") as file:
        for formatted_line in header_data:
            file.write(formatted_line)
        for formatted_line in formatted_data:
            file.write(formatted_line + "\n")
        file.write("\n")
        file.write(footer_data)
        file.write("\n")
        file.write("\n")
    if Debug:
        print(f"\n************* {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.gjf *************")
        for formatted_line in header_data:
            print(formatted_line)
        for formatted_line in formatted_data:
            print(formatted_line)
        print(footer_data)
        print("**************************************************************************\n")
    messages.append(f"\t>>> {MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.gjf has been created.")
    help_check_exit(messages, HelpList)
    return None


# shファイルの作成関数
def write_sh_file(messages, HelpList, EG_or_SP, Charge, Operator, GJF_Name, basis_function):
    """
    shファイルの作成関数
    :param basis_function:
    :param GJF_Name:
    :param Operator:
    :param Charge:
    :param EG_or_SP:
    :param HelpList:
    :param messages:
    """
    SH_Name = f"G-{Operator}_{Charge}_{EG_or_SP}_{basis_function[1]}.sh"
    print(f"\nCreating {SH_Name}...")
    # shファイルの作成
    Sh_txt = Stereotyped.Sh_txt.splitlines()
    Sh_txt[12] = f"g16 {GJF_Name}\n"
    with open(f"{SH_Name}", "w") as file:
        for line in Sh_txt:
            file.write(line + "\n")
    messages.append(f"\t>>> {SH_Name} has been created.")
    help_check_exit(messages, HelpList)
    qsub_temp = f"qsub {SH_Name}"
    return qsub_temp


# ジョブを実行する関数
def job_execution(qsub_temp, Debug):
    """
    ジョブを実行する関数
    :param Debug:
    :param qsub_temp:
    """
    print(f"\n{Color.BLUE}Jobs are submitting...{Color.RESET}")
    if qsub_temp == "":
        print("\t>>> Any job was not submitted. Calculations with the conditions might be finished.")
    else:
        subprocess.run(qsub_temp.split())
        jobID = get_job_id()
        RunningList = Running_JobIDList()
        Flag, wait_job_count = check_jobs(RunningList, [jobID])
        start_time = datetime.datetime.now()
        if Flag:
            formated_ST = start_time.strftime("%m/%d %H:%M:%S")
            print(f"\n'{int(len([jobID]))}' calculations was submitted!! at {formated_ST}")
            print(f"\n"
                  f"Wait until jobID {[jobID][-1]}!!\n"
                  f"start ID: {RunningList[0]}\n"
                  f"\n")
        else:
            pass

        while Flag:
            Flag, wait_job_count = check_jobs(RunningList, [jobID])
            formated_NOW, elapsed_time = getElapsedTime(start_time)
            Now_time = datetime.datetime.now()
            formated_end_time = (
                (Now_time + datetime.timedelta(minutes=wait_job_count)).strftime("%m/%d %H:%M:%S"))
            sys.stdout.write(
                "\033[1F\033[G%s" %
                f"\t{formated_NOW} ({elapsed_time} min. passed): Job {RunningList[0]} is in progress.            \n"
                f"\tNext Check >>> {wait_job_count} minute later! forecast: {formated_end_time}    ")
            sys.stdout.flush()
            time.sleep(wait_job_count * 60)
            RunningList = Running_JobIDList()
            Flag, wait_job_count = check_jobs(RunningList, [jobID])
        else:
            print(f"{Color.GREEN}\n\nCalculation cycles until JobID {jobID} were finished.{Color.RESET}")
    if Debug:
        pass
    else:
        rmWildCards(f"./*.chk")
    return None


# jobIDを取得する関数
def get_job_id():
    """
    jobIDを取得する関数
    """
    result = subprocess.run(["qstat"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    lines = result.stdout.strip().splitlines()
    lastJobID = int(lines[-1].strip().split()[0])
    return lastJobID


# 現在実行中のジョブIDのリストを取得する関数
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


# ジョブの状態をチェックする関数
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


# 経過時間を取得する関数
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


# ワイルドカードを削除する関数
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


# Frequencyにマイナスがないかを確認する関数
def check_negative_frequency(LogFileName, messages, HelpList, Debug, Addition_Information):
    """
    Frequencyにマイナスがないかを確認する関数
    :param Addition_Information:
    :param Debug:
    :param LogFileName:
    :param HelpList:
    :param messages:
    """
    count = 0
    print(f"\nChecking if there are any negative frequencies in {LogFileName}...")
    if Debug:
        Addition_Information.append(f"\nfrequencies in {LogFileName}\n")
    with open(LogFileName, "r") as file:
        data = file.read()
    Temp = data.split("Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering")[1]
    for line in Temp.strip().splitlines():
        if "Frequencies" in line:
            elements = line.split()
            for element in elements:
                try:
                    element = float(element)
                    count = count + 1
                    if Debug:
                        print(f"\t>>> Frequency {count}: {element}")
                        Addition_Information.append(f"Frequency {count}: {element}\n")
                    if element < 0:
                        messages.append(f"\t{Color.RED}>>> Negative frequency found {count}: {element}{Color.RESET}")
                        HelpList.append(True)
                except ValueError:
                    pass
    help_check_exit(messages, HelpList)
    print(f"\t{Color.GREEN}>>> No negative frequencies found.{Color.RESET}")
    return None


# 計算結果からエネルギーを取得し、再配置エネルギーを計算する関数
def Calculate_ReorgEnergy(MolecularName, EnergyList, Debug, basis_function):
    """
    エネルギーを取得し、再配置エネルギーを計算する関数
    :param basis_function:
    :param MolecularName:
    :param Debug:
    :param EnergyList:
    """
    print("\nGetting the energy...")
    # 0_EGのエネルギーを取得
    get_energy(MolecularName, "0", "EG", EnergyList, basis_function)
    # +1_EGのエネルギーを取得
    get_energy(MolecularName, "+1", "EG", EnergyList, basis_function)
    # +0_SPのエネルギーを取得
    get_energy(MolecularName, "+0", "SP", EnergyList, basis_function)
    # +1_SPのエネルギーを取得
    get_energy(MolecularName, "+1", "SP", EnergyList, basis_function)

    # エネルギーリストの表示
    # カラムのヘッダー
    headers = ["Charge", "EG_or_SP", "Energy"]
    # 各カラムの幅を求める
    column_widths = {header: max(len(str(row[header])) for row in EnergyList + [dict(zip(headers, headers))])
                     for header in headers}
    # ヘッダー行の表示
    header_row = " | ".join(f"{header:<{column_widths[header]}}" for header in headers)
    print("\t" + "-" * len(header_row))
    print("\t" + header_row)
    print("\t" + "-" * len(header_row))
    # データ行の表示
    for row in EnergyList:
        print("\t" + " | ".join(f"{str(row[header]):<{column_widths[header]}}" for header in headers))
    # 下線の表示
    print("\t" + "-" * len(header_row) + "\n")

    # 再配置エネルギーの計算
    result_meV = reorg_energy(EnergyList, Debug)
    return EnergyList, result_meV


# エネルギーを取得する関数
def get_energy(MolecularName, Charge, EG_or_SP, Energy_List, basis_function):
    """
    エネルギーを取得する関数
    :param basis_function:
    :param EG_or_SP:
    :param Charge:
    :param MolecularName:
    :param Energy_List:
    """
    with open(f"{MolecularName}_{Charge}_{EG_or_SP}_{basis_function[1]}.log", "r") as file:
        data = file.read()
    Temp = data.split("SCF Done:")[-1]
    Temp = Temp.split("=")[1]
    Temp = Temp.split("A.U. after")[0].strip()
    aaa = {
        "Charge": Charge,
        "EG_or_SP": EG_or_SP,
        "Energy": Temp
    }
    Energy_List.append(aaa)
    return None


# 再配置エネルギーを計算する関数
def reorg_energy(EnergyList, Debug):
    """
    再配置エネルギーを計算する関数
    :param Debug:
    :param EnergyList:
    """
    print("Calculating the reorganization energy...")
    # 必要なエネルギー値を抽出
    energy_0_EG = float(
        next(
            item['Energy'] for item in EnergyList
            if item['Charge'] == '0' and item['EG_or_SP'] == 'EG'
        )
    )
    energy_0_SP = float(
        next(
            item['Energy'] for item in EnergyList
            if item['Charge'] == '+0' and item['EG_or_SP'] == 'SP'
        )
    )
    energy_1_EG = float(
        next(
            item['Energy'] for item in EnergyList
            if item['Charge'] == '+1' and item['EG_or_SP'] == 'EG'
        )
    )
    energy_1_SP = float(
        next(
            item['Energy'] for item in EnergyList
            if item['Charge'] == '+1' and item['EG_or_SP'] == 'SP'
        )
    )

    # 式に従って計算
    result_Hartree = ((energy_0_SP - energy_0_EG) + (energy_1_SP - energy_1_EG))
    if Debug:
        print(f"\t>>> λ1 (0_SP - 0_EG): {energy_0_SP - energy_0_EG}")
        print(f"\t>>> λ2 (1_SP - 1_EG): {energy_1_SP - energy_1_EG}")
        print("\t>>> 再配置エネルギー[Hartree]:", result_Hartree, "[Hartree]")
    # HartreeをmeVに変換
    result_meV = hartree_to_meV(result_Hartree)
    print(f"\t{Color.GREEN}>>> 再配置エネルギー[meV]:", result_meV, f"[meV]{Color.RESET}")
    return result_meV


# HartreeをmeVに変換する関数
def hartree_to_meV(Hartree):
    """
    HartreeをmeVに変換する関数
    :param Hartree:
    """
    meV = Hartree * 1000 * Constants.Hartree_to_eV
    return meV


# 計算結果をファイルに保存する関数
def save_results(MolecularName, EnergyList, ReorgEnergy, Addition_Information, Operator, basis_function, args):
    """
    計算結果をファイルに保存する関数
    :param args:
    :param basis_function:
    :param Operator:
    :param Addition_Information:
    :param MolecularName:
    :param EnergyList:
    :param ReorgEnergy:
    """
    print("\nSaving the results...")
    # テキストファイルに整形して書き込み
    if args.function:
        Text_File_Name = f"{MolecularName}_{basis_function[1]}_ReorgEnergy.txt"
    else:
        Text_File_Name = f"{MolecularName}_ReorgEnergy.txt"
    with open(Text_File_Name, "w") as file:
        file.write(f"{'Execution date':<20}{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                   f"{'Molecular Name':<20}{MolecularName}\n"
                   f"{'Operator name':<20}{Operator}\n"
                   f"{'Basis function':<20}{basis_function[0]}\n"
                   f"\n"
                   f"{'ReorgEnergy':<20}{ReorgEnergy} [meV]\n"
                   f"\n")

        file.write("=" * 40 + "\n")
        file.write(f"{'Charge':<10}{'EG/SP':<10}{'Energy [Hartree]':<20}\n")
        file.write("=" * 40 + "\n")
        for item in EnergyList:
            file.write(f"{item['Charge']:<10}{item['EG_or_SP']:<10}{item['Energy']:<20}\n")
        file.write("=" * 40 + "\n")

        # 追加情報
        for Addition in Addition_Information:
            file.write(Addition)
        file.write("\n")
    print(f"\t>>> Results saved to {Color.GREEN}{Text_File_Name}{Color.RESET}")
    return None


if __name__ == '__main__':
    main()
