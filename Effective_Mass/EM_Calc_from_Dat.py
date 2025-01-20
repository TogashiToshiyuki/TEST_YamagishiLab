#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import datetime
import functools
import glob
import math
import os
import time

import matplotlib.pyplot as plt
import numpy as np

print = functools.partial(print, flush=True)


def main():
    print(StandardPhrases.ProgramAbst)
    args, before = get_args()
    EM = EffectiveMass(args)


def get_args():
    before = time.time()
    help_text = StandardPhrases.ProgramAbst

    parser = argparse.ArgumentParser(description=help_text, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--debug', '-d', '-D',
                        help="Start in debug mode",
                        action='store_true')
    args = parser.parse_args()

    if args.debug:
        print(f"{Color.RED}*************** Caution!!! Debug Started!!! ***************{Color.RESET}\n")

    return args, before


class EffectiveMass:
    def __init__(self, args):
        self.args = args
        self.debug = args.debug
        self.messages = []
        self.HelpList = []
        self.dat_files = []
        os.makedirs("./BandInfo", exist_ok=True)
        os.makedirs("./Figures/BandStructures", exist_ok=True)
        os.makedirs("./Figures/EvsAngle", exist_ok=True)
        os.makedirs("./Figures/TIvsAngle", exist_ok=True)
        os.makedirs("./Figures/AllData", exist_ok=True)
        os.makedirs("./Figures/MassvsAngles", exist_ok=True)

        # data filesを読み込むor作成するかを数字で選択する
        while True:
            make_dat_flag = input("Do you want to create a new data file? (y/n): ")
            if make_dat_flag == "y":
                self.dat_files = self.make_dat_files()
                break
            elif make_dat_flag == "n":
                self.dat_files = glob.glob("*.dat")
                break

    def make_dat_files(self):
        DatList_temp = []
        while True:
            while True:
                material_name = input(f"\nPlease enter {Color.GREEN}the material name{Color.RESET}: ")
                flag = input(f"Material name: {Color.GREEN}{material_name}{Color.RESET}. "
                             f"Is this correct? (y/n): ")
                if flag == "y":
                    break

            while True:
                structure = input(f"\nPlease enter {Color.GREEN}the structure{Color.RESET}\n"
                                  f"1: B12\n"
                                  f"2: B3\n"
                                  f">>> ")
                if structure == "1":
                    structure = "B12"
                    break
                elif structure == "2":
                    structure = "B3"
                    break

            while True:
                HOMO_LUMO = input(f"\nPlease select HOMO or LUMO.\n"
                                  f"1: HOMO\n"
                                  f"2: LUMO\n"
                                  f">>> ")
                if HOMO_LUMO == "1":
                    HOMO_LUMO = "HOMO"
                    break
                elif HOMO_LUMO == "2":
                    HOMO_LUMO = "LUMO"
                    break

            Dcol = self.get_float_input("the column distance (nm)", "Column distance")
            Dtrv = self.get_float_input("the transverse distance (nm)", "Transverse distance")
            T12 = self.get_float_input("T12", "T12")
            T13 = self.get_float_input("T13", "T13")
            T23 = self.get_float_input("T23", "T23")
            T34 = self.get_float_input("T34", "T34")
            T35 = self.get_float_input("T35", "T35")

            print(f"\nCreating Dat File...")
            title = f"{material_name}-{structure}"
            print(f"\t>>> {title}-HOMO.dat: ", end="")
            with open(f"./BandInfo/{title}-HOMO.dat", "w") as Fhomo:
                Fhomo.write(f"{title}-HOMO\n")
                Fhomo.write(f"{Constants.n}\n")
                Fhomo.write(f"{Dcol}\n")
                Fhomo.write(f"{Dtrv}\n")
                Fhomo.write("\n")
                Fhomo.write(f"{T12}\n")
                Fhomo.write(f"{T13}\n")
                Fhomo.write(f"{T23}\n")
                Fhomo.write(f"{T34}\n")
                Fhomo.write(f"{T35}\n")
                Fhomo.write("\n")
            DatList_temp.append(f"./BandInfo/{title}-HOMO.dat")
            print(f"{Color.GREEN}Complete!{Color.RED}")

            # Datファイル作成を終了するかを確認する
            Flag = input(f"\nFinish create new dat data? (y/n):")
            if Flag == "y":
                break
        return DatList_temp

    @staticmethod
    def get_float_input(message_1, message_2):
        while True:
            value = input(f"\nPlease enter {Color.GREEN}{message_1}{Color.RESET}: ")
            try:
                value = float(value)
            except ValueError:
                print(f"{Color.RED}Error:{Color.RESET} Please enter a number.")
                continue
            flag = input(f"{message_2}: {Color.GREEN}{value}{Color.RESET}. "
                         f"Is this correct? (y/n): ")
            if flag == "y":
                break
        return value


class Color:
    """
    ANSI escape code for color
    """

    def __init__(self):
        self._Color = "Color"

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


class StandardPhrases:
    def __init__(self):
        self._StandardPhrases = "StandardPhrases"

    ProgramAbst = ("\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n"
                   "*   This program collects and processes calculation results of                *\n"
                   "*       optimized aggregate structures of organic semiconductors              *\n"
                   "*       with HerringBone (HB) structures.                                     *\n"
                   "*   It performs the following tasks:                                          *\n"
                   "*       - Calculates effective masses                                         *\n"
                   "*       - Plots band diagrams                                                 *\n"
                   "*       - Generates summary slides of the results                             *\n"
                   "*                                                                             *\n"
                   "*   To run this program, please ensure the following files are available:     *\n"
                   "*       - “MoleculeName”_3molp1_t”N”d_results                                 *\n"
                   "*       - “MoleculeName”_3molp2_t”N”d_results                                 *\n"
                   "*       - “MoleculeName”_3molp3_t”N”d_results                                 *\n"
                   "*                                                                             *\n"
                   "*   Replace “MoleculeName” with the actual name of your molecule,             *\n"
                   "*       and “N” with the appropriate value used in your calculations.         *\n"
                   "*                                                                             *\n"
                   "*   When you start the program, you will be prompted to enter your name,      *\n"
                   "*       - which will be used as the slide creator’s name.                     *\n"
                   "*                                                                             *\n"
                   "*                                 created by Toshiyuki Togashi 2024/11/27     *\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n")
    AbnormalEnd = "\n************* Programme DID NOT terminate successfully. *************\n"
    HelpText = ("\n"
                "The required file may not exist.\n"
                "Please check.")


class Constants:
    # ハートリーからeVに変換するための定数
    Hartree_to_eV = 27.21138602
    # 気体定数[R] (J/(mol・K))
    R = 8.314462618
    # 真空中の光速[m/s]
    C = 299792458
    # 電気素量[C]
    E = 1.602176634e-19
    # meV/(m/s)^2
    ElMass = 0.51099895 * 10 ** 9 / C ** 2
    # meVs
    h_bar = 6.582119569 * 10 ** (-13)
    # 近似のための分割数
    n = 100
    # プロットパターン
    Pattern = "2"

    # スライドサイズ
    # 4:3 (default) 9144000x6858000
    # 16:9 12193200x6858000
    SLIDE_WIDTH, SLIDE_HEIGHT = 9144000, 6858000
    # スライド中心のX、Y座標（左上が原点）
    IMG_CENTER_X, IMG_CENTER_Y = SLIDE_WIDTH / 2, SLIDE_HEIGHT / 2
    # スライドのアスペクト比
    SLIDE_ASPECT_RATIO = SLIDE_WIDTH / SLIDE_HEIGHT


if __name__ == '__main__':
    main()
