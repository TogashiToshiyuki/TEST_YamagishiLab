#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import functools
import glob
import math
import os
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import differential_evolution

print = functools.partial(print, flush=True)


def main():
    print(StandardPhrases.ProgramAbst)
    args, before = get_args()
    EM = EffectiveMass(args)

    # data filesを読み込むor作成するかを数字で選択する
    while True:
        make_dat_flag = input("Do you want to create a new data file? (y/n): ")
        if make_dat_flag == "y":
            DatList = EM.make_dat_files()
            break
        elif make_dat_flag == "n":
            DatList = glob.glob("./BandInfo/*.dat")
            if not DatList:
                print(f"Dat file does not exist, please place dat file in BandInfo directory.")
                EM.HelpList.append(True)
                EM.help_check_exit()
                break
            else:
                print(f"Dat Files:")
                for dat in DatList:
                    print(f"\t{dat}")
                break

    # Make band structure figures
    print("\n**********\nMaking band structures...")
    Pattern = Constants.Pattern
    FailFiles = []
    MassValues = []
    plotPNGName = ""
    for Dat in DatList:
        print(f"\n>>> Creating band structure figures for {Color.GREEN}{Dat}{Color.RESET}...")
        Params, Fail = EM.getParameters(Dat)
        if Fail:
            FailFiles.append(Dat)
        else:
            if EM.debug:
                EM.displayDef(Params)
            dev, Kcol, Ktrv, Energy_plus_array, Energy_minus_array, Mass, Masses = EM.calcEffMass(Params)
            plotPNGName = EM.plotBandDisp(Params["Comment"], dev, Kcol, Ktrv, Energy_plus_array, Energy_minus_array,
                                          Masses, Pattern)


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
            print(f"{Color.RED}{StandardPhrases.AbnormalEnd}{Color.RESET}")
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

    def getParameters(self, path):
        if not os.path.isfile(path):
            print(f"\t\t{Color.RED}>>> Error: {path} DOES NOT exist in the current directory.{Color.RESET}")
            self.HelpList.append(True)
            self.help_check_exit()
        else:
            pass
        with open(path, "r") as f:
            lines = f.readlines()
        if lines[4].isspace() and lines[10].isspace():
            Comment = lines[0].strip()
            dev = int(lines[1].strip())
            Dcol = round(float(lines[2].strip()), 3)
            Dtrv = round(float(lines[3].strip()), 3)
            TI12 = round(float(lines[5].strip()), 3)
            TI13 = round(float(lines[6].strip()), 3)
            TI23 = round(float(lines[7].strip()), 3)
            TI34 = round(float(lines[8].strip()), 3)
            TI35 = round(float(lines[9].strip()), 3)
            Fail = False
        else:
            print(f"\t\t{Color.RED}>>> Error: {path} is not a valid file.{Color.RESET}")
            Fail = True
            Comment, dev, Dcol, Dtrv, TI12, TI13, TI23, TI34, TI35 = "", 0, 0, 0, 0, 0, 0, 0, 0
        Params = {"Comment": Comment, "dev": dev, "Dcol": Dcol, "Dtrv": Dtrv, "TI12": TI12, "TI13": TI13, "TI23": TI23,
                  "TI34": TI34, "TI35": TI35}
        return Params, Fail

    @staticmethod
    def displayDef(Params):
        Dcol = Params["Dcol"]
        Dtrv = Params["Dtrv"]
        TI12 = Params["TI12"]
        TI13 = Params["TI13"]
        TI23 = Params["TI23"]
        TI34 = Params["TI34"]
        TI35 = Params["TI35"]
        Comment = Params["Comment"]

        print(f"\t*************** Definition of Variables ***************\n"
              f"\n"
              f"\t                       –Dtransv–\n"
              f"\t(Mol1) ––––––––––––––––––––––––––––––––––––– (Mol4) \n"
              f"\t  |     .                                  .   |\n"
              f"\t  |         . TI13               TI34  .       |\n"
              f"\t  |             .                 .            |\n"
              f"\t  |                 .        .                 |   |\n"
              f"\t TI12                 (Mol3)                   |  Dcol\n"
              f"\t  |                 .        .                 |   |\n"
              f"\t  |             .                 .            |\n"
              f"\t  |        .  TI23              TI35   .       |\n"
              f"\t  |    .                                   .   |\n"
              f"\t(Mol2) ––––––––––––––––––––––––––––––––––––– (Mol5) \n"
              f"\n"
              f"\t  Data from '{Comment}'\n"
              f"\t  Distance in col. = {Dcol} Å\n"
              f"\t  Distance in transv. = {Dtrv} Å\n"
              f"\t  Transfer integrals:\n"
              f"\t       TI12 (Mol1-Mol2) = {TI12} meV\n"
              f"\t       TI13 (Mol1-Mol3) = {TI13} meV\n"
              f"\t       TI23 (Mol2-Mol3) = {TI23} meV\n"
              f"\t       TI34 (Mol3-Mol4) = {TI34} meV\n"
              f"\t       TI35 (Mol3-Mol5) = {TI35} meV\n"
              f"\t*******************************************************\n")
        return None

    def calcEffMass(self, Params):
        print(f"\tCalculating effective masses...")
        Comment = Params["Comment"]
        dev = Params["dev"]
        Dcol = Params["Dcol"]
        Dtrv = Params["Dtrv"]
        TI12 = Params["TI12"]
        TI13 = Params["TI13"]
        TI23 = Params["TI23"]
        TI34 = Params["TI34"]
        TI35 = Params["TI35"]

        # k空間の範囲を定義
        Kcol_range = [-math.pi / Dcol, math.pi / Dcol]
        Ktrv_range = [-math.pi / Dtrv, math.pi / Dtrv]

        # エネルギーを最適化する関数を定義
        def energy_func(k):
            Kc, Kt = k
            results = self.calcEnergy(Dcol, Dtrv, TI12, TI13, TI23, TI34, TI35, Kc, Kt)
            if "HOMO" in Comment:
                return -results[0]  # 最大化するためにマイナスを付ける
            elif "LUMO" in Comment:
                return results[1]  # 最小化
            else:
                print("Error: Comment should contain HOMO or LUMO")
                return None

        # 制約条件を設定（kの範囲内に収める）
        bounds = [Kcol_range, Ktrv_range]

        # 最適化を実行
        res = differential_evolution(energy_func, bounds)

        if not res.success:
            print(f"\t{Color.RED}>>> Optimization failed.{Color.RESET}")
            return None, None, None, None, None, None, None

        # 最適化結果からkベクトルを取得
        Kc_be, Kt_be = res.x

        # エネルギー面を可視化
        self.plot_EnergySurface(Dcol, Dtrv, TI12, TI13, TI23, TI34, TI35, Comment, Kc_be, Kt_be)

        # バンド端のエネルギーを計算
        results = self.calcEnergy(Dcol, Dtrv, TI12, TI13, TI23, TI34, TI35, Kc_be, Kt_be)
        if "HOMO" in Comment:
            Energy = results[0]
            d2E_dKc2 = results[2]
            d2E_dKt2 = results[3]
            d2E_dKcKt = results[4]
        elif "LUMO" in Comment:
            Energy = results[1]
            d2E_dKc2 = results[5]
            d2E_dKt2 = results[6]
            d2E_dKcKt = results[7]
        else:
            Energy, d2E_dKc2, d2E_dKt2, d2E_dKcKt = 0, 0, 0, 0

        print(f"\t\t>>> Energy of the band edge: {Energy} meV at ({round(Kc_be, 6)}, {round(Kt_be, 6)})")

        # 有効質量テンソルを構築
        Mass_tensor = np.array([[d2E_dKc2, d2E_dKcKt],
                                [d2E_dKcKt, d2E_dKt2]])

        # 有効質量の計算
        try:
            # テンソルの対角化
            eigenvalues, eigenvectors = np.linalg.eig(Mass_tensor)
            Mass = Constants.h_bar ** 2 / eigenvalues / Constants.ElMass * 1e20
            Masses = [Mass[0], Mass[1], "from tensor"]
            if self.debug:
                print(f"\t\t>>> d2E_dKc2: {d2E_dKc2}")
                print(f"\t\t>>> d2E_dKt2: {d2E_dKt2}")
                print(f"\t\t>>> d2E_dKcKt: {d2E_dKcKt}")
            print(f"\t>>> Mass in {Color.BLUE}Column{Color.RESET}: {Mass[0]} m0 "
                  f"{Color.GREEN}from tensor{Color.RESET}")
            print(f"\t>>> Mass in {Color.BLUE}Transverse{Color.RESET}: {Mass[1]} m0 "
                  f"{Color.GREEN}from tensor{Color.RESET}")
        except Exception as e:
            print(f"{Color.RED}\t>>> Error: {e}{Color.RESET}")
            Mass_col = Constants.h_bar ** 2 / d2E_dKc2 / Constants.ElMass * 1e20
            Mass_trv = Constants.h_bar ** 2 / d2E_dKt2 / Constants.ElMass * 1e20
            if self.debug:
                print(f"\t\t>>> d2E_dKc2: {d2E_dKc2}")
                print(f"\t\t>>> d2E_dKt2: {d2E_dKt2}")
                print(f"\t\t>>> d2E_dKcKt: {d2E_dKcKt}")
            print(f"\t>>> Mass in {Color.BLUE}Column{Color.RESET}: {Mass_col} m0 "
                  f"from {Color.RED}NOT{Color.RESET} tensor")
            print(f"\t>>> Mass in {Color.BLUE}Transverse{Color.RESET}: {Mass_trv} m0 "
                  f"from {Color.RED}NOT{Color.RESET} tensor")
            Mass = np.array([[], []])
            Masses = [Mass_col, Mass_trv, "not from tensor"]

        # バンド図のためのデータを返す
        Kcol_for_plot = round(math.pi / Dcol / dev, 10)
        Ktrv_for_plot = round(math.pi / Dtrv / dev, 10)
        Energy_plus_array, Energy_minus_array = np.zeros((dev + 1, dev + 1)), np.zeros((dev + 1, dev + 1))
        for i in range(dev + 1):
            for j in range(dev + 1):
                Kc = Kcol_for_plot * i
                Kt = Ktrv_for_plot * j
                Energy_plus, Energy_minus, *_ = self.calcEnergy(Dcol, Dtrv, TI12, TI13, TI23, TI34, TI35, Kc, Kt)
                Energy_plus_array[i][j] = Energy_plus
                Energy_minus_array[i][j] = Energy_minus

        return dev, Kcol_for_plot, Ktrv_for_plot, Energy_plus_array, Energy_minus_array, Mass, Masses

    @staticmethod
    def calcEnergy(column, transv, TI12, TI13, TI23, TI34, TI35, Kcol, Ktrv):
        B11 = 2 * TI12 * math.cos(Kcol * column)
        B12R = ((TI13 + TI35) * math.cos(Kcol * (column / 2) + Ktrv * (transv / 2))
                + (TI23 + TI34) * math.cos(Kcol * (column / 2) - Ktrv * (transv / 2)))
        B12I = ((TI35 - TI13) * math.sin(Kcol * (column / 2) + Ktrv * (transv / 2))
                + (TI23 - TI34) * math.sin(Kcol * (column / 2) - Ktrv * (transv / 2)))
        B12 = math.sqrt(B12R ** 2 + B12I ** 2)
        # Avoid division by zero
        epsilon = 1e-20
        B12 = max(B12, epsilon)

        Energy_plus = B11 + B12
        Energy_minus = B11 - B12

        # Kcolで一階微分
        # dB11_dKc = -2 * TI12 * column * math.sin(column * Kcol)
        dB12R_dKc = -1 * (column / 2) * (
                (TI13 + TI35) * math.sin(((column / 2) * Kcol) + ((transv / 2) * Ktrv)) +
                (TI23 + TI34) * math.sin(((column / 2) * Kcol) - ((transv / 2) * Ktrv))
        )
        dB12I_dKc = (column / 2) * (
                (TI35 - TI13) * math.cos(((column / 2) * Kcol) + ((transv / 2) * Ktrv)) +
                (TI23 - TI34) * math.cos(((column / 2) * Kcol) - ((transv / 2) * Ktrv))
        )
        # dB12_dKc = (1 / B12) * (dB12R_dKc * B12R + dB12I_dKc * B12I)
        # dEplus_dKc = dB11_dKc + dB12_dKc
        # dEminus_dKc = dB11_dKc - dB12_dKc

        # Ktrvで一階微分
        # dB11_dKt = 0
        dB12R_dKt = -1 * (transv / 2) * (
                (TI13 + TI35) * math.sin(((column / 2) * Kcol) + ((transv / 2) * Ktrv)) -
                (TI23 + TI34) * math.sin(((column / 2) * Kcol) - ((transv / 2) * Ktrv))
        )
        dB12I_dKt = (transv / 2) * (
                (TI35 - TI13) * math.cos(((column / 2) * Kcol) + ((transv / 2) * Ktrv)) -
                (TI23 - TI34) * math.cos(((column / 2) * Kcol) - ((transv / 2) * Ktrv))
        )
        # dB12_dKt = (1 / B12) * (dB12R_dKt * B12R + dB12I_dKt * B12I)
        # dEplus_dKt = dB11_dKt + dB12_dKt
        # dEminus_dKt = dB11_dKt - dB12_dKt

        # Kcolで二階微分
        d2B11_dKc2 = -2 * TI12 * (column ** 2) * math.cos(column * Kcol)
        d2B12R_dKc2 = -(column / 2) ** 2 * (
                (TI13 + TI35) * math.cos((column / 2) * Kcol + (transv / 2) * Ktrv) +
                (TI23 + TI34) * math.cos((column / 2) * Kcol - (transv / 2) * Ktrv)
        )
        d2B12I_dKc2 = -(column / 2) ** 2 * (
                (TI35 - TI13) * math.sin((column / 2) * Kcol + (transv / 2) * Ktrv) +
                (TI23 - TI34) * math.sin((column / 2) * Kcol - (transv / 2) * Ktrv)
        )
        d2B12_dKc2 = (1 / B12) * (
                (dB12R_dKc ** 2) + (dB12I_dKc ** 2) + B12R * d2B12R_dKc2 + B12I * d2B12I_dKc2
        ) - ((1 / B12) ** 3) * (
                             (B12R * dB12R_dKc + B12I * dB12I_dKc) ** 2
                     )
        d2Eplus_dKc2 = d2B11_dKc2 + d2B12_dKc2
        d2Eminus_dKc2 = d2B11_dKc2 - d2B12_dKc2

        # Ktrvで二階微分
        d2B11_dKt2 = 0
        d2B12R_dKt2 = -(transv / 2) ** 2 * (
                (TI13 + TI35) * math.cos((column / 2) * Kcol + (transv / 2) * Ktrv) +
                (TI23 + TI34) * math.cos((column / 2) * Kcol - (transv / 2) * Ktrv)
        )
        d2B12I_dKt2 = -(transv / 2) ** 2 * (
                (TI35 - TI13) * math.sin((column / 2) * Kcol + (transv / 2) * Ktrv) +
                (TI23 - TI34) * math.sin((column / 2) * Kcol - (transv / 2) * Ktrv)
        )
        d2B12_dKt2 = (1 / B12) * ((dB12R_dKt ** 2) + (dB12I_dKt ** 2) + B12R * d2B12R_dKt2 + B12I * d2B12I_dKt2) - (
                (1 / B12) ** 3) * (
                             (B12R * dB12R_dKt + B12I * dB12I_dKt) ** 2
                     )
        d2Eplus_dKt2 = d2B11_dKt2 + d2B12_dKt2
        d2Eminus_dKt2 = d2B11_dKt2 - d2B12_dKt2

        # KcolとKtrvでの混合微分
        d2B11_dKcKt = 0
        d2B12R_dKcdKt = -1 * (column / 2) * (transv / 2) * (
                (TI13 + TI35) * math.cos((column / 2) * Kcol + (transv / 2) * Ktrv) -
                (TI23 + TI34) * math.cos((column / 2) * Kcol - (transv / 2) * Ktrv)
        )
        d2B12I_dKcdKt = -1 * (column / 2) * (transv / 2) * (
                (TI35 - TI13) * math.sin((column / 2) * Kcol + (transv / 2) * Ktrv) -
                (TI23 - TI34) * math.sin((column / 2) * Kcol - (transv / 2) * Ktrv)
        )
        d2B12_dKcKt = (1 / B12) * (
                dB12R_dKc * dB12R_dKt + B12R * d2B12R_dKcdKt + dB12I_dKc * dB12I_dKt + B12I * d2B12I_dKcdKt) - (
                              ((1 / B12) ** 3) * (B12R * dB12R_dKc + B12I * dB12I_dKc) *
                              (B12R * dB12R_dKt + B12I * dB12I_dKt)
                      )
        dEplus_dKcKt = d2B11_dKcKt + d2B12_dKcKt
        dEminus_dKcKt = d2B11_dKcKt - d2B12_dKcKt

        return (Energy_plus, Energy_minus,
                d2Eplus_dKc2, d2Eplus_dKt2, dEplus_dKcKt, d2Eminus_dKc2, d2Eminus_dKt2, dEminus_dKcKt)

    def plot_EnergySurface(self, Dcol, Dtrv, TI12, TI13, TI23, TI34, TI35, Comment, Kc_be, Kt_be):
        if self.debug:
            # k空間を細かく分割
            Kc_values = np.linspace(-math.pi / Dcol, math.pi / Dcol, 200)
            Kt_values = np.linspace(-math.pi / Dtrv, math.pi / Dtrv, 200)
            Kc_grid, Kt_grid = np.meshgrid(Kc_values, Kt_values)
            Energy_grid = np.zeros_like(Kc_grid)

            for i in range(len(Kc_values)):
                for j in range(len(Kt_values)):
                    Kc = Kc_values[i]
                    Kt = Kt_values[j]
                    results = self.calcEnergy(Dcol, Dtrv, TI12, TI13, TI23, TI34, TI35, Kc, Kt)
                    if "HOMO" in Comment:
                        Energy_grid[j, i] = results[0]
                    elif "LUMO" in Comment:
                        Energy_grid[j, i] = results[1]

            plt.figure()
            plt.contourf(Kc_grid, Kt_grid, Energy_grid, levels=50, cmap='viridis')
            plt.colorbar(label='Energy (meV)')
            plt.xlabel('Kc')
            plt.ylabel('Kt')
            plt.title(f'Energy Surface for {Comment}')

            # 極値点をプロット
            plt.plot(Kc_be, Kt_be, '.', markersize=10, label='Band Edge', color='red')
            plt.legend()
            plt.show()
            plt.close()
        else:
            pass
        return None

    def plotBandDisp(self, Comment, dev, Kcol, Ktrv, Energy_plus_array, Energy_minus_array, Masses, Pattern):
        # Make plot for Energy dispersion
        ZeroTo = np.arange(0, dev + 1, 1)
        ToZero = np.flipud(ZeroTo)
        Zero = np.zeros(dev + 1)
        Top = np.full(dev + 1, dev)

        if "1" in Pattern:
            xdev = []
            y1 = []
            y2 = []
            xUnitLen = Ktrv
            x_dev_temp, y1_temp, y2_temp = self.EnergyPlotParts(dev, Energy_plus_array, Energy_minus_array,
                                                                xUnitLen, Top, ToZero)
            xdev.extend(x_dev_temp)
            y1.extend(y1_temp)
            y2.extend(y2_temp)

            xUnitLen = Kcol
            x_dev_temp, y1_temp, y2_temp = self.EnergyPlotParts(dev, Energy_plus_array, Energy_minus_array,
                                                                xUnitLen, ToZero, Zero)
            xdev.extend(x_dev_temp)
            y1.extend(y1_temp)
            y2.extend(y2_temp)

            xUnitLen = Ktrv
            x_dev_temp, y1_temp, y2_temp = self.EnergyPlotParts(dev, Energy_plus_array, Energy_minus_array,
                                                                xUnitLen, Zero, ZeroTo)
            xdev.extend(x_dev_temp)
            y1.extend(y1_temp)
            y2.extend(y2_temp)

            xUnitLen = Kcol
            x_dev_temp, y1_temp, y2_temp = self.EnergyPlotParts(dev, Energy_plus_array, Energy_minus_array,
                                                                xUnitLen, ZeroTo, Top)
            xdev.extend(x_dev_temp)
            y1.extend(y1_temp)
            y2.extend(y2_temp)

            xUnitLen = math.sqrt(Kcol ** 2 + Ktrv ** 2)
            x_dev_temp, y1_temp, y2_temp = self.EnergyPlotParts(dev, Energy_plus_array, Energy_minus_array,
                                                                xUnitLen, ToZero, ToZero)
            xdev.extend(x_dev_temp)
            y1.extend(y1_temp)
            y2.extend(y2_temp)

            x, val = [], 0
            for i in range(len(xdev)):
                if i == 0:
                    val = 0 + xdev[i]
                    x.append(val)
                else:
                    val = val + xdev[i]
                    x.append(val)

            point1 = x[0]  # C point: pi/a,pi/b
            point2 = x[dev + 1]  # X point: pi/a,0
            point3 = x[2 * (dev + 1)]  # Gamma point: 0.0, 0.0
            point4 = x[3 * (dev + 1)]  # Y point: 0,pi/b
            point5 = x[4 * (dev + 1)]  # C point: pi/a, pi/b
            point6 = x[-1]  # Gamma point: 0.0, 0.0
            xtickvals = [point1, point2, point3, point4, point5, point6]
            xticktexts = ["C", "X", "gamma", "Y", "C", "gamma"]

        elif "2" in Pattern:
            xdev = []
            y1 = []
            y2 = []

            xUnitLen = Kcol
            x_dev_temp, y1_temp, y2_temp = self.EnergyPlotParts(dev, Energy_plus_array, Energy_minus_array,
                                                                xUnitLen, ToZero, Zero)
            xdev.extend(x_dev_temp)
            y1.extend(y1_temp)
            y2.extend(y2_temp)

            xUnitLen = Ktrv
            x_dev_temp, y1_temp, y2_temp = self.EnergyPlotParts(dev, Energy_plus_array, Energy_minus_array,
                                                                xUnitLen, Zero, ZeroTo)
            xdev.extend(x_dev_temp)
            y1.extend(y1_temp)
            y2.extend(y2_temp)

            xUnitLen = Kcol
            x_dev_temp, y1_temp, y2_temp = self.EnergyPlotParts(dev, Energy_plus_array, Energy_minus_array,
                                                                xUnitLen, ZeroTo, Top)
            xdev.extend(x_dev_temp)
            y1.extend(y1_temp)
            y2.extend(y2_temp)

            xUnitLen = Ktrv
            x_dev_temp, y1_temp, y2_temp = self.EnergyPlotParts(dev, Energy_plus_array, Energy_minus_array,
                                                                xUnitLen, Top, ToZero)
            xdev.extend(x_dev_temp)
            y1.extend(y1_temp)
            y2.extend(y2_temp)

            x, val = [], 0
            for i in range(len(xdev)):
                if i == 0:
                    val = 0 + xdev[i]
                    x.append(val)
                else:
                    val = val + xdev[i]
                    x.append(val)

            # print(len(x))
            point1 = x[0]  # X point: pi/a,0
            point2 = x[dev + 1]  # Gamma point: 0.0, 0.0
            point3 = x[2 * (dev + 1)]  # Y point: 0,pi/b
            point4 = x[3 * (dev + 1)]  # C point: pi/a, pi/b
            point5 = x[-1]  # X point: pi/a,0
            xtickvals = [point1, point2, point3, point4, point5]
            xticktexts = ["X", "gamma", "Y", "C", "X"]

        else:
            print(f"\t{Color.RED}>>> Error: Pattern should be 1 or 2.{Color.RESET}")
            self.HelpList.append(True)
            self.help_check_exit()
            exit()

        plt.figure()
        plt.plot(x, y1, color="r")
        plt.plot(x, y2, color="r")
        for xtickval in xtickvals:
            plt.axvline(xtickval, color="k", linewidth=0.75)
        plt.xticks(xtickvals, xticktexts)
        plt.tick_params(axis="both", direction="in", labelsize=14)
        plt.xlim(xtickvals[0], xtickvals[-1])
        plt.ylabel("E (meV)", fontsize=18)
        plt.title(f"{Comment}: mc={round(Masses[0], 2)}m0, mt={round(Masses[1], 2)}m0 '{Masses[2]}'", fontsize=10)
        plt.tight_layout()
        if self.debug:
            plt.show()
            pass
        else:
            pass

        plotPNGname = Comment[:Comment.find(",")].strip()
        plotPNGname = plotPNGname[:-2]
        plt.savefig(f"./Figures/BandStructures/{plotPNGname}-{Pattern.strip()}.png", format="png",
                    dpi=500, bbox_inches="tight")
        plt.close()

        with open(f"./BandInfo/{plotPNGname}-{Pattern.strip()}.txt", "w") as f:
            for i in range(len(x)):
                f.write(f"{x[i]}\t{y1[i]}\t{y2[i]}\n")

        return plotPNGname

    @staticmethod
    def EnergyPlotParts(dev, E_plus_array, E_minus_array, xUnitLen, D_col, D_trv):
        x_dev_temp = np.full(dev + 1, xUnitLen)
        x_dev_temp[0] = 0
        y1_temp, y2_temp = [], []
        for i in range(dev + 1):
            y1_temp.append(E_plus_array[int(D_col[i]), int(D_trv[i])])
            y2_temp.append(E_minus_array[int(D_col[i]), int(D_trv[i])])
        return x_dev_temp, y1_temp, y2_temp


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
