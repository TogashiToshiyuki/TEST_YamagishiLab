#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import math
import os

import numpy as np
import sympy as sp


def main():
    args = get_args()
    EM = EffectiveMass(args)

    # Make band information files
    print("\n**********\nMaking band information files...")
    if len(EM.p1Files) != 0 and len(EM.p2Files) != 0:
        DatList_p12 = EM.mkBandInfo_p1p2()
    else:
        DatList_p12 = []
    if len(EM.p3Files) != 0:
        DatList_p3 = EM.mkBandInfo_p3()
    else:
        DatList_p3 = []
    DatList = DatList_p12 + DatList_p3
    DatList.sort()
    print(f"\n\t>>> {Color.GREEN}The band information files are created.{Color.RESET}")

    # Make band structure figures
    print("\n**********\nMaking band structures...")
    Pattern = "2"
    FailFiles = []
    MassValues = []
    for Dat in DatList:
        print(f"\n>>> Creating band structure figures for {Color.GREEN}{Dat}{Color.RESET}...")
        Params, Fail = EM.getParameters(Dat)
        if Fail:
            FailFiles.append(Dat)
        else:
            if EM.debug:
                EM.displayDef(Params)
            EM.calcEffMass(Params, Dat)


def get_args():
    help_text = "This is a help text"

    parser = argparse.ArgumentParser(description=help_text, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('MaterName',
                        help='Material name')
    parser.add_argument('--debug', '-d', '-D',
                        help="Start in debug mode",
                        action='store_true')
    args = parser.parse_args()
    return args


class EffectiveMass:
    def __init__(self, args):
        self.args = args
        self.MaterName = args.MaterName
        self.debug = args.debug
        self.messages = []
        self.HelpList = []
        self.p1Files, self.p2Files, self.p3Files, self.AllFiles, self.Angles, self.Tilt_Angle = self.File_Set_Check()
        os.makedirs("./BandInfo", exist_ok=True)
        os.makedirs("./Figures/BandStructures", exist_ok=True)

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

    def File_Set_Check(self):
        """
        Check if the required files exist.
        :return:
        """
        print(f"{Color.GREEN}Checking the required files...{Color.RESET}")
        p1Files, p2Files, p3Files = self.getMinTIFileName()
        AllFiles = self.getAllFileNames()
        self.help_check_exit()
        print(f"{Color.GREEN}All required files are available.{Color.RESET}")

        print("\nTI files:")
        for p1File in p1Files:
            print(f"\t{p1File}")
        for p2File in p2Files:
            print(f"\t{p2File}")
        for p3File in p3Files:
            print(f"\t{p3File}")

        print("\nAll files:")
        for AllFile in AllFiles:
            print(f"\t{AllFile}")

        Angles = []
        for AllFile in AllFiles:
            Angles_temp = self.getAngle(AllFile)
            Angles = Angles + Angles_temp
        Angles = list(set(Angles))
        Angles.sort()
        if self.debug:
            print(f"\nAngles: {Angles}")

        print(f"\nGetting the tilt angles...")
        Tilt_temp = []
        for p1File in p1Files:
            Tilt_temp.append(p1File.split("/")[1].split("_")[2])
        for p2File in p2Files:
            Tilt_temp.append(p2File.split("/")[1].split("_")[2])
        for p3File in p3Files:
            Tilt_temp.append(p3File.split("/")[1].split("_")[2])
        Tilt_temp = list(set(Tilt_temp))
        if len(Tilt_temp) != 1:
            print(f"{Color.RED}\t>>> Error: Tilt angles are not consistent.{Color.RESET}")
            self.HelpList.append(True)
        else:
            print(f"\t>>> Tilt angles: {Color.GREEN}{Tilt_temp[0]}{Color.RESET}")
        return p1Files, p2Files, p3Files, AllFiles, Angles, Tilt_temp[0]

    def getMinTIFileName(self):
        Dirs = os.listdir("./")
        p1Files, p2Files, p3Files = [], [], []
        for Dir in Dirs:
            if self.MaterName in Dir and "_results" in Dir:
                Files_temp = os.listdir(Dir)
                for File in Files_temp:
                    if "3molp1" in File and "min-TIs-" in File:
                        p1Files.append("./" + Dir + "/" + File)
                    elif "3molp2" in File and "min-TIs-" in File:
                        p2Files.append("./" + Dir + "/" + File)
                    elif "3molp3" in File and "min-TIs-" in File:
                        p3Files.append("./" + Dir + "/" + File)
            else:
                pass

        if len(p1Files) == 0:
            print(f"{Color.RED}\t>>> Error: No files found for 3molp1.{Color.RESET}")
            self.HelpList.append(True)
        elif len(p1Files) != 3:
            print(f"{Color.RED}\t>>> Error: 3molp1 files are not complete.{Color.RESET}")
            self.HelpList.append(True)
        else:
            print(f"\t>>> 3molp1 files: {Color.GREEN}All Available.{Color.RESET}")

        if len(p2Files) == 0:
            print(f"{Color.RED}\t>>> Error: No files found for 3molp2.{Color.RESET}")
            self.HelpList.append(True)
        elif len(p2Files) != 3:
            print(f"{Color.RED}\t>>> Error: 3molp2 files are not complete.{Color.RESET}")
            self.HelpList.append(True)
        else:
            print(f"\t>>> 3molp2 files: {Color.GREEN}All Available.{Color.RESET}")

        if len(p3Files) == 0:
            print(f"{Color.RED}\t>>> Error: No files found for 3molp3.{Color.RESET}")
            self.HelpList.append(True)
        elif len(p3Files) != 3:
            print(f"{Color.RED}\t>>> Error: 3molp3 files are not complete.{Color.RESET}")
            self.HelpList.append(True)
        else:
            print(f"\t>>> 3molp3 files: {Color.GREEN}All Available.{Color.RESET}")

        return p1Files, p2Files, p3Files

    def getAllFileNames(self):
        Dirs = os.listdir("./")
        AllFiles = []
        for Dir in Dirs:
            if self.MaterName in Dir and "_results" in Dir and os.path.isdir(Dir):
                Files_temp = os.listdir(Dir)
                for File in Files_temp:
                    if "3mol" in File and "_all.txt" in File:
                        AllFiles.append("./" + Dir + "/" + File)
                    else:
                        pass
            else:
                pass

        if len(AllFiles) == 0:
            print(f"{Color.RED}\t>>> Error: No files found for 3mol_all.txt{Color.RESET}")
            self.HelpList.append(True)
        else:
            print(f"\t>>> 3mol_all.txt files: {Color.GREEN}All Available.{Color.RESET}")

        return AllFiles

    @staticmethod
    def getAngle(AllFile):
        with open(AllFile, "r") as f:
            lines = f.readlines()
        Angles_temp = []
        for line in lines:
            Contents = line.strip().split()
            try:
                Angles_temp.append(float(Contents[1]))
            except (ValueError, IndexError):
                pass
        return Angles_temp

    def mkBandInfo_p1p2(self):

        p1Data_12, p1Data_23, p1Data_31 = [], [], []
        for p1File in self.p1Files:
            if "-12.txt" in p1File:
                p1Data_12 = self.TextFileToData(p1File)
            elif "-23.txt" in p1File:
                p1Data_23 = self.TextFileToData(p1File)
            elif "-31.txt" in p1File:
                p1Data_31 = self.TextFileToData(p1File)

        p2Data_12, p2Data_23, p2Data_31 = [], [], []
        for p2File in self.p2Files:
            if "-12.txt" in p2File:
                p2Data_12 = self.TextFileToData(p2File)
            elif "-23.txt" in p2File:
                p2Data_23 = self.TextFileToData(p2File)
            elif "-31.txt" in p2File:
                p2Data_31 = self.TextFileToData(p2File)

        DatList_temp = []

        for Angle in self.Angles:
            name1, name2 = [], []
            for line in p1Data_12:
                DataList = line.strip().split()
                if float(DataList[1]) == Angle:
                    name1 = DataList[0].split("_")
                    Dcol = float(DataList[2].strip())
                    Dtrv1 = float(DataList[3].strip())
                    T12_HOMO = float(DataList[15].strip())
                    T12_LUMO = float(DataList[13])

            for line in p1Data_23:
                DataList = line.split()
                if float(DataList[1]) == Angle:
                    T23_HOMO = float(DataList[15].strip())
                    T23_LUMO = float(DataList[13].strip())
                else:
                    pass

            for line in p1Data_31:
                DataList = line.split()
                if float(DataList[1]) == Angle:
                    T13_HOMO = float(DataList[15].strip())
                    T13_LUMO = float(DataList[13].strip())
                else:
                    pass

            for line in p2Data_12:
                DataList = line.strip().split()
                if float(DataList[1]) == Angle:
                    name2 = DataList[0].split("_")
                    Dtrv2 = float(DataList[3].strip())
                else:
                    pass

            for line in p2Data_23:
                DataList = line.split()
                if float(DataList[1]) == Angle:
                    T35_HOMO = float(DataList[15].strip())
                    T35_LUMO = float(DataList[13].strip())
                else:
                    pass

            for line in p2Data_31:
                DataList = line.split()
                if float(DataList[1]) == Angle:
                    T34_HOMO = float(DataList[15].strip())
                    T34_LUMO = float(DataList[13].strip())
                else:
                    pass

            # 同じ物質のデータかどうかを確認
            if name1[0] == name2[0]:
                name = name1
            else:
                self.messages.append(f"{Color.RED}"
                                     f"\t>>> Error: Data for the Different Material might exist."
                                     f"{Color.RESET}")
                self.HelpList.append(True)
            self.help_check_exit()

            Dtrv = Dtrv1 + Dtrv2
            title = f"{name[0]}_{self.Tilt_Angle}-B12-{int(Angle)}d"
            with open(f"./BandInfo/{title}-HOMO.dat", "w") as Fhomo:
                Fhomo.write(f"{title}-HOMO\n")
                Fhomo.write(f"{Constants.n}\n")
                Fhomo.write(f"{Dcol}\n")
                Fhomo.write(f"{Dtrv}\n")
                Fhomo.write("\n")
                Fhomo.write(f"{T12_HOMO}\n")
                Fhomo.write(f"{T13_HOMO}\n")
                Fhomo.write(f"{T23_HOMO}\n")
                Fhomo.write(f"{T34_HOMO}\n")
                Fhomo.write(f"{T35_HOMO}\n")
                Fhomo.write("\n")
            DatList_temp.append(f"./BandInfo/{title}-HOMO.dat")
            print(f"\t>>> {title}-HOMO.dat: Complete!")

            with open(f"./BandInfo/{title}-LUMO.dat", "w") as Flumo:
                Flumo.write(f"{title}-LUMO\n")
                Flumo.write(f"{Constants.n}\n")
                Flumo.write(f"{Dcol}\n")
                Flumo.write(f"{Dtrv}\n")
                Flumo.write("\n")
                Flumo.write(f"{T12_LUMO}\n")
                Flumo.write(f"{T13_LUMO}\n")
                Flumo.write(f"{T23_LUMO}\n")
                Flumo.write(f"{T34_LUMO}\n")
                Flumo.write(f"{T35_LUMO}\n")
                Flumo.write("\n")
            DatList_temp.append(f"./BandInfo/{title}-LUMO.dat")
            print(f"\t>>> {title}-LUMO.dat: Complete!")

        return DatList_temp

    def mkBandInfo_p3(self):

        p3Data_12, p3Data_23, p3Data_13 = [], [], []
        for p3 in self.p3Files:
            if "-12.txt" in p3:
                p3Data_12 = self.TextFileToData(p3)
            elif "-23.txt" in p3:
                p3Data_23 = self.TextFileToData(p3)
            elif "-31.txt" in p3:
                p3Data_13 = self.TextFileToData(p3)
            else:
                pass

        DatList_temp = []

        for Angle in self.Angles:
            for line in p3Data_12:
                DataList = line.strip().split()
                if float(DataList[1]) == Angle:
                    name = DataList[0].split("_")
                    Dcol = float(DataList[2].strip())
                    Dtrv = float(DataList[3].strip()) * 2
                    T12_HOMO = float(DataList[15])
                    T12_LUMO = float(DataList[13])
                else:
                    pass

            for line in p3Data_23:
                DataList = line.strip().split()
                if float(DataList[1]) == Angle:
                    T23_HOMO = float(DataList[15])
                    T23_LUMO = float(DataList[13])
                else:
                    pass

            for line in p3Data_13:
                DataList = line.strip().split()
                if float(DataList[1]) == Angle:
                    T13_HOMO = float(DataList[15])
                    T13_LUMO = float(DataList[13])
                else:
                    pass

            title = f"{name[0]}_{self.Tilt_Angle}-B3-{int(Angle)}d"
            with open(f"./BandInfo/{title}-HOMO.dat", "w") as Fhomo:
                Fhomo.write(f"{title}-HOMO\n")
                Fhomo.write(f"{Constants.n}\n")
                Fhomo.write(f"{Dcol}\n")
                Fhomo.write(f"{Dtrv}\n")
                Fhomo.write("\n")
                Fhomo.write(f"{T12_HOMO}\n")
                Fhomo.write(f"{T13_HOMO}\n")
                Fhomo.write(f"{T23_HOMO}\n")
                Fhomo.write(f"{T13_HOMO}\n")
                Fhomo.write(f"{T23_HOMO}\n")
                Fhomo.write("\n")
            DatList_temp.append(f"./BandInfo/{title}-HOMO.dat")
            print(f"\t>>> {title}-HOMO.dat: Complete!")

            with open(f"./BandInfo/{title}-LUMO.dat", "w") as Flumo:
                Flumo.write(f"{title}-LUMO\n")
                Flumo.write(f"{Constants.n}\n")
                Flumo.write(f"{Dcol}\n")
                Flumo.write(f"{Dtrv}\n")
                Flumo.write("\n")
                Flumo.write(f"{T12_LUMO}\n")
                Flumo.write(f"{T13_LUMO}\n")
                Flumo.write(f"{T23_LUMO}\n")
                Flumo.write(f"{T13_LUMO}\n")
                Flumo.write(f"{T23_LUMO}\n")
                Flumo.write("\n")
            DatList_temp.append(f"./BandInfo/{title}-LUMO.dat")
            print(f"\t>>> {title}-LUMO.dat: Complete!")

        return DatList_temp

    @staticmethod
    def TextFileToData(path):
        with open(path, "r") as f:
            lines = f.readlines()
        del lines[0:2]
        return lines

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

    def calcEffMass(self, Params, Dat):
        print(f"\tCalculating effective masses for {Color.GREEN}{Dat}{Color.RESET}...")
        Comment = Params["Comment"]
        dev = Params["dev"]
        Dcol = Params["Dcol"]
        Dtrv = Params["Dtrv"]
        TI12 = Params["TI12"]
        TI13 = Params["TI13"]
        TI23 = Params["TI23"]
        TI34 = Params["TI34"]
        TI35 = Params["TI35"]

        # Calculate 2D Energy Dispersion
        Kcol = round(math.pi / Dcol / dev, 10)
        Ktrv = round(math.pi / Dtrv / dev, 10)
        if self.debug:
            print(f"\t>>> Kcol: {Kcol}\n\t>>> Ktrv: {Ktrv}")
        Energy_plus_array, Energy_minus_array = np.zeros((dev + 1, dev + 1)), np.zeros((dev + 1, dev + 1))
        CE = CalculateEnergy()
        for i in range(dev + 1):
            for j in range(dev + 1):
                Energy_plus = CE.Energy_Plus(i, j, TI12, TI13, TI23, TI34, TI35, Kcol, Ktrv)
                Energy_minus = CE.Energy_Minus(i, j, TI12, TI13, TI23, TI34, TI35, Kcol, Ktrv)
                Energy_plus_array[i][j] = Energy_plus
                Energy_minus_array[i][j] = Energy_minus
        if "HOMO" in Comment:
            BandEdge_tmp = np.max(Energy_plus_array)
            BandEdge_index = np.unravel_index(np.argmax(Energy_plus_array), Energy_plus_array.shape)
        elif "LUMO" in Comment:
            BandEdge_tmp = np.min(Energy_minus_array)
            BandEdge_index = np.unravel_index(np.argmin(Energy_minus_array), Energy_minus_array.shape)
        else:
            BandEdge_tmp, BandEdge_index = 0, [0, 0]
            print(f"{Color.RED}\t>>> Error: There is NO information specifying HOMO or LUMO in the file.\n"
                  f"\t>>> The first line (Comment line) should contain HOMO or LUMO.{Color.RESET}")
            self.HelpList.append(True)
        self.help_check_exit()
        if self.debug:
            print(f"Energy of the band edge: {BandEdge_tmp} meV "
                  f"at ({BandEdge_index[0]}, {BandEdge_index[1]}) "
                  f"({BandEdge_index[0] * Kcol}, {BandEdge_index[1] * Ktrv})")

        # バンドエッジ近傍のエネルギーを再計算
        points = 5
        Kcol_for_M = np.zeros(2 * points + 1)
        Ktrv_for_M = np.zeros(2 * points + 1)
        for i in range(2 * points + 1):
            Kcol_for_M[i] = Kcol * (BandEdge_index[0] + i - points)
            Ktrv_for_M[i] = Ktrv * (BandEdge_index[1] + i - points)
        E4M = np.zeros((2 * points + 1, 2 * points + 1))
        for i in range(len(Kcol_for_M)):
            for h in range(len(Ktrv_for_M)):
                Energy_plus = CE.Energy_Plus(Dcol, Dtrv, TI12, TI13, TI23, TI34, TI35, Kcol_for_M[i], Ktrv_for_M[i])
                Energy_minus = CE.Energy_Minus(Dcol, Dtrv, TI12, TI13, TI23, TI34, TI35, Kcol_for_M[i], Ktrv_for_M[i])
                if "HOMO" in Comment:
                    E4M[i][h] = Energy_plus
                elif "LUMO" in Comment:
                    E4M[i][h] = Energy_minus

    @staticmethod
    def calcEnergy(column, transv, TI12, TI13, TI23, TI34, TI35, Kcol, Ktrv):
        B11 = 2 * TI12 * math.cos(Kcol * column)
        B12R = ((TI13 + TI35) * math.cos(Kcol * (column / 2) + Ktrv * (transv / 2))
                + (TI23 + TI34) * math.cos(Kcol * (column / 2) - Ktrv * (transv / 2)))
        B12I = ((TI35 - TI13) * math.sin(Kcol * (column / 2) + Ktrv * (transv / 2))
                + (TI23 - TI34) * math.sin(Kcol * (column / 2) - Ktrv * (transv / 2)))
        B12 = math.sqrt(B12R ** 2 + B12I ** 2)
        Energy_plus = B11 + B12
        Energy_minus = B11 - B12
        return Energy_plus, Energy_minus


class CalculateEnergy:
    def __init__(self):
        D_column, D_transv = sp.symbols("D_column D_transv")
        ti12, ti13, ti23, ti34, ti35 = sp.symbols("TI12 TI13 TI23 TI34 TI35")
        kcol, ktrv = sp.symbols("Kcol Ktrv")

        B11 = 2 * ti12 * sp.cos(kcol * D_column)
        B12R = ((ti13 + ti35) * sp.cos(kcol * (D_column / 2) + ktrv * (D_transv / 2))
                + (ti23 + ti34) * sp.cos(kcol * (D_column / 2) - ktrv * (D_transv / 2)))
        B12I = ((ti35 - ti13) * sp.sin(kcol * (D_column / 2) + ktrv * (D_transv / 2))
                + (ti23 - ti34) * sp.sin(kcol * (D_column / 2) - ktrv * (D_transv / 2)))
        B12 = sp.sqrt(B12R ** 2 + B12I ** 2)

        self.Energy_plus = B11 + B12
        self.Energy_minus = B11 - B12

        diff = sp.diff(self.Energy_plus, kcol)

    def Energy_Plus(self, column, transv, TI12, TI13, TI23, TI34, TI35, Kcol, Ktrv):
        Energy_plus = self.Energy_plus.subs([('D_column', column),
                                             ('D_transv', transv),
                                             ('TI12', TI12),
                                             ('TI13', TI13),
                                             ('TI23', TI23),
                                             ('TI34', TI34),
                                             ('TI35', TI35),
                                             ('Kcol', Kcol),
                                             ('Ktrv', Ktrv)])

        return Energy_plus

    def Energy_Minus(self, column, transv, TI12, TI13, TI23, TI34, TI35, Kcol, Ktrv):
        Energy_minus = self.Energy_minus.subs([('D_column', column),
                                               ('D_transv', transv),
                                               ('TI12', TI12),
                                               ('TI13', TI13),
                                               ('TI23', TI23),
                                               ('TI34', TI34),
                                               ('TI35', TI35),
                                               ('Kcol', Kcol),
                                               ('Ktrv', Ktrv)])

        return Energy_minus


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
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n"
                   "*   To Search the energetically stable aggregation                      *\n"
                   "*            for BrickWork (BW) aggregation                             *\n"
                   "*                                                                       *\n"
                   "*   The programme requires the following arguments.                     *\n"
                   "*     - xyz file of the molecule from which the calculations are made.  *\n"
                   "*                                                                       *\n"
                   "*   The following files are also required.                              *\n"
                   "*     - CalcSetting_BW.txt                                              *\n"
                   "*                                                                       *\n"
                   "*    Option                                                             *\n"
                   "*      - Tilt angle can be added to the molecule                        *\n"
                   "*                          created by Toshiyuki Togashi 2024/11/12      *\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n")
    AbnormalEnd = "\n************* Programme DID NOT terminate successfully. *************\n"
    HelpText = ("\n"
                "The required file may not exist.\n"
                "Please check.")


if __name__ == '__main__':
    main()
