#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import functools
import os
import time

print = functools.partial(print, flush=True)


def main():
    args, before = get_args()
    EM = EffectiveMass(args)

    # Make band information files
    print("\n**********\nMaking band information files...")


def get_args():
    before = time.time()
    help_text = StandardPhrases.ProgramAbst

    parser = argparse.ArgumentParser(description=help_text, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('MaterName',
                        help='Material name')
    parser.add_argument('--debug', '-d', '-D',
                        help="Start in debug mode",
                        action='store_true')
    parser.add_argument('--gif', '-g', '-G',
                        help="Create GJF File from All File",
                        action='store_true')
    args = parser.parse_args()

    if args.debug:
        print(f"{Color.RED}*************** Caution!!! Debug Started!!! ***************{Color.RESET}\n")

    return args, before


class EffectiveMass:
    def __init__(self, args):
        self.args = args
        self.MaterName = args.MaterName
        self.debug = args.debug
        self.gif = args.gif
        self.messages = []
        self.HelpList = []
        self.p1Files, self.AllFiles, self.Other = self.File_Set_Check()
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

    def File_Set_Check(self):
        print(f"{Color.GREEN}Checking the required files...{Color.RESET}")
        p1Files = self.getMinTIFileName()
        AllFiles = self.getAllFileNames()
        self.help_check_exit()
        print(f"{Color.GREEN}All required files are available.{Color.RESET}")

        print("\nTI files:")
        for p1File in p1Files:
            print(f"\t{p1File}")

        print("\nAll files:")
        for AllFile in AllFiles:
            print(f"\t{AllFile}")

        Other = []
        for AllFile in AllFiles:
            Other_temp = self.getOther(AllFile)
            Other = Other + Other_temp
        Other = list(set(Other))
        Other.sort()
        if self.debug:
            print(f"\nOther: {Other}")

        return p1Files, AllFiles, Other

    def getMinTIFileName(self):
        Dirs = os.listdir("./")
        p1Files = []
        for Dir in Dirs:
            if self.MaterName in Dir and "_results" in Dir:
                Files_temp = os.listdir(Dir)
                for File in Files_temp:
                    if "3molp1" in File and "min-TIs-" in File:
                        p1Files.append("./" + Dir + "/" + File)

        if len(p1Files) == 0:
            print(f"{Color.RED}\t>>> Error: No files found for 3molp1.{Color.RESET}")
            self.HelpList.append(True)
        elif len(p1Files) != 3:
            print(f"{Color.RED}\t>>> Error: 3molp1 files are not complete.{Color.RESET}")
            self.HelpList.append(True)
        else:
            print(f"\t>>> 3molp1 files: {Color.GREEN}All Available.{Color.RESET}")

        return p1Files

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
    def getOther(AllFile):
        with open(AllFile, "r") as f:
            lines = f.readlines()
        Other_temp = []
        for line in lines:
            Contents = line.strip().split()
            try:
                Other_temp.append(float(Contents[1]))
            except (ValueError, IndexError):
                pass
        return Other_temp


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


if __name__ == "__main__":
    main()
