import argparse
import datetime
import functools
import glob
import os
import subprocess
import sys
import time

print = functools.partial(print, flush=True)


class Basefunctions:
    function = {
        1: ["b3lyp/6-31g(d)", "b3lyp_6-31Gd"],
        2: ["b3lyp/6-311+g(d,p)", "b3lyp_6-311+Gdp"],
        3: ["exit", "none"],  # 終了
    }


def main():
    # program start
    # argument parser
    args, before = arg_parser()

    # ReorgEnergy class
    re = ReorgEnergy(args, before)

    # Calculate
    MyJobID = re.calculate()

    # Show log
    re.show_log(MyJobID)

    # program end
    print(f"{Color.GREEN}"
          f"************************* ALL PROCESSES END *************************"
          f"{Color.RESET}\n")
    return


def arg_parser():
    before = time.time()
    parser = argparse.ArgumentParser(description=StandardPhrases.ProgramAbst,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('MaterName', help="Mater Name")
    parser.add_argument('--debug', '-d', '-D',
                        help="Debug mode",
                        action='store_false')
    parser.add_argument('--function', '-f', '-F',
                        help="Change the basis function",
                        action='store_true')
    args = parser.parse_args()

    # View Program Overview
    print(StandardPhrases.ProgramAbst)

    if args.debug:
        print(f"{Color.RED}*************** Caution!!! Debug Started!!! ***************{Color.RESET}")
        print(f"In debug mode, the following files are not deleted\n"
              f"\t - sh files\n"
              f"\t - chk files\n"
              f"\t - log files\n"
              f"\t - gjf files\n")
    return args, before


class ReorgEnergy:
    def __init__(self, args, before):
        self._ReorgEnergy = "ReorgEnergy"
        self.args = args
        self.before = before
        self.MaterName = args.MaterName[:-4]
        self.debug = args.debug
        self.function = args.function
        self.messages, self.HelpList = [], []
        self.Info_for_Debug = []
        self.StartTime = datetime.datetime.now()

        # Retrieve the operator name
        Operator = input(f"\n{Color.GREEN}Retrieve the operator name.{Color.RESET}\n"
                         "\t>>> The name entered here will be used to identify the operator.\n"
                         f"\tPlease enter the operator.\n"
                         f"\t{Color.GREEN}>>> {Color.RESET}")
        if Operator == "":
            Operator = "ONE"
        else:
            pass
            print("")
        self.Operator = Operator

        # Check the required files
        self.File_Set_Check()

        # Set the basis function
        print("\nSetting the basis function...")
        if self.function:
            while True:
                print(f"\n"
                      f"{Color.BLUE}Selected basis function{Color.RESET}")
                for key, value in Basefunctions.function.items():
                    print(f"\t{key}: {value[0]}")
                function_set = input("\t>>> ")

                if function_set.isdecimal() and 1 <= int(function_set) <= len(Basefunctions.function):
                    basis_function = Basefunctions.function.get(int(function_set))
                    if basis_function[0] == "exit":
                        print(f"\t{Color.RED}>>> Exit the program.{Color.RESET}")
                        exit()
                    else:
                        pass
                    break
                else:
                    print(f"{Color.RED}Error: Invalid input.{Color.RESET}")
        else:
            basis_function = ["b3lyp/6-31g(d)", "b3lyp_6-31Gd"]
        print(f"\t>>> The basis function is set to {basis_function[1]}.")
        self.basis_function = basis_function

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

        print(f"\t{self.MaterName}.gjf: ", end="")
        if os.path.exists(f"{self.MaterName}.gjf"):
            print(f"{Color.GREEN}Found{Color.RESET}")
        else:
            print(f"{Color.RED}Not Found{Color.RESET}")
            self.HelpList.append(True)
            self.messages.append(f"{Color.RED}\t>>> Some required files are not available.\n"
                                 f"\t>>> Please upload the required files and try again.{Color.RESET}")
        self.help_check_exit()
        print(f"\n\t{Color.GREEN}>>> All required files are available.{Color.RESET}")
        return None

    def job_submission(self, qsubList, which, dirpath):
        """
        Submit the job
        :param qsubList:
        :param which:
        :param dirpath:
        :return:
        """
        MyJobIDList = []
        print("\n**********\nJobs are submitting...")
        if len(qsubList) == 0:
            self.messages.append(f"\t>>> Job was not submitted.\n"
                                 f"\t>>> {Color.GREEN}Calculations with the conditions might be finished.{Color.RESET}")
            self.help_check_exit()
        else:
            for qsub in qsubList:
                qsub = qsub.split()
                subprocess.run(qsub, cwd=dirpath)

            Wait_minutes = 1
            MyJobIDList = self.My_JobIDList(qsubList)
            MyJobIDList.sort()
            RunningJobIDList = self.Running_JobIDList()

            Flag, wait_job_count = self.check_jobs(RunningJobIDList, MyJobIDList)
            start_time = datetime.datetime.now()

            if Flag:
                formated_ST = start_time.strftime("%m/%d %H:%M:%S")
                print(f"{Color.GREEN}\t>>> '{int(len(MyJobIDList))}' calculations for ReorgEnergy was submitted!!"
                      f" {Color.RESET}at {formated_ST}")
                print(f"\n"
                      f"{Color.GREEN}Wait until jobID {MyJobIDList[-1]}!!{Color.RESET}\n"
                      f"start ID: {RunningJobIDList[0]}\n"
                      f"\n")
            else:
                pass

            while Flag:
                Flag, wait_job_count = self.check_jobs(RunningJobIDList, MyJobIDList)
                formated_NOW, elapsed_time = self.getElapsedTime(start_time)
                Now_Time = datetime.datetime.now()
                formated_end_time = (
                    (Now_Time + datetime.timedelta(minutes=(wait_job_count * Wait_minutes))).strftime("%m/%d %H:%M:%S"))
                sys.stdout.write(
                    "\033[1F\033[G%s" %
                    f"\t{formated_NOW} ({elapsed_time} min. passed): Job {RunningJobIDList[0]} is in progress.       \n"
                    f"\tNext Check >>> {Wait_minutes * wait_job_count} minute later! forecast: {formated_end_time}    ")
                sys.stdout.flush()
                time.sleep(Wait_minutes * wait_job_count * 60)
                RunningJobIDList = self.Running_JobIDList()
                Flag, wait_job_count = self.check_jobs(RunningJobIDList, MyJobIDList)
            else:
                print(f"{Color.GREEN}\n\n"
                      f"Calculation cycles for {which} until JobID {MyJobIDList[-1]} were finished.{Color.RESET}")
        return MyJobIDList[-1]

    @staticmethod
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

    @staticmethod
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

    @staticmethod
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

    @staticmethod
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

    def calculate(self):
        print(f"{Color.GREEN}Calculating...{Color.RESET}")

        if self.debug:
            Debug_arg = "--debug"
        else:
            Debug_arg = ""

        with open(f"G-Reorg-{self.Operator}-{self.MaterName}-{self.basis_function[1]}.sh", "w") as file:
            file.write(f"#!/bin/sh\n")
            file.write(f"\n")
            file.write(f"#$ -S /bin/sh\n")
            file.write(f"#$ -cwd\n")
            file.write(f"#$ -V\n")
            file.write(f"#$ -pe gau 12\n")
            file.write(f"#$ -q all.q\n")
            file.write(f"\n")
            file.write(f"module load gaussian/g16\n")
            file.write(f"export GAUSS_SCRDIR=/scr/$JOB_ID\n")
            file.write(f"mkdir /scr/$JOB_ID\n")
            file.write(f"\n")
            file.write(f"BG_ReorgEnergy {self.MaterName} {self.basis_function[1]} {Debug_arg} "
                       f"| tee Reorg-{self.Operator}-{self.MaterName}-{self.basis_function[1]}.log\n")
            file.write(f"rm -rf /scr/$JOB_ID\n")
            file.write(f"\n")

        # Execute the shell script
        qsubList = [f"qsub G-Reorg-{self.Operator}-{self.MaterName}-{self.basis_function[1]}.sh"]
        JobID = self.job_submission(qsubList, "Reorg_BG", os.getcwd())

        print(f"\n{Color.GREEN}All calculations are finished!!{Color.RESET}")
        return JobID

    def show_log(self, JobID):
        print(f"\nChecking whether the calculations are finished successfully...")
        if not os.path.getsize(f"G-Reorg-{self.Operator}-{self.MaterName}-{self.basis_function[1]}.sh.e{JobID}") == 0:
            print(f"\t>>> {Color.RED}Calculations are not finished successfully.{Color.RESET}\n"
                  f"\t>>> Please check the calculation conditions.\n")

            print("********************************************************************************")
            with open(f"G-Reorg-{self.Operator}-{self.MaterName}-{self.basis_function[1]}.sh.e{JobID}", "r") as file:
                print(file.read())
            print("********************************************************************************")
            self.HelpList.append(True)
            self.messages.append(f"\n\t>>> The calculation was not successful.\n"
                                 f"\t>>> Please check the calculation conditions.")
            self.help_check_exit()
        else:
            print(f"\t>>> {Color.GREEN}Calculations are finished successfully!!{Color.RESET}")

        print(f"{Color.GREEN}Log file is displayed.{Color.RESET}")
        with open(f"{self.MaterName}_ReorgEnergy_{self.basis_function[1]}.txt", "r") as file:
            print(file.read())

        if not self.debug:
            self.rmWildCards("*.chk")
            self.rmWildCards("*.sh*")
        return None

    @staticmethod
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


class StandardPhrases:
    def __init__(self):
        self._StandardPhrases = "StandardPhrases"

    ProgramAbst = ("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n"
                   "*   Programme for calculating relocation energy.            *\n"
                   "*                                                           *\n"
                   "*   This programme requires the following files.            *\n"
                   "*     - A gjf file of the molecule you want to calculate.   *\n"
                   "*                                                           *\n"
                   "*   This programme requires the following arguments.        *\n"
                   "*     - Molecule name.                                      *\n"
                   "*                created by Toshiyuki Togashi 2024/10/10    *\n"
                   "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n")

    AbnormalEnd = "\n************* Programme DID NOT terminate successfully. *************\n"


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
