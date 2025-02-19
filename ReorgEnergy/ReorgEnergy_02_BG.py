import argparse
import datetime
import functools
import os
import platform
import subprocess
import sys
import time

print = functools.partial(print, flush=True)


class Basefunctions:
    function = {
        "b3lyp_6-31Gd": "b3lyp/6-31g(d)",
        "b3lyp_6-311+Gdp": "b3lyp/6-311+g(d,p)",
        "b3lyp_cc-pVTZ": "b3lyp/cc-pvtz",
    }


def main():
    before = time.time()
    print("**************************************************")
    print("*              ReorgEnergy for NITT              *")
    print("*        Created by T. Togashi (24/12/10)        *")
    print("**************************************************")
    args = arg_parser()
    print(f"Material Name: {args.MaterName}")
    print(f"Function: {args.Function}")
    print(f"Start Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    re = ReorgEnergy(args)

    re.calc_from_gjf("EG", "0")  # 0荷、0価の構造
    re.calc_from_log("EG", "+1", f"{args.MaterName}_0_EG_{args.Function}.log")  # +1荷、+1価の構造
    re.calc_from_log("EG", "-1", f"{args.MaterName}_0_EG_{args.Function}.log")  # -1荷、-1価の構造
    re.calc_from_log("SP", "+0", f"{args.MaterName}_+1_EG_{args.Function}.log")  # 0荷、+1価の構造
    re.calc_from_log("SP", "-0", f"{args.MaterName}_-1_EG_{args.Function}.log")  # 0荷、-1価の構造
    re.calc_from_log("SP", "+1", f"{args.MaterName}_0_EG_{args.Function}.log")  # +1荷、0価の構造
    re.calc_from_log("SP", "-1", f"{args.MaterName}_0_EG_{args.Function}.log")  # -1荷、0価の構造

    result_meV = re.calc_reorg_energy()

    re.save_result(result_meV)

    after = time.time()
    elapsed_time = after - before
    formatted_time = time.strftime("%H h %M m %S s", time.gmtime(elapsed_time))
    print(f"")
    print(f"Start Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Elapsed Time: {formatted_time}")


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('MaterName')
    parser.add_argument('Function')
    parser.add_argument('--debug', '-d',
                        action='store_true')
    parser.add_argument('-g', '--g09')
    args = parser.parse_args()

    return args


class ReorgEnergy:
    def __init__(self, args):
        self.args = args
        self.MaterName = args.MaterName
        self.debug = args.debug
        self.Function_Name = args.Function
        self.function = Basefunctions.function[args.Function]
        self.bond_info = self.get_band_info()
        if args.g09:
            self.gaussian_command = 'g09'
        else:
            self.gaussian_command = 'g16'
        self.EnergyList = []
        self.Freq_0_EG, self.Freq_plus1_EG, self.Freq_minus1_EG = [], [], []

    def get_band_info(self):
        print("\n結合情報を取得しています...")
        bond_info = []
        with open(f"{self.MaterName}.gjf", 'r') as file:
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

        print("結合情報を取得しました。\n")

        return "\n".join(bond_info)

    @staticmethod
    def format_coordinate(coord):
        """
        座標をフォーマットする関数
        :param coord:
        :return:
        """
        x, y, z = map(float, coord)
        x = f'{x:.10f}'.rjust(15, ' ')
        y = f'{y:.10f}'.rjust(15, ' ')
        z = f'{z:.10f}'.rjust(15, ' ')
        return f'{x} {y} {z}'

    def calc_from_gjf(self, EG_or_SP, Charge):
        if os.path.exists(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.log"):
            print(f"{Charge}荷、{EG_or_SP}の計算は既に終了しています。")
        else:
            print(f"{Charge}荷、{EG_or_SP}のgjfファイルを作成しています...")
            with open(f"{self.MaterName}.gjf", "r") as file:
                data = file.read()
            List = data.splitlines()[6:]
            data = []
            for element in List:
                try:
                    coordinate = [float(element.split()[1]), float(element.split()[2]), float(element.split()[3])]
                    formatted_coordinate = self.format_coordinate(coordinate)
                    formatted_line = f" {element.split()[0]}                 {formatted_coordinate}"
                    data.append(formatted_line)
                except IndexError:
                    break

            gjf = self.write_gjf_file(EG_or_SP, Charge, data)
            print(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.gjfを作成しました。")

            print(f"{Charge}荷、{EG_or_SP}を計算しています...")
            self.run_gaussian(gjf)
            print(f"{Charge}荷、{EG_or_SP}の計算が終了しました。")

        self.Check_normal_termination(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.log", True)
        if EG_or_SP == "EG":
            self.Check_negative_frequency(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.log", Charge)
        print("")
        return None

    def calc_from_log(self, EG_or_SP, Charge, Base_Log_File):
        if os.path.exists(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.log"):
            print(f"{Charge}荷、{EG_or_SP}の計算は既に終了しています。")
        else:
            print(f"{Charge}荷、{EG_or_SP}のgjfファイルを作成しています...")
            data = self.molecular_log_open(Base_Log_File)
            formatted_data = self.format_atomic_data(data)
            gjf = self.write_gjf_file(EG_or_SP, Charge, formatted_data)
            print(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.gjfを作成しました。")

            print(f"{Charge}荷、{EG_or_SP}を計算しています...")
            self.run_gaussian(gjf)
            print(f"{Charge}荷、{EG_or_SP}の計算が終了しました。")

        self.Check_normal_termination(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.log", True)
        if EG_or_SP == "EG":
            self.Check_negative_frequency(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.log", Charge)
        print("")
        return None

    def write_gjf_file(self, EG_or_SP, Charge, formatted_data):
        # 仮のヘッダーを作成
        header_data = StandardPhrases.Header_Template.splitlines()
        header_data[0] = header_data[0] + "\n"
        header_data[1] = header_data[1] + "\n"
        header_data[2] = f"%chk={self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.chk" + "\n"
        # EGかSPかで分岐
        if EG_or_SP == "EG":
            header_data[3] = f"# opt freq=noraman {self.function} geom=connectivity" + "\n"
        elif EG_or_SP == "SP":
            header_data[3] = f"# {self.function} geom=connectivity" + "\n"
        else:
            print("Error: EG or SP is not specified.")
            sys.exit()
        header_data[4] = header_data[4] + "\n"
        header_data[5] = header_data[5] + "\n"
        header_data[6] = header_data[6] + "\n"
        # Chargeによって分岐
        if Charge == "0" or Charge == "+0" or Charge == "-0":
            header_data[7] = "0 1\n"
        elif Charge == "+1":
            header_data[7] = "1 2\n"
        elif Charge == "-1":
            header_data[7] = "-1 2\n"
        # GJFファイルへの書き込み
        with open(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.gjf", "w") as file:
            for formatted_line in header_data:
                file.write(formatted_line)
            for formatted_line in formatted_data:
                file.write(formatted_line + "\n")
            file.write("\n")
            file.write(self.bond_info)
            file.write("\n")
            file.write("\n")
        if self.debug:
            print(f"\n************* {self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.gjf *************")
            with open(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.gjf", "r") as file:
                data = file.read()
            print(data)
            print("**************************************************************************\n")
        return f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.gjf"

    def run_gaussian(self, gjf):
        command_list = [self.gaussian_command, gjf]
        res = subprocess.run(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

        # check error
        if res.returncode:
            if platform.system() == 'Windows':
                print(f"Failed to execute {' '.join(command_list)}")
            else:
                print(res.stderr.strip())

        return None

    @staticmethod
    def Check_normal_termination(Log_File, Flag):
        if Flag:
            print("計算が正常に終了したかを確認しています...")
        with open(Log_File, "r") as file:
            data = file.read()
        if "Normal termination" in data:
            if Flag:
                print("計算は正常に終了しました。")
        else:
            if Flag:
                print(f"{Color.RED}計算が正常に終了しませんでした。{Color.RESET}")
                print(f"{Color.RED}エラーが発生している可能性があります。{Color.RESET}")
                print(f"{Color.RED}ログファイルを確認してください。{Color.RESET}")
                print(f"{Color.RED}プログラムを終了します。{Color.RESET}")
            exit()
        return None

    def Check_negative_frequency(self, Log_File, Charge):
        print("振動数が負の値を持っていないかを確認しています...")
        count, messages, HelpList = 0, [], []
        with open(Log_File, "r") as file:
            data = file.read()
        Temp = data.split("Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering")[1]
        for line in Temp.strip().splitlines():
            if "Frequencies" in line:
                elements = line.split()
                for element in elements:
                    try:
                        element = float(element)
                        count = count + 1
                        if Charge == "0":
                            self.Freq_0_EG.append([count, element])
                        elif Charge == "+1":
                            self.Freq_plus1_EG.append([count, element])
                        elif Charge == "-1":
                            self.Freq_minus1_EG.append([count, element])
                        else:
                            pass
                        if self.debug:
                            print(f"\t>>> 周波数 {count}: {element}")
                        if element < 0:
                            messages.append(
                                f"\t{Color.RED}>>> Negative frequency found {count}: {element}{Color.RESET}")
                            HelpList.append(True)
                    except ValueError:
                        pass
        if True in HelpList:
            print(f"{Color.RED}EG, {Charge}価での計算結果にて負の振動数が見つかりました。{Color.RESET}")
            sys.stderr.write(f"{Color.RED}EG, {Charge}価での計算結果にて負の振動数が見つかりました。{Color.RESET}\n")
            for message in messages:
                print(message)
                sys.stderr.write(f"{message}\n")
            print(f"{Color.RED}ログファイルを確認してください。{Color.RESET}")
            print(f"{Color.RED}プログラムを終了します。{Color.RESET}")
        else:
            print("負の振動数は見つかりませんでした。")
        return None

    def molecular_log_open(self, Log_File):
        with open(Log_File, "r") as file:
            data = file.read()
        self.Check_normal_termination(Log_File, False)
        return data

    def format_atomic_data(self, data):
        """
        データを加工する関数
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
            formatted_coordinate = self.format_coordinate(coordinate)

            # 元素記号を取得
            symbol = None
            try:
                symbol = PeriodicTable.atomic_symbols.get(atomic_number)[0]
            except KeyError:
                print(f"予期しない原子を検出しました。")
                print(f"原子番号: {atomic_number}")
                print(f"ファイルを確認するか、プログラムに原子を追加してください。")
                exit()

            # フォーマットされた行を作成
            formatted_line = f" {symbol:<2}  {formatted_coordinate}"
            formatted_data.append(formatted_line)

        return formatted_data

    def calc_reorg_energy(self):
        HelpList = []
        print("logファイルを読み込んでいます...")
        self.get_energy("0", "EG")
        self.get_energy("+1", "EG")
        self.get_energy("-1", "EG")
        self.get_energy("+0", "SP")
        self.get_energy("-0", "SP")
        self.get_energy("+1", "SP")
        self.get_energy("-1", "SP")

        energy_0_EG = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '0' and item['EG_or_SP'] == 'EG'
            )
        )
        energy_plus0_SP = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '+0' and item['EG_or_SP'] == 'SP'
            )
        )
        energy_plus1_EG = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '+1' and item['EG_or_SP'] == 'EG'
            )
        )
        energy_plus1_SP = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '+1' and item['EG_or_SP'] == 'SP'
            )
        )
        energy_minus0_SP = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '-0' and item['EG_or_SP'] == 'SP'
            )
        )
        energy_minus1_EG = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '-1' and item['EG_or_SP'] == 'EG'
            )
        )
        energy_minus1_SP = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '-1' and item['EG_or_SP'] == 'SP'
            )
        )

        # カラムのヘッダー
        headers = ["Charge", "EG_or_SP", "Energy"]
        # 各カラムの幅を求める
        column_widths = {header: max(len(str(row[header])) for row in self.EnergyList + [dict(zip(headers, headers))])
                         for header in headers}
        # ヘッダー行の表示
        header_row = " | ".join(f"{header:<{column_widths[header]}}" for header in headers)
        print("-" * len(header_row))
        print(header_row)
        print("-" * len(header_row))
        # データ行の表示
        for row in self.EnergyList:
            print(" | ".join(f"{str(row[header]):<{column_widths[header]}}" for header in headers))
        # 下線の表示
        print("-" * len(header_row))

        # 式に従って計算
        anion_lambda1 = energy_plus0_SP - energy_0_EG
        anion_lambda2 = energy_plus1_SP - energy_plus1_EG
        result_Hartree_anion = (anion_lambda1 + anion_lambda2)
        result_meV_anion = self.hartree_to_meV(result_Hartree_anion)
        if anion_lambda1 < 0 or anion_lambda2 < 0:
            print(f"{Color.RED}minus value of λ1 or/and λ2 is found in calculations for anion.{Color.RESET}")
            print(f"anion λ1: {anion_lambda1}")
            print(f"anion λ2: {anion_lambda2}")
            print(f"ReorgEnergy [Hartree]: {result_Hartree_anion}")
            print(f"ReorgEnergy [meV]: {result_meV_anion}")

            sys.stderr.write(f"{Color.RED}minus value of λ1 or/and λ2 is found "
                             f"in calculations for anion.{Color.RESET}\n")
            sys.stderr.write(f"anion λ1: {anion_lambda1}\n")
            sys.stderr.write(f"anion λ2: {anion_lambda2}\n")
            sys.stderr.write(f"ReorgEnergy [Hartree]: {result_Hartree_anion}\n")
            sys.stderr.write(f"ReorgEnergy [meV]: {result_meV_anion}\n")
            HelpList.append(True)

        cation_lambda1 = energy_minus0_SP - energy_0_EG
        cation_lambda2 = energy_minus1_SP - energy_minus1_EG
        result_Hartree_cation = (cation_lambda1 + cation_lambda2)
        result_meV_cation = self.hartree_to_meV(result_Hartree_cation)
        if cation_lambda1 < 0 or cation_lambda2 < 0:
            print(f"{Color.RED}minus value of λ1 or/and λ2 is found in calculations for cation.{Color.RESET}")
            print(f"cation λ1: {cation_lambda1}")
            print(f"cation λ2: {cation_lambda2}")
            print(f"ReorgEnergy [Hartree]: {result_Hartree_cation}")
            print(f"ReorgEnergy [meV]: {result_meV_cation}")

            sys.stderr.write(f"{Color.RED}minus value of λ1 or/and λ2 is found "
                             f"in calculations for cation.{Color.RESET}\n")
            sys.stderr.write(f"cation λ1: {cation_lambda1}\n")
            sys.stderr.write(f"cation λ2: {cation_lambda2}\n")
            sys.stderr.write(f"ReorgEnergy [Hartree]: {result_Hartree_cation}\n")
            sys.stderr.write(f"ReorgEnergy [meV]: {result_meV_cation}\n")
            HelpList.append(True)

        if True in HelpList:
            sys.exit(1)

        print("logファイルを読み込みました。")

        # 結果の表示
        print("\n-- Reorganization Energy --")

        print("Anion:")
        print(f"λ1 (+0_SP - 0_EG): {round(energy_plus0_SP - energy_0_EG, 8)}[Hartree]")
        print(f"λ2 (+1_SP - +1_EG): {round(energy_plus1_SP - energy_plus1_EG, 8)}[Hartree]")
        print("再配置エネルギー[Hartree]:", round(result_Hartree_anion, 8), "[Hartree]")
        print("再配置エネルギー[meV]:", round(result_meV_anion, 8), "[meV]")

        print("\nCation:")
        print(f"λ1 (-0_SP - 0_EG): {round(energy_minus0_SP - energy_0_EG, 8)}[Hartree]")
        print(f"λ2 (-1_SP - -1_EG): {round(energy_minus1_SP - energy_minus1_EG, 8)}[Hartree]")
        print("再配置エネルギー[Hartree]:", round(result_Hartree_cation, 8), "[Hartree]")
        print("再配置エネルギー[meV]:", round(result_meV_cation, 8), "[meV]")

        return [result_meV_anion, result_meV_cation]

    def get_energy(self, Charge, EG_or_SP):
        with open(f"{self.MaterName}_{Charge}_{EG_or_SP}_{self.Function_Name}.log", "r") as file:
            data = file.read()
        Temp = data.split("SCF Done")[-1]
        Temp = Temp.split("=")[1]
        Temp = Temp.split("A.U. after")[0].strip()

        DataList = {"Charge": Charge,
                    "EG_or_SP": EG_or_SP,
                    "Energy": Temp}
        self.EnergyList.append(DataList)
        return None

    @staticmethod
    def hartree_to_meV(Hartree):
        """
        HartreeをmeVに変換する関数
        :param Hartree:
        """
        meV = Hartree * 1000 * Constants.Hartree_to_eV
        return meV

    def save_result(self, result_meV):
        print("\n結果を保存しています...")

        energy_0_EG = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '0' and item['EG_or_SP'] == 'EG'
            )
        )
        energy_plus0_SP = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '+0' and item['EG_or_SP'] == 'SP'
            )
        )
        energy_plus1_EG = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '+1' and item['EG_or_SP'] == 'EG'
            )
        )
        energy_plus1_SP = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '+1' and item['EG_or_SP'] == 'SP'
            )
        )
        energy_minus0_SP = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '-0' and item['EG_or_SP'] == 'SP'
            )
        )
        energy_minus1_EG = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '-1' and item['EG_or_SP'] == 'EG'
            )
        )
        energy_minus1_SP = float(
            next(
                item['Energy'] for item in self.EnergyList
                if item['Charge'] == '-1' and item['EG_or_SP'] == 'SP'
            )
        )

        min_freq_0_EG = min([freq[1] for freq in self.Freq_0_EG])
        min_freq_plus1_EG = min([freq[1] for freq in self.Freq_plus1_EG])
        min_freq_minus1_EG = min([freq[1] for freq in self.Freq_minus1_EG])

        with open(f"{self.MaterName}_ReorgEnergy_{self.Function_Name}.txt", "w") as file:
            file.write(f"Execution date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            file.write(f"Material Name: {self.MaterName}\n")
            file.write(f"Function: {self.function}\n")
            file.write(f"\n")
            file.write(f"Reorganization Energy\n")
            file.write(f"Anion: {round(result_meV[0], 8)} [meV]\n")
            file.write(f"Cation: {round(result_meV[1], 8)} [meV]\n")

            file.write("\n" + "=" * 40 + "\n")
            file.write(f"{'Charge':<10}{'EG/SP':<10}{'Energy [Hartree]':<20}\n")
            file.write("=" * 40 + "\n")
            for item in self.EnergyList:
                file.write(f"{item['Charge']:<10}{item['EG_or_SP']:<10}{item['Energy']:<20}\n")
            file.write("=" * 40 + "\n")

            file.write(f"-Anion-\n")
            file.write(f"λ1 (-0_SP - 0_EG): {round(energy_minus0_SP - energy_0_EG, 8)} [Hartree]\n")
            file.write(f"λ2 (-1_SP - -1_EG): {round(energy_minus1_SP - energy_minus1_EG, 8)} [Hartree]\n")
            file.write(f"\n")
            file.write(f"-Cation-\n")
            file.write(f"λ1 (+0_SP - 0_EG): {round(energy_plus0_SP - energy_0_EG, 8)} [Hartree]\n")
            file.write(f"λ2 (+1_SP - +1_EG): {round(energy_plus1_SP - energy_plus1_EG, 8)} [Hartree]\n")
            file.write(f"\n")
            file.write(f"Minimum Frequency of 0_EG: {min_freq_0_EG}\n")
            file.write(f"Minimum Frequency of +1_EG: {min_freq_plus1_EG}\n")
            file.write(f"Minimum Frequency of -1_EG: {min_freq_minus1_EG}\n")

        print("結果を保存しました。")
        return None


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


class StandardPhrases:
    def __init__(self):
        self._StandardPhrases = "StandardPhrases"

    Header_Template = ("%nprocshared=12\n"
                       "%mem=32GB\n"
                       "%chk=Template.chk\n"
                       "# opt freq=noraman b3lyp/6-31g(d) geom=connectivity\n"
                       "\n"
                       "Title Card Required\n"
                       "\n"
                       "0 1\n")


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


class Constants:
    # ハートリーからeVに変換するための定数
    Hartree_to_eV = 27.21138602
    # 気体定数[R] (J/(mol・K))
    R = 8.314462618
    # 真空中の光速[m/s]
    C = 299792458
    # 電気素量[C]
    E = 1.602176634e-19


if __name__ == "__main__":
    main()
