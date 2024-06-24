import re
import os
import json
from math import pi, exp
import math
from numpy import array
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# FILE_DIR = r'D:\shujufenxi\Y6\y6\y6ea'
# V_Absolute_Value_File = r'D:\wq-data-statistic\Y6\ouheshuju\N35.990687'


# 数字_A数字.log --> ENUT 文件
# 数字oni_A数字.log --> EANION 文件
# 1_A1.log --> 1_oni_A1.log
# 每一个 1_A1.log 都有一个对应的 1oni_A1.log


class FindSCF:
    """
    封装成类，方便调用
    """

    def __init__(self, file_dir, v_absolute_value_file):
        self.file_dir = file_dir
        self.v_absolute_value_file = v_absolute_value_file
        self.ALL_LOG_FILE_LIST = [f for f in os.listdir(file_dir) if
                                  os.path.isfile(os.path.join(file_dir, f)) and f.endswith('.log')]

        self.file_ea_count = 0
        self.final_value_list_ea = []
        self.final_E_value_list = []
        self.pair_list = []
        self.result_dict_homo = {}
        self.result_dict_lumo = {}
        self.csv_filename = None

    def get_ea_file_list(self):
        enut_file_list = []
        enut_pattern = re.compile(r'\d+\.log')

        for log_file_name in self.ALL_LOG_FILE_LIST:
            if enut_pattern.match(log_file_name):
                enut_file_list.append(log_file_name)
        print('所有的enut文件B：', enut_file_list, '\n', 'enut文件B数量：', len(enut_file_list))
        # 查找enut file 对应的 eanion file
        # 最终生成字典 {enut: enut_file, eanion: eanion_file}
        self.file_ea_count = len(enut_file_list)

        final_file_list = []
        for enut_file_name in enut_file_list:
            eanion_file_name = enut_file_name.replace('.log', '_oni.log')
            eacat_file_name = enut_file_name.replace('.log', '_cat.log')
            if eacat_file_name not in self.ALL_LOG_FILE_LIST:
                print('没有找到对应的eanion file', eacat_file_name)
            if eanion_file_name in self.ALL_LOG_FILE_LIST:
                final_file_list.append({
                    'serial_number': int(eanion_file_name.split('_oni')[0]),
                    'enut': enut_file_name,
                    'eanion': eanion_file_name,
                    'eacat': eacat_file_name
                })
            else:
                print('没有找到对应的eanion file:', eanion_file_name)
        self.final_file_list_ea = final_file_list

    def get_scf_value(self):
        for file_dict in self.final_file_list_ea:
            enut_file_path = os.path.join(self.file_dir, file_dict['enut'])
            eanion_file_path = os.path.join(self.file_dir, file_dict['eanion'])
            eacat_file_path = os.path.join(self.file_dir, file_dict['eacat'])
            value_dict = {
                'serial_number': file_dict['serial_number'],  # 序号，好像是婉清需要的
                'enut_file': file_dict['enut'],
                'enut_value': None,
                'eanion_file': file_dict['eanion'],
                'eanion_value': None,
                'eacat_file': file_dict['eacat'],
                'eacat_value': None,
                'e_total_value': None,
                'i_total_value': None,
            }
            # 先读取enut file
            with open(enut_file_path, 'r') as enut_file:
                enut_lines = enut_file.readlines()
                for line in enut_lines:
                    if 'SCF Done:' in line:
                        line_str_list = line.split(' ')
                        value = float(line_str_list[7])
                        value_dict['enut_value'] = value
                        break
                    if line == enut_lines[-1]:
                        print('not find enut file')
            # 再读取eanion file

            with open(eanion_file_path, 'r') as eanion_file:
                eanion_lines = eanion_file.readlines()
                # i=i+1
                for line in eanion_lines:
                    if 'SCF Done:' in line:
                        line_str_list = line.split(' ')
                        value = float(line_str_list[7])
                        print(line)
                        print(value_dict['eanion_file'])
                        value_dict['eanion_value'] = value
                        break
                    if line == enut_lines[-1]:
                        print('not find ea file')
            # 读取cations file
            with open(eacat_file_path, 'r') as eacat_file:
                eacat_lines = eacat_file.readlines()
                for line in eacat_lines:
                    if 'SCF Done:' in line:
                        line_str_list = line.split(' ')
                        value = float(line_str_list[7])
                        value_dict['eacat_value'] = value
                        break
                    if line == enut_lines[-1]:
                        print('not find cat file')
            # 计算e_total
            print(value_dict['eanion_value'])
            value_dict['e_total_value'] = abs(value_dict['enut_value'] - value_dict['eanion_value']) * 27.2114
            # 计算i_total
            value_dict['i_total_value'] = abs(value_dict['eacat_value'] - value_dict['enut_value']) * 27.2114
            self.final_value_list_ea.append(value_dict)
            # print('fianl_value_list', self.final_value_list_ea)


    def calculate_delta_E2(self):
        for i, value_list in enumerate(self.pair_list, start=1):
            data = {'serial_number': i,
                    'pair_list': value_list,
                    'delta_ea': None,
                    'delta_ip': None,
                    'A': None,
                    'B': None
                    }
            print(type(self.final_value_list_ea))
            print('self.final_value_list_ea::', self.final_value_list_ea)
            for value_dict in self.final_value_list_ea:
                if value_dict['serial_number'] == value_list[0]:
                    data['A'] = value_dict
                if value_dict['serial_number'] == value_list[1]:
                    data['B'] = value_dict
            data['delta_ea'] = abs(data['A']['e_total_value'] - data['B']['e_total_value'])
            data['delta_ip'] = abs(data['A']['i_total_value'] - data['B']['i_total_value'])
            self.final_E_value_list.append(data)


    def get_v_absolute_value_for_more_pair(self):
        #因为计算耦合的分子对数超过了用于计算载流子迁移的分子对数，所以需要去掉不需要的对
        #具体做法是先获取pair的值
        #然后匹配
        pair_list = self.get_molpair(self.csv_filename, 8)
        print(pair_list)
        # 初始化状态
        in_running_block = False
        current_block = {}
        result_dict = {}
        result_dict_lumo = {}
        result_dict_homo = {}

        # 定义匹配 "RUNNING 数字" 格式的正则表达式
        pattern = re.compile(r'RUNNING (\d+)')

        # 读取文件内容
        with open(self.v_absolute_value_file, 'r') as file:
            for line in file:
                # 跳过无关数据
                if not in_running_block and not line.startswith('RUNNING'):
                    continue

                # 如果进入了 RUNNING 块
                if line.startswith('RUNNING'):
                    in_running_block = True

                    # 如果当前块不为空，将其添加到字典中
                    if current_block:
                        result_dict[int(current_block['number'])] = current_block['data']
                        current_block = {}

                    # 提取序号并将当前行添加到当前块中
                    match = pattern.match(line)
                    if match:
                        current_block['number'] = match.group(1)
                    current_block['data'] = [line.strip()]
                else:
                    # 将当前行添加到当前块中
                    current_block['data'].append(line.strip())

        # 将最后一个块添加到字典中
        if current_block:
            result_dict[int(current_block['number'])] = current_block['data']
        # 遍历 result_dict，改为 {serial_number: [v_absolute_value]}
        v_pattern = re.compile(r'J_eff\s+([-+]?\d*\.\d+(?:[eE][-+]?\d+)?)\s+eV')
        for key in result_dict.keys():
            # print(result_dict[key])
            first_match = 1
            for line in result_dict[key]:
                match = v_pattern.match(line)
                if match and first_match == 1:
                    result_dict_lumo[key] = float(match.group(1))
                    first_match = 0
                elif match and first_match == 0:
                    result_dict_homo[key] = float(match.group(1))
                    break

        #self.result_dict_homo = result_dict_homo
        #self.result_dict_lumo = result_dict_lumo
        # print(result_dict_homo)

        # 将 v_absolute_value 添加到 self.final_E_value_list 中
        for data in self.final_E_value_list:
            print('serial_number of data:', data['serial_number'])
            print('get_v_absolute_value, data in self.final_E_value_list: ', data)
            print('result_dict_lumo:', result_dict_lumo)
            for i, pair in enumerate(pair_list, start= 1):
                if data['pair_list'] == pair:
                    data['v_absolute_value_lumo'] = abs(result_dict_lumo[i])
                    data['v_absolute_value_homo'] = abs(result_dict_homo[i])
                    break
                else:
                    print('error find v value')

    def save_to_json(self,file_name):
        # 将结果写入Json文件
        with open(file_name, 'w') as file:
            json.dump(self.final_E_value_list, file, indent=4)
        print('数据提取的结果已保存到 result.json 文件中，请检查')

    def calculate_value_k(self, is_lumo):
        # 计算 k 值
        """
        h = 6.582119514 * 10 ** -16
        KbT = 0.025852
        lamda = 0.327166199682
        k = (2 pi / h) * ( v_absolute_value ** 2 / (4 * pi * lamda * KbT) ** 0.5 ) * exp(- (delta_E + lamda) ** 2 / (4 * lamda * KbT))
        @return:
        """
        h = 6.582119514 * 10 ** -16
        KbT = 0.025852
        lamda = 0.2633
        total_V_value = 0
        for data in self.final_E_value_list:
            print("data", data)
            print(data['v_absolute_value_lumo'])
            if is_lumo:
                v_absolute_value = data['v_absolute_value_lumo']
                delta_E = data['delta_ea']
            else:
                v_absolute_value = data['v_absolute_value_homo']
                delta_E = data['delta_ip']

            k = (2 * pi / h) * (v_absolute_value ** 2 / (4 * pi * lamda * KbT) ** 0.5) * exp(
                - (delta_E + lamda) ** 2 / (4 * lamda * KbT))
            # print('k值：', k, type(k))
            data['k_value'] = k
            total_V_value += v_absolute_value
        print('average V', total_V_value / len(self.final_E_value_list))

    def statistic_value_k(self):
        """
        统计 k 值在 小于 10^8, 在 10^8- 10^9 之间,  在 10^9 - 10^10 之间， 10^10 ~ 10^11 之间 、 10^11 - 10^12 之间, 和 大于 10^12 的数量以及百分比
        @return:
        """
        total_count = len(self.final_E_value_list)
        count_less_1e2 = 0
        count_1e2_1e3 = 0
        count_1e3_1e4 = 0
        count_1e4_1e5 = 0
        count_1e5_1e6 = 0

        count_1e6_1e7 = 0
        count_1e7_1e8 = 0
        count_1e8_1e9 = 0
        count_1e9_1e10 = 0
        count_1e10_1e11 = 0
        count_1e11_1e12 = 0
        count_greater_equal_1e12 = 0
        total_k_value = 0
        for data in self.final_E_value_list:
            total_k_value = total_k_value + data['k_value']
            if data['k_value'] < 1e2:
                count_less_1e2 += 1
                print(data['serial_number'])
            elif 1e2 <= data['k_value'] < 1e3:
                count_1e2_1e3 += 1
            elif 1e3 <= data['k_value'] < 1e4:
                count_1e3_1e4 += 1
            elif 1e4 <= data['k_value'] < 1e5:
                count_1e4_1e5 += 1
            elif 1e5 <= data['k_value'] < 1e6:
                count_1e5_1e6 += 1
            elif 1e6 <= data['k_value'] < 1e7:
                count_1e6_1e7 += 1
            elif 1e7 <= data['k_value'] < 1e8:
                count_1e7_1e8 += 1
            elif 1e8 <= data['k_value'] < 1e9:
                count_1e8_1e9 += 1
            elif 1e9 <= data['k_value'] < 1e10:
                count_1e9_1e10 += 1
            elif 1e10 <= data['k_value'] < 1e11:
                count_1e10_1e11 += 1
            elif 1e11 <= data['k_value'] < 1e12:
                count_1e11_1e12 += 1
            else:
                count_greater_equal_1e12 += 1
        print('k值的数量-total_count', total_count)
        print('k值小于 1e2 的数量：', count_less_1e2, '占比：', round(count_less_1e2 / total_count * 100, 2), '%')
        print('k值在 1e2 - 1e3 之间的数量：', count_1e2_1e3, '占比：', round(count_1e2_1e3 / total_count * 100, 2), '%')
        print('k值在 1e3 - 1e4 之间的数量：', count_1e3_1e4, '占比：', round(count_1e3_1e4 / total_count * 100, 2), '%')
        print('k值在 1e4 - 1e5 之间的数量：', count_1e4_1e5, '占比：', round(count_1e4_1e5 / total_count * 100, 2), '%')
        print('k值在 1e5 - 1e6 之间的数量：', count_1e5_1e6, '占比：', round(count_1e5_1e6 / total_count * 100, 2), '%')
        print('k值在 1e6 - 1e7 之间的数量：', count_1e6_1e7, '占比：', round(count_1e6_1e7 / total_count * 100, 2), '%')
        print('k值在 1e7 - 1e8 之间的数量：', count_1e7_1e8, '占比：', round(count_1e7_1e8 / total_count * 100, 2), '%')
        print('k值在 1e8 - 1e9 之间的数量：', count_1e8_1e9, '占比：', round(count_1e8_1e9 / total_count * 100, 2), '%')
        print('k值在 1e9 - 1e10 之间的数量：', count_1e9_1e10, '占比：', round(count_1e9_1e10 / total_count * 100, 2),
              '%')
        print('k值在 1e10 - 1e11 之间的数量：', count_1e10_1e11, '占比：', round(count_1e10_1e11 / total_count * 100, 2),
              '%')
        print('k值在 1e11 - 1e12 之间的数量：', count_1e11_1e12, '占比：', round(count_1e11_1e12 / total_count * 100, 2),
              '%')
        print('k值大于等于 1e12 的数量：', count_greater_equal_1e12, '占比：',
              round(count_greater_equal_1e12 / total_count * 100, 2), '%')
        print('average k value', total_k_value/total_count)

    def get_molpair(self, filename, close_contact_num):
        self.csv_filename = filename
        mol_close_contact = pd.read_csv(filename, names=range(1, 201))
        mol_close_contact.index += 1

        faceon_pair = []
        for mol_line in range(1, 201):
            for mol_cul in range(1, 201):
                if (mol_close_contact.at[mol_line, mol_cul] >= close_contact_num):
                    faceon_pair.append([mol_line, mol_cul])
       # self.pair_list = faceon_pair
        return faceon_pair

    def statistic_v_absolute_value(self, is_lumo):
        between_0_001_count = 0
        between_001_002_count = 0
        between_002_003_count = 0
        between_003_004_count = 0
        between_004_005_count = 0
        between_005_006_count = 0
        between_006_007_count = 0
        between_007_008_count = 0
        between_008_009_count = 0
        between_009_01_count = 0
        bigger_01_count = 0
        #data_cord = {}
        for data in self.final_E_value_list:
            if is_lumo:
                v_value = abs(data['v_absolute_value_lumo'])
            else:
                v_value = abs(data['v_absolute_value_homo'])

            if 0 <= v_value <= 0.01:
                between_0_001_count += 1
            elif 0.01 < v_value <= 0.02:
                between_001_002_count += 1
            elif 0.02 < v_value <= 0.03:
                between_002_003_count += 1
            elif 0.03 < v_value <= 0.04:
                between_003_004_count += 1
            elif 0.04 < v_value <= 0.05:
                between_004_005_count += 1
            elif 0.05 < v_value <= 0.06:
                between_005_006_count += 1
            elif 0.06 < v_value <= 0.07:
                between_006_007_count += 1
            elif 0.07 < v_value <= 0.08:
                between_007_008_count += 1
                '''
                data_cord.setdefault(data['serial_number'], data['pair_list'])
                print(data['serial_number'])
                print(data['pair_list'])
                print(v_value)
                '''
            elif 0.08 < v_value <= 0.09:
                between_008_009_count += 1
                '''
                data_cord.setdefault(data['serial_number'], data['pair_list'])
                print('0.08',data['serial_number'])
                print(data['pair_list'])
                print(v_value)
                '''
            elif 0.09 < v_value <= 0.1:
                between_009_01_count += 1
                '''
                data_cord.setdefault(data['serial_number'], data['pair_list'])
                print('0.09',data['serial_number'])
                print(data['pair_list'])
                print(v_value)
                '''
            else:
                bigger_01_count += 1
                '''
                data_cord.setdefault(data['serial_number'], data['pair_list'])
                print('0.1+',data['serial_number'])
                print(data['pair_list'])
                print(v_value)
                '''
        '''
        print(data_cord)
        data_pair = []
        data_mol = set()
        for key, value in data_cord.items():
            data_pair.append(str(key))
            data_mol.update(set(str(val) for val in value))
        print('--data_pair',' '.join(data_pair))
        print('--data_mol', ' '.join(data_mol))
        '''
        total_count = len(self.final_E_value_list)
        # 检查所有的数量是否正确
        if total_count != between_0_001_count + between_001_002_count + between_002_003_count + between_003_004_count + between_004_005_count + between_005_006_count + between_006_007_count + between_007_008_count + between_008_009_count + between_009_01_count + bigger_01_count:
            print(' statistic_v_absolute_value 计算的数量不正确，请检查')
            return
        print('between 0 001', between_0_001_count)
        print('between 001 002', between_001_002_count)
        print('between 002 003', between_002_003_count)
        print('between 003 004', between_003_004_count)
        print('between 004 005', between_004_005_count)
        print('between 005 006', between_005_006_count)
        print('between 006 007', between_006_007_count)
        print('between 007 008', between_007_008_count)
        print('between 008 009', between_008_009_count)
        print('between 009 010', between_009_01_count)
        print('between biger 01', bigger_01_count)

    def exec(self, filename):
        # self.get_A_file_list()
        self.get_ea_file_list()

        self.get_scf_value()
        print('exec, final_value_list_ea: ', self.final_value_list_ea)
        self.pair_list = self.get_molpair(filename, 8)
        self.calculate_delta_E2()
        self.get_v_absolute_value_for_more_pair()
        self.calculate_value_k(1)
        self.save_to_json('n3300_lumo0521.json')
        self.statistic_v_absolute_value(1)
        self.statistic_value_k()


def get_transport(mol_name):

    if str(mol_name) == 'y6150':
        FILE_DIR_y6 = r'D:\shujufenxi\Y6\disorder\y6\y6150'
        V_Absolute_Value_File_y6 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\y61500v.replace_result'
        filename_y6 = r'D:\shujufenxi\Y6\disorder\y6\v\data_15000000.csv'
        t1_y6 = FindSCF(FILE_DIR_y6, V_Absolute_Value_File_y6)
        t1_y6.exec(filename_y6)
    elif str(mol_name) == 'n3300':
        FILE_DIR_n3 = r'D:\shujufenxi\Y6\disorder\n32\n33000' #n32900
        V_Absolute_Value_File_n3 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\n3300.replace_result'
        #filename_n3 = r'D:\shujufenxi\Y6\disorder\n32\v\data_29900000.csv'
        filename_n3 = r'D:\shujufenxi\Y6\finaldata\n3zhenguishuju\data_30000000.csv'
        t1_n3 = FindSCF(FILE_DIR_n3,V_Absolute_Value_File_n3)
        t1_n3.exec(filename_n3)
    elif str(mol_name) == 'n4300':
        FILE_DIR_n4 = r'D:\shujufenxi\Y6\disorder\n42\n4299000'
        V_Absolute_Value_File_n4 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\n42990.replace_result'
        filename_n4 = r'D:\shujufenxi\Y6\disorder\n42\v\data_29900000.csv'
        t1_n4 = FindSCF(FILE_DIR_n4, V_Absolute_Value_File_n4)
        t1_n4.exec(filename_n4)
    else:
        print('cant recoginatize')

get_transport('n3300')


