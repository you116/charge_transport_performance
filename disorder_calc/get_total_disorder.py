import re
import os
import json
from math import pi, exp
import math
from numpy import array
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#FILE_DIR= 'D:\shujufenxi\Y6\disorder\y6\y6110'
class Get_total_disorder:
    """

    """
    def __init__(self):
        self.ALL_LOG_FILE_LIST = [f for f in os.listdir(FILE_DIR) if os.path.isfile(os.path.join(FILE_DIR, f)) and f.endswith('.log')]
        self.final_file_list_ea = []
        self.final_value_list_ea = []

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
                print('没有找到对应的eanion file',eacat_file_name)
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

        return final_file_list
    def get_scf_value(self, file_dict_list):
        final_value_list = []
        for file_dict in file_dict_list:
            enut_file_path = os.path.join(FILE_DIR, file_dict['enut'])
            eanion_file_path = os.path.join(FILE_DIR, file_dict['eanion'])
            eacat_file_path = os.path.join(FILE_DIR, file_dict['eacat'])
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
            i=0
            with open(eanion_file_path, 'r') as eanion_file:
                eanion_lines = eanion_file.readlines()
                i=i+1
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
            #读取cations file
            with open(eacat_file_path,'r') as eacat_file:
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
            value_dict['e_total_value'] = (value_dict['enut_value'] - value_dict['eanion_value'])*27.2114
            #计算i_total
            value_dict['i_total_value'] = (value_dict['eacat_value'] - value_dict['enut_value'])*27.2114
            final_value_list.append(value_dict)
        return final_value_list

    def save_to_json(self, filename):
        filename = str(filename)
        with open(filename, 'w') as file:
            json.dump(self.final_value_list_ea, file, indent=4)
        print('数据提取到了')

    def exec(self, out_filename):
        self.get_ea_file_list()
        self.final_value_list_ea = self.get_scf_value(self.final_file_list_ea)
        self.save_to_json(out_filename)

#执行部分，FILE_DIR存放IP或EA值得文件目录，将提取的值放入.json文件中
for value in ['y6110', 'y6115', 'y6120', 'y6125', 'y6130', 'y6135','y6140', 'y6145']:
    FILE_DIR = 'D:\shujufenxi\Y6\disorder\\y63\\'+ value
    t1 = Get_total_disorder()
    json_filename = 'D:\shujufenxi\Y6\disorder\\y63\\'+ value+'.json'
    t1.exec(json_filename)




