import re
import os
import json
from math import pi, exp
import math
from numpy import array
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class Get_net_picture:

    def __init__(self, v_absolute_value_file, is_lumo = 1):
        self.result_dict_homo = {}
        self.result_dict_lumo = {}
        self.v_absolute_value_file = v_absolute_value_file
        self.is_lumo = is_lumo
        self.pair_list = []

    def get_molpair(self, filename):
        mol_close_contact = pd.read_csv(filename, names=range(1, 201))
        mol_close_contact.index += 1

        faceon_pair = []
        for mol_line in range(1, 201):
            for mol_cul in range(1, 201):
                if (mol_close_contact.at[mol_line, mol_cul] >= 8):
                    faceon_pair.append([mol_line, mol_cul])
        self.pair_list = faceon_pair


    def get_v_absolute_value(self):
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
            print(result_dict[key])
            first_match = 1
            for line in result_dict[key]:
                match = v_pattern.match(line)
                if match and first_match == 1:
                    result_dict_lumo[key] = float(match.group(1))
                    first_match = 0
                elif match and first_match == 0:
                    result_dict_homo[key] = float(match.group(1))
                    break

        self.result_dict_homo = result_dict_homo
        self.result_dict_lumo = result_dict_lumo
        print(result_dict_homo)

    def draw_net_picture_no_ea_value(self):
        data3 = np.zeros((200, 200))
        data4 = pd.DataFrame(data3, index=range(1, 201), columns=range(1, 201))
        result_dict = {}
        #print(len(self.final_E_value_list))
        #print(len(self.pair_list))
        if self.is_lumo:
            result_dict = self.result_dict_lumo
        else:
            result_dict = self.result_dict_homo

        for i, lists in enumerate(self.pair_list, start=1):
            data4.at[lists[0], lists[1]] = result_dict[i]
            data4.at[lists[1], lists[0]] = result_dict[i]
        data4.to_csv('name.csv')
        len_vslist = []
        v_max_list = []
        len_net_list = []
        num_list = []
        average_num = []
        x_list = np.arange(0.0, 0.052, 0.002)
        # x_list = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]
        for wq_x in x_list:
            data5 = abs(data4) > wq_x
            print("----mibe----", data5)
            average_num.append(data5.values.sum() / 200)
            V = list(data5.index[data5.any()])
            net = []
            while V:
                seed2 = [V[0]]
                seed = []
                while len(seed) != len(seed2):
                    # columns_l =
                    seed = seed2
                    print("----one", data5.loc[seed])
                    if len(seed) == 1:
                        columns_l = data5.loc[seed[0]]
                        #
                        print(columns_l)
                    else:
                        print("----two", data5.loc[seed])
                        columns_l = data5.loc[seed].any()
                    #
                    # print('---c',columns_list)

                    line_data = data5[columns_l]
                    print("======two", line_data)
                    columns_list = list(line_data.index)
                    print('---5', columns_list)
                    # print('----l', np.where(line_data.any()))
                    if isinstance(line_data, pd.Series):
                        seed2 = list(line_data[line_data].index)
                    else:
                        e = line_data.any()
                        seed2 = list(e[e].index)
                print(seed2)
                net.append(set(seed2 + columns_list))
                for h in set(seed2 + columns_list):
                    V.remove(h)
                print('net', net)
                print(V)
                print(seed)
            Vs = []
            Vmax = []
            for x in net:
                Vs.extend(x)
                Vmax.append(len(x))

            len_vslist.append(len(Vs))
            v_max_list.append(max(Vmax))
            len_net_list.append(len(net))
            num_list.append(
                (200 - len(Vs)) + len(net)
            )
        return {
            'len_vs_list': len_vslist,
            'v_max_list': v_max_list,
            'len_net_list': len_net_list,
            'num_list': num_list,
            'average_num': average_num
        }

    def exct(self, csv_file):
        self.get_molpair(csv_file)
        self.get_v_absolute_value()
        y_value = self.draw_net_picture_no_ea_value()
        return y_value

V_Absolute_Value_File_n3 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\n3300.replace_result'
#V_Absolute_Value_File_y6 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\y614975v.'
V_Absolute_Value_File_n4 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\n4300v.replace_result'
V_Absolute_Value_File_y62 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\y614950v.replace_result'
#V_Absolute_Value_File_n327 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\n3275.46737'
V_Absolute_Value_File_n329 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\n32990.replace_result'
V_Absolute_Value_File_y63 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\y61500v.replace_result'
V_Absolute_Value_File_n429 = r'D:\wq-data-statistic\Y6\ouheshuju\wb97xd\\n42990.replace_result'
filename_y6 = r'D:\shujufenxi\Y6\disorder\y6\v\data_14975000.csv'
filename_n3 = r'D:\shujufenxi\Y6\disorder\n32\v\data_30000000.csv'
filename_n4 = r'D:\shujufenxi\Y6\disorder\n42\v\data_30000000.csv'
filename_y62 = r'D:\shujufenxi\Y6\disorder\y6\v\data_14950000.csv'
filename_n327 = r'D:\shujufenxi\Y6\disorder\n3\v\n3300\data_30000000.csv'
filename_n329 = r'D:\shujufenxi\Y6\disorder\n32\v\data_29900000.csv'
filename_y63 = r'D:\shujufenxi\Y6\disorder\y6\v\data_15000000.csv'
filename_n429 = r'D:\shujufenxi\Y6\disorder\n42\v\data_29900000.csv'


is_lumo = 0

#t1_y6 = Get_net_picture(V_Absolute_Value_File_y6, is_lumo)
#y6_y_value = t1_y6.exct(filename_y6)
t1_n4 = Get_net_picture(V_Absolute_Value_File_n4, is_lumo)
n4_y_value = t1_n4.exct(filename_n4)
t1_y62 = Get_net_picture(V_Absolute_Value_File_y62, is_lumo)
y62_y_value = t1_y62.exct(filename_y62)
t1_n3 = Get_net_picture(V_Absolute_Value_File_n3, is_lumo)
n3_y_value = t1_n3.exct(filename_n3)
#t1_n327 = Get_net_picture(V_Absolute_Value_File_n327, is_lumo)
#n327_y_value = t1_n327.exct(filename_n327)
t1_n329 = Get_net_picture(V_Absolute_Value_File_n329, is_lumo)
n329_y_value = t1_n329.exct(filename_n329)
t1_y63 = Get_net_picture(V_Absolute_Value_File_y63, is_lumo)
y63_y_value = t1_y63.exct(filename_y63)
t1_n429 = Get_net_picture(V_Absolute_Value_File_n429, is_lumo)
n429_y_value = t1_n429.exct(filename_n429)


x_list = np.arange(0.0, 0.052, 0.002)
'''
num_y62 = y62_y_value['num_list']
#num_y6 = y6_y_value['num_list']
num_n4 = n4_y_value['num_list']
num_n3 = n3_y_value['num_list']
#num_n327 = n327_y_value['num_list']
num_n329 = n329_y_value['num_list']
num_y63 = y63_y_value['num_list']
num_n429 = n429_y_value['num_list']

max_y62 = y62_y_value['v_max_list']
#max_y6 = y6_y_value['v_max_list']
max_n4 = n4_y_value['v_max_list']
max_n3 = n3_y_value['v_max_list']
#max_n327 = n327_y_value['v_max_list']
max_n329 = n329_y_value['v_max_list']
max_y63 = y63_y_value['v_max_list']
max_n429 = n429_y_value['v_max_list']
average_y62 = y62_y_value['average_num']
#average_y6 = y6_y_value['average_num']
average_n4 = n4_y_value['average_num']
average_n3 = n3_y_value['average_num']
#average_n327 = n327_y_value['average_num']
average_n329 = n329_y_value['average_num']
average_y63 = y63_y_value['average_num']
average_n429 = n429_y_value['average_num']
'''
num_y6_total = [(x+y)/2 for x, y in zip(y62_y_value['num_list'],y63_y_value['num_list'])]
num_n3_total = [(x+y)/2 for x, y in zip(n3_y_value['num_list'],n329_y_value['num_list'])]
num_n4_total = [(x+y)/2 for x, y in zip(n4_y_value['num_list'], n429_y_value['num_list'])]

print('num_y6', num_y6_total)
print('num_n3', num_n3_total)
print('num_n4', num_n4_total)

max_y6_total = [(x+y)/2 for x, y in zip(y62_y_value['v_max_list'], y63_y_value['v_max_list'])]
max_n4_total = [(x+y)/2 for x, y in zip(n4_y_value['v_max_list'], n429_y_value['v_max_list'])]
max_n3_total = [(x+y)/2 for x, y in zip(n3_y_value['v_max_list'], n329_y_value['v_max_list'])]

print('max_y6', max_y6_total)
print('max_n3', max_n3_total)
print('max_n4', max_n4_total)

average_y6_total = [(x + y)/2 for x, y in zip(y62_y_value['average_num'], y63_y_value['average_num'])]
average_n4_total = [(x+y)/2 for x, y in zip(n4_y_value['average_num'], n429_y_value['average_num'])]
average_n3_total = [(x+y)/2 for x, y in zip(n3_y_value['average_num'], n329_y_value['average_num'])]

print('average_y6', average_y6_total)
print('average_n3', average_n3_total)
print('average_n4', average_n4_total)

# 将数据写入文件
import json
with open('0423_draw_data_LUMO0.json', 'w') as f:
    data_write = {
        'is_lumo': is_lumo,
        'num_y6': num_y6_total,
        'num_n3': num_n3_total,
        'num_n4': num_n4_total,
        'max_y6': max_y6_total,
        'max_n3': max_n3_total,
        'max_n4': max_n4_total,
        'average_y6': average_y6_total,
        'average_n3': average_n3_total,
        'average_n4': average_n4_total
    }
    json.dump(data_write, f)

    print('文件写入完成')


'''
plt.figure()
plt.plot(x_list, num_y6, 'bo-', label='Y6')
plt.plot(x_list, num_y62, 'ro-', label='y62-50')
plt.plot(x_list, num_n4, 'go-', label='N4')
plt.plot(x_list, num_n429, 'yo-', label='N3')
plt.plot(x_list, num_y63, 'mo-', label='y63-50')
plt.plot(x_list, num_n3, 'ko-', label='n3')
plt.xlabel('V(meV)')
plt.ylabel('number of cluster')
plt.xlabel('V(meV)')
plt.legend()

plt.figure()
plt.plot(x_list, max_y6, 'bo-', label='Y6')
plt.plot(x_list, max_y62, 'ro-', label='y62-50')
plt.plot(x_list, max_n4, 'go-', label='N4')
plt.plot(x_list, max_n429, 'yo-', label='N3')
plt.plot(x_list, max_y63, 'mo-', label='y63-50')
plt.plot(x_list, max_n3, 'ko-', label='n3')
plt.xlabel('V(meV)')
plt.ylabel('number of cluster')
plt.xlabel('V(meV)')
plt.legend()
plt.show()



t1_y62 = Get_net_picture(V_Absolute_Value_File_y62, is_lumo)
y62_y_value = t1_y62.exct(filename_y62)
num_y6_total = y62_y_value['num_list']
max_y6_total = y62_y_value['v_max_list']
print('num_y6_total', num_y6_total)
print('v_max_list', max_y6_total)
pd.DataFrame(num_y6_total).to_csv('num_y.csv')
pd.DataFrame(max_y6_total).to_csv('max.csv')


plt.figure()
plt.plot(x_list, num_y6_total, 'bo-', label='Y6')
#plt.plot(x_list, num_y62, 'ro-', label='y62-50')
#plt.plot(x_list, num_n4_total, 'go-', label='N4')
#plt.plot(x_list, num_n3_total, 'yo-', label='N3')
#plt.plot(x_list, num_y63, 'ko-', label='y63')
plt.xlabel('V(meV)')
plt.ylabel('number of cluster')
plt.xlabel('V(meV)')
plt.legend()

plt.figure()
plt.plot(x_list, max_y6_total, 'bo-', label='Y6')
#plt.plot(x_list, max_y62, 'ro-', label='y62-50')
#plt.plot(x_list, max_n4_total, 'go-', label='N4')
#plt.plot(x_list, max_n3_total, 'yo-', label='N3')
#plt.plot(x_list, max_y63, 'ko-', label='y63')
plt.xlabel('V(meV)')
plt.ylabel('max cluster size')
plt.xlabel('V(meV)')
plt.legend()
'''
plt.figure()
plt.plot(x_list, average_y6_total, 'bo-', label='Y6')
#plt.plot(x_list, average_y62, 'ro-', label='y62-50')
plt.plot(x_list, average_n4_total, 'go-', label='N4')
plt.plot(x_list, average_n3_total, 'yo-', label='N3')
#plt.plot(x_list, average_y63, 'ko-', label='y63')
plt.xlabel('V(meV)')
plt.ylabel('Nc(per molecule)')
plt.xlabel('V(meV)')
plt.legend()

plt.show()
#plt.savefig()




