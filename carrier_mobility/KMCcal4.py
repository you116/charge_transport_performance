import datetime
import os
import random
import numpy as np
import json
import matplotlib.pyplot as plt

def gen_neighbor_list():
    with open('data/n3300_lumo.json') as file:
        json_data = json.load(file)
        data_zero = np.zeros((200, 200))
        for data in json_data:
            data_zero[data['pair_list'][0] - 1, data['pair_list'][1] - 1] = data['k_value'] / 10 ** 12
            data_zero[data['pair_list'][1] - 1, data['pair_list'][0] - 1] = data['k_value'] / 10 ** 12
    return data_zero


def calcul_gai_lv_local(mol_data, nonzero_index_list, gai_lv_value):
    """
    计算概率分布， 最终跳跃到哪个分子？
    return  目标分子编号
    """
    jump_to_mol_index = None
    # 当前 mol_data 中的非0元素值
    value_list = []
    for index in nonzero_index_list[0]:
        value_list.append(mol_data[index])

    # 对value_list进行累加替换， 例如 [0.1, 0.2, 0.3] -> [0.1, 0.3, 0.6]， 后一个元素等于前一个元素与自身之和
    for i in range(1, len(value_list)):
        value_list[i] += value_list[i - 1]

    for i in range(len(value_list)):
        value_list[i] = value_list[i] / mol_data.sum()
    # 比较 gai_lv_value 与 value_list 中的元素，看它落在哪个区间
    for i in range(len(value_list)):
        if gai_lv_value < value_list[i]:
            jump_to_mol_index = nonzero_index_list[0][i]
            break
    return jump_to_mol_index


def calcul_carrier_mobility_process(neighbor_array, critical_time=10**6):
    """
    计算载流子迁移过程
    neighbor_array: 邻接矩阵， 每一行记录一个分子的电子跳跃到相邻分子的抽象概率，第一行为第一个分子，以此类推。
    若某一行第 i 列有值，则表示第 i 个分子是第一个分子的邻居，值为跳跃概率。
    每次跳跃的所需时间为该分子的邻居跳跃概率之和 分之一， 即此行数值之和分之一。
    当时间达到微秒级别时，返回最终跳跃位置。
    @return: res, start_index_copy, start_mol_index, total_time
    """
    # 随机选择一个分子作为跳跃起始点
    start_mol_index = start_index_copy = random.randint(0, 199)
    total_time = 0
    while total_time < critical_time:
        # 当前的分子 / 行
        current_mol = neighbor_array[start_mol_index]
        # 当前所需的跳跃时间
        jump_time = current_mol.sum()
        # 找出当前分子数组中的所有非0元素以及对应的列索引
        nonzero_index_list = np.nonzero(current_mol)
        # 生成一个概率随机数，0-1之间的浮点数
        gai_lv = round(random.uniform(0, 1), 6)
        # 计算概率分布
        # 从 nonzero_index_list 中随机选择一个非0元素作为下一个分子
        start_mol_index = calcul_gai_lv_local(current_mol, nonzero_index_list, gai_lv)
        if start_mol_index is None:
            print('error, start_mol_index is None')
            return False, start_index_copy, start_mol_index, total_time
        total_time += jump_time
    if start_mol_index == start_index_copy:
        print('return to the start point')
        return False, start_index_copy, start_mol_index, total_time

    return True, start_index_copy, start_mol_index, total_time


def calcul_mol_distance(mol_index_A, mol_index_B, coordinate_data_list):
    """
    计算分子距离, 使用 “分子中的环” 作为输入数据

    :param mol_index_A: 分子 A 的编号
    :param mol_index_B 分子 B 的编号
    :param coordinate_data_list  分子坐标数据
    """
    mol_A_center_point = None   # 分子 A 的中心点 {'x': 0.0, 'y': 0.0, 'z': 0.0}
    mol_B_center_point = None   # 分子 B 的中心点 {'x': 0.0, 'y': 0.0, 'z': 0.0}

    for mol_data in coordinate_data_list:
        mol_id = mol_data.get('mol_id', None)
        if mol_id-1 == mol_index_A:
            mol_A_ring_list = mol_data.get('ring_list', [])
            for ring_data in mol_A_ring_list:
                ring_index = ring_data.get('ring_index', None)
                if ring_index == 6:
                    mol_A_center_point = ring_data.get('center_point', None)
                    break
        elif mol_id-1 == mol_index_B:
            mol_B_ring_list = mol_data.get('ring_list', [])
            for ring_data in mol_B_ring_list:
                ring_index = ring_data.get('ring_index', None)
                if ring_index == 6:
                    mol_B_center_point = ring_data.get('center_point', None)
                    break

        if mol_A_center_point is not None and mol_B_center_point is not None:
            break

    if mol_A_center_point is None or mol_A_center_point is None:
        print('没有找到分子坐标数据！')
        return None

    # 计算两个分子中心点的距离
    distance = np.sqrt((mol_A_center_point['x'] - mol_B_center_point['x']) ** 2 +
                       (mol_A_center_point['y'] - mol_B_center_point['y']) ** 2 +
                       (mol_A_center_point['z'] - mol_B_center_point['z']) ** 2)
    return distance


def simulation(critical_time=10**6, simulation_times=5, coordinate_data_file_path='data/6_分子中的环.json', result_save_dir='data'):
    """
    :param critical_time: 临界时间
    :param simulation_times: 模拟次数
    :param coordinate_data_file_path: 分子坐标数据文件路径
    :param result_save_dir: 模拟结果保存路径

    结果是返回多次模拟 simulation_times 的平均值
    """
    neighbor_array = gen_neighbor_list()
    with open(coordinate_data_file_path) as file:
        coordinate_data_list = json.load(file)
    final_data = {
        'average_distance': None,
        'average_distance_square': None,
        'average_spend_time': None,
        'detail_data_list': []
    }
    detail_list = []
    i = 0
    while i < simulation_times:
        is_success, start_mol_id, end_mol_id, spend_time = calcul_carrier_mobility_process(neighbor_array, critical_time)
        if is_success:
            distance = calcul_mol_distance(start_mol_id, end_mol_id, coordinate_data_list)
            detail_list.append({
                'start_mol_id': int(start_mol_id),
                'end_mol_id': int(end_mol_id),
                'spend_time': float(spend_time),
                'distance': float(distance),
                'distance_square': float(distance ** 2)  # 距离的平方
            })
            i += 1
        else:
            print(f'第 {i+1} 次模拟失败！ 不要担心，继续模拟...')

    # 计算平均值
    distance_list = []
    distance_square_list = []
    spend_time_list = []
    for detail in detail_list:
        distance_list.append(detail['distance'])
        distance_square_list.append(detail['distance_square'])
        spend_time_list.append(detail['spend_time'])
    final_data['average_distance'] = float(np.mean(distance_list))
    final_data['average_distance_square'] = float(np.mean(distance_square_list))
    final_data['average_spend_time'] = float(np.mean(spend_time_list))

    # 将模拟结果写入文件,
    # 只有当 result_save_dir 不为空时 才保存到文件，
    # 只有当 result_save_dir 不为空时， 才记录 detail_list，否则将之 留空，以减少内存占用
    if result_save_dir:
        final_data['detail_data_list'] = detail_list
        file_name = f'{simulation_times}次_{critical_time}微秒{datetime.datetime.now().strftime("%Y%m%d%H%M%S")}.json'
        file_path = os.path.join(result_save_dir, file_name)
        with open(file_path, 'w') as file:
            json.dump(final_data, file, indent=4)
        print('模拟完成！文件名为：', file_name)

    return final_data


def analyze_and_draw(simulation_times, critical_time_list=None, coordinate_data_file_path='data/6_分子中的环.json', result_save_dir='data'):
    """
    分析和绘图
    :param critical_time_list: 临界时间列表
    :param simulation_times: 模拟次数
    :param coordinate_data_file_path: 分子坐标数据文件路径
    :param result_save_dir: 模拟结果保存路径
    """
    if critical_time_list is None:
        critical_time_list = [10 ** 6, 2 * (10 ** 6), 3 * (10 ** 6)]
    print(f'共有 {len(critical_time_list)} 个临界时间！')
    final_data_list = []
    for critical_time in critical_time_list:
        final_data = simulation(critical_time, simulation_times, coordinate_data_file_path, result_save_dir)
        final_data_list.append(final_data)
        print(f'临界时间为 {critical_time} 的模拟完成')

    # 绘图
    print('开始绘图...')
    x = critical_time_list
    y = []
    for final_data in final_data_list:
        y.append(final_data['average_distance_square'])
    plt.plot(x, y)
    plt.xlabel('time_spend')
    plt.ylabel('average_distance_square')
    plt.title('average_distance vs critical_time')
    plt.show()


def draw_Scatter_plot(simulation_times, critical_time_list, coordinate_data_file_path='data/6_分子中的环.json', result_save_dir='data'):
    """
    绘制散点图
    :param simulation_times: 模拟次数, 每次一个点
    :param critical_time_list: 临界时间，为 x 轴
    :param coordinate_data_file_path: 分子坐标数据文件路径
    :param result_save_dir: 模拟结果保存路径
    """
    final_data_list = []
    for critical_time in critical_time_list:
        final_data = simulation(critical_time, simulation_times, coordinate_data_file_path, result_save_dir)
        final_data_list.append(final_data)
        print(f'临界时间为 {critical_time} 的模拟完成')

    x_list = []
    y_list = []
    for final_data in final_data_list:
        detail_data_list = final_data.get('detail_data_list', [])
        for detail_data in detail_data_list:
            x_list.append(detail_data['spend_time'])
            y_list.append(detail_data['distance_square'])

    plt.scatter(x_list, y_list)
    plt.xlabel('time_spend')
    plt.ylabel('distance ** 2')
    plt.title('distance vs critical_time')
    plt.show()


if __name__ == '__main__':
    # analyze_and_draw(simulation_times=2, result_save_dir=None)

    draw_Scatter_plot(simulation_times=2, critical_time_list=[10 ** 6, 2 * (10 ** 6), 3 * (10 ** 6)], result_save_dir='data')
