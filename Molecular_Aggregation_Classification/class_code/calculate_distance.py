"""
计算成对的环的距离
"""
import json
import os

import numpy as np
#


class CalculateDistance:

    def __init__(self, ring_file_path=None, pair_file_path=None, ring_data=None, pair_data=None, frame_num=None, save_file_dir=None, is_save_process_file=False):
        if is_save_process_file and not save_file_dir:
            raise ValueError('保存文件的目录不能为空！')
        self.ring_data = ring_data
        self.pair_data = pair_data
        self.frame_num = frame_num
        self.save_file_dir = save_file_dir
        self.is_save_process_file = is_save_process_file
        self.ring_file_path = ring_file_path
        self.pair_file_path = pair_file_path
        self.ring_data_list = self.read_ring_file()
        self.pair_list = self.read_pair_file()
        self.further_pair_data = None   # 进一步整理成对的环数据
        self.calculate_distance_data = None     # 计算的距离数据
        self.kmeans_data = None    # KMeans 算法需要的数据

    def read_pair_file(self):
        """
        读取成对情况文件
        :return:
        """
        if self.pair_data:
            return self.pair_data.get('partners_data')
        else:
            with open(self.pair_file_path, 'r') as f:
                pair_list = json.load(f).get('partners_data')
            return pair_list

    def read_ring_file(self):
        if self.ring_data:
            return self.ring_data
        else:
            with open(self.ring_file_path, 'r') as f:
                ring_data_list = json.load(f)
            return ring_data_list

    def save_to_file(self):
        further_pair_file_path = os.path.join(self.save_file_dir, f'{self.frame_num}_进一步整理成对的环.json')
        with open(further_pair_file_path, 'w') as f:
            json.dump(self.further_pair_data, f, indent=4)
        calculate_file_path = os.path.join(self.save_file_dir, f'{self.frame_num}_计算的距离数据.json')
        with open(calculate_file_path, 'w') as f:
            json.dump(self.calculate_distance_data, f, indent=4)
        kmeans_file_path = os.path.join(self.save_file_dir, f'{self.frame_num}_KMeans数据.json')
        with open(kmeans_file_path, 'w') as f:
            json.dump(self.kmeans_data, f, indent=4)

    def build_data(self):
        """
        构建数据
        :return:
        """
        final_data_list = []
        for pair_data in self.pair_list:
            main_mol_id = pair_data.get('mol_id')    # 当前的分子编号，称为 “主分子”
            partner_list = pair_data.get('partner_list')    # 与它成对的分子编号列表, 称为 “伙伴分子”

            # 找到主分子的环列表
            main_mol_ring_list = None   # 主分子的环列表, 里面有12个环的数据
            for ring_data in self.ring_data_list:
                if ring_data.get('mol_id') == main_mol_id:
                    main_mol_ring_list = ring_data.get('ring_list')
                    break
            if not main_mol_ring_list:
                print(f'没有找到 主分子 {main_mol_id} 的环列表数据！')
                return None

            # 找到伙伴分子的环列表
            partner_ring_data_list = []    # 伙伴分子的环列表, 里面有每个伙伴分子的12个环的数据
            for partner_mol_id in partner_list:
                partner_mol_ring_data = {
                    'partner_mol_id': partner_mol_id,
                    'ring_list': None
                }
                for ring_data in self.ring_data_list:
                    if ring_data.get('mol_id') == partner_mol_id:
                        partner_mol_ring_data['ring_list'] = ring_data.get('ring_list')
                        break
                if not partner_mol_ring_data['ring_list']:
                    print(f'没有找到 伙伴分子 {partner_mol_id} 的环列表数据！')
                    return None
                partner_ring_data_list.append(partner_mol_ring_data)

            final_data = {
                'main_mol_id': main_mol_id,   # 主分子的编号
                'main_mol_ring_list': main_mol_ring_list,    # 主分子的环列表
                'partner_mol_id_list': partner_list,    # 伙伴分子的编号列表
                'partner_ring_data_list': partner_ring_data_list    # 伙伴分子的环列表
            }
            final_data_list.append(final_data)

        self.further_pair_data = final_data_list
        return final_data_list

    @staticmethod
    def calcul_two_point_distance(point1, point2):
        """
        计算两个点的距离
        :param point1: {'x': 1, 'y': 2, 'z': 3}     # 点1的坐标
        :param point2: {'x': 4, 'y': 5, 'z': 6}     # 点2的坐标
        :return: number   # 两个点的距离
        """
        x1, y1, z1 = point1['x'], point1['y'], point1['z']
        x2, y2, z2 = point2['x'], point2['y'], point2['z']
        A = np.array([x1, y1, z1])
        B = np.array([x2, y2, z2])
        distance = np.linalg.norm(B-A)
        return distance

    def calculate(self):
        # if calculate_file_path:
        #     with open(calculate_file_path, 'r') as f:
        #         calculate_data_list = json.load(f)
        # else:
        #     calculate_data_list = self.build_data()

        final_data_list = []
        for calcul_data in self.further_pair_data:

            main_mol_id = calcul_data.get('main_mol_id')    # 主分子的编号
            main_mol_ring_list = calcul_data.get('main_mol_ring_list')  # 主分子的环列表
            partner_mol_id_list = calcul_data.get('partner_mol_id_list')    # 伙伴分子的编号列表
            partner_ring_data_list = calcul_data.get('partner_ring_data_list')  # 伙伴分子的环列表

            distance_data = {
                'main_mol_id': main_mol_id,
                'partner_mol_id_list': partner_mol_id_list,
                'distance_list': []  # 距离列表
            }

            # 遍历主分子的每个环，计算它与伙伴分子的每个环的距离
            for main_ring in main_mol_ring_list:
                main_ring_point = main_ring.get('center_point')    # 主分子的环的中心点坐标
                for partner_data in partner_ring_data_list:
                    partner_mol_id = partner_data.get('partner_mol_id')
                    partner_ring_list = partner_data.get('ring_list')
                    distance_detail = {
                        'main_mol_id': main_mol_id,
                        'main_ring_index': main_ring.get('ring_index'),    # 主分子的环的编号，从 1 开始, 代表当前的距离数据是主分子的第几个环
                        'partner_mol_id': partner_mol_id,
                        'distance': []    # 距离列表
                    }
                    for partner_ring in partner_ring_list:
                        partner_ring_point = partner_ring.get('center_point')
                        distance = self.calcul_two_point_distance(main_ring_point, partner_ring_point)
                        distance_detail['distance'].append(distance)
                    distance_data['distance_list'].append(distance_detail)
            final_data_list.append(distance_data)
        self.calculate_distance_data = final_data_list
        return final_data_list

    def conv_distance_data(self, distance_data_file_path=None):
        """
        将计算的距离数据转换成 KMeans 聚类算法需要的数据格式
        kmeans 的基本单元是分子对
        最终的数据格式是：
        human: Json 格式，便于人眼查看
        { pair: [1, 31], distance: [ [ 123, 234, 2312, 3423 ], [ 123, 234, 2312, 3423 ] ] }
        pair: [1, 31] 代表分子对的编号, 只能有两个
        distance 是一个二维数组，代表两个分子的每个环的距离
        kmeans: 纯数组格式
        [ [ 123, 234, 2312, 3423 ], [ 123, 234, 2312, 3423 ] ]
        """
        final_data = {
            'human_data': [],   # 人工标注的数据, 便于人眼查看
            'kmeans_data': []   # kmeans 算法需要的数据， 纯数组
        }
        # if distance_data_file_path:
        #     with open(distance_data_file_path, 'r') as f:
        #         distance_data_list = json.load(f)
        # else:
        #     distance_data_list = self.calculate()

        for distance_data in self.calculate_distance_data:
            main_mol_id = distance_data.get('main_mol_id')
            partner_mol_id_list = distance_data.get('partner_mol_id_list')
            distance_list = distance_data.get('distance_list')

            for partner_id in partner_mol_id_list:
                human_data = {
                    'pair': [main_mol_id, partner_id],
                    'distance': []
                }
                kmeans_data = []
                for distance_detail in distance_list:
                    if distance_detail.get('partner_mol_id') == partner_id:
                        human_data['distance'].append(distance_detail.get('distance'))
                        kmeans_data.append(distance_detail.get('distance'))
                final_data['human_data'].append(human_data)
                final_data['kmeans_data'].append(kmeans_data)

        self.kmeans_data = final_data
        return final_data

    def exec(self):
        """
        执行
        :return:
        """
        self.build_data()
        self.calculate()
        self.conv_distance_data()
        if self.is_save_process_file:
            self.save_to_file()


if __name__ == '__main__':
    # 读取环的文件路径
    Ring_File_Path = r'D:\Code\wq-data-statistic\Molecular_Aggregation_Classification\process_file\1_分子中的环.json'
    # 读取成对情况文件路径
    Pair_File_Path = r'D:\Code\wq-data-statistic\Molecular_Aggregation_Classification\process_file\1_成对情况.json'
    cd = CalculateDistance(ring_file_path=Ring_File_Path, pair_file_path=Pair_File_Path, frame_num=99, save_file_dir=r'D:\Code\wq-data-statistic\Molecular_Aggregation_Classification\process_file', is_save_process_file=True)
    cd.exec()
    print('Done!')
