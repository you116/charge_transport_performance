"""
计算每个分子中的环的中心点坐标
按照图示，每个分子有 12 个环，从左到右 分别记为  ring_1, ring_2, ..., ring_12
手动录入环的原子编号，如 ring_1_atom_id_list = [1, 2, 3, 4, 5, 6, 7, 8, 9]
"""
import json
import os

import numpy as np

from config_public import RING_ATOM_ID_LIST


class Ring(object):
    """
    环
    """

    def __init__(self, molecular_coordinate_data=None, molecular_coordinate_file_path=None, frame_num=None, save_file_dir=None, is_save_process_file=False):
        if is_save_process_file and not save_file_dir:
            raise ValueError('保存文件的目录不能为空！')
        self.save_file_dir = save_file_dir
        self.is_save_process_file = is_save_process_file
        self.molecular_coordinate_data = molecular_coordinate_data
        self.molecular_coordinate_file_path = molecular_coordinate_file_path
        self.frame_num = frame_num
        self.origin_data = self.read_origin_data()
        self.final_ring_data = None   # 最终的环数据

    def read_origin_data(self):
        """
        读取原始数据，将原始数据保存到 self.origin_data 中
        :return:
        """
        if self.molecular_coordinate_data:
            return self.molecular_coordinate_data
        else:
            with open(self.molecular_coordinate_file_path, 'r') as f:
                origin_data = json.load(f)
            return origin_data

    @classmethod
    def calcul_center_point(cls, two_array):
        """
        求二维数组的中心点坐标，算数平均值，几何中心/质心
        :param two_array: 二维数组，一个环中所有原子的坐标，可能是 5 个，也可能是 6 个
        :return:
        """
        points = np.array(two_array)
        center_point = np.mean(points, axis=0)
        return center_point

    def ring_of_mol(self):
        """
        遍历每个分子，找出分子中的每个环，计算环的中心点坐标
        :return:
        """
        final_mol_data_list = []    # 最终的分子数据，在原始数据的基础上，添加了环的中心点坐标和环的原子编号

        for mol in self.origin_data:
            mol_atom_list = mol['atom_list']    # 分子中的原子列表
            ring_index = 1    # 环的编号，从 1 开始
            mol['ring_list'] = []    # 分子中的环列表

            for ring in RING_ATOM_ID_LIST:
                current_ring_atom_list = []    # 当前环中的原子列表
                for atom in mol_atom_list:
                    if atom['atom_id'] in ring:
                        current_ring_atom_list.append([float(atom['x']), float(atom['y']), float(atom['z'])])
                center_point = self.calcul_center_point(current_ring_atom_list)
                ring_data = {
                    'ring_index': ring_index,
                    'ring_atom_id_list': ring,
                    'center_point': {'x': center_point[0], 'y': center_point[1], 'z': center_point[2]}
                }
                mol['ring_list'].append(ring_data)
                ring_index += 1

            final_mol_data_list.append(mol)
        self.final_ring_data = final_mol_data_list
        return final_mol_data_list

    def save_to_file(self):
        """
        将 ring_mol_data_list 保存到文件中
        :return:
        """
        ring_file_path = os.path.join(self.save_file_dir, f'{self.frame_num}_分子中的环.json')

        with open(ring_file_path, 'w') as f:
            json.dump(self.final_ring_data, f, indent=4)

    def exec(self):
        self.ring_of_mol()
        if self.is_save_process_file:
            self.save_to_file()
        else:
            pass


if __name__ == '__main__':
    r = Ring()
    r.save_to_file()
