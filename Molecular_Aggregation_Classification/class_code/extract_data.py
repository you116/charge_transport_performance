"""
提取数据
./molecular_data/1_all_mol_data.txt   每个分子中所有原子的坐标数据
./molecular_data/in_pair.csv        判断分子成对情况
"""
import json
import os


class DataExtract(object):

    def __init__(self, coordinate_file_path=None, in_pair_file_path=None, frame_num=None, save_file_dir=None, is_save_process_file=False):

        if is_save_process_file and not save_file_dir:
            raise ValueError('保存文件的目录不能为空！')

        # 原子坐标数据，以原子为单位
        #   {
        #         'mol_id': line_list[0],
        #         'atom_id': line_list[1],
        #         'x': line_list[2],
        #         'y': line_list[3],
        #         'z': line_list[4],
        #     }
        self.frame_num = frame_num  # 数据帧编号
        self.save_file_dir = save_file_dir  # 保存文件的目录
        self.is_save_process_file = is_save_process_file
        self.coordinate_file_path = coordinate_file_path
        self.in_pair_file_path = in_pair_file_path
        self.all_atom_coordinate_data = None    # 所有原子坐标数据，以原子为单位
        # 分子坐标数据，以分子为单位
        #   {
        #       mol_id: 1,
        #       atom_list: [
        #           {
        #               'atom_id': 1,
        #               'x': 123,
        #               'y': 879374,
        #               'z': 2343,
        #           },
        #           ...
        #       ]
        self.molecular_coordinate_data = None   # 按分子分类的原子坐标数据
        # 成对情况数据
        #   {
        #       pair_count: 123,                                # 有多少对
        #       meet_stand_count: 123,                     # 所有达到成对要求的分子数量
        #       pairs: [ [1, 2], [3, 4], ... ],         # 成对的分子
        #       meet_stand_mol_ids: [1, 2, 3, 4, ...]              # 所有达到成对要求的分子
        #   }
        self.in_pair_data = None    # 成对情况数据

    def build_all_atom_coordinate_data(self):
        """
        读取坐标数据
        :return:
        """
        with open(self.coordinate_file_path, 'r') as f:
            lines = f.readlines()[1:]

        final_data = []
        for line in lines:
            line = line.strip()
            line_list = line.split(',')
            atom_data = {
                'mol_id': int(line_list[0]),
                'atom_id': int(line_list[1]),
                'x': line_list[2],
                'y': line_list[3],
                'z': line_list[4],
            }
            final_data.append(atom_data)
        self.all_atom_coordinate_data = final_data

    def build_molecular_coordinate_data(self):
        """
        从 all_atom_coordinate_data 中构建 molecular_coordinate_data
        :return:
        """
        if not self.all_atom_coordinate_data:
            self.build_all_atom_coordinate_data()

        final_data = []
        mol_id_set = set()
        for atom_data in self.all_atom_coordinate_data:
            mol_id_set.add(atom_data['mol_id'])

        for mol_id in mol_id_set:
            mol_data = {
                'mol_id': mol_id,
                'atom_list': []
            }
            for atom_data in self.all_atom_coordinate_data:
                if atom_data['mol_id'] == mol_id:
                    mol_data['atom_list'].append(atom_data)
            final_data.append(mol_data)
        self.molecular_coordinate_data = final_data

    def save_to_file(self):
        """
        将 all_atom_coordinate_data 和 molecular_coordinate_data 保存到 Json 文件中
        """
        all_atom_coord_file_path = os.path.join(self.save_file_dir, f'{self.frame_num}_所有原子坐标.json')
        atom_coord_by_mol_file_path = os.path.join(self.save_file_dir, f'{self.frame_num}_按分子分类的原子坐标.json')
        in_pair_file_path = os.path.join(self.save_file_dir, f'{self.frame_num}_成对情况.json')

        # 检测保存文件的目录是否存在，不存在则创建
        if not os.path.exists(self.save_file_dir):
            os.makedirs(self.save_file_dir)

        with open(all_atom_coord_file_path, 'w') as f:
            json.dump(self.all_atom_coordinate_data, f, indent=4)

        with open(atom_coord_by_mol_file_path, 'w') as f:
            json.dump(self.molecular_coordinate_data, f, indent=4)

        with open(in_pair_file_path, 'w') as f:
            json.dump(self.in_pair_data, f, indent=4)

    def build_in_pair_data(self):
        """
        读取成对情况数据 从 IN_PAIR_FILE_PATH 中
        将 IN_PAIR_FILE_PATH 文件读成一个二维数组
        遍历二维数组，当 第 x 行 第 y 列 的值 >= 25 时，表示 分子x 和 y 成对
        """
        with open(self.in_pair_file_path, 'r') as f:
            lines = f.readlines()

        final_data = []
        for line in lines:
            line = line.strip()
            line_list = line.split(',')
            final_data.append(line_list)

        # 遍历 final_data
        #   当 final_data[i][j] >= 25 时，表示 分子 i 和 j 成对
        detail_data = []
        all_data = []
        for i in range(len(final_data)):
            for j in range(len(final_data[i])):
                if i == j:
                    continue
                if float(final_data[i][j]) >= 25:
                    # print(f'分子{i+1} 和 分子{j+1} 成对, 值为 {final_data[i][j]}')
                    detail_data.append([i+1, j+1])
                    all_data.append(i+1)
                    all_data.append(j+1)

        self.in_pair_data = {
            'pair_count': len(detail_data),     # 有多少对
            'meet_stand_count': len(set(all_data)),     # 所有达到成对要求的分子数量
            'pairs': detail_data,   # 成对的分子
            'meet_stand_mol_ids': list(set(all_data)),  # 所有达到成对要求的分子
            'partners_data': [],    # 每个分子的配对分子， { 'mol_id': 1, 'partner_list': [2, 3, 4, ...] }
        }
        # partner 记录键值形式的每一个分子的配对分子， { 'mol_id': 1, 'partner_list': [2, 3, 4, ...] }
        # 去重，将 [1, 2] 和 [2, 1] 视为同一对
        final_pairs = []
        for pair in self.in_pair_data['pairs']:
            if [pair[1], pair[0]] not in final_pairs:
                final_pairs.append(pair)
        self.in_pair_data['pairs'] = final_pairs

        partner = {}
        for mol_id in self.in_pair_data['meet_stand_mol_ids']:
            partner[mol_id] = []
        for pair in self.in_pair_data['pairs']:
            partner[pair[0]].append(pair[1])
            partner[pair[1]].append(pair[0])
        for mol_id in partner:
            self.in_pair_data['partners_data'].append({
                'mol_id': mol_id,
                'partner_list': partner[mol_id]
            })

        # partners_data 去重， mol_id: 3, partner_list: [4, 5, 6], 则在 partners_data 中， mol_id: 4, partner_list: [3, 5, 6] 中删除 3
        for partner_data in self.in_pair_data['partners_data']:
            for partner_id in partner_data['partner_list']:
                for partner_data_2 in self.in_pair_data['partners_data']:
                    if partner_data_2['mol_id'] == partner_id:
                        if partner_data['mol_id'] in partner_data_2['partner_list']:
                            partner_data_2['partner_list'].remove(partner_data['mol_id'])

        # 最终删掉 partners_data 中 partner_list 为空的数据
        final_partners_data = []
        for partner_data in self.in_pair_data['partners_data']:
            if partner_data['partner_list']:
                final_partners_data.append(partner_data)
        self.in_pair_data['partners_data'] = final_partners_data

    def exec(self):
        self.build_all_atom_coordinate_data()
        self.build_molecular_coordinate_data()
        self.build_in_pair_data()
        if self.is_save_process_file:
            self.save_to_file()
        else:
            pass


if __name__ == '__main__':
    de = DataExtract(
        coordinate_file_path=r'/Molecular_Aggregation_Classification/molecular_data/1_all_mol_data.txt',
        in_pair_file_path=r'/Molecular_Aggregation_Classification/molecular_data/1_in_pairs.csv',
        frame_num=22,
        save_file_dir=r'/Molecular_Aggregation_Classification/test_save_dir',
        is_save_process_file=True
    )
    de.save_to_file()


