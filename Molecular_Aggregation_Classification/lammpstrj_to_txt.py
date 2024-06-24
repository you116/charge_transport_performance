import numpy as np
import math
import pandas as pd


class Exec_dumpfile:
    '''
    打开dump file文件，对数据进行处理，超出周期性边界条件的原子归位
    数据清洗，需要下一步的计算有，临近元素计算，生成分子对坐标的gjf文件，临近分子的V计算，
    absolute E计算
    提前改变的参数，盒子大小，省略时原子的位置，
    self.cord_data最终清理后的数据，包括周期性边界条件的处理，以及原子序号与元素的映射
    不同分子计算mol的方式不同，具体改变calc_mol_pair方法
    '''
    def __init__(self):
        self.cord_data = None
        self.atom_num_per_mol = None
        self.pair_list = None
        self.cord_data_no_c_chain = None

    def cal_distance(self,vec1, vec2):
        vec1 = np.array(vec1)
        vec2 = np.array(vec2)
        vec21 = vec1 - vec2
        distance = math.sqrt(np.dot(vec21, vec21))
        return distance
    def cord_map(self, data2):
        map_dict =  {1:'C',2:'C',3:'C',4:'N',5:'S',6:'C',7:'C',8:'N',9:'C',10:'C',11:'S',12:'C',13:'C',14:'S',15:'C',16:'C',17:'C',18:'C',19:'C',20:'C',21:'C',22:'C',23:'C',24:'C',25:'C',26:'O',27:'F',28:'C',29:'C',30:'N',31:'C',32:'N',33:'H',34:'H',35:'C',36:'C'}
        data2['atom_type'] = data2['type'].map(map_dict)
        return data2
    def extract_cord(self, filename, mid_box):
        with open(filename, encoding='utf_8') as file:
            data = pd.read_table(file, sep='\\s+', skiprows=9,
                                 names=['atomid', 'mol', 'type', 'charge', 'x', 'y', 'z'])
        data.sort_values('atomid', inplace=True, ascending=True)
        self.atom_num_per_mol = int(data.iloc[-1,0]/200)
        data['mol_atom_id'] = data['atomid'] - self.atom_num_per_mol * (data['mol'] - 1)
        data['mol_atom_id'] = data['mol_atom_id'].apply(int)
        data2 = data.set_index(['mol','mol_atom_id'], drop=False)
        #获取box_size
        with open(filename, encoding='utf_8') as file:
            data11 = pd.read_table(file, skiprows=0, nrows=9)
        data12 = data11.iloc[5, 0].split(' ')
        box_size = float(data12[1])-float(data12[0])

        #数据整理，周期性边界条件的处理
        for mol in range(1, 201):
            mid_var_x = -1
            mid_var_y = -1
            mid_var_z = -1
            # print(data2.loc[(mol,1)])
            for atomid in range(1, self.atom_num_per_mol+1):
                value = data2.loc[(mol, atomid)] - data2.loc[(mol, 1)]
                # print(value)
                if (value.x <= -50) | (value.x >= 50):
                    mid_var_x = 1
                if (value.y <= -50) | (value.y >= 50):
                    mid_var_y = 1
                if (value.z <= -50) | (value.z >= 50):
                    mid_var_z = 1
            for atomid in range(1, self.atom_num_per_mol+1):
                if (mid_var_x == 1):
                    if (data2.loc[(mol, atomid), 'x'] <= mid_box):
                        data2.loc[(mol, atomid), 'x'] = data2.loc[(mol, atomid), 'x'] + box_size
                if (mid_var_y == 1):
                    if (data2.loc[(mol, atomid), 'y'] <= mid_box):
                        data2.loc[(mol, atomid), 'y'] = data2.loc[(mol, atomid), 'y'] + box_size
                if (mid_var_z == 1):
                    if (data2.loc[(mol, atomid), 'z'] <= mid_box):
                        data2.loc[(mol, atomid), 'z'] = data2.loc[(mol, atomid), 'z'] + box_size
        self.cord_data = self.cord_map(data2)
for num in range(30000000,30050000,50000):
    t1 = Exec_dumpfile()
    t1.extract_cord('D:\shujufenxi\Y6\\finaldata\\n3lammps3\\N3dump'+str(num)+'.lammpstrj', 85)
    data = t1.cord_data[['x','y','z']]
    data.to_csv('D:\wq-data-statistic\Molecular_Aggregation_Classification\\N3\\'+'_all_mol_data'+'.txt')
    print('finished calculation', num)

