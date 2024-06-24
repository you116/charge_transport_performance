import pandas as pd
import numpy as np
import math
"""
在进行分子对提取并将其包装成高斯输入文件用于电子耦合的计算以及电子亲合能的计算
lammps输出文件中 atom_type 有关信息需要将其替换为具体的元素类型

在提取坐标时需要考虑周期性边界条件，此处的处理是如果分子上原子间的距离超过周期性盒子的一半，则视分子在周期性边界条件上，将其坐标归正
如果支持周期性边界条件的处理，则无需变动

从lammps中获取到分子坐标后，将其准备成高斯文件时，需要知道内存、核数、泛函基组等等信息，以及结尾时分子的连接信息。

考虑将烷基链缩短，需要C-C键换成C-H键
"""
#分子骨架上的原子序号
bone_atom_list = list(range(1, 32)) + list(range(36, 70)) + [78, 79, 80, 81]

#准备含有连接信息的坐标文件放入single.gjf和pair.gjf中
MAP_DICT = {1:'C',2:'C',3:'C',4:'N',5:'S',6:'C',7:'C',8:'N',9:'C',10:'C',11:'S',12:'C',13:'C',14:'S',15:'C',16:'C',17:'C',18:'C',19:'C',20:'C',21:'C',22:'C',23:'C',24:'C',25:'C',26:'O',27:'F',28:'C',29:'C',30:'N',31:'C',32:'N',33:'H',34:'H',35:'C',36:'C'}

#是否去除烷基链，
drop_alkyl_chain = True
#C-H对应的原子序号放入dict中
c_h_dict = {32:194, 33:197, 34:126, 35:157}
#去除烷基链后剩余的原子序号
list_left = list(range(1,82))

#lammps坐标文件
cord_filename_path = r'D:\daimazhengli\pair_extract\N4dump30000000.lammpstrj'

#含有分子上原子间距离信息文件----只计算EA、IP时不用提供
distance_filename_path = r'/'

#计算耦合时提供的文件有single.gjf, pair.gjf
single_file_path = r'D:\daimazhengli\pair_extract\single.gjf'
pair_file_path = r'D:\daimazhengli\pair_extract\pair.gjf'


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
        self.total_mol_num = None
        self.pair_list = None
        self.cord_data_no_c_chain = None

    def cal_distance(self,vec1, vec2):
        vec1 = np.array(vec1)
        vec2 = np.array(vec2)
        vec21 = vec1 - vec2
        distance = math.sqrt(np.dot(vec21, vec21))
        return distance
    def cord_map(self, data2):
        #map_dict =  {1:'C',2:'C',3:'C',4:'N',5:'S',6:'C',7:'C',8:'N',9:'C',10:'C',11:'S',12:'C',13:'C',14:'S',15:'C',16:'C',17:'C',18:'C',19:'C',20:'C',21:'C',22:'C',23:'C',24:'C',25:'C',26:'O',27:'F',28:'C',29:'C',30:'N',31:'C',32:'N',33:'H',34:'H',35:'C',36:'C'}
        data2['atom_type'] = data2['type'].map(MAP_DICT)
        return data2
    def extract_cord(self, filename):
        with open(filename, encoding='utf_8') as file:
            data = pd.read_table(file, sep='\\s+', skiprows=9,
                                 names=['atomid', 'mol', 'type', 'charge', 'x', 'y', 'z'])
        data.sort_values(by=['mol', 'atomid'], inplace=True, ascending=True)
        total_mol_num = data['mol'].max()
        self.total_mol_num  = total_mol_num
        self.atom_num_per_mol = int(data.iloc[-1,0]/total_mol_num)
        data['mol_atom_id'] = data['atomid'] - self.atom_num_per_mol * (data['mol'] - 1)
        data['mol_atom_id'] = data['mol_atom_id'].apply(int)
        data2 = data.set_index(['mol','mol_atom_id'], drop=False)
        #获取box_size
        with open(filename, encoding='utf_8') as file:
            data11 = pd.read_table(file, skiprows=0, nrows=9)
        data12 = data11.iloc[5, 0].split(' ')
        box_size = float(data12[1])-float(data12[0])
        mid_box = float(data12[0])+box_size/2

        #数据整理，周期性边界条件的处理
        for mol in range(1, total_mol_num+1):
            mid_var_x = -1
            mid_var_y = -1
            mid_var_z = -1
            # print(data2.loc[(mol,1)])
            for atomid in range(1, self.atom_num_per_mol+1):
                value = data2.loc[(mol, atomid)] - data2.loc[(mol, 1)]
                # print(value)
                if (value.x <= -box_size/2) | (value.x >= box_size/2):
                    mid_var_x = 1
                if (value.y <= -box_size/2) | (value.y >= box_size/2):
                    mid_var_y = 1
                if (value.z <= -box_size/2) | (value.z >= box_size/2):
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

    def calc_mol_pair(self, output_filename):
        total_mol_num = self.total_mol_num
        mol_close_contact = pd.DataFrame(np.zeros((total_mol_num, total_mol_num)), index=range(1, total_mol_num+1), columns=range(1, total_mol_num+1))
        print('start to get molecular pair')

        flag_j = []
        data2 = self.cord_data
        for mol in range(1, total_mol_num+1):
            for mol2 in range(mol, total_mol_num+1):
                if (mol == mol2): continue
                for atom1 in bone_atom_list:
                    vec1 = [data2.at[(mol, atom1), 'x'], data2.at[(mol, atom1), 'y'], data2.at[(mol, atom1), 'z']]
                    for atom2 in bone_atom_list:
                        if (mol_close_contact.at[mol, mol2] >= 25) & (flag_j == [mol, mol2]):
                            continue
                        vec2 = [data2.at[(mol2, atom2), 'x'], data2.at[(mol2, atom2), 'y'],
                                data2.at[(mol2, atom2), 'z']]
                        distance = self.cal_distance(vec1, vec2)
                        if distance <= 4:
                            mol_close_contact.at[mol, mol2] = mol_close_contact.at[mol, mol2] + 1
                            flag_j = [mol, mol2]
        mol_close_contact.to_csv(output_filename, index=False, header=False)
        print('get molecular pair finished')

    def get_mol_pair(self, filename, min_close_contact_num):
        mol_close_contact = pd.read_csv(filename, names=range(1, self.total_mol_num+1))
        mol_close_contact.index += 1
        faceon_pair = []
        for mol_line in range(1, self.total_mol_num+1):
            for mol_cul in range(1, self.total_mol_num+1):
                if (mol_close_contact.at[mol_line, mol_cul] >= min_close_contact_num):
                    faceon_pair.append([mol_line, mol_cul])
        self.pair_list = faceon_pair
    def change_atomcord(self,c_seq,h_seq, data3):
        data4 = data3.loc[h_seq, ['x', 'y', 'z']] - data3.loc[c_seq, ['x', 'y', 'z']]
        data4['o'] = np.sum(data4.apply(lambda x: x ** 2), axis='columns').apply(math.sqrt)
        data4['p'] = data4['o'].apply(lambda y: 1.09 / y)
        data4['x'] = data4['x'] * data4['p']
        data4['y'] = data4['y'] * data4['p']
        data4['z'] = data4['z'] * data4['p']
        data5 = data4[['x', 'y', 'z']] + data3.loc[c_seq, ['x', 'y', 'z']]
        # print('data4\n', type(data4))

        for mol in range(1, self.total_mol_num+1):
            data3.loc[(h_seq, mol), ['x', 'y', 'z']] = data5.loc[mol, ['x', 'y', 'z']]

    def drop_c_chain(self, c_h_dict):
        '''
        计算EA，和IP时，去掉分子对的烷基链，方便计算
        需要注意的是，去掉多余的链断掉的位置C变成H，键长也发生变化
        :return:
        '''
        #list_left = list(range(1,82))
        data = self.cord_data.set_index(['mol_atom_id', 'mol'],drop=False)
        for key, value in c_h_dict.items():
            self.change_atomcord(key, value, data)
            list_left.append(value)
        data_out = data.loc[list_left]
        data_out.loc[list(c_h_dict.values()), 'atom_type'] = 'H'
        data_out = data_out.swaplevel('mol_atom_id', 'mol').sort_index(level=0)
        self.cord_data_no_c_chain = data_out

    def file_starts(self, name,charge_str='0 1'):
        # chk = '%chk='+name.rstrip('.gjf')+'.chk'
        title = name
        part_one = pd.DataFrame([
            '%nprocshared=40',
            '%mem=120GB',
            # chk,
            '#  wb97xd/6-31g(d,p)  geom=connectivity  iop(3/107=0096400000,3/108=0096400000)',
            ' ',
            title,
            ' ',
            charge_str,
        ], columns=['atom'])
        return part_one
    def file_ends(self):
        with open(single_file_path, encoding='utf_8', mode='r') as fl:
            # data_starts = pd.read_table(fl, nrows=9, header=None)
            data_ends = pd.read_table(fl, names=['atom'])

        return data_ends
    def pair_file_ends(self):
        with open(pair_file_path, encoding='utf_8', mode='r') as fl:
            # data_starts = pd.read_table(fl, nrows=9, header=None)
            data_pair_ends = pd.read_table(fl, skiprows=380, names=['atom'])
        return data_pair_ends

    def get_total_mol(self, filename, charge_str='0 1'):
        total_pair = list(range(1, self.total_mol_num+1))
        for pair in total_pair:
            filename1 = str(pair)+str(filename)
            mol1 = self.cord_data_no_c_chain.loc[pair, ['atom_type', 'x', 'y', 'z']]
            gjf1 = pd.concat([self.file_starts(filename, charge_str), mol1, pd.DataFrame([np.NaN]), self.file_ends()])
            gjf1.to_csv(filename1, sep=' ', index=False, header=False, quotechar=' ')

    def get_pair_mol(self):
        total_pair = self.pair_list
        for i, pair in enumerate(total_pair,start=1):
            mol_pair = self.cord_data_no_c_chain.loc[pair, ['atom_type', 'x', 'y', 'z']]
            filename3 = str(i) + '_pair.gjf'
            gjf3 = pd.concat([self.file_starts(filename3), mol_pair, pd.DataFrame([np.NaN]), self.pair_file_ends()])
            gjf3.to_csv(filename3, sep=' ', index=False, header=False, quotechar=' ')
    def get_all_pair_mol(self, distance_cord, min_close_contact_num):
        self.get_mol_pair(distance_cord, min_close_contact_num)
        faceon_pair = self.pair_list
        data_out = self.cord_data_no_c_chain
        for i, pair in enumerate(faceon_pair, start=1):
            mol1 = data_out.loc[pair[0], ['atom_type', 'x', 'y', 'z']]
            mol2 = data_out.loc[pair[1],['atom_type', 'x' , 'y', 'z']]
            mol_pair = data_out.loc[pair, ['atom_type', 'x', 'y', 'z']]

            filename3 = str(i) + '_pair.gjf'
            filename1 = str(i) + '_A' + str(pair[0]) + '.gjf'
            filename2 = str(i) + '_B' + str(pair[1]) + '.gjf'

            gjf1 = pd.concat([self.file_starts(filename1), mol1, pd.DataFrame([np.NaN]), self.file_ends()])
            gjf2 = pd.concat([self.file_starts(filename2), mol2, pd.DataFrame([np.NaN]), self.file_ends()])
            gjf3 = pd.concat([self.file_starts(filename3), mol_pair, pd.DataFrame([np.NaN]), self.pair_file_ends()])

            gjf1.to_csv(filename1, sep=' ', index=False, header=False, quotechar=' ')
            gjf2.to_csv(filename2, sep=' ', index=False, header=False, quotechar=' ')
            gjf3.to_csv(filename3, sep=' ', index=False, header=False, quotechar=' ')
   # def exec(self, cord_filename=None, distance_filename=None):


#执行
t1 = Exec_dumpfile()
t1.extract_cord(cord_filename_path)

if drop_alkyl_chain:
    t1.drop_c_chain(c_h_dict)
#获取distance_filename的方式，在此处调用t1.calc_mol_pair(output_filename)
#输出用于计算耦合的文件时，请先调用get_mol_pair, 然后调用get_all_pair_mol, 并确保提供了含有分子间距离信息的文件distance_cord以及最小接触的原子的数量min_close_contact_num，示例： t1.get_all_pair_mol(distance_cord, min_close_contact_num)

#获取用于计算EA和IP值的高斯文件时，请调用get_total_mol，示例：t1.get_total_mol('_oni.gjf', charge_str = '-1 2')

t1.get_total_mol('_oni.gjf', charge_str='-1 2')
