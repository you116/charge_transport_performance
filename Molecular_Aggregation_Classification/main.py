"""
模块入口文件
"""
import os

from class_code.analyze import KmeansAnalyzer
from class_code.calculate_distance import CalculateDistance
from class_code.extract_data import DataExtract
from class_code.ring_of_mol import Ring

# 存放原始数据的目录
Origin_Data_Dir_Path = r'D:\wq-data-statistic\Molecular_Aggregation_Classification\n30517'
# 原始数据包含1. 坐标文件，文件格式: .txt; 内容：mol mol_atom_id x y z; 2. 成对信息，文件格式：.csv; 内容：包含成对信息的二维
# 是否保存处理过程文件
Is_Save_Process_File = True
# 保存过程处理文件的目录
Save_Process_File_Dir = r'D:\data_calculate_ZJY\N3\process3'
# 计算结果保存目录
Result_Save_Dir = r'D:\data_calculate_ZJY\N3\result5'
# 总聚类的数量
total_cluster_num = 6
# 是否输出灰度图
print_gray_and_color_picture = True


def check_dir():
    """
    检查 Origin_Data_Dir_Path 是否存在、是否是目录
    Save_Process_File_Dir 若不存在，则创建
    Result_Save_Dir 若不存在，则创建
    :return:
    """
    if not os.path.exists(Origin_Data_Dir_Path):
        raise ValueError(f'原始数据的目录 {Origin_Data_Dir_Path} 不存在！')
    if not os.path.isdir(Origin_Data_Dir_Path):
        raise ValueError(f'原始数据的目录 {Origin_Data_Dir_Path} 不是目录！')
    if Is_Save_Process_File:
        if not os.path.exists(Save_Process_File_Dir):
            os.makedirs(Save_Process_File_Dir)
            print(f'创建保存过程处理文件的目录 {Save_Process_File_Dir} 成功！')
    if not os.path.exists(Result_Save_Dir):
        os.makedirs(Result_Save_Dir)
        print(f'创建计算结果保存目录 {Result_Save_Dir} 成功！')
    print(f'目录检查通过！')


def check_origin_data():
    """
    检查原始数据, 遍历 Origin_Data_Dir_Path 下所有文件。
    必须包含 帧数_all_mol_data.txt、帧数_in_pairs.csv, 帧数为 1 - n 的整数
    :return:
    """
    file_list = os.listdir(Origin_Data_Dir_Path)
    frame_num_list_mol_data = []
    frame_num_list_pairs = []
    for file in file_list:
        if file.endswith('_all_mol_data.txt'):
            frame_num = int(file.split('_')[0])
            frame_num_list_mol_data.append(frame_num)
        elif file.endswith('_in_pairs.csv'):
            frame_num = int(file.split('_')[0])
            frame_num_list_pairs.append(frame_num)
    # mol_data 和 pairs 文件的帧数必须一一对应
    frame_num_list_mol_data.sort()
    frame_num_list_pairs.sort()
    if frame_num_list_mol_data != frame_num_list_pairs:
        raise ValueError(f'帧数_all_mol_data.txt 和 帧数_in_pairs.csv 不一一对应！')
    print(f'原始数据检查通过！')
    return frame_num_list_mol_data


def main():
    check_dir()
    frame_num_list = check_origin_data()
    file_path_list = []
    for frame_num in frame_num_list:
        mol_file_path = os.path.join(Origin_Data_Dir_Path, f'{frame_num}_all_mol_data.txt')
        pair_file_path = os.path.join(Origin_Data_Dir_Path, f'{frame_num}_in_pairs.csv')
        data = {
            'coord_file_path': mol_file_path,
            'pair_file_path': pair_file_path,
            'frame_num': frame_num,
        }
        file_path_list.append(data)

    print(f'原始文件读取完成，共有 {len(frame_num_list)} 帧的数据！')

    # 最后合并一次本批数据，生成一张大的 kmeans 图片
    merge_kmeans_data_list = []

    for file_path_data in file_path_list:
        print(f'正在处理帧数 {file_path_data.get("frame_num")} 的数据...')
        # 基础数据提取
        base_data_extract = DataExtract(
            coordinate_file_path=file_path_data.get('coord_file_path'),
            in_pair_file_path=file_path_data.get('pair_file_path'),
            frame_num=file_path_data.get('frame_num'),
            save_file_dir=Save_Process_File_Dir,
            is_save_process_file=Is_Save_Process_File,
        )
        base_data_extract.exec()
        print(f'{file_path_data.get("frame_num")} 的基础数据处理完成！')

        # 环数据提取
        ring_data_extract = Ring(
            molecular_coordinate_data=base_data_extract.molecular_coordinate_data,
            molecular_coordinate_file_path=None,
            frame_num=file_path_data.get('frame_num'),
            save_file_dir=Save_Process_File_Dir,
            is_save_process_file=Is_Save_Process_File,
        )
        ring_data_extract.exec()
        print(f'{file_path_data.get("frame_num")} 的环数据处理完成！')

        # 距离计算
        distance_calculator = CalculateDistance(
            ring_file_path=None,
            pair_file_path=None,
            ring_data=ring_data_extract.final_ring_data,
            pair_data=base_data_extract.in_pair_data,
            frame_num=file_path_data.get('frame_num'),
            save_file_dir=Save_Process_File_Dir,
            is_save_process_file=Is_Save_Process_File,
        )
        distance_calculator.exec()
        print(f'{file_path_data.get("frame_num")} 的距离计算完成！')
        # merge_kmeans_data_list.append(distance_calculator.kmeans_data.get('kmeans_data'))
        # 分析处理6
        analyzer = KmeansAnalyzer(
            data_file_path=None,
            data_input=distance_calculator.kmeans_data,
            frame_num=file_path_data.get('frame_num'),
            result_save_dir=Result_Save_Dir,
            is_draw_gray_and_color=print_gray_and_color_picture,
        )
        analyzer.exec()
        merge_kmeans_data_list.append(analyzer.final_kmeans_data)

    print(f'每一帧数据处理完成')

    # 最后合并一次本批数据，生成一张大的 kmeans 图片
    KmeansAnalyzer.merge_kmeans(merge_kmeans_data_list, Result_Save_Dir, n_clusters=total_cluster_num)


if __name__ == '__main__':
    main()
