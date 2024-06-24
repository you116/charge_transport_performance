# 分子聚合分类

## 处理顺序
1， extract_data.py 进行最初的数据处理 ，生成 原子坐标.json、按分子分类的原子坐标.json、成对情况.json

2， ring_of_mol.py 处理环相关的数据，生成  分子中的环.json

3， calculate_distance.py   计算距离， 先处理 环数据，生成 “进一步整理成对的环.json”  和 最终的距离数据  “计算的距离数据.json”

---
代码中有详细注释


## 用法 Usage

- 将所有的原始坐标文件，按照示例的格式，放在 一个目录中，必须是 xx_all_mol_data.txt 和 xx_in_pairs.csv 格式，xx 代表帧数编号，为阿拉伯数字。
- 即，每个帧都应有一个坐标文件 xx_all_mol_data.txt 和分子对 xx_in_pairs.csv 文件。 详见 molecular_data/ 目录的示例
- 示例的提取lammps文件的脚本为lammpstrj_to_txt.py

- 在 config_public.py 中配置环上的原子序号，以及根据分子对称性选取的环。
- 在main.py 中，修改文件路径，运行即可。
- 此外在 main.py 中可以配置聚类的数量，