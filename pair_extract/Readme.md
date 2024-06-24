---
计算耦合值
---
1. 调用 t1.calc_mol_pair(output_filename)
   output_filename 为输出文件位置
2. 调用 t1.get_all_pair_mol(distance_cord, min_close_contact_num)
    distance_cord 为步骤1的输出文件，min_close_contact_num 为紧密接触的数量，提取 face_on 分子对时，分子上五元环为25，总分子对为1个或者8个

---
计算EA/IP值
---
IP = Ecation - Eneutral
EA = Eneutral - Eanion
1. 调用 t1.get_total_mol('_oni.gjf', charge_str = '-1 2') 获得 Eanion
2. 调用 t1.get_total_mol('_cat.gjf', charge_str = '1 2') 获得 Ecation
