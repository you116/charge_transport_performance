---
pair_extract文件
---
1. 需要获取分子对----计算分子上原子间的距离矩阵----pair_extract/face_onextract_ouhe.py
2. 获取用于耦合值计算的输出文件
3. 获取用于EA值计算的输出文件

---
molecular_aggregation_classification文件
---
用于聚类计算，需要输入分子对信息

---
net_analysis
--
用于计算电连接网络，需要获取分子对信息，输出用于耦合值计算的高斯文件，批量计算高斯文件，网络计算

---
charge_transport_rate
---
电荷转移速率计算，需要知道分子对信息，重组能大小，每个分子对的EA值等

---
disorder_calc
---
用于计算IP和EA的总无序、动态无序和静态无序，需要获取分子对信息，和EA值

---
carrier_mobility
---
蒙特卡洛方法计算载流子迁移率，需要知道分子对间的电荷转移速率信息。

---
流程计算文件
---
批量执行高斯计算和电荷耦合值







