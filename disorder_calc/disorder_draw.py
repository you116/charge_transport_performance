import json
import matplotlib.pyplot as plt
import  numpy as np
import scipy.stats as stats

"""
从get_total_disorder 中获取json文件，然后提取json文件，画无序的图
"""
import matplotlib.pyplot as plt

#提取文件


total_file = []
e_total_value_list = []
i_total_value_list = []

class disorder_draw:
    def __init__(self):
        self.e_total_value_list = []
        self.i_total_value_list = []
        self.x_ip = None
        self.y_ip = None
        self.x_ea = None
        self.y_ea = None

    def get_value(self,filename):
        with open(filename) as file:
            data = json.load(file)
            for value in data:
                self.e_total_value_list.append(value['e_total_value'])
                self.i_total_value_list.append(value['i_total_value'])

    def plot_scf_value_mean(self):

        data_ea = abs(np.array(self.e_total_value_list))
        mean_ea = np.mean(data_ea)
        std_dev_ea = np.std(data_ea)

        data_ip = abs(np.array(self.i_total_value_list))
        mean_ip = np.mean(data_ip)
        std_dev_ip = np.std(data_ip)

        self.x_ip = np.linspace(mean_ip - 3 * std_dev_ip, mean_ip + 3 * std_dev_ip, 100)
        self.y_ip = stats.norm.pdf(self.x_ip, mean_ip, std_dev_ip)

        print('mean_ea',mean_ea)
        print('mean_ip', mean_ip)
        print('std_dev_ip',std_dev_ip)
        print('std_dev_ea', std_dev_ea)
        self.x_ea = np.linspace(mean_ea - 3 * std_dev_ea, mean_ea + 3 * std_dev_ea, 100)
        self.y_ea = stats.norm.pdf(self.x_ea, mean_ea, std_dev_ea)


'''
for filename in ['12', '30', '80', '100', '120', '140', '160', '180']:
    filename2 = 'D:\shujufenxi\Y6\disorder\y6\static\\'+ filename+'.json'
    t1.get_value(filename2)
'''

t1 = disorder_draw()

t1.get_value('D:\shujufenxi\Y6\disorder\y63\y611120.json')

t1.plot_scf_value_mean()
print(len(t1.e_total_value_list))
print(len(t1.i_total_value_list))

plt.plot(t1.x_ip,t1.y_ip)
plt.show()





