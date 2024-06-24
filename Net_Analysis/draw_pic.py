import json
import numpy as np
from matplotlib import font_manager
import matplotlib.pyplot as plt


x_list = np.arange(0.0, 0.052, 0.002)
# print('X_list: ', x_list)


def read_draw_data():
    with open('0423_draw_data_lum1.json', 'r') as lum1_file:
        lum1_data = json.load(lum1_file)

    with open('0423_draw_data_LUMO0.json', 'r') as lum0_file:
        lum0_data = json.load(lum0_file)

    return lum0_data, lum1_data


def draw_pic(x, y):
    # 设置全局字体为新罗马
    plt.rcParams['font.family'] = 'Times New Roman'

    #plt.figure(figsize=(6, 6)) # 设置画布大小
    plt.plot(x, y, 'bo', label='PDIN')

    # plt.ylim((0, 1))
    plt.xlim((0, 0.052))
    plt.ylabel('', fontsize=24, labelpad=10)  # y轴标签, 20号字体，标签距离坐标轴20个像素
    plt.xlabel('V(eV)', fontsize=24, labelpad=10)  # x轴标签, 20号字体，标签距离坐标轴20个像素
    plt.xticks(fontsize=20)
    # plt.xticks(fontsize=25, rotation=45)    # x轴标签旋转45度
    # 自定义 x 轴刻度显示为两位小数
    plt.gca().xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    # 移除 x 轴第一个刻度值
    plt.gca().xaxis.get_major_ticks()[0].label1.set_visible(False)
    plt.yticks(fontsize=20)
    plt.legend(prop={'size': 10},)
  #  plt.tight_layout()  # 自动调整图像边缘, 使图像完全显示, 不超出边界
    # 设置 x 轴刻度数量为五个
    #plt.locator_params(axis='x', nbins=6)
    plt.show()


def draw_pic_2(x_list, y_list):
    # 设置全局字体为新罗马
    plt.rcParams['font.family'] = 'Times New Roman'

    colors = ['brown', 'green', 'blue']  # 线的颜色列表
    #line_styles = ['o-', 's-', '^-',]  # 线的样式列表 'o-' 圆点  , '-'普通线 , '^-'三角线 , 's-' 方块线,
    #name = ['Y6', 'N4', 'N3']
    #plt.figure(figsize=(7, 6))
    # 绘制多条线
    for i in range(len(x_list)):
        #平滑曲线
        # window_size = 5
        # y_smooth = np.convolve(y_list[i], np.ones(window_size) / window_size, mode='same')

        plt.plot(x_list[i], y_list[i], 'o-', color=colors[i], label='Line {}'.format(i+1))

    plt.ylim((0, 200))
    plt.ylabel('Number of clusters', fontsize=22, labelpad=10)
    plt.xlabel('V(eV)', fontsize=22, labelpad=10)
    plt.xticks(fontsize=20)
    plt.gca().xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    plt.gca().xaxis.get_major_ticks()[0].label1.set_visible(False)
    plt.yticks(fontsize=20)
    plt.legend(['Y6', 'N4', 'N3'], prop={'size': 15})
    plt.tight_layout()
    plt.savefig('wq_pic_2.png')
    plt.show()


if __name__ == '__main__':
    d0, d1 = read_draw_data()
    print('d0 keys', d0.keys())
    print('d1 keys', d1.keys())
    num_y6_lum0 = d0.get('num_y6')
    num_n3_lum0 = d0.get('num_n3')
    num_n4_lum0 = d0.get('num_n4')
    max_y6 = d1.get('max_y6')
    max_n3 = d1.get('max_n3')
    max_n4 = d1.get('max_n4')
    average_y6 = d0.get('average_y6')
    average_n3 = d0.get('average_n3')
    average_n4 = d0.get('average_n4')

    draw_pic_2([x_list, x_list, x_list], [num_y6_lum0, num_n4_lum0, num_n3_lum0])
    #draw_pic_2([x_list, x_list, x_list], [max_y6, max_n4, max_n3])
    #draw_pic_2([x_list, x_list, x_list], [average_y6, average_n4, average_n3])

