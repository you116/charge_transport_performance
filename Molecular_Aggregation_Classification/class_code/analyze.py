"""
数据分析模块，使用 Kmeans 数据
"""
import json
import os
import time
from collections import Counter
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import numpy as np
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
from config_public import RING_ID_FOR_CLUSTER
from config_public import RING_ID_FOR_CLUSTER_T
class KmeansAnalyzer:
    """
    数据分析类, 使用 KMeans、DBSCAN、HDBSCAN 等聚类算法
    """

    def __init__(self, data_file_path=None, data_input=None, frame_num=None, result_save_dir=None, is_draw_gray_and_color=False):
        if not data_file_path and not data_input:
            raise ValueError('文件路径 和 数据列表 不能同时为空！')
        self.pending_data = self.build_data(data_file_path, data_input)  # 待分析的数据
        self.human_data = self.build_human_data(data_file_path, data_input)  # 人工标记的数据
        self.final_kmeans_data = None   # 最终的 KMeans 数据 7x7
        self.frame_num = frame_num
        self.result_save_dir = result_save_dir
        self.n_clusters = None

        self.kmeans_labels = None   # KMeans 聚类标签
        self.kmeans_centers = None  # KMeans 聚类中心点
        self.marked_data = None  # 标记的数据结果
        self.human_result_text = None  # 人工标记的结果文本
        self.is_draw_gray_and_color = is_draw_gray_and_color    # 是否生成灰度图和彩色图

    def build_data(self, data_file_path, data_input):
        if data_file_path:
            with open(data_file_path, 'r') as f:
                kmeans_data = json.load(f).get('kmeans_data')

        else:
            kmeans_data = data_input.get('kmeans_data')
        # print(f'数据构建完成！ 长度为 {len(kmeans_data)} ！， shape 为 {np.array(kmeans_data).shape} ！')
        return kmeans_data

    def zoom_out_kmeans_data_ywq(self):
        """
        根据分子对称性缩小矩阵。 在将分子上的环缩小为点之后，获得环与环之间的距离矩阵，需要根据对称性对其进行缩小
        举例：pending_data 为很多个 12x12 的矩阵，将每个矩阵缩小为 7x7 的矩阵
        每个 12x12 的矩阵表示一对分子 A 的12个环与另一对分子 B 的12个环的距离矩阵
        1，1 表示 A 的第一个环与 B 的第一个环的距离；1，2 表示 A 的第一个环与 B 的第二个环的距离 ...
        目的：分子是对称的，所以 12 个环是对称的，所以只需要 7 个环即可， 每个分子对取距离最小的 5 个环 加上中间的两个环。
        """
        final_data_list = []
        for origin_every_pair_data in self.pending_data:
            # every_pair_data 是一个 12x12 的矩阵
            # 提取出 4 个 7*7 的子矩阵。
            every_pair_data = np.array(origin_every_pair_data)
            ring_id_for_cluster = RING_ID_FOR_CLUSTER
            ring_id_for_cluster_t = RING_ID_FOR_CLUSTER_T
            ring_num_choose = len(ring_id_for_cluster)-1
            """
            s_a1_b1_1 = every_pair_data[:9, :9]   # A 的左半边（1-7）与 B 的左半边（1-7）的距离矩阵
            s_a1_b2_1 = every_pair_data[:9, 3:]  # A 的左半边（1-7）与 B 的右半边（6-12）的距离矩阵
            s_a2_b1_1 = every_pair_data[3:, :9]  # A 的右半边（6-12）与 B 的左半边（1-7）的距离矩阵
            s_a2_b2_1 = every_pair_data[3:, 3:]  # A 的右半边（6-12）与 B 的右半边（6-12）的距禝矩阵


            s_a1_b1 = s_a1_b1_1
            s_a1_b2 = s_a1_b2_1[:, np.arange(9)[::-1]]
            s_a2_b1 = s_a2_b1_1[np.arange(9)[::-1], :]
            s_a2_b2_2 = s_a2_b2_1[np.arange(9)[::-1], :]
            s_a2_b2 = s_a2_b2_2[:, np.arange(9)[::-1]]
            """
            s_a1_b1 = every_pair_data[np.ix_(ring_id_for_cluster, ring_id_for_cluster)]
            s_a1_b2 = every_pair_data[np.ix_(ring_id_for_cluster, ring_id_for_cluster_t)]
            s_a2_b1 = every_pair_data[np.ix_(ring_id_for_cluster_t, ring_id_for_cluster)]
            s_a2_b2 = every_pair_data[np.ix_(ring_id_for_cluster_t, ring_id_for_cluster_t)]
            print('s_', s_a2_b1)
            # 比较四个子矩阵，找出值最小的一个矩阵
            counts = [np.zeros((ring_num_choose, ring_num_choose)), np.zeros((ring_num_choose, ring_num_choose)), np.zeros((ring_num_choose, ring_num_choose)), np.zeros((ring_num_choose, ring_num_choose))]   # 计数器
            matrices = [s_a1_b1, s_a1_b2, s_a2_b1, s_a2_b2]
            for i in range(ring_num_choose):
                for j in range(ring_num_choose):
                    values = [m[i, j] for m in matrices]
                    min_index = values.index(min(values))
                    counts[min_index][i, j] += 1
            sums = [np.sum(c) for c in counts]
            max_index = sums.index(max(sums))   # 拥有小值最多的矩阵，即为最小的矩阵
            final_data_list.append(matrices[max_index])
        self.final_kmeans_data = final_data_list
        print(final_data_list)
        print(f'缩小矩阵完成！ shape 为 {np.array(final_data_list).shape} ！')

    def build_human_data(self, data_file_path, data_input):
        if data_file_path:
            with open(data_file_path, 'r') as f:
                human_data = json.load(f).get('human_data')
        else:
            human_data = data_input.get('human_data')
        # print(f'human 数据构建完成！ 长度为 {len(human_data)} ！')
        return human_data

    def calcul_kmeans(self):
        """
        计算 KMeans 聚类
        :param n_clusters: 聚类数量, 分为几类？
        :return:
        """
        # KMeans 算法需要一个二维输入，形状为 (n_samples, n_features)，所以需要将数据转换为二维，将每个 12x12 的矩阵展平为一个一维数组，然后将这些数组堆叠起来
        data_array = np.array(self.final_kmeans_data)
        data_2d = data_array.reshape((data_array.shape[0], -1))
        # print(f'Kmeans 计算，数据转换为二维完成！ shape 为 {data_2d.shape} ！')
        kmeans = KMeans(n_clusters=self.n_clusters)
        kmeans.fit(data_2d)
        labels = kmeans.labels_
        centers = kmeans.cluster_centers_
        self.kmeans_labels = labels
        self.kmeans_centers = centers
        return labels, centers

    def test_kmeans_n_clusters(self):
        """
        肘部法则  Elbow Method
        轮廓系数    Silhouette Coefficient
        测试Kmeans应该选择的聚类数量，n_clusters=?
        """
        wcss = []
        sil = []
        print('使用肘部法则和轮廓系数测试 KMeans 聚类数量...')
        for i in range(2, 11):
            data_array = np.array(self.final_kmeans_data)
            data_2d = data_array.reshape((data_array.shape[0], -1))
            kmeans = KMeans(n_clusters=i)
            kmeans.fit(data_2d)
            labels = kmeans.predict(data_2d)
            wcss.append(kmeans.inertia_)
            sil.append(silhouette_score(data_2d, labels))
        #
        # print('肘部法则测试结果：', wcss, '长度：', len(wcss))
        # print('轮廓系数测试结果：', sil, '长度：', len(sil))
        # 绘制 WCSS
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.plot(range(2, 11), wcss, marker='o')
        plt.title(f'Frame_{self.frame_num} Elbow Method')
        plt.xlabel(f'Frame_{self.frame_num} Number of clusters')
        plt.ylabel('WCSS')
        # 绘制轮廓系数
        # 绘制轮廓系数
        plt.subplot(1, 2, 2)
        plt.plot(range(2, 11), sil, marker='o')
        plt.title('Silhouette Coefficient')
        plt.xlabel('Number of clusters')
        plt.ylabel('Silhouette Coefficient')
        plt.tight_layout()
        # plt.ion()
        print(f'当前帧数：{self.frame_num}，请根据弹出的示意图，选择合适的 K 值，关闭图像后输入 K 值并回车 ！')
        plt.show()
        # plt.pause(300)
        n_clusters = input('请输入 K 值：')
        # 判断输入是否合法
        while not n_clusters.isdigit() or int(n_clusters) < 2 or int(n_clusters) > 12:
            n_clusters = input('输入不合法，请重新输入 K 值：')
        self.n_clusters = int(n_clusters)
        print(f'帧数：{self.frame_num} KMeans 聚类数量选择完成！ n_clusters = {self.n_clusters} ！')

    def mark_data(self):
        """
        标记数据
        :param labels: kmeans.labels_
        :param centers: kmeans.cluster_centers_
        :return:
        """
        final_data_list = []
        """
        TODO 原数据改变为 7*7 的形状, 搁置标记
        """
        i = 0
        for label in self.kmeans_labels:
            data = {
                'pair': self.human_data[i].get('pair'),
                'class': int(label),
            }
            final_data_list.append(data)
            i += 1
        self.marked_data = final_data_list
        return final_data_list

    def draw_kmeans_result(self):
        """
        绘制 KMeans 结果
        :param labels:
        :param centers:
        :return:
        """
        data_array = np.array(self.final_kmeans_data)
        data_reshape = data_array.reshape((data_array.shape[0], -1))
        tsne = TSNE(n_components=2, random_state=0)
        data_2d = tsne.fit_transform(data_reshape)
        # print(f'数据降维完成！ shape 为 {data_2d.shape} ！')
        fig, ax = plt.subplots(figsize=(10, 10))
        scatter = ax.scatter(data_2d[:, 0], data_2d[:, 1], c=self.kmeans_labels, cmap='viridis')

        # 添加图例
        legend1 = ax.legend(*scatter.legend_elements(),
                            loc="upper right", title="Clusters")
        ax.add_artist(legend1)
        #ax.set_xticks([])
        #ax.set_yticks([])
        plt.savefig(os.path.join(self.result_save_dir, f'{self.frame_num}_KMeans结果图.png'))
        # plt.show()

    def analyze_kmeans_result(self):
        """
        分析 KMeans 结果
        :return:
        """
        # 将聚类中心点还原为 7*7 的矩阵， 共有 n_clusters 个聚类中心点

        centers = self.kmeans_centers.reshape((self.kmeans_centers.shape[0], len(RING_ID_FOR_CLUSTER), len(RING_ID_FOR_CLUSTER)))

        re_text_list = []

        counts = Counter(self.kmeans_labels)
        for i, center in enumerate(centers):
            # print(f'第 {i} 类的聚类中心点为：\n{center} ！')
            re_text_list.append(f'第 {i} 类的数量为：{counts[i]} ！')

        self.human_result_text = re_text_list

    def save_result(self):
        """
        保存结果
        :return:
        """
        result_file_path = os.path.join(self.result_save_dir, f'{self.frame_num}_KMeans结果.json')
        write_data = {
            'marked_data': self.marked_data,
            'human_result_text': self.human_result_text,
            'kmeans_labels': str(self.kmeans_labels),
            'kmeans_centers': str(self.kmeans_centers),
        }
        print(f'保存 KMeans 结果到文件 {result_file_path} 中！\n')
        with open(result_file_path, 'w') as f:
            json.dump(write_data, f, indent=4, ensure_ascii=False)

    """
    24-03-13 ---- 在 kmeans 之前对数据进一步处理 from WQ：
    1.取倒数，2.归一化处理（我搜的是用z-score）3.将归一化后的数据处理成一定范围的图像，不难，有现成的调用
    
    取倒数和归一化 ： reciprocal_normalized 函数。
    画出灰度图： draw_grayscale_image 函数。
    画出彩色图： draw_color_image 函数。
    画出灰度图和彩色图： draw_grayscale_and_color_image 函数。
            ----  03-13 end ----
    """
    def reciprocal_normalized(self):
        """
        取倒数, 归一化处理, 将归一化后的数据处理成灰度图
        """
        # print('取倒数')
        # print('self.final_kmeans_data: 类型：', type(self.final_kmeans_data), "长度：", len(self.final_kmeans_data))

        op_data = np.array(self.final_kmeans_data)
        if (op_data == 0).any():
            raise ValueError('!!!  注意： 数据中包含 0 值，取倒数时将无穷大。 这代表数据可能有错误，请检查！')
        reciprocal_data = 1.0 / op_data
        # print('取倒数之后的数据: 类型：', type(reciprocal_data), "长度：", len(reciprocal_data), 'shape: ', reciprocal_data.shape)

        # 归一化处理
        # print('归一化处理')
        min_value = np.min(reciprocal_data)
        max_value = np.max(reciprocal_data)
        normalized_data = (reciprocal_data - min_value) / (max_value - min_value)
        # print('归一化之后的数据: 类型：', type(normalized_data), "长度：", len(normalized_data), 'shape: ', normalized_data.shape)
        self.final_kmeans_data = normalized_data

    def draw_grayscale_image(self, save_dir):
        """
        画出灰度图, 不显示颜色条
        """
        print('正在生成灰度图...')
        # i = 0
        # for data in self.final_kmeans_data:
        #     plt.imshow(data, cmap='gray_r')
        #     plt.axis('off')
        #     plt.savefig(os.path.join(save_dir, f'{self.frame_num}_{i}_灰度图.png'))
        #     i += 1
        #     plt.close()
        for i in range(len(self.final_kmeans_data)):
            pair = self.human_data[i].get('pair')
            pair_text = 'PE' + str(pair)    # PE PairElement, 成对的元素
            class_text = 'C' + str(self.kmeans_labels[i])   # C Class
            img_file_name = f'F{self.frame_num}_{class_text}_PN{i}_{pair_text}_灰度.png'
            plt.imshow(self.final_kmeans_data[i], cmap='gray_r')
            plt.axis('off')
            plt.savefig(os.path.join(save_dir, img_file_name))
            plt.close()

    def draw_color_image(self, save_dir):
        """
        画出彩色图, 保留颜色条
        """
        print('正在生成彩色图...')


        for i in range(len(self.final_kmeans_data)):
            pair = self.human_data[i].get('pair')
            pair_text = 'PE' + str(pair)    # PE PairElement, 成对的元素
            class_text = 'C' + str(self.kmeans_labels[i])   # C Class
            img_file_name = f'F{self.frame_num}_{class_text}_PN{i}_{pair_text}_彩色.png'
            plt.imshow(self.final_kmeans_data[i], cmap='jet')
            plt.axis('off')
            plt.colorbar()
            plt.savefig(os.path.join(save_dir, img_file_name))
            plt.close()

    def draw_grayscale_and_color_image(self):
        """
        画出灰度图和彩色图
        """
        if self.is_draw_gray_and_color:
            gray_save_dir = os.path.join(self.result_save_dir, 'gray_image')
            color_save_dir = os.path.join(self.result_save_dir, 'color_image')
            if not os.path.exists(gray_save_dir):
                os.makedirs(gray_save_dir)
            if not os.path.exists(color_save_dir):
                os.makedirs(color_save_dir)
            self.draw_grayscale_image(gray_save_dir)
            self.draw_color_image(color_save_dir)
            print(f'灰度图和彩色图生成完成！保存在 {gray_save_dir} 和 {color_save_dir} 中！')

    def exec(self):
        self.zoom_out_kmeans_data_ywq()
        self.reciprocal_normalized()    # 取倒数，归一化处理

        self.test_kmeans_n_clusters()  # 测试 KMeans 聚类数量

        self.calcul_kmeans()
        self.draw_grayscale_and_color_image()  # 画出灰度图
        self.mark_data()
        self.analyze_kmeans_result()
        self.save_result()
        self.draw_kmeans_result()
        print(f'KMeans 分析完成！')
        # print('human_data_len', len(self.human_data))
        # print('pending_data_len', len(self.pending_data))
        # print('final_kmeans_data_len', len(self.final_kmeans_data))

    @classmethod
    def merge_kmeans(cls, data_list, save_result_dir, n_clusters=5):
        """
        合并数据做一次大的 kmeans，用来画图, 仅生成图片
        :param n_clusters:
        :param data_list:
        :param save_result_dir:
        :return:
        """
        print(f'正在合并 {len(data_list)} 个数据，进行 KMeans 分析... ！')
        data_array = np.concatenate(data_list, axis=0)
        print(f'合并完成！ shape 为 {data_array.shape} ！')
        data_2d = data_array.reshape((data_array.shape[0], -1))
        kmeans = KMeans(n_clusters=n_clusters)
        kmeans.fit(data_2d)
        labels = kmeans.labels_

        tsne = TSNE(n_components=2, random_state=0)
        draw_data_2d = tsne.fit_transform(data_2d)
        fig, ax = plt.subplots(figsize=(10, 10))
        scatter = ax.scatter(draw_data_2d[:, 0], draw_data_2d[:, 1], c=labels, cmap='viridis')
        #legend1 = ax.legend(*scatter.legend_elements(),
        #                    loc="upper right", title="Clusters")
        #ax.add_artist(legend1)
        ax.set_xticks([])
        ax.set_yticks([])
        file_path = os.path.join(save_result_dir, f'000合并的KMeans结果图{time.time()}.png')

        plt.savefig(file_path)
        print('合并的 KMeans 结果图已保存到文件中！', file_path)


if __name__ == '__main__':
    Data_File_Path = r'D:\data_calculate_ZJY\N4\process\2_KMeans数据.json'
    analyzer = KmeansAnalyzer(data_file_path=Data_File_Path, frame_num='2',
                              result_save_dir=r'D:\data_calculate_ZJY\N4\result',
                              is_draw_gray_and_color=True)
    # print('human_data_len', len(analyzer.human_data))
    # print('pending_data_len', len(analyzer.pending_data))
    analyzer.exec()
    # print('final_kmeans_data_len', len(analyzer.final_kmeans_data))
    # test_data = [analyzer.pending_data, analyzer.pending_data, analyzer.pending_data]
    # KmeansAnalyzer.merge_kmeans(test_data, r'D:\Code\wq-data-statistic\Molecular_Aggregation_Classification\test_save_dir')