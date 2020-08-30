import os
import numpy as np
from main import *

import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'

from predict import *

def plot_carlibrating_curve():
    def gene_data(cancer_name = 'hepatitis', threshold = 0.1):
        '''
        处理数据，产生数值矩阵
        '''
        y = np.loadtxt('../tem_data/probas_txt/label_mat.txt')
        perfs_m = np.zeros((20, 10))
        for i in range(1, 21):
            i /= 10
            state = 'nor'
            proba = np.load('../tem_data/{}_probas/probas_{}_tau_{}.npy'.format(cancer_name, state, i))
            proba = proba_scaling(proba, 'lin')
            pred = calibrate(proba, threshold)
            perfs = multi_label_performance(pred, y)
            j = int(i * 10 - 1)
            perfs_m[j] = perfs
        return perfs_m

    def wrap_data(data):
        '''
        将数值矩阵转化为namedtuple类型，方便后续处理
        '''
        from collections import namedtuple
        Perf = namedtuple('P', ['PPV','SE','SP','ACC2', 'MCC', 'AIM', 'CVR', 'ACC1', 'ATR', 'AFR'])
        data_tuple = Perf(*data.T)
        return data_tuple

    def plot_perf(data_tuple, perf_type):
        font = {'family' : 'Times New Roman', 'weight' : 'normal', }
        fig, ax = plt.subplots()
        x = np.arange(0.1, 2.1, 0.1)
        ax.set_xlabel('τ')
        ax.set_ylabel('MCC score')
        ax.set_xticks(x)
        ax.grid(True)#添加网格
        for i in range(len(data_tuples)):
            data_tuple = data_tuples[i]
            y = getattr(data_tuple, perf_type)
            ax.plot(x, y, label = '{} = {}'.format(r'$ \alpha $',(i + 1)/10))
            ax.legend()
        file_name = os.path.join(fig_path, 'C_MCC.{}'.format(form))
        plt.savefig(file_name, format = form, dpi = n_dpi)

    data_tuples = []
    for threshold in (x/10 for x in range(1, 4)):
        data = gene_data(cancer_name = 'hepatitis', threshold = threshold)
        data_tuple = wrap_data(data)
        data_tuples.append(data_tuple)
    plot_perf(data_tuples, 'MCC')

def plot_location_distribution():
    with open('../data/location_name.txt') as f:
        x = f.read().split('\n')[:-1]
    #y = [31, 2274, 97, 524, 526, 600, 78, 761, 121, 366, 2864, 1326]
    _, loc = construct_and_return_loc_net()
    y = loc.sum(0).astype(int)
    x.reverse()
    fig, ax = plt.subplots()
    b = ax.barh(x, y)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]

    ax.set_xlim(0, 3500)

    for rect in b:
        w = rect.get_width()
        ax.text(w, rect.get_y()+rect.get_height()/2, '%d' % int(w))

    ax.set_xlabel('number of proteins')
    plt.tight_layout()
    file_name = os.path.join(fig_path, 'location_distribution.{}'.format(form))
    plt.savefig(file_name, format = form, dpi = n_dpi)

def plot_multiplicity():
    #y = [4118, 1736, 504, 98, 15, 2]
    _, loc = construct_and_return_loc_net()
    y = np.bincount(loc.sum(1).astype(int))[1:]
    x = range(1, 7)
    fig, ax = plt.subplots()
    b = ax.bar(x, y)
    for a,b in zip(x,y):  
         ax.text(a, b+0.05, '%.0f' % b, ha='center', va= 'bottom',fontsize=11)  

    ax.set_xlabel('subcellular location multiplicity')
    ax.set_ylabel('number of proteins')
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    file_name = os.path.join(fig_path, 'multiplicity.{}'.format(form))
    plt.savefig(file_name, format = form, dpi = n_dpi)


if __name__ == '__main__':
    fig_path = '../figures'
    form = 'svg'
    form = 'tif'
    form = 'eps'
    n_dpi = None
    plot_carlibrating_curve()
    plot_location_distribution()
    plot_multiplicity()
