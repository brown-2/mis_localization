from main import *
from predict import *
from split_wrap import *
from scipy.linalg import expm
import os
from pathlib import Path

def gene_diff_record(cancer_name):
    dir_path = '../tem_data/{}_records/'.format(cancer_name)
    print(cancer_name)
    for i in range(1, 21):
        i /= 10
        diff = np.load(dir_path + 'diff_{}.npy'.format(i))
        std = diff.std()
        mean = diff.mean()
        print(cancer_name, i, mean, std)


if __name__ == '__main__':
    cancer_names = [
            'hepatitis',
            'leukemia',
            'breast'
            ]
    for cancer_name in cancer_names:
        gene_diff_record(cancer_name)
