from main import *
from predict import *
from split_wrap import *
from scipy.linalg import expm
import os


def gene_txt_results(cancer_name):
    in_path = '../tem_data/{}_probas/'.format(cancer_name)

    id_list = get_annotated_ids()
    id_list = np.array(id_list).reshape(-1,1)

    for i in range(1, 21):
        i /= 10
        nor_proba = np.load('../tem_data/{}_probas/probas_nor_tau_{}.npy'.format(cancer_name, i))
        ill_proba = np.load('../tem_data/{}_probas/probas_ill_tau_{}.npy'.format(cancer_name, i))
        nor_proba = proba_scaling(nor_proba)
        ill_proba = proba_scaling(ill_proba)
        nor_out_path = '../tem_data/probas_txt/{}_nor_probas_{}.txt'.format(cancer_name, i)
        ill_out_path = '../tem_data/probas_txt/{}_ill_probas_{}.txt'.format(cancer_name, i)

        np.savetxt(
                fname = nor_out_path,
                X = np.append(id_list, np.round(nor_proba, 3), 1),
                fmt = '%5s'
                )
        np.savetxt(
                fname = ill_out_path,
                X = np.append(id_list, np.round(ill_proba, 3), 1),
                fmt = '%5s'
                )



if __name__ == '__main__':
    cancer_names = [
            'hepatitis',
            'leukemia',
            'breast'
            ]
    for cancer_name in cancer_names:
        gene_txt_results(cancer_name)
