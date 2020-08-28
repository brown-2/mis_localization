from main import *
from predict import *
from split_wrap import *
from scipy.linalg import expm
import os
from pathlib import Path

def gene_diff_record(cancer_name):
    dir_path = '../tem_data/{}_records/'.format(cancer_name)
    if not Path(dir_path).is_dir():
        Path(dir_path).mkdir()
    _, y = construct_and_return_loc_net()
    for i in tqdm(range(1, 21)):
        i /= 10
        nor_proba = np.load('../tem_data/{}_probas/probas_nor_tau_{}.npy'.format(cancer_name, i))
        ill_proba = np.load('../tem_data/{}_probas/probas_ill_tau_{}.npy'.format(cancer_name, i))
        nor_proba = proba_scaling(nor_proba)
        ill_proba = proba_scaling(ill_proba)

        diff = ill_proba - nor_proba
        np.savetxt(
                fname = dir_path + 'diff_{}.txt'.format(i),
                X = diff,
                fmt = '%.10f'
                )
        np.save(dir_path + 'diff_{}.npy'.format(i), diff)

        sorted_indices = diff.reshape(-1).argsort().tolist()
        sorted_indices.reverse()
        with open('../data/uni_ids_with_loc.txt') as f:
            proteins = f.read().split()
        with open('../data/location_name.txt') as f:
            locations = f.read().split('\n')#位置列表，12个
            locations = np.array(locations)
        with open(dir_path + '{}_record_{}.txt'.format(cancer_name, i), 'w') as f:
            id2gene = gene_mapping('id')
            id2entry = gene_mapping('id2entry')
            for i in sorted_indices:
                row, col = i // 12, i % 12
                protein_id = proteins[row]
                entry = id2entry[protein_id]
                gene_names = id2gene[protein_id]
                f.write(
                        '\t'.join([
                            protein_id, 
                            entry, 
                            #locations[col].ljust(21), 
                            '{}({:.3f}-{:.3f})'.format(locations[col], nor_proba[row][col], ill_proba[row][col]).ljust(40), 
                            str(round(diff[row, col],3)),
                            ' '.join(gene_names).ljust(30)
                            ]) + '\t'
                        )
                original = []
                for j in np.where(y[row])[0]:
                    location = locations[j]
                    nor_p = nor_proba[row][j]
                    ill_p = ill_proba[row][j]
                    original.append('{}({:.3f}-{:.3f})'.format(location, nor_p, ill_p))
                f.write('|'.join(original))
                f.write('\n')

if __name__ == '__main__':
    cancer_names = [
            'hepatitis',
            'leukemia',
            'breast'
            ]
    for cancer_name in cancer_names:
        gene_diff_record(cancer_name)
