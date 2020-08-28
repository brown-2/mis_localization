from main import *
from predict import *
from split_wrap import *
from scipy.linalg import expm
import os


def gene_PCC_file(file_name, pcc_process_mode = 'abs'):
    if 'mas5.csv' not in file_name.split('_')[-1]:
            print('file name format error! Please use xx_normal/ill_xxmas5.csv')
            return

    ncol = file_name.split('_')[-1][:-8]
    ncol = int(ncol)
    l, a = construct_and_return_expr_matrix(file_name = file_name, ncol = ncol)
    a = np.corrcoef(a)
    np.fill_diagonal(a, 0)
    a[np.isnan(a)] = 1
    if pcc_process_mode == 'abs':
        a = np.abs(a)
    elif pcc_process_mode == 'lin_scale':
        a += 1
        a /= a.max()
    elif pcc_process_mode == 'raw':
        pass
    return a

def alternative_proba(normal_name, ill_name, dir_suffix = 'probas', pcc_process_mode = 'abs'):
    cancer_name = normal_name.split('_')[0]

    params = {
            'kernel':'precomputed',
            'C':10,
            'cache_size':4000,
            'class_weight':'balanced',
            'probability':True,
            'decision_function_shape':'ovo',
            'random_state': 100,
            }
    l_loc, y = construct_and_return_loc_net()
    l_loc_new, y_new, d = locative_loc_net(l_loc, y)
    with open('../data/uni_ids.txt') as f:
        l_ppi = f.read().split()

    #加权重
    weight = y.sum(0).max()/y.sum(0)
    weight **= 2
    params['class_weight'] = dict(zip(range(12), weight))

    from pathlib import Path
    dir_path = '../tem_data/{}_{}/'.format(cancer_name, dir_suffix)
    if not Path(dir_path).is_dir():
        Path(dir_path).mkdir()

    #e_nor_01_path = Path('../tem_data/{}_probas/e_nor_0.1.npy'.format(cancer_name))
    #e_ill_01_path = Path('../tem_data/{}_probas/e_ill_0.1.npy'.format(cancer_name))
    e_nor_01_path = Path(dir_path + 'e_nor_0.1.npy')
    e_ill_01_path = Path(dir_path + 'e_ill_0.1.npy')
    if e_nor_01_path.exists():
        # 如果已经有，直接加载
        e_nor_01 = np.load(str(e_nor_01_path))
        e_ill_01 = np.load(str(e_ill_01_path))
    else:
        #如果没有，计算出来，并保存方便以后使用
        PCC_normal = gene_PCC_file(normal_name)
        PCC_ill = gene_PCC_file(ill_name)

        #ECC
        ecc = np.load('../tem_data/ECC.npy')
        ecc_ill = np.load('../tem_data/ECC_{}.npy'.format(cancer_name))
        t_nor = PCC_normal * ecc
        t_ill = PCC_ill * ecc_ill
        l_new, t_nor_new = locative_laplacian(l_ppi, t_nor, d)
        l_new, t_ill_new = locative_laplacian(l_ppi, t_ill, d)
        np.fill_diagonal(t_nor_new, -t_nor_new.sum(0))
        np.fill_diagonal(t_ill_new, -t_ill_new.sum(0))

        e_nor_01 = expm(t_nor_new/10)
        e_ill_01 = expm(t_ill_new/10)
        np.save(str(e_nor_01_path), e_nor_01)
        np.save(str(e_ill_01_path), e_ill_01)

    e_nor_curr = np.eye(len(e_nor_01))
    e_ill_curr = np.eye(len(e_nor_01))
    for i in tqdm(range(1, 21)):
        i /= 10
        # 检查是否已存在

        e_nor_curr = e_nor_curr @ e_nor_01
        e_ill_curr = e_ill_curr @ e_ill_01

        X = get_K_from_e(e_nor_curr)
        X_ill = get_K_from_e(e_ill_curr)
        proba_nor, proba_ill = KFold_predict_locative(X, X_ill, y_new, y, params)
        np.save(
                file = dir_path + 'probas_nor_tau_{}.npy'.format(i),
                arr = proba_nor,
                )
        np.save(
                file = dir_path + 'probas_ill_tau_{}.npy'.format(i),
                arr = proba_ill,
                )
        print('probas of tau {} has been saved.'.format(i))

if __name__ == '__main__':
    normal_names = [
            'hepatitis_normal_6mas5.csv',
            'leukemia_bone_marrow_normal_10mas5.csv',
            'breast_normal_20mas5.csv'
            ]
    ill_names = [
            'hepatitis_ill_21mas5.csv',
            'leukemia_bone_marrow_ill_7mas5.csv',
            'breast_ill_17mas5.csv'
            ]
    for normal_name, ill_name in zip(normal_names, ill_names):
        alternative_proba(normal_name, ill_name)
