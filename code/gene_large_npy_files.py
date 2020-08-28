from main import *
from pathlib import Path

def gene_expr(file_name, pcc_process_mode = 'abs'):
    '''
    返回表达数据，16319行，不同列数
    '''
    if 'mas5.csv' not in file_name.split('_')[-1]:
            print('file name format error! Please use xx_normal/ill_xxmas5.csv')
            return

    ncol = file_name.split('_')[-1][:-8]
    ncol = int(ncol)
    l, expr = construct_and_return_expr_matrix(file_name = file_name, ncol = ncol)
    return expr
def modify_a(a, expr_nor, expr_ill):
    '''
    输入邻接矩阵和表达数据，返回修改过的邻接矩阵
    '''
    #计算原始的pcc并将主对角线元素归零
    pcc_nor = np.corrcoef(expr_nor)
    pcc_ill = np.corrcoef(expr_ill)
    np.fill_diagonal(pcc_nor, 0)
    np.fill_diagonal(pcc_ill, 0)
    nor_non = np.isnan(pcc_nor)
    ill_non = np.isnan(pcc_ill)
    print('nor_non:', nor_non.sum())
    print('ill_non:', ill_non.sum())
    pcc_nor[np.isnan(pcc_nor)] = 1
    pcc_ill[np.isnan(pcc_ill)] = 1

    # 得到差分矩阵
    diff = pcc_ill - pcc_nor
    new_a = a.copy()
    diff_on_edge = diff.copy()
    diff_on_edge[a==0] = 0
    print('diff_on_edge number:', np.count_nonzero(diff_on_edge))

    diff_off_edge = diff.copy()
    diff_off_edge[a==1] = 0
    print('diff_off_edge number:', np.count_nonzero(diff_off_edge))

    # 得到标准差
    on_edge_std = diff[np.logical_and(a==1, np.logical_not(nor_non))].std()
    on_edge_mean = diff[np.logical_and(a==1, np.logical_not(nor_non))].mean()
    off_edge_std = diff[np.logical_and(a==0, np.logical_not(nor_non))].std()
    off_edge_mean = diff[np.logical_and(a==0, np.logical_not(nor_non))].mean()
    print('on_edge_mean', on_edge_mean)
    print('off_edge_mean', off_edge_mean)
    print('on_edge_std', on_edge_std)
    print('off_edge_std', off_edge_std)

    #修改拓扑
    new_a[diff_on_edge < on_edge_mean - on_edge_std * 3] = 0
    print(new_a.sum())
    new_a[diff_off_edge > off_edge_mean + off_edge_std * 3] = 1
    print(new_a.sum())


    return new_a

if __name__ == '__main__':
    _, a = construct_and_return_uni_net()
    ecc = ECC(a)

    np.save('../tem_data/a.npy', a)
    np.save('../tem_data/ECC.npy', ecc)

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
        print(normal_name.split('_')[0])
        expr_nor = gene_expr(normal_name)
        expr_ill = gene_expr(ill_name)

        cancer_name = normal_name.split('_')[0]

        new_a = modify_a(a, expr_nor, expr_ill)
        np.save('../tem_data/a_{}.npy'.format(cancer_name), new_a)

        cancer_name = normal_name.split('_')[0]
        new_a = np.load('../tem_data/a_{}.npy'.format(cancer_name))
        new_ECC = ECC(new_a)
        print('ECC:',new_ECC.sum())
        np.save('../tem_data/ECC_{}.npy'.format(cancer_name), new_ECC)

