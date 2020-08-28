from main import *
def locative_loc_net(l, loc):
    '''
    输入两个对象, 蛋白列表和定位矩阵. l, loc
    返回三个对象, 新的蛋白名称列表, 新的定位矩阵和被拆分蛋白的字典l_new, a_new, d
    '''
    loc_new = np.zeros((loc.sum().astype(int), loc.shape[1]))
    n_loc = loc.sum(1).astype(int)
    d = dict()
    l_new = []
    i_new = 0
    for i in range(len(l)):
        if n_loc[i] == 1:
            l_new.append(l[i])
            loc_new[i_new] = loc[i]
            i_new += 1
        else:
            d[l[i]] = []#记录被拆分的蛋白
            nonzero_index = loc[i].nonzero()[0]
            for j in range(n_loc[i]):
                l_new.append(l[i] + '_{}'.format(j))
                d[l[i]].append(l[i] + '_{}'.format(j))
                loc_new[i_new, nonzero_index[j]] = 1
                i_new += 1
    return l_new, loc_new, d

def locative_laplacian(l, a, d):
    '''
    输入三个对象, 蛋白列表, 系数矩阵和待拆分蛋白的字典
    返回两个对象, 新的蛋白列表和系数矩阵
    '''
    l_new = []
    for i in l:
        if i not in d:
            l_new.append(i)
        else:
            for j in d[i]:
                l_new.append(j)
    mapping = dict(zip(l_new, range(len(l_new))))
    a_new = np.zeros((len(l_new), len(l_new)))
    for x,y in zip(*a.nonzero()):
        p1, p2 = l[x], l[y]
        indices_x, indices_y = [], []
        if p1 in d:
            for p in d[p1]:
                indices_x.append(mapping[p])
        else:
            indices_x.append(mapping[p1])
        if p2 in d:
            for p in d[p2]:
                indices_y.append(mapping[p])
        else:
            indices_y.append(mapping[p2])
        for ix in indices_x:
            for iy in indices_y:
                a_new[ix, iy] = a[x, y]
    return l_new, a_new
def gene_real2virtual_mapping_dict():
    '''
    辅助函数，用来生成映射字典，方便转换索引
    '''
    with open('../data/uni_ids_with_loc.txt') as f:
        l_real = f.read().split()
    with open('../data/locative_uni_ids_with_loc.txt') as f:
        l_virtual = f.read().split()
    with open('../data/split_protein_dict.json') as f:
        d = json.load(f)
    mapping = dict(zip(l_virtual, range(len(l_virtual))))
    d_new = {}
    for i in range(len(l_real)):
        p = l_real[i]
        if p not in d:
            d_new[i] = mapping[p]
        else:
            d_new[i] = []
            for j in d[p]:
                d_new[i].append(mapping[j])
    return d_new
    
def real2virtual_map(indices_real):
    '''
    映射函数，将一组实体蛋白的索引映射成相应的locative蛋白索引
    '''
    d = gene_real2virtual_mapping_dict()
    indices_virtual = []
    for index in indices_real:
        t = d[index]
        if isinstance(t, int):
            indices_virtual.append(t)
        else:
            indices_virtual.extend(t)
    return np.array(indices_virtual)
def real2virtual_map_unique(indices_real):
    d = gene_real2virtual_mapping_dict()
    indices_virtual = []
    for index in indices_real:
        t = d[index]
        if isinstance(t, int):
            indices_virtual.append(t)
        else:
            indices_virtual.append(t[0])
    return np.array(indices_virtual)

def unittest():
    l = [0,1,2,3]
    a = np.array([
        [0,1,0,0],
        [1,0,1,1],
        [0,1,0,0],
        [0,1,0,0]])
    d = {1:[4,5], 2:[6,7]}
    l_new, a_new = locative_laplacian(l, a, d)
    assert l_new == [0,4,5,6,7,3]
    a_test = np.array([ 
        [0,1,1,0,0,0],
        [1,0,0,1,1,1],
        [1,0,0,1,1,1],
        [0,1,1,0,0,0],
        [0,1,1,0,0,0],
        [0,1,1,0,0,0],
        ])
    assert (a_new == a_test).all()






