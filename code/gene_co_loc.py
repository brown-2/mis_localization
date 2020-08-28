from main import *
from pathlib import Path
def gene_micro_co_location():
    l, loc = construct_and_return_loc_net()
    link = np.zeros((len(l), len(l)))
    for i in tqdm(range(len(l))):
        for j in range(i + 1, len(l)):
            link[i, j] = link[j, i] = \
                    np.logical_and(loc[i], loc[j]).sum()/\
                    np.logical_or(loc[i], loc[j]).sum()
    return link
def gene_co_location():
    micro = np.load('../tem_data/co_loc_micro.npy')
    a = np.load('../tem_data/a.npy')
    l = get_uni_list()
    l_loc, _ = construct_and_return_loc_net()
    mapping = dict(zip(l, range(len(l))))
    res = np.zeros(a.shape)
    for i in tqdm(range(len(micro))):
        for j in range(i + 1, len(micro)):
            x, y = l_loc[i], l_loc[j]
            ix, iy = mapping[x], mapping[y]
            res[ix, iy] =res[iy, ix] = micro[i, j]
    return res
def gene_adjacency_between_labels():
    l, a = construct_and_return_uni_net()
    l_loc, loc = construct_and_return_loc_net()
    mapping = dict(zip(l, range(len(l))))
    res = np.zeros(a.shape)
    for i in tqdm(range(len(l_loc))):
        for j in range(i + 1, len(l_loc)):
            x, y = l_loc[i], l_loc[j]
            ix, iy = mapping[x], mapping[y]
            res[ix, iy] =res[iy, ix] = 1
    return res



path = Path('../tem_data/co_loc_micro.npy')
if not path.is_file():
    co_loc_micro = gene_micro_co_location()
    np.save(path, co_loc_micro)
path2 = Path('../tem_data/co_loc.npy')
if not path2.is_file():
    co_loc = gene_co_location()
    np.save(path2, co_loc)

path3 = Path('../tem_data/a_between_labels.npy')
if not path3.is_file():
    a_between_labels = gene_adjacency_between_labels()
    np.save(path3, a_between_labels)

