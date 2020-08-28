import os 
import time
import json
import gzip
import numpy as np
from tqdm import tqdm
from sklearn.svm import SVC
from sklearn.multiclass import OneVsRestClassifier
from pathlib import Path

def get_annotated_ids():
    with open('../data/uni_ids_with_loc.txt') as f:
        return f.read().split()
def construct_and_return_bio_net():
    ppi_file = '../data/BIOGRID-ORGANISM-Homo_sapiens-3.5.179.mitab.txt'
    with open(ppi_file) as f:
        next(f) # skip the head line
        data = f.readlines()

    id_l, interact_l = [], []
    for line in data:
        line = line.split('\t')
        if '0915' in line[11] or '0407' in line[11]:
            id_1 = line[2].split('|')[0].split(':')[1]
            id_2 = line[3].split('|')[0].split(':')[1]
            # The format of line[2] and line[3]: biogrid:100010|... 
            if id_1 == id_2:
                continue

            id_l.append(id_1)
            id_l.append(id_2)
            interact_l.append((id_1, id_2))

    id_l = list(set(id_l))
    id_l.sort()
    index_map = dict(zip(id_l, range(len(id_l))))
    a= np.zeros((len(id_l), len(id_l)))

    for i in interact_l:
        p1, p2 = i
        x, y = index_map[p1], index_map[p2]
        a[x, y] = a[y, x] = 1

    out_path = '../data/bio_ids.txt'
    if not Path(out_path).exists():
        with open(out_path, 'w') as f:
            for row in id_l:
                f.write(row + '\n')
    return id_l, a

def get_uni_list():
    with open('../data/uni_ids.txt') as f:
        l = f.read().split()
    return l

def construct_and_return_uni_net():
    '''
    This function uses file 'bio_ids.txt' and file 'bio-uni.txt'.
    'bio-uni.txt' is downloaded from uniprot using 'bio_ids.txt'
    依赖"../data/bio-uni.txt"文件，
    '''
    uni_id_set = set() 
    d = {}
    with open('../data/bio-uni.txt') as f:
        next(f)
        data = f.read().split('\n')[:-1]
    for line in data:
        uni_id, bio_ids = line.strip().split('\t')
        uni_id_set.add(uni_id)
        bio_ids = bio_ids.split(',') 
        # data has lines such that:
        # P04745    106773,106774,106775
        for bio_id in bio_ids:
            if bio_id in d:
                d[bio_id] += '|' + uni_id
            else:
                d[bio_id] = uni_id

    uni_list = list(uni_id_set)
    uni_list.sort()
    index_map = dict(zip(uni_list, range(len(uni_list))))
    #get uni_list and d: the id map dict

    b = np.zeros((len(uni_list), len(uni_list)))
    bio_list, a = construct_and_return_bio_net()
    indices = list(zip(*np.where(a)))
    for x1, y1 in indices:
        bio1, bio2 = bio_list[x1], bio_list[y1]
        try:
            u1, u2 = d[bio1], d[bio2]
        except KeyError:
            continue
        try:
            x2, y2 = index_map[u1], index_map[u2]
            b[x2, y2] = 1
        except KeyError:
            for u1 in u1.split('|'):
                for u2 in u2.split('|'):
                    x2, y2 = index_map[u1], index_map[u2]
                    b[x2, y2] = 1
    # get b: uni net array
    np.fill_diagonal(b, 0)
    # 消除孤立点
    nonzero_index = b.sum(0).nonzero()[0]
    b = b[nonzero_index, :][:, nonzero_index]
    uni_list = np.array(uni_list)
    uni_list = uni_list[nonzero_index].tolist()
    return uni_list, b

def construct_and_return_loc_net():
    def judge_GO_line(line, GO_and_childrean_list):
        is_GO_line, code_evidence, in_GO_list = False, False, False

        if line.startswith('DR   GO;') and 'C:' in line \
                and ('IDA' in line or 'HDA' in line)\
                and line[9:19] in GO_and_childrean_list:
            return True
        else: return False

    with open('../data/GO_list.txt') as f:
        GO_list = f.read().split()

    with open('../data/uni_ids.txt') as f:
        uni_list = f.read().split()
        s = set(uni_list)

    uniprot_file = '../data/uniprot_sprot_human.dat.gz'
    with gzip.open(uniprot_file) as f:
        data = f.read().decode()
    entry_list = data.split('//\n')
    d = {}

    for entry in entry_list:
        AC, CC_list = None, []
        lines = entry.split('\n')
        for line in lines:
            if not AC:
                if line.startswith('AC'):
                    AC = line.split()[1].replace(';', '')
            elif judge_GO_line(line, GO_list):
                CC = line[9:19]
                CC_list.append(CC)
        if AC in s and CC_list:
            d[AC] = CC_list
    uni_list = list(d)
    uni_list.sort()
    a = np.zeros((len(d), len(GO_list)))
    for i in range(len(d)):
        p = uni_list[i]
        for GO in d[p]:
            a[i][GO_list.index(GO)] = 1
    return uni_list, a




def ECC(a, epsilon = 0):
    res = np.zeros(a.shape)
    sum_list = a.sum(0)
    #for i in tqdm(range(len(a))):
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i, j]:
                molecule = np.logical_and(a[i], a[j]).sum()
                denomanator = min(sum_list[i], sum_list[j]) - 1
                if denomanator == 0:
                    val = epsilon
                else:
                    val = molecule / denomanator
                res[i, j] = res[j, i] = val
    return res

def my_expm(a, n_iter = 15, diag = 'laplace'):
    t = a.copy()
    np.fill_diagonal(t, 0)
    if diag == 'laplace':
        np.fill_diagonal(t, -a.sum(0))
    elif diag == 'max':
        np.fill_diagonal(t, -a.sum(0).max())
    elif diag == 'zero':
        pass
    res = np.eye(len(t))
    n = 2 ** n_iter
    res += t/n
    for i in tqdm(range(n_iter)):
        res = res @ res
    return res



def construct_and_return_dict_relate_probe_and_uniprot(direct = 0):
    '''
    direct = 0 return prob to uniprot dict;
    direct = 1 return uniprot tp prob dict.
    '''
    d0, d1 = {}, {}
    # UniprotID_list = []
    with open('../data/probe2UniprotAC_file.csv') as f:
        next(f)
        for line in f.readlines():
            line = line.strip().replace('"', '')
            _, probID, UniprotID = line.split(',')
            if probID in d0:
                d0[probID] += '|' + UniprotID
            else:
                d0[probID] = UniprotID
            if UniprotID in d1:
                d1[UniprotID] += '|' + probID
            else:
                d1[UniprotID] = probID
            # UniprotID_list.append(UniprotID)
    del d1['NA']
    res = [d0, d1]
    return res[direct]

def construct_and_return_expr_matrix(method = 'mean', ncol = 23, file_name = 'probe_expression_file.csv'):
    '''
    method parameter can be 'mean', 'max' or 'median'
    '''
    d = construct_and_return_dict_relate_probe_and_uniprot(1)
    # d is a dict from uniprotAC to probe ID
    with open('../data/' + file_name)as f:
        n = f.read().count('\n') - 1
        eset_mat = np.zeros((n, ncol))
        probe_list = []
        f.seek(0)
        next(f)
        for i in range(n):
            line = f.readline().replace('"', '')
            probe, *values = line.split(',')
            probe_list.append(probe)
            eset_mat[i] = list(map(float, values))
    index_map1 = dict(zip(probe_list, range(len(probe_list))))
    #import pdb; pdb.set_trace()

    with open('../data/uni_ids.txt') as f:
        uni_list = f.read().split()
    #uni_list = list(d)
    #uni_list.sort()
    uniprot_mat = np.zeros((len(uni_list), ncol))
    method = np.median if method == 'median' else eval('np.ndarray.' + method) 
    #for i in tqdm(range(len(d))):
    for i in range(len(uni_list)):
        prot = uni_list[i]
        try:
            mapped_probes = [x for x in d[prot].split('|') if x in index_map1]
        except KeyError:
            continue
        indices = list(map(index_map1.get, mapped_probes))
        if indices:
            uniprot_mat[i] = method(eset_mat[indices], 0)
    return uni_list, uniprot_mat



def gene_mapping(key = 'gene'):
    map_file = '../data/ID_map.tab'
    with open(map_file) as  f:
        from collections import namedtuple, defaultdict
        head = f.readline().strip().replace(' ', '_')
        Line = namedtuple('Line', head)
        d = defaultdict(list)
        for line in f:
            line = Line(*line.strip().split('\t'))
            if key == 'gene':
                for gene in line.Gene_names.split():
                    d[gene].append(line.Entry)
            elif key == 'id':
                for gene in line.Gene_names.split():
                    d[line.Entry].append(gene)
            elif key == 'id2entry':
                d[line.Entry] = line.Entry_name
    return d
