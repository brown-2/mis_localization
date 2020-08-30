import numpy as np
from main import *
from predict import *
def gene_loc_dict():
	with open('../data/prediction_humploc.txt') as f:
		head, _= f.read().split('//')
	Loc2GO = {}
	for i in head.split('\n'):
		if i:
			name, GO = i.split('\t')
			Loc2GO[name] = GO
	with open('../data/GO_list.txt') as f:
		locs = f.read().split()
	Loc2index = {}
	for loc in Loc2GO:
		if Loc2GO[loc] in locs:
			Loc2index[loc] = locs.index(Loc2GO[loc])
	return Loc2index

def gene_loc_array():
	with open('../data/uni_ids_with_loc.txt') as f:
		data = f.read().strip().split()
		prot2index = dict(zip(data, range(len(data))))

	with open('../data/prediction_humploc.txt') as f:
		_, data= f.read().split('//')
		data = data.strip().split('\n')
	a = np.zeros((len(prot2index), 12))
	Loc2index = gene_loc_dict()
	for line in data:
		ID, locs = line.split('\t')
		try:
			index = prot2index[ID]
		except KeyError:
			continue
		for loc in locs.strip().replace('.', '').split():
			try:
				loc_index = Loc2index[loc]
				a[index, loc_index] = 1
			except KeyError:
				continue
	return a
def print_comparation_record():
    pred = gene_loc_array()
    pred = pred[:, np.arange(6, 11)]
    y = np.loadtxt('../tem_data/probas_txt/label_mat.txt')
    y = y[:, np.arange(6, 11)]
    indices = y.any(1)
    ids = np.array(get_annotated_ids())
    ids = ids[indices]
    for i in ids:
        print(i)
    print(len(ids))




	
if __name__ == '__main__':
	print_comparation_record()
