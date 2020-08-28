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
	pred = pred[indices]
	y = y[indices]

	proba = np.load('../tem_data/hepatitis_probas/probas_nor_tau_1.1.npy')
	proba = proba_scaling(proba, 'lin')
	pred2 = calibrate(proba, 0.3)
	pred2 = pred2[:, np.arange(6, 11)]
	pred2 = pred2[indices]

	def _print():
		print('PPV: {}, SE: {}, SP: {}, ACC2: {}, MCC: {}, AIM: {}, CVR: {}, ACC1: {}, ATR: {}, AFT: {}'.format(*multi))
		print('PPV:',*PPV)
		print('SE:',*SE)
		print('SP:',*SP)
		print('ACC:',*ACC)
		print('MCC:',*MCC)
		print()
	multi = multi_label_performance(pred, y)
	PPV, SE, SP, ACC, MCC = performance(pred, y)
	print('Performance of Hum-ploc:')
	_print()

	multi = multi_label_performance(pred2, y)
	PPV, SE, SP, ACC, MCC = performance(pred2, y)
	print('Performance of our method:')
	_print()
def gene_tfpn(pred, y_test):
	TP = np.logical_and(pred, y_test).sum(0)
	TN = np.logical_and(np.logical_not(pred), np.logical_not(y_test)).sum(0)
	FP = np.logical_and(pred, np.logical_not(y_test)).sum(0)
	FN = np.logical_and(np.logical_not(pred), y_test).sum(0)
	print('TP:', *TP, sep = ',')
	print('TN:', *TN, sep = ',')
	print('FP:', *FP, sep = ',')
	print('FN:', *FN, sep = ',')
	print('Global:')
	print('TP: {}, TN: {}, FP: {}, FN: {}'.format(TP.sum(), TN.sum(), FP.sum(), FN.sum()))
	print()

def print_tfpn():
	y = np.loadtxt('../tem_data/probas_txt/label_mat.txt')
	with open('../data/location_name.txt') as f:
		names = f.read().split()

	proba = np.load('../tem_data/hepatitis_probas/probas_nor_tau_1.1.npy')
	proba = proba_scaling(proba, 'lin')
	pred2 = calibrate(proba, 0.3)
	print(' ', *names, sep =  ',')
	gene_tfpn(pred2, y)

	pred = gene_loc_array()
	pred = pred[:, np.arange(6, 11)]
	y = y[:, np.arange(6, 11)]
	indices = y.any(1)
	pred = pred[indices]
	y = y[indices]

	pred2 = pred2[:, np.arange(6, 11)]
	pred2 = pred2[indices]
	print(' ', *names[6:11], sep =  ',')
	gene_tfpn(pred, y)
	print(' ', *names[6:11], sep =  ',')
	gene_tfpn(pred2, y)
def print_MCC_tpfn():
    cancer_names = ['hepatitis', 'breast', 'leukemia']
    _, y = construct_and_return_loc_net()
    for cancer_name in cancer_names:
        for tau in range(1, 21):
            tau /= 10
            for threshold in [0.1, 0.2, 0.3]:
                print('cancer_name:', cancer_name)
                print('tau:', tau)
                print('threshold:', threshold)
                proba = np.load('../tem_data/{}_probas/probas_nor_tau_{}.npy'.format(cancer_name, tau))
                proba = proba_scaling(proba, 'lin')
                pred = calibrate(proba, threshold)
                gene_tfpn(pred, y)



	
if __name__ == '__main__':
	print_comparation_record()
	print_tfpn()
        #print_MCC_tpfn()
