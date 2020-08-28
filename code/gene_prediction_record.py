from predict import *
y = np.loadtxt('../tem_data/probas_txt/label_mat.txt')

def gene_prediction_performance_of_cancer(cancer_name):
    for i in range(1, 21):
        i /= 10
        state = 'nor'
        proba = np.load('../tem_data//{}_probas/probas_{}_tau_{}.npy'.format(cancer_name, state, i))
        proba = proba_scaling(proba, 'lin')
        print('tau {}'.format(i))
        for threshold in (0.1, 0.2, 0.3):
            print('threshold', threshold)
            pred = calibrate(proba, threshold)
            print_performances(pred, y)
        print()





if __name__ == '__main__':
    '''
    cancer_names = [
            'hepatitis',
            'leukemia',
            'breast'
            ]
    for cancer_name in cancer_names:
        print(cancer_name + ':')
        gene_prediction_performance_of_cancer(cancer_name)
        '''
    import sys
    cancer_name = sys.argv[1]
    gene_prediction_performance_of_cancer(cancer_name)
