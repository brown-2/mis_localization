import os
import math
from time import perf_counter
from scipy.linalg import expm
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold

from sklearn.preprocessing import FunctionTransformer
from sklearn.preprocessing import power_transform, PowerTransformer
from sklearn.preprocessing import KernelCenterer
from sklearn.preprocessing import label_binarize

from sklearn.svm import SVC
from sklearn.pipeline import make_pipeline, Pipeline

from main import *
from split_wrap import *
def get_K_from_e(e):
    '''
    用来切分核，把带标签的部分切出来
    '''
    with open('../data/locative_uni_ids.txt') as f:
        l_ppi = f.read().split()
    index_map = dict(zip(l_ppi, range(len(l_ppi))))
    with open('../data/locative_uni_ids_with_loc.txt') as f:
        l_loc = f.read().split()
    #l_loc, _ = construct_and_return_loc_net()
    
    indices = np.array(list(map(index_map.get, l_loc)))
    X = e[indices, :][:, indices]
    return X

def split_dataset(X, y, test_size = 0.1, random_state = 100):
    '''
    用来划分训练集和测试集
    '''
    indices_train, indices_test, X_train, X_test, y_train, y_test = \
            train_test_split(np.arange(len(X)), X, y, test_size = test_size,\
            random_state = random_state)
    X_train = X_train[:, indices_train]
    X_test = X_test[:, indices_train]
    return X_train, y_train, X_test, y_test

def train_test(X_train, y_train, X_test, params):
    clf = SVC(**params)
    clf.fit(X_train, np.argmax(y_train, 1))
    '''
    print(np.bincount(clf.predict(X_train).astype(int)))
    print(clf.predict_proba(X_train).max(0))
    print(clf.predict_proba(X_train).sum(0))
    '''
    print('score:', clf.score(X_train, np.argmax(y_train, 1)))
    proba = clf.predict_proba(X_test)
    return proba

def performance(pred, y_test):
    TP = np.logical_and(pred, y_test).sum(0)
    TN = np.logical_and(np.logical_not(pred), np.logical_not(y_test)).sum(0)
    FP = np.logical_and(pred, np.logical_not(y_test)).sum(0)
    FN = np.logical_and(np.logical_not(pred), y_test).sum(0)
    PPV = TP / (TP + FP)
    NPV = TN / (TN + FN)
    SE = TP / (TP + FN)
    SE[np.isnan(SE)] = 0
    SP = TN / (TN + FP)
    ACC = (TP + TN) / (TP + TN + FP + FN)
    MCC = (TP * TN - FP * FN)/\
            ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5
    MCC[np.isnan(MCC)] = 0
    return np.round(PPV, 3), np.round(SE, 3), np.round(SP, 3), np.round(ACC, 3), np.round(MCC, 3)
def multi_label_performance(pred, y_test):
    TP_m = np.logical_and(pred, y_test)
    TN_m = np.logical_and(np.logical_not(pred), np.logical_not(y_test))
    FP_m = np.logical_and(pred, np.logical_not(y_test))
    FN_m = np.logical_and(np.logical_not(pred), y_test)

    AIM_t = TP_m.sum(1)/pred.sum(1)
    AIM_t[np.isnan(AIM_t)] = 0
    AIM = AIM_t.mean()

    CVR = (TP_m.sum(1)/y_test.sum(1)).mean()
    ACC1 = (TP_m.sum(1)/np.logical_or(pred, y_test).sum(1)).mean()
    ATR = 1 - np.any(pred - y_test, 1).mean()
    AFR = (np.logical_or(pred, y_test).sum(1) - TP_m.sum(1)).sum()/\
            (pred.shape[0] * pred.shape[1])

    TP, TN, FP, FN = TP_m.sum(), TN_m.sum(), FP_m.sum(), FN_m.sum()
    PPV = TP / (TP + FP)
    SE = TP / (TP + FN)
    SP = TN / (TN + FP)
    ACC2 = (TP + TN) / (TP + TN + FP + FN)
    MCC = (TP * TN - FP * FN)/\
            ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5
    return [np.round(x, 3) for x in (PPV, SE, SP, ACC2, MCC, AIM, CVR, ACC1, ATR, AFR)]


def KFold_predict(X, y, params):
    kf = KFold(n_splits = 10, random_state = 100)
    y_pred = np.zeros(y.shape)
    y_pred_proba = np.zeros(y.shape)
    for indices_train, indices_test in kf.split(X):
        X_train = X[indices_train, :][:, indices_train]
        X_test = X[indices_test, :][:, indices_train]
        y_train = y[indices_train]
        y_test = y[indices_test]
        pred, proba = train_test(X_train, y_train, X_test, C)
        #pred[np.arange(len(pred)), proba.argmax(1)] = 1
        y_pred[indices_test] = pred
        y_pred_proba[indices_test] = proba
    return y_pred, y_pred_proba

def KFold_predict_locative(X, X_ill, y_new, y, params):
    kf = KFold(n_splits = 10, random_state = 100)
    y_nor_proba = np.zeros(y.shape)
    y_ill_proba = np.zeros(y.shape)
    for indices_train, indices_test in kf.split(y):
        indices_train_locative= real2virtual_map(indices_train)
        X_train = X[indices_train_locative, :][:, indices_train_locative]
        y_train = y_new[indices_train_locative]
        clf = SVC(**params)
        clf.fit(X_train, y_train.argmax(1))

        indices_test_locative= real2virtual_map_unique(indices_test)
        X_test = X[indices_test_locative, :][:, indices_train_locative]
        X_ill_test = X_ill[indices_test_locative, :][:, indices_train_locative]

        proba_nor = clf.predict_proba(X_test)
        proba_ill = clf.predict_proba(X_ill_test)
        y_nor_proba[indices_test] = proba_nor
        y_ill_proba[indices_test] = proba_ill
    return y_nor_proba, y_ill_proba

def print_performances(y_pred, y):
    Multi = multi_label_performance(y_pred, y)
    PPV, SE, SP, ACC, MCC = performance(y_pred, y)
    print(*Multi)
    print(*PPV)
    print(*SE)
    print(*SP)
    print(*ACC)
    print(*MCC)

def calibrate(proba, C):
    pred = np.zeros(proba.shape)
    thresholds = proba.max(1) * (1 - C) + proba.min(1) * C
    for i in range(len(proba)):
        threshold = thresholds[i]
        pred[i][proba[i] >= threshold] = 1
    return pred

def proba_scaling(proba, scaling_type = 'lin'):
    '''
    scaling_type 参数可以是 'lin' 或 'nor'
    '''
    if scaling_type == 'lin':
        proba_lin = proba.copy()

        # 减去最小值
        p_min = proba_lin.min(0)
        for j in range(proba.shape[1]):
            proba_lin[:,j] = proba_lin[:,j] - p_min[j]

        # 除以最大值 
        p_max = proba_lin.max(0)
        for j in range(proba.shape[1]):
            proba_lin[:,j] = proba_lin[:,j] / p_max[j]

        # 行归一化
        sum_row = proba_lin.sum(1)
        for i in range(proba.shape[0]):
            proba_lin[i] /= sum_row[i]
        return proba_lin

    elif scaling_type == 'nor':
        # 列规范化
        proba_nor = proba.copy()
        p_mean = proba_nor.mean(0)
        p_std = proba_nor.std(0)
        for j in range(proba.shape[1]):
            proba_nor[:,j]  = (proba_nor[:,j] - p_mean[j])/p_std[j]

        # 行归一化
        sum_row = proba_nor.sum(1)
        for i in range(proba.shape[0]):
            proba_nor[i] /= sum_row[i]
        return proba_nor
    else:
        print('scaling_type parament must be lin or nor')
        raise TypeError

if __name__ == '__main__':
    params = {
            'kernel':'precomputed',
            'C':1000,
            'cache_size':4000,
            'class_weight':'balanced',
            'probability':False,
            'decision_function_shape':'ovr',
            }
    e = np.load('../tem_data/e_raw_locative.npy')
    X = get_K_from_e(e)
    l, y = construct_and_return_loc_net()

    #加权重
    weight = y.sum(0).max()/y.sum(0)
    #weight **= 2
    params['class_weight'] = dict(zip(range(12), weight))
    print(weight.astype(int))

    l_new, y_new, d = locative_loc_net(l, y)
    #proba = KFold_predict_locative(X, y_new, y, params)
    proba = KFold_predict_dec(X, y_new, y, params)
    p_mean = proba.mean(0)
    p_std = proba.std(0)
    p_max = proba.max(0)
    proba_nor = proba.copy()
    proba_lin = proba.copy()
    for j in range(proba.shape[1]):
        proba_nor[:,j]  = (proba_nor[:,j] - p_mean[j])/p_std[j]
        proba_lin[:,j] = proba_lin[:,j]/p_max[j]
    pred = calibrate(proba, 0.5)
    pred_nor = calibrate(proba_nor, 0.5)
    pred_lin = calibrate(proba_lin, 0.5)
