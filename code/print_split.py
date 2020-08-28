from main import *
import numpy as np
from sklearn.model_selection import KFold
kf = KFold(n_splits = 10, random_state = 0)
with open('../data/uni_ids_with_loc.txt') as f:
    data = f.read().split()
print('start_index', 'end_index')
for train_index, test_index in kf.split(data):
    print(test_index[0], test_index[-1])

