from main import *
import numpy as np
with open('../data/location_name.txt') as f:
    locations = f.read().split('\n')[:-1]
    if len(locations) != 12:
        print('Error')
    locations = np.array(locations)

l, y = construct_and_return_loc_net()
y = (y!=0)
with open('../tem_data/loc_annotations.txt', 'w') as f:
    for i in range(len(l)):
        f.write(l[i] + '\t' + '|'.join(locations[y[i]]) + '\n')


