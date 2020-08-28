from main import *
from split_wrap import *

def BioNum():
    uniprot_file = '../data/uniprot_sprot_human.dat.gz'
    with gzip.open(uniprot_file) as f:
        data = f.read().decode()
    entry_list = data.split('//\n')
    return len(entry_list)

l_bio, a_bio = construct_and_return_bio_net()
print('The number of biogrid proteins:', len(l_bio))
print('The number of biogrid edges:', a_bio.sum()//2)

n = BioNum()
print('The number of uniprot proteins:', n)

l_uni, a_uni = construct_and_return_uni_net()
print('The number of common proteins:', len(l_uni))
print('The number of final edges:', a_uni.sum()//2)


l_loc, y_loc = construct_and_return_loc_net()
print('The number of labeled proteins:', len(l_loc))

l_new, loc_new, d = locative_loc_net(l_loc, y_loc)
print('The number of locative labeled proteins:', len(l_new))

l_new, a_new = locative_laplacian(l_uni, a_uni, d)
print('The number of locative uniprot proteins:', len(l_new))
print('The number of locative edges:', a_new.sum()//2)

