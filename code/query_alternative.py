from main import *
def preprocess():
    map_file = '../data/ID_map.tab'
    with open(map_file, 'r+') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i]
            line = line.split('\t')[:-1]
            line[-1] += '\n'
            line = '\t'.join(line)
            lines[i] = line
    with open(map_file, 'w') as f:
        f.writelines(lines)

def get_diff(proteins, threshold = 0):
    l = get_annotated_ids()
    index_mapping = dict(zip(l, range(len(l))))
    d = gene_mapping('id')
    for cancer_name in ['hepatitis', 'breast', 'leukemia']:
        for i in range(1, 21):
            i /= 10
            path = '../tem_data/{}_records/diff_{}.npy'.format(cancer_name, i)
            a = np.load(path)
            for protein in proteins:
                index = index_mapping[protein]
                diffs = a[index]
                if (np.abs(diffs)>threshold).any():
                    print('in file:', path)
                    print(*d[protein], *np.round(diffs, 3))
                



if __name__ == '__main__':
    import sys
    query_file = sys.argv[1]
    with open(query_file) as  f:
        proteins = f.read().replace(',', ' ').split()

    d = gene_mapping('gene')
    found = []
    not_founds = []
    for i in proteins:
        if i in d:
            found.extend(d[i])
        else:
            not_founds.append(i)
    #print('Following proteins not in dict:',*not_founds)
    print('Following proteins found:')
    get_diff(found, 0.5)
