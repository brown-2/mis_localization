name = 'covdesc'
with open(name) as f:
    data = f.read().split('\n')
    data = [x for x in data if x != '' and x !=  '-']
for i in range(len(data)):
    x,y,*tem = data[i].split('\t', 2)
    y = y.replace('(', '').replace(')', '')
    state, H, num = y.split()
    x += '_' + H + '_' + num + '.CEL.gz'
    data[i] = x + '\t' + state
with open('t.txt','w') as f:
    f.write('\t' + 'state' + '\n')
    for row in data:
        print(row)
        f.write(row + '\n')

