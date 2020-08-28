with open('covdesc') as f: 
    data = f.read().split('\n')
res = ''
fix = '.CEL.gz'
oneline = data[0] + fix
for i in range(1, len(data)):
    term = data[i]
    if term.startswith('GSM'):
        term += fix
    if ' ' in term:
        term = term.replace(' ', '_')
    if oneline != '':
        oneline += '\t'
    oneline += term
    try:
        if data[i + 1].startswith('GSM'):
            res += '\n' + oneline
            oneline = ''
    except IndexError:
            res += '\n' + oneline
with open('covdesc.txt', 'w') as f:
    f.write('\tstate\ttype\ttitle')
    f.write(res)


