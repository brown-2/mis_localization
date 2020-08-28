with open('covdesc') as f: 
    data = f.read()
fix = '.CEL.gz'
data = data.replace('benign breast abnormalities', 'benign_breast_abnormalities').replace(\
        'ectopic (gastrointestinal and brain) cancers', 'ectopic_cancers').replace(\
        'malignant breast cancer', 'malignant_breast_cancer').replace(\
        'Pre-Surgery (malignant)', 'Pre-Surgery_(malignant)').replace(\
        'Post-Surgery (malignant)', 'Post-Surgery_(malignant)').replace(\
        'training set', 'training_set').replace(\
        'validation set', 'validation_set')

res = ''
for line in data.split('\n'):
    line = line[:9] + fix + line[9:]
    res += line + '\n'
with open('covdesc.txt', 'w') as f:
    f.write('\tstate\ttype\ttitle\n')
    f.write(res)

