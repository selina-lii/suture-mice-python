
def loadmat(filename):
    import mat73
    import scipy.io
    try:
        data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    except:
        data = mat73.loadmat(filename, use_attrdict=True)
    return data

import numpy as np

mouseDir = 'C:\\Users\\selinali\\lab\\Mice'
mice = ['T01', 'T02', 'T03', 'T04', 'Sut1', 'Sut2', 'Sut3', 'Sut4']
import glob
names = ['B005', 'B04', 'B12', 'B24', 'B48', 'B72', 'B96',
              'S005', 'S04', 'S12', 'S24', 'S48', 'S72', 'S96',
              'U005', 'U04', 'U12', 'U24', 'U48', 'U72', 'U96', 'U120', 'U144', 'U168',
              'R005', 'R04', 'R12', 'R24', 'R48']
n_neurons = np.empty(len(mice),dtype=np.ndarray)
n_ses = np.empty(len(mice),dtype=np.ndarray)
sum_n_ses = 0
for i, mouse in enumerate(mice):
    allfiles = glob.glob(mouseDir + '\\' + mouse + '\\*_Suite2p_*.mat')
    if len(allfiles)>0:
        n_neurons[i] = np.empty(len(allfiles),dtype=np.ndarray)
    for j, file in enumerate(allfiles):
        data = loadmat(file)['suite2pData']
        nses = len(data.startIdx)
        sum_n_ses = sum_n_ses + nses
        n_neurons[i][j] = [len(data.cellIdx)]*nses
        print(n_neurons[i][j])

# the above code takes a long time to run but its data has been saved locally: #
np.save('n_neurons.npy',n_neurons,allow_pickle=True)
n_neurons=np.load('n_neurons.npy',allow_pickle=True)

n_allNeurons = int(max(max([max(v) for v in n_neurons if v is not None])) * 1.7)

inhabitants_of_table = np.empty([sum_n_ses, n_allNeurons], dtype=bool)

import random
cursor = 0
cross_day_idxs=np.empty(len(n_neurons),dtype=np.ndarray)
for i, n_neu_mice in enumerate(n_neurons):
    if n_neu_mice is not None:
        cross_day_idxs[i]=np.empty(len(n_neu_mice),dtype=np.ndarray)
        for j, n_neu_run in enumerate(n_neu_mice):
            cross_day_idxs[i][j] = sorted(random.sample(range(n_allNeurons),n_neu_run[0]))
            for _ in range(len(n_neu_run)):
                inhabitants_of_table[cursor][cross_day_idxs[i][j]]=True
                cursor=cursor+1

