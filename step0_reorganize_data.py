import numpy as np
import pickle
import glob
import magicpixie
import sys
import time
log_time = time.strftime("%Y%m%d%H%M")

###################### log #########################
old_stdout = sys.stdout
log_file = open('log_'+log_time+'.txt','w')
sys.stdout = log_file

#####################################################
def loadmat(filename):
    import mat73
    import scipy.io
    try:
        data = mat73.loadmat(filename, use_attrdict=True)
    except:
        #not using this as for now.. instead converted non-7.3 files into v7.3 for sut3/4
        data = scipy.io.loadmat(filename, struct_as_record=True, squeeze_me=True)
    return data

def run_names(data):
    runs = []
    hours = ['0025', '005', '01', '04', '06', '12', '24', '48', '72', '96', '120', '144', '168']
    condition = ['B', 'S', 'U', 'R']
    for data_col,cond in zip(data,condition):
        for d,hour in zip(data_col,hours):
            if d==1:
                runs.append(cond+hour)
    return runs


#mouseDir = 'D:\\333_Galaxy Maotuan\\I love studying\\2022 winter\\lab\\papers\\2022-2-18'
mouseDir = 'C:\\Users\\selinali\\lab\\Mice'
mice = ['Sut4', 'Sut1', 'Sut2', 'Sut3', 'T01', 'T02', 'T03', 'T04']
# original order - mice = ['T01', 'T02', 'T03', 'T04', 'Sut1', 'Sut2', 'Sut3', 'Sut4']


with open('mice_days.pickle', 'rb') as file:
    mice_days=pickle.load(file)
    reo=[7,4,5,6,0,1,2,3]
    mice_days=mice_days[reo,:,:]
    mice_n_files = np.sum(np.sum(mice_days,axis=2),axis=1)

for i, mouse in enumerate(mice):
    allfiles = glob.glob(mouseDir + '\\' + mouse + '\\*_Suite2p_*.mat')
    print(mouse)
    if len(allfiles)>0:
        assert mice_n_files[i] == len(allfiles)
        print(mouse)
        names = run_names(mice_days[i])
        for j, file in enumerate(allfiles):
            if mouse=='Sut1' and j>14:
                continue
            if not glob.glob("%s_%s_%s_" % (mouse, str(j).zfill(2), names[j]) + '*.pickle'):
                print(names[j] + ' ' + str(j))
                data = loadmat(file)['suite2pData']
                magicpixie.pixiedust(data=data,name=names[j],idx=j,cut_sessions=True if i>3 else False)

        '''
        n_ses = len(data.startIdx)
        sum_n_ses = sum_n_ses + n_ses
        n_neurons[i][j] = [len(data.cellIdx)]*n_ses
        print(n_neurons[i][j])
        '''

sys.stdout = old_stdout
log_file.close()

from magicpixie import *
acrossdayidxs = np.swapaxes(loadmat('Sut4CellIndex.mat')['cellindexalign'], 0,1)
acrossdayidxs=np.delete(acrossdayidxs,7,0)
neurons_across_days = np.empty([acrossdayidxs.shape[0]*2,acrossdayidxs.shape[1]],dtype=Neuron)

allfiles = glob.glob('Sut4_*.pickle')
cnt = 0
for filename, dayidxs in zip(allfiles, acrossdayidxs):
    print(filename)
    with open(filename, 'rb') as file:
        run=pickle.load(file)
    file.close()
    neuronidxs = [x-1 for x in dayidxs if x > 0]
    temp = np.empty([len(neuronidxs),2],dtype=Neuron)
    print(len(neuronidxs)==len(run.sessions[0].neurons))
    for session in run.sessions[0:2]:
        for i, idx in enumerate(dayidxs):
            idx=idx-1
            if idx>=0:
                neurons_across_days[cnt, i] = session.neurons[idx]
                #neuron.plotTrace()
                #neuron.plotTrialsMat()
        cnt=cnt+1

