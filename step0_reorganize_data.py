import numpy as np
import pickle
import glob
import magicpixie
import sys
import time
log_time = time.strftime("%Y%m%d%H%M")

###################### log #########################
#old_stdout = sys.stdout
#log_file = open('log_'+thismouse+'_'+str(thisrun)+'_'+log_time+'.txt','w')
#sys.stdout = log_file

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
#mice = ['T03','T02','Sut3']
mice = ['T01', 'T02', 'T03', 'T04', 'Sut1', 'Sut2', 'Sut3', 'Sut4']
plotFolder = 'C:\\Users\\selinali\\lab\\RDE20141\\2022-7-7-plots'

with open('mice_days.pickle', 'rb') as file:
    mice_days=pickle.load(file)
    reo=[2,1,6]
    mice_days=mice_days[reo,:,:]
    mice_n_files = np.sum(np.sum(mice_days,axis=2),axis=1)

for i, mouse in enumerate(mice):
    #if mouse==thismouse:
        allfiles = glob.glob(mouseDir + '\\' + mouse + '\\*_Suite2p_*.mat')
        print(mouse)
        if len(allfiles)>0:
            assert mice_n_files[i] == len(allfiles)
            print(mouse)
            names = run_names(mice_days[i])
            for j, file in enumerate(allfiles):
                print(names[j] + ' ' + str(j))
                '''                filename = glob.glob("%s_%s_%s_" % (mouse, str(j).zfill(2), names[j]) + '*.pickle')[0]
                                with open(filename,'wb') as file:
                                    run=pickle.load(file)
                                    roiMask = session.neurons[0].roiMask
                                    roiMask.circs=None
                                    roiMask.med=None
                                    roiMask.fills=None
                                    roiMask.box=None
                                    for session in run.sessions:
                                        session.ghostNeuron=magicpixie.Neuron(idx=-1, idx_glob=-1, sessions=session, runs=j, roiMask=roiMask)
                '''
                ######################################### 7-12 #############################################################
                #if j==thisrun:
                try:
                    filename=glob.glob("%s_%s_%s_" % (mouse, str(j).zfill(2), names[j]) + '*.pickle')[0]
                except:
                    data = loadmat(file)['suite2pData']
                    magicpixie.pixiedust(data=data,name=names[j],idx=j)
                    filename = glob.glob("%s_%s_%s_" % (mouse, str(j).zfill(2), names[j]) + '*.pickle')[0]
                magicpixie.sparkles(filename=filename,plotFolder=plotFolder,cut_sessions=True if reo[i]>3 else False, roi=True,trace=True,trialsMat=True)
                ######################################### 7-12 #############################################################
                '''
                n_ses = len(data.startIdx)
                sum_n_ses = sum_n_ses + n_ses
                n_neurons[i][j] = [len(data.cellIdx)]*n_ses
                print(n_neurons[i][j])'''

#sys.stdout = old_stdout
#log_file.close()

