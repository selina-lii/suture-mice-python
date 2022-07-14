'''from scipy.io.wavfile import write
write('test.wav', 1563, data.spks[0,:])'''


plotFolder = 'C:\\Users\\selinali\\lab\\RDE20141\\2022-7-7-plots'


def loadmat(filename):
    import mat73
    import scipy.io
    try:
        data = mat73.loadmat(filename, use_attrdict=True)
    except:
        #not using this as for now.. instead converted non-7.3 files into v7.3 for sut3/4
        data = scipy.io.loadmat(filename, struct_as_record=True, squeeze_me=True)
    return data


import glob
from magicpixie import *
#data = np.swapaxes(loadmat('Sut4CellIndex.mat')['cellindexalign'], 0,1)
data = loadmat('Sut4CellIndex.mat')['cellindexalign']
data=np.delete(data,7,1) #7,0
acrossdayidxs=np.zeros_like(data,dtype=int)
for i,dd in enumerate(data):
    for j,d in enumerate(dd):
        acrossdayidxs[i,j]=d-1 if d>0 else -1

'''neurons_across_days = np.empty([acrossdayidxs.shape[0]*2,acrossdayidxs.shape[1]],dtype=Neuron)

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
        cnt=cnt+1'''

import cv2

class plotType:
    name=None
    framerate=None
    compressionRatio=1
    perSes=True
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)

crossdayFolder='C:\\Users\\selinali\\lab\\RDE20141\\across_day_Sut4_' + statName

def plotCrossdayWrapper():
    plottypes=np.empty(3,dtype=plotType)
    plottypes[0]=plotType(name='trace',framerate=2,compressionRatio=3)
    plottypes[1]=plotType(name='trialmat',framerate=2)
    plottypes[2]=plotType(name='roi',framerate=6,perSes=False)
    for plottype in plottypes:
        if plottype.perSes:
            for
        plotCrossday(,,,)

def plotCrossday():
    statName=plottype.name
    framerate=plottype.framerate
    compressionRatio=plottype.compressionRatio
    ses=

    img_example=cv2.imread(plotFolder+'\\Sut4_00_B005_191117_ses1_'+statName+'\\Sut4_00_B005_191117_ses1_#0_'+statName+'.jpg')
    w=img_example.shape[1]//compressionRatio
    h=img_example.shape[0]//compressionRatio
    blankscreen=np.zeros((w, h, 3), dtype = "uint8")
    blankscreen.fill(255)

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')

    mouse='T01'
    folders=glob.glob(plotFolder + '\\' + mouse + '*ses0*' + statName)
    for cd_i,reorder in enumerate(acrossdayidxs[0:100]):
        out = cv2.VideoWriter(crossdayFolder+'\\'+str(cd_i)+'.mp4',fourcc,framerate,[w,h])
        for i, f in zip(reorder, folders):
            if i==-1:
                img = blankscreen
            else:
                flast = f.split('\\')[-1].rsplit('_', 1)[0]
                img = cv2.imread(f + '\\' + flast + '_#' + str(i) + '_'+statName+'.jpg')
            img=cv2.resize(img,size)
            out.write(img)
        out.release()


#organize this
# do spontaneous data
# across days
#mark up overlap

#should run have image size????
#still have to do database for our metadata and katie's....
#create a blank refimg w ghost neuron or separate function..

import os, os.path
for fol in os.listdir(plotFolder):
    print(fol+'  '+str(len(glob.glob(plotFolder+'\\'+fol+'\\'+'*.jpg'))))
# rerun T03 later days
# Sut1 resize trialsmat
#build a mongodb.
    #all these should be method classes, with data stored in the database
    # two databases: 1.meta(including filenames) 2.dff
    # plotting info is still hardcoded
    # add user input ports -> GUI
    # after putting in a database we combine sessions for plots much more easily.
    # neuron.overview() will generate a single info card/profile per neuron that's comparable between neurons
    # neuron_crossday probably will have database ids instead of real objects. most of the other classes will be either
    # pure methods class or have object pointers to the database

import pymongo
myclient=pymongo.MongoClient("mongodb://localhost:27017/")
mydb=myclient["sut_mice"]
mycol=mydb["mouse"]
print(mydb.list_collection_names())
print(myclient.list_database_names())

