from dataclasses import dataclass


class Grid:
    r=None
    c=None
    rlabels=None
    clabels=None
    reorder=None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def split(self, matrix, pool=True):
        splits = np.array([np.hsplit(x, self.c[:-1]) for x in np.split(matrix, self.r[:-1])], dtype=object)
        if pool:
            for r in range(splits.shape[0]):
                for c in range(splits.shape[1]):
                    splits[r][c] = splits[r][c].flatten()
        return splits


class StimObject:
    label=None
    type=None
    color=None
    index=None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        return "StimObject:%s %s°(%s)" % (self.type, str(self.label), str(self.index))

    def __str__(self):
        return "%s %s°" % (self.type, str(self.label))

    def __eq__(self, other):
        return self.label == other.label

    def __lt__(self, other):
        return self.label < other.label


class StimEvent:
    on=None
    off=None
    trialOn=None
    trialOff=None
    stim=None
    ses=None
    pos_ses=None
    pos_glob=None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        return "StimEvent: %s (on=%s:off=%s)" % (self.stim.__repr__(), str(self.on), str(self.off))

    def __eq__(self, other):
        return self.stim == other.stim

    def __lt__(self, other):
        return self.stim < other.stim

    def trialLen(self):
        return self.trialOff - self.trialOn

    def stimLen(self):
        return self.off - self.on

    def setTrialOnOff(self, shifton, shiftoff, stimLen):
        self.trialOn = self.on-shifton
        self.trialOff = self.on+stimLen+shiftoff

    def zero(self,shift):
        self.on = self.on - shift
        self.off = self.off - shift
        if self.trialOn is not None:
            self.trialOn = self.trialOn - shift
        if self.trialOff is not None:
            self.trialOff = self.trialOff - shift


class StimTrace:
    stimEvents=None
    end=None
    stimObjects=None
    stimLen=None
    minTrialLen=None

    reorder=None
    grid=None
    trialLen=None

    def __init__(self, stimEvents=None,end=None,stimObjects=None):
        self.stimEvents = stimEvents
        self.end = end
        self.stimObjects = stimObjects
        self.stimLen = max([x.stimLen() for x in self.stimEvents])
        self.minTrialLen = min([x1.on - x.on for (x1, x) in zip(self.stimEvents[1:], self.stimEvents[:-1])])

        self.reorder = None
        self.grid = None
        self.trialLen = None

    def trace(self):
        trace = [0] * self.end
        for event in self.stimEvents:
            trace[event.on:event.off] = 1
        import matplotlib.pyplot as plt
        plt.clf()
        plt.plot(trace)
        plt.show()
        plt.figure().close()
        return trace

    def cutouts(self, imagingTrace, flat=False):
        # keep or cut off selected segments
        # flat or layered
        heatmap = np.empty([len(self.stimEvents), self.trialLen])
        for i,event in enumerate(self.stimEvents):
            heatmap[i] = imagingTrace[event.trialOn:event.trialOff]
        return heatmap if not flat else heatmap.flatten()

    def reorderByStim(self,matrix):
        return matrix[self.reorder]

    def stimObjectsLabels(self,statName,values):
        label = statName+': '
        for obj, val in zip(self.stimObjects, values):
            label = label + str(obj) + ' = ' + "%.2f" % val + '; \n'
        return label

    def setGrid(self,shiftOn,shiftOff):
        for event in self.stimEvents:
            event.setTrialOnOff(shiftOn, shiftOff, self.stimLen)
        if event.trialOff > self.end:
            self.stimEvents = self.stimEvents[:-1]
            print('Last trial exceeds total number of frames by %s frames with current shiftOff = %s. '
                  'Deleting the last trial.' % (str(event.trialOff - self.end),str(shiftOff)))
        unq, unq_cnt = np.unique(sorted([x.stim for x in self.stimEvents]), return_counts=True)
        self.grid = Grid(r=np.cumsum(unq_cnt),
                         c=np.cumsum([shiftOn, self.stimLen, shiftOff]),
                         rlabels=[str(x) for x in unq],
                         clabels=['stimOn', 'stimOff'])
        self.trialLen = self.stimEvents[0].trialLen()
        self.reorder = np.argsort(a=self.stimEvents, kind='stable')


class Session:  # or neuron??
    '''convert between suite2pdata and internal modules'''


    def __init__(self, index=None, stimTrace=None, neurons=None, end_glob=None, start_glob=None, end=None,
                 neurons_vis_driven=None, type=None):
        self.index = index
        self.type = type
        self.stimTrace = stimTrace
        self.end_glob = end_glob
        self.start_glob = start_glob
        self.neurons = neurons
        self.neurons_vis_driven = neurons_vis_driven
        if end is not None:
            self.end = end
        else:
            self.end = end_glob - start_glob


class Neuron:
    dff = None
    idx = None
    idx_glob = None
    sessions = None
    runs = None

    # qc and stuff
    snr = None
    spks = None
    auc = None

    # plotting properties
    roiMask = None
    trialsMat = None
    avgRes = None

    # statistics
    visdriven = None
    js_dists = None

    plot_label = None
    plot_title = None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __str__(self):
        return str(self.idx_glob)

    def setPlotTitle(self,runInfo):
        self.plot_title = runInfo

    #dff = random.choices(range(-1, 5), k=14000)
    def plotTrace(self,save_folder,stretch=4000):
        import matplotlib.pyplot as plt
        n_splits=len(self.dff)//stretch+1
        splits=np.array_split(self.dff,n_splits)
        plt.clf()
        fig = plt.figure(figsize=(n_splits*stretch//200, n_splits*stretch//1000))
        for i,split in enumerate(splits):
            plt.subplot(n_splits, 1, i+1)
            plt.tight_layout()
            plt.plot(split)

        statName='trace'
        filename = self.plot_title + '# ' + str(self.idx) + '_' + statName + '.jpg'
        save_folder = save_folder + '\\' + self.plot_title + '_' + statName
        import os
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        plt.savefig(save_folder + '\\' + filename)
        plt.figure().clear()
        return

    def plotTrialsMat(self,grid,runName,xlabel,save_folder):
        import matplotlib.pyplot as plt
        ASPECTRATIO = 4

        plt.clf()
        plt.imshow(self.trialsMat, cmap=plt.get_cmap('turbo'))
        plt.colorbar()

        width = self.trialsMat.shape[0]
        height = self.trialsMat.shape[1]
        plt_title = self.plot_title

        plt.gca().set_aspect(width / height * ASPECTRATIO)
        #plt.gca().invert_yaxis()
        plt.title(plt_title)
        plt.text(-100, height-50, xlabel)

        # gridlines
        for (pos,gridLabel) in zip(grid.c,grid.clabels):
            plt.axvline(pos, c='w', linestyle='--')
            plt.text(pos - 10, self.trialsMat.shape[0], gridLabel)
        for (pos, gridlabel) in zip(grid.r,grid.rlabels):
            plt.axhline(pos, c='w', linestyle='--')
            plt.text(0, pos, gridlabel, c='w')

        statName = 'trialsMat'
        filename = self.plot_title + '# ' + str(self.idx) + '_' + statName + '.jpg'
        save_folder = save_folder + '\\' + self.plot_title + '_' + statName
        import os
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        plt.savefig(save_folder + '\\' + filename)
        plt.figure().clear()

    def plotAvgResponse(self):
        # TODO
        return

    def plotROI(self):
        ref_img = self.sessions.refImg
        # TODO
        return


class ImagePixel:
    def __init__(self, x=None, y=None, overlap=None, trace=None, neurons=None):
        self.x = x
        self.y = y
        self.overlap = overlap
        self.trace = trace
        self.neurons = neurons


class RoiMask:
    def __init__(self, fills=None, circs=None, med=None):
        self.fills = fills
        self.circs = circs
        self.med = med


class NeuronAcrossTime:
    def __init__(self, neurons):
        self.neurons = neurons


class Run:
    sessions=None
    stimTrace=None
    mouse=None
    idx=None
    date=None
    paradigm=None # or type, e.g.long term learning / suture-resuture
    magnification=None #order of imaging
    notes=None #e.g. 'magnification on'
    bad_day=None
    cond=None
    day=None
    dprime=None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __str__(self):
        return "%s_%s_%s_%s" % (str(self.idx), self.name, self.mouse, self.date)

    def name(self):
        return self.cond + str(self.day)


class VisDriven:
    def __init__(self, stat_name=None, threshold=None):
        self.stat_name=stat_name
        self.threshold=threshold

    def js_dist(self, trialMat, grid, baselineIdxs=None, nonbaseIdxs=None):
        if baselineIdxs is None:
            baselineIdxs = [0, 0]
        if nonbaseIdxs is None:
            nonbaseIdxs = [1, 2]

        import scipy.spatial.distance
        pools = grid.split(trialMat)
        js_dists = np.empty([len(grid.r), len(baselineIdxs)])
        for ori in range(pools.shape[0]):
            for (cond, (bcond, nbcond)) in enumerate(zip(baselineIdxs, nonbaseIdxs)):
                bp = pools[ori, bcond]
                nbp = pools[ori, nbcond]
                maxr = max(max(bp), max(nbp))
                minr = min(min(bp), min(nbp))

                bhist = np.histogram(bp, bins=10, range=[minr, maxr])[0]
                bhist = bhist / sum(bhist)
                nbhist = np.histogram(nbp, bins=10, range=[minr, maxr])[0]
                nbhist = nbhist / sum(nbhist)

                js_dists[ori, cond] = scipy.spatial.distance.jensenshannon(bhist, nbhist, base=2)
        return js_dists

    def js_distsProc(self, neuron, maxlevel=None, exci=False):
        js_dists=neuron.js_dists[:,0] if exci else neuron.js_dists
        if maxlevel is None:
            return js_dists
        elif maxlevel == 'cond':
            return np.amax(js_dists, axis=1)
        elif maxlevel == 'all':
            return np.max(js_dists)
        else:
            print('invalid maxlevel: can only be *cond* or *all*')

    def js_summary(self, neuron, stimTrace):
        return stimTrace.stimObjectsLabels('js_dist', neuron.js_distsProc(maxlevel='cond'))

class SpontBehavior:
    def stats(self,neuron):
        import scipy.stats
        dff = neuron.dff
        neuron.rmse     = np.sqrt(np.mean(dff ** 2)) #rmse
        neuron.kurtosis = scipy.stats.kurtosis(dff) #kurtosis
        neuron.stdev    = np.std(dff) #stdev
        neuron.sampen   = scipy.stats.entropy(dff) #sampen

class Foo:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)




########### load one run's data #######################################

def loadmat(filename):
    import mat73
    import scipy.io
    try:
        data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    except:
        data = mat73.loadmat(filename, use_attrdict=True)
    return data

import numpy as np


mouseDir = 'D:\\333_Galaxy Maotuan\\I love studying\\2022 winter\\lab\\papers\\2022-2-18'
mice = ['T01', 'T02', 'T03', 'Sut3', 'Sut4']
mouse = mice[0]
import glob

allfiles = glob.glob(mouseDir + '\\' + mouse + '\\*_Suite2p_*.mat')




file = allfiles[0]
data = loadmat(file)['suite2pData']
name='S005'
idx=1
################ SESSIONS ##########################
import math

freq = math.ceil(data.ops.fs)
Stim = data.Stim
sesStarts = [int(x - 1) for x in data.startIdx]
sesEnds = [int(x - 1) for x in data.endIdx]
sesPoles = sesStarts[1:]

numVisstim = Stim.numVisstim
stimSessions = np.nonzero(numVisstim)[0]
n_StimSessions = len(stimSessions)
n_sessions = len(sesStarts)

############## STIM ###################################################

orientations = np.empty(len(Stim.orientationsUsed), dtype=StimObject)
for i, ori in enumerate(Stim.orientationsUsed):
    orientations[i] = StimObject(type='ori', label=ori, index=i)

stimEvents = np.empty(len(Stim.condition), dtype=StimEvent)
for i, (on, off, type) in enumerate(zip(Stim.trialonsets, Stim.trialoffsets, Stim.condition)):
    type = type - 1  # conversion between matlab:1 and python:0
    stimEvents[i] = StimEvent(on=on, off=off, stim=orientations[type], pos_glob=i)

numVisstim = numVisstim[stimSessions]
stimPoles = [None] * (n_StimSessions - 1)
for n in range(n_StimSessions - 1):
    stimPoles[n] = sum(numVisstim[0:n + 1])

stimEvents_bySes = np.hsplit(stimEvents, stimPoles)


######################### NEURON ##############################
def getRoiMask(stat):
    med = ImagePixel(x=stat.med[0], y=stat.med[1])
    fills = np.empty(len(stat.xpix), dtype=ImagePixel)
    for ii, (x, y, o) in enumerate(zip(stat.xpix, stat.ypix, stat.overlap)):
        fills[ii] = ImagePixel(x=x, y=y, overlap=o)
    circs = np.empty(len(stat.xext), dtype=ImagePixel)
    for ii, (x, y) in enumerate(zip(stat.xext, stat.yext)):
        circs[ii] = ImagePixel(x=x, y=y)
    return RoiMask(fills=fills, circs=circs, med=med)


neurons_bySes = np.empty([data.dFF.shape[0], n_sessions], dtype=Neuron)
neuronIdxs = data.iscell[:, 0].astype(bool)
probIsCell = data.iscell[:, 1][neuronIdxs]
stats = data.stat[neuronIdxs]
for i, (dff, stat, spks, snr, cellidx, auc) in enumerate(
        zip(data.dFF, stats, data.spks, data.snr, data.cellIdx, data.AUC)):
    dff = np.hsplit(dff, sesPoles)
    spks = np.hsplit(spks, sesPoles)
    roiMask = getRoiMask(stat)
    for j in range(n_sessions):
        neurons_bySes[i, j] = Neuron(dff=dff[j], idx=i, idx_glob=cellidx, sessions=j, runs=idx, roiMask=roiMask,
                                     snr=snr,
                                     auc=auc[j], spks=spks[j])
neurons_bySes = np.swapaxes(neurons_bySes, 0, 1)

################ getting session stuff ############
sessions = np.empty(n_sessions, dtype=Session)

for i in range(n_sessions):
    sesStart = sesStarts[i]
    sesEnd = sesEnds[i]
    sessions[i] = Session(index=i,
                          neurons=neurons_bySes[i],
                          end_glob=sesEnd,
                          start_glob=sesStart,
                          type='spont')
    for neuron in sessions[i].neurons:
        neuron.sessions = sessions[i]

################# add stimTrace into ses ############
shiftOn = freq * 3
shiftOff = freq * 3
for i, stimevents_ses in enumerate(stimEvents_bySes):
    sesIdx = stimSessions[i]
    sesStart = sesStarts[sesIdx]
    sesEnd = sesEnds[sesIdx]
    for j, event in enumerate(stimevents_ses):
        event.ses = sesIdx
        event.pos_ses = j
        event.zero(sesStart)
    stimTrace = StimTrace(stimEvents=stimevents_ses,
                           end=sesEnd - sesStart,
                           stimObjects=orientations)
    stimTrace.setGrid(shiftOn, shiftOff)
    sessions[sesIdx].stimTrace = stimTrace
    sessions[sesIdx].type = 'stim'

############## create RUN ############################
# TODO: get name table for each mouse
run = Run(sessions=sessions,
          mouse=data.nidaqAligned.mouse, idx=idx, name=name, date=str(data.nidaqAligned.date))
saveRunFolder = 'D:\\333_Galaxy Maotuan\\I love studying\\2022 winter\\lab\\papers\\2022-6-28'

#################### js_dist #####################


plotFolder = 'D:\\333_Galaxy Maotuan\\I love studying\\2022 winter\\lab\\papers\\2022-5-11'
js = VisDriven(stat_name='Jensen-Shannon distance',threshold=.3)
js_dists=np.empty([len(sessions[1].neurons),n_StimSessions])
for i, i_ss in enumerate(stimSessions):
    sesInfo = str(run)+'_ses'+str(i_ss)
    st = sessions[i_ss].stimTrace
    st.setGrid(shiftOn, shiftOff)
    vis_driven_neurons=[]
    for j,neu in enumerate(sessions[i_ss].neurons):
        neu.setPlotTitle(sesInfo)
        neu.trialsMat = st.reorderByStim(st.cutouts(neu.dff))
        neu.js_dists = js.js_dist(trialMat=neu.trialsMat, grid=st.grid)
        js_dists[j,i]=js.js_distsProc(neuron=neu, maxlevel='all', exci=True)
        if js_dists[j,i]>0.3:
            neu.visdriven=True
            vis_driven_neurons.append(neu)
            '''neu.plotTrialsMat(grid=st.grid,
                              runName=str(run),
                              xlabel=VisDriven.js_summary(neuron=neu, stimTrace=st),
                              save_folder=plotFolder+'\\'+str(run))'''
        else:
            neu.visdriven=False
        sessions[i_ss].neurons_vis_driven = vis_driven_neurons

plotFolder = 'D:\\333_Galaxy Maotuan\\I love studying\\2022 winter\\lab\\papers\\2022-5-11'
spontNeurons=sessions[0].neurons
sesInfo = str(run)+'_ses0'
neu.setPlotTitle(sesInfo)
for neu in sessions[1].neurons_vis_driven:
    spontNeu = spontNeurons[neu.idx]
    spontNeu.setPlotTitle(sesInfo)
    spontNeu.plotTrace(save_folder=plotFolder,stretch=4000)

    # plot traces of vis driven neurons in spontaneous sessions
    # calculate their stats
    # put stats under session


############## all done ################
class Foo:
    test=None
    test2=None
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

import sqlite3
connection = sqlite3.connect('test.db')

foo=Foo(test='test',test2='2')
for key in vars(foo).keys():
    print(key)


# session identifier
# run identifier
'''
driven thresholds
crosscor
slide of background on what 2p imaging is
aknowledgement, spaceholders for plots
homeostatic plasticity
what a neuron looks like across all stimulations
10-15min ~12 slides


'''

import pickle

with open('test.pickle', 'ab') as file:
    pickle.dump(sessions, file)
    pickle.dump(item, file, pickle.HIGHEST_PROTOCOL)


with open('test.pickle', 'rb') as file:
    pickle.load(file)





