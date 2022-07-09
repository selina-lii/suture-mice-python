import numpy as np
import pickle
import math
import matplotlib.pyplot as plt
import os

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

    def __init__(self, stimEvents=None, end=None, stimObjects=None):
        self.stimEvents = stimEvents
        self.end = end
        self.stimObjects = stimObjects

        self.stimLen = max([x.stimLen() for x in self.stimEvents])
        self.minTrialLen = min([x1.on - x.on for (x1, x) in zip(self.stimEvents[1:], self.stimEvents[:-1])])

        self.reorder = None
        self.grid = None
        self.trialLen = None

    def trace(self):
        trace = np.zeros(self.end)
        fig = plt.figure(figsize=(self.end/500,3))
        for event in self.stimEvents:
            trace[event.on:event.off] = 1
        plt.clf()
        plt.plot(trace)
        plt.show()
        #plt.figure().close()
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
        for i,(obj, val) in enumerate(zip(self.stimObjects, values)):
            label = label + str(obj) + ' = ' + "%.2f" % val + ';'
            if (i+1)%4==0:
                label = label + '\n'
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

        if self.trialLen > self.minTrialLen:
            print('At least one trial window overlaps previous/next trial by %s frames with current shiftOff = %s. '
                  'Consider decreasing the number of frames required for each trial.' % (str(self.trialLen - self.minTrialLen), str(shiftOff)))

        self.reorder = np.argsort(a=self.stimEvents, kind='stable')


class Session:  # or neuron??
    '''convert between suite2pdata and internal modules'''

    def __init__(self, idx=None, name=None, type=None, stimTrace=None, neurons=None, end_glob=None, start_glob=None, end=None,
                 neurons_vis_driven=None,refImg=None,run=None):
        self.idx = idx
        self.name = name
        self.type = type
        self.stimTrace = stimTrace
        self.end_glob = end_glob
        self.start_glob = start_glob
        self.neurons = neurons
        self.neurons_vis_driven = neurons_vis_driven
        self.refImg = refImg
        self.run = run
        if end is not None:
            self.end = end
        else:
            self.end = end_glob - start_glob

    def __str__(self):
        return str(self.run)+'_ses'+str(self.idx)

    def isStim(self):
        return self.type=='stim'

class Neuron:
    dff = None
    idx = None
    idx_glob = None
    session = None
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

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __str__(self):
        return str(self.idx_glob)

    def plotTitle(self):
        # 'T03_1_B04_ses2 #3'
        return str(self.session)+' #'+str(self.idx)

    def plotFilename(self,folder,plotType):
        # 'T03_1_B04_ses2_3_trialsMat.jpg'
        return folder+'//'+str(self.session)+'_#'+str(self.idx)+'_'+plotType+'.jpg'

    def plotFolder(self,parentFolder,plotType):
        # 'C://some_path//plot_folder//T03_1_B04_ses2_trialsMat'
        return parentFolder+'//'+str(self.session)+'_'+plotType

    def plotTrace(self,parentFolder,stretch=2500):
        import matplotlib.pyplot as plt
        m=stretch/4e5
        n_splits=len(self.dff)//stretch+1
        splits=np.array_split(self.dff,n_splits)
        plt.clf()
        fig = plt.figure(figsize=(50, (n_splits+1)*stretch//1000))
        plt.margins(x=0)
        plt.title(self.plotTitle(),fontsize=16)
        for i,split in enumerate(splits):
            ax = plt.subplot(n_splits+1, 1, i+1)
            ax.margins(m)
            plt.tight_layout()
            plt.plot(split)
        ax = plt.subplot(n_splits+1, 1, i+2)
        ax.margins(m)
        plt.tight_layout()
        plt.plot(self.dff,color='k')
        save_folder = self.plotFolder(parentFolder=parentFolder,plotType='trace')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        plt.savefig(self.plotFilename(folder=save_folder, plotType='trace'))
        #plt.savefig(self.plotFilename(folder='C:\\Users\\selinali\\PycharmProjects\\scientificProject',plotType='trace'))
        plt.close()

    def plotTrialsMat(self,grid,xlabel,parentFolder):
        import matplotlib.pyplot as plt

        plt.clf()
        height = self.trialsMat.shape[0]
        width = self.trialsMat.shape[1]
        fig = plt.figure(figsize=(height/6, width/30))

        plt.imshow(self.trialsMat, cmap=plt.get_cmap('turbo'))
        plt.colorbar(shrink=0.6)
        plt.margins(0)

        ASPECTRATIO=5
        plt.gca().set_aspect(height / width * ASPECTRATIO)

        # gridlines
        for (pos,gridLabel) in zip(grid.c,grid.clabels):
            plt.axvline(pos+0.5, c='w', linestyle='--')
            plt.text(pos-12, 2, gridLabel,c='w')
        for (pos, gridlabel) in zip(grid.r,grid.rlabels):
            plt.axhline(pos+0.5, c='w', linestyle='--')
            plt.text(0, pos-0.5, gridlabel, c='w')

        plt.title(self.plotTitle())
        plt.text(0, 0, xlabel)
        plt.tight_layout()

        save_folder = self.plotFolder(parentFolder=parentFolder,plotType='trialsMat')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        plt.savefig(self.plotFilename(folder=save_folder,plotType='trialsmat'),bbox_inches='tight',pad_inches=0.2)
        # plt.savefig(self.plotFilename(folder='C:\\Users\\selinali\\PycharmProjects\\scientificProject',plotType='trialsmat'),bbox_inches='tight',pad_inches=0.2)

        plt.close()

    def plotAvgResponse(self):
        # TODO
        return

    def plotROI(self):
        ref_img = self.session.refImg
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
    n_stimSes=None

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        self.n_stimSes=len([ses.type=='stim' for ses in self.sessions])

    def __str__(self):
        return "%s_%s_%s_%s" % (self.mouse, str(self.idx).zfill(2), self.name, self.date)

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
        return stimTrace.stimObjectsLabels('js_dist', self.js_distsProc(neuron=neuron, maxlevel='cond'))


class SpontBehavior:
    def stats(self,neuron):
        import scipy.stats
        dff = neuron.dff
        neuron.rmse     = np.sqrt(np.mean(dff ** 2)) #rmse
        neuron.kurtosis = scipy.stats.kurtosis(dff) #kurtosis
        neuron.stdev    = np.std(dff) #stdev
        neuron.sampen   = scipy.stats.entropy(dff) #sampen


def intarr(arr,mat2py=False):
    if mat2py:
        return [int(x-1) for x in arr]
    return [int(x) for x in arr]

def pixiedust(data,name,idx,cut_sessions=False):

    ################ SESSIONS/get fields from data ##########################

    freq = math.ceil(data.ops.fs)
    Stim = data.Stim
    sesStarts = intarr(data.startIdx,mat2py=True)
    sesEnds = intarr(data.endIdx,mat2py=True)
    sesPoles = sesStarts[1:]
    n_visstim = intarr(Stim.numVisstim)
    n_sessions = len(sesStarts)

    trialOnsets = intarr(Stim.trialonsets)
    trialoffsets = intarr(Stim.trialoffsets)
    condition =  intarr(Stim.condition,mat2py=True)
    refImg = data.ops.refImg
    oris_bySes = Stim.orientations
    oris_used = intarr(Stim.orientationsUsed)
    ################ getting session stuff ############
    sessions = np.empty(n_sessions, dtype=Session)

    for i in range(n_sessions):
        sesStart = sesStarts[i]
        sesEnd = sesEnds[i]
        type = 'spont' if n_visstim[i]==0 else 'stim'
        sessions[i] = Session(idx=i,
                              end_glob=sesEnd,
                              start_glob=sesStart,
                              type=type,
                              refImg=refImg)

    ############## STIM ###################################################
    orientations = np.empty(len(Stim.orientationsUsed), dtype=StimObject)
    for i, ori in enumerate(oris_used):
        orientations[i] = StimObject(type='ori', label=ori, index=i)

    stimEvents = np.empty(len(Stim.condition), dtype=StimEvent)
    for i, (on, off, type) in enumerate(zip(trialOnsets, trialoffsets, condition)):
        stimEvents[i] = StimEvent(on=on, off=off, stim=orientations[type], pos_glob=i)

    ################# add stimTrace into ses ############
    shiftOn = freq * 3
    shiftOff = freq * 3
    cur = -1
    poles = np.unique(np.cumsum(n_visstim))

    def findStimObjectsSes(oris_ses,orientations):
        if len(oris_ses)==len(orientations):
            return orientations
        else:
            oris_list = []
            for ori in orientations:
                if ori.label in oris_ses:
                    oris_list.append(ori)
            return np.array(oris_list)

    for i, (session,oris_ses) in enumerate(zip(sessions, oris_bySes)):
        if session.isStim():
            cur = cur + 1
            stimEvents_ses = stimEvents[poles[cur]:poles[cur+1]]
            for j, event in enumerate(stimEvents_ses):
                event.ses = i
                event.pos_ses = j
                event.zero(session.start_glob)
            oris_ses = [int(x[0]) for x in oris_ses] #TODO for some reason????
            oris_objs_ses = findStimObjectsSes(oris_ses,orientations)
            stimTrace = StimTrace(stimEvents=stimEvents_ses,
                                   end=session.end,
                                   stimObjects=oris_objs_ses)
            stimTrace.setGrid(shiftOn, shiftOff)
            session.stimTrace = stimTrace

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
    stats = [data.stat[i] for i in neuronIdxs]
    for i, (dff, stat, spks, snr, cellidx, auc) in enumerate(
            zip(data.dFF, stats, data.spks, data.snr, data.cellIdx, data.AUC)):
        dff = np.hsplit(dff, sesPoles)
        spks = np.hsplit(spks, sesPoles)
        roiMask = getRoiMask(stat)
        for j in range(n_sessions):
            neurons_bySes[i, j] = Neuron(dff=dff[j], idx=i, idx_glob=int(cellidx), sessions=j, runs=idx, roiMask=roiMask,
                                         snr=snr,
                                         auc=auc[j], spks=spks[j])
    neurons_bySes = np.swapaxes(neurons_bySes, 0, 1)

    ################ getting session stuff ############
    for i,session in enumerate(sessions):
        session.neurons = neurons_bySes[i]
        for neuron in session.neurons:
            neuron.session = session
            # or i?

    ############## create RUN ############################
    # TODO: get name table for each mouse
    run = Run(sessions=sessions, mouse=data.nidaqAligned.mouse, idx=idx, name=name, date=str(int(data.nidaqAligned.date)))
    saveRunFolder = 'D:\\333_Galaxy Maotuan\\I love studying\\2022 winter\\lab\\papers\\2022-6-28'
    for session in sessions:
        session.run=run

    #################### js_dist #####################
    js = VisDriven(stat_name='Jensen-Shannon distance',threshold=.3)
    js_dists=np.empty([len(sessions[0].neurons),run.n_stimSes])
    for i, session in enumerate(sessions):
        if session.isStim():
            st = session.stimTrace
            vis_driven_neurons=[]
            for j,neu in enumerate(session.neurons):
                neu.trialsMat = st.reorderByStim(st.cutouts(neu.dff))
                neu.js_dists = js.js_dist(trialMat=neu.trialsMat, grid=st.grid)
                js_dists[j,i]=js.js_distsProc(neuron=neu, maxlevel='all', exci=True)
                if js_dists[j,i]>0.3:
                    neu.visdriven=True
                    vis_driven_neurons.append(neu)
                else:
                    neu.visdriven=False
                session.neurons_vis_driven = vis_driven_neurons

    plotFolder = 'C:\\Users\\selinali\\lab\\RDE20141\\2022-7-7-plots'

    for i, session in enumerate(sessions[0:2] if cut_sessions else sessions):
        if session.isStim():
            st = session.stimTrace
            for neu in session.neurons:
                neu.plotTrialsMat(grid=st.grid,
                                  xlabel=js.js_summary(neuron=neu, stimTrace=st),
                                  parentFolder=plotFolder)
        for neu in session.neurons:
            neu.plotTrace(parentFolder=plotFolder)

    filename=str(run)+'.pickle'
    with open(filename, 'wb') as file:
        pickle.dump(run, file)

#TODO:
#   1. show more info from jsdistsummary (which condition was selected)...??
#   2. for plot trace:
#       1. overlay faded stimTrace in the background and peak active orientations like we did in matlab, if it's a stim session
#       2. trace ticks into seconds, define stretch based on frame rate
#   5. generate all trials mats and traces for all neurons
#   8. positions of neurons across days
#       2. make a function that draws ROI masks
#       3. loop that with plottrialsmat and plottrace for all neurons
#   3. put all neurons into cross day then pull out ROI (and all other plots) for each 'neuron across day'
#       write a function that select a subset of plots based on input

#   7. ofc we can run clustering on spontaneous neurons but
#   8. in the long term - how can we get videos of each neuron into this data by segmenting the original suite2p video
#   9. i wanna see if there's any form of wave propagating or if there's any wave like patterns like the ones i saw in whole brain fmri
#   10. write the same function that pulls out all the files, but with an actual photo viewer oh wow *sparkles*
#   11. so we'll have the whole 2p video loaded
#   then if you click on some neuron it's gonna show a bubble with its index, stats and. a randomly generated name (lol fr
#   and a button that says 'plot all properties' or drop down list 'plot __' like trials mat or trace..
#   then another window pops up that gives you the plot on demand
#   if you select two neurons you can see their correlation
#   or compare plots in a panel called 'compare activity'
#   but ofc eventually there'll be a network plotted on all these neurons based on their strength of connectivity
#   for the plotting, neurons would need to have colors given by mapping principles (adjacent neurons should have different colors)
#   it'll be best to connect this to suite2p but that might be harder..
#   lol i'm basically here dreaming up an entire lego castle.
#   12. baseline level of activity for two populations - spont and driven
#   pull out all their traces, average and get mean & std
#   13. honestly I wanna generate soundwaves for neuronal activity cause that might work better than looking a 15,000 frames long line...
#   and like if you wanna compare just have two tones overlayed to each other

    plotFolder = 'D:\\333_Galaxy Maotuan\\I love studying\\2022 winter\\lab\\papers\\2022-5-11'
    spontNeurons=sessions[0].neurons
    sesInfo = str(run)+'_ses0'
    for neu in sessions[1].neurons_vis_driven:
        spontNeu = spontNeurons[neu.idx]
        spontNeu.setPlotTitle(sesInfo)
        spontNeu.plotTrace(save_folder=plotFolder,stretch=4000)

        # plot traces of vis driven neurons in spontaneous sessions
        # calculate their stats
        # put stats under session

    ############## all done ################

    '''
    # session identifier
    # run identifier

    driven thresholds
    crosscor
    slide of background on what 2p imaging is
    aknowledgement, spaceholders for plots
    homeostatic plasticity
    what a neuron looks like across all stimulations
    10-15min ~12 slides
    '''
    
#crossdayidx = loadmat('C:\\Users\\selinali\\lab\\RDE20141\\Sut4CellIndex.mat')['cellindexalign']
import glob
allfiles = glob.glob('Sut4' + '*.pickle')
for filename in allfiles:
    with open('Sut1_0_B005_190924.pickle', 'rb') as file:
        x=pickle.load(file)


from scipy.io.wavfile import write
write('test.wav', 1563, data.spks[0,:])