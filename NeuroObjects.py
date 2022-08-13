
class StimTrace():
    def stimBlocks(self,ax):
        ymin,ymax=ax.get_ylim()
        for event in self.stimEvents:
            ax.plot([event.on,event.off], [ymax,ymax], c=event.stim.color, marker='|',linewidth=20)
            ax.vlines([event.on,event.off], ymin, ymax,colors=['k','k'],linestyles='dashed',alpha=0.3)
        return ax

    def trace(self, highlight_oris=None,scale=1):
        trace = np.zeros(self.end)
        if highlight_oris is None:
            for event in self.stimEvents:
                trace[event.on:event.off] = scale
        else:
            for event in self.stimEvents:
                    trace[event.on:event.off] = scale*2 if event.stim in highlight_oris else scale
        return trace


class Neuron():
    # statistics
    visdriven = None
    js_dists = None
    active_oris = None

    plot_label = None

    filename_trace = None
    filename_roi = None
    filename_trialsMat = None

    def plotTitle(self):
        # 'T03_1_B04_ses2 #3'
        return str(self.session)+' #'+str(self.idx)

    def plotFilename(self,folder,plotType,includeSession=True):
        # 'T03_1_B04_ses2_3_trialsMat.jpg'
        if not includeSession:
            return folder+'\\'+str(self.run)+'_#'+str(self.idx)+'_'+plotType+'.jpg'
        return folder+'\\'+str(self.session)+'_#'+str(self.idx)+'_'+plotType+'.jpg'

    def plotFolder(self,parentFolder,plotType,includeSession=True):
        # 'C:\\some_path\\plot_folder\\T03_1_B04_ses2_trialsMat'
        if not includeSession:
            return parentFolder+'\\'+str(self.run)+'_'+plotType
        return parentFolder+'\\'+str(self.session)+'_'+plotType

    def plotTrace(self,parentFolder,stretch=2500):
        plotType='trace'
        filename=self.precheck(parentFolder,plotType)
        if filename is None:
            return
        eps=0.01
        ymin=round(min(self.dff)-eps,2)
        ymax=round(max(self.dff)+eps,2)
        mar=stretch/4e5
        n_splits=len(self.dff)//stretch+1
        w=50
        h=(n_splits+1)*stretch//1000
        n_tks=round(stretch/100)
        isStim=self.session.isStim()

        splits=np.array_split(self.dff,n_splits)

        if isStim:
            st = self.session.stimTrace
            splits_stim = np.array_split(st.trace(scale=ymax/3,
                                                  highlight_oris=None),#self.active_oris
                                         n_splits)
        plt.clf()
        fig = plt.figure(figsize=(w, h))

        def frame2sec(tk,shift):
            return '%d'%(tk+shift) +'('+ '%.2f'%((tk+shift)/self.session.freq) +'s)'

        shift=0
        for i,split in enumerate(splits):
            ax = plt.subplot(n_splits+1, 1, i+1)
            ax.margins(mar)
            ax.set_ylim(ymin=ymin,ymax=ymax)
            plt.plot(split)
            tks = np.arange(0, len(split), 200)
            ax.set_xticks(tks)
            ax.set_xticklabels([frame2sec(tk,shift) for tk in tks])
            if isStim:
                plt.plot(splits_stim[i],c='k',alpha=0.3)
            shift = shift + len(split)

        ax = plt.subplot(n_splits+1, 1, i+2)
        ax.margins(mar)
        plt.plot(self.dff,color='k')
        if isStim:
            ax = self.session.stimTrace.stimBlocks(ax=ax)

        plt.title(self.plotTitle(), fontsize=16)
        plt.tight_layout()
        self.savePlot(plt.gca(), filename,plotType)

    def statsTag(self,xlabel):
        plotType='trialsmat'
        filename=self.precheck(parentFolder,plotType)
        xlabel=    something
        plt.text(10, -2*xlabel.count('\n'), xlabel)
        if filename is None:
            plt.title(self.plotTitle())

    def plotAvgResponse(self,parentFolder):
        pass

    def plotROI(self,parentFolder):
        plotType='roi'
        filename=self.precheck(parentFolder,plotType,includeSession=False)
        if filename is None:
            return
        ax = self.roiMask.plot()
        ax.set_title(self.plotTitle())
        self.savePlot(ax, filename,plotType)

    def precheck(self,parentFolder,plotType,includeSession=True):
        save_folder = self.plotFolder(parentFolder=parentFolder,plotType=plotType,includeSession=includeSession)
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        filename=self.plotFilename(folder=save_folder,plotType=plotType,includeSession=includeSession)
        return filename if not os.path.exists(filename) else None

class RoiMask():
    circs=None
    box=None
    med=None
    h=None
    w=None

    def __init__(self, **kwargs):
        margin=2
        self.box=[[min(xs)-margin,max(xs)+margin],
                  [min(ys)-margin,max(ys)+margin]]
        self.h,self.w=self.refImg.shape
        self.h=self.h/70
        self.w=self.w/70

    def plot(self):
        plt.clf()
        plt.figure(figsize=(self.w, self.h))
        plt.imshow(self.refImg, cmap='pink')
        flpink='#FF1493'
        flyellow='#67DCD5'
        plt.vlines(x=self.box[0],ymin=self.box[1][0],ymax=self.box[1][1],colors=[flpink,flpink],linewidth=0.5)
        plt.hlines(y=self.box[1],xmin=self.box[0][0],xmax=self.box[0][1],colors=[flpink,flpink],linewidth=0.5)
        plt.plot([pix.x for pix in self.circs],[pix.y for pix in self.circs],c=flyellow,alpha=0.1)
        plt.plot(self.med.x, self.med.y, c='r',marker='o',markersize=0.5)
        return plt.gca()


class Run:
    #config and qc
    dprime=None
    bad_day=None
    paradigm=None # or type, e.g.long term learning / suture-resuture
    magnification=None #order of imaging
    framerate=None
    notes=None #e.g. 'magnification on'

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def name(self):
        return self.cond + str(self.day)

def js_distsProc(neuron, maxlevel=None, exci=False):
    js_dists=neuron.js_dists[:,0] if exci else neuron.js_dists
    if maxlevel is None:
        return js_dists
    elif maxlevel == 'cond':
        return np.amax(js_dists, axis=1)
    elif maxlevel == 'all':
        return np.max(js_dists)
    else:
        print('invalid maxlevel: can only be *cond* or *all*')

def js_summary(neuron, stimTrace):
    return stimTrace.stimObjectsLabels('js_dist', js_distsProc(neuron=neuron, maxlevel='cond'))
