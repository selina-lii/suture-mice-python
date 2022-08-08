import math
import operator
import os
import numpy as np
import glob
import matplotlib.pyplot as plt
from math import log2,inf,ceil
import time
from datetime import datetime
import shutil
from scipy.special import rel_entr
from scipy.spatial.distance import jensenshannon
from EntropyHub import SampEn
from scipy.stats import kurtosis,sem
from LyxTools import *

############################# ADD DATA ##############################################

def get_name_tag(id_mouse,name_mouse,name_ses,fp_suite2pData):
    fn_suite2pData = fp_suite2pData.split('\\')[-1].split('_')

    id_mouse=id_mouse
    name_mouse=name_mouse
    cond_code = name_ses[0]
    hour_code = name_ses[1:]
    year = fn_suite2pData[1][0:2]
    month = fn_suite2pData[1][2:4]
    day = fn_suite2pData[1][4:6]

    nametag = DBinterface(id_mouse=id_mouse, name_mouse=name_mouse,
                           cond_code=cond_code, hour_code=hour_code,year=year, month=month, day=day)

    return nametag

def loop(db, config, save_refimg=False, save_dff=False, add_mice=False, add_sessions=False, add_runs=False, add_stims=False, add_neurons=False):
    for id_mouse, name_mouse in enumerate(config.name_mice):
        if not config.select_mouse or (config.select_mouse and name_mouse == config.name_mice_selected):

            print(name_mouse)
            fps_mouse_suite2pData = glob.glob("%s\\%s\\*dff*.mat" % (config.inputdir, name_mouse))

            if len(fps_mouse_suite2pData) > 0:
                n_files = config.n_mice_days[id_mouse]
                assert n_files == len(fps_mouse_suite2pData)
                name_sess = ses_names(config.gui_mice_days[id_mouse],config.conditions,
                                      config.hours)

                if add_mice:
                    conds = np.array([x[0] for x in mouse['name_sess']])
                    _, cond_poles = np.unique(conds, return_index=True)
                    cond_poles.sort()
                    conds = conds[cond_poles]

                    mouse = DBinterface(_id=id_mouse,
                                        name=name_mouse,
                                        n_sess=n_files,
                                        name_sess=name_sess,
                                        conds=conds,
                                        cond_poles=cond_poles[1:])
                    insert(mouse, db.mouse, id_mouse)

                if save_dff:
                    print(
                        'As for 7-24-22, extraction of dff was done in MATLAB and stored in C:\\Users\\selinali\\lab\\sut\\dff.')
                    pass

                read_suite2pData = save_refimg or add_sessions or add_runs or add_neurons

                for id_ses, fp_suite2pData in enumerate(fps_mouse_suite2pData):
                    if read_suite2pData:
                        data = loadmat(fp_suite2pData)['suite2pData']
                        print(name_mouse + ' ' + str(id_ses))

                        _id_ses = "%d%02d" % (id_mouse, id_ses)

                        name_ses = name_sess[id_ses]
                        name_tag=get_name_tag(id_mouse,name_mouse,name_ses,fp_suite2pData)

                        n_run = len(data.startIdx)
                        n_visstim = intarr(data.Stim.numVisstim)

                        if add_sessions:
                            fp_refimg = "%s\\refimg\\%s_%s_ref.png" % (config.workdir, name_mouse, name_ses)
                            fp_dff = "%s\\dff\\%s_%s_%s%s%s_dff.mat" % (
                                config.workdir, name_mouse, name_ses, name_tag.year, name_tag.month, name_tag.day)
                            n_neu = len(data.cellIdx)

                            run_id_spont = [];run_id_stim = []
                            for i, val in enumerate(n_visstim):
                                if val == 0:
                                    run_id_spont.append(i)
                                else:
                                    run_id_stim.append(i)

                            try:
                                magnification = float(data.config.Magnification)
                            except:
                                magnification = None
                            framerate = float(data.ops.fs)

                            ses = DBinterface(_id=_id_ses, n_run=n_run, n_neu=n_neu, run_id_stim=run_id_stim,
                                              run_id_spont=run_id_spont,
                                              magnification=magnification, framerate=framerate,
                                              fp_suite2pData=fp_suite2pData, fp_dff=fp_dff, fp_refimg=fp_refimg,
                                              id_ses=id_ses
                                              )
                            ses.properties.update(name_tag.properties)
                            insert(ses, db.session, _id_ses)

                        if save_refimg:
                            fp_refimg = "%s\\refimg\\%s_%s_ref.png" % (config.workdir, name_mouse, name_ses)
                            plt.imsave(fp_refimg, data.ops.refImg, cmap='pink')

                        if add_neurons:
                            neuronIdxs = np.where(data.iscell[:, 0].astype(bool))[0]
                            stats = [data.stat[i] for i in neuronIdxs]
                            snrs = data.snr
                            cellIdxs = intarr(data.cellIdx, mat2py=True)

                            for id_neu, (stat, snr, id_neu_glob) in enumerate(zip(stats, snrs, cellIdxs)):
                                _id_neu = "%d%02d%d%04d" % (id_mouse, id_ses, 0, id_neu)
                                snr = snr.item()
                                try:
                                    roi_pix_x = intarr(stat.xext)
                                    roi_pix_y = intarr(stat.yext)
                                except:
                                    roi_pix_x = intarr(stat.xcirc)
                                    roi_pix_y = intarr(stat.ycirc)
                                roi_med = [int(stat.med[1]), int(stat.med[0])]
                                neuron = DBinterface(_id=_id_neu, id_mouse=id_mouse, _id_ses=_id_ses, id_neu=id_neu,
                                                     snr=snr, id_neu_glob=id_neu_glob,
                                                     roi_pix_x=roi_pix_x, roi_pix_y=roi_pix_y, roi_med=roi_med)
                                insert(neuron, db.neuron, _id_neu)

                        if add_runs or add_stims:
                            #degrees_of_angle = intarr(data.Stim.orientationsUsed)
                            ons = intarr(data.Stim.trialonsets)
                            offs = intarr(data.Stim.trialoffsets)
                            stim_trace = intarr(data.Stim.condition, mat2py=True)
                            run_starts = intarr(data.startIdx, mat2py=True)
                            run_ends = intarr(data.endIdx, mat2py=True)
                            id_stim_glob = 0

                            # ref back to ses
                            stim_labels_ses, n_stim_ori_ses = np.unique(stim_trace, return_counts=True)
                            stim_labels_ses=intarr(stim_labels_ses)
                            n_stim_ori_ses=intarr(n_stim_ori_ses)
                            sf_id(db.session, _id_ses, 'stim_labels', stim_labels_ses)
                            sf_id(db.session, _id_ses, 'n_stim_ori', n_stim_ori_ses)
                            n_stim_ori_ses = [0]*len(config.stim_labels)
                            # for each run

                            for id_run in range(n_run):
                                _id_run = "%d%02d%d" % (id_mouse, id_ses, id_run)
                                type = 'spont' if n_visstim[id_run] == 0 else 'stim'

                                # config and qc
                                end_glob = run_ends[id_run]
                                start_glob = run_starts[id_run]
                                end = end_glob - start_glob

                                # files

                                # counts
                                n_neu = len(data.cellIdx)
                                run = DBinterface(_id=_id_run, type=type, end_glob=end_glob, start_glob=start_glob, end=end,
                                                  n_neu=n_neu,id_run=id_run,id_ses=_id_ses)
                                run.properties.update(name_tag.properties)

                                if type=='stim':
                                    #TODO: fp_trialsmat to be tested
                                    fp_trialsmat='%s%s\\%s_%s_%s'%(config.workdir,config.trialsmatdir,name_mouse,name_ses,id_run)

                                    n_stim_ori_run = [0]*len(config.stim_labels) #TODO: kinda a logical loophole in our data..
                                    id_stim = 0
                                    n_stim = n_visstim[id_run]
                                    n_stim_countdown = n_stim
                                    start_run_glob = run_starts[id_run]
                                    max_stim_len = 0
                                    off_prev = None
                                    iti_max = 0
                                    iti_min = run_ends[-1]  # equivalent to 'inf'

                                    while n_stim_countdown > 0:
                                        _id_stim = "%d%02d%d%03d" % (id_mouse, id_ses, id_run, id_stim)

                                        on_glob = ons[id_stim_glob]
                                        on = on_glob - start_run_glob

                                        off_glob = offs[id_stim_glob]
                                        off = off_glob - start_run_glob

                                        stim_len = off - on
                                        if stim_len > max_stim_len: max_stim_len = stim_len
                                        if off_prev is not None:
                                            gap = on - off_prev
                                            if gap > iti_max: iti_max = gap
                                            if gap < iti_min: iti_min = gap

                                        id_stim_label = stim_trace[id_stim_glob]

                                        id_stim_ori = n_stim_ori_run[id_stim_label]
                                        id_stim_ori_glob = n_stim_ori_ses[id_stim_label]

                                        if add_stims:
                                            stim_event = DBinterface(_id=_id_stim,
                                                                     on_glob=on_glob, off_glob=off_glob, on=on, off=off,
                                                                     stim_len=stim_len,
                                                                     id_stim_glob=id_stim_glob,
                                                                     id_stim_ori=id_stim_ori, id_stim_ori_glob=id_stim_ori_glob,
                                                                     id_stim_label=id_stim_label,
                                                                     id_stim=id_stim, id_run=_id_run, id_ses=_id_ses,
                                                                     id_mouse=id_mouse, name_mouse=name_mouse
                                                                     )
                                            insert(stim_event, db.stim_event, _id_stim)

                                        n_stim_countdown = n_stim_countdown - 1
                                        id_stim_glob = id_stim_glob + 1
                                        id_stim = id_stim + 1
                                        n_stim_ori_ses[id_stim_label] = n_stim_ori_ses[id_stim_label] + 1
                                        n_stim_ori_run[id_stim_label] = n_stim_ori_run[id_stim_label] + 1
                                        off_prev = off

                                    # only keeping the labels that actually appeared in this run
                                    stim_labels_run = []
                                    n_stim_ori_run_list = []
                                    stim_labels_low_n = []
                                    for label,n_stim_ori in enumerate(n_stim_ori_run):
                                        if n_stim_ori != 0:
                                            stim_labels_run.append(label)
                                            n_stim_ori_run_list.append(n_stim_ori)
                                            if n_stim_ori <= config.n_stim_ori_threshold:
                                                stim_labels_low_n.append(label)

                                    run.properties.update(n_stim = n_stim,
                                                          stim_labels=stim_labels_run,
                                                          n_stim_ori=n_stim_ori_run_list,
                                                          stim_labels_low_n=stim_labels_low_n,
                                                          max_stim_len=max_stim_len,
                                                          iti_max=iti_max,iti_min=iti_min)
                                if add_runs:
                                    insert(run, db.run, _id_run)

def add_grids(stim_labels):
    for run in db.run.find({'type':'stim'}):
        dropped_last_trial = False
        stim_len = run['max_stim_len']

        _id_run = run['_id']
        _id_grid = _id_run
        pre_on = run['iti_min'] // 2
        post_on = run['iti_min'] // 2 + run['max_stim_len']
        ons = []
        offs = []

        stims = list(db.stim_event.find({'id_run': run['_id']}).sort('on'))
        for stim in stims:
            on = stim['on'] - pre_on;ons.append(on)
            off = stim['on'] + post_on;offs.append(off)
        if off > run['end']:
            stims.pop(-1);ons.pop(-1);offs.pop(-1)
            print(_id_grid + ':' + 'Last trial exceeds total number of frames by %s frames and is flagged for exclusion' %
                (off - run['end']))
            dropped_last_trial = True

        ons_reor = []
        offs_reor = []

        for i, label in enumerate(run['stim_labels']):
            stims = list(db.stim_event.find({'$and': [{'id_run': _id_run}, {'id_stim_label': label}]}).sort("id_stim_ori"))
            for stim in stims:
                if dropped_last_trial and stim['_id']==run['n_stim']-1:
                    run['n_stim_ori'][i] = run['n_stim_ori'][i] - 1
                    run['n_stim'] = run['n_stim'] - 1
                else:
                    on = stim['on'] - pre_on;ons_reor.append(on)
                    off = stim['on'] + post_on;offs_reor.append(off)

        row_poles = np.cumsum(run['n_stim_ori'])[:-1].tolist()
        row_end = run['n_stim']
        row_labels = [stim_labels[i] for i in run['stim_labels']]
        col_poles = [pre_on, pre_on + stim_len]
        col_end = pre_on + post_on
        col_labels = ['stim_on', 'stim_off']

        grid = DBinterface(_id=_id_grid,
                           ons=ons, offs=offs, ons_reor=ons_reor, offs_reor=offs_reor,
                           row_poles=row_poles, row_end=row_end, row_labels=row_labels,
                           col_poles=col_poles, col_end=col_end, col_labels=col_labels,
                           dropped_last_trial=dropped_last_trial,
                           mouse_id=run['id_mouse'], ses_id=run['id_ses'], run_id=run['id_run'])
        insert(grid, db.grid, _id_grid)

#TODO #should i have a function for this?
def add_dependencies():
    for mouse in db.mouse.find({},{'_id':1}):
        runs=list(db.run.find({'id_ses':str(mouse['_id'])+'00'},{'start_glob':1,'end_glob':1,'run_id_stim':1,'run_id_spont':1}))

        ends=[x['end_glob'] for x in runs]
        starts=[x['start_glob'] for x in runs]
        sf_id(db.mouse,mouse['_id'],'ends',ends)
        sf_id(db.mouse,mouse['_id'],'starts',starts)

    for mouse in db.mouse.find({}, {'_id': 1}):
        ses = db.session.find_one({'id_mouse': mouse['_id']},{'run_id_stim': 1, 'run_id_spont': 1})
        run_id_stim=ses['run_id_stim']
        run_id_spont=ses['run_id_spont']
        sf_id(db.mouse, mouse['_id'], 'run_id_stim', run_id_stim)
        sf_id(db.mouse, mouse['_id'], 'run_id_spont', run_id_spont)


#TODO
def add_empty_neu_runs():
    #stim:
    for ses in db.session.find({},{'run_id_stim':1,'n_neu':1}):
        _id_ses=ses['_id']
        runs=ses['run_id_stim']
        for id_neu in range(ses['n_neu']):
            for id_run in runs:
                _id='%s%d%04d'%(_id_ses,id_run,id_neu)
                id_mouse=int(_id_ses[0])
                _id_run=_id[0:3]
                neu=DBinterface(_id=_id,id_mouse=id_mouse,_id_ses=_id_ses,_id_run=_id_run,id_run=id_run,id_neu=id_neu)
                neu.insert(db.neu_run)

def add_basic_neu_runs_spont():
    #spont:
    for ses in db.session.find({},{'run_id_spont':1,'fp_dff':1}):
        dff_ses=loadmat(ses['fp_dff'])['dFF']
        _id_ses=ses['_id']
        for id_run in ses['run_id_spont']:
            _id_run='%s%01d'%(_id_ses,id_run)
            run=db.run.find_one({'_id':_id_run},{'start_glob':1,'end_glob':1})
            dff_run=dff_ses[:,run['start_glob']:run['end_glob']]
            for id_neu,dff in enumerate(dff_run):
                _id='%s%04d'%(_id_run,id_neu)
                id_mouse=int(_id_ses[0])

                rms=float(np.sqrt(np.mean(np.square(dff))))
                kurt=kurtosis(dff)
                #samp_en=float(SampEn(dff)[0][1]) #TODO: not sure... taking m=2 rn
                mean=float(np.mean(dff))
                std=float(np.std(dff))
                neu=DBinterface(_id=_id,mean=mean,std=std,min=float(min(dff)),max=float(max(dff)),
                                rms=rms,kurt=kurt,#samp_en=samp_en,
                                id_mouse=id_mouse,_id_ses=_id_ses,_id_run=_id_run,id_run=id_run,id_neu=id_neu)
                neu.insert(db.neu_run_spont)

#TODO: this is never,never,never gonna finish at 10sec/neuron
def add_sampen():
    for ses in db.session.find({},{'run_id_spont':1,'fp_dff':1}):
        dff_ses=loadmat(ses['fp_dff'])['dFF']
        _id_ses=ses['_id']
        for id_run in ses['run_id_spont']:
            _id_run='%s%01d'%(_id_ses,id_run)
            run=db.run.find_one({'_id':_id_run},{'start_glob':1,'end_glob':1})
            dff_run=dff_ses[:,run['start_glob']:run['end_glob']]
            for id_neu,dff in enumerate(dff_run):
                _id='%s%04d'%(_id_run,id_neu)
                samp_en=float(SampEn(dff)[0][1]) #TODO: not sure... taking m=2 rn
                sf_id(db.neu_run_spont,_id,'samp_en',samp_en)

def get_sampen():
    for ses in db.session.find({'id_mouse':2,'id_ses':{'$gt':5}},{'run_id_spont':1,'fp_dff':1}):
        dff_ses=loadmat(ses['fp_dff'])['dFF']
        _id_ses=ses['_id']
        for id_run in ses['run_id_spont']:
            _id_run='%s%01d'%(_id_ses,id_run)
            run=db.run.find_one({'_id':_id_run},{'start_glob':1,'end_glob':1})
            dff_run=dff_ses[:,run['start_glob']:run['end_glob']]
            for id_neu,dff in enumerate(dff_run):
                _id = '%s%04d' % (_id_run, id_neu)
                if fd_id(db.temp_sampen,_id) is None:
                    samp_en=[list(x) for x in SampEn(dff,tau=0.2*np.std(dff))]
                    try:
                        db.temp_sampen.insert_one({'_id':_id,'samp_en':samp_en})
                    except:
                        pass

############################# DEFINE VISUALLY DRIVEN NEURONS ##############################################
def tm(dff,grid):
    trialsmat = np.empty([grid['row_end'], grid['col_end']], )
    for j, (on, off) in enumerate(zip(grid['ons_reor'], grid['offs_reor'])):
        trialsmat[j] = dff[on:off]
    return trialsmat

# kl or js
def dist_distrib(dff, grid, run, baselineIdxs, nonbaseIdxs, dist_type='kl_raw'):
    trialmat = tm(dff,grid)
    pools = np.array([np.hsplit(x, grid['col_poles']) for x in np.split(trialmat, grid['row_poles'])], dtype=object)
    for r in range(pools.shape[0]):
        for c in range(pools.shape[1]):
            pools[r][c] = pools[r][c].flatten()

    js_dists = np.empty([len(baselineIdxs), len(run['stim_labels'])])
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

            if dist_type=='kl':
                kl=[p*(log2(p/q)) if p!=0 and q!=0 else 0 for p,q in zip(nbhist,bhist)]
                js_dists[cond, ori]=max(np.sum(kl),0)

            elif dist_type=='kl_raw':
                js_dists[cond, ori] = sum(rel_entr(bhist, nbhist))

            elif dist_type=='js':
                dist_mean=[(b+nb)/2 for nb,b in zip(nbhist,bhist)]
                kl_l=[p*(log2(p/q)) if p!=0 else 0 for p,q in zip(nbhist,dist_mean)]
                kl_r=[p*(log2(p/q)) if p!=0 else 0 for p,q in zip(bhist,dist_mean)]
                js_dists[cond, ori] =np.sqrt((sum(kl_l)+sum(kl_r))/2)
                #js_dists[cond, ori] = jensenshannon(bhist, nbhist)

            else:
                KeyError('unsupported type of distribution distance')
    return js_dists

def visual_driven_loop(db,dist_type):
    baselineIdxs = [0, 0, 1]
    nonbaseIdxs = [1, 2, 2]
    #names = ['b1', 's', 'b2']
    for ses in db.session.find():
        _id_ses = ses['_id']
        dffs_ses = loadmat(ses['fp_dff'])['dFF']
        for id_run in ses['run_id_stim']:
            _id_run = _id_ses + str(id_run)
            _id_grid = _id_run
            print(_id_run)
            run = fd_id(db.run,_id_run)
            dffs_run = dffs_ses[:, run['start_glob']:run['end_glob']]
            grid = fd_id(db.grid,_id_grid)
            for neu_id, dff in enumerate(dffs_run):
                vd = dist_distrib(dff, grid, run, baselineIdxs, nonbaseIdxs,dist_type=dist_type)
                #sf_id(db.neu_run, "%s%d%04d"%(_id_ses,id_run,neu_id), dist_type, vd.tolist())

def summary_of_array(arr,keep_rows=None,keep_cols=None,max=False,binarize=False,thres=None,
                     scale_rows=None,scale_rows_factor=None,scale_cols=None,scale_cols_factor=None):
    #TODO: is it okay not to copy arr into a new var.. gonna be deep copy & change the original value
    arrmax_row=None
    arrmax_col=None
    arrmax=None

    arr_bin=None
    arrmax_row_bin=None
    arrmax_col_bin=None
    arrmax_bin=None

    if isinstance(arr,list):
        arr=np.asarray(arr)

    if keep_rows is not None:
        arr = arr[keep_rows]
    if keep_cols is not None:
        arr = arr[:, keep_cols]

    if scale_rows is not None:
        if isinstance(scale_rows_factor,list):
            assert len(scale_rows)==scale_rows_factor
            for r,rf in zip(scale_rows,scale_rows_factor):
                arr[r]=arr[r]*rf
        else:
            arr[scale_rows]=arr[scale_rows]*scale_rows_factor
    if scale_cols is not None:
        if isinstance(scale_cols_factor,list):
            assert len(scale_cols)==scale_cols_factor
            for c,cf in zip(scale_cols,scale_cols_factor):
                arr[:,c]=arr[:,c]*cf
        else:
            arr[:,scale_cols]=arr[:,scale_cols]*scale_cols_factor

    if max:
        arrmax_col=np.max(arr,axis=0)
        arrmax_row=np.max(arr,axis=1)
        arrmax = np.max(arrmax_row)

    if binarize:  # there must be a threshold
        assert thres is not None
        arr_bin = np.where(arr > thres, True, False)
        if max:
            arrmax_col_bin = np.logical_or.reduce(arr_bin, 0)
            arrmax_row_bin = np.logical_or.reduce(arr_bin, 1)
            arrmax_bin = arrmax > thres

    return arr.tolist(),arrmax_row.tolist(),arrmax_col.tolist(),float(arrmax),\
           arr_bin.tolist(),arrmax_row_bin.tolist(),arrmax_col_bin.tolist(),bool(arrmax_bin)

def vd_proc(jss,stim_labels_low_n,scale_down_factor,thres,selected_conds=None):
    return summary_of_array(jss, keep_rows=selected_conds, max=True, binarize=True, thres=thres,
                     scale_cols=stim_labels_low_n,
                     scale_cols_factor=scale_down_factor)

#have a functional system to calculate statistics of visual drivenness and put them into neu_run not neuron,
# run it
# and another function to pull out plots into a separate folder (functional filename) so we can see if the ones that are determined
# to be visually driven are actually visually driven

#TODO: put this into.. idk the right shape since we have neu_run now
def js_proc_loop(threshold,scale_down_factor):
    #selected_conds=[0,1]
    for ses in db.session.find({},{'run_id_stim':1,'stim_labels':1,'id_mouse':1}):
        print(ses['_id'])
        thres=threshold[ses['id_mouse']]
        stim_labels_low_n_runs=[]
        for id_run in ses['run_id_stim']:
            stim_labels_low_n=db.run.find_one({'_id': '%s%01d' % (ses['_id'], id_run)}, {'stim_labels_low_n': 1})['stim_labels_low_n']
            stim_labels_low_n=[ses['stim_labels'].index(i) for i in stim_labels_low_n]
            stim_labels_low_n_runs.append(stim_labels_low_n)
        for neu in db.neu_run.find({'id_ses':ses['_id']},{'js':1}):
            for id_run,stim_labels_low_n in zip(ses['run_id_stim'],stim_labels_low_n_runs):
                vd_corrected,vd_max_cond,vd_max_ori,vd_max,is_visdriven_mat,is_visdriven_cond,is_visdriven_ori,is_visdriven=\
                    vd_proc(neu['js'],stim_labels_low_n,scale_down_factor,thres)#selected_conds
                sf_id(db.neu_run, neu['_id'], 'js_corrected', js_corrected)
                sf_id(db.neu_run, neu['_id'], 'js_max_cond', js_max_cond)
                sf_id(db.neu_run, neu['_id'], 'js_max_ori', js_max_ori)
                sf_id(db.neu_run, neu['_id'], 'js_max', js_max)
                sf_id(db.neu_run, neu['_id'], 'is_visdriven_mat', is_visdriven_mat)
                sf_id(db.neu_run, neu['_id'], 'is_visdriven_cond', is_visdriven_cond)
                sf_id(db.neu_run, neu['_id'], 'is_visdriven_ori', is_visdriven_ori)
                sf_id(db.neu_run, neu['_id'], 'is_visdriven', is_visdriven)

############################# PLOTTING V1.0 ############################################
def plot_mouse_stat(data,title,xticklabels,filename,err=None):
    plt.clf()
    plt.figure(figsize=(len(data), 7.5))
    l = plt.plot(data)[0]
    if err is not None:
        plt.fill_between(l.get_xdata(),data+err,data-err,color=l.get_color(),alpha=0.2)
    plt.title(title)
    ax = plt.gca()
    ax.set_xticks(list(range(len(data))))
    ax.set_xticklabels(xticklabels, fontsize=8)
    plt.savefig(filename, bbox_inches='tight', pad_inches=0.2)
    plt.close()

def n_visdriven(db):
    n_vd_neu=np.zeros([2,7,30],dtype=int)
    for mouse in db.mouse.find():
        id_mouse = mouse['_id']
        for ses in db.session.find({'id_mouse':id_mouse},{'run_id_stim':1,'id_ses':1}):
            for i,id_run in enumerate(ses['run_id_stim']):
                n_vd_neu[i,id_mouse,ses['id_ses']]=len(list(db.neuron.find({'id_ses':ses['_id'],'is_visdriven_%d'%id_run:True},{'is_visdriven_%d':1})))

    for i,_ in enumerate(n_vd_neu):
        for j,each in enumerate(_):
            mouse=fd_id(db.mouse,j)
            if mouse is not None and len(mouse['run_id_stim'])>i:
                mouse_name = mouse['name']
                id_run = mouse['run_id_stim'][i]

                n=mouse['n_sess']
                names=mouse['name_sess']
                d=each[:n]
                tt="# visually driven neurons - %s run%d" % (mouse_name, id_run)
                fn='%s_run%d_nvdn.jpg' % (mouse_name, id_run)
                plot_mouse_stat(d,n,tt,names,fn)
    return n_vd_neu


def average_response_spontaneous(db):
    avgact_vd_neu_in_spont=np.zeros([2,7,30])

    for mouse in db.mouse.find():
        id_mouse=mouse['_id']
        for j in range(mouse['n_sess']):
            ses = db.session.find_one({'_id': '%d%02d' % (mouse['_id'],j)}, {'run_id_stim': 1,'fp_dff' : 1})
            neus = list(db.neuron.find({'id_ses':ses['_id']},{'_id':0,'is_visdriven_1':1,'is_visdriven_2':1,'is_visdriven_3':1}))
            dff = loadmat(ses['fp_dff'])['dFF']
            for i,id_run in enumerate(ses['run_id_stim']):
                run = db.run.find_one({'_id':'%s%d'%(ses['_id'],id_run)},{'start_glob':1,'end_glob':1})
                ids=[]
                for k,neu in enumerate(neus):
                    if neu['is_visdriven_%d'%id_run] is True:
                        ids.append(k)
                if len(ids)!=0:
                    dff_run=dff[ids,run['start_glob']:run['end_glob']]
                    avgact_vd_neu_in_spont[i,id_mouse,j]=np.mean(dff_run)

    for mouse in db.mouse.find():
        n = mouse['n_sess']
        m_n = mouse['name']
        names = mouse['name_sess']
        for i,id_run in enumerate(mouse['run_id_stim']):
            tt = "average response of driven neurons in spontaneous run - %s run%d" % (m_n, id_run)
            fn = 'avgact_%s_run%d.jpg' % (m_n, id_run)
            d = avgact_vd_neu_in_spont[i,mouse['_id'],:n]
            plot_mouse_stat(d, n, tt, names, fn)
    return avgact_vd_neu_in_spont

#assigning crossday index for each neuron
def crossday(db,crossdaydir):
    for mouse in db.mouse.find():
        id_mouse=mouse['_id']
        run_id_spont=mouse['run_id_spont']
        run_id_stim = mouse['run_id_stim']
        fp_crossday=crossdaydir+'\\'+mouse['name']+'.mat'
        try:
            ids_cd=loadmat(fp_crossday)
            ids_cd=ids_cd.get(list(ids_cd.keys())[-1])
        except:
            continue
        for id_cdneu,neu_cd in enumerate(ids_cd):
            neu_cd=intarr(neu_cd,mat2py=True)
            for id_ses,id_neu in enumerate(neu_cd):
                if id_neu!=-1:
                    #_id_neu = '%d%02d0%04d' % (id_mouse, id_ses, id_neu)
                    for id_run in run_id_spont:
                        _id_neu = '%d%02d%d%04d' % (id_mouse, id_ses, id_run, id_neu)
                        sf_id(db.neu_run_spont, _id_neu, 'id_cd', id_cdneu)
                    #print('%s: %s'%(_id_neu,id_cdneu) )
                    #sf_id(db.neuron,_id_neu,'id_cd',id_cdneu)

def crossday_imshow():
    for mouse in db.mouse.find():
        id_mouse=mouse['_id']
        fp_crossday=config.workdir+config.crossdaydir+'\\'+mouse['name']+'.mat'
        try:
            ids_cd=loadmat(fp_crossday)
            ids_cd=ids_cd.get(list(ids_cd.keys())[-1])
        except:
            continue
        plt.clf()
        n=ids_cd.shape[0]//100+1
        nx=20
        ny=n//nx+1
        fig, axs = plt.subplots(ny,nx,figsize=(nx*ids_cd.shape[1]//12, ny*6))
        fig.suptitle('Horizontally stacked subplots')
        for ax in axs.flatten():
            ax.set_visible(False)
        for ax,slice in zip(axs.flatten(),np.array_split((ids_cd>0)*1,n)):
            ax.imshow(slice)
            ax.set_visible(True)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
        fig.suptitle('%s crossday cell consistency' % (mouse['name']))
        plt.tight_layout()
        savePlot_fig(fig, 'cd check_%d.jpg'%id_mouse)

#TODO: this is still a script not a function
def crossday_meanact():
    stat='mean' #'rms'
    for mouse in db.mouse.find({'_id':{'$gt':4}}):
        id_mouse=mouse['_id']
        print(mouse['name'])
        n_sess=mouse['n_sess']
        #number of days active to be a valid crossday neuron
        thres_valid=0
        try:
            ids_cd=loadmat(mouse['fp_crossday'])
            ids_cd=ids_cd.get(list(ids_cd.keys())[-1])
            if thres_valid==0:
                valid_ids=list(range(ids_cd.shape[0]))
            else:
                valid_ids=np.where(np.sum(ids_cd>0,axis=1)>thres_valid)[0].tolist()
        except:
            continue
        for id_run in mouse['run_id_spont']:
            print('run%d'%id_run)
            stats = np.empty([n_sess, len(valid_ids)], dtype=float)
            stats[:] = np.nan
            for i,id_cd in enumerate(valid_ids):
                sess=[]
                stats_cd=[]
                findquery = {'id_mouse': id_mouse, 'id_run': id_run, 'id_cd': id_cd, 'is_nonphys': {'$exists': 0}}
                projection = {'_id': 0, stat: 1, 'id_ses': 1}
                for neu in db.neu_run_spont.find(findquery,projection):
                    sess.append(neu['id_ses'])
                    stats_cd.append(neu[stat])
                stats[sess,i]=stats_cd
                if id_cd%1000==0:
                    print(id_cd)

            plt.clf()
            fig = plt.figure(figsize=(20,15))
            plt.plot(stats, alpha=0.4)
            plt.title('%s run %d crossday %s' % (mouse['name'],id_run,stat))
            plt.tight_layout()
            savePlot(fig, 'c\\cd %s run %d mouse %d.jpg'%(stat,id_run,id_mouse))

def add_neurons_crossday():
    for mouse in db.mouse.find():
        try:
            ids_cd=loadmat(mouse['fp_crossday'])
            ids_cd=ids_cd.get(list(ids_cd.keys())[-1])
        except:
            continue
        print(mouse['name'])
        id_mouse=mouse['_id']
        n_sess=mouse['n_sess']
        cond_poles=mouse['cond_poles']

        for id_cdneu, id_neus in enumerate(ids_cd):
            _id='%d%04d'%(id_mouse,id_cdneu)

            vec=(id_neus!=0)
            is_on_n_sess=np.count_nonzero(vec)
            is_on_percent_sess=is_on_n_sess/n_sess*100
            is_on_all_conds=np.all([np.any(x) for x in np.split(vec,cond_poles)]).tolist()
            id_neus=id_neus.tolist()

            DBinterface(_id=_id,id_mouse=id_mouse,id_cdneu=id_cdneu,id_neus=id_neus,
                        is_on_all_conds=is_on_all_conds,
                        is_on_n_sess=is_on_n_sess,
                        is_on_percent_sess=is_on_percent_sess).insert(db.neu_cd)


def crossday_mean_acitivty_subselect():
    stats=['mean','std','rms']
    projection = {'_id': 0, 'id_ses': 1}
    for stat in stats:
        projection[stat]=1

    for mouse in db.mouse.find():
        print(mouse['name'])
        id_mouse=mouse['_id']
        n_sess = mouse['n_sess']
        for i_run,id_run in enumerate(mouse['run_id_spont']):
            print('run%d' % id_run)
            stats_mean = np.empty([len(stats),n_sess])
            stats_err = np.empty([len(stats),n_sess])
            neus_cd=list(db.neu_cd.find({'is_on_all_conds':True,'id_mouse':id_mouse},{'_id':0,'id_cdneu':1}))
            valid_ids=[x['id_cdneu'] for x in neus_cd]

            stat_result = np.empty([len(stats), n_sess, len(valid_ids)], dtype=float)
            stat_result[:] = np.nan
            for i,id_cd in enumerate(valid_ids):
                findquery = {'id_mouse': id_mouse, 'id_run': id_run, 'id_cd': id_cd, 'is_nonphys': {'$exists': 0}}
                for neu in db.neu_run_spont.find(findquery, projection):
                    for j,stat in enumerate(stats):
                        stat_result[j, neu['id_ses'], i] = neu[stat]
            stats_mean[:,:]=np.nanmean(stat_result,axis=2)
            stats_err[:,:]=sem(stat_result,axis=2,nan_policy='omit')
            print(np.count_nonzero(~np.isnan(stat_result[0]),axis=1))
            if id_mouse==0 and id_run==2:
                stats_mean[:,4]=np.nan #because the run is not right rn
                stats_err[:,4]=np.nan

            for j, stat in enumerate(stats):
                plt.clf()
                fig = plt.figure(figsize=(12*(n_sess/16),9))
                mean=stats_mean[j]
                err=stats_err[j]
                l = plt.plot(mean)[0]
                plt.fill_between(l.get_xdata(), mean - err, mean + err, color=l.get_color(),
                                 alpha=0.2, edgecolor='none')
                plt.title('%s crossday %s for spontaneous run#%d \n(subselect neurons that are active for at least 1 day on all all conditons)' % (mouse['name'], stat, i_run))
                plt.tight_layout()
                savePlot(fig, '%s\\cd %s run %d mouse %d subselect.jpg'%(stat,stat,id_run,id_mouse))


def crossday_mean_acitivty_across_thresholds():
    thresholds=[0.5,0.6,0.7,0.8,0.9]
    stats=['mean','std','rms']
    projection = {'_id': 0, 'id_ses': 1}
    for stat in stats:
        projection[stat]=1

    for mouse in db.mouse.find():
        print(mouse['name'])
        id_mouse=mouse['_id']
        n_sess=mouse['n_sess']
        try:
            ids_cd=loadmat(mouse['fp_crossday'])
            ids_cd=ids_cd.get(list(ids_cd.keys())[-1])
        except:
            continue
        for i_run,id_run in enumerate(mouse['run_id_spont']):
            print('run%d' % id_run)
            stats_mean = np.empty([len(stats),len(thresholds),n_sess])
            stats_err = np.empty([len(stats),len(thresholds),n_sess])
            n_valid=[]
            for k,thres in enumerate(thresholds):
                valid_ids=np.where(np.sum(ids_cd>0,axis=1)>thres*n_sess)[0].tolist()
                n_valid.append(len(valid_ids))
                stat_result = np.empty([len(stats), n_sess, len(valid_ids)], dtype=float)
                stat_result[:] = np.nan
                for i,id_cd in enumerate(valid_ids):
                    findquery = {'id_mouse': id_mouse, 'id_run': id_run, 'id_cd': id_cd, 'is_nonphys': {'$exists': 0}}
                    for neu in db.neu_run_spont.find(findquery, projection):
                        for j,stat in enumerate(stats):
                            stat_result[j, neu['id_ses'], i] = neu[stat]
                stats_mean[:,k,:]=np.nanmean(stat_result,axis=2)
                stats_err[:,k,:]=sem(stat_result,axis=2,nan_policy='omit')
                print(np.count_nonzero(~np.isnan(stat_result[0]),axis=1))
            if id_mouse==0 and id_run==2:
                stats_mean[:,:,4]=np.nan #because the run is not right rn
                stats_err[:,:,4]=np.nan

            for j, stat in enumerate(stats):
                plt.clf()
                fig = plt.figure(figsize=(12*(n_sess/16),9))
                for mean,err in zip(stats_mean[j],stats_err[j]):
                    l = plt.plot(mean)[0]
                    plt.fill_between(l.get_xdata(), mean - err, mean + err, color=l.get_color(),
                                     alpha=0.2, edgecolor='none')
                plt.legend(['>%d%% (n=%s)' % (t*100,n) for t,n in zip(thresholds,n_valid)],fontsize=8)
                plt.title('%s crossday %s for spontaneous run#%d' % (mouse['name'], stat, i_run))
                plt.tight_layout()
                savePlot(fig, '%s\\cd %s run %d mouse %d.jpg'%(stat,stat,id_run,id_mouse))


'''def find_crossday(db,crossdaydir):
    for mouse in db.mouse.find():
        for
        db.neuron.find({'id_mouse'},{'_id':1})
'''
def plot_trialsmat_loop(db,config,id_mouse_sel=None):
    if id_mouse_sel is not None:
        query = {'id_mouse':id_mouse_sel}
    else:
        query = {}
    for ses in db.session.find(query,{'fp_dff':1,'run_id_stim':1,'n_neu':1,'id_mouse':1}):
        id_mouse=ses['id_mouse']
        asp = config.plt_tm_cb_asp[id_mouse]
        pad = config.plt_tm_cb_pad[id_mouse]
        dff_all=loadmat(ses['fp_dff'])['dFF']
        for id_run in ses['run_id_stim']:
            _id_run = ses['_id']+str(id_run)
            run = db.run.find_one({'_id':_id_run},{'start_glob':1,'end_glob':1,'fp_trialsmat':1})
            fpout = run['fp_trialsmat']
            if not os.path.isdir(fpout):
                os.mkdir(fpout)
            grid=db.grid.find_one({'_id':_id_run},{'ons':0,'offs':0})
            dff_r=dff_all[:,run['start_glob']:run['end_glob']]
            for id_neu,dff in enumerate(dff_r):
                _id_neu='%s%04d'%(_id_run, id_neu)
                plotTrialsMat(dff, grid, fpout, _id_neu, asp, pad)

def plotTrialsMat(dff, grid, folder, _id_neu, cb_asp, cb_pad):
    trialsmat = tm(dff, grid)
    h = trialsmat.shape[0]
    w = trialsmat.shape[1]

    plt.clf()
    fig = plt.figure(figsize=(w/16, h/8))

    plt.imshow(trialsmat, cmap=plt.get_cmap('turbo'))
    plt.colorbar(aspect=cb_asp,shrink=0.5,orientation='horizontal',pad=cb_pad)
    plt.margins(0)
    plt.gca().set_aspect(3)

    # gridlines
    for (pos,gridLabel) in zip(grid['col_poles'],grid['col_labels']):
        plt.axvline(pos+0.5, c='w', linestyle='--')
        plt.text(pos-10, 2, gridLabel, c='w')
    grid['row_poles'].append(grid['row_end'])
    for (pos, gridlabel) in zip(grid['row_poles'],grid['row_labels']):
        plt.axhline(pos+0.5, c='w', linestyle='--')
        plt.text(0, pos-0.5, str(gridlabel)+'°', c='w')

    plt.tight_layout()
    savePlot(plt.gca(), folder +'\\' + _id_neu + '.jpg')

#TODO: average response

#TODO: might be a module
def db_of_plots():
    fp_in = 'C:\\Users\\selinali\\lab\\RDE20141\\2022-7-7-plots\\trace\\'
    fp_db=config.workdir + '\\2022-7-22-plots\\db_trace\\'
    for file in glob.glob(fp_in + '\\*\\*.jpg'):
        shutil.move(file, fp_db + os.path.basename(file))

def test_threshold_for_js():
    for mouse in db.mouse.find({}, {'_id': 1, 'run_id_stim': 1, 'name': 1, 'js_m+1std': 1}):
        id_mouse = mouse['_id']
        for id_run, thres in zip(mouse['run_id_stim'], mouse['js_m+1std']):
            fp_out = config.workdir + '\\2022-7-22-plots\\test_tm2\\%s_run%d' % (mouse['name'], id_run)
            mkdir(fp_out)
            for type in ['above_thres', 'below_thres']:
                fp_out_type = fp_out + '\\' + type
                mkdir(fp_out_type)
                if type == 'above_thres':
                    query = {'id_mouse': id_mouse, 'id_run': id_run, 'js_max': {'$gt': thres}, 'is_visdriven': False}
                else:
                    query = {'id_mouse': id_mouse, 'id_run': id_run, 'js_max': {'$lt': thres, '$gt': thres - .05},
                             'is_visdriven': False}
                ll = list(db.neu_run.find(query, {'_id': 1, 'js_max': 1}).sort('js_max', -1))
                for i, neu in enumerate(ll):
                    shutil.copy2(config.workdir + '\\2022-7-22-plots\\db_trialsmat\\' + neu['_id'] + '.jpg',
                                 fp_out_type + '\\%d_%s_%.2f.jpg' % (i, neu['_id'], neu['js_max']))

def add_AUC():
    for ses in db.session.find({'id_mouse':6,'id_ses':{'$gt':3}}, {'fp_suite2pData': 1}):
        print(ses['_id'])
        data = loadmat(ses['fp_suite2pData'])['suite2pData']
        data = data.AUC
        for id_neu, auc in enumerate(data):
            _id='%s0%04d'%(ses['_id'],id_neu)
            sf_id(db.neuron,_id,'auc',auc.tolist())


def histograms_of_stim(fp=None):
    if fp is None:
        fp='C:\\Users\\selinali\lab\\sut\\2022-7-22-plots\\stim_stats'
    stats = ['js_max']
    for run in db.run.find({'type':'stim'},{'name_mouse':1,'cond_code':1,'hour_code':1}):
        _id_run=run['_id']
        ll=list(db.neu_run.find({'_id_run':_id_run},{'js_max':1}))
        for stat in stats:
            mkdir(fp+'\\'+stat)
            ll_s=[x[stat] for x in ll]
            mean=np.mean(ll_s)
            std=np.std(ll_s)

            plt.clf()
            plt.hist(ll_s,bins=15)
            plt.axvline(mean+std,color='r')
            plt.title('%s %s%s histogram of %s for stim run%s'%(run['name_mouse'],run['cond_code'],run['hour_code'],stat,_id_run[-1]))
            plt.xlabel('mean=%.6f, std=%.6f' % (mean,std))
            fn='%s_%s_%s_%s_%s%s'%(stat,_id_run[-1],run['name_mouse'],_id_run[1:3],run['cond_code'],run['hour_code'])
            plt.tight_layout()
            savePlot(plt.gca(), fp+'\\'+stat+'\\'+fn + '.jpg')

types=['b',     'l',    'u',    'u',    'u',    'u']
thres_u=[0.6,   10,     inf,    2,      200,    3]
thres_l=[0,     0,      -2,     -inf,   -inf,   -inf]
stats = ['mean','max', 'min',   'std',  'kurt', 'rms']
def histograms_of_spont(fp=None):
    if fp is None:
        fp='C:\\Users\\selinali\lab\\sut\\2022-7-22-plots\\spont_stats'
    for run in db.run.find({'type':'spont'},{'name_mouse':1,'cond_code':1,'hour_code':1}):
        _id_run=run['_id']
        ll=list(db.neu_run_spont.find({'_id_run':_id_run}))
        for stat,out_u,out_l in zip(stats,thres_u,thres_l):
            mkdir(fp+'\\'+stat)
            ll_s=[x[stat] for x in ll]
            ll_s_phys=[]
            outs=[]
            cnt=0
            for x in ll_s:
                if x<out_l or x>out_u:
                    cnt=cnt+1
                    x=int(x)
                    if x==0 or x==-1 or x==1:
                        pass
                    else:
                        outs.append(int(x))
                else:
                    ll_s_phys.append(x)

            mean = np.mean(ll_s_phys)
            std = np.std(ll_s_phys)
            median = np.median(ll_s_phys)

            plt.clf()
            plt.hist(ll_s_phys,bins=20)
            plt.axvline(mean + std, color='r')
            plt.axvline(mean, color='y')
            plt.axvline(mean - std, color='r')
            plt.axvline(median,color='k')
            plt.text(median,0,'median',color='k',rotation=-90)
            plt.text(mean, 0, 'mean', color='y',rotation=-90)
            plt.text(median+std, 0, '+1std', color='r',rotation=-90)
            plt.text(median-std, 0, '-1std', color='r',rotation=-90)
            tex='mean=%.6f\nmedian=%.6f:\n+1std=%.6f\n-1std=%.6f\n'%(mean,median,mean+std,mean-std)
            plt.text(plt.gca().get_xlim()[1]*0.6,plt.gca().get_ylim()[1]*0.6,tex)
            plt.title('%s %s%s histogram of %s for spontaneous run%s'%(run['name_mouse'],run['cond_code'],run['hour_code'],stat,_id_run[-1]))
            plt.xlabel('#outliers%d:%s'%(cnt,outs))
            plt.tight_layout()
            fn='%s_%s_%s_%s_%s%s'%(stat,_id_run[-1],run['name_mouse'],_id_run[1:3],run['cond_code'],run['hour_code'])
            savePlot(plt.gca(), fp+'\\'+stat+'\\'+fn + '.jpg')


def mark_nonphysiological_neurons():
    '''    thres_u = [10, 20, inf, 2, 200, 3]
        thres_l = [-2, 0, -2, -inf, -inf, -inf]
    '''

    thres_u = [0.5, 100, inf, 10, 500, 10]
    thres_l = [0, -5, -10, -inf, -inf, -inf]
    stats = ['mean', 'max', 'min', 'std', 'kurt', 'rms']
    for neu in db.neu_run_spont.find({}, {'mean': 1, 'std': 1, 'min': 1, 'max': 1, 'rms': 1, 'kurt': 1}):
        for stat, out_u, out_l in zip(stats, thres_u, thres_l):
            x = neu[stat]
            if x < out_l or x > out_u:
                sf_id(db.neu_run_spont, neu['_id'], 'is_nonphys', True)


def get_keys_collection(col):
    return col.aggregate([
        {"$project": {"arrayofkeyvalue": {"$objectToArray": "$$ROOT"}}},
        {"$unwind": "$arrayofkeyvalue"},
        {"$group": {"_id": None, "allkeys": {"$addToSet": "$arrayofkeyvalue.k"}}}
    ]).next()["allkeys"]

def add_indexs_collection(col):
    start_time = time.time()
    keys=get_keys_collection(col)
    print(keys)
    for k in keys:
        col.create_index([(k, 1)])
    end_time = time.time()
    print('time elapsed:%.2f'%(end_time - start_time))

def add_auc_test_query_pattern():
    #comparing runtime between two search schema...
    start_time = time.time()
    #projection={'_id':0,'auc':1,'_id_ses':1,'id_neu':1}
    projection={'auc':1}
    for neu in db.neuron.find({},projection):
        _id=neu['_id']
        for id_run,auc in enumerate(neu['auc']):
            findquery={'_id':'%s%d%s' %(_id[:3],id_run,_id[4:])}
                #time elapsed:170.34 / 3min
            #findquery={'_id_run':'%s%d' %(neu['_id_ses'],id_run),'id_neu':neu['id_neu']}
                #time elapsed:479.65 / 8min
            db.neu_run2.update_one(findquery,{'$set':{'auc':auc}})
    end_time = time.time()
    print('time elapsed:%.2f'%(end_time - start_time))

def std_of_auc_across_days():
    stats=['auc']
    projection = {'_id': 0}
    for stat in stats:
        projection[stat]=1
    stat_all=[]
    for mouse in db.mouse.find({'_id':{'$gt':0}}):
        mouse_name = mouse['name']
        name_sess = mouse['name_sess']
        n_sess = mouse['n_sess']
        id_mouse = mouse['_id']
        print(mouse['name'])
        for id_run in range(mouse['n_runs']):
            stats_mean = np.empty([len(stats), n_sess], dtype=float)
            stats_err = np.empty([len(stats), n_sess], dtype=float)
            for ses in db.session.find({'id_mouse':id_mouse},{'id_ses':1}):
                id_ses=ses['id_ses']

                if id_mouse == 0 and id_run>1 and id_ses == 4:
                    stats_mean[:, id_ses] = np.nan
                    stats_err[:, id_ses] = np.nan
                    continue

                neus=list(db.neu_run2.find({'_id_run': '%d%02d%d' % (id_mouse,id_ses, id_run),'is_nonphys':{'$exists':0}}, projection))
                for s,stat in enumerate(stats):
                    stat_vals=[x[stat] for x in neus]
                    stats_mean[s,id_ses]=np.mean(stat_vals)
                    stats_err[s,id_ses]=sem(stat_vals)

            for s,stat in enumerate(stats):
                tt="average %s for %s run %d" % (stat,mouse_name, id_run)
                fn='%s %s run%d.jpg' % (stat,mouse_name, id_run)
                plot_mouse_stat(stats_mean[s],tt,name_sess,fn,err=stats_err[s])
    return stat_all

def average_response():
    pass

start_time = time.time()
print(datetime.now())

db = pymongo.MongoClient("mongodb://localhost:27017/").sut_mice2

# create collections
#db['mouse'];db['session'];db['run'];db['neuron'];db['stim_event'];db['grid'];db['backup'];

# macro info
#CHANGE workdir, inputdir AND EVERYTHING ELSE YOU NEED TO CHANGE TO GET IT WORKING :D
config = get_config(db)

# add all data by looping through the .mat files once
#loop(db,config,save_refimg=False,save_dff=False,add_mice=False,add_sessions=False,add_runs=True,add_stims=False,add_neurons=False)

# add grids for js and stuff
#add_grids(config.stim_labels)
#add_empty_neu_runs()

#visual_driven_loop(db)
#n_vd_neu = js_proc_loop(config.js_thresholds,config.js_scale_down_factor)

#avgact_vd_neu_in_spont = average_response_spontaneous(db)

#crossday(db,config.workdir+config.crossdaydir)

#plot_trialsmat_loop(db,config)

#visual_driven_loop(db,'js')

#add_AUC()

#add_basic_neu_runs_spont()

#will take  forever to run ;(
#get_sampen()
#histograms_of_spont()
#histograms_of_stim()

#mark_nonphysiological_neurons()
#crossday_meanact()
#crossday_mean_acitivty_across_thresholds()
#std_of_auc_across_days()

crossday_mean_acitivty_subselect()

end_time = time.time()
print(datetime.now())
print('time elapsed:%.2f'%(end_time - start_time))
assert 0==1
############################### scripts

for i,neu in enumerate(list(db.neu_run.find({'id_mouse':0,'js_max':{'$gt':0.3},'is_visdriven':False},{'_id':1}).sort('js_max',-1))):
    shutil.copy2(config.workdir+'\\2022-7-22-plots\\db_trialsmat\\'+neu['_id']+'.jpg',
                 config.workdir+'\\2022-7-22-plots\\test_tm_subthres\\%d_%s.jpg'%(i,neu['_id']))

#0.38606758684929166
for neu in db.neu_run.find({},{'_id':1}):
    #sf_id(db.neu_run,neu['_id'],'_id_ses',neu['id_ses'])
    sf_id(db.neu_run,neu['_id'],'id_run',int(int(neu['_id'][3])))

for mouse in db.mouse.find({},{'_id':1,'run_id_stim':1}):
    mm=[]
    sd=[]
    m1sd=[]
    for id_run in mouse['run_id_stim']:
        ll = list(db.neu_run.find({'id_mouse': mouse['_id'], 'id_run': id_run}, {'js_max': 1}))
        ll_js=[x['js_max'] for x in ll]
        mm.append(np.mean(ll_js))
        sd.append(np.std(ll_js))
    m1sd=[x+y for x,y in zip(mm,sd)]
    sf_id(db.mouse,mouse['_id'],'js_mean',mm)
    sf_id(db.mouse, mouse['_id'], 'js_std', sd)
    sf_id(db.mouse, mouse['_id'], 'js_m+1std', m1sd)

