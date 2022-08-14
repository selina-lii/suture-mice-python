import os
import numpy as np
import glob

from math import log2,inf,ceil
import time
from datetime import datetime
import shutil
from scipy.special import rel_entr
from scipy.spatial.distance import jensenshannon
from EntropyHub import SampEn
from scipy.stats import kurtosis,sem
from LyxTools import *
from plotting_functions import *
from subprocess import call

############################# ADD DATA ##############################################
def initiate_database():
    # create collections
    collections=['mouse','session','run','neuron','stim_event','grid','neu_run']
    for col in collections:
        db[col]
    set_config()

# TODO will be a GUI
def set_config():
    mice_names = ['T01', 'T02', 'T03', 'Sut1', 'Sut2', 'Sut3', 'Sut4']
    hours = ['0025', '005', '01', '04', '06', '08', '12', '24', '48', '72', '96', '120', '144', '168']
    conditions = ['B', 'S', 'U', 'R']

    dir_input = 'C:\\Users\\selinali\\lab\\Mice'
    dir_work = 'C:\\Users\\selinali\\lab\\sut'
    dir_tm = dir_work+ '\\2022-7-22-plots\\trialsmat'

    stim_labels = [0, 45, 90, 135, 180, 225, 270, 315]
    stim_colors = ['#AC041B', '#F8693D', '#EDB120', '#21C40F', '#67DCD5', '#2A46FA', '#A81CBF', '#FF1493']
    n_stim_ori_threshold = 6

    config = DBinterface(_id=0,
                         mice_names=mice_names,
                         hours=hours,
                         conditions=conditions,
                         dir_input=dir_input,
                         dir_work=dir_work,
                         dir_tm=dir_tm,
                         stim_labels=stim_labels,
                         stim_colors=stim_colors,
                         n_stim_ori_threshold=n_stim_ori_threshold).insert(db.config)
    return config

def add_mice():
    for _id, name in enumerate(config.mice_names):
        print(name)
        fps_mouse_suite2pData = glob.glob("%s\\%s\\*dff*.mat" % (config.dir_input, name))

        if len(fps_mouse_suite2pData) > 0:
            n_sess = config.n_mice_days[_id]
            assert n_sess == len(fps_mouse_suite2pData)
            DBinterface(_id=_id,name=name,n_sess=n_sess).insert(db.mouse)

    call("gui_mice_days.py", shell=True)

    for mouse in db.mouse.find():
        conds = np.array([x[0] for x in mouse['name_sess']])
        _, cond_poles = np.unique(conds, return_index=True)
        cond_poles.sort()
        conds = conds[cond_poles]
        sf_id(db.mouse,mouse['_id'],'conds',conds)
        sf_id(db.mouse,mouse['_id'],'cond_poles',cond_poles[1:])

def add_neurons():
    for mouse in db.mouse.find():
        id_mouse=mouse['_id']
        n_runs=mouse['n_runs']
        for ses in db.session.find({'id_mouse':id_mouse}):
            _id_ses=ses['_id']
            stats = loadmat(ses['fp_stats'])['stats']
            data = loadmat(ses['fp_s2p_useful'])
            dff = loadmat(ses['fp_dff'])['dff']

            id_neurons = np.where(data.iscell[:, 0].astype(bool))[0]
            stats = [stats[i] for i in id_neurons]
            snrs = data.snr
            aucs = data.AUC
            cellIdxs = intarr(data.cellIdx, mat2py=True)

            for id_neu, (stat, snr, id_neu_glob,auc_neu) in enumerate(zip(stats, snrs, cellIdxs, aucs)):
                _id_neu = "%s0%04d" % (_id_ses, id_neu)
                snr = snr.item()
                try:
                    roi_pix_x = intarr(stat.xext)
                    roi_pix_y = intarr(stat.yext)
                except:
                    roi_pix_x = intarr(stat.xcirc)
                    roi_pix_y = intarr(stat.ycirc)
                roi_med = [int(stat.med[1]), int(stat.med[0])]
                DBinterface(_id=_id_neu, id_mouse=id_mouse, _id_ses=_id_ses, id_neu=id_neu, id_neu_glob=id_neu_glob,
                            roi_pix_x=roi_pix_x, roi_pix_y=roi_pix_y, roi_med=roi_med, snr=snr).insert(db.neuron)

                for id_run,auc in enumerate(auc_neu):
                    pass





#TODO still messy
def mainloop(save_refimg=False, add_sessions=False, add_runs=False, add_stims=False, update_mouse=False):
    for mouse in db.mouse.find():
        name_mouse=mouse['name']
        id_mouse=mouse['_id']
        name_sess=mouse['name_sess']
        fps_mouse_suite2pData = glob.glob("%s\\%s\\*dff*.mat" % (config.dir_input, name_mouse))

        for id_ses, (name_ses,fp_suite2pData) in enumerate(zip(name_sess,fps_mouse_suite2pData)):
                print('%d %s'%(id_ses,name_ses))
                _id_ses = "%d%02d" % (id_mouse, id_ses)

                name_tag=get_name_tag(id_mouse,name_mouse,name_ses,fp_suite2pData)
                data = loadmat(fp_suite2pData)['suite2pData']

                n_run = len(data.startIdx)
                n_visstim = intarr(data.Stim.numVisstim)
                run_starts = intarr(data.startIdx, mat2py=True)
                run_ends = intarr(data.endIdx, mat2py=True)
                run_id_spont = np.where(n_visstim == 0)[0]
                run_id_stim = np.nonzero(n_visstim)

                if update_mouse and id_ses==0:
                    d=dict(run_starts=run_starts,run_ends=run_ends,run_id_spont=run_id_spont,run_id_stim=run_id_stim)
                    sf_id(db.mouse,id_mouse,d)

                if save_refimg:
                    fp_refimg = "%s\\refimg\\%s_%s_ref.png" % (config.workdir, name_mouse, name_ses)
                    plt.imsave(fp_refimg, data.ops.refImg, cmap='pink')

                if add_sessions:
                    n_neu = len(data.cellIdx)
                    framerate = float(data.ops.fs)

                    fp_s2p_useful = "%s\\s2p_useful\\%s_%s_%s_s2p_useful.mat" % (
                        config.workdir, name_mouse, name_ses, name_tag.date_of_acq)
                    fp_refimg = "%s\\refimg\\%s_%s_ref.png" % (config.workdir, name_mouse, name_ses)
                    fp_dff = "%s\\dff\\%s_%s_%s_dff.mat" % (
                        config.workdir, name_mouse, name_ses, name_tag.date_of_acq)
                    try:
                        magnification = float(data.config.Magnification)
                    except:
                        magnification = None #TODO some files were missing magnification?

                    stim_trace = intarr(data.Stim.condition, mat2py=True)
                    stim_labels, n_stim_ori = np.unique(stim_trace, return_counts=True)

                    ses = DBinterface(_id=_id_ses, n_run=n_run, n_neu=n_neu,
                                      run_id_stim=run_id_stim,  run_id_spont=run_id_spont,
                                      magnification=magnification, framerate=framerate,
                                      fp_suite2pData=fp_suite2pData, fp_s2p_useful=fp_s2p_useful,
                                      fp_dff=fp_dff, fp_refimg=fp_refimg,
                                      stim_labels=intarr(stim_labels), n_stim_ori=intarr(n_stim_ori),
                                      id_ses=id_ses, id_mouse=id_mouse)
                    ses.properties.update(name_tag.properties)
                    ses.insert(db.session)

                if add_runs:
                    for id_run in run_id_spont:
                        _id_run = "%d%02d%d" % (id_mouse, id_ses, id_run)
                        end = run_ends[id_run]
                        start = run_starts[id_run]
                        end = end - start
                        run = DBinterface(_id=_id_run, type='spont', end_ses=end, start_ses=start, end=end,
                                          id_run=id_run,_id_ses=_id_ses)
                        run.properties.update(name_tag.properties)
                        run.insert(db.run)

                if add_runs or add_stims:
                    poles=np.cumsum(n_visstim[run_id_stim])[:-1]
                    ons_runs = np.split(data.Stim.trialonsets,poles)
                    offs_runs = np.split(data.Stim.trialoffsets,poles)
                    stim_trace_runs = np.split(intarr(data.Stim.condition, mat2py=True),poles)
                    id_stim_ses=0

                    for id_run,ons,offs,stim_trace in zip(run_id_stim,ons_runs,offs_runs,stim_trace_runs):
                        _id_run = "%d%02d%d" % (id_mouse, id_ses, id_run)
                        end = run_ends[id_run]
                        start = run_starts[id_run]
                        length = end - start

                        n_stim_ori = [0]*len(config.stim_labels) #TODO: kinda a logical loophole in our data..
                        n_stim = n_visstim[id_run]
                        start_run_ses = run_starts[id_run]
                        max_stim_len = 0
                        off_prev = None
                        iti_max = 0
                        iti_min = run_ends[-1]  # equivalent to 'inf'

                        for id_stim,on_ses,off_ses,label in enumerate(zip(ons,offs,stim_trace)):
                            _id_stim = "%d%02d%d%03d" % (id_mouse, id_ses, id_run, id_stim)
                            on = on_ses - start_run_ses
                            off = off_ses - start_run_ses
                            stim_len = off - on

                            if stim_len > max_stim_len: max_stim_len = stim_len
                            if off_prev is not None:
                                gap = on - off_prev
                                if gap > iti_max: iti_max = gap
                                if gap < iti_min: iti_min = gap

                            n_stim_ori[label] += 1

                            if add_stims:
                                DBinterface(_id=_id_stim,
                                         on_ses=on_ses, off_ses=off_ses, on=on, off=off,
                                         stim_len=stim_len, label=label,
                                         id_stim=id_stim, id_stim_ses=id_stim_ses, _id_run=_id_run, _id_ses=_id_ses,
                                         id_mouse=id_mouse).insert(db.stim_event)

                            off_prev = off
                            id_stim_ses+=1

                        if add_runs:
                            # only keeping the labels that actually appeared in this run
                            thres=config.n_stim_ori_threshold
                            stim_labels = [label for label, n in enumerate(n_stim_ori) if n!=0]
                            n_stim_ori = np.asarray(n_stim_ori)[stim_labels].tolist()
                            stim_labels_low_n = [label for label, n in zip(stim_labels,n_stim_ori) if n<thres]

                            DBinterface(_id=_id_run, type='stim', id_run=id_run,id_ses=_id_ses,
                                        end=end, start=start, length=length,
                                        n_stim = n_stim, n_stim_ori=n_stim_ori, max_stim_len=max_stim_len,
                                        stim_labels=stim_labels, stim_labels_low_n=stim_labels_low_n,
                                        label_first_stim=stim_trace[0], label_last_stim=stim_trace[end],
                                        iti_max=iti_max, iti_min=iti_min).insert(db.run)

#dropping first and last trial
def add_grids(stim_labels):
    for run in db.run.find({'type':'stim'}):
        _id_run = run['_id']
        stim_len = run['max_stim_len']
        pre_on = run['iti_min'] // 2
        post_on = run['iti_min'] // 2 + run['max_stim_len']

        run['n_stim'] -= 2
        run['n_stim_ori'][run['label_first_stim']] -= 1
        run['n_stim_ori'][run['label_last_stim']] -= 1

        query= {'_id_run': _id_run, 'is_edge':{'$exists': 0}}

        stims = list(db.stim_event.find(query).sort('on'))
        ons=[x['on'] - pre_on for x in stims]
        offs=[x['on'] + post_on for x in stims]

        reorder = [t[0] for t in sorted(enumerate(stims), key=lambda k: (k[1]['id_stim_label'], k[1]['on']))]
        ons_reor=np.asarray(ons)[reorder].tolist()
        offs_reor=np.asarray(offs)[reorder].tolist()

        row_poles = np.cumsum(run['n_stim_ori'])[:-1].tolist()
        row_end = run['n_stim']
        row_labels = [stim_labels[i] for i in run['stim_labels']]

        col_poles = [pre_on, pre_on + stim_len]
        col_end = pre_on + post_on
        col_labels = ['stim_on', 'stim_off']

        DBinterface(_id=_id_run,
                   ons=ons, offs=offs, ons_reor=ons_reor, offs_reor=offs_reor,
                   row_poles=row_poles, row_end=row_end, row_labels=row_labels,
                   col_poles=col_poles, col_end=col_end, col_labels=col_labels).insert(db.grid)

#TODO
def add_empty_neu_runs():
    for ses in db.session.find({},{'run_id_stim':1,'run_id_spont':1,'n_neu':1}):
        _id_ses=ses['_id']
        runs=ses['run_id_stim']
        for id_neu in range(ses['n_neu']):
            for id_run in runs:
                _id='%s%d%04d'%(_id_ses,id_run,id_neu)
                id_mouse=int(_id_ses[0])
                _id_run=_id[0:3]
                neu=DBinterface(_id=_id,id_mouse=id_mouse,_id_ses=_id_ses,_id_run=_id_run,id_run=id_run,id_neu=id_neu)
                neu.insert(db.neu_run)

def add_basic_neu_runs():
    for ses in db.session.find({},{'run_id_spont':1,'fp_dff':1}):
        dff_ses=loadmat(ses['fp_dff'])['dFF']
        _id_ses=ses['_id']
        for id_run in range(ses['n_runs']):
            _id_run='%s%01d'%(_id_ses,id_run)
            run=db.run.find_one({'_id':_id_run},{'start_ses':1,'end_ses':1})
            dff_run=dff_ses[:,run['start_ses']:run['end_ses']]
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
            run=db.run.find_one({'_id':_id_run},{'start_ses':1,'end_ses':1})
            dff_run=dff_ses[:,run['start_ses']:run['end_ses']]
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
            run=db.run.find_one({'_id':_id_run},{'start_ses':1,'end_ses':1})
            dff_run=dff_ses[:,run['start_ses']:run['end_ses']]
            for id_neu,dff in enumerate(dff_run):
                _id = '%s%04d' % (_id_run, id_neu)
                if fd_id(db.temp_sampen,_id) is None:
                    samp_en=[list(x) for x in SampEn(dff,tau=0.2*np.std(dff))]
                    try:
                        db.temp_sampen.insert_one({'_id':_id,'samp_en':samp_en})
                    except:
                        pass

############################# DEFINE VISUALLY DRIVEN NEURONS ##############################################

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
    for ses in db.session.find():
        _id_ses = ses['_id']
        dffs_ses = loadmat(ses['fp_dff'])['dFF']
        for id_run in ses['run_id_stim']:
            _id_run = _id_ses + str(id_run)
            _id_grid = _id_run
            print(_id_run)
            run = fd_id(db.run,_id_run)
            dffs_run = dffs_ses[:, run['start_ses']:run['end_ses']]
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



def add_auc():
    start_time = time.time()
    projection={'auc':1}
    for neu in db.neuron.find({},projection):
        _id=neu['_id']
        for id_run,auc in enumerate(neu['auc']):
            findquery={'_id':'%s%d%s' %(_id[:3],id_run,_id[4:])}
            db.neu_run2.update_one(findquery,{'$set':{'auc':auc}})
    end_time = time.time()
    print('time elapsed:%.2f'%(end_time - start_time))


def average_response():
    pass

def cd_pull_files():
    meta='C:\\Users\\selinali\\lab\\sut\\2022-7-22-plots\\trialsmat_cd'
    for mouse in db.mouse.find({'n_neu_cd':{'$exists':1}}):
        id_mouse = mouse['_id']
        name=mouse['name']
        n_neu_cd=mouse['n_neu_cd']
        mkdir('%s\\%s'%(meta,name))
        for id_run in mouse['run_id_stim']:
            mkdir('%s\\%s\\%d'%(meta,name,id_run))
            for id_cdneu in range(n_neu_cd):
                tgt_folder='%s\\%s\\%d\\%d'%(meta,name,id_run,id_cdneu)
                mkdir(tgt_folder)
                id_neus=db.neu_cd.find_one({'id_mouse':id_mouse,'id_cdneu':id_cdneu},{'id_neus':1})['id_neus']
                for id_ses,id_neu in enumerate(id_neus):
                    if id_neu!=0:
                        id='%d%02d%d%04d'%(id_mouse,id_ses,id_run,id_neu-1)
                        src_folder='C:\\Users\\selinali\\lab\\sut\\2022-7-22-plots\\db_trialsmat'
                        try:
                            shutil.copyfile('%s\\%s.jpg'%(src_folder,id),'%s\\%s.jpg'%(tgt_folder,id))
                        except:
                            pass

def find_cd_dff_lims():
    for mouse in db.mouse.find({'keep':True}):
        print(mouse['name'])
        id_mouse = mouse['_id']
        n_runs = mouse['n_runs']
        for neu_cd in db.neu_cd.find({'id_mouse':id_mouse,'is_empty':{'$exists':0}},{'_id':0,'id_cdneu':1}):
            id_cd=neu_cd['id_cdneu']
            print(id_cd)
            for id_run in range(n_runs):
                neus = list(db.neu_run2.find({'id_mouse': id_mouse, 'id_run': id_run, 'id_cd': id_cd,'is_nonphys':{'$exists':0}},
                                             {'max': 1, 'min': 1}))
                if len(neus)!=0:
                    ymax=max([x['max'] for x in neus])
                    ymin=min([x['min'] for x in neus])
                    _ids=[x['_id'] for x in neus]
                    sf_ids(db.neu_run2,_ids,dict(dff_lims_cd=[ymin,ymax]))

def set_id_cd_for_all_neu_runs():
    id_mouse=5
    for neu_cd in db.neu_cd.find({'id_mouse':id_mouse}, {'_id': 0, 'id_cdneu': 1, 'id_neus': 1}):
        id_cd = neu_cd['id_cdneu']
        if id_cd%500==0:
            print(id_cd)
        for id_ses, id_neu in enumerate(neu_cd['id_neus']):
            if id_neu != 0:
                id_neu -= 1
                db.neu_run2.update_many(
                    {'id_mouse': id_mouse, 'id_ses': id_ses, 'id_neu': id_neu},
                    {'$set': {'id_cd': id_cd}})

start_time = time.time()
print(datetime.now())

db = pymongo.MongoClient("mongodb://localhost:27017/").sut_mice2
config = get_config(db)

#cd_visdriven_on_last_baseline(db)
#cd_pull_files()
# find_cd_dff_lims()
plot_trace_loop(config.workdir+config.db_trace,db,2)

end_time = time.time()
print(datetime.now())
print('time elapsed:%.2f'%(end_time - start_time))