import os.path

import matplotlib.pyplot as plt
from math import inf
from LyxTools import loadmat,savePlot,mkdir,savePlot_fig
import numpy as np
from scipy.stats import sem

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

def average_response_spontaneous(db):
    avgact_vd_neu_in_spont=np.zeros([2,7,30])

    for mouse in db.mouse.find():
        id_mouse=mouse['_id']
        for j in range(mouse['n_sess']):
            ses = db.session.find_one({'_id': '%d%02d' % (mouse['_id'],j)}, {'run_id_stim': 1,'fp_dff' : 1})
            neus = list(db.neuron.find({'id_ses':ses['_id']},{'_id':0,'is_visdriven_1':1,'is_visdriven_2':1,'is_visdriven_3':1}))
            dff = loadmat(ses['fp_dff'])['dFF']
            for i,id_run in enumerate(ses['run_id_stim']):
                run = db.run.find_one({'_id':'%s%d'%(ses['_id'],id_run)},{'start_ses':1,'end_ses':1})
                ids=[]
                for k,neu in enumerate(neus):
                    if neu['is_visdriven_%d'%id_run] is True:
                        ids.append(k)
                if len(ids)!=0:
                    dff_run=dff[ids,run['start_ses']:run['end_ses']]
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
            run = db.run.find_one({'_id':_id_run},{'start_ses':1,'end_ses':1,'fp_trialsmat':1})
            fpout = run['fp_trialsmat']
            mkdir(fpout)
            grid=db.grid.find_one({'_id':_id_run},{'ons':0,'offs':0})
            dff_r=dff_all[:,run['start_ses']:run['end_ses']]
            for id_neu,dff in enumerate(dff_r):
                _id_neu='%s%04d'%(_id_run, id_neu)
                plotTrialsMat(dff, grid, fpout, _id_neu, asp, pad)

def tm(dff,grid):
    trialsmat = np.empty([grid['row_end'], grid['col_end']], )
    for j, (on, off) in enumerate(zip(grid['ons_reor'], grid['offs_reor'])):
        trialsmat[j] = dff[on:off]
    return trialsmat

def plotTrialsMat(dff, grid, folder, _id_neu, cb_asp, cb_pad):
    trialsmat = tm(dff, grid)
    h = trialsmat.shape[0]
    w = trialsmat.shape[1]

    plt.clf()
    plt.figure(figsize=(w/16, h/8))

    plt.imshow(trialsmat, cmap=plt.get_cmap('turbo'))
    plt.colorbar(aspect=cb_asp,shrink=0.5,orientation='horizontal',pad=cb_pad)
    plt.margins(0)
    plt.gca().set_aspect(3)

    # gridlines
    for (pos,gridLabel) in zip(grid['col_poles'],grid['col_labels']):
        plt.axvline(pos-0.5, c='w', linestyle='--')
        plt.text(pos-10, 2, gridLabel, c='w')
    grid['row_poles'].append(grid['row_end'])
    for (pos, gridlabel) in zip(grid['row_poles'],grid['row_labels']):
        plt.axhline(pos-0.5, c='w', linestyle='--')
        plt.text(0, pos-0.5, str(gridlabel)+'°', c='w')

    plt.tight_layout()
    savePlot(plt.gca(), folder +'\\' + _id_neu + '.jpg')

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

def cd_visdriven_on_last_baseline(db):
    stats=['mean','std','rms','auc']
    projection = {'_id': 0, 'id_ses': 1}
    for stat in stats:
        projection[stat]=1

    for stat in stats:
        for mouse in db.mouse.find():
            print(mouse['name'])
            id_mouse=mouse['_id']
            n_sess = mouse['n_sess']
            try:
                ids_cd=loadmat(mouse['fp_crossday'])
                ids_cd = ids_cd.get(list(ids_cd.keys())[-1])
            except:
                continue

            id_ses=mouse['cond_poles'][0]-1
            _id_ses='%d%02d'%(id_mouse, id_ses)
            mouse['cond_poles'].insert(0, 0) # for plotting

            for i_run,id_run in enumerate(mouse['run_id_spont']):
                print('run%d' % id_run)
                ids_cd_select=[]
                vd_neus = [x['_id'] for x in list(db.neu_run2.find({'is_visdriven': True, 'id_mouse':id_mouse, 'id_ses': {'$lt':id_ses}, 'id_run': id_run + 1},
                                    {'_id': 1}))]
                #vd_neus = [x['id_neu'] for x in list(db.neu_run2.find({'is_visdriven': True, '_id_run': '%s%d' % (_id_ses, id_run + 1)},
                #                    {'_id': 0, 'id_neu': 1}))]
                stats = np.empty([n_sess, len(vd_neus)], dtype=float)
                stats[:] = np.nan
                for i,_id_neu in enumerate(vd_neus):
                    id_cd=db.neuron.find_one({'_id': '%s0%s' % (_id_neu[:3],_id_neu[4:])},{'_id':0,'id_cd':1})
                    if not 'id_cd' in id_cd.keys():
                        print('neu# %s not found in crossday' % _id_neu)
                        continue
                    if id_cd not in ids_cd_select:
                        id_cd=id_cd['id_cd']
                        ids_cd_select.append(id_cd)
                        sess=[]
                        stats_cd=[]
                        findquery = {'id_mouse': id_mouse, 'id_run': id_run, 'id_cd': id_cd, 'is_nonphys': {'$exists': 0}}
                        projection = {'_id': 0, stat: 1, 'id_ses': 1}
                        for neu in db.neu_run2.find(findquery,projection):
                            sess.append(neu['id_ses'])
                            stats_cd.append(neu[stat])
                        stats[sess,i]=stats_cd
                        if id_mouse == 0 and id_run==2:
                            stats[4,:]=np.nan

                mean=np.nanmean(stats*2,axis=1)
                err=sem(stats*2,axis=1,nan_policy='omit')
                print(np.count_nonzero(~np.isnan(stats[0])))

                plt.clf()
                fig = plt.figure(figsize=(20*(n_sess/16),15))
                l = plt.plot(mean,color='b',linewidth=5)[0]
                plt.fill_between(l.get_xdata(), mean - err, mean + err, color=l.get_color(),
                                 alpha=0.4, edgecolor='none')
                plt.plot(stats, alpha=0.5)
                plt.title('%s run %d crossday %s (subselect vis driven neurons on the last baseline day)' % (mouse['name'],id_run,stat))
                for x,label in zip(mouse['cond_poles'],mouse['conds']):
                    plt.axvline(x=x, label=label,color='r',linestyle='--')
                    plt.text(x-0.3,plt.gca().get_ylim()[1]*0.95,label,fontsize=40)
                ax = plt.gca()
                ax.set_xticks(list(range(n_sess)))
                ax.set_xticklabels(mouse['name_sess'], fontsize=15)
                plt.tight_layout()
                savePlot(fig, 'cd %s run %d mouse %d.jpg'%(stat,id_run,id_mouse))


def plot_trace_loop(folder,db,id_mouse):
    plt.figure(figsize=(15, 3.6))
    for mouse in db.mouse.find({'_id':id_mouse}):
        print(mouse['name'])
        for ses in db.session.find({'id_mouse':id_mouse,'id_ses':{'$gt':11}},{'fp_dff':1,'framerate':1,'n_neu':1}):
            _id_ses=ses['_id']
            print(_id_ses)
            dff_ses=loadmat(ses['fp_dff'])['dFF']
            for id_run in mouse['run_id_spont']:
                _id_run = '%s%01d' % (_id_ses, id_run)
                run = db.run.find_one({'_id': _id_run}, {'start_ses': 1, 'end_ses': 1})
                dff_run = dff_ses[:, run['start_ses']:run['end_ses']]
                for id_neu,dff in enumerate(dff_run):
                        _id = '%s%d%04d' % (_id_ses, id_run, id_neu)
                        filename = folder + '\\' + _id+'.png'
                        if not os.path.isfile(filename):
                            neu = db.neu_run2.find_one({'_id': _id},{'dff_lims_cd':1})
                            try:
                                plt.clf()
                                plt.ylim(neu['dff_lims_cd'][0], neu['dff_lims_cd'][1])
                                plt.plot(dff, linewidth=0.9)
                                plt.margins(0.01)
                                plt.savefig(filename, bbox_inches='tight')
                            except:
                                print(_id)



def stimBlocks(ax):
    ymin,ymax=ax.get_ylim()
    for event in stimEvents:
        ax.plot([event.on,event.off], [ymax,ymax], c=event.stim.color, marker='|',linewidth=20)
        ax.vlines([event.on,event.off], ymin, ymax,colors=['k','k'],linestyles='dashed',alpha=0.3)
    return ax



























