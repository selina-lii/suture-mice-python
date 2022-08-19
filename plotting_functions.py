import os.path
import matplotlib.pyplot as plt
from math import inf
from LyxTools import *
import numpy as np
from scipy.stats import sem


def average_response_spontaneous_wrapper(db):
    avgact_vd_neu_in_spont = np.zeros([2, 7, 30])
    stats=mouse_loop(db,average_response_spontaneous,)
    for mouse in db.mouse.find():
        n = mouse['n_sess']
        m_n = mouse['name']
        names = mouse['name_sess']
        for i, id_run in enumerate(mouse['run_id_stim']):
            tt = "average response of driven neurons in spontaneous run - %s run%d" % (m_n, id_run)
            fn = 'avgact_%s_run%d.jpg' % (m_n, id_run)
            d = avgact_vd_neu_in_spont[i, mouse['_id'], :n]
            plot_mouse_stat(d, n, tt, names, fn)

def average_response_spontaneous(db,ses):

    for run in db.run.find():
        means=get_neu_runs(db,run,proj('mean'))
        return np.mean(means)


    for mouse in db.mouse.find():
        id_mouse = mouse['_id']
        for j in range(mouse['n_sess']):
            ses = db.session.find_one({'_id': '%d%02d' % (mouse['_id'], j)}, {'run_id_stim': 1, 'fp_dff': 1})
            neus = list(db.neuron.find({'id_ses': ses['_id']},
                                       {'_id': 0, 'is_visdriven_1': 1, 'is_visdriven_2': 1, 'is_visdriven_3': 1}))
            dff = loadmat(ses['fp_dff'])['dFF']
            for i, id_run in enumerate(ses['run_id_stim']):
                run = db.run.find_one({'_id': '%s%d' % (ses['_id'], id_run)}, {'start_ses': 1, 'end_ses': 1})
                ids = []
                for k, neu in enumerate(neus):
                    if neu['is_visdriven_%d' % id_run] is True:
                        ids.append(k)
                if len(ids) != 0:
                    dff_run = dff[ids, run['start_ses']:run['end_ses']]
                    avgact_vd_neu_in_spont[i, id_mouse, j] = np.mean(dff_run)

    return avgact_vd_neu_in_spont

def plot_trialsmat_loop(db, config, id_mouse_sel=None):
    if id_mouse_sel is not None:
        query = {'id_mouse': id_mouse_sel}
    else:
        query = {}
    for ses in db.session.find(query, {'fp_dff': 1, 'run_id_stim': 1, 'n_neu': 1, 'id_mouse': 1}):
        id_mouse = ses['id_mouse']
        asp = config.plt_tm_cb_asp[id_mouse]
        pad = config.plt_tm_cb_pad[id_mouse]
        dff_all = loadmat(ses['fp_dff'])['dFF']
        for id_run in ses['run_id_stim']:
            _id_run = ses['_id'] + str(id_run)
            run = db.run.find_one({'_id': _id_run}, {'start_ses': 1, 'end_ses': 1, 'fp_trialsmat': 1})
            fpout = run['fp_trialsmat']
            mkdir(fpout)
            grid = db.grid.find_one({'_id': _id_run}, {'ons': 0, 'offs': 0})
            dff_r = dff_all[:, run['start_ses']:run['end_ses']]
            for id_neu, dff in enumerate(dff_r):
                _id_neu = '%s%04d' % (_id_run, id_neu)
                plotTrialsMat(dff, grid, fpout, _id_neu, asp, pad)

def tm(dff, grid):
    trialsmat = np.empty([grid['row_end'], grid['col_end']], )
    for j, (on, off) in enumerate(zip(grid['ons_reor'], grid['offs_reor'])):
        trialsmat[j] = dff[on:off]
    return trialsmat


def plotTrialsMat(dff, grid, folder, _id_neu, cb_asp, cb_pad):
    trialsmat = tm(dff, grid)
    h = trialsmat.shape[0]
    w = trialsmat.shape[1]

    plt.clf()
    plt.figure(figsize=(w / 16, h / 8))

    plt.imshow(trialsmat, cmap=plt.get_cmap('turbo'))
    plt.colorbar(aspect=cb_asp, shrink=0.5, orientation='horizontal', pad=cb_pad)
    plt.margins(0)
    plt.gca().set_aspect(3)

    # gridlines
    for (pos, gridLabel) in zip(grid['col_poles'], grid['col_labels']):
        plt.axvline(pos - 0.5, c='w', linestyle='--')
        plt.text(pos - 10, 2, gridLabel, c='w')
    grid['row_poles'].append(grid['row_end'])
    for (pos, gridlabel) in zip(grid['row_poles'], grid['row_labels']):
        plt.axhline(pos - 0.5, c='w', linestyle='--')
        plt.text(0, pos - 0.5, str(gridlabel) + '°', c='w')

    plt.tight_layout()
    savePlot(plt.gca(), folder + '\\' + _id_neu + '.jpg')


def histograms_of_stim(fp=None):
    if fp is None:
        fp = 'C:\\Users\\selinali\lab\\sut\\2022-7-22-plots\\stim_stats'
    stats = ['js_max']
    for run in db.run.find({'type': 'stim'}, {'name_mouse': 1, 'cond_code': 1, 'hour_code': 1}):
        _id_run = run['_id']
        ll = list(db.neu_run.find({'_id_run': _id_run}, {'js_max': 1}))
        for stat in stats:
            mkdir(fp + '\\' + stat)
            ll_s = [x[stat] for x in ll]
            mean = np.mean(ll_s)
            std = np.std(ll_s)

            plt.clf()
            plt.hist(ll_s, bins=15)
            plt.axvline(mean + std, color='r')
            plt.title('%s %s%s histogram of %s for stim run%s' % (
            run['name_mouse'], run['cond_code'], run['hour_code'], stat, _id_run[-1]))
            plt.xlabel('mean=%.6f, std=%.6f' % (mean, std))
            fn = '%s_%s_%s_%s_%s%s' % (
            stat, _id_run[-1], run['name_mouse'], _id_run[1:3], run['cond_code'], run['hour_code'])
            plt.tight_layout()
            savePlot(plt.gca(), fp + '\\' + stat + '\\' + fn + '.jpg')


types = ['b', 'l', 'u', 'u', 'u', 'u']
thres_u = [0.6, 10, inf, 2, 200, 3]
thres_l = [0, 0, -2, -inf, -inf, -inf]
stats = ['mean', 'max', 'min', 'std', 'kurt', 'rms']


def histograms_of_spont(fp=None):
    if fp is None:
        fp = 'C:\\Users\\selinali\lab\\sut\\2022-7-22-plots\\spont_stats'
    for run in db.run.find({'type': 'spont'}, {'name_mouse': 1, 'cond_code': 1, 'hour_code': 1}):
        _id_run = run['_id']
        ll = list(db.neu_run_spont.find({'_id_run': _id_run}))
        for stat, out_u, out_l in zip(stats, thres_u, thres_l):
            mkdir(fp + '\\' + stat)
            ll_s = [x[stat] for x in ll]
            ll_s_phys = []
            outs = []
            cnt = 0
            for x in ll_s:
                if x < out_l or x > out_u:
                    cnt = cnt + 1
                    x = int(x)
                    if x == 0 or x == -1 or x == 1:
                        pass
                    else:
                        outs.append(int(x))
                else:
                    ll_s_phys.append(x)

            mean = np.mean(ll_s_phys)
            std = np.std(ll_s_phys)
            median = np.median(ll_s_phys)

            plt.clf()
            plt.hist(ll_s_phys, bins=20)
            plt.axvline(mean + std, color='r')
            plt.axvline(mean, color='y')
            plt.axvline(mean - std, color='r')
            plt.axvline(median, color='k')
            plt.text(median, 0, 'median', color='k', rotation=-90)
            plt.text(mean, 0, 'mean', color='y', rotation=-90)
            plt.text(median + std, 0, '+1std', color='r', rotation=-90)
            plt.text(median - std, 0, '-1std', color='r', rotation=-90)
            tex = 'mean=%.6f\nmedian=%.6f:\n+1std=%.6f\n-1std=%.6f\n' % (mean, median, mean + std, mean - std)
            plt.text(plt.gca().get_xlim()[1] * 0.6, plt.gca().get_ylim()[1] * 0.6, tex)
            plt.title('%s %s%s histogram of %s for spontaneous run%s' % (
            run['name_mouse'], run['cond_code'], run['hour_code'], stat, _id_run[-1]))
            plt.xlabel('#outliers%d:%s' % (cnt, outs))
            plt.tight_layout()
            fn = '%s_%s_%s_%s_%s%s' % (
            stat, _id_run[-1], run['name_mouse'], _id_run[1:3], run['cond_code'], run['hour_code'])
            savePlot(plt.gca(), fp + '\\' + stat + '\\' + fn + '.jpg')


def cd_std_auc_wrapper(db):
    mouse_loop(db, cd_std_auc, cd=True, stat='auc')

def cd_std_auc(mouse):
    for id_run in range(mouse.n_runs):
        stats_mean = np.empty([len(stats), mouse.n_sess], dtype=float)
        stats_err = np.empty([len(stats), mouse.n_sess], dtype=float)
        for ses in db.session.find({'id_mouse': mouse._id}, {'id_ses': 1}):
            id_ses = ses['id_ses']

            if mouse._id == 0 and id_run > 1 and id_ses == 4:
                stats_mean[:, id_ses] = np.nan
                stats_err[:, id_ses] = np.nan
                continue

            neus = list(
                db.neu_run2.find('_id_run': '%d%02d%d' % (mouse._id, id_ses, id_run), 'is_nonphys': {'$exists': 0}},
                                 projection))
            for s, stat in enumerate(stats):
                stat_vals = [x[stat] for x in neus]
                stats_mean[s, id_ses] = np.mean(stat_vals)
                stats_err[s, id_ses] = sem(stat_vals)

        for s, stat in enumerate(stats):
            tt = "average %s for %s run %d" % (stat, mouse.name, id_run)
            fn = '%s %s run%d.jpg' % (stat, mouse.name, id_run)
            plot_mouse_stat(stats_mean[s], tt, mouse.name_sess, fn, err=stats_err[s])



def cd_compare_conditions_mean(db):


def crossday_mean_acitivty_across_thresholds_wrapper(db):
    stats = ['mean', 'std', 'rms']
    for stat in stats:
        mouse_loop(db, crossday_mean_acitivty_across_thresholds, cd=True, stat=stat,outdir=db.config.find_one()['_testpath'])

def crossday_mean_acitivty_across_thresholds():
    thresholds = [0.5, 0.6, 0.7, 0.8, 0.9]

    id_runs=mouse.run_id_spont
    for i, id_run in enumerate(id_runs):
        print('run%d' % id_run)
        ids_cd_selected=[x['_id'] for x in list(db.neu_cd.find(dict(is_on_all_conds=True, id_mouse=id_mouse), {'_id': 1}))]
        cd_stats=get_cd_stats(db, mouse, id_run, stat, ids_cd_selected)
        mean, err, _ = meanerr(cd_stats)
        title = ' crossday %s for %s spontaneous run %d' % (mouse['name'], stat, id_run)
        filename = '%s\\%d%d_cd_meanact_all_conds.jpg' % (outdir,id_run,mouse._id)
        legend=['>%d%% (n=%s)' % (t * 100, n) for t, n in zip(thresholds, n_valid)]
        cd_plot_line(mouse, title, filename, mean=mean, err=err,legend=legend)


        for i, id_run in enumerate(mouse['run_id_spont']):
            print('run%d' % id_run)
            stats_mean = np.empty([len(stats), len(thresholds), n_sess])
            stats_err = np.empty([len(stats), len(thresholds), n_sess])
            n_valid = []
            for k, thres in enumerate(thresholds):
                valid_ids = np.where(np.sum(ids_cd > 0, axis=1) > thres * n_sess)[0].tolist()
                n_valid.append(len(valid_ids))
                stat_result = np.empty([len(stats), n_sess, len(valid_ids)], dtype=float)
                stat_result[:] = np.nan
                for i, id_cd in enumerate(valid_ids):
                    findquery = {'id_mouse': id_mouse, 'id_run': id_run, 'id_cd': id_cd, 'is_nonphys': {'$exists': 0}}
                    for neu in db.neu_run_spont.find(findquery, projection):
                        for j, stat in enumerate(stats):
                            stat_result[j, neu['id_ses'], i] = neu[stat]
                stats_mean[:, k, :] = np.nanmean(stat_result, axis=2)
                stats_err[:, k, :] = sem(stat_result, axis=2, nan_policy='omit')
                print(np.count_nonzero(~np.isnan(stat_result[0]), axis=1))


def crossday_mean_acitivty_subselect_all_conds_wrapper(db):
    stats = ['mean', 'std', 'rms', 'auc']
    for stat in stats:
        mouse_loop(db, crossday_mean_acitivty_subselect_all_conds, cd=True, stat=stat,outdir=db.config.find_one()['_testpath'])

def crossday_mean_acitivty_subselect_all_conds():
    id_runs=mouse.run_id_spont
    for i_run, id_run in enumerate(id_runs):
        print('run%d' % id_run)
        ids_cd_selected=[x['_id'] for x in list(db.neu_cd.find(dict(is_on_all_conds=True, id_mouse=id_mouse), {'_id': 1}))]
        cd_stats=get_cd_stats(db, mouse, id_run, stat, ids_cd_selected)
        mean, err, _ = meanerr(cd_stats)
        title = '%s run %d crossday %s (subselect neurons that are present on at least one day in all conditions)' % (
                mouse.name, id_run, stat)
        filename = '%s\\%d%d_cd_meanact_all_conds.jpg' % (outdir,id_run,mouse._id)
        cd_plot_line(mouse, title, filename, mean=mean, err=err)


def crossday_meanact_wrapper(db):
    mouse_loop(db, crossday_meanact, cd=True, outdir=db.config.find_one()['_testpath'])

def crossday_meanact(outdir):
    id_runs = mouse.run_id_spont
    for i_run, id_run in enumerate(id_runs):
        print('run%d' % id_run)
        cd_stats=get_cd_stats(db,mouse,id_run,stat)
        title = '%s run %d crossday mean activity' % (
                mouse.name, id_run)
        filename = '%s\\%d%d_cd_meanact.jpg' % (outdir,id_run,mouse._id)
        cd_plot_line(mouse, title, filename, data=cd_stats,alpha_data=0.4)


def cd_visdriven_on_last_baseline_wrapper(db):
    stats = ['mean', 'std', 'rms', 'auc']
    for stat in stats:
        mouse_loop(db, cd_visdriven_on_last_baseline, cd=True, stat=stat,outdir=db.config.find_one()['_testpath'])

def cd_visdriven_on_last_baseline(db, mouse, stat, outdir):
    id_ses=mouse.cond_poles[0]-1
    id_runs=mouse.run_id_spont

    for i_run, id_run in enumerate(id_runs):
        print('run%d' % id_run)
        vd_neus = get_neu_runs(db, mouse._id, id_ses, id_run + 1,also=dict(is_visdriven=True),fields='id_cd') #TODO +1: oops this is hard coding..
        cd_stats=get_cd_stats(db, mouse, id_run, stat, neu_to_cd(vd_neus))
        mean, err, _ = meanerr(cd_stats)
        title = '%s run %d crossday %s (subselect vis driven neurons on the last baseline day)' % (
                mouse.name, id_run, stat)
        filename = '%s\\%d%d_cd_%s.jpg' % (outdir,id_run, mouse._id,stat)
        cd_plot_line(mouse, title, filename, data=cd_stats, mean=mean, err=err)

def converter(**kwargs):
    l = list(kwargs.values())
    if len(kwargs)==1:
        l=list(kwargs.values())[0]
        if len(l)==3:

    else:
        print(1)


def cd_mean_between_conds(db,mouse,stat,outdir):
    ids_cd = load_ids_cd(mouse)
    if ids_cd is None:
        return

    ses_queries = [dict(cond_code='B', hour={'$gt': 24}),
                   dict(cond_code='S', hour={'$gt': 5, '$lt': 24}),
                   dict(cond_code='S', hour={'$gt': 24})]
    find_query=dict(id_mouse=mouse._id)
    labels=['late baseline','early suture','late suture']
    ids_cd_selected=[]
    for ses_q in ses_queries:
        find_query.update(ses_q)
        id_sess=[x['id_ses'] for x in db.session.find(find_query,{'id_ses':1})]
        ids_cd_selected.append([i if sum(x[id_sess])!=0 for x in ids_cd])


    id_runs=mouse.run_id_spont
    for i_run, id_run in enumerate(id_runs):
        print('run%d' % id_run)
        for cd_neu in mouse.neu_cd.find(dict(id_mouse=mouse._id,id_ses={'$in':id_sess})):
            vd_neus = get_neu_runs(db, mouse._id, id_sess, id_run + 1,proj('id_cd')) #TODO +1: oops this is hard coding..
            cd_stats=get_cd_stats(db, mouse, id_run, stat, neu_to_cd(vd_neus))
            mean, err, _ = meanerr(cd_stats)
            title = '%s run %d crossday %s (subselect vis driven neurons on the last baseline day)' % (
                    mouse.name, id_run, stat)
            filename = '%s\\%d%d_cd_%s.jpg' % (outdir,id_run, mouse._id,stat)
        cd_plot_bar(title, filename, data=cd_stats, mean=mean, err=err)


#pick neu-runs
def get_neu_runs(db, id_mouse=None, id_ses=None, id_run=None, also=None,fields=None):
    
    if type(id_ses) is list:
        id_ses = {'$in': id_ses}
    find_query=dict(is_nonphys={'$exists':False}, id_mouse=id_mouse, id_ses=id_ses, id_run=id_run)
    find_query.update(also)
    
    if isinstance(fields,string):
        return [x[field] for x in list(db.neu_run2.find(dict(_id_run=run._id), proj(fields,_id=False)))]
    elif isinstance(fields,dict):
        return list(db.neu_run2.find(find_query, fields))

# bridges within database
'''neu->neu_cd'''
def neu_to_cd(neus_selected):
    ids_cd_selected = []
    for i, neu in enumerate(neus_selected):
        try:
            ids_cd_selected.append(neu['id_cd'])
        except:
            print('neu# %s not found in crossday' % neu['_id'])
    ids_cd_selected=np.unique(ids_cd_selected).tolist()
    return ids_cd_selected
'''neu_cd->stat'''
def get_cd_stats(db, mouse, id_run, stat, ids_selected=None,sess_selected=None):
    n_neu_cd=mouse.nneu_cd if ids_selected is None else len(ids_selected)
    stats = np.empty([mouse.n_sess, n_neu_cd], dtype=float)
    stats[:]=np.nan
    findquery = dict(id_mouse=mouse._id, id_run=id_run, is_nonphys={'$exists': 0})
    if sess_selected is not None:
        findquery['id_ses']={'$in':sess_selected}
    projection = proj(['id_ses', stat], inc__id=False)

    if ids_selected is None:
        ids_selected=list(range(n_neu_cd))

    for i,id_cd in enumerate(ids_selected):
        findquery['id_cd'] = id_cd
        neus = list(db.neu_run2.find(findquery, projection))
        stats[[x['id_ses'] for x in neus], i] = [x[stat] for x in neus]
    return stats

def get_neu_run_stats(db,run,stat,ids_selected=None):
    if ids_selected is not None:
        stats=[]
        for id_neu in ids_selected:
            stats.append(db.neu_run2.find(dict(_id='%s%04d'%(run._id,id_neu)),dict(stat=stat))[stat])
        return stats
    return [x[stat] for x in list(db.neu_run2.find(dict(_id_run=run._id),dict(stat=stat)))]

# loops
def mouse_loop(db, func, cd=False,**kwargs):
    rr=[]
    for mouse in db.mouse.find():
        mouse = DBinterface(mouse)
        print(mouse.name)
        if cd:
            r = func(db,mouse,**kwargs)
        else:
            r = run_loop(db, func, mouse, **kwargs)
        rr.append(r)
    return rr

def run_loop(db, func, mouse, **kwargs):
    rr=[]
    for ses in db.session.find(dict(id_mouse=mouse._id)):
        ses = DBinterface(ses)
        print(ses._id)
        r=func(db=db, ses=ses, **kwargs)
        rr.append()
    return rr

# plottings
# TODO the name 'show_all' is really bad
def cd_plot_line(mouse, title, filename, data=None, mean=None, err=None, legend=None, show_cond_lines=True, alpha_data=0.5):
    plt.figure(figsize=(1.25 * mouse.n_sess, 15))
    ax=plt.gca()

    if mean is not None:
        plot_meanerr(ax,mean,err)
    if legend is not None:
        plt.legend(legend, fontsize=8)
    if data is not None:
        plt.plot(data, alpha=alpha_data)
    if show_cond_lines:
        plot_cond_poles(ax,mouse,fontsize=20)

    set_xticks_sess_names(ax,mouse)
    ax.tight_layout()
    ax.title(title,fontsize=40)
    savePlot(ax, filename)

def cd_plot_bar( title, filename, mean, bar_labels,data=None, alpha_data=0.5):
    plt.figure(figsize=(len(mean), 12))
    ax = plt.gca()

    ax.bar(bar_labels, mean)
    if data is not None:
        plt.plot(data, alpha=alpha_data,color='k',marker='o')

    ax.tight_layout()
    ax.title(title,fontsize=40)
    savePlot(ax, filename)



def plot_mouse_stat(mouse,data, title, filename, err=None):
    plt.clf()
    plt.figure(figsize=(len(data), 7.5))
    l = plt.plot(data)[0]
    if err is not None:
        plt.fill_between(l.get_xdata(), data + err, data - err, color=l.get_color(), alpha=0.2)
    plt.title(title)
    set_xticks_sess_names(plt.gca(),mouse,fontsize=8)
    plt.savefig(filename, bbox_inches='tight', pad_inches=0.2)
    plt.close()

def set_xticks_sess_names(ax,mouse,fontsize):
    ax.set_xticks(list(range(mouse.n_sess)))
    ax.set_xticklabels(mouse.name_sess, fontsize=fontsize)
    return ax

def plot_cond_poles(ax,mouse):
    cond_poles = [x for x in mouse.cond_poles]
    cond_poles.insert(0, 0)
    for x, label in zip(cond_poles, mouse.conds):
        ax.axvline(x=x, label=label, color='r', linestyle='--')
        ax.text(x - 0.3, ax.get_ylim()[1] * 0.95, get_cond_name(label), fontsize=40)
    return ax

def plot_meanerr(ax,mean,err=None):
    if len(mean.shape) == 1:
        l = ax.plot(mean, color='b', linewidth=5)[0]
        if err is not None:
            ax.fill_between(l.get_xdata(), mean - err, mean + err, color=l.get_color(),
                             alpha=0.2, edgecolor='none')
    else:
        if err is not None:
            for m, e in zip(mean, err):
                l = ax.plot(mean, color='b', linewidth=2)[0]
                ax.fill_between(l.get_xdata(), m - e, m + e, color=l.get_color(),
                                 alpha=0.4, edgecolor='none')
        else:
            ax.plot(mean, color='b', linewidth=2)
    return ax

def plot_traces(db, ses, id_runs, folder):
    plt.figure(figsize=(15, 3.6))
    dff_ses = loadmat(ses.fp_dff)['dFF']
    for id_run in id_runs:
        run = db.run.find_one(dict(_id='%s%01d' % (ses._id, id_run)),
                              dict(start_ses=1, end_ses=1))
        dff_run = dff_ses[:, run['start_ses']:run['end_ses']]
        for id_neu, dff in enumerate(dff_run):
            _id = '%s%d%04d' % (ses._id, id_run, id_neu)
            filename = folder + '\\' + _id + '.png'
            if not os.path.isfile(filename):
                neu = db.neu_run2.find_one(id(_id), {'dff_lims_cd': 1})
                try:
                    plt.clf()
                    plt.ylim(neu['dff_lims_cd'][0], neu['dff_lims_cd'][1])
                    plt.plot(dff, linewidth=0.9)
                    plt.margins(0.01)
                    plt.savefig(filename, bbox_inches='tight')
                except:
                    print(_id)

def stimBlocks(ax,stim_events,color_scheme):
    ymin, ymax = ax.get_ylim()
    for event in stim_events:
        ax.plot([event.on, event.off], [ymax, ymax], c=color_scheme[event.label], marker='|', linewidth=20)
        ax.vlines([event.on, event.off], ymin, ymax, colors=['k', 'k'], linestyles='dashed', alpha=0.3)
    return ax

# level 1-basic stats
def meanerr(data): # scaling if you want to inflate the trend so you get a better look..?
    mean = np.nanmean(data, axis=1)
    err = sem(data, axis=1, nan_policy='omit')
    ndata = len(data[0])
    return mean, err, ndata

#helpers
def get_cond_name(cond_code):
    if cond_code == 'B':
        return 'baseline'
    elif cond_code == 'S':
        return 'suture'
    elif cond_code == 'U':
        return 'unsuture'
    elif cond_code == 'R':
        return 'resuture'
    else:
        NameError('invalid condition code')

    #projections
    
def proj(fields, _id=True, reverse=False):
    p = dict()
    for fn in fields:
        p[fn] = 1
    if not _id:
        p['_id'] = 0
    if reverse:
        for key, val in p.items():
            p[key] = 1 - val
    return p


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


def crossday_imshow():
    for mouse in db.mouse.find():
        id_mouse=mouse['_id']
        fp_crossday=config.workdir+config.crossdaydir+'\\'+mouse['name']+'.mat'
        ids_cd=load_ids_cd(mouse)
        if ids_cd is None:
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
