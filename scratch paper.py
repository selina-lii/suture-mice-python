

def cd_select(thresholds,stats,):
    thresholds=[0.5,0.6,0.7,0.8,0.9]
    stats=['mean','std','rms']
    fieldquery = {'_id': 0, 'id_ses': 1}
    for stat in stats:
        fieldquery[stat]=1


def plot_crossday():

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
                    for neu in db.neu_run_spont.find(findquery, fieldquery):
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
