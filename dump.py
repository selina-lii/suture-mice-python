'''
def insert_dff_and_spks(db):
    for i, mouse in enumerate(mice):
        if True:#mice is 'Sut1':
            allfiles = glob.glob(mouseDir + '\\' + mouse + '\\*_Suite2p_*.mat')
            print(mouse)
            if len(allfiles) > 0:
                assert mice_n_files[i] == len(allfiles)
                names = run_names(mice_days[i])
                for j, file in enumerate(allfiles):
                    if True:#j==12:
                        name=names[j]
                        print(name + ' ' + str(j))
                        data = loadmat(file)['suite2pData']
                        sesPoles=intarr(data.startIdx[1:])

                        mouse_id = i
                        run_id = j

                        #timeseries
                        dffss = np.hsplit(data.dFF,sesPoles)
                        spksss =np.hsplit(data.spks,sesPoles)

                        run = db.run.find_one({"$and": [{"mouse_id": mouse_id}, {"run_id": run_id}]})
                        _run_id=run['_id']

                        for ses_id,(dffs,spkss) in enumerate(zip(dffss,spksss)):
                            ses = db.session.find_one({"$and": [{"mouse_id": mouse_id}, {"run_id": run_id}, {"ses_id": ses_id}]})
                            _ses_id = ses['_id']
                            for neu_id,(dff,spks) in enumerate(zip(dffs,spkss)):
                                query={"$and": [{"mouse_id": mouse_id}, {"run_id": run_id}, {"neu_id": neu_id}]}
                                neu = db.neuron.find_one(query)
                                _neu_id = neu['_id']
                                dff=[x.item() for x in dff]
                                spks=[x.item() for x in spks]
                                #TODO: js_dists and vis_driven
                                dff_doc = DBinterface(mouse_id=mouse_id, run_id=run_id,ses_id=ses_id,neu_id=neu_id,
                                                      dff=dff, len=len(dff),
                                                      _neu=_neu_id)
                                _dff_id = dff_doc.insert(db.dff)
                                spks_doc = DBinterface(mouse_id=mouse_id, run_id=run_id, ses_id=ses_id, neu_id=neu_id,
                                                       spks=spks, len=len(spks),
                                                       _neu=_neu_id)
                                _spks_id = spks_doc.insert(db.spks)
                                db.neuron.update_one(query,{"$set": {"_dff":_dff_id.inserted_id, "_spks":_spks_id.inserted_id}})
'''
import os.path
import shutil

'''def set_parent(db):
    # mouse<->run
    for mouse_name in mice:
        if True:#mouse_name == 'Sut4':
            mouse=db.mouse.find_one({"mouse_name":mouse_name})
            if mouse is not None:
                _id_mouse=mouse['_id']
                print(_id_mouse)
                query={"mouse_name": mouse_name}
                db.run.update_many(query, {"$set": {"_parent": _id_mouse}})
                runs=list(db.run.find(query).sort('run_id'))
                _id_runs=[]
                for run in runs:
                    _id_runs.append(run['_id'])
                db.mouse.update_one({"mouse_name":mouse_name},{"$set":{"_children": _id_runs}})
    # run<->ses
    for run in db.run.find({}):
        _id_run=run['_id']
        print(_id_run)
        query={"$and":[{"mouse_name":run['mouse_name']},{"run_id":run['run_id']}]}
        db.session.update_many(query, {"$set": {"_parent": _id_run}})
        sessions=list(db.session.find(query).sort('ses_id'))
        _id_sess = []
        for ses in sessions:
            _id_sess.append(ses['_id'])
        db.run.update_one(query, {"$set": {"_children": _id_sess}})
    # ses<->stimEvent
    for ses in db.session.find({}):
        _id_ses=ses['_id']
        query={"$and":[{"mouse_name":ses['mouse_name']},{"run_id":ses['run_id']},{"ses_id":ses['ses_id']}]}
        stimEvents=list(db.stim_event.find(query).sort('stim_id')) #why do you have to do list here but not with the other one
        if len(stimEvents) is not 0:
            db.stim_event.update_many(query, {"$set": {"_parent": _id_ses}})
            _id_stim_events = []
            for event in stimEvents:
                _id_stim_events.append(event['_id'])
            db.session.update_one(query, {"$set": {"_stim_events": _id_stim_events}})
    print('ses<->stimEvent:done')
    # stimEvent<->stimObject ???
    # run<->neuron
'''
'''def custom__id():
    for neu in db.neuron.find():
        neu['_id']="%01d%02d%01d%04d"%(neu['mouse_id'],neu['run_id'],0,neu['neu_id'])
    for stim in db.stim_event.find():
        stim['_id']="%01d%02d%01d%03d"%(stim['mouse_id'],stim['run_id'],stim['ses_id'],stim['stim_id'])
        ses['min_trial_len']=
        ses['oris']=
        ses['n_oris']=
        ses['n_oris_each']='''
'''def messed_up(): #damn... this is stupid i forgot to zero stim id
    db.stim_event.delete_many({})
    for each in db.stim_event_backup.find():
        db.stim_event.insert_one(each)
    for run in db.run.find():
        _id_run=run['_id']
        n_stim=0
        for j in range(run['n_ses']-1):
            _id_ses=run['_id']+str(j+1)
            ses=db.session.find_one({'_id':_id_ses})
            if ses['is_stim'] and (db.stim_event.find_one({'_id':_id_ses+'000'}) is None):
                for stim in db.stim_event.find({"$and":[{"mouse_id":int(_id_ses[0])},{"run_id":int(_id_ses[1:3])},
                                                        {"ses_id":int(_id_ses[3])}]}):
                    query= {'_id': stim['_id']}
                    stim_id_ses=int(stim['_id'][-3:])-n_stim #zero
                    stim['_id']=_id_ses+'%03d'%stim_id_ses
                    stim['stim_id_glob'] = stim['stim_id']
                    stim['stim_id_ses'] = stim_id_ses
                    db.stim_event.delete_one(query)
                    db.stim_event.insert_one(stim)
            elif ses['is_stim']: #preserve order
                n_stim=n_stim+ses['n_stim']
                for stim in db.stim_event.find({"ses_id": _id_ses[3]}):
                    query = {'_id': stim['_id']}
                    stim['stim_id_ses']=stim['stim_id_glob']=stim['stim_id']
                    stim.pop('stim_id')
                    db.stim_event.delete_one(query)
                    db.stim_event.insert_one(stim)

                    #stim['stim_id_ori_ses']=
                    #stim['stim_id_ori_glob']= #set in the following few functions
'''
'''def change_ses_and_run_id():
    for stim in db.stim_event.find({}):
        db.stim_event.update_one({'_id': stim['_id']},{"$set": {'ses_id': stim['_id'][:4]}})
    for stim in db.stim_event.find({}):
        db.stim_event.update_one({'_id': stim['_id']}, {"$set": {'run_id': stim['_id'][:3]}})
    for neu in db.neuron.find({}):
        db.neuron.update_one({'_id': neu['_id']},{"$set": {'run_id': neu['_id'][:3]}})
'''
'''
def set_list_and_n_of_oris_in_session():
    for ses in db.session.find({"is_stim":True}):
        _id_ses=ses['_id']
        stims=list(db.stim_event.find({"ses_id":_id_ses}).sort("label_id"))
        labels=[]
        for stim in stims:
            labels.append(stim['label'])
        unq_labels,counts=np.unique(labels,return_counts=True)
        unq_labels=unq_labels.tolist()
        counts=counts.tolist()
        db.session.update_one({'_id':_id_ses},{"$set":{'stim_labels':unq_labels}})
        db.session.update_one({'_id': _id_ses}, {"$set": {'n_stim_ori':counts}})
        for label in unq_labels:
            stims = list(db.stim_event.find({"$and":[{"ses_id": _id_ses}, {"label": label}]}).sort("start_ses"))
            for id,stim in enumerate(stims):
                db.stim_event.update_one({'_id': stim['_id']}, {"$set": {'stim_id_ori_ses': id}})

    for run in db.run.find({}):
        _id_run = run['_id']
        stims = list(db.stim_event.find({"run_id": _id_run}).sort("label_id"))
        labels = []
        for stim in stims:
            labels.append(stim['label'])
        unq_labels, counts = np.unique(labels, return_counts=True)
        unq_labels = unq_labels.tolist()
        counts = counts.tolist()
        db.run.update_one({'_id': _id_run}, {"$set": {'stim_labels': unq_labels}})
        db.run.update_one({'_id': _id_run}, {"$set": {'n_stim_ori': counts}})
        for label in unq_labels:
            stims = list(db.stim_event.find({"$and": [{"run_id": _id_run}, {"label": label}]}).sort("start_glob"))
            for id, stim in enumerate(stims):
                db.stim_event.update_one({'_id': stim['_id']}, {"$set": {'stim_id_ori_run': id}})
'''
'''
def set_stim_gaps({"is_stim":True}):
    for ses in db.session.find():
        n_stim = ses['n_stim']
        for i in range(n_stim - 1):
            _id_p = "%s%03d" % (ses['_id'], i)
            prev = db.stim_event.find_one({"_id": _id_p})
            _id_n = "%s%03d" % (ses['_id'], i + 1)
            next = db.stim_event.find_one({"_id": _id_n})
            gap = next['start_ses'] - prev['end_ses']
            db.stim_event.update_one({'_id': _id_n}, {"$set": {"gap_to_prev": gap}})
            db.stim_event.update_one({'_id': _id_p}, {"$set": {"gap_to_next": gap}})
            
def set_session_min_trial_len():
    for ses in db.session.find():
        _id_ses = ses['_id']
        sortres = list(db.stim_event.find({"$and": [{"mouse_id": int(_id_ses[0])}, {"run_id": int(_id_ses[1:3])},
                                                    {"ses_id": int(_id_ses[3])}]}).sort("gap_to_next"))
        db.session.update_one({'_id': _id_ses}, {"$set": {'iti_min': sortres[1]["gap_to_next"]}})
        db.session.update_one({'_id': _id_ses}, {"$set": {'iti_max': sortres[-1]["gap_to_next"]}})
        
def max_stim_len_for_session():
    for ses in db.session.find({"is_stim": True}):
        _id_ses = ses['_id']
        stims = db.stim_event.find({"ses_id": _id_ses}, {'len': 1})
        maxlen = 0
        for stim in stims:
            maxlen = max(maxlen, stim['len'])
        db.session.update_one({'_id': _id_ses}, {"$set": {'max_stim_len': maxlen}})

    for run in db.run.find({}):
        _id_run = run['_id']
        stims = db.stim_event.find({"run_id": _id_run}, {'len': 1})
        maxlen = 0
        for stim in stims:
            maxlen = max(maxlen, stim['len'])
        db.run.update_one({'_id': _id_run}, {"$set": {'max_stim_len': maxlen}})

def update_fp_dff_run():
    for run in db.run.find():
        fp_dff = 'C:\\Users\\selinali\\lab\\RDE20141\\2022-7-19-dff\\' + run['fp_dff'].split('\\')[-1]
        db.run.update_one({'_id': run['_id']}, {'$set': {'fp_dff': fp_dff}})

def mark_n_stim_below_thres():
    for ses in db.session.find({'n_stim_ori': {'$exists': 1}}):
        labels = [i for i in ses['n_stim_ori'] if i <= 6]
        if len(labels) != 0:
            print(ses['_id'])
            db.session.update_one({'_id': ses['_id']}, {'$set': {'stim_labels_low_n': labels}})

def insert_neurons(db):  # ayyyyyy the grandest part of the project
    for i, mouse in enumerate(mice):
        allfiles = glob.glob(mouseDir + '\\' + mouse + '\\*_Suite2p_*.mat')
        print(mouse)
        if mouse == 'Sut1' or mouse == 'Sut3' or mouse == 'Sut4':
            if len(allfiles) > 0:
                assert mice_n_files[i] == len(allfiles)
                names = run_names(mice_days[i])
                for j, file in enumerate(allfiles):
                    if j > 10:
                        name = names[j]
                        print(name + ' ' + str(j))
                        data = loadmat(file)['suite2pData']
                        # sesPoles=intarr(data.startIdx[1:])

                        mouse_id = i
                        run_id = j

                        # timeseries
                        # dffss = np.hsplit(data.dFF,sesPoles)
                        # spksss=np.hsplit(data.spks,sesPoles)
                        # might have behavior traces in other exps

                        neuronIdxs = np.where(data.iscell[:, 0].astype(bool))[0]
                        stats = [data.stat[i] for i in neuronIdxs]
                        snrs = data.snr
                        cellIdxs = intarr(data.cellIdx, mat2py=True)

                        run = db.run.find_one({"$and": [{"mouse_id": mouse_id}, {"run_id": run_id}]})
                        _run_id = run['_id']

                        _id_neurons = []
                        for neu_id, (stat, snr, neu_id_glob) in enumerate(zip(stats, snrs, cellIdxs)):
                            _id = "%01d%02d%01d%04d" % (mouse_id, run_id, 0, neu_id)
                            snr = snr.item()
                            try:
                                roi_pix_x = intarr(stat.xext)
                                roi_pix_y = intarr(stat.yext)
                            except:
                                roi_pix_x = intarr(stat.xcirc)
                                roi_pix_y = intarr(stat.ycirc)
                            roi_med = [int(stat.med[1]), int(stat.med[0])]
                            neuron = DBinterface(_id=_id, mouse_id=mouse_id, run_id=run_id,
                                                 snr=snr, neu_id=neu_id, neu_id_glob=neu_id_glob,
                                                 roi_pix_x=roi_pix_x, roi_pix_y=roi_pix_y, roi_med=roi_med,

                                                 )
                            try:
                                neuron.insert(db.neuron)
                            except:
                                pass


def set_grid():
    for ses in db.session.find({'is_stim': True}):
        last_stim = db.stim_event.find_one({'$and': [{'ses_id': ses['_id']}, {'stim_id_ses': ses['n_stim'] - 1}]})
        db.stim_event.update_one({'_id': last_stim['_id']}, {'$set': {'last_stim': True}})

    for ses in db.session.find({'is_stim': True}):
        dropped_last_trial = False
        stim_len_unified = ses['max_stim_len']

        _id = ses['_id']
        pre_on = ses['iti_min'] // 2
        post_on = ses['iti_min'] // 2 + ses['max_stim_len']
        ons = []
        offs = []

        stims = list(db.stim_event.find({'ses_id': ses['_id']}).sort("start_ses"))
        for stim in stims:
            on = stim['start_ses'] - pre_on;
            ons.append(on)
            off = stim['start_ses'] + post_on;
            offs.append(off)
        if off > ses['end']:
            stims.pop(-1);
            ons.pop(-1);
            offs.pop(-1)
            print(_id + ':' + 'Last trial exceeds total number of frames by %s frames and is flagged for exclusion' % (
                    off - ses['end']))
            dropped_last_trial = True

    for ses in db.session.find({'is_stim': True}):
        _id = ses['_id']
        pre_on = ses['iti_min'] // 2
        post_on = ses['iti_min'] // 2 + ses['max_stim_len']
        ons_reor = []
        offs_reor = []
        grid = db.grid.find_one({'_id': _id})

        for i, label in enumerate(ses['stim_labels']):
            stims = list(
                db.stim_event.find({'$and': [{'ses_id': ses['_id']}, {'label': label}]}).sort("stim_id_ori_ses"))
            for stim in stims:
                if grid['dropped_last_trial'] and ('last_stim' in stim):
                    ses['n_stim_ori'][i] = ses['n_stim_ori'][i] - 1
                    ses['n_stim'] = ses['n_stim'] - 1
                else:
                    on = stim['start_ses'] - pre_on;
                    ons_reor.append(on)
                    off = stim['start_ses'] + post_on;
                    offs_reor.append(off)
        db.grid.update_one({'_id': _id}, {'$set': {'ons_reor': ons_reor}})
        db.grid.update_one({'_id': _id}, {'$set': {'offs_reor': offs_reor}})

        row_poles = np.cumsum(ses['n_stim_ori'])[:-1].tolist()
        row_end = ses['n_stim']
        row_labels = ses['stim_labels']
        col_poles = [pre_on, pre_on + stim_len_unified]
        col_end = pre_on + post_on
        col_labels = ['stim_on', 'stim_off']

        grid = DBinterface(_id=_id, mouse_id=ses['mouse_id'], run_id=ses['run_id'], ses_id=ses['ses_id'],
                           ons=ons, offs=offs, ons_reor=ons_reor, offs_reor=offs_reor,
                           row_poles=row_poles, row_end=row_end, row_labels=row_labels,
                           col_poles=col_poles, col_end=col_end, col_labels=col_labels,
                           dropped_last_trial=dropped_last_trial)
        grid.insert(db.grid)


'''



def fp_crossday(config, mouse):
    d = loadmat(config['workdir'] + config['fp_crossday'] + '\\' + mouse['mouse_name'] + '.mat')
    for key, item in d.items():
        if key[0] != '_': return item


def set_fp_crossday():
    for mouse_name in mice:
        if mouse_name != 'Sut2':
            fp_crossday = 'C:\\Users\\selinali\\lab\\crossday_id\\' + mouse_name + 'CellIndex.mat'
            db.mouse.update_one({'mouse_name': mouse_name}, {"$set": {'fp_crossday': fp_crossday}})


def set_crossday_id():
    for mouse in db.mouse.find():
        ids_crossday = fp_crossday(config, mouse)
        for run in db.run.find({'mouse_id': mouse['_id']}):
            for neu in db.neuron.find_one({"_id": run_id}):


def set_session_plot_folder():
    config = db.config.find_one()
    fol = config['workdir'] + config['fp_trialsmat']
    for ses in db.session.find({'is_stim': True}, {'_id': 1}):
        os.mkdir(fol + '\\' + ses['_id'])
'''
def make_neuruns():
    for mouse in db.mouse.find():
        id_mouse=mouse['_id']
        for neu in db.neuron.find({'id_mouse':id_mouse}):
            id_neu=neu['_id']
            for id_run in mouse['run_id_stim']:
                DBinterface(_id='%s%d%s'%(id_neu[0:3],id_run,id_neu[4:]),
                            js=neu['js_%d'%id_run],
                            js_max_cond=neu['js_max_cond_%d'%id_run],
                            js_max_ori=neu['js_max_ori_%d' % id_run],
                            js_max=neu['js_max_%d' % id_run],
                            is_visdriven_cond=neu['is_visdriven_cond_%d' % id_run],
                            is_visdriven_ori=neu['is_visdriven_ori_%d' % id_run],
                            is_visdriven=neu['is_visdriven_%d' % id_run],
                            js_corrected=neu['js_corrected_%d' % id_run],
                            is_visdriven_mat=neu['is_visdriven_mat_%d' % id_run],
                            id_mouse=id_mouse,
                            id_ses=id_neu[:3],
                            id_run=id_run,
                            id_neu=neu['id_neu']
                            ).insert(db.neu_run)
'''
'''


def messed_up_again():#DONE: rename all files using _id not id_neu
    for ses in db.session.find({},{'fp_dff':1,'run_id_stim':1,'n_neu':1,'id_mouse':1,'name_mouse':1}):
        id_mouse=ses['id_mouse']
        for id_run in ses['run_id_stim']:
            _id_run = ses['_id']+str(id_run)
            run = db.run.find_one({'_id':_id_run},{'fp_trialsmat':1})
            fpout = "\\".join(run['fp_trialsmat'].split('\\')[:-1])+'\\'+ses['name_mouse']+'\\'+run['fp_trialsmat'].split('\\')[-1]
            #sf_id(db.run,_id_run,'fp_trialsmat',fpout)
            for id_neu in range(ses['n_neu']):
                _id_neu='%s%04d'%(_id_run, id_neu)
                fp=fpout+'\\'+str(id_neu)+'.jpg'
                fpnew=fpout+'\\'+_id_neu+'.jpg'
                os.rename(fp,fpnew)

def rename_preDB_batch_of_plots():
    #for file in glob.glob('C:\\Users\\selinali\\lab\\RDE20141\\2022-7-7-plots\\rois\\*'):
        #[name_mouse,id_ses,name_ses,date,_]=os.path.basename(file).split('_')
        #fp_new=os.path.dirname(file)+'\\'+'%s_%s'%(name_mouse,name_ses)
        #shutil.move(file,fp_new)
    for file in glob.glob('C:\\Users\\selinali\\lab\\RDE20141\\2022-7-7-plots\\rois\\*\\*.jpg'):
        [name_mouse,id_ses,name_ses,date,id_neu,_]=os.path.basename(file).split('_')
        mouse=db.mouse.find_one({'name':name_mouse},{'_id':1})
        fp_new=os.path.dirname(file)+'\\'+'%s%s0%04d.jpg'%(mouse['_id'],id_ses,int(id_neu[1:]))
        shutil.move(file,fp_new)

    for file in glob.glob('C:\\Users\\selinali\\lab\\RDE20141\\2022-7-7-plots\\trace\\*\\*.jpg'):
        [name_mouse,id_ses,name_ses,date,id_run,id_neu,_]=os.path.basename(file).split('_')
        mouse=db.mouse.find_one({'name':name_mouse},{'_id':1})
        fp_new=os.path.dirname(file)+'\\'+'%s%s%s%04d.jpg'%(mouse['_id'],id_ses,id_run[-1],int(id_neu[1:]))
        shutil.move(file,fp_new)

def tmp_mouse_run_id_stim():
    for mouse in db.mouse.find():
        ses=db.session.find_one({'id_mouse':mouse['_id']})
        sf_id(db.mouse,mouse['_id'],'run_id_stim',ses['run_id_stim'])

'''

def sort_folders():
    for file in glob.glob('C:\\Users\\selinali\\lab\\RDE20141\\2022-7-7-plots\\*trace'):
        shutil.move(file,'C:\\Users\\selinali\\lab\\RDE20141\\2022-7-7-plots\\trace\\'+os.path.basename(file))

def patch_1032():
    for mouse in gf(db.mouse,'name_sess'):
        i=mouse['_id']
        cond_sess=np.array([x[0] for x in mouse['name_sess']])
        _, cond_poles = np.unique(cond_sess, return_index=True)
        cond_poles.sort()
        conds=cond_sess[cond_poles]
        sf_id(db.mouse,i,'conds',conds.tolist())
        sf_id(db.mouse,i,'cond_poles',cond_poles[1:].tolist())

