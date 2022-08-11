
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

                if add_mice:
                    conds = np.array([x[0] for x in mouse['name_sess']])
                    _, cond_poles = np.unique(conds, return_index=True)
                    cond_poles.sort()
                    conds = conds[cond_poles]

                    mouse = DBinterface(_id=id_mouse,
                                        name=name_mouse,
                                        n_sess=n_files,
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
