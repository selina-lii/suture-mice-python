


def sf(col, doc, fn, val):
    col.update_one({'_id':doc['_id']}, {"$set": {fn: val}})

def sf_id(col, _id, fn, val):
    col.update_one({'_id':_id}, {"$set": {fn: val}})


def find(col,_id):
    return col.find_one({'_id':_id})







#API with file organziation

def fp_crossday(db):
    workdir = db.config.find_one()['workdir']
    col=db.mouse
    for doc in col.find():
        fp = workdir + '\\crossday_id\\' + doc['mouse_name'] + '.mat'
        sf(col=col,doc=doc,fn='fp_crossday',val=fp)
        print(fp)

def fp_dff(db):
    import glob
    workdir = db.config.find_one()['workdir']
    col=db.session
    for mouse in db.mouse.find():
        fps = glob.glob("%s\\dff\\%s*.mat" % (workdir, mouse['mouse_name']))
        for i,fp in enumerate(fps):
            _id="%d%02d"%(mouse['_id'],i)
            sf_id(col=col,_id=_id,fn='fp_dff',val=fp)
        print("%s\t\t#files:\t%d" %(mouse['mouse_name'],i+1))

def fp_refimg(db):
    import glob
    workdir = db.config.find_one()['workdir']
    col=db.session
    for mouse in db.mouse.find():
        fps = glob.glob("%s\\dff\\%s*.mat" % (workdir, mouse['mouse_name']))
        for i,fp in enumerate(fps):
            _id="%d%02d"%(mouse['_id'],i)
            sf_id(col=col,_id=_id,fn='fp_dff',val=fp)
        print("%s\t\t#files:\t%d" %(mouse['mouse_name'],i+1))

