# da tools

import matplotlib.pyplot as plt
import pymongo

################################## <MATLAB - PYTHON INTERFACES> #############################################
def loadmat(filename):
    import mat73
    import scipy.io
    try:
        data = mat73.loadmat(filename, use_attrdict=True)
    except:
        # not using this as for now.. instead converted non-7.3 files into v7.3 for sut3/4
        data = scipy.io.loadmat(filename, struct_as_record=True, squeeze_me=True)
    return data

def intarr(arr, mat2py=False):
    if mat2py:
        return [int(x - 1) for x in arr]
    return [int(x) for x in arr]

################################## <DATASET - DB INTERFACES> #############################################
class DBinterface:
    properties = None

    def __init__(self, iterable=(), **kwargs):
        self.__dict__.update(iterable, **kwargs)
        try:
            self.properties = kwargs
        except:
            pass

    def insert(self, col):
        return col.insert_one(self.properties)
    #obj.properties.update(key=val,...)


# where all the macro/exp level info go
def set_config(db):
    name_mice = ['T01', 'T02', 'T03', 'Sut1', 'Sut2', 'Sut3', 'Sut4']
    hours = ['0025', '005', '01', '04', '06', '12', '24', '48', '72', '96', '120', '144', '168']
    conditions = ['B', 'S', 'U', 'R']

    inputdir='C:\\Users\\selinali\\lab\\Mice'
    workdir='C:\\Users\\selinali\\lab\\sut'
    trialsmatdir='2022-7-22-plots\\trialsmat'

    stim_labels = [0, 45, 90, 135, 180, 225, 270, 315]
    stim_colors = ['#AC041B', '#F8693D', '#EDB120', '#21C40F', '#67DCD5', '#2A46FA', '#A81CBF', '#FF1493']

    select_mouse=False
    name_mouse_selected=''

    n_stim_ori_threshold=6
    js_threshold=[0.4]*3+[0.3]*4

    DBinterface(_id=0, name_mice=name_mice, hours=hours, conditions=conditions,
                inputdir=inputdir,workdir=workdir,trialsmatdir=trialsmatdir,
                stim_labels=stim_labels,stim_colors=stim_colors,
                select_mouse=select_mouse,name_mouse_selected=name_mouse_selected,
                n_stim_ori_threshold=n_stim_ori_threshold).insert(db.config)

def get_config(db):
    return DBinterface(db.config.find_one())

def ses_names(data,conditions,hours):
    sess=[]
    for data_col, cond in zip(data, conditions):
        for d, hour in zip(data_col, hours):
            if d == 1:
                sess.append(cond + hour)
    return sess

def insert(doc, col, _id, overwrite=True):
    try:
        doc.insert(col)
    except:
        if overwrite:
            col.delete_one({'_id': _id})
            doc.insert(col)
        else:
            print('%s: document already exists'% _id)


################################## <DB OPERATIONS> #############################################

def startMongo():
    return pymongo.MongoClient("mongodb://localhost:27017/").sut_mice

def backup(db):
    for each in db.stim_event.find():
        db.stim_event_backup.insert_one(each)
    for each in db.stim_event_backup.find():
        db.stim_event.insert_one(each)

    for each in db.neuron.find():
        db.z_neuron_backup.insert_one(each)
    for each in db.neuron_lite.find():
        db.neuron.insert_one(each)

    for each in db.backup.find({}):
        db.session.insert_one(each)
    for each in db.session.find():
        db.run.insert_one(each)

#find document by id
def fd_id(col, _id):
    return col.find_one({'_id':_id})

#delete field
def df(col, fieldname):
    col.update_many({}, {"$unset": {fieldname: 1}})

#set field
def sf(col, doc, fieldname, value):
    col.update_one({'_id':doc['_id']}, {"$set": {fieldname: value}})

#set field by id
def sf_id(col, _id, fn, val):
    col.update_one({'_id':_id}, {"$set": {fn: val}})

def sf_config(db,fieldname,value):
    db.config.update_one({},{'$set':{fieldname:value}})

#rename field
def rnf(col,fieldname,new):
    col.update_many({}, {"$rename": {fieldname: new}})

#clear all documents from collection
def clear_col(db):  # CAREFUL!!!!!!!!!!!!!!!!!!!!
    db.neuron_backup.delete_many({})
    db.stim_event_backup.delete_many({})
    db.session.delete_many({})

#partially clear the database.. for now
def clear_db(db):
    db.session.delete_many({})
    db.run.delete_many({})
    db.stim_event.delete_many({})

################################## <PLOTTING> #############################################
def savePlot(ax,filename):
    fig = ax.get_figure()
    fig.savefig(filename,bbox_inches='tight',pad_inches=0.2)
    fig.clear()
    plt.close(fig)
