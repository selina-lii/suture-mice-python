# da tools

import matplotlib.pyplot as plt
import pymongo
import os
import mat73
################################## <MATLAB - PYTHON INTERFACES> #############################################
import scipy.io


def loadmat(filename):
    try:
        data = mat73.loadmat(filename, use_attrdict=True)
        return data
    except:
        # not using this as for now.. instead converted non-7.3 files into v7.3 for sut3/4
        return scipy.io.loadmat(filename)

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

    def insert(self, col, overwrite=True):
        try:
            col.insert_one(self.properties)
        except:
            if overwrite:
                col.delete_one({'_id': self._id})
                col.insert_one(self.properties)
            else:
                print('%s: document already exists' % self._id)


def get_config(db):
    return DBinterface(db.config.find_one())

def get_name_tag(id_mouse,name_mouse,name_ses,fp_suite2pData):
    fn_suite2pData = fp_suite2pData.split('\\')[-1].split('_')
    cond_code = name_ses[0]
    hour_code = name_ses[1:]
    date_of_acq = fn_suite2pData[1]
    nametag = DBinterface(id_mouse=id_mouse, name_mouse=name_mouse,
                           cond_code=cond_code, hour_code=hour_code,date_of_acq=date_of_acq)
    return nametag

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
def sf(col, doc, dict):
    col.update_one({'_id':doc['_id']}, {"$set": dict})

#set field by id
def sf_id(col, _id, dict):
    col.update_one({'_id':_id}, {"$set": dict})

def sf_ids(col, dict_ids, dict):
    col.update_many({'_id':{'$in':dict_ids}}, {"$set": dict})

def sf_config(db,fieldname,value):
    db.config.update_one({},{'$set':{fieldname:value}})

#loop through one field in all documents in a collection
def gf(col,fieldname):
    return col.find({},{fieldname:1})

#rename field
def rnf(col,fieldname,new):
    col.update_many({}, {"$rename": {fieldname: new}})

#clear all documents from collection
def clear_col(db):  # CAREFUL!!!!!!!!!!!!!!!!!!!!
    db.neu_run_spont.delete_many({})
    db.stim_event_backup.delete_many({})
    db.session.delete_many({})

#partially clear the database.. for now
def clear_db(db):
    db.session.delete_many({})
    db.run.delete_many({})
    db.stim_event.delete_many({})

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

################################## <PLOTTING> #############################################
def savePlot(ax,filename):
    fig = ax.get_figure()
    fig.savefig(filename,bbox_inches='tight',pad_inches=0.2)
    fig.clear()
    plt.close(fig)

def savePlot_fig(fig,filename):
    fig.savefig(filename,bbox_inches='tight',pad_inches=0.2)
    fig.clear()
    plt.close(fig)


def mkdir(fp):
    if not os.path.isdir(fp):
        os.mkdir(fp)

