import tkinter as tk
from tkinter import ttk
import numpy as np
from LyxTools import sf_id
import pymongo

# for clicking days for each mice.

class Example(tk.LabelFrame):
    mice_names = ['T01', 'T02', 'T03', 'Sut1', 'Sut2', 'Sut3', 'Sut4']
    hours = ['0025', '005', '01', '04', '06', '08', '12', '24', '48', '72', '96', '120', '144', '168']
    conditions = ['B', 'S', 'U', 'R']
    mouse_cursor = 0
    buttons = np.empty([len(conditions), len(hours)],dtype=tk.Checkbutton)
    checkvars = np.empty([len(conditions), len(hours)],dtype=tk.IntVar)
    data = np.zeros([len(mice_names), len(conditions), len(hours)], dtype='bool')

    def __init__(self, *args, **kwargs):
        try:
            self.data=np.asarray(db.config.find_one({},{'gui_mice_days':1})['gui_mice_days'])
            #Use when modifying hours or conditions
            '''            mice=list(db.mouse.find({},{'name_sess':1}).sort('_id'))
            for mouse in mice:
                id = mouse['_id']
                print(id)
                name_sess=mouse['name_sess']
                for name in name_sess:
                    cond=name[0]
                    hour=name[1:]
                    self.data[id,self.conditions.index(cond),self.hours.index(hour)]=True'''
        except:
            pass

        tk.LabelFrame.__init__(self, *args, **kwargs)

        for i,cond in enumerate(self.conditions):
            for j,hour in enumerate(self.hours):
                var = tk.IntVar()
                button = ttk.Checkbutton(self, variable=var, onvalue=1, offvalue=0, text=cond + hour, command=self.click(var))
                button.grid(row=i, column=j, sticky="ew")
                button.state(['!alternate'])
                val=int(self.data[self.mouse_cursor, i, j])
                if va    l == 1:
    button.state(['selected'])


var.set(val)
self.checkvars[i, j] = var
self.buttons[i, j] = button
root.title(self.mice_names[self.mouse_cursor])
b = tk.Button(root, text="Submit", command=self.submit)
b.pack()


def click(self,var):
        var.set(1)

    def submit(self):
        if self.mouse_cursor == len(self.mice_names):
            print('end')
            self.save()
            root.destroy()
            return
        else:
            for i, a in enumerate(self.checkvars):
                for j, c in enumerate(a):
                    self.data[self.mouse_cursor, i, j] = c.get()
            print(self.data[self.mouse_cursor])
            self.mouse_cursor = self.mouse_cursor + 1
            if self.mouse_cursor != len(self.mice_names):
                root.title(self.mice_names[self.mouse_cursor])
                for i, a in enumerate(self.checkvars):
                    for j, c in enumerate(a):
                        self.checkvars[i, j].set(int(self.data[self.mouse_cursor, i, j]))

    def save(self):
        db.config.update_one({},{'$set':{'gui_mice_days':self.data.tolist()}})
        for mouse in db.mouse.find({}, {'_id':1}):
            name_sess=[]
            data=self.data[mouse['_id']]
            for data_row, cond in zip(data, self.conditions):
                for d, hour in zip(data_row, self.hours):
                    if d==True:
                        name_sess.append(cond + hour)
            sf_id(db.mouse,mouse['_id'],'name_sess',name_sess)

db = pymongo.MongoClient("mongodb://localhost:27017/").sut_mice2
root = tk.Tk()
Example(root).pack(side="top", fill="both", expand=True, padx=10, pady=10)
root.mainloop()
