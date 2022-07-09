import tkinter as tk
from tkinter import ttk
import numpy as np
import pickle

# for clicking days for each mice.
##this is ... soo messy but i don't wanna tidy it up

class Example(tk.LabelFrame):
    mice = ['T01', 'T02', 'T03', 'T04', 'Sut1', 'Sut2', 'Sut3', 'Sut4']
    mouse_cursor = 0
    hours = ['0025', '005', '01', '04', '06', '12', '24', '48', '72', '96', '120', '144', '168']
    condition = ['B', 'S', 'U', 'R']

    buttons = np.empty([len(condition), len(hours)],dtype=tk.Checkbutton)
    checkvars = np.empty([len(condition), len(hours)],dtype=tk.IntVar)
    data = np.zeros([len(mice),len(condition), len(hours)])

    def __init__(self, *args, **kwargs):
        try:
            with open('mice_days.pickle', 'rb') as file:
                self.data = pickle.load(file)
        except:
            pass

        tk.LabelFrame.__init__(self, *args, **kwargs)

        for i,cond in enumerate(self.condition):
            for j,hour in enumerate(self.hours):
                var = tk.IntVar()
                button = ttk.Checkbutton(self, variable=var, onvalue=1, offvalue=0, text=cond + hour, command=self.click(var))
                button.grid(row=i, column=j, sticky="ew")
                button.state(['!alternate'])
                val=int(self.data[self.mouse_cursor, i, j])
                if val==1:
                    button.state(['selected'])
                var.set(val)
                self.checkvars[i, j]=var
                self.buttons[i, j]=button
        root.title(self.mice[self.mouse_cursor])
        b = tk.Button(root, text="Submit", command=self.submit)
        b.pack()

    def click(self,var):
        var.set(1)

    def submit(self):
        if self.mouse_cursor == len(self.mice):
            print('end')
            self.save()
            return
        else:
            for i, a in enumerate(self.checkvars):
                for j, c in enumerate(a):
                    self.data[self.mouse_cursor, i, j] = c.get()
            print(self.data[self.mouse_cursor])
            self.mouse_cursor = self.mouse_cursor + 1
            if self.mouse_cursor != len(self.mice):
                root.title(self.mice[self.mouse_cursor])
                for i, a in enumerate(self.checkvars):
                    for j, c in enumerate(a):
                        self.checkvars[i, j].set(int(self.data[self.mouse_cursor, i, j]))

    def save(self):
        with open('mice_days.pickle', 'wb') as file:
            pickle.dump(self.data, file)


root = tk.Tk()
Example(root).pack(side="top", fill="both", expand=True, padx=10, pady=10)
root.mainloop()

with open('mice_days.pickle', 'rb') as file:
    data=pickle.load(file)

