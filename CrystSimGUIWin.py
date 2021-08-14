#!/usr/bin/env python
# coding: utf-8

# In[1]:


from IPython import get_ipython
import numpy as np

#import matplotlib
from matplotlib import pyplot as plt
import sys
import os
#import math
from scipy.optimize import least_squares, curve_fit
#matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure

# append pcrystal2 dir to path, if not in the same directory
# sys.path.append('')
import pcrystal2, pcrystal2_mod


# In[2]:


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)




class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        #Display text in tooltip window
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        tip_label = tk.Label(tw, text=self.text, justify="left",
                      background="#ffffe0", borderwidth=1,
                      font=("tahoma", "12", "normal"), fg='black')
        tip_label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)
    

    
    
def DIYmdestroy():
    global DIYmWindow, firstError
    
    DIYmWindow.grab_release()
    if firstError == True:
        if DIYmWindow.winfo_exists():
            DIYmWindow.destroy()
            
    
def DIYmessage(titletext, text):
    global DIYmWindow, firstError
    
    if firstError == True:
        if DIYmWindow.winfo_exists():
            DIYmWindow.destroy()
    
    DIYmWindow = tk.Toplevel()
    DIYmWindow.title(titletext)
    DIYmWindow.geometry("450x150")
    firstError = True
    DIYmWindow.grab_set()
    DIYmframe=tk.Frame(DIYmWindow, bg="#d9dcde", bd=5)
    DIYmframe.place(relx=0.5, rely=0, relwidth=0.99, relheight=0.99, anchor="n")
    DIYmlabel=tk.Label(DIYmframe, text=text, bg="#d9dcde", bd=10, font=("Times", 14), fg='black', anchor='n')
    DIYmlabel.place(relx=0.01, rely=0, relwidth=0.99, relheight=0.5)
    
    DIYmOKbutton = tk.Button(DIYmframe, text="OK", bg="#778899", font=("Times", 16), command=lambda: DIYmdestroy(), fg='black')
    DIYmOKbutton.place(relx=0.375, rely=0.675, relwidth=0.25, relheight=0.25)
    
    
    
    


# In[3]:


from math import pi, exp
import scipy as sp

import mie_aniso

def validate(string): #checks for only float values
    regex = re.compile(r"(\+|\-)?.[0-9.]*$")
    result = regex.match(string)
    return (string == ""
             or (string.count('+') <= 1
                and string.count('-') <= 1
                and string.count('.') <= 1
                and result is not None
                and result.group(0) != ""))
    
def on_validate(P):
    return validate(P)    


def I_trans_tot(alpha_b, L, Ibe = 1):
    """
    alpha_b: backward extinction coefficient
    calculates total transmission with backward extinction
    """
    return Ibe/(np.exp(-2*alpha_b*L) + 1)*2*np.exp(-2*alpha_b*L)
    
def I_trans_spc(alpha_b, alpha_f, L, Ibe = 1):
    return Ibe*np.exp(-L*(alpha_b+alpha_f))

def I_trans_dif(alpha_b, alpha_f, L, Ibe = 1):
    return I_trans_tot(alpha_b, L, Ibe) - I_trans_spc(alpha_b, alpha_f, L, Ibe)

def haze(alpha_b, alpha_f, L):
    return I_trans_dif(alpha_b, alpha_f, L)/I_trans_tot(alpha_b, L)

def test_function(entry,entry2,entry3,entry4,entry5): 
    global Hlabel6, Hlabel7, hazeWindow
    
    if float(entry2)<=0:
        DIYmessage("Error!", "Spherulite size is invalid")
    else:
        r0 = float(entry2) # radius in micron
        dn = float(entry5) # Birefringence
        n_bg = float(entry4) # background refractive index
        LL = float(entry) # sample thickness in micron
        wl = float(entry3) # wavelength in micron
    
        Hlabel6["text"]= "Sample thickness [um]: %s \nSpherulite radius [um]:  %s \nWavelength [um]: %s \nRefractive index [-]:  %s \nBirefringence [-]: %s" % (LL, "{:<10.4g}".format(r0), wl, n_bg, dn)
    
        #print("Sample thickness:", LL)
        #print("Spherulite radius:", r0)
        #print("Wavelength:", wl)
        #print("Refractive index:", n_bg)
        #print("Birefringence:", dn)
    
        n_r = n_bg+dn*0.5 # radial refractive index
        n_tt = n_bg-dn*0.5 # tangential refractive index
    
        f_V = 1 # volume fraction (sum v_i)/V_container
        n0 = f_V/(4./3*r0**3*pi) # nucleus density N/V [1/um^3]
    
    # create single mie object
        sm = mie_aniso.single_mie(r0, wl, n_r**2, n_tt**2, n_bg**2, method='my')
    
    # calculate forward and backward cross section
        (C_f, C_b) = sm.csca_fb()
    
    # total transmission
        Hlabel7["text"]="Total transmission: %s \nDiffuse transmission: %s \nSpecular transmission: %s \nHaze: %s" % ("{:<10.4g}".format(I_trans_tot(C_b*n0, LL)),"{:<10.4g}".format(I_trans_dif(C_b*n0, C_f*n0, LL)),"{:<10.4g}".format(I_trans_spc(C_b*n0, C_f*n0, LL)), "{:<10.4g}".format(haze(C_b*n0, C_f*n0, LL)))
        """print("total transmission:","{:<10.4g}".format(I_trans_tot(C_b*n0, LL)))
        print("diffuse transmission:", round(I_trans_dif(C_b*n0, C_f*n0, LL),4))
        print("specular transmission:", round(I_trans_spc(C_b*n0, C_f*n0, LL),4))
        print("haze:", round(haze(C_b*n0, C_f*n0, LL),4))"""
    
def HazeSim(avgsph_in):
    global Hlabel6, Hlabel7, hazeWindow, Hframe, hazeWindowfirstopen
 
    HEIGHT=1000
    WIDTH=1400
    
    if hazeWindowfirstopen == True:
        if hazeWindow.winfo_exists():
            hazeWindow.destroy()
    
    
    hazeWindow = tk.Toplevel()
    hazeWindow.title("Haze simulation for semicrystalline polymers")
    hazeWindow.geometry("1400x1000")
    hazeWindowfirstopen = True
    
    Hcanvas=tk.Canvas(hazeWindow, height=HEIGHT, width=WIDTH)
    Hcanvas.pack()

    Hframe=tk.Frame(hazeWindow, bg="#7be9ed", bd=5)
    Hframe.place(relx=0.5, rely=0, relwidth=1, relheight=1, anchor="n")

    image_path = resource_path("JM_JoPS_cover_dark_.gif")
    """Hbackground_image=tk.PhotoImage(file = "/home/aliz/Crystal/fcrystal/JM_JoPS_cover_dark_.png")"""
    Hbackground_image=tk.PhotoImage(file = image_path)
    Hbackground_label = tk.Label(Hframe, image=Hbackground_image)
    Hbackground_label.place(x=0, y=0, relwidth=1, relheight=1)
    
    Hlabel=tk.Label(Hframe, text="Haze simulation for semicrystalline polymers", bg="#c4cdcf", fg='black', bd=10, font=("Times", 36))
    Hlabel.place(relx=0.125, rely=0, relwidth=0.75, relheight=0.1)

    Hlabel1=tk.Label(Hframe, text="Sample thickness [um]", bg="#c4cdcf", font=("Times", 16), fg='black')
    Hlabel1.place(relx=0.02, rely=0.2, relwidth=0.15, relheight=0.075)
            
    Hentry=tk.Entry(Hframe, font=("Times", 16), justify="center", bd = 5,validate="key", bg='white', fg='black')
    Hentry.place(relx=0.2, rely=0.2, relwidth=0.1, relheight=0.075)
    vcmd = (Hentry.register(on_validate), '%P') 
    Hentry.config(validatecommand=vcmd)
    Hentry.insert(0, "1000")

    Hlabel2=tk.Label(Hframe, text="Spherulite radius [um]", bg="#c4cdcf", font=("Times", 16), fg='black')
    Hlabel2.place(relx=0.02, rely=0.31, relwidth=0.15, relheight=0.075)

    Hentry2=tk.Entry(Hframe, font=("Times", 16), justify="center", bd = 5,validate="key", bg='white', fg='black')
    Hentry2.place(relx=0.2, rely=0.31, relwidth=0.1, relheight=0.075)
    vcmd = (Hentry2.register(on_validate), '%P') 
    Hentry2.config(validatecommand=vcmd)
    Hentry2.insert(0,round(avgsph_in/2,3))

    Hlabel3=tk.Label(Hframe, text="Wavelength [um]", bg="#c4cdcf", font=("Times", 16), fg='black')
    Hlabel3.place(relx=0.02, rely=0.42, relwidth=0.15, relheight=0.075)
    
    Hentry3=tk.Entry(Hframe, font=("Times", 16), justify="center", bd = 5,validate="key", bg='white', fg='black')
    Hentry3.place(relx=0.2, rely=0.42, relwidth=0.1, relheight=0.075)
    vcmd = (Hentry3.register(on_validate), '%P') 
    Hentry3.config(validatecommand=vcmd)
    Hentry3.insert(0, "0.6328")
            
    Hlabel4=tk.Label(Hframe, text="Refractive index [-]", bg="#c4cdcf",font=("Times", 16), bd = 5, fg='black')
    Hlabel4.place(relx=0.02, rely=0.53, relwidth=0.15, relheight=0.075)
            
    Hentry4=tk.Entry(Hframe, font=("Times", 16), justify="center", bd = 5,validate="key", bg='white', fg='black')
    Hentry4.place(relx=0.2, rely=0.53, relwidth=0.1, relheight=0.075)
    vcmd = (Hentry4.register(on_validate), '%P') 
    Hentry4.config(validatecommand=vcmd)
    Hentry4.insert(0, "1.5")

    Hlabel5=tk.Label(Hframe, text="Birefringence [-]", bg="#c4cdcf",font=("Times", 16), bd = 5, fg='black')
    Hlabel5.place(relx=0.02, rely=0.64, relwidth=0.15, relheight=0.075)
            
    Hentry5=tk.Entry(Hframe, font=("Times", 16), justify="center", bd = 5,validate="key", bg='white', fg='black')
    Hentry5.place(relx=0.2, rely=0.64, relwidth=0.1, relheight=0.075)
    vcmd = (Hentry5.register(on_validate), '%P') 
    Hentry5.config(validatecommand=vcmd)

    Hentry5.insert(0, "0.0027")
    
    Hside_frame=tk.Frame(hazeWindow, bg="#96d9e3")
    Hside_frame.place(relx=0.5, rely=0.2, relwidth=0.35, relheight=0.515, anchor="n")

    Hlabel6=tk.Label(Hside_frame, font=("Times", 22), anchor="nw", justify="left", bd=4, fg='black', bg="#c4cdcf")
    Hlabel6.place(relwidth=1, relheight=0.5)

    Hlabel7=tk.Label(Hside_frame, font=("Times", 22), anchor="nw",bg="#7be9ed", justify="left", bd=4, fg='black')
    Hlabel7.place(relx=0, rely=0.5, relwidth=1, relheight=0.5)
                                                                                                                                       
    Hbutton = tk.Button(Hframe, text="Simulate haze", bg="#7be9ed", font=("Times", 24), command=lambda: test_function(Hentry.get(), Hentry2.get(), Hentry3.get(), Hentry4.get(), Hentry5.get()))
    Hbutton.place(relx=0.02, rely=0.8, relwidth=0.20, relheight=0.05)

    
    

    hazeWindow.mainloop()


# In[4]:


def HL_func(T, G0, KG):
    #calculate G in the entire T(array) temperature range
    #returns G_fit, array of growth speeds

    U=755 #activation energy/R
    Tg=263.15 #glass transition temperature
    Tref=30.0 #reference temperature
    T0m=481.15 #equilibrium melting temperature

    G_fit=[]
    #G_fit=G0*exp(-U/(T-Tg+Tref))*exp(KG*(T0m**2)*(T0m+T)/(2*(T**2)*(T0m-T)))

    """ezt amúgy át lehetne írni a fortranos H-L függvényre
    vagy azt átírni, hogy kezeljen listát :'('"""
    
    for i in range(np.size(T)):
        G = G0*exp(-U/(T[i]-Tg+Tref))*exp(KG*(T0m**2)*(T0m+T[i])/(2*(T[i]**2)*(T0m-T[i])))
        G_fit.append(G)
    return G_fit

def Avrami_func(t, n, k):
    
    x_t=[]
    
    for i in range(np.size(t)):
        x_aktualis = 1-exp(-k*(t[i]**n))
        x_t.append(x_aktualis)
    return x_t

def DFresetArea():
    global DFbutton_loadtxt, DFframe, isTxt, DF_Ttext_area, DF_Gtext_area
    global DF_T_entry1, DF_T_entry2, DF_T_entry3, DF_T_entry4, DF_T_entry5
    global DF_G_entry1, DF_G_entry2, DF_G_entry3, DF_G_entry4, DF_G_entry5
    
    DF_Ttext_area.destroy()
    DF_Gtext_area.destroy()
    
    DF_T_entry1 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_T_entry1.place(relx=0.1, rely=0.175, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_T_entry1.insert(0, 0)
    
    DF_T_entry2 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_T_entry2.place(relx=0.1, rely=0.3, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_T_entry2.insert(0, 0)
    
    DF_T_entry3 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_T_entry3.place(relx=0.1, rely=0.425, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_T_entry3.insert(0, 0)
    
    DF_T_entry4 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_T_entry4.place(relx=0.1, rely=0.55, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_T_entry4.insert(0, 0)
    
    DF_T_entry5 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_T_entry5.place(relx=0.1, rely=0.675, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_T_entry5.insert(0, 0)
    
    DF_G_entry1 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_G_entry1.place(relx=0.55, rely=0.175, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_G_entry1.insert(0, 0)
    
    DF_G_entry2 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_G_entry2.place(relx=0.55, rely=0.3, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_G_entry2.insert(0, 0)
    
    DF_G_entry3 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_G_entry3.place(relx=0.55, rely=0.425, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_G_entry3.insert(0, 0)
    
    DF_G_entry4 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_G_entry4.place(relx=0.55, rely=0.55, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_G_entry4.insert(0, 0)
    
    DF_G_entry5 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_G_entry5.place(relx=0.55, rely=0.675, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_G_entry5.insert(0, 0)
    
    DFbutton_loadtxt.configure(command=lambda: DFloadtxt(), text="Load data from .txt file")
    
    isTxt=False

def DFloadtxt():
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    global DFbutton_loadtxt, DFframe, isTxt, DF_Ttext_area, DF_Gtext_area, T_arr, G_arr
    global DF_T_entry1, DF_T_entry2, DF_T_entry3, DF_T_entry4, DF_T_entry5
    global DF_G_entry1, DF_G_entry2, DF_G_entry3, DF_G_entry4, DF_G_entry5
    
    T_arr=[0]
    G_arr=[0]
    DFfilename = askopenfilename()
    try:
        try:
            data_input = np.loadtxt(DFfilename,skiprows=1)
        except (ValueError, TypeError) as ee:
            DIYmessage("Error!", "Cannot read some data from text file")
        else:
            try:
                T_arr=data_input[:,0]
                G_arr=data_input[:,1]
            except IndexError:
                DIYmessage("Error!", "Cannot read data from text file")
            else:
                DF_T_entry1.destroy()
                DF_T_entry2.destroy()
                DF_T_entry3.destroy()
                DF_T_entry4.destroy()
                DF_T_entry5.destroy()
                DF_G_entry1.destroy()
                DF_G_entry2.destroy()
                DF_G_entry3.destroy()
                DF_G_entry4.destroy()
                DF_G_entry5.destroy()
                DF_Ttext_area = tk.Text(DFframe, width = 35, height = 65, font = ("Times New Roman", 15), bg='white', fg='black') 
                DF_Ttext_area.place(relx=0.1, rely=0.175, relwidth=0.35, relheight=0.55, anchor='nw')
                DF_Gtext_area = tk.Text(DFframe, width = 35, height = 65, font = ("Times New Roman", 15), bg='white', fg='black') 
                DF_Gtext_area.place(relx=0.55, rely=0.175, relwidth=0.35, relheight=0.55, anchor='nw')
                for i in range(np.size(T_arr)):
                    DF_Ttext_area.insert(tk.INSERT,str(T_arr[i])+' \n')
                    T_arr[i] = T_arr[i]+273.15
                for i in range(np.size(G_arr)):
                    DF_Gtext_area.insert(tk.INSERT,str(G_arr[i])+' \n')
                DFbutton_loadtxt.configure(command=lambda: DFresetArea(), text="Type data manually")
                DF_Ttext_area.configure(state='disabled')
                DF_Gtext_area.configure(state='disabled')
                isTxt = True
    except (UnicodeDecodeError, OSError) as eeee:
        DIYmessage("Error!", "Text file is invalid")
        
            
            
    
    
def dataFit():
    global HL_params, DFlabel_results, isTxt, T_arr, G_arr
    global DF_T_entry1, DF_T_entry2, DF_T_entry3, DF_T_entry4, DF_T_entry5
    global DF_G_entry1, DF_G_entry2, DF_G_entry3, DF_G_entry4, DF_G_entry5
    global Tmin, Tmax
    
    p_0 = (20277.70942,-2.67238)
    if not isTxt: #echeck if every value is a valid number
        T_arr=[0]
        G_arr=[0]
        try:
            G_arr = [float(DF_G_entry1.get()),float(DF_G_entry2.get()),float(DF_G_entry3.get()),float(DF_G_entry4.get()),float(DF_G_entry5.get())]
        except (ValueError, TypeError) as ee:
            G_arr = [-1]
        else:
            try:
                T_arr = [float(DF_T_entry1.get())+273.15,float(DF_T_entry2.get())+273.15,float(DF_T_entry3.get())+273.15,float(DF_T_entry4.get())+273.15,float(DF_T_entry5.get())+273.15]
            except (ValueError, TypeError) as ee:
                T_arr = [-1]
            else:
                k=0
                while (k<=np.size(T_arr)-1)&(k<=np.size(G_arr)-1): #a lista utolsó eleme a listaméret-1-edik :'(((((
                    if (T_arr[k]==273.15) or (G_arr[k]==0):
                        T_arr.pop(k)
                        G_arr.pop(k)
                    else:
                        k=k+1
    else:
        try:
            for i in range(np.size(T_arr)):
                #print(T_arr[i])
                T_arr[i]=float(T_arr[i])
        except (ValueError, TypeError) as ee:
            T_arr=[-1]
        try:
            for i in range(np.size(G_arr)):
                G_arr[i]=float(G_arr[i])
        except (ValueError, TypeError) as ee:
            G_arr=[-1]

    if (all(350<i<460 for i in T_arr)):
        if (all(k>0 for k in G_arr)):
            try:
                fitted = curve_fit(HL_func, T_arr, G_arr, p0=p_0) #determine the iterative parameters
                HL_params[0] = fitted[0][0]
                HL_params[1] = fitted[0][1]
            #determine Tmax/Tmin
                go = True
                Tmax = 423.15
                while go:
                    Tmax = Tmax + 0.05
                    if pcrystal2_mod.fcrystal.hl_calc(Tmax, HL_params)<1e-10:
                        go = False
                go = True
                Tmin = 313.15
                while go:
                    Tmin = Tmin - 0.05
                    if pcrystal2_mod.fcrystal.hl_calc(Tmin, HL_params)<1e-10:
                        go = False
                DFlabel_results["text"] ="The fitted parameters: \nG0: %s \nKG: %s \n Valid temperature range:\n %s°C - %s°C" % ("{:<10.4g}".
                    format(HL_params[0]), "{:<10.4g}".format(HL_params[1]), "{:<10.4g}".format(Tmin-273.15), "{:<10.4g}".format(Tmax-273.15))
                #DFlabel_results["text"] ="Valid temperature range:\n %s°C - %s°C" % ("{:<10.4g}".format(Tmin-273.15), "{:<10.4g}".format(Tmax-273.15))
            except (ValueError, RuntimeError, OverflowError) as e:
                DIYmessage("Error!", "Fitting has failed")
        else:
            DIYmessage("Error!", "Invalid spherulite growth speed")
    else:
        DIYmessage("Error!", "Temperature too high, too low or invalid")

    
    
    
def dataFit_window():
    global HL_params, DFlabel_results, DFbutton_loadtxt, DFframe, isTxt, dataFitfirstopen, datafitWindow
    global DF_T_entry1, DF_T_entry2, DF_T_entry3, DF_T_entry4, DF_T_entry5
    global DF_G_entry1, DF_G_entry2, DF_G_entry3, DF_G_entry4, DF_G_entry5
    
    
    isTxt=False
    
    if dataFitfirstopen == True:
        if datafitWindow.winfo_exists():
            datafitWindow.destroy()
    
    datafitWindow = tk.Toplevel()
    datafitWindow.title("Fit experimental data")
    datafitWindow.geometry("700x500")
    dataFitfirstopen = True
    DFframe=tk.Frame(datafitWindow, bg="#c28285", bd=5)
    DFframe.place(relx=0.5, rely=0, relwidth=0.95, relheight=0.95, anchor="n")
    
    DFlabel_T=tk.Label(DFframe, text="Temperatures [°C]", bg="#c4cdcf", bd=10, font=("Times", 12), fg='black')
    DFlabel_T.place(relx=0.1, rely=0.05, relwidth=0.35, relheight=0.075)
    DFlabel_G=tk.Label(DFframe, text="Spherulite growth speeds [m/s]", bg="#c4cdcf", bd=10, font=("Times", 12), fg='black')
    DFlabel_G.place(relx=0.55, rely=0.05, relwidth=0.35, relheight=0.075)
    
    DF_T_entry1 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_T_entry1.place(relx=0.1, rely=0.175, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_T_entry1.insert(0, 0)
    
    DF_T_entry2 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_T_entry2.place(relx=0.1, rely=0.3, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_T_entry2.insert(0, 0)
    
    DF_T_entry3 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_T_entry3.place(relx=0.1, rely=0.425, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_T_entry3.insert(0, 0)
    
    DF_T_entry4 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_T_entry4.place(relx=0.1, rely=0.55, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_T_entry4.insert(0, 0)
    
    DF_T_entry5 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_T_entry5.place(relx=0.1, rely=0.675, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_T_entry5.insert(0, 0)
    
    DF_G_entry1 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_G_entry1.place(relx=0.55, rely=0.175, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_G_entry1.insert(0, 0)
    
    DF_G_entry2 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_G_entry2.place(relx=0.55, rely=0.3, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_G_entry2.insert(0, 0)
    
    DF_G_entry3 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_G_entry3.place(relx=0.55, rely=0.425, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_G_entry3.insert(0, 0)
    
    DF_G_entry4 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_G_entry4.place(relx=0.55, rely=0.55, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_G_entry4.insert(0, 0)
    
    DF_G_entry5 = tk.Entry(DFframe, validate="key", bg='white', fg='black')
    DF_G_entry5.place(relx=0.55, rely=0.675, relwidth=0.35, relheight=0.075, anchor='nw')
    DF_G_entry5.insert(0, 0)
    

    
    DFbutton_calculate = tk.Button(DFframe, text="Calculate", bg="#778899", font=("Times", 16), fg='black', command=lambda: dataFit())
    DFbutton_calculate.place(relx=0.1, rely=0.795, relwidth=0.35, relheight=0.075)
    
    DFbutton_loadtxt = tk.Button(DFframe, text="Load data from .txt file", bg="#778899", font=("Times", 16), fg='black', command=lambda: DFloadtxt())
    DFbutton_loadtxt.place(relx=0.1, rely=0.9, relwidth=0.35, relheight=0.075)
    
    DFlabel_results=tk.Label(DFframe, text="", bg="#d3d3d3", bd=10, font=("Times", 12), fg='black')
    DFlabel_results.place(relx=0.55, rely=0.795, relwidth=0.35, relheight=0.18)
    
    
    
    datafitWindow.mainloop()


# In[5]:


import tkinter as tk
from tkinter import font
import tkinter.scrolledtext as st
import re
from datetime import datetime
from scipy.ndimage import gaussian_filter1d



HEIGHT = 800
WIDTH = 1620



def residual_array(params, t_exp, k_exp, G, plot = False, verbose = 1, min_vf = 0.05, max_vf = 0.99):
    """
    The program calculates the volume fraction curve vs. time for the given parameters.
    Plots the experimental data with the calculated one between min_vf and max_vf (in case of exp. data)
      with the time axis aligned so that the curves are at the same position at 50% volume fraction
    Calculates the difference at each time-point between the calculated and experimental data.
    
    params: [log10(rho_inst), log10(F)], where 
            - rho_inst is the density of instantaneously born seeds [1/um^3]
            - F is the birthspeed of thermal seeds [1/(um^3 * s)]
    t_exp:  time data of the DSC measurement in seconds
    k_exp:  volume fraction data of the DSC measurement. 0 < k_exp < 100
    G:      growth speed [um/s], constant currently
    
    returns the difference vector between experimental and calculated volume fraction
    """
    global n_t, tmaxC 
    global diameter_out, filename_out, dt_string #txt
    global frame1, frame2, frame3, text_area, newWindowtext, text_result #widget
    global isSeparate #status
    #ami kell hogy global legyen mert használja: MINDEN AMIBE KIÍR, diameter_out hogy számolhasson vele, n_t tesztszám miatt, filename!
    
    

    (t1, k1, vols, N_thrm, N_inst, VV) = pcrystal2_mod.fcrystal_cdf(10**params[0], 10**params[1], G, tmaxC, verbose = (verbose > 1), n_test = n_t)
    inds = np.logical_and(k_exp > min_vf*100, k_exp < max_vf*100)
    
    t1_center = t1[np.nonzero(k1>=0.5)[0][0]]
    texp_center = t_exp[np.nonzero(k_exp>=50)[0][0]]
    y = np.interp(t_exp[inds]-texp_center, t1-t1_center, k1)

    vols[:,0]=2*((vols[:,0]*3./4/np.pi)**(1./3.))

    

    
    
    if verbose>2:
        np.savetxt(filename_out+'_gocok'+dt_string+'.txt', vols, delimiter=",")
        np.savetxt(filename_out+'_szimkonv'+dt_string+'.txt', np.c_[t1, k1])
        one_line="Total nucleus density [1/m^3]: "+"{:<10.4g}".format((10**params[0]+N_thrm/VV)*10**18)+' '+"{:<10.2%}".format(N_inst/float(N_thrm+N_inst))+'/ '+"{:<10.2%}".format(N_thrm/float(N_thrm+N_inst))+'\n'
        with open(filename_out+'_szimkonv'+dt_string+'.txt', 'r+') as fp:
            lines = fp.readlines()     # lines is list of line, each element '...\n'
            lines.insert(0, one_line)  # you can use any index if you know the line index
            fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
            fp.writelines(lines)       # write whole lists again to the same file

          

    if plot:
        frame1.destroy()
        frame1=tk.Frame(frame, bg="#ffffff", bd=5)
        frame1.place(relx=0.32, rely=0.15, relwidth=0.33, relheight=0.55, anchor="nw")
        fig = Figure(figsize = (5, 5), dpi = 100) 
        plot1 = fig.add_subplot(111) 
        plot1.plot(t_exp[inds],k_exp[inds]/100,label = 'measured')
        plot1.plot(t_exp[inds],y,label = 'simulated')
        plot1.legend()
        plot1.set_xlabel('Time [s]')
        plot1.set_ylabel('Conversion [-]')
        plot1.set_title('Conversion curves') #kirajzolja a konverziógörbét

        canvas = FigureCanvasTkAgg(fig, master = frame1)   #draw the plot
        canvas.draw() 
        canvas.get_tk_widget().place(relx=0, rely=0, relwidth=0.75, relheight=0.75) 
  
        toolbar = NavigationToolbar2Tk(canvas, frame1) #add the toolbar to the plot
        toolbar.update() 
        canvas.get_tk_widget().pack()
 
       
    if verbose > 0:
        text_area.configure(state='normal')
        #text_area.update()
        text_area.insert(tk.INSERT,"rho_i                     rho_t                         F                       G                    f_i                         f_t         cost"+'\n') 
        #text_area.insert(tk.INSERT,"{:<10.4g}{:<10.4g}{:<10.4g}{:<10.4g}{:<10.2%}{:<10.2%}{:<10.4g}".
        #      format(10**params[0], N_thrm/VV, 10**params[1], G, N_inst/float(N_thrm+N_inst), N_thrm/float(N_thrm+N_inst),0.5*np.sum((k_exp[inds]/100. - y)**2)))
        text_area.insert(tk.INSERT,"{:<10.4g}".format(10**params[0]*10**18)) 
        text_area.insert(tk.INSERT,"             ")
        text_area.insert(tk.INSERT,"{:<10.4g}".format(N_thrm/VV*10**18)) 
        text_area.insert(tk.INSERT,"             ")  
        text_area.insert(tk.INSERT,"{:<10.4g}".format(10**params[1]*10**18)) 
        text_area.insert(tk.INSERT,"            ")
        text_area.insert(tk.INSERT,"{:<10.4g}".format(G*10**-6)) 
        text_area.insert(tk.INSERT,"           ")
        text_area.insert(tk.INSERT,"{:<10.2%}".format(N_inst/float(N_thrm+N_inst))) 
        text_area.insert(tk.INSERT,"           ")
        text_area.insert(tk.INSERT,"{:<10.2%}".format(N_thrm/float(N_thrm+N_inst))) 
        text_area.insert(tk.INSERT,"        ")
        text_area.insert(tk.INSERT,"{:<10.4g}".format(0.5*np.sum((k_exp[inds]/100. - y)**2))) 
        text_area.insert(tk.INSERT,"     ")
        text_area.insert(tk.INSERT,' \n') #kiírja az új ablakba az iterációt (nem folyamatos????)
        text_area.configure(state='disabled')
        
        
    if verbose == 3:
        text_result.configure(state='normal')
      #  text_result.insert(tk.INSERT,"rho_i     rho_t     F         G         f_i       f_t       cost"+'\n') 
      #  text_result.insert(tk.INSERT,"{:<10.4g}{:<10.4g}{:<10.4g}{:<10.4g}{:<10.2%}{:<10.2%}{:<10.4g}".
      #        format(10**params[0], N_thrm/VV, 10**params[1], G, N_inst/float(N_thrm+N_inst), N_thrm/float(N_thrm+N_inst),0.5*np.sum((k_exp[inds]/100. - y)**2)))
        text_result.insert(tk.INSERT,"Total nucleus density [1/m^3]: ") 
        text_result.insert(tk.INSERT,"{:<10.4g}".format((10**params[0]+N_thrm/VV)*10**18)) 
        text_result.insert(tk.INSERT,' \n')
        text_result.insert(tk.INSERT,"Density of instantenously born seeds [1/m^3]: ") 
        text_result.insert(tk.INSERT,"{:<10.4g}".format(10**params[0]*10**18)) 
        text_result.insert(tk.INSERT,' \n')
        if N_thrm/VV*10**18>=1:
            text_result.insert(tk.INSERT,"Density of thermal seeds [1/m^3]: ") 
            text_result.insert(tk.INSERT,"{:<10.4g}".format(N_thrm/VV*10**18)) 
            text_result.insert(tk.INSERT,' \n')
        else:
            text_result.insert(tk.INSERT,"Density of thermal seeds [1/m^3]: ") 
            text_result.insert(tk.INSERT,"0") 
            text_result.insert(tk.INSERT,' \n')
        text_result.insert(tk.INSERT,"Birthspeed of thermal seeds [1/m^3s]: ") 
        text_result.insert(tk.INSERT,"{:<10.4g}".format(10**params[1]*10**18)) 
        text_result.insert(tk.INSERT,' \n')
        text_result.insert(tk.INSERT,"% of instantanenously born seeds: ") 
        text_result.insert(tk.INSERT,"{:<10.2%}".format(N_inst/float(N_thrm+N_inst))) 
        text_result.insert(tk.INSERT,'   ')
        text_result.insert(tk.INSERT,"% of thermal seeds: ") 
        text_result.insert(tk.INSERT,"{:<10.2%}".format(N_thrm/float(N_thrm+N_inst))) 
        text_result.insert(tk.INSERT,' \n')
        text_result.configure(state='disabled')
        #kiírja a fő ablakba a legjobb szimulációt
        
    
    return k_exp[inds]/100. - y

def residual_array_aniso(params, t_exp, k_exp, T_start, coolingrate, HLparams, Tcp, plot = False, verbose = 1, min_vf = 0.05, max_vf = 0.99):
    """
    The program calculates the volume fraction curve vs. time for the given parameters.
    Plots the experimental data with the calculated one between min_vf and max_vf (in case of exp. data)
      with the time axis aligned so that the curves are at the same position at 50% volume fraction
    Calculates the difference at each time-point between the calculated and experimental data.
    
    params: [log10(rho_inst), log10(F)], where 
            - rho_inst is the density of instantaneously born seeds [1/um^3]
            - F is the average birthspeed of continuously born seeds [1/(um^3 * s)]
    t_exp:  time data of the DSC measurement in seconds
    k_exp:  volume fraction data of the DSC measurement. 0 < k_exp < 100
    T_start: starting temperature of the crystallization [K]
    coolingrate: cooling rate [K/s]
    HLparams: fitted parameters (array of 2)
    
    returns the difference vector between experimental and calculated volume fraction
    """
    
    
    global n_t, tmaxC 
    global diameter_out, filename_out, dt_string #txt
    global frame1, frame2, frame3, text_area, newWindowtext, text_result #widget
    global isSeparate #status
    
    (t1, k1, vols, N_thrm, N_inst, VV) = pcrystal2_mod.fcrystal_cdf_aniso(10**params[0], 10**params[1], T_start, coolingrate, HLparams, tmaxC, Tcp, verbose = (verbose > 1), n_test = n_t)
    inds = np.logical_and(k_exp > min_vf*100, k_exp < max_vf*100)
    
    t1_center = t1[np.nonzero(k1>=0.5)[0][0]]
    texp_center = t_exp[np.nonzero(k_exp>=50)[0][0]]


    y = np.interp(t_exp[inds]-texp_center, t1-t1_center, k1)
    vols[:,0]=2*((vols[:,0]*3./4/np.pi)**(1./3.))
    

    
 
    if verbose>2:
        
        np.savetxt(filename_out+'_gocok'+dt_string+'.txt', vols, delimiter=",")
        np.savetxt(filename_out+'_szimkonv'+dt_string+'.txt', np.c_[t1, k1])
        one_line="Total nucleus density [1/m^3]: "+"{:<10.4g}".format((10**params[0]+N_thrm/VV)*10**18)+' '+"{:<10.2%}".format(N_inst/float(N_thrm+N_inst))+'/ '+"{:<10.2%}".format(N_thrm/float(N_thrm+N_inst))+'\n'
        with open(filename_out+'_szimkonv'+dt_string+'.txt', 'r+') as fp:
            lines = fp.readlines()     # lines is list of line, each element '...\n'
            lines.insert(0, one_line)  # you can use any index if you know the line index
            fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
            fp.writelines(lines)       # write whole lists again to the same file

            
            
            
            
    if verbose > 0:
        #text_area.update()
        text_area.configure(state='normal')
        text_area.insert(tk.INSERT,"rho_i                     rho_t                         F                    f_i                         f_t         cost"+'\n') 
        #text_area.insert(tk.INSERT,"{:<10.4g}{:<10.4g}{:<10.4g}{:<10.4g}{:<10.2%}{:<10.2%}{:<10.4g}".
        #      format(10**params[0], N_thrm/VV, 10**params[1], G, N_inst/float(N_thrm+N_inst), N_thrm/float(N_thrm+N_inst),0.5*np.sum((k_exp[inds]/100. - y)**2)))
        text_area.insert(tk.INSERT,"{:<10.4g}".format(10**params[0]*10**18)) 
        text_area.insert(tk.INSERT,"             ")
        text_area.insert(tk.INSERT,"{:<10.4g}".format(N_thrm/VV*10**18)) 
        text_area.insert(tk.INSERT,"             ")  
        text_area.insert(tk.INSERT,"{:<10.4g}".format(10**params[1]*10**18)) 
        text_area.insert(tk.INSERT,"            ")
       # text_area.insert(tk.INSERT,"{:<10.4g}".format(G*10**-6)) 
       # text_area.insert(tk.INSERT,"           ")
        text_area.insert(tk.INSERT,"{:<10.2%}".format(N_inst/float(N_thrm+N_inst))) 
        text_area.insert(tk.INSERT,"           ")
        text_area.insert(tk.INSERT,"{:<10.2%}".format(N_thrm/float(N_thrm+N_inst))) 
        text_area.insert(tk.INSERT,"        ")
        text_area.insert(tk.INSERT,"{:<10.4g}".format(0.5*np.sum((k_exp[inds]/100. - y)**2))) 
        text_area.insert(tk.INSERT,"     ")
        text_area.insert(tk.INSERT,' \n') #kiírja az új ablakba az iterációt (nem folyamatos????)
        text_area.configure(state='disabled')
        
        
    if verbose == 3:
      #  text_result.insert(tk.INSERT,"rho_i     rho_t     F         G         f_i       f_t       cost"+'\n') 
      #  text_result.insert(tk.INSERT,"{:<10.4g}{:<10.4g}{:<10.4g}{:<10.4g}{:<10.2%}{:<10.2%}{:<10.4g}".
      #        format(10**params[0], N_thrm/VV, 10**params[1], G, N_inst/float(N_thrm+N_inst), N_thrm/float(N_thrm+N_inst),0.5*np.sum((k_exp[inds]/100. - y)**2)))

        text_result.configure(state='normal')
        text_result.insert(tk.INSERT,"Total nucleus density [1/m^3]: ") 
        text_result.insert(tk.INSERT,"{:<10.4g}".format((10**params[0]+N_thrm/VV)*10**18)) 
        text_result.insert(tk.INSERT,' \n')
        text_result.insert(tk.INSERT,"Density of instantenously born seeds [1/m^3]: ") 
        text_result.insert(tk.INSERT,"{:<10.4g}".format(10**params[0]*10**18)) 
        text_result.insert(tk.INSERT,' \n')
        if N_thrm/VV*10**18>=1:
            text_result.insert(tk.INSERT,"Density of continuously born seeds [1/m^3]: ") 
            text_result.insert(tk.INSERT,"{:<10.4g}".format(N_thrm/VV*10**18)) 
            text_result.insert(tk.INSERT,' \n')
        else:
            text_result.insert(tk.INSERT,"Density of continuously born seeds [1/m^3]: ") 
            text_result.insert(tk.INSERT,"0") 
            text_result.insert(tk.INSERT,' \n')
        text_result.insert(tk.INSERT,"Birthspeed of continuously born seeds [1/m^3s]: ") 
        text_result.insert(tk.INSERT,"{:<10.4g}".format(10**params[1]*10**18)) 
        text_result.insert(tk.INSERT,' \n')
        text_result.insert(tk.INSERT,"% of instantanenously born seeds: ") 
        text_result.insert(tk.INSERT,"{:<10.2%}".format(N_inst/float(N_thrm+N_inst))) 
        text_result.insert(tk.INSERT,'   ')
        text_result.insert(tk.INSERT,"% of continuously born seeds: ") 
        text_result.insert(tk.INSERT,"{:<10.2%}".format(N_thrm/float(N_thrm+N_inst))) 
        text_result.insert(tk.INSERT,' \n')
        text_result.configure(state='disabled')
        #kiírja a fő ablakba a legjobb szimulációt  
        
        
        frame1.destroy()
        frame1=tk.Frame(frame, bg="#ffffff", bd=5)
        frame1.place(relx=0.32, rely=0.15, relwidth=0.33, relheight=0.55, anchor="nw")
        
        fig = Figure(figsize = (5, 5), dpi = 100) 
        plot1 = fig.add_subplot(111) 
        plot1.plot(t_exp[inds],k_exp[inds]/100,label = 'measured')
        plot1.plot(t_exp[inds],y,label = 'simulated')
        plot3 = plot1.twinx()
        plot3.plot(t_exp[inds],T_start-t_exp[inds]*coolingrate-273.15,label = 'temperature', color='orange')
        plot1.legend(loc="right")
        plot3.legend(loc="lower right")
        plot1.set_xlabel('Time [s]')
        plot1.set_ylabel('Conversion [-]')
        plot3.set_ylabel('Temperature [°C]')
        plot1.set_title('Conversion curves') #kirajzolja a konverziógörbét
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master = frame1)   #draw the plot
        canvas.draw() 
        canvas.get_tk_widget().place(relx=0, rely=0, relwidth=0.75, relheight=0.75) 

        toolbar = NavigationToolbar2Tk(canvas, frame1) #add the toolbar to the plot
        toolbar.update() 
        canvas.get_tk_widget().pack()
    

    return k_exp[inds]/100. - y








def openNewWindow():
    global newWindowtext, text_area, isfirstrun, frame1, frame2, isAniso

    if isfirstrun:
        text_area.configure(state='normal')
        text_area.delete(1.0,tk.END)
    else:
        newWindowtext = tk.Toplevel()
        newWindowtext.title("the current iteration")
        newWindowtext.geometry("920x300")
        text_area = st.ScrolledText(newWindowtext, width = 100, height = 10, font = ("Times New Roman", 15), bg='white', fg='black') 
    
        text_area.grid(column = 0, pady = 10, padx = 10) #megnyitja az iterációt kiíró ablakot ha még nincs
    if isAniso.get()==2:
        text_area.insert(tk.INSERT,"Density of             Density of           Birthspeed         % of               % of               cost"+'\n')
        text_area.insert(tk.INSERT,"instantenously         continuous         of continuous       instantenously  continuous"+'\n')
        text_area.insert(tk.INSERT,"born seeds [1/m^3]    seeds [1/m^3]  seeds [1/m^3s]    born seeds         seeds"+'\n')
    else:
        text_area.insert(tk.INSERT,"Density of             Density of           Birthspeed         Spherulite         % of               % of               cost"+'\n')
        text_area.insert(tk.INSERT,"instantenously         thermal            of thermal         growth speed  instantenously  thermal"+'\n')
        text_area.insert(tk.INSERT,"born seeds [1/m^3]    seeds [1/m^3]  seeds [1/m^3s]    [m/s]             born seeds         seeds"+'\n')
    #text_area.insert(tk.INSERT,"rho_i                     rho_t                         F                       G                    f_i                         f_t         cost"+'\n')
    text_area.update()
    if isfirstrun: #diagramok kiírása ha nem első lefutás
        frame1.destroy()
        frame2.destroy()
        frame1=tk.Frame(frame, bg="#ffffff", bd=5)
        frame1.place(relx=0.32, rely=0.15, relwidth=0.33, relheight=0.55, anchor="nw")
        frame2=tk.Frame(frame, bg="#ffffff", bd=5)
        frame2.place(relx=0.67, rely=0.15, relwidth=0.325, relheight=0.55, anchor="nw")
        
        
    root.update()
    
def percentage_seek(perc, t_temp, k_temp):
    #find the first instance of a certain conversion%
    tmax_seek=t_temp[-1]
    i=0
    while 1:
        if ((k_temp[i])<perc)&((k_temp[i+1])>=perc):
            tmax_seek=t_temp[i]*60.
            break
        else:
            i=i+1
            if i==(np.size(t_temp)-1):
                  break
    return tmax_seek

def Tcp_seek(t_temp, k_temp, T_start, cr, tmax, tmin):
    from scipy import interpolate
    global Tcp_default
    k_smooth = gaussian_filter1d(k_temp, 100)
    k_d2 = np.gradient(np.gradient(k_smooth))
    infls = np.where(np.diff(np.sign(k_d2)))[0]
    #print(infls)
    inds = np.logical_and(infls>tmin,infls<tmax)
    infls = infls[inds]
    #print(infls)
    """inds = np.logical_and((T_start-t_temp[infls]*cr)>Tmin+273.15, (T_start-t_temp[infls]*cr)<Tmax+273.15)
    Tcp_ido = infls[inds]"""
    Tcp = T_start-infls/60*cr
    
    
    if np.size(Tcp)!=1:
        tk_interp_func = interpolate.interp1d(t_temp, k_temp, kind='cubic')
        t_temp_uj = np.linspace(t_temp[0], t_temp[-1], num = np.size(t_temp)*10, endpoint=True)
        k_temp_uj = tk_interp_func(t_temp_uj)
        k_smooth = gaussian_filter1d(k_temp, 100)
        k_d2 = np.gradient(np.gradient(k_smooth))
        infls = np.where(np.diff(np.sign(k_d2)))[0]
        inds = np.logical_and(infls>tmin,infls<tmax)
        infls = infls[inds]
        Tcp = T_start-infls/60*cr
        
        if np.size(Tcp)!=1: 
        
            answer = tk.messagebox.askokcancel("Error!"," Unable to determine peak crystallization temperature\n Press OK to proceed with default Tcp ("+str(round(Tcp_default-273.15,1))+" °C)")
            if answer:
                Tcp = Tcp_default
            else:
                raise ValueError
        else:
            Tcp = float(Tcp)
    else:
        Tcp = float(Tcp)
    return Tcp


def plot_conversioninput(T_st, cr):
    global conversion_input, frame1, isTemp, isAniso
    t_temp = conversion_input[:,0]
    k_temp = conversion_input[:,1]
    inds = np.logical_and(k_temp > 0.05*100, k_temp < 0.99*100)
        
    frame1.destroy()
    frame1=tk.Frame(frame, bg="#ffffff", bd=5)
    frame1.place(relx=0.32, rely=0.15, relwidth=0.33, relheight=0.55, anchor="nw")
        
    fig = Figure(figsize = (5, 5), dpi = 100) 
    plot1 = fig.add_subplot(111) 
    plot1.plot(t_temp[inds]*60,k_temp[inds]/100,label = 'measured')
    plot3 = plot1.twinx()
    if isAniso.get()==2:
        plot3.plot(t_temp[inds]*60,T_st-t_temp[inds]*cr-273.15,label = 'temperature', color='orange')
        plot3.legend(loc="lower right")
        plot3.set_ylabel('Temperature [°C]')
        
    plot1.legend(loc="right")
    plot1.set_xlabel('Time [s]')
    plot1.set_ylabel('Conversion [-]')

    plot1.set_title('Conversion curves') #kirajzolja a konverziógörbét
    fig.tight_layout()
        
    canvas = FigureCanvasTkAgg(fig, master = frame1)  
    canvas.draw() 
    canvas.get_tk_widget().place(relx=0, rely=0, relwidth=0.75, relheight=0.75) 

    toolbar = NavigationToolbar2Tk(canvas, frame1) 
    toolbar.update() 
    canvas.get_tk_widget().pack()

def browse():
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    global conversion_input, filename, filename_out
    global frame1, cr_entry, Tst_entry
    global Tmin, Tmax, T_st, tmaxC, tminC, Tcp
    global isTemp, isSmooth, isAniso
    conversion_input=np.zeros((1,2)) #biztonság kedvéért kinullázza a régi fájlt
    filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
       
    try:
        conversion_input = np.loadtxt(filename,skiprows=3)
        filename_out = filename[:-4]
    except (UnicodeDecodeError, ValueError) as e:
        DIYmessage("Error!", "DSC file is invalid")
        filename='no text'
        
    else:
        label_file=tk.Label(frame, text="The conversion curve file:", bg="#d9dcde", bd=10, font=("Times", 12), fg='black')
        label_file.place(relx=0.05, rely=0.815, relwidth=0.25, relheight=0.05)    
        label_file2=tk.Label(frame, text=filename_out, bg="#d9dcde", bd=10, font=("Times", 8), fg='black')
        label_file2.place(relx=0.05, rely=0.865, relwidth=0.25, relheight=0.05)

        #print(isTemp, isAniso.get())
        if isTemp&(isAniso.get()==2):
            t_temp = conversion_input[:,0] + 273.15
            k_temp = conversion_input[:,1]
            
            ind_limit = np.logical_and(t_temp<=Tmax, t_temp>=Tmin)
            t_temp = t_temp[ind_limit]
            k_temp = k_temp[ind_limit]
                    
            try:
                cr=float(cr_entry.get())
            except ValueError:
                DIYmessage("Error!", "Cooling rate is invalid")
                filename='no text'
            else:
                T_st = t_temp[0]
                for i in range(np.size(t_temp)):
                    t_temp[i] = (T_st-t_temp[i])/cr
                
                tmaxC = percentage_seek(100.,t_temp,k_temp)
                tminC = percentage_seek(0.0001,t_temp,k_temp)
                
                try:
                    Tcp = Tcp_seek(t_temp, k_temp, T_st, cr, tmaxC, tminC)
                except ValueError:
                    filename='no text'
                else:
                    conversion_input = np.zeros((np.size(t_temp),2))
                    conversion_input[:,0] = t_temp
                    conversion_input[:,1] = k_temp
                    plot_conversioninput(T_st, cr)
                
                    isSmooth = False
                
                    CreateToolTip(label_file2, text = ('Tcp = '+str(round(Tcp-273.15, 1))+'°C'))
                
                
        elif (not isTemp)&(isAniso.get()==2):            
            try:
                cr = float(cr_entry.get())
                T_st = float(Tst_entry.get())+273.15
            except ValueError:
                DIYmessage("Error!", "Cooling rate or starting temperature is invalid")
                filename='no text'
            else:
                tmaxC = percentage_seek(100.,conversion_input[:,0],conversion_input[:,1])
                tminC = percentage_seek(0.0001,conversion_input[:,0],conversion_input[:,1])
                try:
                    Tcp = Tcp_seek(conversion_input[:,0], conversion_input[:,1], T_st, cr, tmaxC, tminC)
                except ValueError:
                    filename='no text'
                else:
                    plot_conversioninput(T_st, cr)
                    isSmooth = False
                    CreateToolTip(label_file2, text = ('Tcp = '+str(round(Tcp-273.15, 1))+'°C'))

        elif isAniso.get()==1:
            
            T_st = 1
            cr = 1 #ezek csak azért, hogy a plot_conversioninput lefusson, nem számol velük
            plot_conversioninput(T_st, cr)
            tmaxC = percentage_seek(100.,conversion_input[:,0],conversion_input[:,1])
            tminC = percentage_seek(0.0001,conversion_input[:,0],conversion_input[:,1])
            isSmooth = False

        

        
        
def plot_sph(file_name):
    global isSeparate, phantom_counter, dt_string, frame2
    sph_true_list_inst=[]
    sph_true_list_sp=[]
    sph_true_list=[]
    sph_list=[]
    
    frame2.destroy()
    frame2=tk.Frame(frame, bg="#ffffff", bd=5)
    frame2.place(relx=0.67, rely=0.15, relwidth=0.325, relheight=0.55, anchor="nw")
    
    if (not file_name == 'no text')&isfirstrun:
        try:
            file_name_out = file_name[:-4]
            sph_list=np.loadtxt(file_name_out+'_gocok'+dt_string+'.txt', delimiter=",", skiprows=1)
            sph_true_counter=0

            phantom_counter=0
            if isSeparate.get()==True:
                for k in range(int(np.size(sph_list)/2)): #kiválasztja a nem 0 szferolitokat
                    if sph_list[k,0] != 0:
                        if sph_list[k,1]==0: #születési idő szerint válogat
                            sph_true_list_inst.append(sph_list[k,0])
                        else:
                            sph_true_list_sp.append(sph_list[k,0])
                        sph_true_counter = sph_true_counter+1
                    else:
                        phantom_counter = phantom_counter+1
                avgsph=(sum(sph_true_list_inst)+sum(sph_true_list_sp))/sph_true_counter
            else:
                for k in range(int(np.size(sph_list)/2)): 
                    if sph_list[k,0] != 0:
                        sph_true_list.append(sph_list[k,0])
                        sph_true_counter = sph_true_counter+1
                    else:
                        phantom_counter = phantom_counter+1      
                avgsph=sum(sph_true_list)/sph_true_counter





            fig2 = Figure(figsize = (5, 5), dpi = 100) #kirajzolja a méreteloszlást
            plot2 = fig2.add_subplot(111) 
            if isSeparate.get()==True:
                plot2.hist([sph_true_list_inst, sph_true_list_sp], bins=100, color = ["orange", "skyblue"], label = ['instant', 'continuous'], density=False) #csak a nem 0 szferolitméreteket ábrázolja
                plot2.legend()
            else:
                plot2.hist(sph_true_list, bins=100, color = (0,0,1,0.5)) 
            plot2.set_xlabel('Spherulite diameter [$\mu$m]')
            plot2.set_ylabel('Count')
            plot2.set_title('Spherulite size distribution')

            canvas = FigureCanvasTkAgg(fig2, master = frame2)   #draw the plot
            canvas.draw() 
            canvas.get_tk_widget().place(relx=0, rely=0, relwidth=0.75, relheight=0.75) 

            toolbar = NavigationToolbar2Tk(canvas, frame2) #add the toolbar to the plot
            toolbar.update() 
            canvas.get_tk_widget().pack()

            frame2.update()
            return avgsph
        except OSError:
            DIYmessage("Error!", "Spherulite distribution data not found!")
    elif (file_name == 'no text')&isfirstrun:
        DIYmessage("Error!", "Spherulite distribution data not found!")
    
def simulate(p0, G):
    global d_s, conversion_input, filename, filename_out
    global text_result, label_avgsph2, text_area, frame2
    global sph_true_list, phantom_counter, avg_sph 
    global isfirstrun, dt_string, isSeparate
    
    
    openNewWindow()
    if isfirstrun:
        text_result.configure(state='normal')
        text_result.delete(1.0,tk.END)
    text_result.update()
    t_exp = conversion_input[:,0]*60. #másodperc
    k_exp = conversion_input[:,1]

    
    res_lsq = least_squares(residual_array, p0, args=(t_exp, k_exp, G), bounds=([-9,-np.inf],[1,1]), verbose = 2, method = 'trf', diff_step = d_s)
    text_area.configure(state='normal')
    text_area.insert(tk.INSERT,' \n')  
    text_area.configure(state='disabled')
    #residual_array(res_lsq.x, t_exp, k_exp, G, plot = False, verbose = 4, min_vf = 0.05, max_vf = 0.99)
    residual_array(res_lsq.x, t_exp, k_exp, G, plot = True, verbose = 3, min_vf = 0.05, max_vf = 0.99) 

    text_result.configure(state='normal')
    text_result.insert(tk.INSERT,'Cost: ')
    text_result.insert(tk.INSERT,"{:<10.4g}".format(res_lsq.cost))
    text_result.insert(tk.INSERT,' \n')
    text_result.insert(tk.INSERT,'Number of evaluations performed: ')
    text_result.insert(tk.INSERT,"{:<10.4g}".format(res_lsq.nfev))
    text_result.configure(state='disabled')
    
    one_line='Cost: '+"{:<10.4g}".format(res_lsq.cost)+'\n'
    with open(filename_out+'_szimkonv'+dt_string+'.txt', 'r+') as fp:
        lines = fp.readlines()     # lines is list of line, each element '...\n'
        lines.insert(1, one_line)  # you can use any index if you know the line index
        fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
        fp.writelines(lines)       # write whole lists again to the same file
    

    if not isfirstrun: #ez végtelenül suta itt de az isfirstrun első szimulációnál még false de már kéne a plot
        isfirstrun = True
        avg_sph = plot_sph(filename)
        isfirstrun = False
    else:
        avg_sph = plot_sph(filename)
    label_avgsph2.config(text=str(round(avg_sph,3)))
    
    one_line='Average spherulite size: '+"{:<10.4g}".format(avg_sph)+'\n'
    with open(filename_out+'_gocok'+dt_string+'.txt', 'r+') as fp:
        lines = fp.readlines()     # lines is list of line, each element '...\n'
        lines.insert(0, one_line)  # you can use any index if you know the line index
        fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
        fp.writelines(lines)       # write whole lists again to the same file

    
    return res_lsq.cost
    
def simulate_aniso(p0, HLparams, T_start, cr):
    global d_s, conversion_input, filename, Tcp
    global text_result, label_avgsph2, text_area, label_avgsph3, frame2
    global sph_true_list, phantom_counter, avg_sph 
    global isfirstrun, dt_string
    openNewWindow()
    if isfirstrun:
        text_result.configure(state='normal')
        text_result.delete(1.0,tk.END)
    text_result.update()
    t_exp = conversion_input[:,0]*60.
    k_exp = conversion_input[:,1]

    res_lsq = least_squares(residual_array_aniso, p0, args=(t_exp, k_exp, T_start, cr, HLparams, Tcp), bounds=([-9,-np.inf],[1,1]), verbose = 2, method = 'trf', diff_step = d_s)
    text_area.configure(state='normal')
    text_area.insert(tk.INSERT,' \n') 
    text_area.configure(state='disabled')
    #residual_array_aniso(res_lsq.x, t_exp, k_exp, T_start, cr, HLparams, plot = False, verbose = 4, min_vf = 0.05, max_vf = 0.99)
    residual_array_aniso(res_lsq.x, t_exp, k_exp, T_start, cr, HLparams, Tcp, plot = True, verbose = 3, min_vf = 0.05, max_vf = 0.99)

    text_result.configure(state='normal')
    text_result.insert(tk.INSERT,'Cost: ')
    text_result.insert(tk.INSERT,"{:<10.4g}".format(res_lsq.cost))
    text_result.insert(tk.INSERT,' \n')
    text_result.insert(tk.INSERT,'Number of evaluations performed: ')
    text_result.insert(tk.INSERT,"{:<10.4g}".format(res_lsq.nfev))
    text_result.configure(state='disabled')
    
    one_line='Cost: '+"{:<10.4g}".format(res_lsq.cost)+'\n'
    with open(filename_out+'_szimkonv'+dt_string+'.txt', 'r+') as fp:
        lines = fp.readlines()     # lines is list of line, each element '...\n'
        lines.insert(1, one_line)  # you can use any index if you know the line index
        fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
        fp.writelines(lines)       # write whole lists again to the same file
    
    
    
    #ˇˇˇ ez mind izoterm/anizotermtől független a végéig
    if not isfirstrun: 
        isfirstrun = True
        avg_sph = plot_sph(filename)
        isfirstrun = False
    else:
        avg_sph = plot_sph(filename)
    label_avgsph2.config(text=str(round(avg_sph,3)))
    
    one_line='Average spherulite size: '+"{:<10.4g}".format(avg_sph)+'\n'
    with open(filename_out+'_gocok'+dt_string+'.txt', 'r+') as fp:
        lines = fp.readlines()     # lines is list of line, each element '...\n'
        lines.insert(0, one_line)  # you can use any index if you know the line index
        fp.seek(0)                 # file pointer locates at the beginning to write the whole file again
        fp.writelines(lines)       # write whole lists again to the same file
    
    return res_lsq.cost
    

    
def simulate_start():
    global n_t, T_st, phantom_counter
    global buttons, F_entry, G_entry, rhoi_entry, nt_entry
    global isfirstrun, G_setting_mode, isAniso, custom_nt 
    global dt_string, HL_params
    buttons.config(text="In progress...")
    buttons.update()
    
    cr = 5
    T_start = 408
    G_temp= 3.277e-2
    
    now = datetime.now()
    dt_string = now.strftime("_%d-%m-%Y %H_%M_%S")
    try:
        F_szim = float(F_entry.get())-18.
        rhoi_szim = float(rhoi_entry.get())-18. 
        parameter_szim = [rhoi_szim, F_szim]
    except ValueError:
        F_szim=200
        rhoi_szim=200 #ha F vagy rhoi nem szám, akkor megfogja
        
    
    if isAniso.get()==1:    
        try:
            if G_setting_mode:
                G_temp = (float(G_entry.get()))*10**6
            else:
                temp = float(G_entry.get())+273.15 #celsius-kelvin
                if (temp>Tmin)&(temp<Tmax):
                    G_temp=pcrystal2_mod.fcrystal.hl_calc(temp, HL_params)*10**6 #program mikrométerben számol!!!! még mindig!!!!!
                else:
                    #tk.messagebox.showerror("Error!", "Temperature is invalid")
                    G_temp=-1
        except ValueError:
            G_temp=-1 #ha G nem szám akkor ez megfogja 
    else:
        try:
            if isTemp:
                T_start = T_st
            else:
                T_start = float(Tst_entry.get())+273.15 #celsius-kelvin
            if (T_start>Tmin)&(T_start<Tmax):
                G_temp=1
            else:
                G_temp=-1
            cr = float(cr_entry.get())/60. #felsőalsókorlát?, másodperc
        except ValueError:
            G_temp=-1 #this is to indicate if any parameters are invalid, not only G
        
    if custom_nt:
        try:
            n_t=float(nt_entry.get())
        except ValueError:
            n_t=1
    if filename != 'no text':
        if (n_t<5*10**8) and (n_t>10): #test point 
            if is_integer(n_t): 
                if G_temp>0: 
                    if (rhoi_szim>=-9)&(rhoi_szim<1)&(F_szim<1): #ez a 30 ez nagyon felülről közelítve lamellavastagságból
                        try:
                            if isComplex:
                                complex_simulation(parameter_szim, G_temp, HL_params, T_start, cr)
                            else:
                                if isAniso.get()==1:
                                    simulate(parameter_szim, G_temp)
                                else:
                                    simulate_aniso(parameter_szim, HL_params, T_start, cr)
                            isfirstrun=True
                        except IndexError:
                            DIYmessage("Error!", "Input file is not a valid conversion curve")
                            isfirstrun=True
                    else:
                        DIYmessage("Error!", "Density or birthspeed of seeds is invalid")
                else:
                    DIYmessage("Error!", "Temperature, cooling rate or spherulite growth speed is invalid")
            else:
                DIYmessage("Error!", "Number of test points must be an integer")
        else:
            DIYmessage("Error!", "Number of test points is too high or too low")
            n_t=1e5
        
    else:
        DIYmessage("Error!", "No conversion curve found")
    if (not isComplex)&(phantom_counter>1000):
        DIYmessage("Warning!", "The list of simulated spherulites contains phantom spherulites.\nDue to this, the size and size distribution data may be inaccurate.")
    buttons.config(text="Simulate")


def change_G(a): 
    global G_setting_mode 
    if a==1:
        mb.config(text='Temperature [°C]')
        G_entry.delete(0, 'end')
        G_entry.insert(0, 135)
        G_setting_mode=False
    else:
        mb.config(text='Spherulite growth speed [m/s]')
        G_entry.delete(0, 'end')
        G_entry.insert(0, 0.00000003277)
        G_setting_mode=True
        
        
def resetArea(a):
    global mb, G_entry, Tst_entry, label_cr, label_Tst, cr_entry, label_conv, mb4
    if a==2:
        
        if not label_cr.winfo_exists():
            mb.destroy()
            G_entry.destroy()
            label_conv.destroy()
        
            mb4=tk.Menubutton(frame, text="Conversion curve", font=("Times", 16), bg="#a3b5bf", fg='black')
            mb4.place(relx=0.05, rely=0.22, relwidth=0.12, relheight=0.05)
            mb4.menu =  tk.Menu(mb4, tearoff = 0)
            mb4["menu"] =  mb4.menu
            mb4.menu.add_command(label=str('time-%'),command = lambda:set_isTemp(1) )
            mb4.menu.add_command(label=str('temperature-% (set cooling rate first!)'),command = lambda:set_isTemp(2) )
            CreateToolTip(mb4, text = 'The program processes conversion data where the first column is elapsed time or sample temperature (in case of anisotherm nucleation) \n and the second column is the % of the crystalline fraction.')

            
            
            
            label_Tst = tk.Label(frame, text='Starting temperature [°C]', bg="#d9dcde", bd=10, font=("Times", 12), fg='black')
            label_Tst.place(relx=0.05, rely=0.37, relwidth=0.12, relheight=0.05, anchor='nw')
        
            Tst_entry = tk.Entry(frame, validate="key", bg='white', fg='black')
            Tst_entry.place(relx=0.18, rely=0.37, relwidth=0.12, relheight=0.05, anchor='nw')
            Tst_entry.insert(0, 135)
            if isTemp==True:
                Tst_entry.config(state='disabled')
        
            label_cr=tk.Label(frame, text='Cooling rate [°C/min]', bg="#d9dcde", bd=10, font=("Times", 16), fg='black')
            label_cr.place(relx=0.05, rely=0.295, relwidth=0.12, relheight=0.05)
        
            cr_entry = tk.Entry(frame, validate="key", bg='white', fg='black')
            cr_entry.place(relx=0.18, rely=0.295, relwidth=0.12, relheight=0.05, anchor='nw')
            cr_entry.insert(0, 10)
        
    else:
        
        if not label_conv.winfo_exists():
            label_Tst.destroy()
            Tst_entry.destroy()
            label_cr.destroy()
            cr_entry.destroy()
            mb4.destroy()
        
            label_conv=tk.Label(frame, text="Conversion curve", font=("Times", 16), bg="#d9dcde", fg='black')
            label_conv.place(relx=0.05, rely=0.22, relwidth=0.12, relheight=0.05)
        
            mb=tk.Menubutton(frame, text="Spherulite growth speed [m/s]",  bg="#a3b5bf", fg='black', font=("Times", 10))
            mb.place(relx=0.05, rely=0.37, relwidth=0.12, relheight=0.05)
            mb.menu =  tk.Menu(mb, tearoff = 0)
            mb["menu"] =  mb.menu
            mb.menu.add_command(label=str('Temperature [°C]'), font=("Times", 12), command=lambda: change_G(1) )
            mb.menu.add_command(label=str('Spherulite growth speed [m/s]'), font=("Times", 12), command=lambda: change_G(2) )
            G_entry = tk.Entry(frame, validate="key", bg='white', fg='black')
            G_entry.place(relx=0.18, rely=0.37, relwidth=0.12, relheight=0.05, anchor='nw')
            G_entry.insert(0, 0.00000003277)
        
        
def complex_simulation(parameter_szim, G_temp, HLparams, T_start, cr):
    global n_t, d_s, isAniso, isfirstrun, cc, rhoi_entry, F_entry
    global phantom_counter, text_area
    F_szim = parameter_szim[1]
    rhoi_szim = parameter_szim[0]
    
    
    if cc==1:
        matrix_size = 1 #ezt később lehet paraméterezni
        cost_matrix = np.zeros((matrix_size*2+1,matrix_size*2+1))
        bestloc=[matrix_size,matrix_size]
        for i in range(-matrix_size, matrix_size+1): #a tartomány felső határa már nincs benne
            for j in range(-matrix_size, matrix_size+1):
                if isAniso.get()==1:
                    cost_matrix[i+matrix_size,j+matrix_size]=simulate([rhoi_szim+i, F_szim+j], G_temp)
                else:
                    cost_matrix[i+matrix_size,j+matrix_size]=simulate_aniso([rhoi_szim+i, F_szim+j], HLparams, T_start, cr)
                isfirstrun = True
        bestloc[0] = np.argmin(cost_matrix)//(matrix_size*2+1)-matrix_size
        bestloc[1] = np.argmin(cost_matrix)%(matrix_size*2+1)-matrix_size
        #text_area.configure(state='normal')
      
    elif cc==2:
        cost_matrix = np.zeros(5)
        if isAniso.get()==1:
            cost_matrix[0]=simulate([rhoi_szim, F_szim], G_temp)
            isfirstrun = True
            cost_matrix[1]=simulate([-8.999, F_szim], G_temp)
            isfirstrun = True
            cost_matrix[2]=simulate([rhoi_szim, -11], G_temp)
            isfirstrun = True
            cost_matrix[3]=simulate([rhoi_szim+2, F_szim+2], G_temp)
            isfirstrun = True
            cost_matrix[4]=simulate([rhoi_szim-1, F_szim-1], G_temp)
            isfirstrun = True #valahol visszaállítja?
        else:
            cost_matrix[0]=simulate_aniso([rhoi_szim, F_szim], HLparams, T_start, cr)
            isfirstrun = True
            cost_matrix[1]=simulate_aniso([-8.999, F_szim], HLparams, T_start, cr)
            isfirstrun = True
            cost_matrix[2]=simulate_aniso([rhoi_szim, -11], HLparams, T_start, cr)
            isfirstrun = True
            cost_matrix[3]=simulate_aniso([rhoi_szim+2, F_szim+2], HLparams, T_start, cr)
            isfirstrun = True
            cost_matrix[4]=simulate_aniso([rhoi_szim-1, F_szim-1], HLparams, T_start, cr)
            isfirstrun = True
        bestloc = np.argmin(cost_matrix)
    
    #text_area.insert(tk.INSERT,' \n')
    #text_area.insert(tk.INSERT,cost_matrix)
    #text_area.insert(tk.INSERT,' \n')
    print(cost_matrix)
    print(bestloc) 
    n_t = n_t*10
    d_s = d_s/10
    if cc==1:
        #text_area.insert(tk.INSERT,[rhoi_szim+bestloc[0]+18, F_szim+bestloc[1]+18])
        #text_area.insert(tk.INSERT,' \n')
        print([rhoi_szim+bestloc[0]+18, F_szim+bestloc[1]+18])
        F_entry.delete(0,"end")
        rhoi_entry.delete(0,"end")
        rhoi_entry.insert(0, rhoi_szim+bestloc[0]+18)
        F_entry.insert(0, F_szim+bestloc[1]+18)
        if isAniso.get()==1:
            simulate([rhoi_szim+bestloc[0], F_szim+bestloc[1]], G_temp)
        else:
            simulate_aniso([rhoi_szim+bestloc[0], F_szim+bestloc[1]], HLparams, T_start, cr)
            
    elif cc==2:
        if isAniso.get()==1:
            if bestloc==0:
                simulate([rhoi_szim, F_szim], G_temp)
            elif bestloc==1:
                simulate([-8.999, F_szim], G_temp)
            elif bestloc==2:
                simulate([rhoi_szim, -11], G_temp)
            elif bestloc==3:
                simulate([rhoi_szim+2, F_szim+2], G_temp)
            elif bestloc==4:
                simulate([rhoi_szim-1, F_szim-1], G_temp)
        else:
            if bestloc==0:
                simulate_aniso([rhoi_szim, F_szim], HLparams, T_start, cr)
            elif bestloc==1:
                simulate_aniso([-8.999, F_szim], HLparams, T_start, cr)
            elif bestloc==2:
                simulate_aniso([rhoi_szim, -11], HLparams, T_start, cr)
            elif bestloc==3:
                simulate_aniso([rhoi_szim+2, F_szim+2], HLparams, T_start, cr)
            elif bestloc==4:
                simulate_aniso([rhoi_szim-1, F_szim-1], HLparams, T_start, cr)
    #text_area.configure(state='disabled')
    if phantom_counter>1000:
        DIYmessage("Warning!", "The list of simulated spherulites contains phantom spherulites.\nDue to this, the size and size distribution data may be inaccurate.")
    n_t = n_t/10
    d_s = d_s*10
    
    
def smooth():
    global conversion_input, tmaxC, tminC, isSmooth, filename, button_smooth
    
    k_0=4*math.pi*1e12*(3.28e-8)**2
    p_0=[3,k_0]

    if filename != 'no text':
        ci_temp = np.zeros((int(np.size(conversion_input)/2),2))
        ci_temp[:,0] = conversion_input[:,0]*60
        ci_temp[:,1] = conversion_input[:,1]/100


        t_tenyleges = ci_temp[:,0].tolist()
        k_tenyleges = ci_temp[:,1].tolist()


        k=0
        while (k<=(np.size(t_tenyleges)-1)): #a lista utolsó eleme a listaméret-1-edik :'(((((
            if (t_tenyleges[k]>tmaxC) or (t_tenyleges[k]<tminC):
                t_tenyleges.pop(k)
                k_tenyleges.pop(k)
            else:
                t_tenyleges[k]=t_tenyleges[k]-tminC
                k=k+1

        try:
            fitted = curve_fit(Avrami_func, t_tenyleges, k_tenyleges, p0=p_0)
            Avrami_params = fitted[0]


            k_illesztett = Avrami_func(t_tenyleges, Avrami_params[0], Avrami_params[1])

            j=0
            for k in range(int(np.size(ci_temp)/2)): 
                if (ci_temp[k,0]<tmaxC) and (ci_temp[k,0]>tminC):
                    ci_temp[k,0]=t_tenyleges[j]+tminC
                    ci_temp[k,1]=k_illesztett[j]
                    j=j+1
            
            
            conversion_input[:,0]=ci_temp[:,0]/60
            conversion_input[:,1]=ci_temp[:,1]*100
            isSmooth=True
            button_smooth.config(text="Done!")
        except (ValueError, RuntimeError, OverflowError) as e:
            DIYmessage("Error!", "Unable to smooth out the curve")
    else:
        DIYmessage("Error!", "No conversion curve found")
        
def sim_mode(mooode):
    global n_t, d_s, frame, custom_nt, nt_entry
    if custom_nt:
            nt_entry.destroy()
            root.update()
    if mooode==2:
        n_t = 1e6
        d_s = 0.2
        custom_nt=False
    elif mooode==1:
        n_t = 1e5
        d_s = 0.2 
        custom_nt=False
    elif mooode==3:
        nt_entry = tk.Entry(frame, validate="key", bg='white', fg='black')
        nt_entry.place(relx=0.18, rely=0.63, relwidth=0.12, relheight=0.05, anchor='nw')
        vcmd = (nt_entry.register(on_validate), '%P') #nem enged nem floatot bevenni
        nt_entry.config(validatecommand=vcmd)
        nt_entry.insert(0, "100000")
        custom_nt=True 
       # n_t=int(nt_entry.get())
        d_s = 0.2
        
def complex_mode(mooode):
    global isComplex, cc
    #if isComplex:
        #complex_entry.destroy()
        #root.update()
    if mooode==1:
        isComplex = False
    elif mooode==2:
        isComplex = True
        cc = 1
    elif mooode==3:
        isComplex = True
        cc = 2
        #complex_entry = tk.Entry(frame, validate="key")
        #complex_entry.place(relx=0.18, rely=0.675, relwidth=0.12, relheight=0.05, anchor='nw')
        #vcmd = (nt_entry.register(on_validate), '%P')
        #complex_entry.config(validatecommand=vcmd)
        #complex_entry.insert(0, "2")
        
def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()

def set_isTemp(a):
    global Tst_entry, isTemp, isAniso
    if a==1:
        isTemp=False
        if isAniso.get()==2:
            Tst_entry.destroy()
            Tst_entry = tk.Entry(frame, validate="key", bg='white', fg='black')
            Tst_entry.place(relx=0.18, rely=0.37, relwidth=0.12, relheight=0.05, anchor='nw')
            Tst_entry.insert(0, 135)
        #print(isTemp)
    else:
        isTemp=True
        if isAniso.get()==2:
            Tst_entry.config(state='disabled')

    
    

n_t = 1e5
d_s = 0.2
G_temp = 3.277e-8
cr = 5
T_start = 408
Tcp_default = 398.15
custom_nt=False
isfirstrun=False
isComplex=False
isSmooth=False
G_setting_mode=True
filename = 'no text' #DSC ellenőrzésére
avg_sph=0
phantom_counter=0 #kezdőértékek, hogy ne adjon ki hibákat ha enélkül elindítom

hazeWindowfirstopen = False
dataFitfirstopen = False
firstError = False 

p_0 = (20277.70942,-2.67238)
fitted = curve_fit(HL_func, [398.15, 403.15, 408.15, 413.15, 418.15], [2.77e-7, 9.86e-8, 3.28e-8, 1.21e-8, 4.77e-9], p0=p_0)
HL_params = fitted[0]

U=755 
Tg=263.15 
Tref=30.0 
T0m=481.15

HL_params = [HL_params[0], HL_params[1], U, Tg, Tref, T0m]



goo = True
Tmax = 423.15
while goo:
    Tmax = Tmax + 0.05
    if pcrystal2_mod.fcrystal.hl_calc(Tmax, HL_params)<1e-10:
        goo = False
goo = True
Tmin = 313.15
while goo:
    Tmin = Tmin - 0.05
    if pcrystal2_mod.fcrystal.hl_calc(Tmin, HL_params)<1e-10:
        goo = False


        
        
root = tk.Tk()
root.title('Morphological structure estimation of semicrystalline polymers')

canvas=tk.Canvas(root, height=HEIGHT, width=WIDTH)
canvas.pack()

frame=tk.Frame(root, bg="#e75480", bd=5)
frame.place(relx=0.5, rely=0, relwidth=0.95, relheight=0.95, anchor="n")

image_path = resource_path("_03_Tc135_tc30_10X.gif")
"""background_image=tk.PhotoImage(file = "/home/aliz/Crystal/fcrystal/03_Tc135_tc30_10X.png")"""
background_image=tk.PhotoImage(file = image_path)
background_label = tk.Label(frame, image=background_image)
background_label.place(x=0, y=0, relwidth=1, relheight=1)



frame1=tk.Frame(frame, bg="#ffffff", bd=5)
frame1.place(relx=0.32, rely=0.15, relwidth=0.33, relheight=0.55, anchor="nw")
    
frame2=tk.Frame(frame, bg="#ffffff", bd=5)
frame2.place(relx=0.67, rely=0.15, relwidth=0.325, relheight=0.55, anchor="nw")

isSeparate = tk.BooleanVar()
c1 = tk.Checkbutton(frame, text='Plot instantaneous and continuous seeds separately', variable=isSeparate, font=("Times", 12), bg="#a3b5bf", fg='black', command = lambda:plot_sph(filename))
c1.place(relx=0.67, rely=0.705, relwidth=0.325, relheight=0.03, anchor="nw")


label_results=tk.Label(frame, text="The calculated values:", bg="#d9dcde", bd=10, font=("Times", 16), fg='black')
label_results.place(relx=0.325, rely=0.8, relwidth=0.425, relheight=0.05)

frame3=tk.Frame(frame, bg="#ffffff", bd=5)
frame3.place(relx=0.325, rely=0.85, relwidth=0.425, relheight=0.135, anchor="nw")

text_result=tk.Text(frame3, height=6, width=120, bg='white', fg='black')
text_result.pack()

label_avgsph=tk.Label(frame, text='Average spherulite size ['+'um]:', bg="#d9dcde", bd=10, font=("Times", 16), fg='black')
label_avgsph.place(relx=0.8, rely=0.8, relwidth=0.18, relheight=0.05)

label_avgsph2=tk.Label(frame, text=str(avg_sph), bg="#d9dcde", bd=10, font=("Times", 16), fg='black')
label_avgsph2.place(relx=0.8, rely=0.85, relwidth=0.18, relheight=0.07)


label_F=tk.Label(frame, text="log10(F [1/m^3s])", bg="#d9dcde", bd=10, font=("Times", 16), fg='black')
label_F.place(relx=0.05, rely=0.445, relwidth=0.12, relheight=0.05)

F_entry = tk.Entry(frame, validate="key", bg='white', fg='black')
F_entry.place(relx=0.18, rely=0.445, relwidth=0.12, relheight=0.05, anchor='nw')
vcmd = (F_entry.register(on_validate), '%P') #nem enged nem floatot bevenni
F_entry.config(validatecommand=vcmd)
F_entry.insert(0, "12")

CreateToolTip(label_F, text = 'F is the birthspeed of continuously born seeds [1/m^3s]')

label_rhoi=tk.Label(frame, text="log10(rhoi [1/m^3])", bg="#d9dcde", bd=10, font=("Times", 16), fg='black')
label_rhoi.place(relx=0.05, rely=0.52, relwidth=0.12, relheight=0.05)

rhoi_entry = tk.Entry(frame, validate="key", bg='white', fg='black')
rhoi_entry.place(relx=0.18, rely=0.52, relwidth=0.12, relheight=0.05, anchor='nw')
vcmd = (rhoi_entry.register(on_validate), '%P')
rhoi_entry.config(validatecommand=vcmd)
rhoi_entry.insert(0, "12")

CreateToolTip(label_rhoi, text = 'rhoi is the density of instantaneously born seeds [1/m^3]')


label=tk.Label(frame, text="Morphological structure estimation of semicrystalline polymers", bg="#d9dcde", bd=10, font=("Times", 24), fg='black', anchor='n')
label.place(relx=0.15, rely=0, relwidth=0.75, relheight=0.075)


isTemp=False
"""mb4=tk.Menubutton(frame, text="Conversion curve", font=("Times", 16), bg="#a3b5bf", fg='black')
mb4.place(relx=0.05, rely=0.22, relwidth=0.12, relheight=0.05)
mb4.menu =  tk.Menu(mb4, tearoff = 0)
mb4["menu"] =  mb4.menu
mb4.menu.add_command(label=str('time-%'),command = lambda:set_isTemp(1) )
mb4.menu.add_command(label=str('temperature-% (set cooling rate first!)'),command = lambda:set_isTemp(2) )"""

label_conv=tk.Label(frame, text="Conversion curve", font=("Times", 16), bg="#d9dcde", fg='black')
label_conv.place(relx=0.05, rely=0.22, relwidth=0.12, relheight=0.05)
#CreateToolTip(label_conv, text = 'The program processes conversion data where the first column is elapsed time or sample temperature (in case of anisotherm nucleation) \n and the second column is the % of the crystalline fraction.')


button = tk.Button(frame, text="Browse", bg="#778899", font=("Times", 16), command=lambda: browse(), fg='black')
button.place(relx=0.18, rely=0.22, relwidth=0.12, relheight=0.05)

buttons = tk.Button(frame, text="Simulate", bg="#778899", font=("Times", 16), command=lambda: simulate_start(), fg='black')
buttons.place(relx=0.1, rely=0.935, relwidth=0.15, relheight=0.05)

button_haze = tk.Button(frame, text="Simulate haze", bg="#778899", font=("Times", 16), command=lambda: HazeSim(avg_sph), fg='black')
button_haze.place(relx=0.8, rely=0.935, relwidth=0.18, relheight=0.05)


button_datafit = tk.Button(frame, text="Fit growth speed", bg="#778899", font=("Times", 16), command=lambda: dataFit_window(), fg='black')
button_datafit.place(relx=0.18, rely=0.14, relwidth=0.12, relheight=0.05)
CreateToolTip(button_datafit, text = 'In this menu, you can enter the measurement data regarding the growth speed of spherulites in the polymer you are using. \n The program will perform a curve fitting to determine the neccessary parameters.')

"""button_smooth = tk.Button(frame, text="Smooth out curve", bg="#778899", font=("Times", 12), command=lambda: smooth(), fg='black')
button_smooth.place(relx=0.32, rely=0.705, relwidth=0.12, relheight=0.03)"""

mb=tk.Menubutton(frame, text="Spherulite growth speed [m/s]", bg="#a3b5bf", fg='black', font=("Times", 10))
mb.place(relx=0.05, rely=0.37, relwidth=0.12, relheight=0.05)
mb.menu =  tk.Menu(mb, tearoff = 0)
mb["menu"] =  mb.menu
mb.menu.add_command(label=str('Temperature [°C]'), font=("Times", 12), command=lambda: change_G(1) )
mb.menu.add_command(label=str('Spherulite growth speed [m/s]'), font=("Times", 12), command=lambda: change_G(2) )


"""label_G=tk.Label(frame, text="G [m/s]", bg="#a3b5bf", bd=10, font=("Times", 16))
label_G.place(relx=0.05, rely=0.425, relwidth=0.12, relheight=0.05)

CreateToolTip(label_G, text = 'G is the growth speed of the spherulites [m/s]')"""

G_entry = tk.Entry(frame, validate="key", bg='white', fg='black')
G_entry.place(relx=0.18, rely=0.37, relwidth=0.12, relheight=0.05, anchor='nw')
G_entry.insert(0, 0.00000003277)
    

"""mb2=tk.Menubutton(frame, text="Sampling mode", bg="#a3b5bf", fg='black')
mb2.place(relx=0.05, rely=0.6, relwidth=0.12, relheight=0.05)
mb2.menu =  tk.Menu(mb2, tearoff = 0)
mb2["menu"] =  mb2.menu
mb2.menu.add_command(label=str('Faster (10^5 test points)'),command=lambda: sim_mode(1) )
mb2.menu.add_command(label=str('Slower (10^6 test points)'),command=lambda: sim_mode(2) )
mb2.menu.add_command(label=str('Custom number of test points'),command=lambda: sim_mode(3) )"""
nt_label = tk.Label(frame, text="Sampling mode", font=("Times", 12), bg="#d9dcde", fg='black')
nt_label.place(relx=0.05, rely=0.6, relwidth=0.12, relheight=0.025)

nt_radio=tk.IntVar() 
nt_r1 = tk.Radiobutton(frame, text='Faster (10^5 test points)', variable=nt_radio, value=1, bg="#a3b5bf", font=("Times", 10), command=lambda: sim_mode(1), fg='black')
nt_r1.place(relx=0.05, rely=0.63, relwidth=0.12, relheight=0.025)
nt_r2 = tk.Radiobutton(frame, text='Slower (10^6 test points)', variable=nt_radio, value=2, bg="#a3b5bf", font=("Times", 10), command=lambda: sim_mode(2), fg='black')
nt_r2.place(relx=0.05, rely=0.66, relwidth=0.12, relheight=0.025)
nt_r3 = tk.Radiobutton(frame, text='Custom number of test points', variable=nt_radio, value=3, bg="#a3b5bf", font=("Times", 10),  command=lambda: sim_mode(3), fg='black')
nt_r3.place(relx=0.05, rely=0.69, relwidth=0.12, relheight=0.025)
nt_r1.select()

CreateToolTip(nt_label, text = 'A higher number of test points results in a more accurate simulation. \nIn case of mainly non-instantaneous nucleation, 10^6 or higher number of testpoints is recommended for accurate results.')

"""mb3=tk.Menubutton(frame, text="Simulation complexity",  bg="#a3b5bf", fg='black')
mb3.place(relx=0.05, rely=0.675, relwidth=0.12, relheight=0.05)
mb3.menu =  tk.Menu(mb3, tearoff = 0)
mb3["menu"] =  mb3.menu
mb3.menu.add_command(label=str('Simple'),command=lambda: complex_mode(1) )
mb3.menu.add_command(label=str('Complex'),command=lambda: complex_mode(2) )
mb3.menu.add_command(label=str('Complex 2.0'),command=lambda: complex_mode(3) )
CreateToolTip(mb3, text = 'Complex simulations test several starting values around the given ones \n and then choose the best fitting to perform a more sensitive simulation with.')"""


complex_label = tk.Label(frame, text="Simulation mode", font=("Times", 12), bg="#d9dcde", fg='black')
complex_label.place(relx=0.05, rely=0.74, relwidth=0.12, relheight=0.025)

complex_radio=tk.IntVar() 
complex_r1 = tk.Radiobutton(frame, text='Simple', variable=complex_radio, value=1, bg="#a3b5bf", font=("Times", 12), command=lambda: complex_mode(1), fg='black')
complex_r1.place(relx=0.05, rely=0.77, relwidth=0.075, relheight=0.025)
complex_r2 = tk.Radiobutton(frame, text='Advanced', variable=complex_radio, value=2, bg="#a3b5bf", font=("Times", 12), command=lambda: complex_mode(2), fg='black')
complex_r2.place(relx=0.13, rely=0.77, relwidth=0.075, relheight=0.025)
"""complex_r3 = tk.Radiobutton(frame, text='Advanced 2.0', variable=complex_radio, value=3, bg="#a3b5bf", font=("Times", 12), command=lambda: complex_mode(3), fg='black')
complex_r3.place(relx=0.21, rely=0.77, relwidth=0.075, relheight=0.025)"""
complex_r1.select()

CreateToolTip(complex_label, text = 'Advanced simulations test several starting values around the given ones \n and then choose the best fitting to perform a more sensitive simulation with.')

isAniso=tk.IntVar() 
r1 = tk.Radiobutton(frame, text='isotherm', variable=isAniso, value=1, bg="#a3b5bf", command=lambda: resetArea(1), font=("Times", 12), fg='black')
r1.place(relx=0.05, rely=0.14, relwidth=0.075, relheight=0.025)
r2 = tk.Radiobutton(frame, text='anisotherm', variable=isAniso, value=2, bg="#a3b5bf", command=lambda: resetArea(2), font=("Times", 12), fg='black')
r2.place(relx=0.05, rely=0.17, relwidth=0.075, relheight=0.025)
r1.select()


label_cr=tk.Label(frame, text='Cooling rate [°C/min]', bg="#d9dcde", bd=10, font=("Times", 16), fg='black')
label_cr.place(relx=0.05, rely=0.295, relwidth=0.12, relheight=0.05)
label_cr.destroy() #ez azért kell, hogy a label_cr-re legyen reference de az exists 0-t adjon rá defaultban



root.mainloop()


# ### 

# In[ ]:





# In[ ]:





# In[ ]:




