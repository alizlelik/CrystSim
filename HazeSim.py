#!/usr/bin/env python
# coding: utf-8

# # Init

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from math import pi
import scipy as sp

import mie_aniso


# In[5]:


#%matplotlib notebook


# In[6]:


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


# Single parameter

# In[7]:


import tkinter as tk
from tkinter import font

HEIGHT = 5000
WIDTH = 7500

def test_function(entry,entry2,entry3,entry4,entry5): 
    
    r0 = float(entry2) # radius in micron
    dn = float(entry5) # Birefringence
    n_bg = float(entry4) # background refractive index
    LL = float(entry) # sample thickness in micron
    wl = float(entry3) # wavelength in micron
    
    label6["text"]= "Sample thickness [um]: %s \nSpherulite radius [um]:  %s \nWavelength [um]: %s \nRefractive index [-]:  %s \nBirefringence [-]: %s" % (LL, r0, wl, n_bg, dn)
    
    print("Sample thickness:", LL)
    print("Spherulite radius:", r0)
    print("Wavelength:", wl)
    print("Refractive index:", n_bg)
    print("Birefringence:", dn)
    
    #r0 = 9.28 # radius in micron
    #wl = 0.6328 # wavelength in micron
    #LL = 1000 # sample thickness in micron
    
    #n_bg = 1.5 # background refractive index
    #dn = 0.002 # Birefringence
    
    n_r = n_bg+dn*0.5 # radial refractive index
    n_t = n_bg-dn*0.5 # tangential refractive index
    
    f_V = 0.4 # volume fraction (sum v_i)/V_container
    n0 = f_V/(4./3*r0**3*pi) # nucleus density N/V [1/um^3]
    
    # create single mie object
    sm = mie_aniso.single_mie(r0, wl, n_r**2, n_t**2, n_bg**2)
    
    # calculate forward and backward cross section
    (C_f, C_b) = sm.csca_fb()
    
    # total transmission
    label7["text"]="Total transmission: %s \nDiffuse transmission: %s \nSpecular transmission: %s \nHaze: %s" % (I_trans_tot(C_b*n0, LL),I_trans_dif(C_b*n0, C_f*n0, LL),I_trans_spc(C_b*n0, C_f*n0, LL), haze(C_b*n0, C_f*n0, LL))
    print("total transmission:", I_trans_tot(C_b*n0, LL))
    print("diffuse transmission:", I_trans_dif(C_b*n0, C_f*n0, LL))
    print("specular transmission:", I_trans_spc(C_b*n0, C_f*n0, LL))
    print("haze:", haze(C_b*n0, C_f*n0, LL))
    
root = tk.Tk()

canvas=tk.Canvas(root, height=HEIGHT, width=WIDTH)
canvas.pack()

#background_image = tk.PhotoImage(file="hatter.png")
#background_label = tk.Label(root, image=background_image)
#background_label.place(relwidth=1, relheight=1)

frame=tk.Frame(root, bg="gray", bd=5)
frame.place(relx=0.5, rely=0, relwidth=1, relheight=1, anchor="n")

label=tk.Label(frame, text="Haze simulation for semicrystalline polymers", bg="#c4cdcf", bd=10, font=("Times", 36))
label.place(relx=0.25, rely=0, relwidth=0.5, relheight=0.15)

label1=tk.Label(frame, text="Sample thickness", bg="#c4cdcf", font=("Times", 24))
label1.place(relx=0.1, rely=0.2, relwidth=0.1, relheight=0.1)
            
entry=tk.Entry(frame, font=("Times", 24), justify="center", bd = 5)
entry.place(relx=0.2, rely=0.2, relwidth=0.1, relheight=0.1)

label2=tk.Label(frame, text="Spherulite radius", bg="#c4cdcf", font=("Times", 24))
label2.place(relx=0.1, rely=0.31, relwidth=0.1, relheight=0.1)

entry2=tk.Entry(frame, font=("Times", 24), justify="center", bd = 5)
entry2.place(relx=0.2, rely=0.31, relwidth=0.1, relheight=0.1)

label3=tk.Label(frame, text="Wavelength", bg="#c4cdcf", font=("Times", 24))
label3.place(relx=0.1, rely=0.42, relwidth=0.1, relheight=0.1)
            
entry3=tk.Entry(frame, font=("Times", 24), justify="center", bd = 5)
entry3.place(relx=0.2, rely=0.42, relwidth=0.1, relheight=0.1)
            
label4=tk.Label(frame, text="Refractive index", bg="#c4cdcf",font=("Times", 24), bd = 5)
label4.place(relx=0.1, rely=0.53, relwidth=0.1, relheight=0.1)
            
entry4=tk.Entry(frame, font=("Times", 24), justify="center", bd = 5)
entry4.place(relx=0.2, rely=0.53, relwidth=0.1, relheight=0.1)

label5=tk.Label(frame, text="Birefringence", bg="#c4cdcf",font=("Times", 24), bd = 5)
label5.place(relx=0.1, rely=0.64, relwidth=0.1, relheight=0.1)
            
entry5=tk.Entry(frame, font=("Times", 24), justify="center", bd = 5)
entry5.place(relx=0.2, rely=0.64, relwidth=0.1, relheight=0.1)

button = tk.Button(frame, text="Simulate", bg="#7be9ed", font=("Times", 24), command=lambda: test_function(entry.get(), entry2.get(), entry3.get(), entry4.get(), entry5.get()))
button.place(relx=0.1, rely=0.8, relwidth=0.1, relheight=0.1)

#label5=tk.Label(frame, text="Sample width", bg="#c4cdcf",font=("Times", 20))
#label5.place(relx=0.1, rely=0.53, relwidth=0.1, relheight=0.1)

side_frame=tk.Frame(root, bg="#96d9e3")
side_frame.place(relx=0.61, rely=0.2, relwidth=0.5, relheight=0.7, anchor="n")

label6=tk.Label(side_frame, font=("Times", 40), anchor="nw", justify="left", bd=4)
label6.place(relwidth=1, relheight=0.5)

label7=tk.Label(side_frame, font=("Times", 40), anchor="nw",bg="#7be9ed", justify="left", bd=4)
label7.place(relx=0, rely=0.5, relwidth=1, relheight=0.5)

root.mainloop()


# In[ ]:





# In[ ]:




