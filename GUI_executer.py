# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer Graphical User Interface function

modififications to previosly versions must be specified here using WAVE codes:
    
W-23.01-1.00.3-00
A-01.01-1.00.3-00
V-22.01-1.00.3-00
E-05.01-1.00.3-00
"""

from tkinter import *
from tkinter.filedialog import *
from tkinter.font import *

import webbrowser
import tkinter.ttk as tka
import subprocess

from queue import Queue, Empty
from tkinter import messagebox
from tkinter import filedialog
from tkinter import simpledialog
from tkinter import ttk

def open_button():
    global path_open, path_temp, path_level, path_meteo, path_senso, path_fetch, folder_path
    path_open = askopenfilename(defaultextension='.set', filetypes=[('SET files','*.set')])
    
    
    if path_open:
        with open(str(path_open),'r') as reader:
            path_temp = reader.readline()
            path_temp = path_temp.replace('\n','')
            path_meteo = reader.readline()
            path_meteo = path_meteo.replace('\n','')
            path_senso = reader.readline()
            path_senso = path_senso.replace('\n','')
            height_wind.delete(0,'end')
            height_wind.insert(END,float(reader.readline()))
            radiation_on.set(int(reader.readline())) 
            reference_level.delete(0,'end')
            reference_level.insert(END,float(reader.readline()))
            contri_wind_var.set(float(reader.readline()))
            latitude_var.set(float(reader.readline()))
            dpi.set(int(reader.readline()))
            typechoose.set(int(reader.readline()))
            if int(typechoose.get()) == 1:
                vari_len.delete(0,'end')
                vari_len.insert(END,float(reader.readline())) 
            elif int(typechoose.get()) == 2:
                path_fetch = reader.readline()
                path_fetch = path_fetch.replace('\n','')
            typelevel.set(int(reader.readline()))
            if int(typelevel.get()) == 1:
                unif_level.delete(0,'end')
                unif_level.insert(END,float(reader.readline())) 
            elif int(typelevel.get()) == 2:
                path_level = reader.readline()
                path_level = path_level.replace('\n','')

            meta.delete(0,'end')
            meta.insert(END,float(reader.readline())) 
            
            filterchoose.set(int(reader.readline()))
            if int(filterchoose.get()) == 2:
                
                lv1.config(state='normal')
                hv1.config(state='normal')                
                
                lv1.delete(0,'end')
                hv1.delete(0,'end')

                lv1.insert(END,float(reader.readline()))
                hv1.insert(END,float(reader.readline()))
            else:
                for i in range(2):
                    nextline = next(reader)


            winsizechoose.set(int(reader.readline()))
            if int(winsizechoose.get()) == 2:  
                winsize.config(state='normal')                          
                winsize.delete(0,'end')
                winsize.insert(END,float(reader.readline()))
            else:
                nextline = next(reader)

               
            win = reader.readline()
            for i in range(4):
                if win.strip() == welch_windows[i]:
                    windows_var.set(welch_windows[i])
                
            mon = reader.readline()
            for i in range(3):
                if mon.strip() == wavelet_windows[i]:
                    mother_var.set(wavelet_windows[i])

            isotherm_button.set(int(reader.readline()))

            num = reader.readline().strip()
            isoa.delete(0,'end')
            if num == '-999':
                isotherma.set(0)
                isoa.config(state='disable')
            else:
                isotherma.set(1)
                isoa.config(state='normal')
                isoa.insert(END,float(num))

            num = reader.readline().strip()
            isob.delete(0,'end')
            if num == '-999':
                isothermb.set(0)
                isob.config(state='disable')
            else:
                isothermb.set(1)
                isob.config(state='normal')
                isob.insert(END,float(num))

            num = reader.readline().strip()
            isoc.delete(0,'end')
            if num == '-999':
                isothermc.set(0)
                isoc.config(state='disable')
            else:
                isothermc.set(1)
                isoc.config(state='normal')
                isoc.insert(END,float(num))
                
            num = reader.readline().strip()
            isod.delete(0,'end')
            if num == '-999':
                isothermd.set(0)
                isod.config(state='disable')
            else:
                isothermd.set(1)
                isod.config(state='normal')                
                isod.insert(END,float(num))
                
            c1 = reader.readline()
            c2 = reader.readline()
            for i in range(6):
                if c1.strip() == comparison_options[i]:
                    compa1_var.set(comparison_options[i])            
                if c2.strip() == comparison_options[i]:
                    compa2_var.set(comparison_options[i])             
 

            sensor_button.set(int(reader.readline()))
            num = reader.readline().strip()
            sena.delete(0,'end')
            if num == '-999':
                sensora.set(0)
                sena.config(state='disable')
            else:
                sensora.set(1)
                sena.config(state='normal')
                sena.insert(END,int(num))

            num = reader.readline().strip()
            senb.delete(0,'end')
            if num == '-999':
                sensorb.set(0)
                senb.config(state='disable')
            else:
                sensorb.set(1)
                senb.config(state='normal')
                senb.insert(END,int(num))

            num = reader.readline().strip()
            senc.delete(0,'end')
            if num == '-999':
                sensorc.set(0)
                senc.config(state='disable')
            else:
                sensorc.set(1)
                senc.config(state='normal')
                senc.insert(END,int(num))
                
            num = reader.readline().strip()
            send.delete(0,'end')
            if num == '-999':
                sensord.set(0)
                send.config(state='disable')
            else:
                sensord.set(1)
                send.config(state='normal')                
                send.insert(END,int(num))    
                
            smooth_button.set(int(reader.readline()))  
            folder_path = reader.readline()
            
def AboutCallBack():
   msg = messagebox.showinfo( "About", " Interwave Analyzer - Version 1.00.3 \n Copyright (C) 2019 Rafael de Carvalho Bueno \n All rights reserved \n \n Developed by \n Rafael de Carvalho Bueno \n\n Co-developer \n Tobias Bleninger \n\n improvements and betterments by \n Andreas Lorke \n\n Report problems and improvements to email adresss below \n rafael.bueno.itt@gmail.com\n \n for mor information, see: \n www.bit.ly/interwave_analyzer \n ")    

def OpenUrl(url):
    webbrowser.open_new(url)


def save_settings(temporary):


    if temporary == 'temporary':
        data_file = 'temporary.txt' 
        space     = '\n'
    else:   
        f = filedialog.asksaveasfile(defaultextension=".set")
        data_file = str(f.name)
        space     =  '\n'

    with open(data_file, 'w') as data:

        data.write(str(path_temp)+space)
        data.write(str(path_meteo)+space)
        data.write(str(path_senso)+space)
        data.write(str(height_wind.get())+'\n')
        data.write(str(radiation_on.get())+'\n')
        data.write(str(reference_level.get())+'\n')
        data.write(str(contri_wind_var.get())+'\n')
        data.write(str(latitude_var.get())+'\n')
        data.write(str(dpi.get())+'\n')
        data.write(str(int(typechoose.get()))+'\n')
        if int(typechoose.get()) == 1:
            data.write(str(vari_len.get())+'\n')
        elif int(typechoose.get()) == 2:
            data.write(str(path_fetch)+space)
        
        data.write(str(int(typelevel.get()))+'\n')
        if int(typelevel.get()) == 1:
            data.write(str(unif_level.get())+'\n')
        elif int(typelevel.get()) == 2:
            data.write(str(path_level)+space)
        data.write(str(meta.get())+'\n')
        
        
        data.write(str(int(filterchoose.get()))+'\n')
        if int(filterchoose.get()) == 2:
            data.write(str(lv1.get())+'\n')
            data.write(str(hv1.get())+'\n')
        else:
            data.write(str(-999)+'\n')
            data.write(str(-999)+'\n')
            
        data.write(str(int(winsizechoose.get()))+'\n')
        if int(winsizechoose.get()) == 2:
            data.write(str(winsize.get())+'\n')
        else:
            data.write(str(-999)+'\n')

        data.write(str(windows_var.get())+'\n')
        data.write(str(mother_var.get())+'\n')
        data.write(str(isotherm_button.get())+'\n')
        if int(isotherma.get()) == 1:
            data.write(str(isoa.get())+'\n')
        else:
            data.write(str(-999)+'\n')
            
        if int(isothermb.get()) == 1:
            data.write(str(isob.get())+'\n')
        else:
            data.write(str(-999)+'\n')
            
        if int(isothermc.get()) == 1:
            data.write(str(isoc.get())+'\n')
        else:
            data.write(str(-999)+'\n')
            
        if int(isothermd.get()) == 1:
            data.write(str(isod.get())+'\n')
        else:
            data.write(str(-999)+'\n')        
        data.write(str(compa1_var.get())+'\n')
        data.write(str(compa2_var.get())+'\n')
        data.write(str(sensor_button.get())+'\n')
        if int(sensora.get()) == 1:
            data.write(str(sena.get())+'\n')
        else:
            data.write(str(-999)+'\n')
            
        if int(sensorb.get()) == 1:
            data.write(str(senb.get())+'\n')
        else:
            data.write(str(-999)+'\n')
            
        if int(sensorc.get()) == 1:
            data.write(str(senc.get())+'\n')
        else:
            data.write(str(-999)+'\n')
            
        if int(sensord.get()) == 1:
            data.write(str(send.get())+'\n')
        else:
            data.write(str(-999)+'\n')        
        data.write(str(smooth_button.get())+'\n')
        data.write(str(folder_path))
        


class StdRedirector():
    def __init__(self, text_widget):
        self.text_space = text_widget

    def write(self, string):
        self.text_space.config(state=NORMAL)
        self.text_space.insert("end", string)
        self.text_space.see("end")
        self.text_space.config(state=DISABLED)
       
def export_data():
    
    import bueno_iwm as iwm
    save_settings('temporary')  
    iwm.main()
    
def temperature_function():
    global path_temp
    path_temp = askopenfilename(defaultextension='.tem', filetypes=[('TEM files','*.tem')])

             
def level_function():
    global path_level
    path_level = askopenfilename(defaultextension='.niv', filetypes=[('NIV files','*.niv')])

   
def meteo_function():
    global path_meteo
    path_meteo = askopenfilename(defaultextension='.met', filetypes=[('MET files','*.met')])
    
     
def senso_function():
    global path_senso
    path_senso = askopenfilename(defaultextension='.sen', filetypes=[('SEN files','*.sen')])
    
 
def fetch_function():
    global path_fetch
    path_fetch = askopenfilename(defaultextension='.fet', filetypes=[('FET files','*.fet')])
    


def output_folder():


    global folder_path
    folder_path = filedialog.askdirectory()

def selected_len():
    global fetch_type
    fetch_type = int(typechoose.get())

    if int(typechoose.get())==2:
        file_len.config(state='normal')
        vari_len.config(state='disable')
    elif int(typechoose.get())==1: 
        file_len.config(state='disable')
        vari_len.config(state='normal')


def selected_level():
    global level_type, unif_level
    level_type = int(typelevel.get())

    if int(typelevel.get())==2:
        file_level.config(state='normal')
        unif_level.config(state='disable')
    elif int(typechoose.get())==1: 
        file_level.config(state='disable')
        unif_level.config(state='normal')
        
def selected_filter():
    global filter_type
    filter_type = int(filterchoose.get())

    if int(filterchoose.get())==1:
        lv1.config(state='disable')
        hv1.config(state='disable')

    elif int(filterchoose.get())==2: 
        lv1.config(state='normal')
        hv1.config(state='normal')

def selected_windowsize():
    global windsize_type
    windsize_type = int(winsizechoose.get())

    if int(winsizechoose.get())==1:
        winsize.config(state='disable')

    elif int(winsizechoose.get())==2: 
        winsize.config(state='normal')

def selected_isotherms():
    global isotherms_type
    isotherms_type = int(isotherm_button.get())

    if int(isotherm_button.get())==0:
        iso1.config(state='disable')
        iso2.config(state='disable')
        iso3.config(state='disable')
        iso4.config(state='disable')
    elif int(isotherm_button.get())==1: 
        iso1.config(state='normal')
        iso2.config(state='normal')
        iso3.config(state='normal')
        iso4.config(state='normal')

def selected_isoa():
    global isoa_type
    isoa_type = int(isotherma.get())
    
    if int(isotherma.get()) == 0:
        isoa.config(state='disable')
    elif int(isotherma.get()) ==1:
        isoa.config(state='normal')
        
def selected_isob():
    global isob_type
    isob_type = int(isothermb.get())
    
    if int(isothermb.get()) == 0:
        isob.config(state='disable')
    elif int(isothermb.get()) ==1:
        isob.config(state='normal')
        
def selected_isoc():
    global isoc_type
    isoc_type = int(isothermc.get())
    
    if int(isothermc.get()) == 0:
        isoc.config(state='disable')
    elif int(isothermc.get()) ==1:
        isoc.config(state='normal')
        
def selected_isod():
    global isod_type
    isod_type = int(isothermd.get())
    
    if int(isothermd.get()) == 0:
        isod.config(state='disable')
    elif int(isothermd.get()) ==1:
        isod.config(state='normal')
       

def selected_sensor():
    global sensor_type
    sensor_type = int(sensor_button.get())

    if int(sensor_button.get())==0:
        sen1.config(state='disable')
        sen2.config(state='disable')
        sen3.config(state='disable')
        sen4.config(state='disable')
    elif int(sensor_button.get())==1: 
        sen1.config(state='normal')
        sen2.config(state='normal')
        sen3.config(state='normal')
        sen4.config(state='normal')

def selected_sena():
    global sena_type
    sena_type = int(sensora.get())
    
    if int(sensora.get()) == 0:
        sena.config(state='disable')
    elif int(sensora.get()) ==1:
        sena.config(state='normal')
        
def selected_senb():
    global senb_type
    senb_type = int(sensorb.get())
    
    if int(sensorb.get()) == 0:
        senb.config(state='disable')
    elif int(sensorb.get()) ==1:
        senb.config(state='normal')
        
def selected_senc():
    global senc_type
    senc_type = int(sensorc.get())
    
    if int(sensorc.get()) == 0:
        senc.config(state='disable')
    elif int(sensorc.get()) ==1:
        senc.config(state='normal')
        
        
def selected_send():
    global send_type
    send_type = int(sensord.get())
    
    if int(sensord.get()) == 0:
        send.config(state='disable')
    elif int(sensord.get()) ==1:
        send.config(state='normal')
   
def window_destroy():
    answer = messagebox.askyesno("exit",'do you really want to exit ?')
    if answer:
        window.destroy()
     
        
window = Tk()


window.geometry("800x800")
window.iconbitmap("interwave_icon.ico")
window.title("Interwave Analyzer") 


menubar = Menu(window)
filemenu = Menu(menubar, tearoff = 0)
filemenu.add_command(label = "Open", command = open_button)
filemenu.add_command(label = "Save as...", command = lambda: save_settings(0))

filemenu.add_separator()

filemenu.add_command(label = "Exit", command = window_destroy)
menubar.add_cascade(label = "File", menu = filemenu)
editmenu = Menu(menubar, tearoff=0)


editmenu.add_separator()

helpmenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label = "Help", menu = helpmenu)
url = 'https://sites.google.com/view/interwaveanalyzer/interwave-analyzer'
helpmenu.add_command(label = "Manual", command = lambda aurl=url:OpenUrl(aurl))
helpmenu.add_command(label = "About", command = AboutCallBack)

tabControl = ttk.Notebook(window)
tabControl.grid(row=1, column=0, columnspan=50, rowspan=15,sticky='NESW')
tab1 = Frame(tabControl)
tabControl.add(tab1,text='Input data')

tab2 = Frame(tabControl)
tabControl.add(tab2,text='Spectral analysis and definitions')

tab3 = Frame(tabControl)
tabControl.add(tab3,text='Isotherm Analysis')


tab4 = Frame(tabControl)
tabControl.add(tab4,text='Output and run')

Label(tab1,anchor="w",font="Verdana 8 bold", text="Input files").grid(row=1,column=0,pady=8,sticky='w')
# ---------------------Temperature time serie----------------------------------

Label(tab1,anchor="w", text="Temperature data:").grid(row=2,column=0,pady=4,sticky='w')
Button(tab1,text='Open File',command=temperature_function).grid(row=2,column=1,pady=4,sticky='w')


# ------------------------------------------------------------------
Label(tab1,anchor="w", text="Meteorological data:").grid(row=4,column=0,pady=4,sticky='w')
Button(tab1,text='Open File',command=meteo_function).grid(row=4,column=1,pady=4,sticky='w')

radiation_on= IntVar()
Checkbutton(tab1,text='Solar Radiation', variable=radiation_on, onvalue=1, offvalue=0).grid(row=3,column=0,pady=4,sticky='w') 
radiation_on.set(1)


Label(tab1, text="Height of the wind measurements (meters):").grid(row=5,column=0,pady=4,sticky='w')
height_wind = Entry(tab1, bd =3)
height_wind.insert(END,10)
height_wind.grid(row=5,column=1,pady=4)

Label(tab1, text="Wind direction contribution (°):").grid(row=6,column=0,pady=4,sticky='w')
contri_wind_var = StringVar(tab1)
contri_wind_var.set(20)
contri_wind = Spinbox(tab1, from_=0, to=180, textvariable=contri_wind_var)
contri_wind.grid(row=6,column=1,pady=4)

Label(tab1, text="Latitude (°):").grid(row=7,column=0,sticky='w',pady=4)
latitude_var = StringVar(tab1)
latitude_var.set(0)
latitude = Spinbox(tab1, from_=-180, to=180, textvariable=latitude_var)
latitude.grid(row=7,column=1,pady=4)

tka.Separator(tab1, orient=HORIZONTAL).grid(column=0, columnspan= 5, row=8, padx=10, pady=10, sticky='we')


# ---------------------------- Wind -------------------------------------------

Label(tab1,anchor="w",font="Verdana 8 bold", text="Basin length (based on wind direction):",width=50).grid(row=9,column=0,pady=8,sticky='w')

typechoose = StringVar()
typechoose.set(1)

rad1 = Radiobutton(tab1,text='Uniform', variable=typechoose, value=1, command=selected_len).grid(row=10,column=0,pady=4,sticky='w') 
rad2 = Radiobutton(tab1,text='File', variable=typechoose, value=2, command=selected_len).grid(row=11,column=0,pady=2,sticky='w')
 
Label(tab1, text="Basin length (meters):").grid(row=10,column=0,pady=4)
vari_len = Entry(tab1, bd =3)
vari_len.grid(row=10,column=1)
vari_len.config(state='normal')


file_len = Button(tab1,text='Open File',command=fetch_function)
file_len.grid(row=11,column=1,sticky='w')
file_len.config(state='disable')

# ------------------------------ Sensor level ---------------------------------
tka.Separator(tab1, orient=HORIZONTAL).grid(column=0, columnspan= 5, row=12, padx=10, pady=10, sticky='we')

Label(tab1,anchor="w", text="Sensor level:").grid(row=13,column=0,pady=4,sticky='w')
Button(tab1,text='Open File',command=senso_function).grid(row=13,column=1,pady=4,sticky='w')

Label(tab1, text="Level of reference (m)").grid(row=14,column=0,pady=4,sticky='w')
reference_level = Entry(tab1, bd =3)
reference_level.insert(END,0)
reference_level.grid(row=14,column=1,pady=4)

# ------------------------------ Water level data -----------------------------
tka.Separator(tab1, orient=HORIZONTAL).grid(column=0, columnspan= 5, row=15, padx=10, pady=10, sticky='we')

Label(tab1,anchor="w",font="Verdana 8 bold", text="Water level data:",width=50).grid(row=16,column=0,pady=8,sticky='w')

typelevel = StringVar()
typelevel.set(1)

Radiobutton(tab1,text='Uniform', variable=typelevel, value=1, command=selected_level).grid(row=17,column=0,pady=4,sticky='w') 
Radiobutton(tab1,text='File', variable=typelevel, value=2, command=selected_level).grid(row=18,column=0,pady=2,sticky='w')
 
Label(tab1, text="Water level (meters):").grid(row=17,column=0,pady=4)
unif_level = Entry(tab1, bd =3)
unif_level.grid(row=17,column=1)
unif_level.config(state='normal')


file_level = Button(tab1,text='Open File',command=level_function)
file_level.grid(row=18,column=1,sticky='w')
file_level.config(state='disable')


# --------------- Spectral Analysys and Definitions   (tab 2) -----------------

Label(tab2,anchor="w",font="Verdana 8 bold", text="Density structure",width=50).grid(row=1,column=0,sticky='w',pady=8)

Label(tab2, text="Metalimnion threshold (kg/m³/m):").grid(row=2,column=0,sticky='w',pady=4)
meta = Entry(tab2, bd =3)
meta.insert(END,0.1)
meta.grid(row=2,column=1,pady=4)

tka.Separator(tab2, orient=HORIZONTAL).grid(column=0, columnspan= 5, row=3, padx=10, pady=8, sticky='we')

# -----------------------------------------------------------------------------

Label(tab2,anchor="w",font="Verdana 8 bold", text="Band pass filter (Cutoff Frequency)",width=50).grid(row=4,column=0,pady=8,sticky='w')

filterchoose = StringVar()
rad1 = Radiobutton(tab2,text='Defined by the Internal wave model', variable=filterchoose, value=1, command=selected_filter).grid(row=5,column=0,pady=4,sticky='w') 
rad2 = Radiobutton(tab2,text='Defined manually: ', variable=filterchoose, value=2, command=selected_filter).grid(row=6,column=0,pady=2,sticky='w')
 
filterchoose.set(1)

Label(tab2, text="High-frequency band period (hour)").grid(row=7,column=0,pady=3,padx=100,sticky='w')
lv1 = Entry(tab2, bd =3)
lv1.grid(row=7,column=1)
lv1.config(state='disable')

Label(tab2, text="Low-frequency band period (hour)").grid(row=8,column=0,pady=3,padx=100,sticky='w')
hv1 = Entry(tab2, bd =3)
hv1.grid(row=8,column=1)
hv1.config(state='disable')

# -----------------------------------------------------------------------------

tka.Separator(tab2, orient=HORIZONTAL).grid(column=0, columnspan= 5, row=11, padx=10, pady=8, sticky='we')

Label(tab2,anchor="w",font="Verdana 8 bold", text="Spectral analysis parameter",width=50).grid(row=12,column=0,pady=8,sticky='w')

Label(tab2, text="Fourier window function:").grid(row=13,column=0,pady=4,padx=100,sticky='w')

welch_windows = ["Hamming", "Hann", "Blackman","Flattop"]
windows_var = StringVar(tab2)
windows_var.set(welch_windows[0])

windows_options = OptionMenu(*(tab2,windows_var)+tuple(welch_windows))
windows_options.grid(row=13,column=1,pady=4,sticky='w')


Label(tab2, text="Wavelet mother function:").grid(row=14,column=0,pady=4,padx=100,sticky='w')

wavelet_windows = ["Morlet", "Paul", "DOG"]
mother_var = StringVar(tab2)
mother_var.set(wavelet_windows[0])

mother_options = OptionMenu(*(tab2,mother_var)+tuple(wavelet_windows))
mother_options.grid(row=14,column=1,pady=4,sticky='w')

# -----------------------------------------------------------------------------


tka.Separator(tab2, orient=HORIZONTAL).grid(column=0, columnspan= 5, row=15, padx=10, pady=8, sticky='we')

Label(tab2,anchor="w",font="Verdana 8 bold", text="Window Size (Welch-averaging)",width=50).grid(row=16,column=0,pady=8,sticky='w')

winsizechoose = StringVar()
rad1 = Radiobutton(tab2,text='Defined by the Internal wave model (10 times the V1H1 mode)', variable=winsizechoose, value=1, command=selected_windowsize).grid(row=17,column=0,pady=4,sticky='w') 
rad2 = Radiobutton(tab2,text='Defined manually (days): ', variable=winsizechoose, value=2, command=selected_windowsize).grid(row=18,column=0,pady=2,sticky='w')
 
winsizechoose.set(1)


winsize = Entry(tab2, bd =3)
winsize.grid(row=18,column=1)
winsize.config(state='disable')

# ----------------------------- Isotherm Analsysis (tab3) ---------------------


isotherm_button= IntVar()
che_iso = Checkbutton(tab3,font="Verdana 8 bold",text='Isotherms analysis', variable=isotherm_button, onvalue=1, offvalue=0, command=selected_isotherms).grid(row=1,column=0,pady=4,sticky='w') 
isotherm_button.set(1)

isotherma = IntVar()
isothermb = IntVar()
isothermc = IntVar()
isothermd = IntVar()

iso1 = Checkbutton(tab3,text='Isotherm 1 (°C)', variable=isotherma, onvalue=1, offvalue=0, command=selected_isoa)
iso2 = Checkbutton(tab3,text='Isotherm 2 (°C)', variable=isothermb, onvalue=1, offvalue=0, command=selected_isob)
iso3 = Checkbutton(tab3,text='Isotherm 3 (°C)', variable=isothermc, onvalue=1, offvalue=0, command=selected_isoc) 
iso4 = Checkbutton(tab3,text='Isotherm 4 (°C)', variable=isothermd, onvalue=1, offvalue=0, command=selected_isod)

iso1.grid(row=2,column=0,pady=2,padx=100,sticky='w') 
iso2.grid(row=3,column=0,pady=2,padx=100,sticky='w') 
iso3.grid(row=4,column=0,pady=2,padx=100,sticky='w') 
iso4.grid(row=5,column=0,pady=2,padx=100,sticky='w') 

isoa = Entry(tab3, bd =3)
isoa.grid(row=2,column=1)
isoa.config(state='disable')

isob = Entry(tab3, bd =3)
isob.grid(row=3,column=1)
isob.config(state='disable')

isoc = Entry(tab3, bd =3)
isoc.grid(row=4,column=1)
isoc.config(state='disable')

isod = Entry(tab3, bd =3)
isod.grid(row=5,column=1)
isod.config(state='disable')



Label(tab3, text="Comparison between isotherms :").grid(row=6,column=0,pady=8,padx=100,sticky='w')


Label(tab3, text="Comparison 1:").grid(row=7,column=0,pady=3,padx=100,sticky='w')
Label(tab3, text="Comparison 2:").grid(row=8,column=0,pady=3,padx=100,sticky='w')

comparison_options = ["None", "Iso. 1-2", "Iso. 1-3", "Iso. 1-4", "Iso. 2-3", "Iso. 2-4", "Iso. 3-4"]
compa1_var = StringVar(tab3)
compa1_var.set(comparison_options[0])

compa1_options = OptionMenu(*(tab3,compa1_var)+tuple(comparison_options))
compa1_options.grid(row=7,column=1,pady=3,sticky='w')

compa2_var = StringVar(tab3)
compa2_var.set(comparison_options[0])

compa2_options = OptionMenu(*(tab3,compa2_var)+tuple(comparison_options))
compa2_options.grid(row=8,column=1,pady=3,sticky='w')
# ------------------------------------------------------------------------------

sensor_button= IntVar()
che_sen = Checkbutton(tab3,font="Verdana 8 bold",text='Sensor analysis', variable=sensor_button, onvalue=1, offvalue=0, command=selected_sensor).grid(row=9,column=0,pady=8,sticky='w') 
sensor_button.set(1)

sensora = IntVar()
sensorb = IntVar()
sensorc = IntVar()
sensord = IntVar()

sensora.set(1)
sensorb.set(1)
sensorc.set(1)


sen1 = Checkbutton(tab3,text='Sensor 1 (-)', variable=sensora, onvalue=1, offvalue=0, command=selected_sena)
sen2 = Checkbutton(tab3,text='Sensor 2 (-)', variable=sensorb, onvalue=1, offvalue=0, command=selected_senb)
sen3 = Checkbutton(tab3,text='Sensor 3 (-)', variable=sensorc, onvalue=1, offvalue=0, command=selected_senc) 
sen4 = Checkbutton(tab3,text='Sensor 4 (-)', variable=sensord, onvalue=1, offvalue=0, command=selected_send)

sen1.grid(row=10,column=0,pady=2,padx=100,sticky='w') 
sen2.grid(row=11,column=0,pady=2,padx=100,sticky='w') 
sen3.grid(row=12,column=0,pady=2,padx=100,sticky='w') 
sen4.grid(row=13,column=0,pady=2,padx=100,sticky='w') 

sena = Entry(tab3, bd =3)
sena.grid(row=10,column=1)
sena.config(state='normal')
sena.insert(END,1)

senb = Entry(tab3, bd =3)
senb.grid(row=11,column=1)
senb.config(state='normal')
senb.insert(END,2)

senc = Entry(tab3, bd =3)
senc.grid(row=12,column=1)
senc.config(state='normal')
senc.insert(END,3)

send = Entry(tab3, bd =3)
send.grid(row=13,column=1)
send.config(state='disable')


smooth_button= IntVar()
Checkbutton(tab3,font="Verdana 8 bold",text='Smooth results (recommended for large data analysis)', variable=smooth_button, onvalue=1, offvalue=0).grid(row=14,column=0,pady=8,sticky='w') 
smooth_button.set(0)


Label(tab4,anchor="w", text="Output folder:").grid(row=1,column=0,pady=4,sticky='w')
Button(tab4,text='Open File',command=output_folder).grid(row=1,column=1,pady=4,sticky='w')

Label(tab4, text="Output images quality (dpi):").grid(row=2,column=0,sticky='w',pady=4)
dpi = StringVar(tab4)
dpi.set(800)
latitude = Spinbox(tab4, from_=50, to=1000, textvariable=dpi)
latitude.grid(row=2,column=1,pady=4)

Button(tab4,font="Verdana 8 bold",text='Run',command=export_data, height = 4, width = 20).grid(row=3,column=1,pady=8,sticky='w')


window.config(menu = menubar)
window.mainloop()