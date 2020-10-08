# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer load data module (input module)

modififications to previosly versions must be specified here using WAVE codes:
    
W-23.03-1.00.3-00
A-01.03-1.00.3-00
V-22.03-1.00.3-00
E-05.03-1.00.3-00
"""
import numpy as np


def file_len(fname, col='off'):
#
#   Internal Function: Find the file length
#   Function to return the length of data
#
    fname = fname.replace('\n','')
    with open(fname, errors='replace') as f:
        line = f.readline()

        for i, l in enumerate(f):
            pass
        
    if col == 'off':
        return i+1
    
    elif col ==  'on':
        
        return i+1, len(line.split())-5

def serial_number (y,mo,d,h,m):
#
#   Internal Function: Serial number to define date
#   Function to return a serial number
#     
    return m + h*60 + d*24*60 + mo*30*24*60 + y*12*30*24*60

def find_nearest(array,value):
#
#   Internal Function: Find index of nearest value
#   Function to return the index (nearest value)  
# 
    return (np.abs(array-value)).argmin()

        
def level(sen, qt):
#
#   Read Function: Sensor distance
#   Function to read the type of sensor and the distance to the reference 
#    
    sensor_level = np.zeros(qt,float)
    speci        = np.zeros(qt,int)
        
    with open(sen, errors='replace') as f:
        
        i = 0
        next(f)
        
        for line in f:          
            line = line.strip()
            camp = line.split('\t')
        
            sensor_level[i] = float(camp[1])
            speci[i]        =   int(camp[2])
            
            i = i + 1

    return sensor_level, speci

def temper_read(nam,lin,qt):
#
#   Read Function: Temperature data
#   Function to read the temperature data  
#     
    yea  = np.zeros(lin,float)
    mon  = np.zeros(lin,float)
    day  = np.zeros(lin,float)
    hou  = np.zeros(lin,float)
    miu  = np.zeros(lin,float)

    temp     = np.zeros((lin,qt),float)
    serial_t = np.zeros(lin,float)
    
    comp_date = ["" for x in range(lin)]

    nam = nam.replace('\n','')
    with open(nam, errors='replace') as f:

        next(f)
        t=0
        for line in f:
            
            line = line.strip()   
            camp = line.split('\t')
                       
            yea[t] = float(camp[0]) 
            mon[t] = float(camp[1])
            day[t] = float(camp[2])
            hou[t] = float(camp[3])
            miu[t] = float(camp[4])                
                            
            comp_date[t] = camp[0]+'/'+camp[1]+'/'+camp[2]+'/'+camp[3]+'/'+camp[4]
            
            serial_t[t] = serial_number(yea[t],mon[t],day[t],hou[t],miu[t])
            
            for z in range(qt):
                temp[t][z] = float(camp[z+5])
    
            t=t+1
        dt = (serial_t[2]-serial_t[1])/60
            
    return comp_date, temp, serial_t, dt
        
               
def cot_read(nac,ean_serie,ean_cota):
#
#   Read Function: Water level data
#   Function to read the elevation data
# 
#   Elevation values does not need to be exacly in the region of sensors,
#   the function will find the best values to match the temperature data
#
    global cot,serial_c
        
    if ean_serie == 1:
        cot = ean_cota
    
    elif ean_serie == 2:
        nac = nac.replace('\n','')
       
        lic = file_len(nac)
    
        cye = np.zeros(lic,float)
        cmo = np.zeros(lic,float)
        cda = np.zeros(lic,float)
        cho = np.zeros(lic,float)
        cmi = np.zeros(lic,float)
        cot = np.zeros(lic,float)
    
        serial_c = np.zeros(lic,float)

        with open(nac, errors='replace') as f:

            next(f)
            t=0
            for line in f:   
                line = line.strip()
                camp = line.split('\t')
        
                cye[t] = float(camp[0]) 
                cmo[t] = float(camp[1])
                cda[t] = float(camp[2])      
                cho[t] = float(camp[3])
                cmi[t] = float(camp[4])
            
                serial_c[t] = serial_number(cye[t],cmo[t],cda[t],cho[t],cmi[t])
        
                cot[t] = float(camp[5])
      
                t=t+1
   
     
def win_read(win,radon):
#
#   Read Function: Wind data
#   Function to read information from meteorological station
# 
    global serial_w, wot, wdi, rad
    win = win.replace('\n','')
    wic = file_len(win)

    wye = np.zeros(wic,float)
    wmo = np.zeros(wic,float)
    wda = np.zeros(wic,float)
    who = np.zeros(wic,float)
    wmi = np.zeros(wic,float)
    wot = np.zeros(wic,float)
    wdi = np.zeros(wic,float)
    rad = np.zeros(wic,float)
    
    serial_w = np.zeros(wic,float)

    with open(win) as f:

        next(f)
        t=0
        for line in f:   
            line = line.strip()
            camp = line.split('\t')
            
            wye[t] = float(camp[0]) 
            wmo[t] = float(camp[1])
            wda[t] = float(camp[2])

            who[t] = float(camp[3])
            wmi[t] = float(camp[4])
            
            wot[t] = float(camp[5])   
            wdi[t] = float(camp[6])
            
            if radon == 1:
                rad[t] = float(camp[7])
            
            serial_w[t] = serial_number(wye[t],wmo[t],wda[t],who[t],wmi[t])
            
            t=t+1    

def windfetch (wd,name_fetch):
#
#   Read Function: Wind-fetch (lake length)
#   Function to read information from wind-fetch file (lake length)
#    
    name_fetch = name_fetch.replace('\n','')
    
    len_fetch = file_len(name_fetch)
    dists = np.zeros(len_fetch,float)
    angle = np.zeros(len_fetch,float)
    
    with open(name_fetch,errors='replace') as f:

        next(f)
        t=0
        for line in f:   
            line = line.strip()
            camp = line.split('\t')
            angle[t]    = float(camp[0])
            dists[t]    = float(camp[1])  
            t = t + 1
            
    return angle, dists
               

def serial_cota(serial_t,nac,lin,qt,temp,sen,ean_serie,ean_cota,z0):
#
#   Serial Function: Water level data
#   Function to match the elevation with temperature values and to adjust 
#   the sensors depth with water level 
# 
    cot_read(nac,ean_serie,ean_cota)
    
    ean   = np.zeros(lin,float)
    h     = np.zeros((lin,qt),float)
    
    tempa = temp
    sensor, condition = level(sen, qt)

    for t in range(lin):
        
        if ean_serie == 1:
            ean[t] = cot
        else:
            indice = find_nearest(serial_c,serial_t[t])
            ean[t] = cot[indice]

        for z in range(qt):
            
            if condition[z] == 1:
                h[t][z] = ean[t] - sensor[z]
            elif condition[z] == 2:
                h[t][z] = z0 + sensor[z]
            if h[t][z] > h[t][z-1] and z>0:
                
                aux        = h[t][z]   # check it !! for me it is wrong
                h[t][z]    = h[t][z-1]
                h[t][z-1]  = aux
                
                aux            = tempa[t][z]
                tempa[t][z]    = tempa[t][z-1]
                tempa[t][z-1]  = aux
    
        
    return ean, h, tempa
 
def serial_wind(serial_t,win,radon,lin):
#
#   Serial Function: Wind data
#   Function to match the wind measurements with temperature values 
#     
    win_read(win,radon)
    
    iw = np.zeros(lin,float)
    dw = np.zeros(lin,float)
    ra = np.zeros(lin,float)
        
    for t in range(lin):
        indice = find_nearest(serial_w,serial_t[t])
        iw[t] = wot[indice]
        dw[t] = wdi[indice]
        ra[t] = rad[indice]
               
    return iw,dw,ra


