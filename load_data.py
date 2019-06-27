# -*- coding: utf-8 -*-
"""
@Created on Thu Jul 20 19:07 2017
@ Author: BUENO, R. (RAFAEL BUENO)


    
"""
import numpy as np


#-----------------------------------------------------------------------------
#    subfunctions
#-----------------------------------------------------------------------------


def file_len(fname, col='off'):
    # find the sample number of a file and the (for col='on'), the number of columns
    # one line of reader is skiped 
    # 
    fname = fname.replace('\n','')
    with open(fname) as f:
        line = f.readline()

        for i, l in enumerate(f):
            pass
        
    if col == 'off':
        return i+1
    
    elif col ==  'on':
        
        return i+1, len(line.split())-5

def serial_number (y,mo,d,h,m):
    # write a serial number for a date and time ('y'ear,'mo'nth,'d'ay,
    #                                                           'h'our,'m'in)
    num = m + h*60 + d*24*60 + mo*30*24*60 + y*12*30*24*60
    return num

def find_nearest(array,value):
    # find the nearest value of a 'value' from a 'array'
    #
    idx = (np.abs(array-value)).argmin()
    return idx

#-----------------------------------------------------------------------------
#    read functions
#-----------------------------------------------------------------------------
def temper_read(nam,lin,qt):

    #
    # function to read the temperature data from files of vossoroca res.
    # qt  = sensors number (temperature chain)
    # dt  = temporal resolution
    # nam = file name of temperature data
    # lin = number of data in time 
    
    yea  = np.zeros(lin,float)
    mon  = np.zeros(lin,float)
    day  = np.zeros(lin,float)
    hou  = np.zeros(lin,float)
    miu  = np.zeros(lin,float)

    
    temp = np.zeros((lin,qt),float)
    
    serial_t = np.zeros(lin,float)
    
    comp_date = ["" for x in range(lin)]

    nam = nam.replace('\n','')
    with open(nam) as f:

        next(f)
        t=0
        for line in f:
            
            line = line.strip()   
            camp = line.split('\t')          # teste.csv
            
            
            yea[t] = float(camp[0]) 
            mon[t] = float(camp[1])
            day[t] = float(camp[2])
            hou[t] = float(camp[3])
            miu[t] = float(camp[4])                
            
                
            comp_date[t] = camp[0]+'/'+camp[1]+'/'+camp[2]+'/'+camp[3]+'/'+camp[4]
            
            serial_t[t] = serial_number(yea[t],mon[t],day[t],hou[t],miu[t])
            
            for z in range(qt):
                temp[t][z] = float(camp[z+5]) # camp[0:1] = date information
    
            t=t+1
        dt = (serial_t[2]-serial_t[1])/60     # dt (hours)
            
    return comp_date, temp, serial_t, dt
        
        

       
def cot_read(nac):
    #
    # function to read the elevation data from files of vossoroca res.
    # cot[t] - float 1d dimension - values of water level in meters
    #          cot[t] is referenced at 800 meters (last sensor - 1 meters 
    #                                   above the bottom of the reservoir)
    #
    # ps: elevation values does not need to be exacly in the region of 
    #     sensors, the function will find the best values to match the temp
    #     data.
    #
    global cot,serial_c
    
    nac = nac.replace('\n','')
    lic = file_len(nac)
    
    cye = np.zeros(lic,float)
    cmo = np.zeros(lic,float)
    cda = np.zeros(lic,float)
    cho = np.zeros(lic,float)
    cmi = np.zeros(lic,float)
    cot = np.zeros(lic,float)
    
    serial_c = np.zeros(lic,float)

    with open(nac) as f:

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
    


def serial_cota(serial_t,nac,lin,qt,temp,sen):
    #
    # Function to match the temperature values with elevation values and
    # to adjust the depth of the sensors with the water level
    #
    # h    [t:lin][z:qt] - depth of the sensors (ref=800m - last sensor) 
    # tempa[t:lin][z:qt] - Temp data adjusted with water level (°C)
    #    
    cot_read(nac)
    
    ean   = np.zeros(lin,float)
    h     = np.zeros((lin,qt),float)
    
    tempa = temp
    sensor, condition = level(sen, qt)

    for t in range(lin):
        
        indice = find_nearest(serial_c,serial_t[t])
        ean[t] = cot[indice]

        for z in range(qt):
            
            if condition[z] == 1:
                h[t][z] = ean[t] - sensor[z]
            elif condition[z] == 2:
                h[t][z] = 799 + sensor[z]
            if h[t][z] > h[t][z-1] and z>0:
                
                aux        = h[t][z]
                h[t][z]    = h[t][z-1]
                h[t][z-1]  = aux
                
                aux            = tempa[t][z]
                tempa[t][z]    = tempa[t][z-1]
                tempa[t][z-1]  = aux
            

    return ean, h, tempa
 
def serial_wind(serial_t,win,lin):
    #
    # Function to match the temperature values with the meteorological data
    #
    # iw[t:lin] - Wind intensity (m/s) 
    # dw[t:lin] - Wind direction (°)
    # ra[t:lin] - Solar radiation (W/m²)
    
    
    
    win_read(win)
    
    iw = np.zeros(lin,float)
    dw = np.zeros(lin,float)
    ra = np.zeros(lin,float)
    
    
    for t in range(lin):
        indice = find_nearest(serial_w,serial_t[t])
        iw[t] = wot[indice]
        dw[t] = wdi[indice]
        ra[t] = rad[indice]
        
    
        
    return iw,dw,ra


def windfetch (wd,name_fetch):
    
    name_fetch = name_fetch.replace('\n','')
    
    len_fetch = file_len(name_fetch)
    dists = np.zeros(len_fetch,float)
    angle = np.zeros(len_fetch,float)
    
    with open(name_fetch) as f:

        next(f)
        t=0
        for line in f:   
            line = line.strip()
            camp = line.split('\t')
            angle[t]    = float(camp[0])
            dists[t]    = float(camp[1])  
            t = t + 1
            
    return angle, dists
            
        
def win_read(win):
    #
    # Function to read the meteorological data
    #
    # the data file need to be formated to run this fuction. (\t) 
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
            rad[t] = float(camp[7])
            
            serial_w[t] = serial_number(wye[t],wmo[t],wda[t],who[t],wmi[t])
            
            t=t+1 
        
def level(sen, qt):
    
    sensor_level = np.zeros(qt,float)
    speci        = np.zeros(qt,int)
        
    with open(sen) as f:
        
        i = 0
        next(f)
        
        for line in f:          
            line = line.strip()
            camp = line.split('\t')
        
            sensor_level[i] = float(camp[1])
            speci[i]        =   int(camp[2])
            
            i = i + 1

    return sensor_level, speci