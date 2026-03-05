# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer - Load variables module

Interwave Analyzer - Version 2 (2026) 
Load variables module version: 2.260305

-------------------------------------------------------------------------------

de Carvalho Bueno, R; Bleninger, T. B.; Lorke, A. 
Internal wave analyzer for thermally stratified lakes 
Environmental Modelling & Software, Elsevier, 2020 


Developed by Rafael de Carvalho Bueno 
https://buenorc.github.io/ 

Improvements and betterments by 
Andreas Lorke & Tobias Bleninger 

Report problems and improvements to email adresss below 
decarvalhobueno@gmail.com

for more information, see: 
https://buenorc.github.io/pages/interwave.html

"""
import os
import numpy as np

import iwpars as pars   # package of addition parameters
import iwmod  as mod    # package of functions

# Internal Function - Find the file length
def file_len(fname, col='off'):

    fname = fname.replace('\n','')
    with open(fname, errors='replace') as f:
        line = f.readline()

        for i, l in enumerate(f):
            pass
        
    if col == 'off':
        return i+1
    
    elif col ==  'on':
        
        return i+1, len(line.split())-5  # return the length of data

# Internal Function - Serial number to define date (min)
def serial_number (y,mo,d,h,m):
 
    return m + h*60 + d*24*60 + mo*30*24*60 + y*12*30*24*60 # return a serial number

# Internal Function - Find index of nearest value
def find_nearest(array,value):

    return (np.abs(array-value)).argmin() # return the index (nearest value)  

# Read Function - Sensor distance
def level(sen, qt):
    
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

    return sensor_level, speci # read the type of sensor and the distance to the reference 

# Read Function - Temperature data
def temper_read(nam,lin,qt):
    
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
        dt = (serial_t[2]-serial_t[1])/60 # hour
            
    return comp_date, temp, serial_t, dt # read the temperature data  
        
# Read Function - Water level data (module to read the elevation data)
def cot_read(nac,ean_serie,ean_cota):
    
#   Elevation values does not need to be exacly in the region of sensors,
#   the function will find the best values to match the temperature data

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
        
        
   
# Read Function - Wind data (read information from meteorological station)
def win_read(win,radon):

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

# Read Function - Loads file with format {H(m) | Dist(m) | Ref(m)}
def loadLen(filepath):

    depths = []
    dists = []
    refs = []

    with open(filepath, errors='replace') as f:
        next(f)  # skip header

        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')

            depths.append(float(parts[0]))
            dists.append(float(parts[1]))
            refs.append(float(parts[2]))

    return {
        "depths": np.array(depths),
        "dists": np.array(dists),
        "refs": np.array(refs)
    }

# Read Function - Basin data depending on final type_length
def loadData(type_length, fna, fna_trans, len_basin, change_basin, maxDepth):

    longData = None
    transData = None

    # Type 1: uniform swimming pool case
    if type_length == 1:
    
        # Create generic vertical structure
        dz = 0.2  # 0.2 m resolution 
        depths = np.linspace(dz, maxDepth + dz, 4)
    
        longData = {
            "depths": depths,
            "dists": np.full_like(depths, len_basin, dtype=float),
            "refs": np.zeros_like(depths, dtype=float),
            "orientation": 270 
        }

    # Type 2: single file (longitudinal)
    elif type_length == 2:

        longData = loadLen(fna)
        longData["orientation"] = 270

    # Type 3: longitudinal + transverse
    elif type_length == 3:

        # Load main
        longData = loadLen(fna)

        # Apply user-defined orientation
        orientation_main = change_basin if change_basin else 270
        longData = mod.basinOrientation(longData, orientation_main)

        # Load transverse
        transData = loadLen(fna_trans)

        # Transverse = +90° from main
        orientation_trans = (orientation_main + 90) % 360
        transData = mod.basinOrientation(transData, orientation_trans)

    return longData, transData

# Internal function - Determines final type_length and associated files
def lengthType(type_length, fna, additional_params):

    fna_trans = None
    warnTransv = False
    missFile = None

    change_basin = 270  # default orientation

    if type_length == 2:

        path_bathy = pars.extractPathBathy(additional_params)

        if path_bathy:
            base_dir = os.path.dirname(fna)
            candidate_path = os.path.join(base_dir, path_bathy)

            if os.path.isfile(candidate_path):
                fna_trans = candidate_path
                type_length = 3
                change_basin = pars.extractChangeBasin(additional_params)
            else:
                warnTransv = True
                missFile = path_bathy

    return type_length, fna, fna_trans, change_basin, warnTransv, missFile
               
# Serial Function - Water level data
def serial_cota(serial_t,nac,lin,qt,temp,sen,ean_serie,ean_cota,z0):

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
 
# Serial Function - Wind data
def serial_wind(serial_t,win,radon,lin):
  
    win_read(win,radon)
    
    iw = np.zeros(lin,float)
    dw = np.zeros(lin,float)
    ra = np.zeros(lin,float)
        
    for t in range(lin):
        indice = find_nearest(serial_w,serial_t[t])
        iw[t] = wot[indice]
        dw[t] = wdi[indice]
        ra[t] = rad[indice]
               
    return iw,dw,ra #   match the wind measurements with temperature values 


