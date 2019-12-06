# -*- coding: utf-8 -*-
"""
@
@ INTERWAVE ANALYZER 
@
@ Program dedicated to interpret the dynamic of thermal stratified basins.
@
@ Created on July 2017
@ Author: BUENO, R. (RAFAEL DE CARVALHO BUENO)
            
"""

# -----------------------------------------------------------------------------
# PACKAGES AND MODULES 
# -----------------------------------------------------------------------------
#
# From Interwave Analyzer Library
#

import graph_package as graph            # package of graphs
import internal_module as mod            # package of functions
import modeling_internalwave   as miw    # package of internal wave model
import load_data as load                 # package of inputs (tab delimiter)
#
# From Python Library
#
import os 
import sys
import time
import datetime
import numpy as np

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader

from tkinter import *
from tkinter.filedialog import *
from tkinter.font import *

class StdoutRedirector(object):

    def __init__(self, text_area):
        self.text_area = text_area

    def write(self, str):
        self.text_area.insert(END, str)
        self.text_area.see(END)

     
def main():
   
    from reportlab.pdfgen import canvas
    
    start_time = time.time()
    
    old_stdout = sys.stdout
    
    root = Tk()
    root.configure(background='white')
    root.title("Interwave Analyzer Running") 
    root.geometry('600x600')
    
    
    outputPanel = Text(root, wrap='word', height=30, width=100)
    outputPanel.grid(column=0, row=0, columnspan = 2, sticky='NSWE', padx=5, pady=5)

    sys.stdout = StdoutRedirector(outputPanel)
    

    print ("> Interwave Analyzer is starting the data processing... ")
    root.update()
    print ("--------------------------------------------------------------------------------------")
    root.update()
    print ("> Interwave Analyzer, version 1.00.0        June 25  2019")
    root.update()  
    print ("> ")
    root.update() 
    print ("--------------------------------------------------------------------------------------")
    root.update()
    print ("> ")
    root.update() 
    print ("> Part I       Reading information from GUI... ")
    root.update()
# -----------------------------------------------------------------------------
#                         INPUTS AND CONFIGURATIONS
# -----------------------------------------------------------------------------

    with open('temporary.txt') as reader:

        nam = reader.readline()
        nac = reader.readline()
        win = reader.readline()
        sen_nam = reader.readline()
        sen_nam = sen_nam.replace('\n','') 
        rw  = float(reader.readline())
        linang = float(reader.readline())
        lat = float(reader.readline())
        depi = int(reader.readline())
        
        
        type_length = float(reader.readline())
        if type_length == 1:
            len_basin = float(reader.readline())
        elif type_length == 2:
            fna = reader.readline()
            
        minval = float(reader.readline()) 
        
# -------------------- frequency cutoff (bandpass filter)  --------------------

# if you whish to  filter the variation according to internal seiche model 
# leave filter_process = on 

# if filter_process is '2', filter defined will be used
# if                is '1', period of V1H1 model is used     
        
        filter_process = int(reader.readline())
        if filter_process == 2:
            low_per  = float(reader.readline())
            high_per = float(reader.readline())
        else:
            for i in range(2):
                next(reader)
              
        window = reader.readline()               
        mother = reader.readline()
        window = window.replace('\n','') 
        mother = mother.replace('\n','') 
    

        if window == 'Hamming':
            window = 'hamming'
        elif window == 'Hann':
            window = 'hann'
        elif window == 'Blackman':
            window =  'blackman'
        elif window == 'Flattop':
            window =  'flattop'
            
    

#-------------------------- Isotherms definition ------------------------------
            
        turn_iso = int(reader.readline())
        
        tan = np.zeros(4,float)
        tau = np.zeros(4,float)
        
        for i in range(4):
            numaux = float(reader.readline())
            
            if turn_iso == 0:
                tau[i] = -999
                tan[i] = 0
            else:
                tau[i] = numaux
                tan[i] = 1
                
        c1 = reader.readline()
        c2 = reader.readline()       
        
        ana1 = mod.comparison_definition(c1)
        ana2 = mod.comparison_definition(c2)

#---------------------------- Sensors definition ------------------------------

        turn_temp = int(reader.readline())  # if '1', the sensor analysis is 'on'
    
        sen = np.zeros(4,float)
        seu = np.zeros(4,int)
        
        for i in range(4):
            seu[i] = int(reader.readline())
            if turn_temp == 0 or seu[i]==-999:
                seu[i] = -999
                sen[i] = 0
            else:
                sen[i] = 1
 
        largelen = int(reader.readline()) # if '1', the data is smoothed to analyze spectrums
        output_path = reader.readline()   # folder that will be used to results
        output_path = output_path+'/'
    
    os.remove('temporary.txt')

 
    root.update()
# -----------------------------------------------------------------------------
# LOAD DATA FROM FILE
# -----------------------------------------------------------------------------

    print ("> Part II      Loading information from files... ")
    root.update()   
       
    lin, qt = load.file_len(nam,col='on') 
        
    date, temp, serialt, dt = load.temper_read(nam,lin,qt)
    ean, h, tempa           = load.serial_cota(serialt,nac,lin,qt,temp,sen_nam)
    wind, dw, ra            = load.serial_wind(serialt,win,lin)

    
    if type_length == 1: 
        dists = len_basin 
    elif type_length == 2:
        angle, dists  = load.windfetch(dw,fna)

    aux_iso = 1    # value to do the corehence analysis with wind and solar rad.


    print ("> Part III     Defining constants and creating variables... ")
    root.update()  
# -----------------------------------------------------------------------------
# DEFINING CONSTANTS 
# -----------------------------------------------------------------------------


    sal = 0              # salinity
    pre = 0              # pressure

    dj       = 0.1       # spatial resolution for countourf graph (m) 
    dci      = 1         # decimal of countouf graph (generally is refered to dj)   

    windlim  = 45        # parametrization of wind contribution (\pm 45° ~ 90°) 


    ddep   = 0.5         # periods for check the model sensitivity
    drho   = 0.1

    lh_sen = 100         # 0.01 for +- 0.5m from the reference thickness 
    lp_sen = 20          # 0.01 for +- drho kg/m³ from the reference density


# -----------------------------------------------------------------------------
# DEFINING ARRAYS 
# -----------------------------------------------------------------------------    

# general parameters
    glin = np.zeros((lin),float)     # reduced gravity (m/s²)
    wast = np.zeros((lin),float)     # friction velocity of the wind  (m/s) 
    s    = np.zeros((lin),float)     # instability frequency (Hz)
    n    = np.zeros((lin),float)     # stability frequency (rad*Hz)
    riw  = np.zeros((lin),float)     # wind richardson number (-) 
    frw  = np.zeros((lin),float)     # wind froude number (-)
    hast = np.zeros((lin),float)     # turbulent fluctuation at the pycnocline (m)                                        
    umal = np.zeros((lin),float)     # maximum shear velocity (m/(s.m))

    thermo_temp = np.zeros((lin),float) # temperature of the thermocline (°C)

# considering a thin metalimnion:
    ht     = np.zeros((lin),float)     # thermocline depth (m) from 800 m   
    he     = np.zeros((lin),float)     # epilimnion thickness (m)
    hh     = np.zeros((lin),float)     # hypolimnion thickness (m)     
    pe     = np.zeros((lin),float)     # average density of the epilimnion (kg/m³)
    ph     = np.zeros((lin),float)     # average density of the hypolimnion (kg/m³)

# 3 layer systems
    h1      = np.zeros((lin),float)     # 1st layer thickness (m)
    h2      = np.zeros((lin),float)     # 2nd layer thickness (m)
    h3      = np.zeros((lin),float)     # 3rd layer thickness (m)
    p1      = np.zeros((lin),float)     # average density of 1st layer (kg/m³)
    p2      = np.zeros((lin),float)     # average density of 2nd layer (kg/m³)
    p3      = np.zeros((lin),float)     # average density of 3rd layer (kg/m³)

    pu     = np.zeros((lin),float)     # "surface" density (1m above) (kg/m³)
    pd     = np.zeros((lin),float)     # "bottom" density (2m from btm) (kg/m³)
    pm     = np.zeros((lin),float)     # density difference (pd - pu) (kg/m³)
#    ze     = np.zeros((lin),float)     # depth of the epilimnion end (m)
#    zh     = np.zeros((lin),float)     # depth of the hypolimnion end (m)
    hH     = np.zeros((lin),float)     # ratio of he/H (-).
    
# parameters
    n        = np.zeros((lin),float)     # stability frequency (rad/Hz)
    riw      = np.zeros((lin),float)     # wind richardson number (-) 
    frw      = np.zeros((lin),float)     # wind froude number (-)
    hast     = np.zeros((lin),float)     # The vertical scale of turbulent 
    strs     = np.zeros((lin),float)     # Winds stress (N/m²)
    wedd     = np.zeros((lin),float)     # Wedderburn Number (-) 
    wedd_inv = np.zeros((lin),float)     # Inverse of the Wedderburn Number (dn/dh)

# isotherms
    isoa = np.zeros((lin),float)     # isotherm 'a' (m) 
    isob = np.zeros((lin),float)     # isotherm 'b' (m) 
    isoc = np.zeros((lin),float)     # isotherm 'c' (m) 
    isod = np.zeros((lin),float)     # isotherm 'd' (m) 

# classification

    genera = np.zeros((lin),float)   # class. of internal seiche generation
    iw_up = np.zeros((lin),float)    # superior limit for internal seiche formation 
    iw_dw = np.zeros((lin),float)    # inf limit for internal wave formation
    
# parameters in 2d (vector[time][depth])
    bvf2d = np.zeros((lin,qt-1),float)  # Brunt-Vaisalla frequency (rad*Hz)
    riw2d = np.zeros((lin,qt-1),float)  # Richardson Number (-)
    hzmid = np.zeros((lin,qt-1),float)  # z depth at middle point (m at ref 800 m)
    ins2d = np.zeros((lin,qt-1),float)  # Instability frequency (dw*/dz) (Hz)
    gli2d = np.zeros((lin,qt-1),float)  # Reduced gravity (m)


# model sensibility 
    v1h1_she = np.zeros((lh_sen),float)
    v1h1_shh = np.zeros((lh_sen),float)
    v1h1_spe = np.zeros((lp_sen),float)
    v1h1_sph = np.zeros((lp_sen),float)
    v2h1_sh2 = np.zeros((lh_sen),float)
    v2h1_sp2 = np.zeros((lp_sen),float)
    v3h1_sh3 = np.zeros((lh_sen),float)
    v3h1_sp3 = np.zeros((lp_sen),float)

    v1mode   = np.zeros((lin),float) 
    v2mode   = np.zeros((lin),float) 

# wind parametrization
    sdur  = np.zeros((lin),float) # wind stress with duration filter
    sdir  = np.zeros((lin),float) # wind stress with duration and direction

    newdw = np.zeros((lin),float) # wind direction for friction velocity of win
    newdi = np.zeros((lin),float)
    newdw_low = np.zeros((lin),float)
    
    ridu = np.zeros((lin),float) # Ri with duration filter
    ridi = np.zeros((lin),float) # Ri


# -----------------------------------------------------------------------------
#                  COMPUTATION OF PARAMETERS AND INDICES 
# -----------------------------------------------------------------------------  

    
    print ("> ")
    root.update()
    load_time = time.time()
    print ("> Execution time for part I, II, and III: ")
    root.update() 
    print ("> "+str(round(load_time-start_time,4))+' seconds')
    root.update()   
# --------------------- Compute the wind fetch [min, ave, max] ----------------

    iw           = mod.velocityten(wind, rw)          # convert wind for 10 m high
    dw_mean      = mod.wind_average(dw,iw)

    dw_min, dw_max = mod.wind_angle(dw_mean,linang)
    

# if basin is selected as integer, it is considered a variation of \pm 5% 
    if type_length == 1: 
        ls_fetch = [len_basin, 0.95*len_basin, 1.05*len_basin]
    elif type_length == 2:
        ls_fetch = [dists[mod.find_nearest(angle,dw_mean)], \
            dists[mod.find_nearest(angle,mod.keepout(dw_min))], \
            dists[mod.find_nearest(angle,mod.keepout(dw_max))]]


    ls_min = np.min(ls_fetch)
    ls_max = np.max(ls_fetch)
    ls_ave = np.average(ls_fetch)    

    ls_fetch = [ls_min,ls_ave,ls_max]   

# -----------------------------------------------------------------------------

    consecutive_dire, consecutive_wind, ver_dire, ver_wind = 0,0,0,0
    auxisa, auxisb, auxisc, auxisd = None, None, None, None
    print ("--------------------------------------------------------------------------------------")
    root.update()
    print('> ')
    root.update()
    print ("> Part IV      Starting computation... ")
    root.update()  
    print('> ')
    root.update()
# Output of functions:
#    
# vector aux - original vector = one dimension reduction (espatial reduction)
#
# auxtem          temp - [array float 2d] temperature     (°C)
# auxh            h    - [array float 2d] elevation       (m from 800m)
# auxean          ean  - [array float 1d] water level     (m from 800m)
# auxiw           iw   - [array float 1d] wind velocity   (m/s)
# auxdw           dw   - [array float 1d] wind direction  (°)
# auxra           ra   - [array float 1d] solar radiation (W/m²)

    for t in range(lin):

    
        auxtem = mod.vector_aux(tempa,t,qt)
        auxh   = mod.vector_aux(h,t,qt)
    
        auxiw  = iw[t]

        auxean = ean[t]

    
        max_ean = np.nanmax(auxean)  # elevation of the second sensor
        min_ean = 800
    
        if (turn_iso == 1):
            if (tau[0] != -999):
                isoa[t] = mod.isotherms(tau[0],qt,auxh,auxtem,max_ean,min_ean,auxisa)
                auxisa = isoa[t]
            if (tau[1] != -999):
                isob[t] = mod.isotherms(tau[1],qt,auxh,auxtem,max_ean,min_ean,auxisb)
                auxisb = isob[t]
            if (tau[2] != -999):
                isoc[t] = mod.isotherms(tau[2],qt,auxh,auxtem,max_ean,min_ean,auxisc)
                auxisc = isoc[t]
            if (tau[3] != -999):
                isod[t] = mod.isotherms(tau[3],qt,auxh,auxtem,max_ean,min_ean,auxisd)
                auxisd = isod[t]
    
    
    
    # 2layer structure
        he[t], hh[t], pe[t], ph[t], glin[t], n[t], pu[t], pd[t] = mod.structure2layer(qt, auxh, auxtem, sal, pre, auxean)      
    
        ht[t] = ean[t] - he[t]
        hH[t] = he[t]/(he[t]+hh[t])
        if (hH[t] > 0.5):
            hH[t] = 1 - hH[t]
        pm[t] = (pd[t]+pu[t])/2
    
        strs[t], wast[t], s[t], riw[t], frw[t], hast[t], umal[t] = mod.wind_parameters(auxiw, rw, pe[t], he[t], n[t], glin[t], auxean)
    

    
    # 2d parameters 
    #
    # output: p, n2d, hmid, s2d, and riw2d --> vector of dim 1 [z]   
        p, n2, hmid = mod.nonwind_para2d(qt,auxh,auxean,auxtem,sal,pre)
        gli2d       = mod.glinha_2d(p,qt)
        s2, riw2    = mod.wind_para2d(auxiw,rw,qt,auxh,pe[t],auxean,n2,hmid,p,gli2d)
    
        for z in range(qt-1):
            bvf2d[t][z] = n2[z]
            riw2d[t][z] = riw2[z]
            hzmid[t][z] = hmid[z]
            ins2d[t][z] = s2[z]
    
    
    # wedderburn function (he = ze for a two layer system)
        wedd[t] = mod.wedderburn(glin[t],he[t],wast[t],ls_fetch[1])  # ls is the fetch
        wedd_inv[t] = 1/wedd[t]
    
    # wind parameterization
        tbsiw = 2*ls_fetch[1]/np.sqrt(glin[t]*he[t]*hh[t]/(he[t]+hh[t])) # Merian
        newdi[t] = None
        newdw_low[t] = None
    
        if(wedd[t] < 20):
            newdw_low[t] = dw[t]
    
        if(wedd[t]<100):

            newdw[t] = dw[t]
        
            if(consecutive_wind > 0):
                wmin, wmax = mod.wind_angle(dw[t-1],windlim)
                if( dw[t] >= wmin and dw[t] <= wmax):
                    if(consecutive_dire==0):
                        tindex = t
                        consecutive_dire = consecutive_dire + 1

                    else:
                        dexmin, dexmax = mod.wind_angle(dw[tindex],windlim)
                        if(dw[t] >= wmin and dw[t] <= wmax):
                            consecutive_dire = consecutive_dire + 1
                            newdi[t] = dw[t]
                        else:
                            if(consecutive_dire > ver_dire):
                                ver_dire = consecutive_dire
                            consecutive_dire = 0
                            
                else:
                    if(consecutive_dire > ver_dire):
                        ver_dire = consecutive_dire
                
                    consecutive_dire = 0
        
            consecutive_wind = consecutive_wind + 1 
        
            if (consecutive_wind > ver_wind):
#                indiwind = t                 # last sample of the main wind 
                ver_wind = consecutive_wind  # number of samples of main wind
            
            if(consecutive_dire > ver_dire):
                ver_dire = consecutive_dire
        
        else:

            consecutive_dire = 0
            consecutive_wind = 0
            newdw[t] = None
        
    # 3layer structure    
        h1[t],h2[t],h3[t],p1[t],p2[t],p3[t] = \
        mod.structure3layer(qt, auxh, auxtem, sal, pre, minval, auxean)
    
        if h1[t] == 999:
            h1[t] = he[t]
            p1[t] = pe[t]
            h2[t] = 0.001*he[t]     
            p2[t] = (pe[t]+ph[t])/2
            h3[t] = hh[t]
            p3[t] = ph[t]
    
    # Generation of internal seiches according to Ri, hh, he, wind fetch 
    #
    # outputs: classification along time series [array int 1d] = 1, 2, 3 or 4
    #
    #    1 = Stratification is broken down by mixing
    #    2 = Large interface displacement occurs and is accompanied by interface shear and Kelvin-Helmholtz billows (strat is broken)
    #    3 = Internal seiche is dominant
    #    4 = Internal seiche with short amplitude.
    #
        genera[t] = mod.class_generation(riw[t],hh[t],he[t],ls_fetch[1])
        iw_dw[t], iw_up[t] = mod.iw_generation(riw[t],hh[t],he[t],ls_fetch[1])
    
        _,v1mode[t],_ = miw.disp_xmodel2(pe[t],ph[t],he[t],hh[t],ls_fetch,1)
        _,v2mode[t],_ = miw.disp_xmodel3(p1[t],p2[t],p3[t],h1[t],h2[t],h3[t],ls_fetch,2,1)
    
        if(t>0 and abs(v2mode[t] - v2mode[t-1])>5*60*60):
    # it is difficult accurs a change in this order (5h in a short dt)
            v2mode[t] = v2mode[t-1]
  


    iso = [isoa, isob, isoc, isod]

# generation periods (total period is divided by 3 shorter periods)
# h: he/H and Wedderburn number for each subperiod (confidence interval with 95%)
                                                   
    group = [genera[i:i+int(lin/3)] for i in range(0, len(genera), int(lin/3))]
    hH_gp = [hH[i:i+int(lin/3)] for i in range(0, len(hH), int(lin/3))]
    wi_gp = [wedd_inv[i:i+int(lin/3)] for i in range(0, len(wedd_inv), int(lin/3))]
    

    

    h_lim1 = mod.ci(hH_gp[0])
    W_lim1 = mod.ci(wi_gp[0])
    h_lim2 = mod.ci(hH_gp[1])
    W_lim2 = mod.ci(wi_gp[1])
    h_lim3 = mod.ci(hH_gp[2])
    W_lim3 = mod.ci(wi_gp[2])
    


#date - formating data
    dx    = [datetime.datetime.strptime(d, '%Y/%m/%d/%H/%M') for d in date]
    dx_gp = [dx[i:i+int(lin/3)] for i in range(0, len(dx), int(lin/3))]
    
           
# defining the subperiods
    P1    = [dx_gp[0][0], dx_gp[0][-1]]
    P2    = [dx_gp[1][0], dx_gp[1][-1]]
    P3    = [dx_gp[2][0], dx_gp[2][-1]]

# average of some variables
    m_pe   = np.nanmean(pe)
    m_ph   = np.nanmean(ph)
    m_he   = np.nanmean(he)
    m_hh   = np.nanmean(hh)
    m_h1   = np.nanmean(h1)
    m_h2   = np.nanmean(h2)
    m_h3   = np.nanmean(h3)
    m_p1   = np.nanmean(p1)
    m_p2   = np.nanmean(p2)
    m_p3   = np.nanmean(p3)
    m_glin = np.nanmean(glin)
    m_wast = np.nanmean(wast)
    m_n    = np.nanmean(n)
    m_riw  = np.nanmean(riw)
    m_ht   = np.nanmean(ht)
    m_newind = np.nanmean(newdw)

    print (">         Parameters were defined")
    root.update() 

# -----------------------------------------------------------------------------
#  correction of parameters according to wind filtering 
# -----------------------------------------------------------------------------

    dura = ver_wind*dt*60*60   # duration of consecutive wind events (seconds)
    dire = ver_dire*dt*60*60   # duration of wind events considering homogeneous 
                           # wind events (seconds)

    fdura = min([1,np.sqrt(4*dura/tbsiw)])   # parametrization (wind duration)
    fdire = min([1,np.sqrt(4*dire/tbsiw)])   # parametrization (wind direc + dura)

# loop in time to define the new richardson based on wind event 
    for t in range(lin):
    
        sdur[t]  = (fdura**2)*strs[t]
        sdir[t]  = (fdire**2)*strs[t]    

        if(fdura == 0):
            fdura = 1
        if(fdire == 0):
            fdire = 1

        ridu[t] = riw[t]/fdura
        ridi[t] = riw[t]/fdire
    
#
# averaged water density and layer thickness for higher vertical modes
#
    np2    = (m_p2+m_p1)/2 
    np3    = (m_p3+m_p2)/2
    nh1    = m_h1/2
    nh2    = (m_h1+m_h2)/2
    nh3    = (m_h2+m_h3)/2
    nh4    = m_h3/2

# ---------------------------------------------------------------------------- 
#  Internal seiche modeling (V1H1)    
# ---------------------------------------------------------------------------- 
# > internal wave period (in seconds)
#
#   pzvnhm[i] - z-momentum n-layered model for VnHm (nonhydrostatic model)
#   pxvnhm[i] - x-momentum n-layered model for VnHm (hydrostatic model)

# approximated internal wave velocity V1H1 (Merian's equation)
    cp = np.sqrt(m_glin*m_he*(1-m_he/(m_he+m_hh)))  

# theoretical periods

    pzv1h1 = np.array(miw.disp_zmodel(m_pe,m_ph,m_he,m_hh,ls_fetch,1))

    pxv1h1 = np.array(miw.disp_xmodel2(m_pe,m_ph,m_he,m_hh,ls_fetch,1))
    pxv1h2 = np.array(miw.disp_xmodel2(m_pe,m_ph,m_he,m_hh,ls_fetch,2))
    pxv1h3 = np.array(miw.disp_xmodel2(m_pe,m_ph,m_he,m_hh,ls_fetch,3))

    pxv2h1 = np.array(miw.disp_xmodel3(m_p1,m_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,2,1))
    pxv2h2 = np.array(miw.disp_xmodel3(m_p1,m_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,2,2))
    pxv2h3 = np.array(miw.disp_xmodel3(m_p1,m_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,2,3))

    pxv3h1 = np.array(miw.disp_xmodel4(m_p1,np2,np3,m_p3,nh1,nh2,nh3,nh4,ls_fetch,3,1))
    pxv3h2 = np.array(miw.disp_xmodel4(m_p1,np2,np3,m_p3,nh1,nh2,nh3,nh4,ls_fetch,3,2))
    pxv3h3 = np.array(miw.disp_xmodel4(m_p1,np2,np3,m_p3,nh1,nh2,nh3,nh4,ls_fetch,3,3))

#  internal wave frequency (without and with Earth rotation effect)
    fz_11  = 1/pzv1h1

    fx_11, cfx11  = 1/pxv1h1, 1/miw.coriolis_effect(lat,pxv1h1)
    fx_21, cfx21  = 1/pxv2h1, 1/miw.coriolis_effect(lat,pxv2h1)
    fx_31, cfx31  = 1/pxv3h1, 1/miw.coriolis_effect(lat,pxv3h1)

# model sensitivity using the hydrostatic model
    aux_he = m_he - ddep
    xhe    = np.linspace(aux_he, m_he+ddep, lh_sen)

    for i in range(lh_sen):
        _,aux,_     = np.array(miw.disp_xmodel2(m_pe,m_ph,aux_he,m_hh,ls_fetch,1))
        v1h1_she[i] = aux/60/60
        aux_he      = aux_he + 0.01


    aux_hh = m_hh - ddep
    xhh    = np.linspace(aux_hh, m_hh+ddep, lh_sen)

    for i in range(lh_sen):
        _,aux,_     = np.array(miw.disp_xmodel2(m_pe,m_ph,m_he,aux_hh,ls_fetch,1))
        v1h1_shh[i] = aux/60/60
        aux_hh      = aux_hh + 0.01


    aux_pe = m_pe - drho
    xpe    = np.linspace(aux_pe, m_pe+drho, lp_sen)

    for i in range(lp_sen):
        _,aux,_     = np.array(miw.disp_xmodel2(aux_pe,m_ph,m_he,m_hh,ls_fetch,1))
        v1h1_spe[i] = aux/60/60
        aux_pe      = aux_pe + 0.01


    aux_ph = m_ph - drho
    xph    = np.linspace(aux_ph, m_ph+drho, lp_sen)

    for i in range(lp_sen):
        _,aux,_     = np.array(miw.disp_xmodel2(m_pe,aux_ph,m_he,m_hh,ls_fetch,1))
        v1h1_sph[i] = aux/60/60
        aux_ph      = aux_ph + 0.01


    aux_h2 = m_h2 - ddep
    xh2    = np.linspace(aux_h2, m_h2+ddep, lh_sen)

    for i in range(lh_sen):
        _,aux,_     = \
        np.array(miw.disp_xmodel3(m_p1,m_p2,m_p3,m_h1,aux_h2,m_h3,ls_fetch,2,1))
        v2h1_sh2[i] = aux /60/60
        aux_h2      = aux_h2 + 0.01


    aux_p2 = m_p2 - drho
    xp2    = np.linspace(aux_p2, m_p2+drho, lp_sen)

    for i in range(lp_sen):
        _,aux,_     = np.array(miw.disp_xmodel3(m_p1,aux_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,2,1))
        v2h1_sp2[i] = aux /60/60
        aux_p2      = aux_p2 + 0.01


    aux_nh3 = nh3 - ddep
    xh3     = np.linspace(aux_nh3, nh3+ddep, lh_sen)

    for i in range(lh_sen):
        _,aux,_     = np.array(miw.disp_xmodel4(m_p1,np2,np3,m_p3,nh1,nh2,aux_nh3,nh4,ls_fetch,3,1))
        v3h1_sh3[i] = aux /60/60
        aux_nh3      = aux_nh3 + 0.01


    aux_np3 = np3 - drho
    xp3     = np.linspace(aux_np3, np3+drho, lp_sen)

    for i in range(lp_sen):
        _,aux,_     = np.array(miw.disp_xmodel4(m_p1,np2,aux_np3,m_p3,nh1,nh2,nh3,nh4,ls_fetch,3,1))
        v3h1_sp3[i] = aux /60/60
        aux_np3      = aux_np3 + 0.01

    print (">         Internal wave periods were estimated")
    root.update() 


# ----------------------------------------------------------------------------
# Theory of Generation of Basin-scale internal waves - probability (%)
# Output: couter (key)
#   
#   key i - class i = chances (%) for each subperiod P1, P2 and P3  
#
    genera_0 = [np.count_nonzero(group[0]==1),np.count_nonzero(group[0]==2),np.count_nonzero(group[0]==3),np.count_nonzero(group[0]==4)]
    genera_1 = [np.count_nonzero(group[1]==1),np.count_nonzero(group[1]==2),np.count_nonzero(group[1]==3),np.count_nonzero(group[1]==4)]
    genera_2 = [np.count_nonzero(group[2]==1),np.count_nonzero(group[2]==2),np.count_nonzero(group[2]==3),np.count_nonzero(group[2]==4)]

    for i in range(4):
        genera_0[i] = genera_0[i]*100/len(group[0]) 
        genera_1[i] = genera_1[i]*100/len(group[1])
        genera_2[i] = genera_2[i]*100/len(group[2]) 
      

    print (">         Periods were classified according to lake mixing")
    root.update() 
# spectral analysis of isotherms and data sensors -----------------------------

# output variables:
#
# (iso/sensor)time_(type)      - [array 1d] array signal time (hour)
# (iso/sensor)per_(type)       - [array 1d] scale 's' (hour)
# (iso/sensor)power_(type)     - [array 2d] wavelet power (unit of type**2)
# (iso/sensor)global_ws_(type) - [array 1d] global wavelet (unit of type**2) 
# sig(isotherm or type)        - [array 1d] significance level (95%)
#
# isotherm - a, b, c, and d 
# sensor   - each line represents one sensor data [array 2d] 
# type     - temp(temperature), sol(solar radiation), win(wind intensity)
#

# ------------------------------ spectral depth -------------------------------

    sensor_filtered = []
    timee           = []
    per             = []
    power_sensor    = []
    freqe           = []
    welch_sensor    = []
    

    order    = 6                 # define the order of the bandpass filter 
    if filter_process == 1 :  # define the range of the bandpass filter
        low1, high1 = 1/((pxv1h1[2]/60)/60), 1/((pxv1h1[0]/60)/60) # banded through internal wave model V1H1 (lower mode)

    else:
        low1, high1 = 1/high_per, 1/low_per
    
    if turn_temp == 1:
        
        
        for index in range(qt):   # spectral analysis for all sensors
        
            new = mod.vector_time (temp,index,lin)
            fs = 1/dt
            
            try:
                filtered_band = mod.butter_bandpass_filter(new, low1, high1, fs, order)
            except ValueError:
                filtered_band = None
            
            aux_time, aux_per, aux_power = mod.wave_spectral ( new, dt, mother)
            aux_freq, aux_wl_per, aux_welch = mod.welch_method(new,0,window, dt)
        
            sensor_filtered.append(filtered_band)

            timee.append(aux_time)
            per.append(aux_per)
            power_sensor.append(aux_power)
            freqe.append(aux_freq) 
            welch_sensor.append(aux_welch)
        
        
# -----------------------------------------------------------------------------
# spectral analysis: Wavelet analysis, Fourier analysis (Welch method), butter
#                    bandpass filter with order 3
#                        
# the window size for Fourier analysis (< windsize ~ > error for higher modes)
#
# isotherm correction: 
#    iso[isotherm]_corrected (m)  -  variation in meters of the isotherm 
#                                    corrected through the surface variation 
#
# -----------------------------------------------------------------------------
        
        
    windsize = 5*pxv2h1[1]             # seconds (5 times the wave periodicity of mode 3)
    # normalized water level variation 

    if filter_process == 1 :        # define the range of the bandpass filter
        lowcut, highcut  = 1/((pxv1h1[2]/60)/60), 1/((pxv1h1[0]/60)/60)
    else:
        lowcut, highcut  =  1/high_per, 1/low_per

    freq      = []
    power     = []
    time_temp = []
    per_temp  = []
    band      = []
    welch     = []

    iso_corrected = np.zeros(4,float)


    if (turn_iso == 1):
        ean_norm  = ean  - np.mean(ean)
        for i in range(4):
            if(tau[i] != -999):
                aux_time, aux_per, aux_power = mod.wave_spectral (iso[i], dt, mother)
            
                aux_freq, _, aux_welch = mod.welch_method(iso[i],windsize,window, dt)
                fs = 1/dt
                
                try:
                    aux_band = mod.butter_bandpass_filter(iso[i], lowcut, highcut, fs, order)
                except ValueError:
                    aux_band = None
                
                iso_corrected[i] = np.mean(abs((iso[i] - np.mean(iso[i])) - ean_norm))
    
                time_temp.append(aux_time)
                per_temp.append(aux_per)
                power.append(aux_power)
                freq.append(aux_freq) 
                welch.append(aux_welch)  
                band.append(aux_band)
            else:
                time_temp.append(0)
                per_temp.append(0)
                power.append(0)
                freq.append(0) 
                welch.append(0)  
                band.append(0)            
    
# spectral analysis of the solar radiation

    time_sol, per_sol, power_sol     = mod.wave_spectral ( ra, dt, mother)
    freq_sol, wl_aper_sol, welch_sol = mod.welch_method (ra, windsize, window, dt)

#spectral analysis of wind intensity

    time_win, per_win, power_win     = mod.wave_spectral ( iw, dt, mother)
    freq_win, wl_aper_win, welch_win = mod.welch_method (iw, windsize, window, dt)

# Coherence (isotherms and meteorological data)
#
# ------------------- for depth ----------------------------------------------
#
    if turn_temp == 1:
    
        d1 = mod.depths(seu[0], h, lin)
        d2 = mod.depths(seu[1], h, lin)
        d3 = mod.depths(seu[2], h, lin)
        d4 = mod.depths(seu[3], h, lin)
    
        depth = [d1, d2, d3, d4]

        s_filtered = []
        phws        = []
        caws        = []
        faws        = []
        c95ws       = []
    
        for i in range(4):
        
            if sen[i] == 1:
                new = mod.vector_time (temp,seu[i],lin)
                a1, a2, a3, a4 = mod.coherence_shift(new, iw, windsize, dt)
            
                phws.append(a1)
                caws.append(a2)
                faws.append(a3)
                c95ws.append(a4)
            
                s_filtered.append(sensor_filtered[seu[i]])

                if sen[0] == 1:
                    phr1, car1, far1, c951r = mod.coherence_shift(new, ra, windsize, dt)
                else:
                    phr1, car1, far1, c951r = mod.zero()
            
            
        
            else:
                phws.append(None)
                caws.append(None)
                faws.append(None)
                c95ws.append(None)
            
                s_filtered.append(None)

        phij  = []
        cohij = []
        fij   = []
        c95ij = []        

        for i in range(4):
            for j in range(4):
                if (i != j and j-i==1):                
                    if(sen[i] == 1 and sen[j] == 1):
                
                        ni = mod.vector_time (temp,seu[i],lin)
                        nj = mod.vector_time (temp,seu[j],lin)
                
                        ph, coh, f, c95 = mod.coherence_shift(ni,nj,windsize,dt)
                
                        phij.append(ph)
                        cohij.append(coh)
                        fij.append(f)
                        c95ij.append(c95)
                
                    else:
                        phij.append(None)
                        cohij.append(None)
                        fij.append(None)
                        c95ij.append(None)    

# ------------------- for isotherms -------------------------------------------
#
# the coherence with solar radiation is done just for one isotherm, 
# in which is defined by tau_int
#
    if (turn_iso == 1):

        phw, caw, faw, c95aw = mod.coherence_shift(iso[aux_iso], iw, windsize, dt)
        phr, car, far, c95ar = mod.coherence_shift(iso[aux_iso], ra, windsize, dt)   


        phiso = []
        coiso = []
        friso = []
        cliso = []

    # coherence and phase analysis between isotherms
        for i in range(4):
            for j in range(4):
                if ( i != j and j > i):
                    if(tau[i] != -999) and (tau[j] != -999):
                    
                        phaux, coaux, fiaux, c9aux = mod.coherence_shift(iso[i],iso[j],windsize,dt)
                        phiso.append(phaux)
                        coiso.append(coaux)
                        friso.append(fiaux)
                        cliso.append(c9aux)
                    else: 
                        phiso.append(None)
                        coiso.append(None)
                        friso.append(None)
                        cliso.append(None)

# ---------- contourf correction and mean temp at thermocline depth -----------
# output variables:
#   templ - temperature 
#   hl    - depth
#   riwl  - richardson number
#  


    templ, hl = graph.correction_contourf(h,temp,dj,lin,dci,qt)
    riwl, hlm = graph.correction_contourf(hzmid,riw2d,dj,lin,dci,qt-1) 

    tt, ithermo = mod.value_nearest(hl,m_ht) 

# outputs:
# tt: thermocline depth according to countourf interpolation
# ithermo: indice on 'hl' 
# -----------------------------------------------------------------------------

    for t in range(lin):   
        thermo_temp[t] = templ[t][ithermo]

# spectral analysis of the thermocline 

    tthermo, pthermo, powerthermo = mod.wave_spectral (thermo_temp, dt, mother)
    
    freqthermo, wlthermo, welchthermo = mod.welch_method(thermo_temp,windsize,window, dt)
    freq_ean, wl_aper_ean, welch_ean  = mod.welch_method(ean,windsize,window, dt)

    print (">         Spectral analysis was computed")
    root.update()
    
    
    print ("> ")
    anal_time = time.time()
    print ("> Execution time for part IV: ")
    root.update() 
    print ("> "+str(round(anal_time-load_time,4))+' seconds')
    root.update()   

    print ("--------------------------------------------------------------------------------------")
    root.update()
    print('> ')
    root.update()
    
# -----------------------------------------------------------------------------
    print ("> Part V        Generating text files... ")
    root.update()
    
    os.makedirs(output_path+'textfiles')

    if turn_iso == 1:
        for i in range(4):
            if(tau[i]!=-999):
                np.savetxt(output_path+'textfiles/spectral_isotherms'+str(tau[i])+'.txt',np.column_stack((freq[i],welch[i])), fmt='%0.8f %0.15f')                
                np.savetxt(output_path+'textfiles/isotherms'+str(tau[i])+'.txt',np.column_stack((time_temp[i],iso[i])), fmt='%0.8f %0.5f')
                
    if turn_temp == 1:
        for i in range(4):
            if(sen[i]==1):
                np.savetxt(output_path+'textfiles/spectral_sensor'+str(sen[i])+'.txt',np.column_stack((freqe[i],welch_sensor[i])), fmt='%0.8f %0.15f')
    

# -----------------------------------------------------------------------------    
    print ("> Part VI       Plotting graphs and making reports... ")
    root.update()  
    print('> ')
    root.update()
    if(depi > 400):
        print (">         Attention: Figures of high quality (DPI > 400)")
        root.update()    
        print (">         You can reduce DPI for a faster plotting")
        root.update() 
    print (">         This may take few minutes")
    root.update()


# -----------------------------------------------------------------------------
#                        Field observations - graphs   
# -----------------------------------------------------------------------------
    show_on  = 0  

    if turn_iso == 1:
        for i in range(4):
            if(tau[i]!=-999):
                plt.figure( figsize=(8,6))
        
                ax1 = plt.subplot2grid((2,2),(0,0),colspan=2)
                ax2 = plt.subplot2grid((2,2),(1,0),colspan=2)
    
                ax1.set_title('(a)',loc='left')
                ax2.set_title('(b)',loc='left')

                graph.isotherm(dx,time_temp[i],iso[i],'navy',tau[i],ax1)
                graph.wavelet(dx,per_temp[i], power[i], show_on, ax2)
                
                plt.savefig(output_path+'iso'+str(i)+'.png', dpi = depi)
    

    plt.figure(figsize=(10,7.5))
    
    ax1 = plt.subplot2grid((2,3),(0,0),colspan=2)
    ax2 = plt.subplot2grid((2,3),(1,0),colspan=2)
    ax3 = plt.subplot2grid((2,3),(1,2),sharey=ax2)

    ax1.set_title('(a) thermocline analysis',loc='left')
    ax2.set_title('(b)',loc='left')
    ax3.set_title('(c)',loc='left')

    graph.thermocline(dx, tthermo,thermo_temp,'navy',round(m_ht, 2),ax1)
    graph.wavelet(dx, pthermo, powerthermo, show_on, ax2)
    graph.psd_thermocline(wlthermo,welchthermo,'navy',int(tt),ax3)

    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.savefig(output_path+'thermocline_analysis.png', dpi = depi)



    plt.figure(figsize=(8,6))

    ax1 = plt.subplot2grid((2,1),(0,0))
    ax2 = plt.subplot2grid((2,1),(1,0))

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')

    graph.wind(dx, time_win, dw, iw, ax1)
    graph.wind_direction(dx, time_win, dw, newdw,newdw_low,newdi,  ax2)


    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.savefig(output_path+'temporal_analysis.png',dpi = depi)


    plt.figure(figsize=(8,6))

    ax1 = plt.subplot2grid((2,1),(0,0))
    ax2 = plt.subplot2grid((2,1),(1,0),sharex=ax1)

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')

    graph.windstress(dx, time_win, strs, ax1)
    graph.wedderburn(dx, time_win, wedd, ax2) 

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.savefig(output_path+'wedderburn_stress.png',dpi = depi)


    plt.figure(figsize=(5,6))
    ax1 = plt.subplot2grid((1,1),(0,0))

    graph.psd_wind(wl_aper_win, welch_win, ax1)
    ax2 = ax1.twiny() 
    graph.psd_sola(wl_aper_sol, welch_sol, ax2)

    plt.savefig(output_path+'meteo_spectra.png', dpi = depi)

    if turn_iso == 1:

        
        plt.figure(figsize=(8,6))
        
        ax1 = plt.subplot2grid((1,2),(0,0))
        ax2 = plt.subplot2grid((1,2),(0,1),sharey=ax1)
        
        ax1.set_title('(a)',loc='left')
        ax2.set_title('(b)',loc='left')

        graph.psd_isoz(tau, freq, welch, fz_11, largelen, ax1)
        
        graph.coherence (caw,car,faw,far,tau,aux_iso, largelen, ax2)

        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.savefig(output_path+'psd_nonhydro.png', dpi = depi)


        plt.figure(figsize=(8,6))

        ax1 = plt.subplot2grid((1,2),(0,0))
        ax2 = plt.subplot2grid((1,2),(0,1),sharey=ax1)

        ax1.set_title('(a) ',loc='left')
        ax2.set_title('(b) Coriolis correction',loc='left')

        graph.psd_isox(tau, freq, welch, fx_11, fx_21, fx_31,largelen, ax1)

    # this has the earth rotation effect
        graph.psd_isox(tau, freq, welch, cfx11, cfx21, cfx31,largelen, ax2)

        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.savefig(output_path+'psd_hydro_coriois.png', dpi = depi)


        plt.figure(figsize=(5,6))

        ax1 = plt.subplot2grid((1,1),(0,0))
        
        graph.psd_isox(tau, freq, welch, fx_11, fx_21, fx_31,largelen, ax1)
        
        plt.savefig(output_path+'psd_hydro.png', dpi = depi)


    plt.figure( figsize=(8,6))

    ax1 = plt.subplot2grid((3,2), (0,0), rowspan=3)
    ax2 = plt.subplot2grid((3,2), (0,1))
    ax3 = plt.subplot2grid((3,2), (1,1))
    ax4 = plt.subplot2grid((3,2), (2,1))

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b) - P1',loc='left')
    ax3.set_title('(c) - P2',loc='left')
    ax4.set_title('(d) - P3',loc='left')

    met = m_h2

    graph.degeneration(P1,P2,P3,h_lim1,h_lim2,h_lim3,W_lim1,W_lim2,W_lim3,m_pe,m_ph,met,m_he+m_hh,ls_fetch[1],ax1)

    graph.generation(genera_0,ax2)
    graph.generation(genera_1,ax3)
    graph.generation(genera_2,ax4)

    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.savefig(output_path+'degenera.png', dpi = depi)





    plt.figure( figsize=(8,6))

    ax1 = plt.subplot2grid((3,2), (0,0), rowspan=3)
    ax2 = plt.subplot2grid((3,2), (0,1), rowspan=3)


    ax1.set_title('(a) - Original Chart',loc='left')
    ax2.set_title('(b) - Zooming',loc='left')


    met = m_h2

    graph.degeneration_evolution(P1,P2,P3,hH_gp,wi_gp,m_pe,m_ph,met,m_he+m_hh,ls_fetch[1],'no',ax1)
    graph.degeneration_evolution(P1,P2,P3,hH_gp,wi_gp,m_pe,m_ph,met,m_he+m_hh,ls_fetch[1],'yes',ax2)

    plt.savefig(output_path+'degenera_evolution.png', dpi = depi)



    plt.figure( figsize=(8,6))

    ax1 = plt.subplot2grid((3,1), (0,0))
    ax2 = plt.subplot2grid((3,1), (1,0), rowspan=3, sharex=ax1)

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')

    graph.richardson(dx,time_win,riw,iw_dw, iw_up, ax1)
    graph.richardson2d(dx,time_win,hlm,riwl,ht,ean,ax2)

    plt.setp(ax1.get_xticklabels(), visible=False)
    
    plt.savefig(output_path+'richardson.png',dpi = depi)


    plt.figure(figsize=(10,7.5))

    ax1 = plt.subplot2grid((3,2),(0,0))
    ax2 = plt.subplot2grid((3,2),(0,1),sharex=ax1)
    ax3 = plt.subplot2grid((3,2),(1,0),colspan=2, rowspan=3)

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')
    ax3.set_title('(c)',loc='left')

    graph.radiation(dx,time_win,ra,ax1)
    graph.density(dx,time_win,pu,pd,pm,ax2)
    graph.tempstructure(dx,time_win,hl,ht,templ,ean,ax3)

    plt.savefig(output_path+'structure.png',dpi = depi)

    plt.figure( figsize=(10,7.5))

    ax1 = plt.subplot2grid((3,2),(0,0))
    ax2 = plt.subplot2grid((3,2),(0,1),sharex=ax1)
    ax3 = plt.subplot2grid((3,2),(1,0),colspan=2, rowspan=3)

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')
    ax3.set_title('(c)',loc='left')

    graph.radiation(dx,time_win,ra,ax1)
    graph.density(dx,time_win,pu,pd,pm,ax2)
    graph.tempstructure_zoom(dx,time_win,hl,ht,templ,ean,ax3)
    plt.savefig(output_path+'structure2.png',dpi = depi)


    plt.figure( figsize=(8,8))

    ax1 = plt.subplot2grid((2,2),(0,0))
    ax2 = plt.subplot2grid((2,2),(0,1))
    ax3 = plt.subplot2grid((2,2),(1,0))
    ax4 = plt.subplot2grid((2,2),(1,1))

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')
    ax3.set_title('(c)',loc='left')
    ax4.set_title('(d)',loc='left')

    graph.depth_sensitivity(xhe,v1h1_she,r'$h_e$',0,ax1)
    graph.densi_sensitivity(xpe,v1h1_spe,r'$\rho_e$',0,ax2)
    graph.depth_sensitivity(xhh,v2h1_sh2,r'$h_h$',1,ax3)
    graph.densi_sensitivity(xph,v2h1_sp2,r'$\rho_h$',1,ax4)

    plt.savefig(output_path+'sensitivity_v1h1.png', dpi = depi)


    plt.figure(figsize=(8,8))

    ax1 = plt.subplot2grid((2,2),(0,0))
    ax2 = plt.subplot2grid((2,2),(0,1))
    ax3 = plt.subplot2grid((2,2),(1,0))
    ax4 = plt.subplot2grid((2,2),(1,1))

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')
    ax3.set_title('(c)',loc='left')
    ax4.set_title('(d)',loc='left')

    graph.depth_sensitivity(xh2,v2h1_sh2,r'$h_2$',0,ax1)
    graph.densi_sensitivity(xp2,v2h1_sp2,r'$\rho_2$',0,ax2)
    graph.depth_sensitivity(xh3,v3h1_sh3,r'$h_3$',1,ax3)
    graph.densi_sensitivity(xp3,v3h1_sp3,r'$\rho_3$',1,ax4)

    plt.savefig(output_path+'sensitivity_multiv.png', dpi = depi)

    
    plt.figure(figsize=(8,3))

    ax1 = plt.subplot2grid((1,1),(0,0))
    graph.ri_compb(dx,time_win, riw, ridi, iw_dw, iw_up, ax1)

    plt.savefig(output_path+'ri_filtering.png',dpi = depi)

    if turn_iso == 1:
        plt.figure(figsize=(8,3))
        
        ax1 = plt.subplot2grid((1,1),(0,0))
        graph.multi_isotherms(dx, iso, tau, time_temp, ax1)
        
        plt.savefig(output_path+'isotherms.png',dpi = depi)


        if(ana1 != -999 and ana2 != -999):
            fig, ax1 = plt.subplots(figsize=(8,5))
            graph.coherence_iso(tau,coiso,friso,ana1,ana2, ax1)   
        
            ax2 = ax1.twinx()
        
            graph.phase_iso(tau,phiso,friso,cliso,ana1,ana2, ax2)   
        

            plt.savefig(output_path+'coherence.png',dpi = depi)



    plt.figure(figsize=(8,3))
    ax1 = plt.subplot2grid((1,1),(0,0))

    graph.wavelet_resona(dx,per_win, power_win, v1mode, v2mode, ax1)
    plt.savefig(output_path+'wind_resonance.png',dpi = depi)


    if turn_iso == 1:

        
        plt.figure(figsize=(8,3))
        ax1 = plt.subplot2grid((1,1),(0,0))
        if filter_process == 1:
            ax1.set_title('Band-pass filter - Bandwidth: '+str(round(pxv1h1[2]/60/60,1))+' h to '+str(round(pxv1h1[0]/60/60,1))+' h',loc='left')
        
        elif filter_process == 2:
            ax1.set_title('Band-pass filter - Bandwidth: '+str(round(low_per,1))+' h to '+str(round(high_per,1))+' h',loc='left')
    
        graph.temp_bandpass(dx, band, tau, ax1)


        
        plt.setp(ax4.get_xticklabels(), visible=False)
        plt.savefig(output_path+'iso_bandpass.png',dpi = depi)


    fig=plt.figure()
    graph.wind_rose(iw,dw)

    plt.savefig(output_path+'windrose.png',dpi = depi)


    if turn_temp == 1:

        for jd in range(4):
        
            if seu[jd] != -999:
                plt.figure( figsize=(8,6))
                ax1 = plt.subplot2grid((2,2),(0,0),colspan=2)
                ax2 = plt.subplot2grid((2,2),(1,0),colspan=2,sharex=ax1)
    
                ax1.set_title('(a)',loc='left')
                ax2.set_title('(b)',loc='left')
    
                ids = seu[jd]
                dep = mod.vector_time (temp,ids,lin)
    
                graph.isodepth(dx, timee[ids],dep,'navy',depth[jd],ax1)
                graph.wavelet_depth(dx, per[ids], power_sensor[ids], timee[ids], ax2)

                plt.savefig(output_path+'depth'+str(int(jd))+'.png', dpi =depi)
    
    
        plt.figure( figsize=(8,6))

        ax1 = plt.subplot2grid((1,2),(0,0))
        ax2 = plt.subplot2grid((1,2),(0,1),sharey=ax1)

        ax1.set_title('(a)',loc='left')
        ax2.set_title('(b)',loc='left')

        graph.psd_depthz(sen,depth, freqe, welch_sensor, fz_11, largelen, ax1)

        graph.coh_depthwind(caws,car1,faws,far1,depth, largelen,sen, ax2)

            
        plt.setp(ax2.get_yticklabels(), visible=False)

        plt.savefig(output_path+'psd_nonhydro_depth.png', dpi = depi)


        plt.figure( figsize=(8,6))

        ax1 = plt.subplot2grid((1,2),(0,0))
        ax2 = plt.subplot2grid((1,2),(0,1),sharey=ax1)

        ax1.set_title('(a)',loc='left')
        ax2.set_title('(b)',loc='left')
    
        graph.psd_isox_depth(sen,depth, freqe, welch_sensor, fx_11, fx_21, fx_31,largelen, ax1)
        graph.psd_isox_depth(sen,depth, freqe, welch_sensor, cfx11, cfx21, cfx31,largelen, ax2)

        plt.setp(ax2.get_yticklabels(), visible=False)

        plt.savefig(output_path+'psd_coriolis_depth.png', dpi = depi)


        plt.figure( figsize=(5,6))
        ax1 = plt.subplot2grid((1,1),(0,0))

        graph.psd_isox_depth(sen,depth, freqe, welch_sensor, fx_11, fx_21, fx_31,largelen, ax1)


        plt.savefig(output_path+'psd_hydro_depth.png', dpi = depi)


        fig, ax1 = plt.subplots(figsize=(8,5))
        graph.coherence_depth(depth,cohij,fij,ax1)
    
        ax2 = ax1.twinx()
        graph.phase_depth(depth, c95ij, phij,fij,ax2)

        plt.savefig(output_path+'phase_depth.png',dpi = depi)


        plt.figure( figsize=(8,3))
        ax1 = plt.subplot2grid((1,1),(0,0))
        
        graph.thermal_variation(dx, depth,seu,temp,time_win, 0, ax1)

        plt.savefig(output_path+'temperature_depth.png',dpi = depi)
    
    
        plt.figure(figsize=(8,3))

        ax1 = plt.subplot2grid((1,1),(0,0))

        if filter_process == 1:
            ax1.set_title('Band-pass filter - Bandwidth: '+str(round(pxv1h1[2]/60/60,1))+' h to '+str(round(pxv1h1[0]/60/60,1))+' h',loc='left')
        
        elif filter_process == 2:
            ax1.set_title('Band-pass filter - Bandwidth: '+str(round(low_per,1))+' h to '+str(round(high_per,1))+' h',loc='left')
        
        graph.depth_bandpass(dx, s_filtered,depth, time_win,sen, ax1)
        plt.savefig(output_path+'depth_bandpass.png',dpi = depi)

  
    
    print (">         Making the Interwave Analyzer's report")
    root.update()
# -----------------------------------------------------------------------------
#  reportlab 
# -----------------------------------------------------------------------------

    ini = date[0].split('/') 
    fin = date[t].split('/')

    canvas = canvas.Canvas(output_path+'form.pdf', pagesize=A4)
    canvas.setLineWidth(.3)
    canvas.setFont('Helvetica',10)
    if (turn_temp == 1 and turn_iso == 1):
            canvas.drawString(520,800,'Page 1/5')
    else:
        if turn_temp == 1 or turn_iso == 1:
            canvas.drawString(520,800,'Page 1/4')
        else:
            canvas.drawString(520,800,'Page 1/3')


    logo =ImageReader('0interwave.png')
    canvas.drawImage(logo, 30, 733,width=120, height=70)

    canvas.drawString(250, 760, 'OFFICIAL REPORT')

    canvas.drawString(30, 703, 'Analyzed Period:  ') 
    canvas.drawString(120, 703, str(ini[2]) +' / '+str(ini[1])+' / '+str(ini[0])+' at '+str(ini[3])+'h  to '+ str(fin[2]) +' / '+str(fin[1])+' / '+str(fin[0])+' at '+str(fin[3])+'h')
    canvas.line(115,700,550,700)

    canvas.drawString(30, 660, 'Basin Length (wind direction:   '+str(round(dw_mean,0))+'° \xb1 '+str(round(linang,0))+'°)')
    canvas.drawString(30, 645, 'min: '+str(round(ls_fetch[0],0))+' m; max: '+str(round(ls_fetch[2],0))+' m; ave: '     +str(round(ls_fetch[1],0))+' m')
    canvas.drawString(30, 625, 'Layers Thickness (Total depth: '+str(round(m_he+m_hh,2))+' m):')
    canvas.drawString(30, 610, 'Two-layer')
    canvas.drawString(110, 610, 'Epilimnion')
    canvas.drawString(180, 610, str(round(m_he,2))+' m')
    canvas.drawString(110, 595, 'Hypolimnion')
    canvas.drawString(180, 595, str(round(m_hh,2))+' m')
    canvas.drawString(30, 580, 'Three-layer')
    canvas.drawString(110, 580, 'Epilimnion')
    canvas.drawString(180, 580, str(round(m_h1,2))+' m')
    canvas.drawString(110, 565, 'Metalimnion')
    canvas.drawString(180, 565, str(round(m_h2,2))+' m')
    canvas.drawString(110, 550, 'Hypolimnion')     
    canvas.drawString(180, 550, str(round(m_h3,2))+' m')              
    canvas.drawString(30, 530, 'Layers Density and stratification:')
    canvas.drawString(30, 515, 'Two-layer ')
    canvas.drawString(110, 515, 'Epilimnion')
    canvas.drawString(180, 515, str(round(m_pe,2))+' kg/m³')
    canvas.drawString(110, 500, 'Hypolimnion')
    canvas.drawString(180, 500, str(round(m_ph,2))+' kg/m³')
    canvas.drawString(30, 485, 'Three-layer')
    canvas.drawString(110, 485, 'Epilimnion')
    canvas.drawString(180, 485, str(round(m_p1,2))+' kg/m³')
    canvas.drawString(110, 470, 'Metalimnion')
    canvas.drawString(180, 470, str(round(m_p2,2))+' kg/m³')
    canvas.drawString(110, 455, 'Hypolimnion')
    canvas.drawString(180, 455,  str(round(m_p3,2))+' kg/m³')

    g1 =ImageReader(output_path+'windrose.png')
    canvas.drawImage(g1, 25, 190,width=250, height=250)

    g2 =ImageReader(output_path+'structure.png')
    canvas.drawImage(g2, 250, 450,width=300, height=225)


    canvas.drawString(280, 420, 'Wind parameters: ')
    canvas.drawString(280, 400, 'Duration of the strongest wind event:')
    canvas.drawString(460, 400,  str(round(dura/(60*60),2))+' h')
    canvas.drawString(280, 380, 'Just considering homogeneous direction:')
    canvas.drawString(280, 365,  str(round(dire/(60*60),2))+' h blowing '+str(round(m_newind,0))+'° North')
    canvas.drawString(280, 345, 'Reduction factor:')
    canvas.drawString(380, 345, 'Duration factor:')
    canvas.drawString(460, 345,  str(round(fdura,3)))
    canvas.drawString(380, 330, 'Direction factor:')
    canvas.drawString(460, 330,  str(round(fdire,3)))
    canvas.drawString(280, 310, 'Mean friction velocity of the wind:')
    canvas.drawString(440, 310,  str(round(m_wast,5))+' m/s')
    canvas.drawString(280, 295, 'Minimum friction velocity of the wind:')
    canvas.drawString(440, 295,  str(round(np.nanmin(wast),5))+' m/s')
    canvas.drawString(280, 280, 'Maximum friction velocity of the wind:')
    canvas.drawString(440, 280,  str(round(np.nanmax(wast),5))+' m/s')


    g3 =ImageReader(output_path+'temporal_analysis.png')
    canvas.drawImage(g3, 250, 50,width=300, height=225)

    canvas.drawString(30, 170, 'Paramters of Stability:')
    canvas.drawString(30, 150, 'Reduced Gravity:')

    canvas.drawString(150, 150,  str(round(m_glin,4))+' \xb1 '+str(round(mod.ciout(glin),4))+' m/s²')
    canvas.drawString(30, 135, 'Brunt-Vaisalla:')
    canvas.drawString(150, 135,  str(round(m_n,4))+' \xb1 '+str(round(mod.ciout(n),4))+' Hz')
    canvas.drawString(30, 120, 'Averaged Richardson number:')
    canvas.drawString(180, 120,  str(round(m_riw,0)))

    canvas.drawString(30, 90, 'Daily averaged variables at thermocline depth:')
    canvas.drawString(30, 70,  'Richardson number:')
    canvas.drawString(140, 70,  str(round(np.nanmean(mod.average(riw,24/dt)),0))+' \xb1 '+str(round(np.std(mod.average(riw,24/dt))*1.96/np.sqrt(dt*lin/24),0))) # 1.96 = 95%

    canvas.drawString(30, 55, 'Wedderburn number:')
    canvas.drawString(140,55,  str(round(np.nanmean(mod.average(wedd,24/dt)),0))+' \xb1 '+str(round(np.std(mod.average(wedd,24/dt))*1.96/np.sqrt(dt*lin/24),0)))

#canvas.drawString(40, 100, 'Isoterms analyzed (zero value is a empty iso)')
#canvas.drawString(40, 85, 'Isotherms: '+str(taua)+' °C; '+str(taub)+' °C; '+str(tauc)+' °C; '+str(taud)+' °C')

    canvas.showPage()

    canvas.setLineWidth(.3)
    canvas.setFont('Helvetica',10)
    if (turn_temp == 1 and turn_iso == 1):
        canvas.drawString(520,800,'Page 2/5')
    else:
        if turn_temp == 1 or turn_iso == 1:
            canvas.drawString(520,800,'Page 2/4')
        else:
            canvas.drawString(520,800,'Page 2/3')


    g7 =ImageReader(output_path+'meteo_spectra.png')
    canvas.drawImage(g7, 320, 500,width=250, height=300)


    canvas.drawString(30,780,'Daily averaged Richardson number at thermocline depth:')
    canvas.drawString(30,760,'Richardson number: ')
    canvas.drawString(200,760, str(round(np.nanmean(mod.average(riw,24/dt)),0))+' \xb1 '+str(round(np.std(mod.average(riw,24/dt))*1.96/np.sqrt(dt*lin/24),0)))
    canvas.drawString(30,745,'Filtered by wind duration:')
    canvas.drawString(200,745, str(round(np.nanmean(mod.average(ridu,24/dt)),0)))
    canvas.drawString(30,730,'Filtered by wind direction:')
    canvas.drawString(200,730, str(round(np.nanmean(mod.average(ridi,24/dt)),0)))

    g13 =ImageReader(output_path+'richardson.png')
    canvas.drawImage(g13, 20, 490,width=300, height=225)


    g13 =ImageReader(output_path+'degenera.png')
    canvas.drawImage(g13, 20, 245,width=300, height=225)

    canvas.drawString(340,480,'Generation & Degeneration Theory¹:')

    canvas.drawString(340,460,'Periods')
    canvas.drawString(400,460,'1/W')
    canvas.drawString(460,460,'he/H:')

    canvas.drawString(340,440,'P1')
    canvas.drawString(400,440, str(round(np.nanmax(wi_gp[0]),5)))
    canvas.drawString(460,440, str(round(np.nanmean(hH_gp[0]),5)))

    canvas.drawString(340,425,'P2')
    canvas.drawString(400,425, str(round(np.nanmax(wi_gp[1]),5)))
    canvas.drawString(460,425, str(round(np.nanmean(hH_gp[1]),5)))

    canvas.drawString(340,410,'P3')
    canvas.drawString(400,410, str(round(np.nanmax(wi_gp[2]),5)))
    canvas.drawString(460,410, str(round(np.nanmean(hH_gp[2]),5)))

    canvas.drawString(340,390,'Maximum amplitude of BSIW:')
    canvas.drawString(340,375,'P1')
    canvas.drawString(460,375, str(round(np.nanmax(wi_gp[0])*m_he/2,3))+' m')
    canvas.drawString(340,360,'P2')
    canvas.drawString(460,360, str(round(np.nanmax(wi_gp[1])*m_he/2,3))+' m')
    canvas.drawString(340,345,'P3')
    canvas.drawString(460,345, str(round(np.nanmax(wi_gp[2])*m_he/2,3))+' m')

    canvas.drawString(340,325,'Internalwave Classification:')
    canvas.drawString(340,310,'P1')
    canvas.drawString(400,310, mod.subclass(np.nanmax(wi_gp[0])))
    canvas.drawString(340,295,'P2')
    canvas.drawString(400,295, mod.subclass(np.nanmax(wi_gp[1])))
    canvas.drawString(340,280,'P3')
    canvas.drawString(400,280, mod.subclass(np.nanmax(wi_gp[2])))

    canvas.drawString(340,260,'¹ Strongest BSIW that should be detected')
    canvas.drawString(150,230,'Probable amplitude of BSIW according to Wedderburn number:')


    canvas.drawString(30,210,'Wedderburn < 100')
    canvas.drawString(30,195,'Periods')
    canvas.drawString(100,195,'Amplitude')
    canvas.drawString(180,195,'Duration Ratio²')

    canvas.drawString(30,180,'P1')
    wig_ave, ratio = mod.amplitude_average(wi_gp[0],m_he,100)
    canvas.drawString(100,180, str(round(wig_ave,3))+' m')
    canvas.drawString(180,180, str(round(ratio,4)))

    canvas.drawString(30,165,'P2')
    wig_ave, ratio = mod.amplitude_average(wi_gp[1],m_he,100)
    canvas.drawString(100,165, str(round(wig_ave,3))+' m')
    canvas.drawString(180,165, str(round(ratio,4)))

    canvas.drawString(30,150,'P3')
    wig_ave, ratio = mod.amplitude_average(wi_gp[2],m_he,100)
    canvas.drawString(100,150, str(round(wig_ave,3))+' m')
    canvas.drawString(180,150, str(round(ratio,4)))

    canvas.drawString(300,210,'Wedderburn < 50')
    canvas.drawString(300,195,'Periods')
    canvas.drawString(370,195,'Amplitude')
    canvas.drawString(450,195,'Duration Ratio²')

    canvas.drawString(300,180,'P1')
    wig_ave, ratio = mod.amplitude_average(wi_gp[0],m_he,50)
    canvas.drawString(370,180, str(round(wig_ave,3))+' m')
    canvas.drawString(450,180, str(round(ratio,4)))

    canvas.drawString(300,165,'P2')
    wig_ave, ratio = mod.amplitude_average(wi_gp[1],m_he,50)
    canvas.drawString(370,165, str(round(wig_ave,3))+' m')
    canvas.drawString(450,165, str(round(ratio,4)))

    canvas.drawString(300,150,'P3')
    wig_ave, ratio = mod.amplitude_average(wi_gp[2],m_he,50)
    canvas.drawString(370,150, str(round(wig_ave,3))+' m')
    canvas.drawString(450,150, str(round(ratio,4)))


    canvas.drawString(30,130,'Wedderburn < 20')
    canvas.drawString(30,110,'Periods')
    canvas.drawString(100,110,'Amplitude')
    canvas.drawString(180,110,'Duration Ratio²')

    canvas.drawString(30,95,'P1')
    wig_ave, ratio = mod.amplitude_average(wi_gp[0],m_he,20)
    canvas.drawString(100,95, str(round(wig_ave,3))+' m')
    canvas.drawString(180,95, str(round(ratio,4)))

    canvas.drawString(30,80,'P2')
    wig_ave, ratio = mod.amplitude_average(wi_gp[1],m_he,20)
    canvas.drawString(100,80, str(round(wig_ave,3))+' m')
    canvas.drawString(180,80, str(round(ratio,4)))

    canvas.drawString(30,65,'P3')
    wig_ave, ratio = mod.amplitude_average(wi_gp[2],m_he,20)
    canvas.drawString(100,65, str(round(wig_ave,3))+' m')
    canvas.drawString(180,65, str(round(ratio,4)))

    canvas.drawString(300,130,'Wedderburn < 3')
    canvas.drawString(300,110,'Periods')
    canvas.drawString(370,110,'Amplitude')
    canvas.drawString(450,110,'Duration Ratio²')

    canvas.drawString(300,95,'P1')
    wig_ave, ratio = mod.amplitude_average(wi_gp[0],m_he,3)
    canvas.drawString(370,95, str(round(wig_ave,3))+' m')
    canvas.drawString(450,95, str(round(ratio,4)))

    canvas.drawString(300,80,'P2')
    wig_ave, ratio = mod.amplitude_average(wi_gp[1],m_he,3)
    canvas.drawString(370,80, str(round(wig_ave,3))+' m')
    canvas.drawString(450,80, str(round(ratio,4)))

    canvas.drawString(300,65,'P3')
    wig_ave, ratio = mod.amplitude_average(wi_gp[2],m_he,3)
    canvas.drawString(370,65, str(round(wig_ave,3))+' m')
    canvas.drawString(450,65, str(round(ratio,4)))

    canvas.drawString(30,40,'² Ratio between duration period of Wedderburn < than the criteria and the total period')



    canvas.showPage()

    canvas.setLineWidth(.3)
    canvas.setFont('Helvetica',10)
    if (turn_temp == 1 and turn_iso == 1):
        canvas.drawString(520,800,'Page 3/5')
    else:
        if turn_temp == 1 or turn_iso == 1:
            canvas.drawString(520,800,'Page 3/4')
        else:
            canvas.drawString(520,800,'Page 3/3')

    canvas.drawString(50,750,'Analysis of the Thermocline and Theoretical results')


    g4 =ImageReader(output_path+'thermocline_analysis.png')
    canvas.drawImage(g4, 0, 380,width=550, height=412.5)
    canvas.drawString(50,390,'¹ Wind criteria to excite internal seiche (wind durastion/internal wave period > 0.25)')

    canvas.drawString(370, 740, 'Internal Seiche Modeling:')

    canvas.drawString(370, 720, 'Merian Equation:')
    canvas.drawString(370, 705, 'Wave speed:')
    canvas.drawString(450, 705,  str(round(cp,4))+' m/s')
    canvas.drawString(370, 690, 'Wave period:')
    canvas.drawString(450, 690,  str(round(tbsiw/60/60,2))+' hours')

    canvas.drawString(370,670, 'Wind duration criteria¹:') #(Tdur/Tbsiw>0.25)
    canvas.drawString(370,650, 'wind duration:')
    canvas.drawString(450,650, str(round(dura/tbsiw,3)))
    canvas.drawString(370,635, 'wind direction:')
    canvas.drawString(450,635, str(round(dire/tbsiw,3)))
    
    canvas.drawString(100, 370, 'Hydrostatic Model for the first three vertical and horizontal modes:')

    canvas.drawString(50, 350, 'V1H1' )
    canvas.drawString(90, 350,  str(round(pxv1h1[1]/(60*60),2))+' h')
    canvas.drawString(140, 350, ' \xb1 ')
    canvas.drawString(160, 350, str(round((pxv1h1[1]-pxv1h1[0])/(60*60),2))+' h')

    canvas.drawString(50, 335, 'V2H1' )
    canvas.drawString(90, 335,  str(round(pxv2h1[1]/(60*60),2))+' h')
    canvas.drawString(140, 335, ' \xb1 ')
    canvas.drawString(160, 335, str(round((pxv2h1[1]-pxv2h1[0])/(60*60),2))+' h')

    canvas.drawString(50, 320, 'V3H1' )
    canvas.drawString(90, 320,  str(round(pxv3h1[1]/(60*60),2))+' h')
    canvas.drawString(140, 320, ' \xb1 ')
    canvas.drawString(160, 320, str(round((pxv3h1[1]-pxv3h1[0])/(60*60),2))+' h')


    canvas.drawString(220, 350, 'V1H2' )
    canvas.drawString(260, 350,  str(round(pxv1h2[1]/(60*60),2))+' h')
    canvas.drawString(310, 350, ' \xb1 ')
    canvas.drawString(330, 350, str(round((pxv1h2[1]-pxv1h2[0])/(60*60),2))+' h')

    canvas.drawString(220, 335, 'V2H2' )
    canvas.drawString(260, 335,  str(round(pxv2h2[1]/(60*60),2))+' h')
    canvas.drawString(310, 335, ' \xb1 ')
    canvas.drawString(330, 335, str(round((pxv2h2[1]-pxv2h2[0])/(60*60),2))+' h')

    canvas.drawString(220, 320, 'V3H2' )
    canvas.drawString(260, 320,  str(round(pxv3h2[1]/(60*60),2))+' h')
    canvas.drawString(310, 320, ' \xb1 ')
    canvas.drawString(330, 320, str(round((pxv3h2[1]-pxv3h2[0])/(60*60),2))+' h')


    canvas.drawString(390, 350, 'V1H3' )
    canvas.drawString(430, 350,  str(round(pxv1h3[1]/(60*60),2))+' h')
    canvas.drawString(470, 350, ' \xb1 ')
    canvas.drawString(490, 350, str(round((pxv1h3[1]-pxv1h3[0])/(60*60),2))+' h')

    canvas.drawString(390, 335, 'V2H3' )
    canvas.drawString(430, 335,  str(round(pxv2h3[1]/(60*60),2))+' h')
    canvas.drawString(470, 335, ' \xb1 ')
    canvas.drawString(490, 335, str(round((pxv2h3[1]-pxv2h3[0])/(60*60),2))+' h')

    canvas.drawString(390, 320, 'V3H3' )
    canvas.drawString(430, 320,  str(round(pxv3h3[1]/(60*60),2))+' h')
    canvas.drawString(470, 320, ' \xb1 ')
    canvas.drawString(490, 320, str(round((pxv3h3[1]-pxv3h3[0])/(60*60),2))+' h')



    g5 =ImageReader(output_path+'sensitivity_v1h1.png')
    canvas.drawImage(g5, 20, 20,width=270, height=270)

    g6 =ImageReader(output_path+'sensitivity_multiv.png')
    canvas.drawImage(g6, 280, 20,width=270, height=270)

    canvas.drawString(160, 290, 'Sensibility Analysis for the Two-layers Internal Wave Model:')

    if (turn_iso == 1):
        canvas.showPage()
        canvas.setLineWidth(.3)
        canvas.setFont('Helvetica',10)
        if (turn_temp == 1):
            canvas.drawString(520,800,'Page 4/5')
        else:
            canvas.drawString(520,800,'Page 4/4')

        ghh =ImageReader(output_path+'psd_hydro.png')
        canvas.drawImage(ghh, 40, 270,width=280, height=336)

        gnh =ImageReader(output_path+'isotherms.png')
        canvas.drawImage(gnh, 130, 631,width=450, height=168.75)


        canvas.drawString(350, 590, 'Normalized isotherms fluctuation¹:')
        if (tau[0]!=-999):
            canvas.drawString(350, 570, str(tau[0])+' °C') 
            canvas.drawString(410, 570, str(round(iso_corrected[0],2))+' m' ) 
        if (tau[1]!=-999):
            canvas.drawString(350, 555, str(tau[1])+' °C') 
            canvas.drawString(410, 555, str(round(iso_corrected[1],2))+' m' ) 
        if (tau[2]!=-999):
            canvas.drawString(350, 540, str(tau[2])+' °C') 
            canvas.drawString(410, 540, str(round(iso_corrected[2],2))+' m' ) 
        if (tau[3]!=-999):
            canvas.drawString(350, 525, str(tau[3])+' °C') 
            canvas.drawString(410, 525, str(round(iso_corrected[3],2))+' m' ) 
    
 
        if(tau[0]!=-999):
            g9 =ImageReader(output_path+'iso0.png')
            canvas.drawImage(g9, 300, 220,width=250, height=230)

        if(tau[1]!=-999):
            g10 =ImageReader(output_path+'iso1.png')
            canvas.drawImage(g10, 30, 0,width=250, height=230)
    
        if(tau[2]!=-999):
            g11 =ImageReader(output_path+'iso2.png')
            if(tau[1]==-999):
                canvas.drawImage(g11, 300, 0,width=250, height=230)
            else:
                canvas.drawImage(g11, 300, 220,width=250, height=230)
    
        canvas.drawString(350, 505, '¹ normalized to neglect the contribution') 
        canvas.drawString(350, 490, 'of the water surface') 
  
    if (turn_temp == 1):
        canvas.showPage()
        canvas.setLineWidth(.3)
        canvas.setFont('Helvetica',10)
        if (turn_iso == 1):
            canvas.drawString(520,800,'Page 5/5')
        else:
            canvas.drawString(520,800,'Page 4/4')

        ghh =ImageReader(output_path+'psd_hydro_depth.png')
        canvas.drawImage(ghh, 40, 270,width=280, height=336)
    
        gnh =ImageReader(output_path+'temperature_depth.png')
        canvas.drawImage(gnh, 130, 631,width=450, height=168.75)

 
        if(sen[0]!=0):
            g9 =ImageReader(output_path+'depth0.png')
            canvas.drawImage(g9, 300, 370,width=250, height=230)

        if(sen[1]!=0):
            g10 =ImageReader(output_path+'depth1.png')
            canvas.drawImage(g10, 30, 0,width=250, height=230)
    
        if(sen[2]!=0):
            g11 =ImageReader(output_path+'depth2.png')
            if(sen[0]!=0):
                canvas.drawImage(g11, 300, 150,width=250, height=230)
            else:
                canvas.drawImage(g11, 300, 370,width=250, height=230)


  
    canvas.showPage()
    canvas.save()


    print ("> ")
    root.update()
    plot_time = time.time()
    print ("> Execution time for part V and VI: ")
    root.update() 
    print ("> "+str(round(plot_time-anal_time,4))+' seconds')
    root.update()  
    print ("> ")
    root.update()
    print ("> ")
    root.update()

    print ("--------------------------------------------------------------------------------------")
    root.update()
    print ("> ")
    root.update()

    print ("> FINISHED            Interwave Analyzer ")
    root.update() 
    print ("> Check the following path for results:")
    root.update()
    print ("> "+output_path)
    root.update() 
    print ("> ")
    root.update()
    print ("> ")
    root.update()
    print ("> For additional information:")
    root.update()
    print ("> https://sites.google.com/view/rafaelbueno/interwave-analyzer")
    root.update()
    print ("> ")
    root.update()
    print ("> ")
    root.update()
    root.update()
    
    root.mainloop()
    sys.stdout = old_stdout
