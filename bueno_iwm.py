# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer main module (heart module)

modififications to previosly versions must be specified here using WAVE codes:
    
W-23.02-1.00.3-00
A-01.02-1.00.3-01
V-22.02-1.00.3-01
E-05.02-1.00.3-01
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
import info_warning as war
#
# From Python Library
#
import os 
import sys
import time
import datetime
import numpy as np

from tkinter import *
from matplotlib import pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader

import matplotlib
matplotlib.use('Agg')

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
    root.geometry('800x800')
    
    root.iconbitmap("interwave_icon.ico")

    outputPanel = Text(root, wrap='word', height=30, width=100)
    outputPanel.grid(column=0, row=0, columnspan = 2, sticky='NSWE', padx=5, pady=5)

    sys.stdout = StdoutRedirector(outputPanel)
    

    print ("> Interwave Analyzer is starting the data processing... ")
    root.update()
    print ("--------------------------------------------------------------------------------------")
    root.update()
    print ("> Interwave Analyzer, version 1.00.2        July   2020")
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
# READ INFORMATION FROM GUI
# -----------------------------------------------------------------------------

    with open('temporary.txt') as reader:

        nam = reader.readline()
        win = reader.readline()
        sen_nam = reader.readline()
        sen_nam = sen_nam.replace('\n','') 
        rw  = float(reader.readline())
        rad = int(reader.readline())
        z0  = float(reader.readline())
        linang = float(reader.readline())
        lat = float(reader.readline())
        depi = int(reader.readline())
        
        
        type_length = float(reader.readline())
        if type_length == 1:
            len_basin = float(reader.readline())
        elif type_length == 2:
            fna = reader.readline()

        ean_serie = float(reader.readline())
        if ean_serie == 1:
            ean_cota = float(reader.readline())
            nac      = -999
        elif ean_serie == 2:
            nac = reader.readline() 
            ean_cota = -999
           
        minval = float(reader.readline()) 
        
        filter_process = int(reader.readline())
        if filter_process == 2:
            low_per  = float(reader.readline())
            high_per = float(reader.readline())
        else:
            for i in range(2):
                next(reader)

        autoff_wsuser = int(reader.readline())
        windsizeuser  = float(reader.readline())
        windsizeuser  = windsizeuser*24*60*60         
              
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

        turn_temp = int(reader.readline())  
    
        sen = np.zeros(4,float)
        seu = np.zeros(4,int)
        
        for i in range(4):
            seu[i] = int(reader.readline())
            if turn_temp == 0 or seu[i]==-999:
                seu[i] = -999
                sen[i] = 0
            else:
                sen[i] = 1
 
        largelen = int(reader.readline()) 
        output_path = reader.readline()   
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
    ean, h, tempa           = load.serial_cota(serialt,nac,lin,qt,temp,sen_nam,ean_serie,ean_cota,z0)
    wind, dw, ra            = load.serial_wind(serialt,win,rad,lin)


    if type_length == 1: 
        angle, dists = None, len_basin 
        
    elif type_length == 2:
        angle, dists  = load.windfetch(dw,fna)

    print ("> Part III     Defining constants and creating variables... ")
    root.update()  
# -----------------------------------------------------------------------------
# DEFINING PARAMETER AND CONSTANTS + CREATING THE DIAGNOSIS FILE 
# -----------------------------------------------------------------------------

    aux_iso = 1       # iso number to apply the corehence analysis with wind

    dj       = 0.1    # spatial resolution for countourf graph (m) 
    dci      = 1      # decimal of countouf graph (generally is refered to dj)   

    ddep   = 0.5      # range used to define the model sensitivity
    drho   = 0.1

    lh_sen = 100      # 0.01 for +- 0.5m from the reference thickness 
    lp_sen = 20       # 0.01 for +- drho kg/m³ from the reference density
    
    dig = open(output_path+'diagnosis.txt', 'w')
    dig.write("-------------------------------------------------------------------------------------\n")
    dig.write("Interwave Analyzer diagnosis\n")
    dig.write("-------------------------------------------------------------------------------------\n\n")
    
# -----------------------------------------------------------------------------
# DEFINING ARRAYS 
# -----------------------------------------------------------------------------    

# general parameters:
    glin = np.zeros((lin),float)     # reduced gravity (m/s²)
    wast = np.zeros((lin),float)     # friction velocity of the wind  (m/s) 
    strs = np.zeros((lin),float)     # Winds stress (N/m²)
    
    n   = np.zeros((lin),float)     # stability frequency (rad*Hz)
    riw = np.zeros((lin),float)     # wind richardson number (-)                                        
    wedd     = np.zeros((lin),float) # Wedderburn Number (-) 
    wedd_inv = np.zeros((lin),float) # Inverse of the Wedderburn Number (dn/dh)

    thermo_temp = np.zeros((lin),float) # temperature of the thermocline (°C)

# two-layer system:
    ht = np.zeros((lin),float)     # thermocline depth (m) 
    he = np.zeros((lin),float)     # epilimnion thickness (m)
    hh = np.zeros((lin),float)     # hypolimnion thickness (m)     
    pe = np.zeros((lin),float)     # average density of the epilimnion (kg/m³)
    ph = np.zeros((lin),float)     # average density of the hypolimnion (kg/m³)

    Pbar = np.zeros((lin),float)
    Hbar = np.zeros((lin),float)

    v1mode   = np.zeros((lin),float)    # model sensitivity V1
    v2mode   = np.zeros((lin),float)    # model sensitivity V2

# three-layer system:
    h1 = np.zeros((lin),float)     # 1st layer thickness (m)
    h2 = np.zeros((lin),float)     # 2nd layer thickness (m)
    h3 = np.zeros((lin),float)     # 3rd layer thickness (m)
    p1 = np.zeros((lin),float)     # average density of 1st layer (kg/m³)
    p2 = np.zeros((lin),float)     # average density of 2nd layer (kg/m³)
    p3 = np.zeros((lin),float)     # average density of 3rd layer (kg/m³)

    pu = np.zeros((lin),float)     # "surface" density (1m above) (kg/m³)
    pd = np.zeros((lin),float)     # "bottom" density (2m from btm) (kg/m³)
    hH = np.zeros((lin),float)     # ratio of he/H (-).

# isotherms:
    isoa = np.zeros((lin),float)     # isotherm 'a' (m) 
    isob = np.zeros((lin),float)     # isotherm 'b' (m) 
    isoc = np.zeros((lin),float)     # isotherm 'c' (m) 
    isod = np.zeros((lin),float)     # isotherm 'd' (m) 

# classification:
    genera = np.zeros((lin),float)   # class. of internal seiche generation
    iw_up = np.zeros((lin),float)    # upper limit for BSIW formation 
    iw_dw = np.zeros((lin),float)    # lower limit for BSIW formation
    
# two-dim parameters (vec[time, depth]):
    riw2d = np.zeros((lin,qt-1),float)  # Richardson Number (-)
    hzmid = np.zeros((lin,qt-1),float)  # z depth at middle point (m from the ref. level)

# wind parametrization
    sdur  = np.zeros((lin),float) # wind stress with duration filter
    sdir  = np.zeros((lin),float) # wind stress with duration and direction

    dw_hom  = np.zeros((lin),float)
    dw_lit = np.zeros((lin),float)  # W according to literature criteria 
    dw_spi  = np.zeros((lin),float) # W according to Spigel et al.
    
    wedu = np.zeros((lin),float) # W filtered by wind duration 
    wedi = np.zeros((lin),float) # W filtered by wind direction 

# -----------------------------------------------------------------------------
# COMPUTATION OF PARAMETERS AND INDICES 
# -----------------------------------------------------------------------------  
  
    print ("> ")
    root.update()
    load_time = time.time()
    print ("> Execution time for part I, II, and III: ")
    root.update() 
    print ("> "+str(round(load_time-start_time,4))+' seconds')
    root.update()   
    
# wind fetch - ls_fetch(min, mean, max):
    iw       = mod.velocityten(wind, rw)    # convert wind for 10 m high   
    dw_mean  = mod.wind_average(dw,iw)      # wind direction average
    
    hmean    = np.mean(h,axis=0) 

    ls_fetch = mod.wind_fetch(dw_mean,linang,angle,dists,type_length,dig)
            
    print ("--------------------------------------------------------------------------------------")
    root.update()
    print('> ')
    root.update()
    print ("> Part IV      Starting computation... ")
    root.update()  
    print('> ')
    root.update()

# temporal parameters:

    consecutive_dire, consecutive_wind, ver_dire, ver_wind = 0,0,0,0
    auxisa, auxisb, auxisc, auxisd = None, None, None, None
    war3,warmode1, warmode2 = 0,0,0
    error_thermo = 0

    for t in range(lin):
   
        auxtem = tempa[t,:]
        auxh   = h[t,:]  
        auxiw  = iw[t]        
        auxean = ean[t]    
        max_ean = np.nanmax(auxean)  
        min_ean = z0
        
        auxtem_ordered = mod.sorting_1d(auxtem)
    
        if (turn_iso == 1):
            if (tau[0] != -999):
                isoa[t] = mod.isotherms(tau[0],qt,auxh,auxtem_ordered,max_ean,min_ean,auxisa)
                auxisa = isoa[t]
            if (tau[1] != -999):
                isob[t] = mod.isotherms(tau[1],qt,auxh,auxtem_ordered,max_ean,min_ean,auxisb)
                auxisb = isob[t]
            if (tau[2] != -999):
                isoc[t] = mod.isotherms(tau[2],qt,auxh,auxtem_ordered,max_ean,min_ean,auxisc)
                auxisc = isoc[t]
            if (tau[3] != -999):
                isod[t] = mod.isotherms(tau[3],qt,auxh,auxtem_ordered,max_ean,min_ean,auxisd)
                auxisd = isod[t]
        
    # 2layer structure
        he[t], hh[t], pe[t], ph[t], glin[t], n[t], pu[t], pd[t], prox_thermo = mod.structure2layer(qt, auxh, auxtem_ordered, auxean, z0)      
        error_thermo = prox_thermo + error_thermo

        ht[t] = ean[t] - he[t]
        hH[t] = he[t]/(he[t]+hh[t])
        
        Hbar[t] = np.sqrt((he[t]+hh[t])/(he[t]*hh[t]))
            
        try:    
            Pbar[t] = np.sqrt(ph[t]/(ph[t]-pe[t]))
        except:
            try:
                Pbar[t] = None
            except:
                Pbar[t] = None
        
        if (hH[t] > 0.5):
            hH[t] = 1 - hH[t]
              
        strs[t], wast[t], riw[t] = mod.wind_parameters(auxiw, rw, pe[t], he[t], n[t], glin[t], auxean)
        
    # 2d parameters  
        p, n2, hmid, gli2d    = mod.thermal_stability(qt,auxh,auxean,auxtem)
        riw2 = mod.richardson(auxiw,rw,qt,auxh,pe[t],auxean,n2,hmid,p,gli2d)
    
        riw2d[t][:] = riw2[:]
        hzmid[t][:] = hmid[:]
    
    # wedderburn function (he = ze for a two layer system)
        wedd[t]     = mod.wedderburn(glin[t],he[t],wast[t],ls_fetch[1])  
        wedd_inv[t] = 1/wedd[t]
    
    # wind parameterization
        try:
            tbsiw = 2*ls_fetch[1]/np.sqrt(glin[t]*he[t]*hh[t]/(he[t]+hh[t])) # Merian
        except RuntimeWarning:
            tbsiw = None
            war.merian(dig)
      
        dw_hom[t] = None
        
        max_spigel = ls_fetch[1]*(he[t]+hh[t])/(4*he[t]*hh[t])
        min_spigel  =0.5*np.sqrt((he[t]+hh[t])/(he[t]*hh[t]))
        
        wedd_aux   = 1
        
        if min_spigel < 1:             # verify the minimum W for IW activity
            wedd_aux = min_spigel

        winsiz = int(0.25/2*tbsiw/dt/3600)
        points = t-winsiz
        polong = t-int(2*winsiz)
    
        
        
        if points < 1: points = 0
        if polong < 1: polong = 0
            
        try: wedd_mean = np.nanmean(wedd[points:t+1])
        except: wedd_mean = wedd[t]       
        
        if(wedd_mean < 20 and wedd_mean > wedd_aux ):
            
            try: dw_lit[t] = mod.wind_average(dw[points:t+1],iw[points:t+1])
            except: dw_lit[t] = dw[t]
                
        else: dw_lit[t] = None
    
        if (wedd_mean < max_spigel) and (wedd_mean > min_spigel):
            
            try: dw_spi[t] = mod.wind_average(dw[points:t+1],iw[points:t+1]) 
            except: dw_spi[t] = dw[t]
        
            if(consecutive_wind > 0):
                                
                wmin, wmax = mod.wind_angle(np.mean(dw[polong:polong+winsiz+1]),linang)
                
                dw10 = mod.wind_average(dw[points:t+1],iw[points:t+1]) 
                                
                if wmin > wmax:
                    if( dw10 >= wmin or dw10 <= wmax):

                        consecutive_dire = consecutive_dire + 1
                        dw_hom[t] = 357

                    else:
                        if(consecutive_dire > ver_dire):
                            ver_dire = consecutive_dire
                
                        consecutive_dire = 0                
                else:                    
                    if( dw10 >= wmin and dw10 <= wmax):
                    
                        consecutive_dire = consecutive_dire + 1
                        dw_hom[t] = 357
                        
                    else:
                        if(consecutive_dire > ver_dire):
                            ver_dire = consecutive_dire
                
                        consecutive_dire = 0                              
        
            consecutive_wind = consecutive_wind + 1 
        
            if (consecutive_wind > ver_wind):
                ver_wind = consecutive_wind  # number of samples of main wind
            
            if(consecutive_dire > ver_dire):
                ver_dire = consecutive_dire
        
        else:

            consecutive_dire = 0
            consecutive_wind = 0
            dw_spi[t] = None
        
    # 3layer structure    
        h1[t],h2[t],h3[t],p1[t],p2[t],p3[t] = mod.structure3layer(qt, auxh, auxtem, minval, auxean, z0)
        
        if h1[t] == -999:            
            war3  = war3 + 1 
            h1[t],h2[t],h3[t],p1[t],p2[t],p3[t] = mod.approx_layer(he[t],hh[t],pe[t],ph[t])
            
            

    
    # Generation of internal seiches according to Ri, hh, he, wind fetch 
    #
    # outputs: classification along time series [array int 1d] = 1, 2, 3 or 4
    #
    #    1 = Stratification is broken down by mixing
    #    2 = Large interface displacement occurs and is accompanied by interface shear and Kelvin-Helmholtz billows (strat is broken)
    #    3 = Internal seiche is dominant
    #    4 = Internal seiche with short amplitude.
    
        genera[t] = mod.class_generation(riw[t],hh[t],he[t],ls_fetch[1])
        iw_dw[t], iw_up[t] = mod.iw_generation(wedd[t],hh[t],he[t],ls_fetch[1])

        _,v1mode[t],_ = np.real(miw.disp_zmodel(pe[t],ph[t],he[t],hh[t],ls_fetch,1))
        _,v2mode[t],_ = np.real(miw.disp_xmodel3(p1[t],p2[t],p3[t],h1[t],h2[t],h3[t],ls_fetch,2,1))
    
    # it is difficult occurs a change in this order (5h in a short dt)
        if(t>0 and abs(v2mode[t] - v2mode[t-1])>5*60*60):
            v2mode[t] = v2mode[t-1]
            warmode2 = warmode2 +1

        if(t>0 and abs(v1mode[t] - v1mode[t-1])>5*60*60):
            v1mode[t] = v1mode[t-1]
            warmode1 = warmode1 + 1

    if warmode1 > 0: war.profile_structure(dig,1,100*warmode1/lin)
    if warmode2 > 0: war.profile_structure(dig,2,100*warmode1/lin)

    if war3 > 0: war.metalimnion(dig,100*war3/lin)
    if error_thermo > 0: war.thermocline(dig,100*error_thermo/lin) 
        
    iso = [isoa, isob, isoc, isod]

# period is devided into three shorter sub-periods
# h: he/H and Wedderburn number for each subperiod (confidence interval with 95%)
    
    group = [genera[i:i+int(lin/3)] for i in range(0, len(genera), int(lin/3))]
    hH_gp = [hH[i:i+int(lin/3)] for i in range(0, len(hH), int(lin/3))]
    wi_gp = [wedd_inv[i:i+int(lin/3)] for i in range(0, len(wedd_inv), int(lin/3))]

    
    h_lim1 = mod.ci(hH_gp[0])    # define the aveage of he/H of subperiod P1
    W_lim1 = mod.ci(wi_gp[0])    # deifne the aveage of W    of subperiod P1
    h_lim2 = mod.ci(hH_gp[1])    # deifne the aveage of he/H of subperiod P2
    W_lim2 = mod.ci(wi_gp[1])    # deifne the aveage of W    of subperiod P2
    h_lim3 = mod.ci(hH_gp[2])    # deifne the aveage of he/W of subperiod P3
    W_lim3 = mod.ci(wi_gp[2])    # deifne the aveage of W    of subperiod P3
    
#date - formating data
    dx    = [datetime.datetime.strptime(d, '%Y/%m/%d/%H/%M') for d in date]
    dx_gp = [dx[i:i+int(lin/3)] for i in range(0, len(dx), int(lin/3))]
    
# defining the subperiods
    P1    = [dx_gp[0][0], dx_gp[0][-1]]
    P2    = [dx_gp[1][0], dx_gp[1][-1]]
    P3    = [dx_gp[2][0], dx_gp[2][-1]]

# average of some variables

    mean_temp,low_temp,high_temp = mod.mean_confidence_interval(tempa)
    mean_h    = np.mean(h,axis=0)
    
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
    
    wedd_lim_lower = np.nanmean(ls_fetch[1]*(m_he+m_hh)/(4*m_he*m_hh))
    wedd_lim_upper = 0.5*np.sqrt((m_he+m_hh)/(m_he*m_hh))
    
    try:
        m_dw_spi= np.nanmean(dw_spi)
    except RuntimeWarning:
        war.homogeneous_condition(dig) 
        m_dw_spi = -999

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

        wedu[t] = wedd[t]/fdura
        wedi[t] = wedd[t]/fdire
    
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
       
        
    pv1h1 = np.array(miw.disp_zmodel(m_pe,m_ph,m_he,m_hh,ls_fetch,1))
    pv1h2 = np.array(miw.disp_zmodel(m_pe,m_ph,m_he,m_hh,ls_fetch,2))
    pv1h3 = np.array(miw.disp_zmodel(m_pe,m_ph,m_he,m_hh,ls_fetch,3))

    pv2h1 = np.array(miw.disp_xmodel3(m_p1,m_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,2,1))
    pv2h2 = np.array(miw.disp_xmodel3(m_p1,m_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,2,2))
    pv2h3 = np.array(miw.disp_xmodel3(m_p1,m_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,2,3))

    pv3h1 = np.array(miw.disp_xmodel4(m_p1,np2,np3,m_p3,nh1,nh2,nh3,nh4,ls_fetch,3,1))
    pv3h2 = np.array(miw.disp_xmodel4(m_p1,np2,np3,m_p3,nh1,nh2,nh3,nh4,ls_fetch,3,2))
    pv3h3 = np.array(miw.disp_xmodel4(m_p1,np2,np3,m_p3,nh1,nh2,nh3,nh4,ls_fetch,3,3))

#  internal wave frequency (without and with Earth rotation effect)
       
    omega_e = 7.2921*10**-5
    lat_rad = np.radians(lat)

    fc      = 2*omega_e*np.sin(lat_rad)
    fo      = abs(fc/(2*np.pi))

    if fc != 0:
        Bu = m_n**2*(m_he+m_hh)**2/(fc**2*ls_fetch[1]**2)
    else:
        Bu = 0

    f11, cf11  = 1/pv1h1, 1/miw.coriolis_effect(fc,pv1h1)
    f21, cf21  = 1/pv2h1, 1/miw.coriolis_effect(fc,pv2h1)
    f31, cf31  = 1/pv3h1, 1/miw.coriolis_effect(fc,pv3h1)

# model sensitivity using the hydrostatic model  

    sens_met = []
    delta = minval/2
    delta_variation = np.linspace(delta,delta+minval,15)
    
    for i in delta_variation:
        hus,hms,hbs,pus,pms,pbs = mod.structure3layer(qt, mean_h, mean_temp, i, np.mean(ean), z0)
        
        try:
            _,sens,_  =miw.disp_xmodel3(pus,pms,pbs,hus,hms,hbs,ls_fetch,2,1)
            sens_met.append(sens)
        except:
            sens_met.append(None)
    
             
    xpe, v1h1_spe = miw.sensitivity_2layer(m_pe,drho,lp_sen,m_pe,m_ph,m_he,m_hh,ls_fetch,1)
    xph, v1h1_sph = miw.sensitivity_2layer(m_ph,drho,lp_sen,m_pe,m_ph,m_he,m_hh,ls_fetch,2)
    xhe, v1h1_she = miw.sensitivity_2layer(m_he,ddep,lh_sen,m_pe,m_ph,m_he,m_hh,ls_fetch,3)  
    xhh, v1h1_shh = miw.sensitivity_2layer(m_hh,ddep,lh_sen,m_pe,m_ph,m_he,m_hh,ls_fetch,4)

    xp1, v2h1_sp1 = miw.sensitivity_3layer(m_p1,drho,lp_sen,m_p1,m_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,typ=1)
    xp2, v2h1_sp2 = miw.sensitivity_3layer(m_p2,drho,lp_sen,m_p1,m_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,typ=2)
    xh1, v2h1_sh1 = miw.sensitivity_3layer(m_h1,ddep,lp_sen,m_p1,m_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,typ=3)
    xh2, v2h1_sh2 = miw.sensitivity_3layer(m_h2,drho,lp_sen,m_p1,m_p2,m_p3,m_h1,m_h2,m_h3,ls_fetch,typ=4)

    rho_bar, dep_bar, Prho, Pdep = miw.sensitivity_dimension(ls_fetch[1], m_pe, m_ph, m_he, m_hh)
 
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

    if int(autoff_wsuser) == 1:
        if pv1h1[1] != None:   
            windsizeuser = 10*pv1h1[1]
        else:
            windsizeuser = 5*24*60*60
            war.welch_average(dig,100*warmode1/lin)              

    sensor_filtered = []
    timee           = []
    per             = []
    power_sensor    = []
    freqe           = []
    welch_sensor    = []
    wr_sensor       = []
    conf_sensor     = []
    
    if filter_process == 1 :  # define the range of the bandpass filter
        low1, high1 = 1/((pv1h1[2]/60)/60), 1/((pv1h1[0]/60)/60) # banded through internal wave model V1H1 (lower mode)

    else:
        low1, high1 = 1/high_per, 1/low_per
    
    if turn_temp == 1:
        
        
        for index in range(qt):   # spectral analysis for all sensors (multi-layer psd needs that information)
        
            new = temp[:,index]
            fs = 1/dt
            
            try:
                filtered_band = mod.butter_bandpass_filter(new, low1, high1, fs)
            except ValueError:
                filtered_band = None
                war.bandpass(dig,'sensor')
            
 
            aux_time, aux_per, aux_power = mod.wave_spectral ( new, dt, mother)
            aux_freq, aux_wl_per, aux_welch, aux_wr, aux_conf = mod.welch_method(new,windsizeuser,window, dt)
        
            sensor_filtered.append(filtered_band)

            timee.append(aux_time)
            per.append(aux_per)
            power_sensor.append(aux_power)
            freqe.append(aux_freq) 
            welch_sensor.append(aux_welch)
            wr_sensor.append(aux_wr)
            conf_sensor.append([aux_conf])

        
        
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
    

    if filter_process == 1 :        # define the range of the bandpass filter
        lowcut, highcut  = 1/((pv1h1[2]/60)/60), 1/((pv1h1[0]/60)/60)
    else:
        lowcut, highcut  =  1/high_per, 1/low_per

    freq      = []
    power     = []
    time_temp = []
    per_temp  = []
    band      = []
    welch     = []
    wr        = []
    conf      = []

    iso_corrected = np.zeros(4,float)


    if (turn_iso == 1):
        ean_norm  = ean  - np.mean(ean)
        
        for i in range(4):
            
            if np.isnan(iso[i]).all() == True and tau[i]!=-999:
                war.isotherm_boundary(dig, str(tau[i]))
                tau[i] = -999        
            
            if(tau[i] != -999):
                
                aux_time, aux_per, aux_power = mod.wave_spectral (iso[i], dt, mother)
                aux_freq, _, aux_welch, aux_wr, aux_conf = mod.welch_method(iso[i],windsizeuser,window, dt)
                    
                
                fs = 1/dt
                
                try:
                    aux_band = mod.butter_bandpass_filter(iso[i], lowcut, highcut, fs)
                except ValueError:
                    aux_band = None
                    war.bandpass(dig,'isotherm')
                
                iso_corrected[i] = np.mean(abs((iso[i] - np.mean(iso[i])) - ean_norm))
    
                time_temp.append(aux_time)
                per_temp.append(aux_per)
                power.append(aux_power)
                freq.append(aux_freq) 
                welch.append(aux_welch)  
                band.append(aux_band)
                wr.append(aux_wr)
                conf.append(aux_conf)
            else:
                time_temp.append(0)
                per_temp.append(0)
                power.append(0)
                freq.append(0) 
                welch.append(0)  
                band.append(0)
                wr.append(0)
                conf.append([0])
    
    banda = []
    for i in range(4):
        banda.append(np.nanmax(abs(band[i])))

    amax = max(np.array(banda))
    aind = tau[np.where(banda == amax)]
         
       
# spectral analysis of the solar radiation
    if rad == 1:
        solar_ws = 5*24*60*60
        time_sol, per_sol, power_sol     = mod.wave_spectral ( ra, dt, mother)
        freq_sol, wl_aper_sol, welch_sol, wr_sol, conf_sol = mod.welch_method (ra, solar_ws, window, dt)

#spectral analysis of wind intensity

    time_win, per_win, power_win     = mod.wave_spectral ( iw, dt, mother)
    freq_win, wl_aper_win, welch_win, wr_win, conf_win = mod.welch_method (iw, windsizeuser, window, dt)

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
                new = temp[:,seu[i]]
                a1, a2, a3, a4 = mod.coherence_shift(new, iw, windsizeuser, dt)
            
                phws.append(a1)
                caws.append(a2)
                faws.append(a3)
                c95ws.append(a4)
            
                s_filtered.append(sensor_filtered[seu[i]])
                    
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
                
                        ni = temp[:,seu[i]]
                        nj = temp[:,seu[j]]
                
                        ph_aux, coh, f, c95 = mod.coherence_shift(ni,nj,windsizeuser,dt)
                
                        phij.append(ph_aux)
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

        phw, caw, faw, c95aw = mod.coherence_shift(iso[aux_iso], iw, windsizeuser, dt)
        
        if rad == 1:
            phr, car, far, c95ar = mod.coherence_shift(iso[aux_iso], ra, windsizeuser, dt)   


        phiso = []
        coiso = []
        friso = []
        cliso = []

    # coherence and phase analysis between isotherms
        for i in range(4):
            for j in range(4):
                if ( i != j and j > i):
                    if(tau[i] != -999) and (tau[j] != -999):
                    
                        phaux, coaux, fiaux, c9aux = mod.coherence_shift(iso[i],iso[j],windsizeuser,dt)
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
#   teml_ordered - ordered temperature profile                        
#  

    templ, hl = graph.correction_contourf(h,temp,dj,lin,dci,qt)
    riwl, hlm = graph.correction_contourf(hzmid,riw2d,dj,lin,dci,qt-1) 
    
    tt, ithermo = mod.find_nearest(hl,m_ht) 
    templ_ordered = mod.sorting_2d(templ)
    
    ri_wind = np.nanmin(mod.average(riw,(0.25/2)*(dura/(60*60)/dt)))
    we_wind = np.nanmin(mod.average(wedd,(0.25/2)*(dura/(60*60)/dt)))
       
# outputs:
# tt: thermocline depth according to countourf interpolation
# ithermo: indice on 'hl' 
# -----------------------------------------------------------------------------
    
    
# -----------------------------------------------------------------------------
# Thorpe scale
        
    tho = np.zeros((lin,len(templ[0])),float)
    
    for t in range(lin):
        
        thermo_temp[t] = templ[t][ithermo]
        
        for z in range(qt):
            _, index_thorpe = mod.find_nearest(templ_ordered[t,:],templ[t,z])
            tho[t,z] = hl[z] - hl[index_thorpe]     # Thorpe scale    


# spectral analysis of the thermocline 

    tthermo, pthermo, powerthermo = mod.wave_spectral (thermo_temp, dt, mother)
    
    freqthermo, wlthermo, welchthermo, wrthermo, confthermo = mod.welch_method(thermo_temp,windsizeuser,window, dt)
    freq_ean, wl_aper_ean, welch_ean, wr_ean, conf_ean  = mod.welch_method(ean,windsizeuser,window, dt)

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
                np.savetxt(output_path+'textfiles/spectral_isotherms'+str(tau[i])+'.txt', np.column_stack((freq[i],welch[i])), delimiter='\t', header='freq(Hz)\tPSD(oC²/Hz)', fmt='%0.8f %0.15f')                
                np.savetxt(output_path+'textfiles/isotherms'+str(tau[i])+'.txt',np.column_stack((time_temp[i],iso[i])), delimiter='\t', header='time(hour)\tisotherm(oC)', fmt='%0.8f %0.5f')
                np.savetxt(output_path+'textfiles/rednoise_isotherms'+str(tau[i])+'.txt', np.column_stack((wr[i],conf[i])), delimiter='\t', header='freq(Hz)\tPSD(oC²/Hz)', fmt='%0.8f %0.15f')
                
    if turn_temp == 1:
        for i in range(4):
            if(sen[i]==1):
                np.savetxt(output_path+'textfiles/spectral_sensor'+str(seu[i])+'.txt',np.column_stack((freqe[seu[i]],welch_sensor[seu[i]])), delimiter='\t', header='freq(Hz)\tPSD(m²/Hz)', fmt='%0.8f %0.15f')
               
    np.savetxt(output_path+'textfiles/wind.txt',np.column_stack((time_win,dw,iw,strs)), delimiter='\t', header='time(hour)\tdirection(o)\tspeed(m/s)\tstress(N/m²)', fmt='%0.8f %0.1f %0.3f %0.6f')
    np.savetxt(output_path+'textfiles/spectral_wind.txt',np.column_stack((wl_aper_win,welch_win)), delimiter='\t', header='period(hour)\tPSD wind((m/s)²/Hz)', fmt='%0.3f %0.5f')
    np.savetxt(output_path+'textfiles/stability.txt',np.column_stack((time_win,riw,wedd)), delimiter='\t', header='time(hour)\t Ri(-)\tWedderburn(-)', fmt='%0.3f %0.8f %0.8f')
    np.savetxt(output_path+'textfiles/thermocline.txt',np.column_stack((time_win,he)), delimiter='\t', header='time(hour)\tthermocline depth(m)', fmt='%0.3f %0.4f')
    np.savetxt(output_path+'textfiles/metalimnion_thickness.txt',np.column_stack((time_win,h2)), delimiter='\t', header='time(hour)\tmetalimnion thickness(m)', fmt='%0.3f %0.4f')
    np.savetxt(output_path+'textfiles/internal_periods.txt',np.column_stack((time_win,v1mode, v2mode)), delimiter='\t', header='time(hour)\tV1H1 period(s)\tV2H1 period(s)', fmt='%0.3f %0.4f %0.4f')
    np.savetxt(output_path+'textfiles/mean_profile.txt',np.column_stack((mean_h,mean_temp)), delimiter='\t', header='mab\ttemperature (dC)', fmt='%0.3f %0.4f')
    
    if rad == 1:
        np.savetxt(output_path+'textfiles/spectral_solar.txt',np.column_stack((wl_aper_sol,welch_sol)), delimiter='\t', header='period(hour)\tPSD sw((W/m²)²/Hz)', fmt='%0.3f %0.5f')

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
                graph.wavelet_iso(dx,per_temp[i], power[i], show_on, ax2)
                
                plt.savefig(output_path+'iso'+str(i)+'.png', dpi = depi)
    

    plt.figure(figsize=(8,6))
    
    ax1 = plt.subplot2grid((2,3),(0,0),colspan=2)
    ax2 = plt.subplot2grid((2,3),(1,0),colspan=2)
    ax3 = plt.subplot2grid((2,3),(1,2),sharey=ax2)

    ax1.set_title('(a) thermocline analysis',loc='left')
    ax2.set_title('(b)',loc='left')
    ax3.set_title('(c)',loc='left')

    graph.thermocline(dx, tthermo,thermo_temp,'navy',round(m_ht, 2),ax1)
    graph.wavelet_iso(dx, pthermo, powerthermo, show_on, ax2)
    graph.psd_thermocline(wlthermo,welchthermo,'navy',int(tt),wrthermo,confthermo,ax3)

    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.savefig(output_path+'thermocline_analysis.png', dpi = depi)



    plt.figure(figsize=(8,3))

    ax1 = plt.subplot2grid((1,1),(0,0))
    graph.wind_direction(dx, time_win, dw, dw_spi,dw_lit,dw_hom,wedd_lim_upper,wedd_lim_lower,ax1)

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


    plt.figure(figsize=(6,5))
    ax1 = plt.subplot2grid((1,1),(0,0))

    graph.psd_wind(wl_aper_win, welch_win, wr_win, conf_win, ax1)
    
    if rad ==1:
        ax2 = ax1.twinx() 
        graph.psd_sola(wl_aper_sol, welch_sol, wr_sol, conf_sol, ax2)

    plt.savefig(output_path+'meteo_spectra.png', dpi = depi)


    plt.figure(figsize=(5,6))
    ax1 = plt.subplot2grid((1,1),(0,0))

    graph.psd_level(freq_ean, welch_ean, wr_ean, conf_ean, ax1)

    plt.savefig(output_path+'level_spectra.png', dpi = depi)


    if turn_iso == 1:

        
        plt.figure(figsize=(8,6))
        
        ax1 = plt.subplot2grid((1,2),(0,0))
        ax2 = plt.subplot2grid((1,2),(0,1),sharey=ax1)
        
        ax1.set_title('(a)',loc='left')
        ax2.set_title('(b)',loc='left')

        graph.psd_iso(tau, freq, welch, f11, f21, f31, largelen,m_n,fo,wr,conf, ax1)
        
        graph.coherence (caw,faw,tau,aux_iso, largelen, ax2)

        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.savefig(output_path+'psd_nonhydro.png', dpi = depi)


        plt.figure(figsize=(8,6))

        ax1 = plt.subplot2grid((1,2),(0,0))
        ax2 = plt.subplot2grid((1,2),(0,1),sharey=ax1)

        ax1.set_title('(a) ',loc='left')
        ax2.set_title('(b) Coriolis correction',loc='left')

        graph.psd_iso(tau, freq, welch, f11, f21, f31,largelen,m_n,fo,wr,conf, ax1)

        # this has the earth rotation effect
        graph.psd_iso(tau, freq, welch, cf11, cf21, cf31,largelen,m_n,fo,wr,conf, ax2)

        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.savefig(output_path+'psd_hydro_coriois.png', dpi = depi)


        plt.figure(figsize=(8,9))

        ax1 = plt.subplot2grid((2,1),(0,0))
        ax2 = plt.subplot2grid((2,1),(1,0),sharex=ax1)

        ax1.set_title('(a) ',loc='left')
        ax2.set_title('(b) ',loc='left')

        graph.psd_iso_conv(tau, freq, welch, f11, f21, f31,m_n,fo,wr,conf,'psd', ax2)
        graph.psd_iso_conv(tau, freq, welch, cf11, cf21, cf31,m_n,fo,wr,conf,'var', ax1)

        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.savefig(output_path+'psd.png', dpi = depi)


    plt.figure( figsize=(8,6))

    ax1 = plt.subplot2grid((3,2), (0,0), rowspan=3)
    ax2 = plt.subplot2grid((3,2), (0,1))
    ax3 = plt.subplot2grid((3,2), (1,1))
    ax4 = plt.subplot2grid((3,2), (2,1))

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b) - P1',loc='left')
    ax3.set_title('(c) - P2',loc='left')
    ax4.set_title('(d) - P3',loc='left')

    graph.degeneration(P1,P2,P3,h_lim1,h_lim2,h_lim3,W_lim1,W_lim2,W_lim3,m_pe,m_ph,m_h2,m_he+m_hh,ls_fetch[1],ax1)

    graph.generation(genera_0,ax2)
    graph.generation(genera_1,ax3)
    graph.generation(genera_2,ax4)

    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.savefig(output_path+'degenera.png', dpi = depi)


    plt.figure(figsize=(8,6))  
    ax1 = plt.subplot2grid((1,2),(0,0))  
    ax2 = plt.subplot2grid((1,2),(0,1))                  

    ax1.set_title('(a)',loc='left')  
    ax2.set_title('(b)',loc='left')
    
    mode1 = pv1h1[1]/3600/2
    average_wedda   = np.nanmin(mod.average(wedi,mode1*0.25/2/dt))
    parh  = m_he/(m_he+m_hh) 
    

    g   = 12.156*(parh**3)-15.714*(parh**2)+2.8426*parh+2.0846
    f   = g*np.exp((parh**2)/0.25)
    paramet = 2*f**2*np.log(amax/(0.1*m_he))

    graph.classification_genera(average_wedda, we_wind, parh,amax/m_he,aind, ax1)
    graph.bueno_parameterization(average_wedda, we_wind, paramet, aind, ax2)

    plt.tight_layout()
    plt.savefig(output_path+'classification_evolution.png', dpi = depi)



    plt.figure( figsize=(8,6))

    ax1 = plt.subplot2grid((3,2), (0,0), rowspan=3)
    ax2 = plt.subplot2grid((3,2), (0,1), rowspan=3)


    ax1.set_title('(a) - Original Chart',loc='left')
    ax2.set_title('(b) - Zooming',loc='left')


    graph.degeneration_evolution(P1,P2,P3,hH_gp,wi_gp,m_pe,m_ph,m_h2,m_he+m_hh,ls_fetch[1],'no',ax1)
    graph.degeneration_evolution(P1,P2,P3,hH_gp,wi_gp,m_pe,m_ph,m_h2,m_he+m_hh,ls_fetch[1],'yes',ax2)

    plt.savefig(output_path+'degenera_evolution.png', dpi = depi)



    plt.figure( figsize=(8,6))

    ax1 = plt.subplot2grid((3,1), (0,0))
    ax2 = plt.subplot2grid((3,1), (1,0), rowspan=3, sharex=ax1)

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')

    graph.wedd_limit(dx,time_win,wedd,iw_dw, iw_up, ax1)
    graph.richardson2d(dx,time_win,hlm,riwl,ht,ean,z0,ax2)
    


    plt.setp(ax1.get_xticklabels(), visible=False)
    
    plt.savefig(output_path+'richardson.png',dpi = depi)


    plt.figure(figsize=(8,6))

    ax1 = plt.subplot2grid((3,1),(0,0))
    ax2 = plt.subplot2grid((3,1),(1,0), rowspan=3, sharex=ax1)

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')

    graph.density(dx,time_win,pu,pd,pe,ph,ax1)
    graph.tempstructure_zoom(dx,time_win,hl,ht,templ,ean,z0,ax2)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.savefig(output_path+'structure_thermo.png',dpi = depi)

    plt.figure(figsize=(8,2.3))

    ax1 = plt.subplot2grid((1,1),(0,0))

    graph.wind(dx, time_win, dw, iw, ax1)
    plt.savefig(output_path+'wind.png',dpi = depi)


    plt.figure( figsize=(8,4))

    ax1 = plt.subplot2grid((1,1),(0,0))

    graph.tempstructure(dx,time_win,hl,ht,templ,ean,z0,ax1)
    plt.savefig(output_path+'structure.png',dpi = depi)



    plt.figure( figsize=(8,3.2))
    ax1 = plt.subplot2grid((1,1),(0,0))

    graph.thorpe_scale(dx,time_win,hl,tho,ean,z0,ax1)
    plt.savefig(output_path+'thorpe_scale.png',dpi = depi)


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
    graph.depth_sensitivity(xhh,v1h1_shh,r'$h_h$',1,ax3)
    graph.densi_sensitivity(xph,v1h1_sph,r'$\rho_h$',1,ax4)

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

    graph.depth_sensitivity(xh1,v2h1_sh1,r'$h_1$',0,ax1)
    graph.densi_sensitivity(xp1,v2h1_sp1,r'$\rho_1$',0,ax2)
    graph.depth_sensitivity(xh2,v2h1_sh2,r'$h_2$',1,ax3)
    graph.densi_sensitivity(xp2,v2h1_sp2,r'$\rho_2$',1,ax4)

    plt.savefig(output_path+'sensitivity_multiv.png', dpi = depi)

    plt.figure(figsize=(8,4))

    ax1 = plt.subplot2grid((1,2),(0,0))
    ax2 = plt.subplot2grid((1,2),(0,1))

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')

    graph.parabar_sensitivity(dep_bar,Pdep,'dep',pv1h1[0],Hbar,pv1h1[2],ax1)
    graph.parabar_sensitivity(rho_bar,Prho,'rho',pv1h1[0],Pbar,pv1h1[2],ax2)

    plt.savefig(output_path+'sensitivity_barparameters.png', dpi = depi)

    plt.figure(figsize=(9,4))

    ax1 = plt.subplot2grid((1,2),(0,0))
    ax2 = plt.subplot2grid((1,2),(0,1))

    ax1.set_title('(a)',loc='left')
    ax2.set_title('(b)',loc='left')

    graph.averaged_profile(mean_temp,low_temp,high_temp,mean_h, z0, ax1)
    graph.sensibility_metalimnion(delta_variation,sens_met,minval,ax2)

    plt.savefig(output_path+'mean_temperature.png', dpi = depi)

    plt.figure(figsize=(8,3))

    ax1 = plt.subplot2grid((1,1),(0,0))
    graph.wedd_compb(dx,time_win, wedd, wedi, iw_dw, iw_up, ax1)

    plt.savefig(output_path+'wedd_filtering.png',dpi = depi)

    if turn_iso == 1:
        plt.figure(figsize=(8,4))
        
        ax1 = plt.subplot2grid((1,1),(0,0))
        graph.multi_isotherms(dx, iso, tau, time_temp, z0, ax1)
        
        plt.savefig(output_path+'isotherms.png',dpi = depi)

        if(ana1 != -999 or ana2 != -999):
            fig, ax1 = plt.subplots(figsize=(8,4))
            graph.coherence_iso(tau,coiso,friso,ana1,ana2, ax1)   
        
            ax2 = ax1.twinx()
        
            graph.phase_iso(tau,phiso,friso,cliso,ana1,ana2, ax2)   
        

            plt.savefig(output_path+'coherence.png',dpi = depi)

    plt.figure(figsize=(8,3))
    ax1 = plt.subplot2grid((1,1),(0,0))

    graph.wavelet_resona(dx,per_win, power_win, v1mode, v2mode, ax1)
    plt.savefig(output_path+'wind_resonance.png',dpi = depi)


    if turn_iso == 1:

        
        plt.figure(figsize=(8,4))
        ax1 = plt.subplot2grid((1,1),(0,0))
        if filter_process == 1:
            ax1.set_title('Band-pass filter - Bandwidth: '+str(round(pv1h1[2]/60/60,1))+' h to '+str(round(pv1h1[0]/60/60,1))+' h',loc='left')
        
        elif filter_process == 2:
            ax1.set_title('Band-pass filter - Bandwidth: '+str(round(low_per,1))+' h to '+str(round(high_per,1))+' h',loc='left')
    
    
        graph.temp_bandpass(dx, band, tau, ax1)


        
        plt.setp(ax4.get_xticklabels(), visible=False)
        plt.savefig(output_path+'iso_bandpass.png',dpi = depi)


    if turn_temp == 1:

        for jd in range(4):
        
            if seu[jd] != -999:
                plt.figure( figsize=(8,6))
                ax1 = plt.subplot2grid((2,2),(0,0),colspan=2)
                ax2 = plt.subplot2grid((2,2),(1,0),colspan=2,sharex=ax1)
    
                ax1.set_title('(a)',loc='left')
                ax2.set_title('(b)',loc='left')
    
                ids = seu[jd]
                dep = temp[:,ids]
    
                graph.temperature(dx, timee[ids],dep,'navy',depth[jd],ax1)
                graph.wavelet_depth(dx, per[ids], power_sensor[ids], timee[ids], ax2)

                plt.savefig(output_path+'depth'+str(int(jd))+'.png', dpi =depi)
    

        plt.figure( figsize=(8,6))
        ax1 = plt.subplot2grid((1,1),(0,0))
        
        graph.psd_multilayer(freqe[0],hmean[::-1],welch_sensor,fo,m_n,z0,ax1)
        plt.savefig(output_path+'psd_multi.png', dpi = depi)

        plt.figure( figsize=(8,6))

        ax1 = plt.subplot2grid((1,2),(0,0))
        ax2 = plt.subplot2grid((1,2),(0,1),sharey=ax1)

        ax1.set_title('(a)',loc='left')
        ax2.set_title('(b)',loc='left')
    
        graph.psd_depth(sen,depth, freqe, welch_sensor, f11, f21, f31,largelen,m_n,fo,wr_sensor, conf_sensor, ax1)
        graph.psd_depth(sen,depth, freqe, welch_sensor, cf11, cf21, cf31,largelen,m_n,fo,wr_sensor,conf_sensor, ax2)

        plt.setp(ax2.get_yticklabels(), visible=False)

        plt.savefig(output_path+'psd_coriolis_depth.png', dpi = depi)


        plt.figure( figsize=(6,6))
        ax1 = plt.subplot2grid((1,1),(0,0))

        graph.psd_depth(sen,depth, freqe, welch_sensor, f11, f21, f31,largelen,m_n,fo,wr_sensor,conf_sensor, ax1)


        plt.savefig(output_path+'psd_hydro_depth.png', dpi = depi)


        fig, ax1 = plt.subplots(figsize=(8,5))
        graph.coherence_depth(depth,cohij,fij,ax1)
    
        ax2 = ax1.twinx()
        graph.phase_depth(depth, c95ij, phij,fij,ax2)

        plt.savefig(output_path+'phase_depth.png',dpi = depi)


        plt.figure( figsize=(8,4))
        ax1 = plt.subplot2grid((1,1),(0,0))
        
        graph.thermal_variation(dx, depth,seu,temp,time_win, ax1)

        plt.savefig(output_path+'temperature_depth.png',dpi = depi)
    
    
        plt.figure(figsize=(8,4))

        ax1 = plt.subplot2grid((1,1),(0,0))

        if filter_process == 1:
            ax1.set_title('Band-pass filter - Bandwidth: '+str(round(pv1h1[2]/60/60,1))+' h to '+str(round(pv1h1[0]/60/60,1))+' h',loc='left')
        
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
    canvas.drawString(520,800,'Page 1/4')


    logo =ImageReader('0interwave.png')
    canvas.drawImage(logo, 30, 733,width=120, height=69.856)

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


    canvas.drawString(30, 425, 'Wind parameters: ')
    canvas.drawString(30, 405, 'Duration of the strongest wind event:')
    canvas.drawString(210, 405,  str(round(dura/(60*60),2))+' h')
    canvas.drawString(30, 385, 'Just considering homogeneous direction:')
    canvas.drawString(30, 370,  str(round(dire/(60*60),2))+' h blowing '+str(round(m_dw_spi,0))+'° (c. nautica)')
    canvas.drawString(30, 350, 'Reduction factor:')
    canvas.drawString(120, 350, 'Duration factor:')
    canvas.drawString(210, 350,  str(round(fdura,3)))
    canvas.drawString(120, 335, 'Direction factor:')
    canvas.drawString(210, 335,  str(round(fdire,3)))    
    
    g2 =ImageReader(output_path+'structure.png')
    canvas.drawImage(g2, 250, 405,width=300, height=150)

    canvas.drawString(30, 300, 'Mean friction velocity of the wind:')
    canvas.drawString(190, 300,  "{:.2e}".format(round(m_wast,5))+' m/s')
    canvas.drawString(30, 285, 'Min friction velocity of the wind:')
    canvas.drawString(190, 285,  "{:.2e}".format(round(np.nanmin(wast),5))+' m/s')
    canvas.drawString(30, 270, 'Max friction velocity of the wind:')
    canvas.drawString(190, 270,  "{:.2e}".format(round(np.nanmax(wast),5))+' m/s')

    canvas.drawString(30, 250, 'Time of wind events favoring BSIW (%)')
    canvas.drawString(30, 230, 'Accoding to W')
    canvas.drawString(150, 230,  str(round(100*(1-np.isnan(dw_spi).sum()/len(dx)),2))+ '%')
    canvas.drawString(30, 215, 'Taking direction:')
    canvas.drawString(150, 215,  str(round(100*(1-np.isnan(dw_hom).sum()/len(dx)),2))+ '%') 

    g2 =ImageReader(output_path+'wind.png')
    canvas.drawImage(g2, 250, 555,width=300, height=86.25)

    g1 =ImageReader(output_path+'mean_temperature.png')
    canvas.drawImage(g1, 250, 251,width=300, height=133.33)

    g3 =ImageReader(output_path+'temporal_analysis.png')
    canvas.drawImage(g3, 250, 127.5,width=300, height=112.5)

    canvas.drawString(30, 170, 'Paramters of Stability:')
    canvas.drawString(30, 150, 'Reduced Gravity:')
    canvas.drawString(120, 150,  "{:.2e}".format(round(m_glin,4))+' \xb1 '+"{:.2e}".format(round(mod.ciout(glin),4))+' m/s²')
    canvas.drawString(30, 135, 'Brunt-Vaisalla:')
    canvas.drawString(120, 135,  "{:.2e}".format(round(m_n,4))+' \xb1 '+"{:.2e}".format(round(mod.ciout(n),4))+' Hz')
    canvas.drawString(30, 120, 'Averaged Richardson number:')
    canvas.drawString(180, 120, "{:.2e}".format(round(m_riw,0)))

    canvas.drawString(30, 90, 'Time averaged (thermocline assuming V1 mode) :')
    canvas.drawString(30, 70,  'Richardson number:')
    canvas.drawString(140, 70,  "{:.2e}".format(round(np.nanmean(mod.average(riw,mode1*0.25/2/dt)),0))+' \xb1 '+"{:.2e}".format(round(np.std(mod.average(riw,mode1*0.25/2/dt))*1.96/np.sqrt(dt*lin/(mode1*0.25/2)),0))) # 1.96 = 95%

    canvas.drawString(30, 55, 'Wedderburn number:')
    canvas.drawString(140,55,  str(round(np.nanmean(mod.average(wedd,mode1*0.25/2/dt)),3))+' \xb1 '+"{:.2e}".format(round(np.std(mod.average(wedd,mode1*0.25/2/dt))*1.96/np.sqrt(dt*lin/(mode1*0.25/2)),3)))

    canvas.drawString(280,90,'Filtered time averaged Wedderburn number:')
    canvas.drawString(280,70,'Filtered by lake mixing criteria:')
    canvas.drawString(460,70, str(round(np.nanmin(mod.average(wedu,mode1*0.25/2/dt)),3)))
    canvas.drawString(280,55,'Filtered by wind homogeneity:')
    canvas.drawString(460,55, str(round(np.nanmin(mod.average(wedi,mode1*0.25/2/dt)),3)))


    canvas.showPage()

    canvas.setLineWidth(.3)
    canvas.setFont('Helvetica',10)
    canvas.drawString(520,800,'Page 2/4')


    g13 =ImageReader(output_path+'degenera.png')
    canvas.drawImage(g13, 20, 240,width=300, height=225)

    g13 =ImageReader(output_path+'thorpe_scale.png')
    canvas.drawImage(g13, 20, 455,width=300, height=120)

    g13 =ImageReader(output_path+'richardson.png')
    canvas.drawImage(g13, 20, 565,width=300, height=225)

    g7 =ImageReader(output_path+'classification_evolution.png')
    canvas.drawImage(g7, 300, 580,width=280, height=210)

    if int(0.5*(dura/(60*60)/dt)) == 0:
        war.average(dig)

 
    canvas.drawString(340,550,'Stability associated to strongest wind event:')
    canvas.drawString(340,530,  'Richardson number:')
    canvas.drawString(460,530,  str(round(np.nanmin(ri_wind),3)))

    canvas.drawString(340,515, 'Wedderburn number:')
    canvas.drawString(460,515,  str(round(np.nanmin(we_wind),3)))
    
    
    
    ampli_spigel = mod.spigel(m_he,m_he+m_hh,average_wedda,we_wind) 
    ampli_bueno  = mod.bueno(m_he,m_hh,average_wedda,we_wind,'amplitude')
    
    canvas.drawString(340,490,'BSIW amplitude according theories:')
    canvas.drawString(340,470,'Spigel and Imberger (1980):')
    canvas.drawString(480,470, str(round(ampli_spigel,2))+' m')

    canvas.drawString(340,455,'Bueno et al. (2020):')
    canvas.drawString(480,455, str(round(ampli_bueno,2))+' m')
    
    canvas.drawString(340,435,'Surface seiche amplitude:')
    canvas.drawString(480,435, str(round(1000*ampli_spigel*m_glin/9.81,2))+' mm')

    canvas.drawString(340,410,'Generation & Degeneration Theory¹:')

    canvas.drawString(340,395,'Periods')
    canvas.drawString(400,395,'1/W')
    canvas.drawString(460,395,'he/H:')

    canvas.drawString(340,375,'P1')
    canvas.drawString(400,375, str(round(np.nanmax(wi_gp[0]),5)))
    canvas.drawString(460,375, str(round(np.nanmean(hH_gp[0]),5)))

    canvas.drawString(340,360,'P2')
    canvas.drawString(400,360, str(round(np.nanmax(wi_gp[1]),5)))
    canvas.drawString(460,360, str(round(np.nanmean(hH_gp[1]),5)))

    canvas.drawString(340,345,'P3')
    canvas.drawString(400,345, str(round(np.nanmax(wi_gp[2]),5)))
    canvas.drawString(460,345, str(round(np.nanmean(hH_gp[2]),5)))

    canvas.drawString(340,325,'Maximum amplitude of BSIW:')
    canvas.drawString(340,310,'P1')
    canvas.drawString(400,310, str(round(np.nanmax(wi_gp[0])*m_he/2,3))+' m')
    canvas.drawString(340,295,'P2')
    canvas.drawString(400,295, str(round(np.nanmax(wi_gp[1])*m_he/2,3))+' m')
    canvas.drawString(340,280,'P3')
    canvas.drawString(400,280, str(round(np.nanmax(wi_gp[2])*m_he/2,3))+' m')

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
    canvas.drawString(520,800,'Page 3/4')


    canvas.drawString(50,750,'Analysis of the Thermocline and Theoretical results')


    g4 =ImageReader(output_path+'thermocline_analysis.png')
    canvas.drawImage(g4, 0, 380,width=550, height=412.5)
    canvas.drawString(50,390,'¹ Wind criteria must be > 0.25 to excite BSIW')

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
    
    canvas.drawString(370,610, 'Burger Number:')
    canvas.drawString(450,610, str(round(Bu,2)))
    
    if Bu <= 1:
        canvas.drawString(370,595, 'Coriolis effect must be considered' )
    else:
        canvas.drawString(370,595, 'Coriolis effect can be neglected')
    
    canvas.drawString(100, 370, 'Hydrostatic Model for the first three vertical and horizontal modes:')

    canvas.drawString(50, 350, 'V1H1' )
    canvas.drawString(90, 350,  str(round(pv1h1[1]/(60*60),2))+' h')
    canvas.drawString(140, 350, ' \xb1 ')
    canvas.drawString(160, 350, str(round((pv1h1[1]-pv1h1[0])/(60*60),2))+' h')

    canvas.drawString(50, 335, 'V2H1' )
    canvas.drawString(90, 335,  str(round(pv2h1[1]/(60*60),2))+' h')
    canvas.drawString(140, 335, ' \xb1 ')
    canvas.drawString(160, 335, str(round((pv2h1[1]-pv2h1[0])/(60*60),2))+' h')

    canvas.drawString(50, 320, 'V3H1' )
    canvas.drawString(90, 320,  str(round(pv3h1[1]/(60*60),2))+' h')
    canvas.drawString(140, 320, ' \xb1 ')
    canvas.drawString(160, 320, str(round((pv3h1[1]-pv3h1[0])/(60*60),2))+' h')


    canvas.drawString(220, 350, 'V1H2' )
    canvas.drawString(260, 350,  str(round(pv1h2[1]/(60*60),2))+' h')
    canvas.drawString(310, 350, ' \xb1 ')
    canvas.drawString(330, 350, str(round((pv1h2[1]-pv1h2[0])/(60*60),2))+' h')

    canvas.drawString(220, 335, 'V2H2' )
    canvas.drawString(260, 335,  str(round(pv2h2[1]/(60*60),2))+' h')
    canvas.drawString(310, 335, ' \xb1 ')
    canvas.drawString(330, 335, str(round((pv2h2[1]-pv2h2[0])/(60*60),2))+' h')

    canvas.drawString(220, 320, 'V3H2' )
    canvas.drawString(260, 320,  str(round(pv3h2[1]/(60*60),2))+' h')
    canvas.drawString(310, 320, ' \xb1 ')
    canvas.drawString(330, 320, str(round((pv3h2[1]-pv3h2[0])/(60*60),2))+' h')


    canvas.drawString(390, 350, 'V1H3' )
    canvas.drawString(430, 350,  str(round(pv1h3[1]/(60*60),2))+' h')
    canvas.drawString(470, 350, ' \xb1 ')
    canvas.drawString(490, 350, str(round((pv1h3[1]-pv1h3[0])/(60*60),2))+' h')

    canvas.drawString(390, 335, 'V2H3' )
    canvas.drawString(430, 335,  str(round(pv2h3[1]/(60*60),2))+' h')
    canvas.drawString(470, 335, ' \xb1 ')
    canvas.drawString(490, 335, str(round((pv2h3[1]-pv2h3[0])/(60*60),2))+' h')

    canvas.drawString(390, 320, 'V3H3' )
    canvas.drawString(430, 320,  str(round(pv3h3[1]/(60*60),2))+' h')
    canvas.drawString(470, 320, ' \xb1 ')
    canvas.drawString(490, 320, str(round((pv3h3[1]-pv3h3[0])/(60*60),2))+' h')



    g5 =ImageReader(output_path+'sensitivity_v1h1.png')
    canvas.drawImage(g5, 20, 20,width=270, height=270)

    g6 =ImageReader(output_path+'sensitivity_multiv.png')
    canvas.drawImage(g6, 280, 20,width=270, height=270)

    canvas.drawString(160, 290, 'Sensibility Analysis for the Two-layers Internal Wave Model:')



    canvas.showPage()
    canvas.setLineWidth(.3)
    canvas.setFont('Helvetica',10)

    canvas.drawString(520,800,'Page 4/4')

    g1 =ImageReader(output_path+'psd_multi.png')
    canvas.drawImage(g1, 20, 540,width=320, height=240)

    g1 =ImageReader(output_path+'meteo_spectra.png')
    canvas.drawImage(g1, 320, 492,width=240, height=200)

    g1 =ImageReader(output_path+'temperature_depth.png')
    canvas.drawImage(g1, 20, 380,width=320, height=160)

    if turn_iso == 1:
        canvas.drawString(40, 360, 'Isotherm analysis:')
        
        g1 =ImageReader(output_path+'isotherms.png')
        canvas.drawImage(g1, 10, 200,width=280, height=140)
  
        g1 =ImageReader(output_path+'iso_bandpass.png')
        canvas.drawImage(g1, 10, 50,width=280, height=140)

        g1 =ImageReader(output_path+'psd_hydro_coriois.png')
        canvas.drawImage(g1, 270, 50,width=300, height=225)      

  
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
    print ("> www.bit.ly/interwave_analyzer")
    root.update()
    print ("> ")
    root.update()
    print ("> ")
    root.update()
    root.update()
    dig.close() 
    
    root.mainloop()
    sys.stdout = old_stdout
