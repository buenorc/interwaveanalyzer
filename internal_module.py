# -*- coding: utf-8 -*-
"""
Module with functions to calculate physical indices/parameter to identify baro-
clinic activity in thermal stratified system

@author: BUENO, R. (Rafael de Carvalho Bueno)
@date  : Jul 12 2017

"""
# pacotes 
from scipy   import signal, interpolate
from wavelib import wavelet
import numpy as np
import warnings
import math


# -----------------------------------------------------------------------------
# Function to define the comparison (coherence and phase shift)
# -----------------------------------------------------------------------------

def comparison_definition(c):

    comparison = ["None\n", "Iso. 1-2\n", "Iso. 1-3\n", "Iso. 1-4\n", "Iso. 2-3\n", "Iso. 2-4\n", "Iso. 3-4\n"]
        
    if c == comparison[0]:
        c = -999
    elif c == comparison[1]:
        c = 12
    elif c == comparison[2]:
        c = 13
    elif c == comparison[3]:
        c = 14
    elif c == comparison[4]:
        c = 23
    elif c == comparison[5]:
        c = 24
    elif c == comparison[6]:
        c = 34
    return c 
            

# -----------------------------------------------------------------------------
# Equations of state
# -----------------------------------------------------------------------------
  
def commission(tau,sal,pre):
    
    #Intergovernmental Oceanographic Commission (2010)
    # rho = f(T(°C),S(g/kg),P(dbar))
    #  2 < sal < 42
    #  0 < pre < 1,000 decibars

    c1 = 9.9984277040408688*10**2
    c2 = 7.3539907257802000*10**0
    c3 = -5.2725024846580537*10**-2
    c4 =  5.1051405427900501*10**-4
    c5 = 2.8372074954162994*10**0
    c6 = -5.7462873738668985*10**-3
    c7 = 2.0165828404011005*10**-3 
    c8 = 1.1506680128760695*10**-2
    c9 = 1.2026027029004581*10**-7
    c10 = 5.5361909365048466*10**-6
    c11 = -2.7563156404651928*10**-8
    c12 = -5.8834769459933364*10**-12
    c14 = 7.2882773179945397*10**-3
    c15 = -4.4270423575705795*10**-5
    c16 = 4.8218167574165732*10**-7
    c17 = 1.9666437776499541*10**-10
    c18 = 2.0192201315731156*10**-3
    c19 = -78386667410747600*10**-6
    c20 = -2.7493971171215844*10**-10
    c21 = 4.6614190290164293*10**-6
    c22 = 1.5182712637288295*10**-9
    c23 = 6.4146293567422886*10**-6
    c24 = -9.5362845886397360*10**-17
    c25 = -9.6237455486277320*10**-18


    a0 = 1 + c14*tau + c15*tau**2 + c16*tau**3 + c17*tau**4 + c18*sal + \
    c19*sal*tau + c20*sal*tau**3 + c21*sal**1.5 + c22*tau**2*sal**1.5
    
    a1 = c23
    a2 = c24*tau**3
    a3 = c25*tau
    
    b0 = c1 + c2*tau + c3*tau**2 + c4*tau**3 + c5*sal + c6*sal*tau + c7*sal**2
    b1 = 0.5*(c8 + c9*tau**2 + c10*sal)
    b2 = c11 + c12*tau**2

    if (pre == 0) :
        v = a0/b0
    else :
        v  = (a0 + a1*pre + a2*pre**2 + a3*pre**3) / \
        (b0 + 2*b1*pre + b2*pre**2)
     
    rho = 1/v
    
    return rho

def density_state(qt,tau,sal,pre):  # water density profile
    
    rho   = np.zeros(qt,float)
    
    for z in range(qt):
        rho[z] = commission(tau[z],sal,pre)
        
    return rho
# -----------------------------------------------------------------------------
# Conversion functions
# -----------------------------------------------------------------------------

def vector_aux (a,i,nj):
    
    b =  np.zeros((nj),float)
    
    for j in range(nj):
        b[j] = a[i][j]
    return b

def vector_time (a,j,ni):
    
    b =  np.zeros((ni),float)
    
    for i in range(ni):
        b[i] = a[i][j]
    return b

def keepout (angle):
    if angle < 90:
        angle = angle + 180 
    else :
        if angle > 270 : 
            angle = angle - 180 
    
    return angle

def ciout(x):
    
    # just the the dev from ci(x)
    nivel = ci(x)[1]-np.mean(ci(x))
    return nivel

def ci(x):
    # 95% -> 1.96 (alpha)
    size    = numberOfNonNans(x)
    sdev    = np.nanstd(x) 
    
    val     = sdev/np.sqrt(size)
    z       = 1.96
    
    dev = val*z
    
    try:
        upper = np.nanmean(x) + dev 
        lower = np.nanmean(x) - dev 
    except RuntimeWarning:
        upper = None
        lower = None
    
    return lower,upper

def numberOfNonNans(data):
    count = 0
    for i in data:
        if not np.isnan(i):
            count += 1
    return count 

def find_nearest(array,value):
    # find the nearest value of a 'value' from a 'array'
    #
    idx = (np.abs(array-value)).argmin()
    return idx


def velocityten(wz,z):
    # convert wind velocity given at z meters to 10 meters
    # wz = wind velocity at z meters 
    
    k   = 0.4 #Von Karman's constant
    exp = math.log(10/z)
    
    l = len(wz)
    w10  = np.zeros(l,float)
       
    for t in range(l):
        
       Cd = 1.3*10**-3
       wind  = wz[t]*(1-math.sqrt(Cd)*exp/k)**-1
        
       if wind >= 5.0 :
            Cd = 1.5*10**-3 
            w10[t]  = wz[t]*(1-math.sqrt(Cd)*exp/k)**-1
       else:
            Cd  = 1.0*10**-3 #drag coefficient
            w10[t]  = wz[t]*(1-math.sqrt(Cd)*exp/k)**-1
    
    return w10

def deriva(pu,pd,zu,zp):
    # function to calculate a discrete derivative
    der = abs(pu-pd)/abs(zu-zp)
    return der

def media(zu,zd):
    # function to measure the mean depth of two depths
    med = abs(zu+zd)/2
    return med

#def mean(d):
    # difference between data and the mean data 
#    mean     = np.mean(d)
 #   diff     = (d - mean)
    
#    return diff

def zero():
    
    return 0,0,0,0

def average(arr,n):
    # n average of arr (1d array)
    n=int(n)
    end = n*int(len(arr)/n)
    return np.nanmean(arr[:end].reshape(-1,n),1)

def value_nearest(array,value):
    # find the nearest value of a 'value' from a 'array'

    idx = (np.abs(array-value)).argmin()
    return array[idx], idx


def near (numbers,myNumber):

    distance = abs(numbers[1] - myNumber)

    idx = 1   
    c=1
    while(c < len(numbers)):
            cdistance = abs(numbers[c] - myNumber)
            if(cdistance < distance):
                    idx = c
                    distance = cdistance           
            c=c+1
    theNumber = numbers[idx]
    
    return theNumber, idx

def interpolation(y1,y2,x1,x2,x):
    
    if(y1 == y2 or x1 == x2):
        return media(y1,y2)
    
    else:
        f = interpolate.interp1d([x1,x2], [y1,y2])
        ynew = f(x)
        return ynew


def wind_average(wd,ws):
   # wind direction average  
    V_east = np.average(np.sin(wd * np.pi/180))
    V_north = np.average(np.cos(wd * np.pi/180))
    mean_wd = np.arctan2(V_east, V_north) * 180/np.pi
    mean_wd = (360 + mean_wd) % 360
    
    return mean_wd
    
def depths(profile, depth, lin):

    d   = []
    aux = np.zeros((lin),float)

    if profile != -999:
        for t in range(lin):
            aux[t] = depth[t][profile]
    else:
        aux = None
        
    d.append(aux)
    return d

def interpol(y1,y2,x1,x2,x):
    
    if(x1==x2):
        return media(y1,y2)
    
    a = (y1-y2)/(x1-x2)
    b = y1 - a*x1
    
    y = a*x + b
    return y
# -----------------------------------------------------------------------------
# Methods to identify thermo
# -----------------------------------------------------------------------------

def isotherms (reqt, qt, h,tau,zmax,zmin,aux_iso):
    
    tau_near, idx = near(tau,reqt)
    
    #if (aux_iso != None and abs(aux_iso-h[idx]) > 4):
    #    return aux_iso

    if(tau_near > reqt):
        if (idx == qt-1):           
            return zmin
        else:
            if (tau[idx] < tau[idx+1]):
                if aux_iso != None:
                    return aux_iso                
                else:
                    return None
  
            else:
                return  interpolation(h[idx],h[idx+1],tau[idx],tau[idx+1],reqt)

    
    else:
        if (idx == 1) :
            return zmax
        else:
            if(tau[idx] > tau[idx-1]):
                return zmax    
            else:
                return interpolation(h[idx-1],h[idx],tau[idx-1],tau[idx],reqt)
    
 
                
# -----------------------------------------------------------------------------
# Methods to identify Thermocline/Halocline/Pycnocline
# -----------------------------------------------------------------------------

def thermo_region (qt,h,tau,sal,pre):
    # identification of the region of the thermocline/pycnocline
    # ps: the thermocline is placed between which sensors  
    #
    # qt : quantidade de sensores (int)
    # h  : depth (1d vector)
    # p  : rho, temperature or salinity (1d vector) 
    
    p = density_state(qt,tau,sal,pre)
    
    z,zt   = 1,1
    dpz = deriva(p[z+1],p[z],h[z+1],h[z])
    
    for z in range(1,qt-2) :   # we are not considering that thermocline could not be 
                             # located between the 6th and 7th sensor.
        newdpz = deriva(p[z+1],p[z],h[z+1],h[z])
        
        if newdpz >= dpz :
            dpz = newdpz
            zt  = z
    
    mid  = media(h[zt],h[zt+1])   # simp. prof da termoclina 
    
    # return: dpz : greater difference
    #         zt  : índice relacionado ao sensor superior do maior gradiente
    #               (1st sensor has indece iqual 0) 
    #         mid : simple thermocline depth 
    #         p   : presurre (vector of the z data) 
    return dpz, zt, mid, p
    
def thermocline_depth (qt,h,tau,sal,pre):
    # identification of the depth of the thermocline/pycnocline
    # ps: the exacly location of the thermocline  
    #
    # qt : quantidade de sensores (int)
    # h  : depth (1d vector)
    # tau  : rho, temperature or salinity (1d vector)  
    #
    # this function use the thermo_region fuction
    dpz, z, mid, p = thermo_region(qt,h,tau,sal,pre)
    
    if(z == 0 or z > qt-3): # Return the normal average when the thermocline is 
        return mid, p       # located between the last or first sensor.           
    
    dpup = deriva(p[z+2],p[z+1],h[z+2],h[z+1])
    dpdn = deriva(p[z],p[z-1],h[z],h[z-1])
    
    tdup = abs(dpz - dpup)
    tddn = abs(dpz - dpdn)
    
    if(tddn==0) :
        tddn = 1
    if(tdup==0) :
        tdup = 0.001
    
    pup  = abs(media(h[z+2],h[z+1])-mid)/tdup
    pdn  = abs(mid-media(h[z],h[z-1]))/tddn    
    pdif = pdn+pup
    
    thermo = h[z+1]*(pup/pdif) + h[z]*(pdn/pdif)
    
    

    
    # return: thermo : thermocline deapth (relation with the surface water)
    #         p      : presurre (vector of the z data) 
    return thermo, p

# -----------------------------------------------------------------------------
# Thickness of the metalimnion
# -----------------------------------------------------------------------------

def metathick (qt,h,tau,sal,pre,minval):
    # function to calculate the bottom and the top of the metalimnion
    # (thickness of the metalimnion layer)
    #     
    # minval: minimum value for gradienty density considering  
    #         the metalimnion interface (float) (kg/m³/m)
    # 
    #
    dpz, zt, mid, p = thermo_region(qt,h,tau,sal,pre)
    
    # interface between metalimnion and heplimnion
    z = zt
    dmin = deriva(p[z+1],p[z],h[z+1],h[z])
    
    while (dmin > minval) and (z > 0):                                            
        z=z-1
        dmin = deriva(p[z+1],p[z],h[z+1],h[z])

    
    aux1 = media(h[z+1],h[z]) - media(h[z+2],h[z+1])
    aux2 = deriva(p[z+1],p[z],h[z+1],h[z])-deriva(p[z+2],p[z+1],h[z+2],h[z+1])
    
    if(aux2 == 0.0):
        ze = media(h[z+1],h[z])
    
    else:
        ze   = media(h[z+1],h[z]) + (minval - dmin)#*aux1/aux2
    
    
    # interface between metalimnion and hypolimnion
    z = zt
    dmin = deriva(p[z+1],p[z],h[z+1],h[z])                                        
    
    while (dmin > minval) and (z < qt-2):    
        z=z+1                                                                     
        dmin = deriva(p[z+1],p[z],h[z+1],h[z])                                    

    
    
    aux1 = media(h[z+1],h[z]) - media(h[z],h[z-1]) 
    aux2 = deriva(p[z+1],p[z],h[z+1],h[z])-deriva(p[z],p[z-1],h[z],h[z-1])
    
    if (aux2 == 0.0):
        zh = media(h[z+1],h[z])                                                 
    else:
        zh   = media(h[z+1],h[z]) + (minval - dmin)*aux1/aux2                 
    
    # return:  
    #    ze: depth of the top of the metalimnion (m in relation with 800)
    #    zh: depth of the botton of the metalimnion (m in relation with 800)      
    return ze,zh
        

# -----------------------------------------------------------------------------
# three layer system - based on thickness of the thermocline
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Average density and thickness of each layer (hpolimnion and epilimnion) 
# -----------------------------------------------------------------------------
    
def density_2layer (qt,h,tau,sal,pre,H):
    # average density and thickness of each layer
    # model with just 2 layers (epilimnion and hipolimnion)
    #
    # qt : quantidade de sensores (int)
    # h  : depth (1d vector)
    # p  : rho, temperature or salinity (1d vector)  
    # H  : maximum depth (depth of the reservoir) (float)
    #
    # this function use the thermocline_depth  
    
    thermo, p = thermocline_depth (qt,h,tau,sal,pre)
    
    hh = thermo - 799   # epilimnion thickness
    he = H - thermo   # hipolimnion thickness
    
    pe,ph,npe,nph=0,0,0,0

    for z in range(qt) :
        if(thermo < h[z]):
            pe   =  pe + p[z]
            npe  =  npe +  1
        else:
            ph   =  ph + p[z]
            nph  = nph +  1

    if npe == 0:
        npe = 1
        pe  = p[0]
        print('Warning 187: epilimnion thickness reached zero for some reason')
    pe = pe/npe     # epilimnion average density       
    ph = ph/nph     # hipolimnion average density
    pu = p[0]       # superficial density (first sensor = 1m from the top)
    pd = p[qt-1]    # bottom density (last sensor = 2m from the bottom)

    return he, hh, pe, ph, pu, pd

def density_3layer (qt,h,tau,sal,pre,minval,H):
    # average density and thickness of each layer
    # model with just 2 layers (epilimnion and hipolimnion)
    #
    # qt : quantidade de sensores (int)
    # h  : depth (1d vector)
    # p  : rho, temperature or salinity (1d vector)  
    # H  : maximum depth (depth of the reservoir) (float)
    #
    
    ze, zh = metathick (qt,h,tau,sal,pre,minval) 
    
    
    if(ze > zh and ze < h[0] and zh > 800):
        
        h1 = H  - ze
        h2 = ze - zh
        h3 = zh - 799
       
        p1, p2, p3, np1, np2, np3 = 0, 0, 0, 0, 0, 0
       
        p = density_state(qt,tau,sal,pre)
            
        for z in range(qt):
            if(ze < h[z]):
               p1   =  p1 + p[z]
               np1  =  np1 +  1
            else:
                
                if(zh < h[z]):
                   p2   =  p2 + p[z]
                   np2  =  np2 +  1
                
                else:
                   p3   =  p3 + p[z]
                   np3  =  np3 +  1           
                   
        p1 = p1/np1

        if np3 == 0:
            p3 = p2
        else:
            p3 = p3/np3
        
        if np2 == 0:
            p2 = (p1+p3)/2
        else:
            p2 = p2/np2
    else:
        
        h1,h2,h3,p1,p2,p3 = 999, 999, 999, 999, 999, 999
    
    return h1,h2,h3,p1,p2,p3
        
        
        
# -----------------------------------------------------------------------------
# Parameters that use wind data 
# -----------------------------------------------------------------------------

def wind_angle(wind_mean,wind_limit):
    
    dw_max = wind_mean + wind_limit
    dw_min = wind_mean - wind_limit

    if (dw_max > 360):
        dw_max = dw_max - 360

    if (dw_min < 0):
        dw_min = 360 + dw_min 
    
    return dw_min, dw_max

def wind_stress(w):
    # wind stress (N/m² or Pa)
    # w: wind velocity at 10 meters above the water surface
    # rw: data of the wind was obtained at rw meters  
    #
    
    if(w > 5):
        Cd = 1.5*10**-3  # drag coefficient
    else:
        Cd = 1.0*10**-3
        if (w == 0.0):
            return -999
            
    
    rho_air = 1.225 # density of the air (kg/m³)
    stress = Cd*rho_air*w**2
    
    return stress

def wind_parameters(w,rw,pe,he,n,glin,H):
    
    # stress = Wind stress (N/m²) 
    # wast   = Friction velocity of the wind  (m/s)
    # s      = Instability frequency caused by the wind (Hz)
    # riw    = Richardson number (non-dim)
    # frw    = Froude number (non-dim) 
    # hast   = The vertical scale of turbulent fluctuation at the pycnocline(m)
    # umal   = maximum shear velocity generated across the 
    #          base of the mixed layer per meter of reservoir length. (m/(s.m))
    # wedd   = Wedderburn number (similiar to Richardson with more he/l) 
    
    k       = 0.4                    #Von Karman's constant 
    stress  = wind_stress(w)    
    
    if stress < 0:
        return None, None, None, None, None, None, None
    
    wast    = math.sqrt(stress/pe) 
    s       = wast/(k*he)    
    

    riw = glin*he/((wast**2))
    frw  = 1/(riw**2)
    hast = wast**2/glin
    umal = wast**2/(2*he*math.sqrt(glin*(H-he)))
    
    return stress, wast, s, riw, frw, hast, umal

def wind_para2d(w,rw,qt,h,pe,H,n2d,hmid,p,glin):
    
    # Function to generate parameters for each point of the grid.
    # Usefull to generate a countourf graph 
    
    # glin = Reduced gravity (m/s²)
    # wast = Friction velocity of the wind  (m/s)
    # s    = Instability frequency caused by the wind (Hz)
    # n    = Bouyancy frequency  (1/s or Hz)
    # riw  = Richardson number (non-dim)

    s2d    = np.zeros(qt-1,float)
    riw2d  = np.zeros(qt-1,float)

    
    for z in range(qt-1):

        win_stress = wind_stress(w)
        
        if win_stress < 0:
            riw2d[z] = None
            s2d[z]   = None
        else:
            k       = 0.4   #Von Karman's constant 
            wast      = math.sqrt(win_stress/media(p[z+1],p[z]))
            s2d[z]    = wast/(k*hmid[z]) 
            riw2d[z] =  (glin[z]*( H - hmid[z]))/((wast**2))                  
    
    return s2d, riw2d


def glinha_2d (p,qt):

    glinha  = np.zeros(qt-1,float)
    
    for z in range(qt-1):

        delta_rho   = p[z+1] - p[z]
        glinha[z]   = 9.81*delta_rho/p[z+1]   
        
    return glinha    
    
    
    
    
def wedderburn(glin,he,wast,ls):
    # Wedderburn Number - improvemenet of Ri (non-dimensional)
    #
    # glin   = reduced gravity (m/s²) 
    # wast   = friction velocity of the wind  (m/s)
    # he     = Depth/thickness of the Epilimnion end (m)
    # ls     = system fetch length (m)  
    
    if wast == 0:
        wedd = None
    else:
        wedd = glin*he**2/(ls*wast**2)
        
    return wedd    
# -----------------------------------------------------------------------------
# Parameters that dont use wind data 
# -----------------------------------------------------------------------------

def structure2layer(qt,h,tau,sal,pre,H):
    
    # glin = Reduced gravity (m/s²)
    # N    = Bouyancy frequency  (1/s or Hz)
    
    
    he,hh,pe,ph,pu, pd = density_2layer(qt,h,tau,sal,pre,H)
    glin   = abs(9.81*(ph-pe)/ph) 
    n    = math.sqrt(glin/he)    
    
    
    return he,hh,pe,ph, glin, n,pu,pd

def structure3layer(qt, h, tau, sal, pre, minval, H):
    
    h1,h2,h3,p1,p2,p3 = density_3layer(qt, h, tau, sal, pre, minval, H)
    
    return h1,h2,h3,p1,p2,p3

def nonwind_para2d(qt,h,H,tau,sal,pre):
    
    # Function to generate Buoyancy freq. for each point of the grid.
    
    # glin = [float] Reduced gravity (m/s²)
    
    n2d    = np.zeros(qt-1,float)
    p      = np.zeros(qt,float)
    hmid   = np.zeros(qt-1,float)
    
    for z in range(qt):
        p[z] = commission(tau[z],sal,pre)
    
    for z in range(qt-1):
        hmid[z]   =  media(h[z],h[z+1]) 
    
        glin      = abs(9.81*(p[z+1]-p[z])/p[z+1])
        if(h[z] != h[z+1]):
            n2d[z]    = math.sqrt(glin/(h[z]-h[z+1]))
        else:
            n2d[z]    = math.sqrt(glin/0.01)
    
    # output:
    #
    # p    = [vec1d] Water density (kg/m³)
    # n2d  = [vec1d] Bouyancy frequency in z (1/s or Hz)
    # hmid = [vec1d] mid point of each point of the grid at 800m reference.
    
    return p, n2d, hmid

# -----------------------------------------------------------------------------
# Frequency scales to degeneration of internal seiches 
# -----------------------------------------------------------------------------

def class_degeneration (fi,sy,ab,V,hh,he,delp,umax,ls,no,co,glin):
    #
    # function to classify the degeneration of internal waves (theoretical)
    # this function uses:
    #   - dampingfreq    : frequency for damping    
    #   - steepiningfreq : frequency to generate non linear internal waves
    #   - kelvinfreq     : frequency to generate K-H billowing 
    #   - boresfreq      : frequency to geberate internal bores
    #    
    # input:  
    #        fi - internal wave frequency considering a two layer system (Hz)   
    #        sy - system (0 to natural fild and 1 to lab experiment)
    #        ab - area of the boundary layer (m²)
    #        V  - volume of the reservoir (m³)
    #        he and hh - thickness of the epi- and hypo-limnion (m)
    #        delp - thickness of the thermocline (m)
    #        umax - maximum shear velocity (m/s)
    #        ls - system fetch length (m)
    #        no - initial amplitude of the internal seiche (m)
    #        co - velocity of the internal wave (m/s)
    #        glin - reduced gravity (m/s²) 


     
    td = 1/dampingfreq(fi,sy,ab,V,hh,he,delp,umax,ls)
    ts = 1/steepingfreq(fi,he,hh,no,co,ls)
    tk = 1/kelvinfreq(glin,no,delp,ls)
    tb = 1/boresfreq(no,fi,he,hh) # mudar
    
    
    if tb < 1/(4*fi):
        return 1 # "internal bore will form."
    else:
        if td < ts :
            if tk < 1/(4*fi) and tk < ts:    
                return 2 # "Kelvin-Helmholtz billows will form."
            else:
                return 3 #"Internal seiches will be damped normally."
        else:
            return 4 #"Non-linear internal wave will form."
    
def boresfreq(no,fi,he,hh):
    # 
    # function to calculate the internal bores frequency (Hz)
    # important: before to use this function, you need to calculate the 
    #            internal wave frequency usng a two layers method.
    #

    aux = 4*no*fi/he
    fb = aux * math.sqrt((math.pow(he,3)+math.pow(hh,3))/\
                         ((hh+he)*math.pow(hh,2)))
    return fb
    
def dampingfreq(fi,sy,ab,V,hh,he,delp,umax,ls):
    # 
    # function to calculate the viscous damping frequency (Hz)
    # important: before to use this function, you need to calculate the 
    #            internal wave frequency usng a two layers method.
    #
    # input:  fi - internal wave frequency considering a two layer system (Hz) 
    #         sy - system (0 to natural fild and 1 to lab experiment)
    #         ab - area of the boundary layer (m²)
    #         V  - volume of the reservoir (m³)
    #         he and hh - thickness of the epi- and hypo-limnion (m)
    #         delp - thickness of the thermocline (m)
    #         umax - maximum shear velocity (m/s)
    #         ls - system fetch length (m)
    #
    # output  fd - frequency per lenght unit (Hz)
    #
    # variables: 
    #         nu - kinematic viscosity of the water at 20°C (m²/s)
    #         H  - total depth - hh + he (m)   
    #
    nu = 1.00*10**-6 
    
    H = hh + he
    
    if sy == 0 :
        # for fild analysis, delp, hh, he is nor required.
        delb   = umax*math.e/(471*math.sqrt(nu*fi))
        alphad = ab*delb/(2*V)
    else:
        # for lab analysis, umax is not required.
        delb   = math.sqrt(nu/(math.pi()*fi))
        alphad = math.pi()*ab*delb/(V) + nu*H/(fi*delp*he*hh) 
    
    
    fd= alphad*fi*ls
    
    return fd

def steepingfreq(fi,he,hh,no,co,ls):
    # 
    # function to calculate the nonlinear steeping frequency-scale (Hz)
    # important: before to use this function, you need to calculate the 
    #            internal wave frequency usng a two layers method.
    #
    # input:  fi - internal wave frequency considering a two layer system (Hz) 
    #         he and hh - thickness of the epi- and hypo-limnion (m)
    #         no - initial amplitude of the internal seiche (m)
    #         co - velocity of the internal wave (m/s)
    #         ls - system fetch length (m)
    #
    # output  fs - frequency times lenght unit (Hz)
    #  
    #
    
    beta = 3*co*(he-hh)/(2*he*hh)
    fs   = beta*no/ls
    
    return fs

def kelvinfreq(glin,no,delp,ls):
    # 
    # function to calculate the Kelvin_helmholtz frequency-scale (Hz)
    # important: before to use this function, you need to calculate the 
    #            internal wave frequency usng a two layers method.
    #
    # input:  glin - reduced gravity (m/s²) 
    #         no - initial amplitude of the internal seiche (m)
    #         delp - thickness of the metalimnion (m)
    #
    # output  fkh - frequency times lenght unit (Hz)
    #
    
    fkh = no*math.sqrt(glin/delp)/ls
    
    return fkh

# -----------------------------------------------------------------------------
# Generation of internal seiches 
# -----------------------------------------------------------------------------
   
def class_generation (riw,hh,he,ls):
    
    # generation of internal seiches according to the Richardson Number
    #
    # outputs: classification along time series [array int 1d] = 1, 2, 3 or 4
    #
    #    1 = Stratification is broken down by mixing
    #    2 = Large interface displacement occurs and is accompanied by 
    #        interface shear and Kelvin-Helmholtz billows (strat is broken)
    #    3 = Internal seiche is dominant
    #    4 = Internal seiche with short amplitude.
    #

    aux1 = math.sqrt((hh+he)/he)
    
    if riw < 1:
        return 1
    else:
        aux2 = ls/(2*he)
        if riw < aux2*aux1:
            return 2
        else:
            if riw < aux1*aux2**2:
                return 3
            else:
                return 4

def iw_generation (riw, hh, he, ls):
    #
    # function to verify the generation of internal waves
    # the lower_gene indicates the lowest Ri for the generation of iw 
    #
    #
    
    aux1 = math.sqrt((hh+he)/hh)
    
    lower_gene = ls/(2*he)*aux1
    upper_gene = ls**2/(4*he**2)*aux1**2
    
    return lower_gene, upper_gene
    
# -----------------------------------------------------------------------------
# Spectral Analysis
# -----------------------------------------------------------------------------    
    


def mask (serie):
    
    xi = np.arange(len(serie))
    
    m = np.isfinite(serie)
    xfiltered = np.interp(xi, xi[m], serie[m])
    
    return xfiltered

def welch_method(serie,windsize,window, dt):
    
    if (windsize == 0):
        janela =  72*60*60 #seconds (window of 3 days = 72hours)
    else:
        janela = windsize 
    

    dt     = 60*60*dt  #seconds
    nperseg = int(janela/dt)
    
    per   = []
    
        
    serie = mask(serie) 
    freq, ff = signal.welch(serie, fs=1.0/dt, window=window, \
                            nperseg=nperseg, detrend='linear',\
                            return_onesided=True, scaling='density',axis=-1)
    
    for x in freq:
        per.append(1.0/(x*60*60))    #hours
        
    return freq, per, ff


def wave_spectral (si, dt, mother): 
    
    # ---------------------  wavelet parameters -------------------------------
 
    pad = 1            # pad the time series with zeroes (recommended)
    dj = 0.10          # this will do 10 sub-octaves per octave
    s0 = 4.*dt         # this says start at a scale of 1 hour
    j1 = 7./dj         # this says do 7 powers-of-two with dj sub-octaves each
   
    # -------------------------------------------------------------------------
      
    si = mask(si)
    l=len(si)
    time  = np.arange(0,l) * dt
    
    variance = np.std(si)**2
    mean     = np.mean(si)
    data     = (si - mean)/np.sqrt(variance)

    warnings.filterwarnings("ignore")

    wav,per,sca,coi = wavelet(data,dt,pad,dj,s0,j1,mother)
    power           = (np.abs(wav))**2   
       
    return time, per, power




def coherence_shift(s1,s2,nfft,dt):

    # coherence

    fcoh, cxy = signal.coherence(s1, s2,  1/dt)
    px, pxy = signal.csd(s1, s2, 1/dt)

    phase = np.angle(pxy)
    phase[phase<-90] += 360
    
    #calculate 95 % confidence level
    edof = (len(s1)/(nfft/2)) * cxy.mean() 
    gamma95 = 1.-(0.01)**(1./(edof-1.))

    conf95 = np.where(cxy>gamma95)
    
    return phase, cxy , fcoh, conf95



def butter_bandpass(lowcut, highcut, fs, order):
    
    nyq = 0.5*fs
    
    low  = lowcut/nyq
    high = highcut / nyq
    

    b, a = signal.butter(order,[low, high], btype='band')

    return b, a




def butter_bandpass_filter(data, lowcut, highcut, fs,  order):

    
    b, a   = butter_bandpass( lowcut, highcut, fs, order=order)
    f = signal.lfilter(b, a, data)
    
    return f


def subclass(iw):
    
    
    if(iw<0.080):
        iwi = 'Small amplitude internal wave'
    elif(iw<0.333):
        iwi = 'Sigmoidal internal wave'
    elif(iw<1.25):
        iwi = 'Non linear internal wave'
    else:
        iwi = 'Strong mixing caused by internal wave'
        
    return iwi

def amplitude_average(iw,he,wedd):
    
    count = 0
    summ  = 0
    lim   = 1/wedd
    for i in range(len(iw)):
        if iw[i] > lim:
            count = count + 1 
            summ  = summ  + iw[i]
    if count == 0:
        average= -999
    else:
        average = (summ/count)*he
    ratio   = count/len(iw)
    
    return average, ratio


    









    
    