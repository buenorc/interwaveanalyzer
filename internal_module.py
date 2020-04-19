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


warnings.simplefilter("error")
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
    
    # Intergovernmental Oceanographic Commission (2010)
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

def depthmean(depth):
    
    hmean = np.zeros(len(depth[0]),float)
    for j in range(len(depth[0])):
        
        hmean[j] = np.mean(depth[:,j])
        
    hmean = hmean[::-1]
    return hmean

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
    
    if z == 10:
        return wz
       
    for t in range(l):
        
        if wz[t] < 5:
            Cd = 0.0010
        else:
            Cd = 0.0015
        
        wind = wz[t]/(1-np.sqrt(Cd)/k*exp) 

        w10[t] = wind
    
    return w10

def deriva(pu,pd,zu,zp):
    # function to calculate a discrete derivative
    der = abs(pu-pd)/abs(zu-zp)
    return der

def media(zu,zd):
    # function to measure the mean depth of two depths
    med = abs(zu+zd)/2
    return med

def sorting(data):
    
    ordered = np.zeros((len(data),len(data[0])),float)

    for t in range(len(data)):
        
        data_aux = vector_aux(data,t,len(data[0]))
        
        temporary = []
        first     = 'first'
        
        for i in range(len(data_aux)):
            
            if np.isnan(data_aux[i]) == False:
                temporary.append(data_aux[i])
            
                if first == 'first':
                    first = 'secondaries'
                    idx   = i
                
        temporary = -np.sort(-np.array(temporary))
    
        auxiliary     = np.full(len(data_aux),np.nan)
        auxiliary[idx:len(temporary)+idx] = temporary
        
        ordered[t][:] = auxiliary
    
    return ordered

def sorting_1d(data):

    temporary = []
    first     = 'first'
        
    for i in range(len(data)):
            
        if np.isnan(data[i]) == False:
            temporary.append(data[i])
            
            if first == 'first':
                first = 'secondaries'
                idx   = i
                
    temporary = -np.sort(-np.array(temporary))
    
    auxiliary     = np.full(len(data),np.nan)
    auxiliary[idx:len(temporary)+idx] = temporary
    
    return auxiliary


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

    distance = abs(numbers[0] - myNumber)

    idx = 0
    
    for c in range(len(numbers)):
            cdistance = abs(numbers[c] - myNumber)
            if(cdistance < distance):
                    idx = c
                    distance = cdistance           
    theNumber = numbers[idx]
    
    return theNumber, idx

def interpolation(y1,y2,x1,x2,x,exc,auxiso):
    
    if(y1 == y2 or x1 == x2):
        return media(y1,y2)
    
    else: 
        f = interpolate.interp1d([x1,x2], [y1,y2])    
        
        try:
            return f(x)
        
        except ValueError:
        
            return auxiso                       

def windbar_fetch(dw_mean,linang,angle,dists):
    
    Ldist = np.zeros(len(angle),float)
    
    dw_min, dw_max = wind_angle(dw_mean,linang)
    
    for i in range(len(angle)):
        if dw_max < dw_min:
            if angle[i] < dw_max or angle[i] > dw_min :
                Ldist[i] = dists[i]
            else:
                Ldist[i] = None
                
        else:
            if angle[i] > dw_min and angle[i] < dw_max:
                Ldist[i] = dists[i]
            else:
                Ldist[i] = None
    
    return Ldist

def wind_average(wd,iw):
   # wind direction average  

    u_east  = np.average(iw*np.sin(np.deg2rad(wd)))
    u_north = np.average(iw*np.cos(np.deg2rad(wd)))

    mean   = (np.arctan2(u_east, u_north)) * 180/np.pi
    mean  = (360+mean)%360

    
    return mean
    
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
    

    if(tau_near > reqt):
        if (idx == qt-1):           
            return None   # we cannot estimat the depth because 
                          # we dont have info bellow that point  
        else:
            return  interpolation(h[idx],h[idx+1],tau[idx],tau[idx+1],reqt,zmax,aux_iso)
    
    else:
        if (idx == 0) :
            return zmax
        else:
            return interpolation(h[idx-1],h[idx],tau[idx-1],tau[idx],reqt,zmax,aux_iso)
    

# -----------------------------------------------------------------------------
# Methods to identify the thermocline
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
    
    # interface between metalimnion and epilimnion
    z = zt
    dmin = deriva(p[z+1],p[z],h[z+1],h[z])
    
    while (dmin > minval) and (z > 0):                                            
        z=z-1
        dmin = deriva(p[z+1],p[z],h[z+1],h[z])

    try:
        aux1 = media(h[z+1],h[z]) - media(h[z+2],h[z+1])
        aux2 = deriva(p[z+1],p[z],h[z+1],h[z])-deriva(p[z+2],p[z+1],h[z+2],h[z+1])
    except IndexError:
        aux1 = 1
        aux2 = 1
    
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
# Average density and thickness of each layer (hypolimnion and epilimnion) 
# -----------------------------------------------------------------------------
    
def density_2layer (qt,h,tau,sal,pre,H,z0):
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
    
    hh = thermo - z0   # epilimnion thickness
    he = H - thermo     # hypolimnion thickness
    
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

    
        
    pe = pe/npe     # epilimnion average density       
    ph = ph/nph     # hipolimnion average density
    pu = p[0]       # superficial density 
    pd = p[qt-1]    # bottom density 
    
    return he, hh, pe, ph, pu, pd

def density_3layer (qt,h,tau,sal,pre,minval,H,z0):
    # average density and thickness of each layer
    # model with just 2 layers (epilimnion and hipolimnion)
    #
    # qt : quantidade de sensores (int)
    # h  : depth (1d vector)
    # p  : rho, temperature or salinity (1d vector)  
    # H  : maximum depth (depth of the reservoir) (float)
    #
    
    ze, zh = metathick (qt,h,tau,sal,pre,minval) 
    
    
    if(ze > zh and ze < h[0] and zh > z0):
        
        h1 = H  - ze
        h2 = ze - zh
        h3 = zh - z0
       
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
            p3 = p[-1]
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

def wind_angle(wind,wind_limit):
    
    dw_max = wind + wind_limit
    dw_min = wind - wind_limit

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
        Cd = 0.0015  # drag coefficient
    else:
        Cd = 0.0010
        
        if (w == 0.0): # low wind to characterize the stability of the lake
            w = 0.01
            
    
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
   
    wast    = math.sqrt(stress/pe) 
    s       = wast/(k*he)    
    
    try:
        riw = glin*he/((wast**2))
        frw  = 1/(riw**2)
        hast = wast**2/glin
        umal = wast**2/(2*he*math.sqrt(glin*(H-he)))
        
    except RuntimeWarning:
        riw,frw,hast,umal = None, None, None, None
    
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
    
    if wast == 0 or glin == 0:
        wedd = None
    else:
        wedd = glin*he**2/(ls*wast**2)
        
    return wedd    
# -----------------------------------------------------------------------------
# Parameters that dont use wind data 
# -----------------------------------------------------------------------------

def structure2layer(qt,h,tau,sal,pre,H,z0):
    
    # glin = Reduced gravity (m/s²)
    # N    = Bouyancy frequency  (1/s or Hz)

    he,hh,pe,ph,pu, pd = density_2layer(qt,h,tau,sal,pre,H,z0)
    
    glin   = abs(9.81*(ph-pe)/ph)  
       
    n    = math.sqrt(glin/he)
    
    return he,hh,pe,ph, glin, n,pu,pd

def structure3layer(qt, h, tau, sal, pre, minval, H, z0):
    
    h1,h2,h3,p1,p2,p3 = density_3layer(qt, h, tau, sal, pre, minval, H, z0)
    
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
        janela =  2*24*60*60 #seconds (window of 2 days)
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



def wavelet_spectral (si, dt, mother): 
    

    mom   = signal.ricker

    l=len(si)
    time  = np.arange(0,l) * dt
    per   = np.arange(1,48)
    
    power = signal.cwt(si,mom,per)
       
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
    

    sos = signal.butter(order,[low, high], analog=False, btype='band', output='sos')

    return sos




def butter_bandpass_filter(data, lowcut, highcut, fs):

    data   = data - np.mean(data)
    sos   = butter_bandpass( lowcut, highcut, fs, 4)
    f = signal.sosfilt(sos, data)
    
    return f


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


def spigel(he,W):
    
    return 0.4996*he/W

def bueno(he,hh,W):
    
    f = 0.76 - 0.52*he/(hh+he)
    k = 0.40
    
    return f*he/(1+np.exp(k*W)) 
    









    
    