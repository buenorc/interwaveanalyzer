# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno

Module with functions to calculate physical indices/parameter to identify 
baroclinic activity in thermal stratified system (brain module)

Internal Functions:  functions that is used just by the software (does not provide any additional information)
External Functions:  functions that can be used externally 

modififications to previosly versions must be specified here using WAVE codes:
    
W-23.04-1.00.0-00
A-01.04-1.00.3-00
V-22.04-1.00.3-00
E-05.04-1.00.3-00

"""

# external modules
import math
import warnings
import numpy.matlib

import numpy as np
import scipy.special as sc

from scipy.stats import t, sem
from scipy import signal, interpolate
from scipy.stats.distributions import chi2

# internal modules
from   wavelib import wavelet
import info_warning as war 


# supress warnings
warnings.simplefilter("error")


def wind_fetch(dw_mean,linang,angle,dists,type_length,dig):
    
    if type_length == 2:
        Lw_cont  = windbar_fetch(dw_mean,linang,angle,dists)
        
        try:
            ls_min = np.nanmin(Lw_cont)
            ls_max = np.nanmax(Lw_cont)
            ls_ave = np.nanmean(Lw_cont)
            ls_fetch = [ls_min,ls_ave,ls_max]  
            
            if np.isnan(Lw_cont).all() == True:
                ls_near  = dists[find_nearest(angle,dw_mean)[1]]
                ls_fetch = [0.95*ls_near,ls_near,1.05*ls_near]
                
                war.coarse_fetch(dig)
            
        except :
            ls_near  = dists[find_nearest(angle,dw_mean)[1]]
            ls_fetch = [0.95*ls_near,ls_near,1.05*ls_near]
            
            war.coarse_fetch(dig) 
    
    else:
        ls_fetch = [0.95*dists, dists, 1.05*dists]
    
    return ls_fetch
 
def comparison_definition(c):
#   
#   Internal Function: Code Definition
#   Function used to identify which pair of isothers will be analyzed 
#       
#   c (input)    label defined in GUI definition 
#   c (output)   label translated to iso map code (define isotherms pairs) 

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

                 
def depths(profile, depth, lin):
#
#   Internal Function: Geometry definitions
#   Function to get the depth that are associated to specified/analyzed sensors
#
    d   = []
    aux = np.zeros((lin),float)

    if profile != -999:
        aux[:] = depth[:,profile]
    else:
        aux = None
        
    d.append(aux)
    return d
          

def commission(tau):
#
#   External Function: Thermal Stratification 
#   Equation of state used to compute the water density from tempearture data
#    
#   Craig (1961) and Millero (1976) 
#   rho(kg/m³) = f(T(°C))

    p0 = 999.842592
    p1 = 6.793952*10**-2
    p2 = 9.095290*10**-3
    p3 = 1.001685*10**-4
    p4 = 1.120083*10**-6
    p5 = 6.536336*10**-9

    rho = p0 + p1*tau - p2*tau**2 + p3*tau**3 - p4*tau**4 + p5*tau**5
    
    return rho


def isotherms (reqt, qt, h,tau,zmax,zmin,aux_iso):
#
#   External Function: Thermal Stratification Analysis
#   Function to define the specified isotherms that will be analyzed
#
    
    tau_near, idx = find_nearest(tau,reqt)


    if(tau_near > reqt):
        if (idx == qt-1):           
            return None   # NaN: When depth cannot be estimated  
                          #      the isotherm is bellow the lake bottom   
        else:
            return  interpolation(h[idx],h[idx+1],tau[idx],tau[idx+1],reqt,aux_iso)
    
    else:
        if (idx == 0) :
            return None   # When the isotherm is not defined (water surface)
        else:
            return interpolation(h[idx-1],h[idx],tau[idx-1],tau[idx],reqt,aux_iso)
    

def thermo_region (qt,h,tau):
#
#   External Function: Thermal Stratification Analysis    
#   Function to define the thermocline depth (standard detemination) 
#
#   qt      input       number of sensors (integer)
#   h       input       depth (1d vector)
#   tau     input       temperature (1d vector)
#   dpz     output      highest gradient difference
#   zt      output      sensor index (1st sensor has index zero)
#   mid     output      simple thermocline depth
#   rho     output      water density (1d vector z-length)
    
    rho = commission(tau)
    
    z,zt   = 1,1
    dpz = deriva(rho[z+1],rho[z],h[z+1],h[z])
    
    for z in range(0,qt-1) :   
            
        new_dpz = deriva(rho[z+1],rho[z],h[z+1],h[z])
        
        if new_dpz > dpz :
            dpz = new_dpz
            zt  = z
    
    mid  = abs(np.mean([h[zt],h[zt+1]]))   
    
    return dpz, zt, mid, rho


    
def thermocline_depth (qt,h,tau):
#
#   External Function: Thermal Stratification Analysis
#   Function to define the thermocline depth (weighted method)
#
#   qt      input       number of sensors (integer)
#   h       input       depth (1d vector)
#   tau     input       temperature (1d vector)
#   thermo  output      thermocline depth (referenced to the water surface)
#   rho     output      water density (1d vector z-length)
    

    error = 0
    dpz, z, mid, rho = thermo_region(qt,h,tau)
    
    
    if z == 0:
        error  = 1
        thermo = mid
    else:
        try:
            
            hplus = (h[z]   - h[z+2])/2
            hminu = (h[z-1] - h[z+1])/2
            
            drho      = (rho[z+1]-rho[z])/(h[z]-h[z+1])
            drho_plus = (rho[z+2]-rho[z+1])/(h[z+1]-h[z+2])
            drho_minu = (rho[z]-rho[z-1])/(h[z-1]-h[z])

            Dplus = hplus/(drho-drho_plus)
            Dminu = hminu/(drho-drho_minu)
   
            thermo = h[z+1]*(Dplus/(Dminu+Dplus)) + h[z]*(Dminu/(Dminu+Dplus))
        
        except:
            error  = 1
            thermo = mid
  

    return thermo, rho, error


def metathick (qt,h,tau,minval,z0):
#
#   External Function: Thermal Stratification Analysis
#   Function to define the metalimnion boundaries (metalimnion thickness)
#     
#   qt      input       number of sensors (integer)
#   h       input       depth (1d vector)
#   tau     input       temperature (1d vector)
#   minval  input       the metalimnion treshold unit kg/m³/m (float number) 
#   ze      output      upper boundary - epilimnion  layer (m from ref. level)
#   zh      output      lower boundary - hypolimnion layer (m from ref. level)
    
    dpz, zt, mid, rho = thermo_region(qt,h,tau)
    
    # 1.1) metalimnion-epilimnion interface:   
    
    z = zt   
                           
    dmin = (rho[z+1] - rho[z])/(h[z] - h[z+1])
    
    try:
        while (dmin > minval) and (z > 0) and (z < qt-2):                                            
            z=z-1
            dmin = (rho[z+1] - rho[z])/(h[z] - h[z+1])

        aux1 =  (h[z]-h[z+2])/2  # now this term is negative as is defined in the paper, but i am not sure !!! (changed to positive)
        rhodown = (rho[z+1]-rho[z])/(h[z] - h[z+1])
        rhoup   = (rho[z+2]-rho[z+1])/(h[z+1] - h[z+2])
        
        aux2 = rhodown - rhoup
  
        if(aux2 == 0.0):
            ze = (h[z]+h[z+1])/2 
        else:
            ze   = (h[z]+h[z+1])/2 + (minval - dmin)*aux1/aux2
 
    except:
        ze = (h[z]+h[z+1])/2
    # 1.2) metalimnion-hypolimnion interface:
    
    z = zt                                  
    dmin =  (rho[z+1] - rho[z])/(h[z] - h[z+1])                                     
    
    try:
        while (dmin > minval) and (z < qt-2) and (z > 0):    
            z=z+1                                                                     
            dmin =  (rho[z+1] - rho[z])/(h[z] - h[z+1])                                   

      
        aux1 = (h[z-1]-h[z+1])/2 
        rhodown = (rho[z+1]-rho[z])/(h[z]-h[z+1])
        rhoup   = (rho[z]-rho[z-1])/(h[z-1]-h[z])
    
        aux2 = rhodown - rhoup
    
        if (aux2 == 0.0):
            zh = (h[z]+h[z-1])/2                                                
        else:
            zh   = (h[z]+h[z-1])/2 + (minval - dmin)*aux1/aux2                 
    except: 
        zh = (h[z]+h[z-1])/2 
        
    ze, zh, error = consistency(ze,zh,h,z0)
  
    return ze,zh,error
    

def thermal_stability(qt,h,H,tau):
#
#   External Function: Thermal Stratification Analysis
#   Function to compute the stratification conditions along depth
#
#   rho     output      Water density along depth (1d N-array)
#   n2d     output      Brunt-Vaisalla frequency (1d N-1-array)       
#   hmid    output      Grid water level for middle layer (from ref. level)
#   glin    output      Reduced gravity (1d N-1-array)
#
    rho    = np.zeros(qt,float)
    n2d    = np.zeros(qt-1,float)
    hmid   = np.zeros(qt-1,float)
    glin   = np.zeros(qt-1,float)
    
    for z in range(qt):
        rho[z] = commission(tau[z])
    
    for z in range(qt-1):
        hmid[z]   =  abs(np.mean([h[z],h[z+1]])) 
    
        glin[z]     = 9.81*abs(rho[z+1] - rho[z])/rho[z+1]

        if(h[z] != h[z+1]):
            n2d[z]    = math.sqrt(glin[z]/(h[z]-h[z+1]))
        else:
            n2d[z]    = math.sqrt(glin[z]/0.01)     # same depth
    
    return rho, n2d, hmid, glin
    
def density_2layer (qt,h,tau,H,z0):
#
#   External Function: Thermal Stratification Analysis
#   Function to determine two-layer structure
#
#   qt      input       number of sensors (integer number)
#   h       input       depth (1d vector)
#   tau     input       temperature (1d vector)
#   H       input       maximum depth of the system (float number)
#   z0      input       reference level
#   he/hh   output      epilimnion/hypolimnion layer thickness
#   pe/ph   output      epilimnion/hypolimnion water density
#   pu/pd   output      superficial/bottom water density
    
    thermo, rho, error = thermocline_depth (qt,h,tau)
    
    hh = thermo - z0   
    he = H - thermo
    
    pe,ph,npe,nph=0,0,0,0

    for z in range(qt) :
        if(thermo < h[z]):
            pe   =  pe + rho[z]
            npe  =  npe +  1
        else:
            ph   =  ph + rho[z]
            nph  = nph +  1

    if npe == 0:
        npe = 1
        pe  = rho[0]
       
    pe = pe/npe       
    ph = ph/nph     
    pu = rho[0]       
    pd = rho[qt-1]    
    
    
    return he, hh, pe, ph, pu, pd, error

def density_3layer (qt,h,tau,minval,H,z0):
#
#   External Function: Thermal Stratification Analysis
#   Function to determine two-layer structure
#
#   qt      input       number of sensors (integer number)
#   h       input       depth (1d vector)
#   tau     input       temperature (1d vector)
#   minval  input       metalimnion threshold 
#   H       input       maximum depth of the system (float number)
#   z0      input       reference level
#   h#      output      layer thickness (layer #)
#   p#      output      water density (layer #)
    
    ze, zh, error = metathick (qt,h,tau,minval,z0)  
    
    try:
        
        h1 = H  - ze
        h2 = ze - zh
        h3 = zh - z0
        
       
        p1, p2, p3, np1, np2, np3 = 0, 0, 0, 0, 0, 0
       
        rho = commission(tau)
            
        for z in range(qt):
            if(ze < h[z]):
               p1   =  p1 + rho[z]
               np1  =  np1 +  1
               
            else:                
                if(zh < h[z]):
                   p2   =  p2 + rho[z]
                   np2  =  np2 +  1
                
                else:
                   p3   =  p3 + rho[z]
                   np3  =  np3 +  1                                     
        p1 = p1/np1

                
        if np3 == 0: p3 = rho[-1]
        else: p3 = p3/np3
        
        
        if np2 == 0: p2 = (p1+p3)/2
        else: p2 = p2/np2
    
    except:
        h1 = (H-z0)/3
        h2 = (H-z0)/3
        h3 = (H-z0)/3        
        
        try:
            p1 = np.nanmean(np.nanmean(rho[0:int(qt/3)]))
            p2 = np.nanmean(np.nanmean(rho[int(qt/3):int(2*qt/3)]))
            p3 = np.nanmean(np.nanmean(rho[int(2*qt/3):qt]))
        except:
            p1 = np.nanmean(rho[0])
            p2 = np.nanmean(rho[int(qt/2)])
            p3 = np.nanmean(rho[-1])
            
    return h1,h2,h3,p1,p2,p3,error

def thickness_decomposition(psi,zph, rhoph):
#
#   External Function: Thermal Stratification Analysis
#   Function to determine the thickness of each layer (model estimator)
#    
    lmode = len(psi[0,:])
    lzph  = len(zph)-1
    
    layer_depth, layer_rho = [], []
    for m in range(lmode):       # loop in modes
        
        laymode, layrho = [], []
        
        rho,i = rhoph[0],1
        for z in range(lzph):    # loop in depth
            
            if psi[z,m]/psi[z+1,m] < 0:
                
                laymode.append(interpolation(zph[z],zph[z+1],psi[z,m],psi[z+1,m],0,0))
                layrho.append(rho/i)
                rho, i = 0, 0
                z = z + 1
                
            i = i + 1
            rho = rhoph[z] + rho
            
        layer_rho.append(layrho)
        layer_depth.append(laymode)
        
    return layer_depth, layer_rho    
        
def structure2layer(qt,h,tau,H,z0):
#
#   External Function: Thermal Stratification Analysis   
#   Function to compute the thermal structure of a two-layer system
#
    he,hh,pe,ph,pu, pd, error = density_2layer(qt,h,tau,H,z0)
    
    glin   = abs(9.81*(ph-pe)/ph)  
       
    n    = math.sqrt(glin/he)
    
    return he,hh,pe,ph, glin, n,pu,pd, error

def structure3layer(qt, h, tau, minval, H, z0):
#
#   External Function: Thermal Stratification Analysis   
#   Function to compute the thermal structure of a three-layer system
#    
    h1,h2,h3,p1,p2,p3,error = density_3layer(qt, h, tau, minval, H, z0)

    return h1,h2,h3,p1,p2,p3, error

def approx_layer(he,hh,pe,ph): 
#
#   External Function: Thermal Stratification Analysis   
#   Function to compute the three-layer the derivative method fails
#  
    h1 = he - 0.05*he
    h2 = 0.05*(he+hh) 
    h3 = hh - 0.05*hh
            
    p1 = pe
    p2 = (pe+ph)/2
    p3 = ph
    
    return h1,h2,h3,p1,p2,p3

def wedderburn(glin,he,wast,ls):
#
#   External Function: Thermal Stability Analysis
#   Function to compute the Wedderburn Number
#
#   glin   input    Reduced gravity (m/s²) 
#   wast   input    Friction velocity of the wind  (m/s)
#   he     input    Thickness of the Epilimnion end (m)
#   ls     input    System fetch length (m)  
#   wedd   output   Wedderburn Number
#    
    if wast == 0 or glin == 0:  # Unstratified/unperturbed ambient 
        wedd = None
    else:
        wedd = glin*he**2/(ls*wast**2)
        
    return wedd   


def consistency(ze,zh,h,z0):
#   
#   External Function: Statistical Analysis        
#   Function to consiste the metalimnion thickness
#        
    error = 0
    
    if np.nanmax(h) < abs(ze) or abs(ze) < np.nanmin(h) or ze < zh:

        ze = (max(h)+zh)/2
        error = 1
        
        if np.nanmax(h) < abs(zh) or abs(zh) < np.nanmin(h) or ze < zh:
            ze = z0 + 2*(max(h)-z0)/3
            zh = z0 + 1*(max(h)-z0)/3

    else:    
        if np.nanmax(h) < abs(zh) or abs(zh) < np.nanmin(h):
            zh = (ze+z0)/2
            error = 1
            
    return ze, zh, error

def chi2inv(p, nfft, nperseg, test=None):
#   
#   External Function: Statistical Analysis        
#   Function to estimate the inverse of cumulative distribution function (percentile)
#     
    if test == None:
        nw2=2*(2.5164*(nfft/nperseg))*1.2
        return chi2.ppf(p, df=nw2)/nw2
    else:
        nw2=(nfft/nperseg)
        return 2*sc.gammaincinv(nw2,p)/nw2  # Inverse incomplete gamma function


def mean_confidence_interval(data, confidence=0.99):
#   
#   External Function: Statistical Analysis        
#   Function to mean and confidence interval
# 

    dof = len(data) -1 
    mean, sigma = np.mean(data,axis=0), sem(data, axis=0)

    h = t.interval(confidence, dof, mean, sigma)

    return mean, h[0], h[1]

def conflevel(Ax,npr,dt,rho,wr,nfft,nperseg):
#   
#   External Function: Statistical Analysis        
#   Function to estimate the confidence levels based on the chi2 test (red noise)
# 
    facchi95=chi2inv(0.95, nfft, nperseg)   # Bernhardt and Kirillin (2013)
    
    fnyq=1/(2*dt)   # Nyquist frequency (half the sampling rate)


    theored=np.zeros((npr))
    for i in range(npr):
        theored[i]=(1-rho**2)/(1-(2*rho*np.cos(np.pi*wr[i]/fnyq))+rho**2)

    theoredun=theored[0];
    theored[0]=0;

    Art = np.mean(theored)
    theored[0]=theoredun
    theored=theored*(Ax/Art)    # Normalisation of the spectrum

    tabtchi=[]
    tabtchi[:]=theored*facchi95  # Chi-square confidence levels
     
    return tabtchi

def rhoAR1(datax):
#   
#   External Function: Statistical Analysis        
#   Function to calculate the lag-1 autocorrelation coefficient 
#   for an AR1 autocorrelation of data.
#
    nrho=len(datax)
    rho=0
    sommesup=0
    sommeinf=0

    moy=np.sum(datax)/nrho
    datam=datax-moy

    for i in range(1,nrho):
        j=i-1
        sommesup=sommesup+(datam[i]*datam[j])
        sommeinf=sommeinf+((datam[j])**2)
        
    rho=sommesup/sommeinf
    return rho

def Rednoise(nt,rho,nsim):
#   
#   External Function: Statistical Analysis        
#   Function to calculates AR1 autocorrelation (Monte-Carlo simulation)
#
    rzero=0
    
    redtab=np.zeros((nt,nsim))
    red=np.zeros(nt)
    i,j=1,1

    srho=np.sqrt(1-rho**2)

    for i in range(nsim):
        
        white=srho*np.random.randn(1)
        red[0]=rho*rzero+white
        
        for j in range(1,nt):
            white=srho*np.random.randn(1)
            red[j]=rho*red[j-1]+white

        redtab[:,i]=red

    return redtab

def RedConf(datax,dt,nsim,nperseg):
#   
#   External Function: Statistical Analysis        
#
    # calculation of the lag-1 autocorrelation coefficient
    rho=rhoAR1(datax)
    
    nt = len(datax)
    # calculation of nsim red noise models
    redtab = Rednoise(nt,rho,nsim)
    
    datan=[]
    i=1
    
    datan=datax-np.mean(datax)

    # spectral analysis of the data
    #nfft = max(256,2**np.ceil(np.log2(abs(len(datan)))))
    nfft = nperseg
    w,po   = signal.welch(datan[:], fs=1/dt, nperseg=nperseg)
    
    # calculation of the area of the data power spectrum
    Ax=np.mean(po)
    

    for i in range(nsim):  # spectral analysis of the nsim red noise signals
        
        red2n=redtab[:,i] -np.mean(redtab[:,i])
        wr,pr = signal.welch(red2n, fs=1/dt, nperseg=nperseg)
        
    npr=len(pr)
    tabtchi = conflevel(Ax,npr,dt,rho,wr,nfft,nperseg)
    
    return wr,tabtchi
    
    
def ciout(x):
#   
#   External Function: Statistical Analysis        
#   confidence interval of x-data (assuming 95% confidence interval )   
#
    nivel = ci(x)[1]-np.mean(ci(x))
    return nivel


def average(arr,n):
#
#   External Functions: Statistical Analysis
#   Function to n-average data (arr)
#    
    n=int(n)
    
    if n == 0:
        return np.nanmean(arr)
    else:
        end = n*int(len(arr)/n)       
        return np.nanmean(arr[:end].reshape(-1,n),1)


def ci(x):
# 
#   External Function: Statistical Analysis
#   Prediction bands (with coverage probability of 95% confidence interval)
#   
    size    = x.size - np.count_nonzero(np.isnan(x)) 
    sdev    = np.nanstd(x) 
    
    val     = sdev/np.sqrt(size)
    z       = 1.96                                  # 95% confidence level
    
    dev = val*z
    
    try:   # If all x is not a number, NaN is returned as prediction bands
        upper = np.nanmean(x) + dev 
        lower = np.nanmean(x) - dev 
    except RuntimeWarning:
        upper = None
        lower = None
    
    return lower,upper



def velocityten(wz,z):
#   
#   External Function: Meteorological Analysis 
#   Convert wind velocity given at z meters to 10 meters
#    
#   wz   (input):    wind velocity at z meters (m/s) 
#   w10  (output):   wind velocity at 10 meters 
    
    k   = 0.4 #Von Karman's constant
    exp = math.log(10/z)
    
    l = len(wz)
    w10  = np.zeros(l,float)
    
    if z == 10:
        return wz
       
    for time in range(l):
        
        if wz[time] < 5:
            Cd = 0.0010
        else:
            Cd = 0.0015
        
        wind = wz[time]/(1-np.sqrt(Cd)/k*exp) 

        w10[time] = wind
    
    return w10


def wind_average(wd,iw):
#
#   External Function: Meteorological Analysis
#   Wwind direction average 
#  

    u_east  = np.average(iw*np.sin(np.deg2rad(wd)))
    u_north = np.average(iw*np.cos(np.deg2rad(wd)))

    mean   = (np.arctan2(u_east, u_north)) * 180/np.pi
    mean  = (360+mean)%360
    
    return mean

   

def windbar_fetch(dw_mean,linang,angle,dists):
#
#   External Function: Meteorological Analysis
#   Wind fetch contribution
#    
#   Specify the wind fetch for each specified angle when the mean wind
#   event is between the treshold angles of wind contribution    
#
#   dw_mean     input       mean wind direction
#   linang      input       angle treshold for wind contribution
#   angle       input       array of    angles specified by user in .fet
#   dists       input       array of distances specified by user in .fet
#   Ldist       output      array of distances when the angle match the specified criteria
#    
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



def wind_angle(wind,linang):
#
#   External Function: Meteorological Analysis
#   Boundaries of wind direction contribution
#    
#   Specify the maximum and minimum direction of wind contribution  
#
#   wind           input        wind direction
#   linang         input        angle treshold for wind contribution
#   dw_min         output       minimum wind direction contribution
#   dw_max         output       maximum wind direction contribution     
    
    dw_max = wind + linang
    dw_min = wind - linang

    if (dw_max > 360):
        dw_max = dw_max - 360

    if (dw_min < 0):
        dw_min = 360 + dw_min 
    
    return dw_min, dw_max



def wind_stress(w):
#
#   External Function: Meteorological Analysis
#   Function to compute the wind stress (N/m² or Pa)
# 
    if(w > 5):
        Cd = 0.0015      # drag coefficient
    else:
        Cd = 0.0010
        
        if (w == 0.0):   # low wind to characterize the system stability
            w = 0.01
              
    rho_air = 1.225      # density of the air (kg/m³)
    stress = Cd*rho_air*w**2
    
    return stress

def wind_parameters(w,rw,pe,he,n,glin,H):
#
#   External Function: Meteorological Analysis
#   Function to general parameters related to system stability (1D)
#     
#   stress      output      Wind stress (N/m²) 
#   wast        output      Friction velocity of the wind  (m/s)
#   riw         output      Richardson number 
    
    stress  = wind_stress(w)     
    wast    = math.sqrt(stress/pe) 
    
    try:
        riw = glin*he/((wast**2))
        
    except RuntimeWarning:  # Nan is used for unstable and non-stratified regions  
        riw=  None
    
    return stress, wast, riw

def richardson(w,rw,qt,h,pe,H,n2d,hmid,p,glin):
#
#   External Function: Meteorological Analysis
#   Function to general parameters related to system stability along depth (2D)
#     
#   riw2d       output      Richardson number (along z-direction)
#
    riw2d   = np.zeros(qt-1,float)
    
    for z in range(qt-1):

        win_stress = wind_stress(w)
        
        wast        = math.sqrt(win_stress/abs(np.mean([p[z+1],p[z]])))     
        riw2d[z]    =  (glin[z]*( H - hmid[z]))/((wast**2))                  
    
    return riw2d    

def deriva(y1,y2,x1,x2):
#   
#   External Function: Mathmatical Analysis
#   Function to calculate a discrete derivative (gradient)
#   
    try:
        der = abs(y1-y2)/abs(x1-x2)
    except RuntimeWarning:
        der = abs(y1-y2)/0.00001
    
    return der

def interpolation(y1,y2,x1,x2,x,auxiso):
#
#   External Function: Mathmatical Analysis
#   Result from linear Interpolation betweem two points (x,y) for x 
#    

    if(y1 == y2 or x1 == x2):

        return np.mean([y1,y2])
    
    else: 
        f = interpolate.interp1d([x1,x2], [y1,y2])    
        
        try:                   # try to fit the function
            return f(x)
        
        except ValueError:     # if x is out of bound, auxiso is returned
                               # auxiso is the last value of f(x) 
            return auxiso    


#def confidence (ff,k):
#
#   External Function: Mathmatical Analysis
#   Connfidence interval for PSD considering conf. level of 95%
#      
#    conf = 0.95   
    
#    chi_val_95 = chi2.isf(q=(1-conf)/2, df=k)
#    y=(ff/len(ff))*(chi_val_95/k)

#    return y

def sorting_2d(data):
#
#   External Function: Conditional Analysis
#   Function to get the a 2d-data sorted 
#    
    ordered = np.zeros((len(data),len(data[0])),float)

    for time in range(len(data)):
        
        data_aux = data[time,:]
        
        temporary = []
        first     = 'on'
        
        for i in range(len(data_aux)):
            
            if np.isnan(data_aux[i]) == False:
                temporary.append(data_aux[i])
            
                if first == 'on':
                    first = 'off'
                    idx   = i
                
        temporary = -np.sort(-np.array(temporary))
    
        auxiliary     = np.full(len(data_aux),np.nan)
        auxiliary[idx:len(temporary)+idx] = temporary
        
        ordered[time][:] = auxiliary
    
    return ordered

def sorting_1d(data):
#
#   External Function: Conditional Analysis
#   Function to get the a 1d-data sorted 
#   
    temporary = []
    first     = 'on'
        
    for i in range(len(data)):
            
        if np.isnan(data[i]) == False:
            temporary.append(data[i])
            
            if first == 'on':
                first = 'off'
                idx   = i
                
    temporary = -np.sort(-np.array(temporary))
    
    auxiliary     = np.full(len(data),np.nan)
    auxiliary[idx:len(temporary)+idx] = temporary
    
    return auxiliary


def find_nearest(array,value):
#
#   External Function: Conditional Analysis 
#   Find the nearest 'value' from a 'array' and returns the value and its index
#
    idx = (np.abs(array-value)).argmin()
    
    return array[idx],idx     # returns the element and its index

def mask (serie):
#
#   External Function: Conditional Analysis
#   Function to filter not a number (NaN) values 
#    
    xi = np.arange(len(serie))
    
    m = np.isfinite(serie)
    xfiltered = np.interp(xi, xi[m], serie[m])
    
    return xfiltered
   
def class_generation (riw,hh,he,ls):
#
#   External Function: Classification Analysis
#   Function to classify the lake mixing/internal wave activity
#
#   Returns     output code
#
#    Regime 1   Stratification is broken down by mixing
#    Regime 2   Large internal displacement accompanied by billows and mixing
#    Regime 3   Internal seiche is dominant
#    Regime 4   Internal seiche with short amplitude.
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



def spigel(he,H,W1,W2):
#
#   External Function: Classification Analysis
#   Function to provide the BSIW amplitude based on Spigel and Imberger theory
#     
    zeta = 0.4996*he/W1
    if zeta > H:    zeta = 0.4996*he/W2
    
    return 0.4996*he/W1

def bueno(he,hh,W1,W2, typ=None):
#
#   External Function: Classification Analysis
#   Function to provide the BSIW amplitude based on Bueno and Bleninger theory
#    
    
    xi = 0.1
    k1 = 6
    k2 = 0.25
    H  = he/(he+hh) 
    
    g   = 12.156*(H**3)-15.714*(H**2)+2.8426*H+2.0846
    f   = g*np.exp((H**2)/k2)
    
    if typ == 'amplitude':
        
        if W1 < 5:
            amp = 0.4996*he/W1
            if amp > he+hh:
                amp = 0.4996*he/W2
        
        return xi*he*np.exp((W1-k1)**2/(2*f**2))
    
    amp = []
    for i in range(len(W1)):
        amp.append(xi*he*np.exp((W1[i]-k1)**2/(2*f**2)))
    
    return np.array(amp)
    

def iw_generation (wedd, hh, he, ls):
#
#   External Function: Classification Analysis
#   Function to identify Regime 3 (Spigel and Imberger, 1980)
#
#   lower_gene      output      lower W number for regime 3
#   upper_gene      putput      highest W number for regime 3
#
    aux1 = math.sqrt((hh+he)/hh)
    
    lower_gene = 0.5*aux1
    upper_gene = ls/(4*he)*aux1**2
    
    return lower_gene, upper_gene
    
def amplitude_average(iw,he,lim):
#
#   External Function: Classification Analysis
#   Function used to estimate the maximum amplitude of BSIW based on W
#   Estimation based on Spigel and Imberger (1980)
#    
#   iw         input       Wedderburn number (-)
#   he         input       epilimnion thickness
#   lim        input       Wedderburn number limit (-)
#   average    output      maximum amplitude according to theory related to the averaged W number     
#   ratio      output      percentage of period with W lower than lim  
#
    count, summ = 0, 0
    lim   = 1/lim
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


def welch_method(serie,size,w, dt):
#
#   External Function: Spectral Analysis
#   Function to compute the PSD based on Welch's method
#
#   serie      input       time-series data
#   size       input       size of the window (in seconds)
#   w          input       window type (char)
#   dt         input       time step (hours)
#   freq       output      frequency (Hz)  
#   per        output      period associated to frequency (hours)  
#   ff         output      power spectral density (x²/Hz)  
#    
    if (size == 0):              # when windowsize is not defined 
        size =  4*24*60*60       # 4 days window (generally used to meteorological data)

    nsim = 10                    # number of Monte Carlo simulations;
    dt   = 60*60*dt              # dt in hour 
    n    = int(size/dt)          # nperseg
    
    if n > len(serie):  # Warning: the specified window size is larger than the time-serie
        n = len(serie)
    
    serie    = mask(serie) 
    serie = serie - np.mean(serie)
    
    freq, ff = signal.welch(serie, fs=1.0/dt, window=w, nperseg=n, detrend='linear', axis=-1)
    
    wr, conf = RedConf(serie,dt,nsim,n)
    
    per  = 1/freq[1:]/60/60           # convert frequency (Hz) to period (h)
        
    return freq[1:], per, ff[1:], wr, conf


def wave_spectral (si, dt, mother): 
#
#   External Function: Spectral Analysis
#   Function to compute the wavelet transfor of signal 'si'
#
#   si         input       time-series data
#   mother     input       mother function type (char)
#   dt         input       time step (hours)
#   time       output      period (hours)  
#   per        output      period associated to frequency (hours)  
#   power      output      power of wavelet (x²)  
#  
#   Wavelet parameters:
# 
    pad = 1         # pad the time series with zeroes (recommended)
    dj = 0.10       # this will do 10 sub-octaves per octave
    s0 = 2.*dt      # this provide a analysis at a scale of 8*dt time-scale
    j1 = 15./dj     # this provides 15 powers-of-two with dj sub-octaves each
#
#      

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


def coherence_shift(s1,s2,window,dt):
#
#   External Function: Spectral Analysis
#   Function to compute coherence between two signals
#
    dt     = 60*60*dt 
    nperseg = int(window/dt)
    
    fcoh, cxy = signal.coherence(s1, s2,  1/dt, nperseg=nperseg)
    px, pxy = signal.csd(s1, s2, 1/dt, nperseg=nperseg)
    
    phase = np.angle(pxy)
    phase[phase<-90] += 360
    
    #calculate 95 % confidence level
    edof = (len(s1)/(nperseg/2)) * cxy.mean() 
    gamma95 = 1.-(0.8)**(1./(edof-1.))

    conf95 = np.where(cxy>gamma95)
    
    return phase, cxy , fcoh, conf95

def butter_bandpass(lowcut, highcut, fs, order):
#
#   External Function: Spectral Analysis
#   Function to compute the coefficeints of the band pass filter
#  
    nyq = 0.5*fs
    
    low  = lowcut/nyq
    high = highcut / nyq
    
    sos = signal.butter(order,[low, high], analog=False, btype='band', output='sos')

    return sos


def butter_bandpass_filter(data, lowcut, highcut, fs):
#
#   External Function: Spectral Analysis
#   Function to band pass filter a signal data
#
    data   = mask(data) 
    data   = data - np.mean(data)
    sos   = butter_bandpass( lowcut, highcut, fs, 1)
    f = signal.sosfilt(sos, data)
    
    return f  