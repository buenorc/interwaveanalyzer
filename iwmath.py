# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer - Modelling functions

Interwave Analyzer - Version 2 (2026) 
Modelling functions module version: 2.260305

-------------------------------------------------------------------------------

de Carvalho Bueno, R; Bleninger, T. B.; Lorke, A. 
Internal wave analyzer for thermally stratified lakess 
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
import numpy as np
import math as ma

from iwmod import thermal_stability


# Decomposition model
def decomposition(L,tau,h, H):

    cond = 0          # when cond is 1 a warnning is raised!
    N    = 100        # resolution of the refined grid

    increment = 3600  # 3600 seconds = 60 minutes = 1 h 
    k    = np.pi/(L)  # wave number


    # depth classification
    a = len(tau)

    tau  = np.concatenate((tau,[-9999]))
    h = np.concatenate((h,[H]))
    H = h[-1]

    refined_depth = np.linspace(0,H,N)
    dz = refined_depth[1] - refined_depth[0]

   
    bv  = np.zeros((N),float)
    nsq = np.zeros((N),float)

    rho, buoy, hmid, glin = thermal_stability(a,h,H,tau)

    buoy = np.concatenate((buoy,[np.nan]))
    buoy = buoy*60  # = [1/min]

    if h[0] < 50:            # if the 1st sensor is not deeper than 50 m  
    
        for i in range(N):   # contour that represents a depth in meters (1 m to 100 m)
            ii = 1           # second index (ii) to evaluate

            while h[ii] < refined_depth[i] :   # 
                ii = ii + 1
        
            if buoy[ii-1] > -1:
                bv[i] = buoy[ii-1]
            else:
                if h[ii-1] < 50:
                    bv[i] = 0.3
                else:
                    bv[i] = 0.08
        
            nsq[i] = bv[i]**2/3600  # buoyancy frequency (N²) to 1/second²
        
    else:       # there is not sensor near the surface (in the first 50 m)

        for i in range(N):    # contour that represents a depth in meters (1 m to 100 m)
            ii = 1            # second index (ii) to evaluate
            
            try:
                while h[ii] > refined_depth[i]:   # 
                    ii = ii + 1
            except:
                    ii = ii - 1
                    
            if buoy[ii-1] > -1:
                bv[i] = buoy[ii-1]
            else:
                if h[ii-1] < 50:
                    bv[i] = 0.3
                else:
                    bv[i] = 0.08
        
            nsq[i] = bv[i]**2/3600  # buoyancy frequency (N²) to 1/second²
        
    
# calculate the first vertical modes w[m,:]
    W = np.ones((N),float)

    W[0] = 0
    W[1] = 1

    finnew = 0
    p      = 0 # internal seiche period


    pnew      = p
    peri      = []
    conv      = []

    w   = np.zeros((5,N),float)
    n   = np.zeros((5,N),float)
    hor = np.zeros((5,N-1),float) 

    for m in range(5):
        e=0.5
        while e >= 0:
            for i in range(2,N):
                f = 2 - k**2*dz**2*(nsq[i-1]*(p**2)/(2*np.pi)**2 - 1 )
                try:
                    W[i] =  -W[i-2] + f * W[i-1]
                except:
                    cond = 1
                    W[i] = np.nan
            
            finold = finnew
            finnew = W[-1]
        
            try:
                e = finold*finnew
            except:
                e = float('inf')
            
            
            pold = pnew
            pnew = p
            p = p + increment

        
        if finnew > 0:
            randoben  = pnew
            randunten = pold
        else:
            randoben  = pold
            randunten = pnew  
    
        finold = finnew
    
    # halfing of the intervalls
        while abs(randunten-randoben)>0.1:

            p=0.5*(randunten+randoben)
        
            for i in range(2,N):
                f = 2 - k**2*dz**2*(nsq[i-1]*(p**2)/(2*np.pi)**2 - 1 )
                try:
                    W[i] =  -W[i-2] + f * W[i-1]  
                except:
                    cond = 1
                    W[i] = np.nan

            finnew = W[-1] 
        
            if finnew < 0:
                randunten = p
            else:
                randoben  = p
        
        try:
            normw=np.sqrt(sum((W*nsq)*np.transpose(W)))
        except:
            normw=1
    
        for i in range(N):
            w[m,i] = W[i]/normw
            n[m,i] = np.sqrt(nsq[i])*W[i]/normw
        
        for i in range(N-1):
            hor[m,i] = w[m,i+1] - w[m,i]
    
        finnew = finold
        
        peri.append(p/3600)  # hour
        p = p + increment

        conv.append((nsq)*w[m,:])
    
    vel = np.transpose(hor)
    

    return vel, conv, refined_depth, peri, cond

# 1D non-hydrostatic analytical model for two-layer system
def disp_zmodel (pe,ph,he,hh,L,m):
    
    gamma = pe/ph
    
    Lmin = np.nanmin(L)
    Lave = np.nanmean(L)
    Lmax = np.nanmax(L)
    
    peri_min = biquadratic(Lmin,he,hh,gamma,m)
    peri_ave = biquadratic(Lave,he,hh,gamma,m)
    peri_max = biquadratic(Lmax,he,hh,gamma,m)
       
    return peri_min, peri_ave, peri_max

# 1D hydrostatic analytical model for three-layer system
def disp_xmodel3(p1,p2,p3,h1,h2,h3,L,vertical,m):
 
    gamma12 = p1/p2
    gamma13 = p1/p3
    gamma23 = p2/p3
       
    A = [[h1, h1,  h1], [h2*gamma12,  h2, h2], [h3*gamma13,  h3*gamma23, h3]]
    
    solv = np.linalg.eigvals(A)
    
    Lmin = np.nanmin(L)
    Lave = np.nanmean(L)
    Lmax = np.nanmax(L)
   
    pv1_min, pv2_min  = eigen3_values(Lmin,solv[0],solv[1],m)
    pv1_ave, pv2_ave  = eigen3_values(Lave,solv[0],solv[1],m)
    pv1_max, pv2_max  = eigen3_values(Lmax,solv[0],solv[1],m)
    
    if(vertical==1):
        return pv1_min,pv1_ave,pv1_max
    else:   
        return pv2_min,pv2_ave,pv2_max



def coriolis_effect(To, lat):
    """
    Estimate the internal seiche period including Coriolis effect.
    
    Parameters
    ----------
    To : float
        Wave period without Coriolis effect [seconds].
    lat : float
        Latitude [degrees].
    
    Returns
    -------
    Tcor : float
        Estimated period including Coriolis effect [seconds].
    """
    Omega = 7.292115e-5  # Earth's rotation rate [s^-1]
    phi = np.radians(lat)
    f = 2 * Omega * np.sin(phi)
    omega0 = 2 * np.pi / To
    omega = np.sqrt(omega0**2 + f**2)
    Tcor = 2 * np.pi / omega
    return Tcor

def biquadratic(L,he,hh,gamma,m):
#
#   Solver Function: Solver for 2-layer model
#
    g = 9.81
    k = m*np.pi/L
    
    th   = np.tanh(k*hh)
    te   = np.tanh(k*he)

    p     = (gamma*te*th + 1)/(k*th)
    q     = -g*(th+te)/th
    r     = -ma.pow(g,2)*(gamma-1)*k*te

    delta = ma.pow(q,2) - 4*p*r
    
    try:
        omega = np.sqrt((-q-np.sqrt(delta))/(2*p))
    except RuntimeWarning:
        return None
    
    peri  = 2*np.pi/omega
    
    return peri 


def eigen3_values(L,lambv1,lambv2,m):
#
#   Solver Function: Solver for 3-layer model
#    
    g     = 9.81   # m/s²

    try:
        peri_v1 =  2*L/(m*np.sqrt(g*lambv1))     # V1 interfacial period (sec) 
    except RuntimeWarning:
        peri_v1 = None
        
    try:
        peri_v2 =  2*L/(m*np.sqrt(g*lambv2))     # V2 interfacial period (sec)
    except RuntimeWarning:
        peri_v2 = None
        
    return peri_v1, peri_v2


def eigen4_values(L,lambv1,lambv2,lambv3,m):
#
#   Solver Function: Solver for 4-layer model
#   
    g     = 9.81   # m/s²

    try:
        peri_v1 =  2*L/(m*np.sqrt(g*lambv1))     # V1 interfacial period (sec) 
    except RuntimeWarning:
        peri_v1 = None
    try:
        peri_v2 =  2*L/(m*np.sqrt(g*lambv2))     # V2 interfacial period (sec)
    except RuntimeWarning:
        peri_v2 = None
    try:
        peri_v3 =  2*L/(m*np.sqrt(g*lambv3))     # V3 interfacial period (sec)
    except RuntimeWarning:
        peri_v3 = None
        
    return peri_v1, peri_v2, peri_v3



def sensitivity_2layer(mean,diff,N,pe,ph,he,hh,fetch,typ):
#
#   Sensitivity Function: Sensitivity analysis to check period variation 
#   based on layer thickness and water density changes (two-layer system)
#      
    x   = np.linspace(mean-diff, mean+diff, N)
    period = np.zeros((N),float)

    for i in range(N):
        
        if   typ == 1:
            _,per,_  = np.array(disp_zmodel(x[i], ph, he, hh,fetch,1))
        elif typ == 2:
            _,per,_  = np.array(disp_zmodel(pe, x[i], he, hh,fetch,1))
        elif typ == 3:
            _,per,_  = np.array(disp_zmodel(pe, ph, x[i], hh,fetch,1))
        elif typ == 4:
            _,per,_  = np.array(disp_zmodel(pe, ph, he, x[i],fetch,1))            
            
        period[i] = per     
        
    return x, period/60/60  

#  Sensitivity analysis to check period variation
def sensitivity_3layer(mean,diff,N,p1,p2,p3,h1,h2,h3,fetch,typ):
     
    aux = mean - diff
    x   = np.linspace(aux, mean+diff, N)
    period = np.zeros((N),float)

    for i in range(N):

        if   typ == 1:
            _,per,_  = np.array(disp_xmodel3(x[i],p2,p3,h1,h2,h3,fetch,2,1))
        elif typ == 2:
            _,per,_  = np.array(disp_xmodel3(p1,x[i],p3,h1,h2,h3,fetch,2,1))
        elif typ == 3:
            _,per,_  = np.array(disp_xmodel3(p1,p2,p3,x[i],h2,h3,fetch,2,1))
        elif typ == 4:
            _,per,_  = np.array(disp_xmodel3(p1,p2,p3,h1,x[i],h3,fetch,2,1))
                
        period[i] = per     
        
    return x, period/60/60  

# Compute sensitivity curves for BSIW period varying
def sensitivity_dimension(L, pe, ph, he, hh):

    N = 50
    g = 9.81

    mean_rho = np.sqrt(ph / (ph - pe))
    mean_dep = np.sqrt((he + hh) / (he * hh))

    drho = 4.0
    ddep = 0.10

    xrho = np.linspace(mean_rho - drho, mean_rho + drho, N)
    xdep = np.linspace(mean_dep - ddep, mean_dep + ddep, N)

    Crho = 2 * L / np.sqrt(g * he * hh / (he + hh))
    Cdep = 2 * L / np.sqrt(g * (ph - pe) / ph)

    yrho = Crho * xrho
    ydep = Cdep * xdep

    return xrho, xdep, yrho / 3600, ydep / 3600  # in hours