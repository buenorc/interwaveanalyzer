# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer modelling function

modififications to previosly versions must be specified here using WAVE codes:
    
W-23.05-1.00.0-00
A-01.05-1.00.3-00
V-22.05-1.00.4-00
E-05.05-1.00.4-00
"""

import numpy as np
import math as ma

from internal_module import thermal_stability



def decomposition(L,tau,h, H):

    N    = 100         # resolution of the refined grid

    increment = 3600  # 3600 seconds = 60 minutes = 1 h 
    k    = np.pi/(L)     # wave number


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

    if h[0] < 50:                 # if the 1st sensor (near the water surface) is not deeper than 50 m  
    
        for i in range(N):         # contour that represents a depth vector in meters (1 m to 100 m)
        
            ii = 1                     # second index (ii) to evaluate

            while h[ii] < refined_depth[i] :   # 
                ii = ii + 1
        
            if buoy[ii-1] > -1:
                bv[i] = buoy[ii-1]
            else:
                if h[ii-1] < 50:
                    bv[i] = 0.3
                else:
                    bv[i] = 0.08
        
            nsq[i] = bv[i]**2/3600  # buoyancy frequency (N²) - 1/min² to 1/second²
        
    else:                              # there is not sensor near the surface (in the first 50 m)

        for i in range(N):         # contour that represents a depth vector in meters (1 m to 100 m)
            ii = 1                     # second index (ii) to evaluate
            while h[ii] > refined_depth[i] :   # 
                ii = ii + 1
        
            if buoy[ii-1] > -1:
                bv[i] = buoy[ii-1]
            else:
                if h[ii-1] < 50:
                    bv[i] = 0.3
                else:
                    bv[i] = 0.08
        
            nsq[i] = bv[i]**2/3600  # buoyancy frequency (N²) - 1/min² to 1/second²
        
    
# calculate the first vertical modes w[m,:]
# find the approximated value of p depending on the mode that is defined

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
                W[i] =  -W[i-2] + f * W[i-1]
            
            
            finold = finnew
            finnew = W[-1]
        
            e = finold*finnew
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
                W[i] =  -W[i-2] + f * W[i-1]  
        

            finnew = W[-1] 
        
            if finnew < 0:
                randunten = p
            else:
                randoben  = p
    
        normw=np.sqrt(sum((W*nsq)*np.transpose(W)))
    
        for i in range(N):
            w[m,i] = W[i]/normw
            n[m,i] = np.sqrt(nsq[i])*W[i]/normw
        
        for i in range(N-1):
            hor[m,i] = w[m,i+1] - w[m,i]
    
        finnew = finold
        
        peri.append(p/3600)  # hour
        p = p + increment

        conv.append((nsq)*w[m,:])
     
    # hortief = np.zeros((a-1),float)
    #
    # for i in range(a-1):
    #     hortief[i] = -1*h[i]
    
    vel = np.transpose(hor)
    

    return vel, conv, refined_depth, peri


def disp_zmodel (pe,ph,he,hh,L,m):
#
#   Model Function: 1D non-hydrostatic analytical model for two-layer system
#   Returns the period of BSIW modes considering basin length variation
#
    gamma = pe/ph
    
    peri_min = biquadratic(L[0],he,hh,gamma,m)
    peri_ave = biquadratic(L[1],he,hh,gamma,m)
    peri_max = biquadratic(L[2],he,hh,gamma,m)
       
    return peri_min, peri_ave, peri_max

def disp_xmodel3(p1,p2,p3,h1,h2,h3,L,vertical,m):
#
#   Model Function: 1D hydrostatic analytical model for three-layer system
#   Returns the period of BSIW modes considering basin length variation
#  
    gamma12 = p1/p2
    gamma13 = p1/p3
    gamma23 = p2/p3
       
    A = [[h1, h1,  h1], [h2*gamma12,  h2, h2], [h3*gamma13,  h3*gamma23, h3]]
    
    solv = np.linalg.eigvals(A)
   
    pv1_min, pv2_min  = eigen3_values(L[0],solv[0],solv[1],m)
    pv1_ave, pv2_ave  = eigen3_values(L[1],solv[0],solv[1],m)
    pv1_max, pv2_max  = eigen3_values(L[2],solv[0],solv[1],m)
    
    if(vertical==1):
        return pv1_min,pv1_ave,pv1_max
    else:   
        return pv2_min,pv2_ave,pv2_max

def disp_xmodel4(p1,p2,p3,p4,h1,h2,h3,h4,L,vertical,m):
#
#   Model Function: 1D hydrostatic analytical model for four-layer system
#   Returns the period of BSIW modes considering basin length variation
#      
    gamma12 = p1/p2
    gamma13 = p1/p3
    gamma23 = p2/p3
    gamma14 = p1/p4
    gamma24 = p2/p4
    gamma34 = p3/p4
    
    A = [[h1, h1,  h1, h1], \
               [h2*gamma12,  h2, h2, h2],\
               [h3*gamma13,  h3*gamma23, h3, h3], \
               [h4*gamma14,  h4*gamma24, h4*gamma34, h4]]
    
    solv = np.linalg.eigvals(A)
    #sorted(solv)
   
    pv1_min, pv2_min, pv3_min = eigen4_values(L[0],solv[0],solv[1],solv[2],m)
    pv1_ave, pv2_ave, pv3_ave = eigen4_values(L[1],solv[0],solv[1],solv[2],m)
    pv1_max, pv2_max, pv3_max = eigen4_values(L[2],solv[0],solv[1],solv[2],m)
    
    if(vertical==1):
        return pv1_min,pv1_ave,pv1_max
    if(vertical==2): 
        return pv2_min,pv2_ave,pv2_max
    else:
        return pv3_min,pv3_ave,pv3_max

def coriolis_effect(fo,to):
#
#   Model Function: Correction model to account Coriolis effect 
#   Returns the period of BSIW modes
#      
    aux  = 4*np.power(np.pi,2) + np.power(to,2)*np.power(fo,2)  
    peri = np.sqrt(4*np.power(np.pi,2)*np.power(to,2)/aux) 
    
    return peri   

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

def sensitivity_3layer(mean,diff,N,p1,p2,p3,h1,h2,h3,fetch,typ):
#
#   Sensitivity Function: Sensitivity analysis to check period variation 
#   based on layer thickness and water density changes (three-layer system)
#      
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


def sensitivity_dimension(L, pe,ph, he, hh):
#
#   Sensitivity Function: Sensitivity analysis to check period variation 
#   Combined with variation on layer thickness and water density, as well as,
#   on lake length 
#    
    N    = 50
    g    = 9.81
    drho = 4.00
    ddep = 0.10
    
    
    mean_rho = np.sqrt(ph/(ph-pe)) 
    mean_dep = np.sqrt((he+hh)/(he*hh))
    
    aux_rho  = mean_rho - drho
    aux_dep  = mean_dep - ddep
    
    xrho   = np.linspace(aux_rho, mean_rho+drho, N)
    xdep   = np.linspace(aux_dep, mean_dep+ddep, N)
    
    yrho   = np.zeros((N),float)
    ydep   = np.zeros((N),float)
    
    Crho = 2*L/np.sqrt(g*he*hh/(he+hh))  
    Cdep = 2*L/np.sqrt(g*(ph-pe)/ph)

    for i in range(N):
        
        yrho[i] = Crho*xrho[i]
        ydep[i] = Cdep*xdep[i]
        
               
    return xrho, xdep, yrho/60/60, ydep/60/60  # period in hours
