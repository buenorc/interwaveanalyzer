# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer modelling function

modififications to previosly versions must be specified here using WAVE codes:
    
W-23.05-1.00.0-00
A-01.05-1.00.3-00
V-22.05-1.00.3-00
E-05.05-1.00.3-00
"""

import numpy as np
import math as ma

from internal_module import thickness_decomposition



def decomposition(rho,z): 
#
#   Model Function: Mode decomposition (numerical modeling)    
#     
    nmodes = int(len(z)/2)
    
    if  int(len(z)/2) > 10:
        nmodes = 10
    
    dz   = np.diff(z)

    hs   = z[:-1]+np.median(np.diff(z))/2    
    
    n2    = np.zeros((len(hs)),float)
    rhoph = np.zeros((len(hs)),float)
    
    for i in range(len(hs)-1):
        
        rhoph[i] = (rho[i]+rho[i+1])/2
        n2[i]  = 9.81*(rho[i]-rho[i+1])/rho[i]/hs[i]

    N2 = n2+(np.arange(len(hs)))*10**-9
    N2 = np.flip(N2)    
    
    # First we are solving for w on dz,2dz,3dz...H-dz
    M = np.shape(N2)[0]-1
    N2mid = N2[:-1]+np.diff(N2)/2.
    
    # matrix for second difference operator
    D = np.diag(-2.*np.ones(M),0)
    D  += np.diag(1.*np.ones(M-1),-1)
    D += np.diag(1.*np.ones(M-1),1)


    for i in range(M):
        for j in range(M):
            D[i,j] = -D[i,j]/dz[i]/dz[i]
    D = np.diag(1./N2mid).dot(D)
    
 
    ce,W = np.linalg.eig(D)
    
    for i in range(M):
        for j in range(M):
            W[i,j] = W[i,j]/np.sqrt(dz[i])
    ce = 1./np.sqrt(ce)
    ind=np.argsort(-ce)
    
    ce=ce[ind[:-2]]
    W=W[:,ind[:-2]]
    
    psi = np.zeros((M+1,M+1-3))
    psi[0,:] = W[0,:]
    psi[1:-1,] = np.diff(W,axis=0)
    psi[-1,:] = -W[-1,:]
    

    A = np.zeros((M-2))

    for i in range(M-2):
        A[i] = np.sqrt(np.sum(psi[:,i]*psi[:,i],axis=0)*dz[i])
    
    psi = psi/A
    
    psi[:,psi[0,:]<0] *= -1
    
    zphi = np.zeros(M+1)
    
    for i in range(M+1):
        zphi[i] = (z[i] + z[i+1])/2
        
    psi_req = np.flip(psi[:,:nmodes],0)
    layer_depth, layer_rho = thickness_decomposition(psi_req,zphi, rhoph)

    return psi_req,zphi,rhoph, nmodes


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
