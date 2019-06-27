# -*- coding: utf-8 -*-
"""
Module to compute internal waves in reservoirs (modeling)

@author: BUENO, R.
@date  : Nov 14 2017

Modified:
    

"""
# pacotes 

import numpy as np
import math as ma


def biquadratic(L,he,hh,gamma,m):
    g = 9.81
    k = m*np.pi/L
    
    th   = np.tanh(k*hh)
    te   = np.tanh(k*he)

    p     = (gamma*te*th + 1)/(k*th)
    q     = -g*(th+te)/th
    r     = -ma.pow(g,2)*(gamma-1)*k*te

    delta = ma.pow(q,2) - 4*p*r
    omega = np.sqrt((-q-np.sqrt(delta))/(2*p))
    
    peri  = 2*np.pi/omega
    
    return peri # seconds

def eigen2_values(L,lamb,m):
    
    g     = 9.81   # m/s²    
    peri =  2*L/(m*np.sqrt(g*lamb))    # interfacial period (sec) 
    
    return peri # seconds

def eigen3_values(L,lambv1,lambv2,m):
    
    g     = 9.81   # m/s²

    peri_v1 =  2*L/(m*np.sqrt(g*lambv1))     # V1 interfacial period (sec) 
    peri_v2 =  2*L/(m*np.sqrt(g*lambv2))     # V2 interfacial period (sec)
    
    return peri_v1, peri_v2

def eigen4_values(L,lambv1,lambv2,lambv3,m):

    g     = 9.81   # m/s²

    peri_v1 =  2*L/(m*np.sqrt(g*lambv1))     # V1 interfacial period (sec) 
    peri_v2 =  2*L/(m*np.sqrt(g*lambv2))     # V2 interfacial period (sec)
    peri_v3 =  2*L/(m*np.sqrt(g*lambv3))     # V3 interfacial period (sec)
    
    return peri_v1, peri_v2, peri_v3
# -----------------------------------------------------------------------------
#  
# -----------------------------------------------------------------------------

def disp_zmodel (pe,ph,he,hh,L,m):
#    
#  dispertion relation of internal waves in a 2 layered system and free surface
#  z-model
#  m : horizontal mode of internal waves
#   
    gamma = pe/ph
    
    peri_min = biquadratic(L[0],he,hh,gamma,m)
    peri_ave = biquadratic(L[1],he,hh,gamma,m)
    peri_max = biquadratic(L[2],he,hh,gamma,m)
    
    
    return peri_min, peri_ave, peri_max

def disp_xmodel2(p1,p2,h1,h2,L,m):
    
    gamma12 = p1/p2   
    A = [[h1, h1],[h2*gamma12,  h2]]
    
    solv = np.linalg.eigvals(A)
    #sorted(solv)

    peri_min = eigen2_values(L[0],solv[0],m)
    peri_ave = eigen2_values(L[1],solv[0],m)
    peri_max = eigen2_values(L[2],solv[0],m)  
    
    return peri_min, peri_ave, peri_max

def disp_xmodel3(p1,p2,p3,h1,h2,h3,L,vertical,m):
    
#    H     = he+hm+hh
    
#    eem   = 1 - pe/pm
#    eeh   = 1 - pe/ph 
#    emh   = 1 - pm/ph
    
#    gamma = eem*he*hm + eeh*he*hh + emh*hm*hh
#    alpha = he*hm*hh*eem*emh
#    delta = ma.pow(gamma,2) - 4*alpha*H

#    lambv1  = 1/(2*H)*(gamma + np.sqrt(delta))     # lambda V1 interfacial
#    lambv2  = 1/(2*H)*(gamma - np.sqrt(delta))     # lambda V2 interfacial
 
    gamma12 = p1/p2
    gamma13 = p1/p3
    gamma23 = p2/p3
    
    A = [[h1, h1,  h1], \
        [h2*gamma12,  h2, h2],\
        [h3*gamma13,  h3*gamma23, h3]]
    
    solv = np.linalg.eigvals(A)
    #sorted(solv)
   
    pv1_min, pv2_min  = eigen3_values(L[0],solv[0],solv[1],m)
    pv1_ave, pv2_ave  = eigen3_values(L[1],solv[0],solv[1],m)
    pv1_max, pv2_max  = eigen3_values(L[2],solv[0],solv[1],m)
    
    if(vertical==1):
        return pv1_min,pv1_ave,pv1_max
    else:   
        return pv2_min,pv2_ave,pv2_max

def disp_xmodel4(p1,p2,p3,p4,h1,h2,h3,h4,L,vertical,m):
    
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

def coriolis_effect(lat,to):
    
    lat     = np.radians(lat)
    omega_e = 2*np.pi/(24*60*60)
    fo      = 2*omega_e*np.sin(lat)
    
    aux = 4*np.power(np.pi,2) + np.power(to,2)*np.power(fo,2)  
    t   = np.sqrt(4*np.power(np.pi,2)*np.power(to,2)/aux)  #second
    
    return t 
    