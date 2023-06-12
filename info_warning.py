# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno

Internalwave Analyzer warning modules (sense module)

modififications to previosly versions must be specified here using WAVE codes:
    
W-23.07-1.00.3-01
A-01.07-1.00.4-01
V-22.07-1.00.4-01
E-05.07-1.00.4-02

"""
def decomp_default(dig,dt):
    dig.write('> Warning: decomposition model has been applied for each time step \n')
    dig.write('> dt='+str(round(int(dt)))+'min \n')
    dig.write('> To run faster, please, increase the time step\n\n\n')

def decomp_changed(dig,dt):
    dig.write('> Warning: the temporal resolution of the decomposition model is higher than the temporal resolution of analysis \n')
    dig.write('> The resolution of the decomposition model has been changed to the same of the resolution data \n')
    dig.write('> dt='+str(round(int(dt)))+'min\n\n\n')
    
def decomp_multiple(dig,dt):
    dig.write('> Warning: the temporal resolution of the decomposition model must be multiple of the temporal resolution of analysis \n')
    dig.write('> The resolution of the decomposition model has been changed to the nearest multiple value\n')
    dig.write('> dt='+str(round(int(dt)))+'min\n\n\n')        
    
def decomp_specified(dig,dt):
    dig.write('> Warning: the temporal resolution of the decomposition was specified by the user \n')
    dig.write('> the temporal resolution is for the decomposition model is specified to '+str(round(int(dt)))+'min \n\n\n')

def coarse_fetch(dig):
    dig.write('> Warning: data from .fet is to coarse to identify the \n')
    dig.write('> variability of the wind fetch \n')
    dig.write('> The wind fetch was estimeted based on the nearest available\n')
    dig.write('> fetch of the mean wind direction considering 95% confidence interval\n\n\n') 

def isotherm_boundary(dig, tau):
    dig.write('> Warning: '+tau+' °C isotherm time series has been deactivated \n')
    dig.write('> it is out of bound  during the whole analyzed period\n')
    dig.write('> Please, spicify a new isotherm that is within the system \n')
    dig.write('> at least during part of the analyzed period\n\n\n') 
    
def three_layer(dig,threrror):
    dig.write('> Warning: Metalimnion borders have not been computed efficiently '+str(int(threrror))+'% of the total analyzed period \n')
    dig.write('> In this case metalimnion thickness is assumed to have a 1/3 thickness of the total water depth \n\n\n')         
  
    
def profile_structure(dig,mode,estimated):
    
    if mode == 1: 
        mode='V1H1'
        stru='two-layer'
    if mode == 2: 
        mode='V2H1'
        stru='three-layer'
    
    dig.write('> Warning: '+mode+' period has been estimated from time averaged profiles  \n')
    dig.write('> '+str(int(estimated))+'% of the total analyzed period \n')
    dig.write('> It may occurs due to the difficulty to Interwave Analyzer identify the '+stru+' structure \n\n\n')         

def merian(dig):
    dig.write('> Warning: Merian equation was not computed correctly\n')
    dig.write('> The Ri will not be filtered considering the duration of the wind event\n\n\n')

def metalimnion(dig,estimated):
    dig.write('> Warning: Metalimnion borders have not been computed '+str(int(estimated))+'% of the total analyzed period \n')
    dig.write('> In this case metalimnion thickness is assumed to be 5% of the total water depth \n')         
    dig.write('> and the metalimnion density, a simple average between epilimnion and hypolimnion \n\n\n') 
    
def thermocline(dig,estimated):
    dig.write('> Warning: Thermocline was estimated using the approximation of the fastest density gradient\n')
    dig.write('> '+str(int(estimated))+'% of the total analyzed period \n')
    dig.write('> It may occurs due to a division by zero in the more sofisticated method  \n\n\n')    
    
def homogeneous_condition(dig):
    dig.write('> Warning: the mean wind direction considering a homogeneous wind event \n')
    dig.write('> under the Spigel and Imberger (1980) conditions has not been computed \n')
    dig.write('> The program could not indentified a homogeneous wind direction\n\n\n')
    
def welch_average(dig,estimated):
    dig.write('> Warning: due to model instability, welch averaging has been changes to 5 days \n\n\n')
    dig.write('> '+str(int(estimated))+'% of the total analyzed period \n')
    dig.write('> It may occurs due to the difficulty to Interwave Analyzer identify the two-layer structure \n\n\n')     

def average(dig):
    dig.write('> Warning: the system did not identify the duration of large wind event \n')
    dig.write('> Interwave Analyzer considered the standard mean of the period \n')
    dig.write('> for parameters averaged by the duration of wind event average \n\n\n')


def bandpass(dig,typ):
    
    if typ == 'sensor':
        dig.write('> Warning: the temperature signals from one or more sensors\n')
        dig.write('> have not been band pass-filtered \n\n\n')            
    if typ == 'isotherm':
        dig.write('> Warning: one or more isotherms have not been band pass-filtered\n\n\n')
        
def spectral(dig,typ):

    if typ == 'sensor':
        dig.write('> Warning: the temperature signals from one or more sensors\n')
        dig.write('> could not be analyzed on spectral methods \n')  
        dig.write('> This problem can arise due to constant array value \n\n\n')            
    if typ == 'isotherm':
        dig.write('> Warning: the temperature signals from one or more sensors\n')
        dig.write('> could not be analyzed on spectral methods \n')  
        dig.write('> This problem can arise due to constant array value \n\n\n')  
    if typ == 'radiation':
        dig.write('> Warning: Radiation is activated but does not vary with time\n')
        dig.write('> Espectral methods cannot be performed \n')  
        dig.write('> Please deactivate solar raidation \n\n\n')  
    if typ == 'wind':
        dig.write('> Warning: Wind does not vary with time\n')
        dig.write('> Espectral methods cannot be performed \n')  
        dig.write('> Please verify wind information \n\n\n')  