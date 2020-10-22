 # -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer ploting function

modififications to previosly versions must be specified here using WAVE codes:
    
W-23.06-1.00.0-00
A-01.06-1.00.3-00
V-22.06-1.00.3-00
E-05.06-1.00.3-00

"""

import math as ma
import numpy as np 
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import internal_module as mod

#from windrose import WindroseAxes
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def equation_sd (H, x, rhoe, rhoh, L, dh):
#
#   Auxiliary Function: Degeneration Regime 1-2  
#   Function to describe the degeneration limit between regime 1 and 2
#   Correlates the dissipative (d) and the steepening (s) timescales    
#     
    nu = 1.1*10**(-6)  
    he = x*H
    
    glin = 9.81*(rhoh - rhoe)/rhoe
    cp   = np.sqrt(glin*(1-x)*he)
    
    aux1 = 2*nu*L/(3*H*cp) 
    aux2 = np.sqrt(np.pi*cp/(2*nu*L))
    aux3 = np.power(1-x,2)
    aux4 = 1/dh
    aux5 = x*(1-2*x)
    
    return aux1*(aux2*aux3 + aux4)/(aux5)
    
def equation_bo (H, x):
#
#   Auxiliar Function: Degeneration Regime for supercritical flow 
#   Function to describe supercritical and subcritical flow degeneration
#   Formation of bores (bo) 
#            
    aux1 = np.power(1-x,2)
    aux2 = np.power(x,3)
    aux3 = np.power(1-x,3)
    
    return np.sqrt(aux1/(aux2 + aux3))

def equation_kh (H, x, dh):
#
#   Auxiliar Function: Degeneration Regime for Kelvin-Helmnholtz billows
#   Function to describe the formation of KH billows (kh)
#                  
    aux1 = 1/x-1
    aux2 = dh/H
    
    return  2*np.sqrt((aux1*aux2))
 
def boundaries (typ, x, y, z,lim):
#
#   Auxiliar Function: Find boundary limit
#   Function to find the upper/lower boundaries
#                  
    if (typ == 'lower'):
        value = min(np.nanmin(x),np.nanmin(y),np.nanmin(z))-0.05
        if(value < lim):
            value = lim
    else:
        value = max(np.nanmax(x),np.nanmax(y),np.nanmax(z))+0.05
        if(value >lim):
            value= lim
            
    return value            

def degeneration(p1,p2,p3,xa,xb,xc,ya,yb,yc,pe,ph,meta,H,ls,ax):
#
#   Plot Function: Averaged degeneration of BSIW (95% c.l.)
#   Graph to show the mean degeneration regime for three analyzed subperiods  
#
    xl = np.arange(0.01, 0.501, 0.001) 
    y1 = equation_sd (H, xl, pe, ph, ls, meta)
    y2 = equation_bo (H, xl)
    y3 = equation_kh (H, xl, meta)
    
    ax.plot(xl, y1, linewidth=1, c='blue', ls='-' )  
    ax.plot(xl, y2, linewidth=1, c='black', ls='-' )
    ax.plot(xl, y3, linewidth=1, c='red', ls='-' )

    cexa, ceya = np.mean(xa), np.mean(ya)  
    cexb, ceyb = np.mean(xb), np.mean(yb)
    cexc, ceyc = np.mean(xc), np.mean(yc)
    
    e1 = mpatches.Ellipse((cexa, ceya), 2*(xa[0]-cexa), 2*(ya[0]-ceya),color='red', label = "P1~"+str(p1[0].day)+'/'+str(p1[0].month)+'-'+str(p1[1].day)+'/'+str(p1[1].month) ,alpha=0.2)
    ax.add_artist(e1)
    
    e2 = mpatches.Ellipse((cexb, ceyb), 2*(xb[0]-cexb), 2*(yb[0]-ceyb),color='navy', label = "P2~"+str(p2[0].day)+'/'+str(p2[0].month)+'-'+str(p2[1].day)+'/'+str(p2[1].month), alpha=0.2)
    ax.add_artist(e2)
    
    e3 = mpatches.Ellipse((cexc, ceyc), 2*(xc[0]-cexc), 2*(yc[0]-ceyc),color='green', label= "P3~"+str(p3[0].day)+'/'+str(p3[0].month)+'-'+str(p3[1].day)+'/'+str(p3[1].month), alpha=0.2)
    ax.add_artist(e3)

    ax.legend(handles=[e1,e2,e3],loc='upper right',prop={'size': 9})
    ax.set_ylabel('1/W')
    ax.set_xlabel(r'$h_e/H$') 
    
    ax.set_xlim([0.1,0.5])
    ax.set_ylim([0.0,2.0])


def degeneration_evolution(p1,p2,p3,x,y,pe,ph,meta,H,ls,zoom,ax):
#
#   Plot Function: Degeneration of BSIW in time
#   Graph to show the degeneration regime along the analyzed period
#
    xl = np.arange(0.01, 0.501, 0.001) 
    y1 = equation_sd (H, xl, pe, ph, ls, meta)
    y2 = equation_bo (H, xl)
    y3 = equation_kh (H, xl, meta)
    
    ax.plot(xl, y1, linewidth=1, c='blue', ls='-' )
    ax.plot(xl, y2, linewidth=1, c='black', ls='-' )
    ax.plot(xl, y3, linewidth=1, c='red', ls='-' )

    window_size = len(x[0])/10
    
    ya_moved = movingaverage(y[0], window_size)
    yb_moved = movingaverage(y[1], window_size)
    yc_moved = movingaverage(y[2], window_size)
    
    xa_moved = movingaverage(x[0], window_size)
    xb_moved = movingaverage(x[1], window_size)
    xc_moved = movingaverage(x[2], window_size)
         
    while (all(ya_moved) == None and all(yb_moved) == None and all(yc_moved)== None) or (window_size == 0):
            
        ya_moved = movingaverage(y[0], window_size)
        yb_moved = movingaverage(y[1], window_size)
        yc_moved = movingaverage(y[2], window_size)  
            
        window_size = window_size - 1 
            
    ax.scatter(xa_moved,ya_moved,marker='o', s=1, color='red', linewidths=1, label='P1')
    ax.scatter(xb_moved,yb_moved,marker='o', s=1, color='navy', linewidths=1, label='P2')
    ax.scatter(xc_moved,yc_moved,marker='o', s=1, color='green', linewidths=1, label='P3')
        
    ax.legend(loc='upper right',prop={'size': 9})
    ax.set_xlabel(r'$h_e/H$') 
    
    minix = boundaries('lower',xa_moved,xb_moved,xc_moved,0.1)
    maxix = boundaries('upper',xa_moved,xb_moved,xc_moved,0.5)
    miniy = boundaries('lower',ya_moved,yb_moved,yc_moved,0.0)
    maxiy = boundaries('upper',ya_moved,yb_moved,yc_moved,2.0)
        
    if (zoom == 'no'):
        
        yd =  (maxiy-miniy)
        xd =  (maxix-minix)

        rect = mpatches.Rectangle( (minix,miniy),xd,yd,linewidth=1, linestyle = '--',fill=None, alpha=1)
        ax.add_patch(rect)
        
        ax.set_ylabel('1/W')
        ax.set_xlim([0.1,0.5])
        ax.set_ylim([0.0,2.0])    
        
    else:
        if window_size ==0: 
            ax.set_xlim([0.1,0.5])
            ax.set_ylim([0.0,2.0])
            
        else:
            try:
                ax.set_xlim([minix,maxix])
            except ValueError:
                ax.set_xlim([0.1,0.5])
            try:
                ax.set_ylim([miniy,maxiy])   
            except ValueError:
                ax.set_ylim([0.0,2.0])                   

def model_plot(freq,col,maxi,ax,typ='log'):
#
#   Auxiliar Plot: Function to define the model results within the PSD plot
#   Generates rectangular boxex on spectrum indicating the BSIW periods
#    
    maxi = maxi*10**3
    d    = freq[0] - freq[2]
    
    if typ == 'log':
        rec    = mpatches.Rectangle((0.0, freq[2]), maxi, d, color=col,alpha=0.5)
        
    else:
        rec    = mpatches.Rectangle((freq[2],0.0), d, maxi, color=col,alpha=0.5)
        
    ax.add_patch(rec)

def significance(sig,wo,c,ax,typ='xax'):
#    
#   Auxiliar Plot: Function to plot the confidence interval for PSD
#   Generates a arrow indicating the significance level of a power spectral density    
#
    if typ == 'xax':
        try:
            ax.plot(sig, wo, linewidth=1, c=c, ls='--')
        except:
            ax.plot(sig[0], wo, linewidth=1, c=c, ls='--')
    else:
        if typ == 'yax_var':
            sig = sig*wo
            
        try:
            ax.plot(wo,sig, linewidth=1, c=c, ls='--')
        except:
            ax.plot(wo,sig[0], linewidth=1, c=c, ls='--')        

def maxi_significance(sig,maximum):
    
    sig_max=maximum
    for i in range(4):
        if len(sig[i]) > 0:
            if sig_max < max(sig[i]):
                sig_max = max(sig[i])
    
    return sig_max
            
   
def psd_iso_conv(tau, freq, psd, v1h1, v2h1, v3h1,n,fo,wo,sig,curve, ax):
#
#   Plot Function: PSD (m²/Hz) and the variance preserving PSD (m²)
#   This shows the PSD in the conventional way (PSD-y) 
#    
    color = ['black','dimgray','salmon','red']
      
    maxima = 0
    minimox = 10**10
    maximox = 10**-3
    

    for i in range(4):
        if(tau[i]!=-999): 
            
            ff = psd[i]
            f  = freq[i]
            s  = sig[i]
            w  = wo[i]
            
            if curve != 'psd':
                ff = ff*f
                significance(s,w,color[i],ax,'yax_var')
            else:
                significance(s,w,color[i],ax,'yax_psd')
                

                      
            ax.plot(f, ff, linewidth=1, c=color[i], ls='-', label=str(tau[i])+"°C")
            
            max_y = valmax(ff, f)
            
            if(max_y > maxima): maxima = max_y
            
            if(10**ma.ceil(ma.log10(min(freq[i]))) < minimox):
                minimox = 10**(ma.ceil(ma.log10(min(freq[i])))-1)
                
            if(10**ma.ceil(ma.log10(max(freq[i]))) > maximox):    
                maximox = 10**(ma.ceil(ma.log10(max(freq[i]))))

    ax.set_xlim([10**-6,maximox])
    
    ax.set_xscale('log')
   
    if curve == 'psd':
        ax.set_yscale('log')
        ax.set_xlabel('frequency (Hz)')
        ax.set_ylabel('PSD isotherms (m²/Hz)') 
    else:
        ax.set_ylabel('variance-preserving PSD (m²)') 
         

    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    
    ax.plot([fo,fo],[0,maxima],linewidth=1, c='black', ls='--')
    ax.plot([n,n],[0,maxima],linewidth=1, c='black', ls='--')
    
    ax.legend(loc='lower left')
    
def psd_iso(tau, freq, psd, v1h1, v2h1, v3h1, largelen,n,fo,wo,sig, ax):
#
#   Plot Function: Power Spectral Density of iostherms variation (m²/Hz)
#   Graph to show PSD of isotherms and model results for the first 3V modes
#    
    color = ['black','dimgray','salmon','red']
      
    maxima = 0
    minimoy = 10**10
    maximoy = 10**-3
    

    for i in range(4):
        if(tau[i]!=-999): 
            
            ff = psd[i]
            f  = freq[i]
            s  = sig[i]
            w  = wo[i]
                      
            if(largelen==1): 
                l0 = int(len(ff)/150)
                ax.plot(smooth(ff, l0), f, linewidth=1, c=color[i], ls='-',label=str(tau[i])+"°C")
            else:
                ax.plot(ff, f, linewidth=1, c=color[i], ls='-', label=str(tau[i])+"°C")
            
            max_x = valmax(ff, f)
            
            if(max_x > maxima): maxima = max_x
            
            if(10**ma.ceil(ma.log10(min(freq[i]))) < minimoy):
                minimoy = 10**(ma.ceil(ma.log10(min(freq[i])))-1)
                
            if(10**ma.ceil(ma.log10(max(freq[i]))) > maximoy):    
                maximoy = 10**(ma.ceil(ma.log10(max(freq[i]))))

            
            significance(s,w,color[i],ax)
            
    maxima=maxima + 100

    ax.set_ylim([10**-6,maximoy])

    ax.set_yscale('log')
    ax.set_xscale('log')    
    
    ax.set_xlabel('PSD isotherms (m²/Hz)') 
    ax.set_ylabel('frequency (Hz)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    model_plot(v1h1,'green',maxima,ax)
    model_plot(v2h1,'blue',maxima,ax)
    model_plot(v3h1,'red',maxima,ax)
    
    maxima = maxi_significance(sig,maxima)
    ax.set_xlim(right=10**ma.ceil(ma.log10(maxima)))
    
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[fo,fo],linewidth=1, c='black', ls='--')
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[n,n],linewidth=1, c='black', ls='--')
    
    ax.legend(loc='lower left')

   

def psd_depth(sen, matrix, freq, psd, v1h1, v2h1, v3h1,largelen,n,fo, wo, sig, ax):
#
#   Plot Function: Power Spectral Density of temperature variation (°C²/Hz)
#   Graph to show PSD of temperature and model results for the first 3V modes
#   
    color = ['black','dimgray','salmon','red']
    
    
    maxima = 0
    minimoy = 10**10
    maximoy = 10**-3
    
    fixed_x = None

    l0 = int(len(psd[0])/150)
    
    for i in range(4):
        
        if(sen[i]==1):
            
            depth = np.nanmean(matrix[i])
            ff    = psd[i]
            f     = freq[i]
            s     = sig[i] 
            w     = wo[i]
            
            if(largelen==1): 
                ax.plot(smooth(ff, l0), f, linewidth=1, c=color[i], ls='-', label=str(round(depth,2))+" m")
            else:
                ax.plot(ff, f, linewidth=1, c=color[i], ls='-', label=str(round(depth,2))+" m")
        
            max_x = valmax(ff,f)
            if(max_x > maxima):
                maxima = max_x
            
            if(10**ma.ceil(ma.log10(min(freq[i]))) < minimoy):
                minimoy = 10**(ma.ceil(ma.log10(min(freq[i])))-1)
                
            if(10**ma.ceil(ma.log10(max(freq[i]))) > maximoy):    
                maximoy = 10**(ma.ceil(ma.log10(max(freq[i]))))

            if fixed_x == None: fixed_x = max(ff)

            significance(s,w,color[i],ax)

    maxima = maxima + 50
    ax.set_ylim([minimoy,maximoy])
 
    ax.set_yscale('log')
    ax.set_xscale('log')    
    
    ax.set_xlabel('PSD isotherms (m²/Hz)') 
    ax.set_ylabel('frequency (Hz)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    model_plot(v1h1,'green',maxima,ax)
    model_plot(v2h1,'blue',maxima,ax)
    model_plot(v3h1,'red',maxima,ax)
    
    maxima = max([maxima,np.nanmax(np.array(sig))])
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[fo,fo],linewidth=1, c='black', ls='--')
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[n,n],linewidth=1, c='black', ls='--')

    ax.set_xlim(right=10**ma.ceil(ma.log10(maxima)))
    ax.legend(loc='lower left')    



def smooth(y, box_pts):
#
#   Auxiliar Function: Smooth data signal
#   Function to smooth time-series
# 
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def smooth_date(t,dt,y,limite=600):
#
#   Auxiliar Function: Smooth signal date
#   Function to smooth the date signal (time)
#   
    new_date = ["" for x in range(limite)]
    
    t_sm     = np.array(dt)
    t_smooth = np.linspace(np.min(t_sm), np.max(t_sm), limite)
    index    = value_nearest_vector(t_sm,t_smooth,limite)
    y_smooth = np.interp(t_smooth, dt, y)
    
    for x in range(limite):
        new_date[x] = t[int(index[x])]
        
    return new_date, t_smooth, y_smooth

def value_nearest_vector(array,value,lim):
#
#   Auxiliar Function: Find nearest value from an array
#   Function to find the index from an array
# 
    idx  = np.zeros(lim,float)
    
    for x in range(lim):
        idx[x] = (np.abs(array-value[x])).argmin()
        
    return idx

def interpolation(y1,y2,x1,x2,x):
#
#   Auxiliar Function: Interpolation function
#   Function to interpolate values from x and y
#
    if(x1==x2):
        return abs(np.mean([y1,y2]))
    
    a = (y1-y2)/(x1-x2)
    b = y1 - a*x1
    
    y = a*x + b
    return y


def correction_contourf(h,num,dj,lin,dci,limit):
#
#   Auxiliar Function: Interpolation function to plot a graph
#   Function to trasnform the elevation 2d array y into a 1d array y'
#
    mx = np.amax(h)                  
    mn = np.amin(h) 

    le = int((ma.ceil(mx)-ma.floor(mn))/dj)                 
    h = np.around(h, decimals=dci)
    nh = np.linspace(ma.ceil(mx),ma.floor(mn), le) 

    countnh = len(nh)
    nnum = np.zeros((lin,le),float)   

    for t in range(lin):    
        i=0
        for z in range(countnh):        
            if nh[z] > h[t][i] or nh[z] < h[t][limit-1]:
                nnum[t][z] = None
            else:
                nnum[t][z] = interpolation(num[t][i],num[t][i+1],h[t][i],h[t][i+1],nh[z])  
            if i < (limit - 1): 
                if nh[z] < h[t][i+1]:
                    i=i+1
    return nnum, nh
    
    
def wavelet_iso(dx,per, power, show, ax):
#
#   Plot Function: Wavelet of isotherm (m²/Hz)
#   Graph to show the wavelet from isotherms data signals
#   
    x = ax.contourf(dx,per,power,100, extend='both', cmap=plt.get_cmap('plasma')) 
    
    ax.set_yscale('log')
    ax.set_ylabel('scale (hours)')
    ax.set_ylim([0,40])
     
    if show == 0: 
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
        plt.colorbar(x, orientation='vertical', label='wavelet (m²/Hz)',aspect=40, ax=ax)

def wavelet_depth(dx,per, power, time, ax):
#
#   Plot Function: Wavelet of isotherm (°C²/Hz)
#   Graph to show the wavelet from isotherms data signals
#  
    x = ax.contourf(dx,per,power,100, extend='both', cmap=plt.get_cmap('plasma')) 
    
    ax.set_yscale('log')
    ax.set_ylabel('scale (hours)')
    ax.set_ylim([0,40])
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    plt.colorbar(x, orientation='vertical', label='wavelet (°C²/Hz)',aspect=40, ax=ax)


def wavelet_resona(dx,per, power, v1, v2, ax):
#
#   Plot Function: Wavelet of wind with model results from seiche modes
#   Graph to show the internal sieche resonance
# 
    x = ax.contourf(dx,per,power,100, extend='both', cmap=plt.get_cmap('plasma')) 
    
    v1 = v1/60/60 
    v2 = v2/60/60
    
    ax.plot(dx, v1, linewidth=1, c='black', ls='-', label='V1H1')
    ax.plot(dx, v2, linewidth=1, c='black', ls='--', label='V2H1')

    ax.legend(loc='upper right', prop={'size': 9})
    
    ax.set_yscale('log')
    ax.set_ylabel('scale (hours)')
    ax.set_ylim([0,40])
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    plt.colorbar(x, orientation='vertical', label='wavelet (m²/(s² Hz))',aspect=40, ax=ax)


def isotherm(date,t,y,color,iso,ax):
#
#   Plot Function: Isotherm variation
#   Graph to show the temporal variation of the analyzed isotherm
# 
    ax.plot(date, y, linewidth=1, c=color, ls='-')
    
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    ax.set_xlim([date[0],date[-1]])
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_ylabel('isotherm '+str(iso)+ ' °C (m)')
    

def temperature(date,t,y,color,depth,ax):
#
#   Plot Function: Temperature variation
#   Graph to show the temporal variation of temperature
#        
    ax.plot(date, y, linewidth=1, c=color, ls='-')
    
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    ax.set_xlim([date[0],date[-1]])
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_ylabel('temp. variation at '+str(round(np.nanmean(depth),2))+' m (°C)')


def multi_isotherms(date, y, iso, time, z0, ax):
#
#   Plot Function: iostherm variation
#   Graph to show the temporal variation of all analyzed isotherms
# 
    color = ['black','dimgray','salmon','red']

    for i in range(4):
        if iso[i]!=-999:

            ax.plot(date, y[i], linewidth=1, c=color[i], ls='-', label = str(iso[i])+'°C')
           
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
    ax.set_xlim([date[0],date[-1]])
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.legend(loc='upper right', prop={'size': 8})
    
    if z0 == 0:
        ax.set_ylabel('height above bottom (m)')   
    else:
        ax.set_ylabel('altitude (m)')
        ax.set_title( 'reference level for bottom: '+str(z0)+' m',loc='left', fontsize= 9)


def temp_bandpass(date, y, iso, ax):
#
#   Plot Function: bandpass filtered isotherms
#   Graph to show the bandpass filter applied on isotherms
#     
    color = ['black','dimgray','salmon','red']
    
    for i in range(4):
        if iso[i]!=-999:
            try:
                ax.plot(date, y[i], linewidth=1, c=color[i], ls='-', label=str(iso[i])+'°C')
            except ValueError:
                pass

    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    

    ax.set_xlim([date[0],date[-1]])
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.legend( loc='upper right', prop={'size': 9})
    ax.set_ylabel('displacement (m)')  


def depth_bandpass(date, y, depth, time,s, ax):
#
#   Plot Function: bandpass filtered thermistors sensors
#   Graph to show the bandpasss filter applied on temperature variation
#        
    color = ['red','maroon','blue','navy']
    
    for i in range(4):
        if(s[i] == 1):
            try:
                ide = np.nanmean(depth[i])
                ax.plot(date, y[i], linewidth=1, c=color[i], ls='-', label=str(round(ide,0))+' m')
            except ValueError:
                pass  
            
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))

    ax.set_xlim([date[0],date[-1]])
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.legend( loc='upper right', prop={'size': 8})
    ax.set_ylabel('thermal variation (°C)')  


def wind_rose(ws,wd):
#
#   Plot Function: Windrose 
#   Graph to show windrose of the analyzed period    
#    
    wd = wd + 180
    wd = np.where(wd > 360, wd-360, wd) # nautica definition to wind direction
    
    
    #ax = WindroseAxes.from_ax()
    #rounda = int(max(ws))+1
    #bins_range = np.linspace(0,rounda,7)
    #ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white', bins=bins_range)
    #ax.set_legend(fontsize=15)

   
def averaged_profile(temp,low,high,h,z0, ax):
#
#   Plot Function: Averaged tempearture profile
#   Graph to show temporal variation of temperature recorded from thermistors   
#   This analysis ignores surface elevation variation
#    
    
    ax.scatter(temp, h, color='red', label='time averaged profile') 

    ax.plot([low,high],[h,h],c='blue',lw=1)
        
    
    ax.plot(temp, h, linewidth=1, c='black', ls='-')      
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    ax.set_xlabel('water temperature (°C)')  
    
    if z0 == 0:
        ax.set_ylabel('height above bottom (m)')   
    else:
        ax.set_ylabel('altitude (m)')
        ax.set_title( 'reference level for bottom: '+str(z0)+' m',loc='right', fontsize= 9)

    
    ax.legend(loc='upper left',prop={'size': 9})
 


def thermal_variation(date,depth, depth_point, temp, time, ax):
#
#   Plot Function: Thermal fluctuation 
#   Graph to show the time average temperature profile 
#    
    ax.plot(date, temp, linewidth=1, c='black', ls='-')
    
    col = ['red','maroon','blue','navy']
    for i in range(4):  
        aux = depth[i]
        if aux[0] is not None:
            
            d = temp[:,depth_point[i]]
            ax.plot(date, d, linewidth=1, c=col[i], ls='-', label=str(round(np.nanmean(depth[i]),1))+' m') 
      
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
    ax.set_xlim([date[0],date[-1]])
    ax.set_ylabel('water temperature (°C)')  
    
    ax.legend(loc='upper right',prop={'size': 9})
    
    
def thermocline(date,t,y,color,thermo,ax):
#
#   Plot Function: Thermocline temperature fluctuation 
#   Graph to show temporal variation of temperature at thermocline depth
#         
    ax.plot(date, y, linewidth=1, c=color, ls='-')
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    ax.set_xlim([date[0],date[-1]])
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_ylabel(str(thermo)+ ' m (°C)')


def psd_thermocline(p,psd,color,depth,wo,sig,ax):
#
#   Plot Function: PSD of thermocline temperature
#   Graph to show PSD of temperature fluctuation at thermocline depth
#   
    ax.plot(psd,p, linewidth=1, c=color, ls='-')
    significance(sig,(1/wo)/3600,'red',ax,typ='xax')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_xlabel('PSD '+ str(depth) +' m (°C²/Hz)')         

def psd_multilayer(f,depth,ff,fi,n2,z0,ax):
#
#   Plot Function: Multi-layer PSD of temperature fluctuation
#   Graph to show PSD of temperature fluctuation along water depth 
#   This analysis ignore surface elevation variation
#     
 
    ff = np.array(ff)
    cf = ax.contourf(f, depth, np.log10(ff[::-1,:]), 100)
    plt.colorbar(cf,orientation='vertical',label='PSD (log(°C²/Hz))',ax=ax)
    
    ax.plot([fi,fi],[min(depth),max(depth)],linewidth=1, color='black', ls='--')
 
    if max(f) > 10**-2:
        ax.plot([n2,n2],[min(depth),max(depth)],linewidth=1, color='black', ls='--')
        
    ax.set_xscale('log')
    ax.set_xlim(left=10**-6)
    ax.set_xlabel('frequency (Hz)')

    if z0 == 0:
        ax.set_ylabel('height above bottom (m)')   
    else:
        ax.set_ylabel('altitude (m)')
        ax.set_title( 'reference level for bottom: '+str(z0)+' m',loc='left',fontsize= 9)



def valmax(psd,f):
#
#   Auxiliar Function: Find maximum spectral energy
#   Function to find the maximum energy of the analyzed spectrum
#
    maxi = 0    
    l    = len(psd)
    
    for i in range(l):    
        if(f[i]>6*10**-6 and psd[i]>maxi):
            maxi = psd[i]
            
    return maxi
        

def psd_wind(period,psd_wind,wo,sig,ax):
#
#   Plot Function: PSD of wind speed
#   Graph to show PSD of wind speed
#   
    ax.plot(period[1:], psd_wind[1:],  linewidth=1, c='navy', ls='-')
    significance(sig,(1/wo)/3600,'navy',ax,typ='yax')
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    ax.tick_params(axis='y', labelcolor='navy' )
    ax.set_ylabel('PSD wind ((m/s)²/Hz)', color='navy') 
    ax.set_xlabel('periodicity (hours)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)

    
def psd_sola(period,psd_sola,wo,sig,ax):
#
#   Plot Function: PSD of solar radiation
#   Graph to show PSD of solar radiation 
#  
    ax.plot(period,psd_sola, linewidth=1, c='red', ls='-')
    significance(sig,(1/wo)/3600,'red',ax,typ='yax')
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    ax.tick_params(axis='y', labelcolor='red')
    
    ax.set_ylabel('PSD solar radiation ((W/m²)²/Hz)', color='red') 
    ax.set_xlabel('periodicity (hours)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)

def psd_level(f, ff, wo, sig, ax):
#
#   Plot Function: PSD of water elevation
#   Graph to show PSD of water elevation
#  
    ax.plot(f, ff, linewidth=1, c='black', ls='-')
    significance(sig,wo,'red',ax,typ='yax')
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    ax.set_ylabel('PSD water elevation (m²/Hz)') 
    ax.set_xlabel('frequency (Hz)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)

def coherence (c1w,f1w,tau,iso, largelen,ax):
#
#   Plot Function: Coherence between isotherms and wind speed
#   Graph to plot the coherence between wind and isotherms
#        
    l10 = int(len(c1w)/100)
    
    if largelen==1:
        ax.plot(smooth(c1w,l10), f1w, linewidth=1, c='navy', ls='-', label='cohe wind')

    else:
        ax.plot(c1w, f1w, linewidth=1, c='navy', ls='-', label='cohe wind')
    
    ax.set_xlim([0,1])
    ax.set_xlabel('coherence iso'+ str(tau[iso]) + '/' + 'met. data')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    ax.legend(loc='lower right')


def code_index(r):
#
#   Auxiliar Function: Code index
#   Auxiliar function to define the combination of isotherms
# 
    if r == 2:
        return 0
    elif r == 3:
        return 1
    elif r == 4:
        return 2
    elif r == 6:
        return 3
    elif r == 8:
        return 4
    elif r == 12:
        return 5

def coherence_iso(t,coh,fre,ana1, ana2, ax):
#
#   Plot Function: Coherence between isotherms
#   Graph to plot the coherence between isotherms
#   
    c = ['black','dimgray','salmon','red','navy','blue']

    ax.grid(axis='y',color='black',ls=":",lw=0.25)
    ax.set_ylim([0,1.01])

    for i in range(4):
        for j in range(4):
            if (i!=j and j>i):
                
                a1 = np.around(float(ana1)/10,0)
                b1 = ana1 - 10*a1
                a2 = np.around(float(ana2)/10,0)
                b2 = ana2 - 10*a2
                
                if((a1 == i+1 and b1 == j+1) or (a2 == i+1 and b2 == j+1)):
                    if(t[i]!=-999 and t[j]!=-999):
                        
                        k = code_index((i+1)*(j+1))
                        ax.plot((1/fre[k])/60/60,coh[k],linewidth=1, ls='-', color=c[k], label=str(t[i])+'/'+str(t[j])+'°C')

    ax.set_ylabel('coherence', color='gray')
    ax.set_xlabel('period (hours)')
    

def coherence_depth(matrix_depth,coh,fre,ax):
#
#   Plot Function: Coherence between thermistors
#   Graph to plot the coherence between temperature fluctuation from sensors
#  
    c = ['red', 'darkgreen', 'darkcyan']
    
    for i in range(3):
        if coh[i] is not None:
            depth1 = np.nanmean(matrix_depth[i])
            depth2 = np.nanmean(matrix_depth[i+1])
            ax.plot(1/fre[i],coh[i],linewidth=1, color=c[i], label=str(round(depth1,1))+'/'+str(round(depth2,1))+' m')

    ax.set_xlim([0,25])
    ax.set_ylabel('coherence',  color='gray')
    ax.set_xlabel('period (hours)')

def phase_iso(t, ph, fre, conf, ana1,ana2, ax):
#
#   Plot Function: Phase-shift between isotherms
#   Graph to plot the phase-shift between isotherms
# 
    c = ['black','dimgray','salmon','red','navy','blue']
    ax.grid(axis='y',color='black',ls=":",lw=0.25)
    ax.set_ylim([0,1.01])

    for i in range(4):
        for j in range(4):
            if (i!=j and j>i):
                a1 = np.around(float(ana1)/10,0)
                b1 = ana1 - 10*a1
                a2 = np.around(float(ana2)/10,0)
                b2 = ana2 - 10*a2
                if((a1 == i+1 and b1 == j+1) or (a2 == i+1 and b2 == j+1)):
                    if(t[i]!=-999 and t[j]!=-999):
                        k = code_index((i+1)*(j+1))
                        f = fre[k]
                        p = ph[k]
                        ax.scatter((1/f[conf[k]])/60/60,abs(p[conf[k]]),marker='x', color=c[k], label=str(t[i])+'/'+str(t[j])+'°C')


    ax.set_ylabel('Phase (degrees)', color='black')
    ax.set_yticks([ 0, np.pi/2, np.pi])
    ax.set_yticklabels([ r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])


    ax.grid(axis='y',color='black',ls=":",lw=0.25)
    ax.legend(loc='upper right',prop={'size': 6})
    ax.set_ylabel('phase (radian)')


def phase_depth(matrix_depth,conf, ph,fr,ax):
#
#   Plot Function: Phase-shift between temperature signals 
#   Graph to plot the phase-shift between temperature signals
# 
    c = ['red', 'darkgreen', 'darkcyan']
    
    for i in range(3):
        if(conf[i] is not None):    
            
            depth1 = np.nanmean(matrix_depth[i])
            depth2 = np.nanmean(matrix_depth[i+1])
            
            f=fr[i]
            p=ph[i]

            ax.scatter(1/f[conf[i]],abs(p[conf[i]]),linewidth=1, marker='o', color=c[i], label=str(round(depth1,1))+'/'+str(round(depth2,1))+' m')
                    
    ax.set_ylabel('Phase (degrees)', color='black')
    ax.set_yticks([ 0, np.pi/2, np.pi])
    ax.set_yticklabels([ r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])


    ax.grid(axis='y',color='black',ls=":",lw=0.25)
    ax.legend(loc='upper right',prop={'size': 6})
    ax.set_ylabel('phase (radian)')


def density(dx,t,pu,pd,pe,ph,ax):
#
#   Plot Function: Water density structure variation
#   Graph to plot the temporal variation of water density (mean and boudanries)
#
    ax.plot(dx, pu, linewidth=1, c='darkred', ls='-')
    ax.plot(dx, pd, linewidth=1, c='darkblue', ls='-')

    ax.plot(dx, pe, linewidth=1, c='red', ls='--')
    ax.plot(dx, ph, linewidth=1, c='blue', ls='--')    
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    ax.invert_yaxis()
    
    ax.set_ylabel(r'$\rho_{water} $'+' (kg/m³)')
    ax.legend(['superficial','bottom','epilimnion','hypolimnion'], loc='upper right',prop={'size': 8})

def radiation(date,t,y,ax):
#
#   Plot Function: Solar radiation
#   Graph to plot solar radiation signal
#    
    ax.plot(date, y, linewidth=1, c='darkgoldenrod', ls='-')
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    ax.set_ylabel('solar rad.(W/m²)')
    
def wind_direction(date,time,y,ys,ylow,ydi,wedd_up,wedd_low,ax):
#
#   Plot Function: Wind direction filtered by Wedderburn number  
#   Graph to plot wind direction signal (considering continuous low W number)
# 
    ax.plot(date, ydi, c='blue', linewidth=3, ls='-', label ='homogeneous wind events')
    
    ax.plot(date, y, linewidth=1, c='lightgray', ls='-')
    ax.plot(date, ys, linewidth=1, c='green', ls='-',label='20 < W < '+str(round(wedd_low,1)))
    ax.plot(date, ylow, linewidth=1, c='red', ls='-', label=str(round(min([wedd_up,1]),1))+'< W < 20')
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.set_ylabel('wind direction '+r'$(^{\circ})$',color='gray') 
    
    ax.legend(loc='upper right', prop={'size': 8})
    ax.set_xlim([date[0],date[-1]])
    ax.set_ylim([0,360])
     
    
def wind(date,time,dire,velo,ax):
#
#   Plot Function: Wind speed
#   Graph to plot wind speed signal 
# 
    l  = len(date)
    dt = time[1]-time[0]

    ax.plot(date, velo, linewidth=1, c='navy', ls='-')

    ax.set_xlim([date[0],date[-1]])
    ax.set_ylim([0,int(2*max(velo))])
    
    ax.tick_params(axis='y', labelcolor='navy' )
    ax.set_ylabel('wind speed (m/s)',color='navy')  
    
    t_smooth, t_numeral, y_smooth = smooth_date(date,time,dire,int(l*dt))
    
    
    ax2 = ax.twinx() 
    ax2.plot(t_smooth, y_smooth, linewidth=1, c='black', ls='-')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax2.tick_params(axis='y', labelcolor='black' )
    ax2.set_ylabel('wind direction '+r'$(^{\circ})$',color='black') 
    ax2.set_ylim([-200,380])
    ax2.set_yticks([0, 180, 360])
    

def windstress (date,time, y, ax):
#
#   Plot Function: Windstress
#   Graph to plot windstress signal 
# 
    ax.plot(date, y, linewidth=1, c='navy', ls='-')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2.00))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    
    ax.set_ylim(bottom=0)
    ax.set_ylabel('wind stress (N/m²)') 


def wedderburn (date,t,y,ax):
#
#   Plot Function: Wedderburn number at thermocline depth
#   Graph to plot Wedderburn number signal at thermocline depth (2-layer system)
#     
    ax.plot(date, y, linewidth=1, c='navy', ls='-')
   
    ax.set_yscale('log')
    ax.set_ylabel('Wedd number') 
    
    ax.set_xlim([date[0],date[-1]])
    
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    

def wedd_limit(dx,time,ri,dw,up,ax):
#
#   Plot Function: Wedderburn number at thermocline (internal wave regime)
#   Graph to plot Ri at thermocline highlighting the occurance of BSIW 
#         
    ax.plot(dx, dw, linewidth=1, c='red', ls='-')
    ax.plot(dx, up, linewidth=1, c='red', ls='-')
    ax.plot(dx, ri, linewidth=1, c='navy', ls='-')
   
    ax.fill_between(dx, dw, up, where=up >= dw, facecolor='red', alpha=0.5, interpolate=True)    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)

    l= len(dx)
    ax.set_xlim([dx[0],dx[l-1]])

    ax.set_yscale('log')
    ax.set_ylabel(r'$W$')

def wedd_compb(dx,time,ri,rid,dw,up,ax):
#
#   Plot Function: Filtered Richardson number at thermocline
#   Graph to plot filtered Ri considering consectutive/aligned wind events 
#         
    ax.plot(dx, dw, linewidth=1, c='red', ls='-')
    ax.plot(dx, up, linewidth=1, c='red', ls='-')
    
    ax.plot(dx, ri, linewidth=1, c='navy', ls='-', label= 'standard W')
    ax.plot(dx, rid, linewidth=1, c='dodgerblue', ls='-', label='filtered W')
  
    ax.fill_between(dx, dw, up, where=up >= dw, facecolor='red', alpha=0.5, interpolate=True)
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.legend(loc='upper right', prop={'size': 9})

    l= len(dx)
    ax.set_xlim([dx[0],dx[l-1]])
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_yscale('log')
    ax.set_ylabel(r'$W$')


def richardson2d(dx,time,h,ri,thermo,cota,z0,ax):
#
#   Plot Function: 2D Richardson number
#   Graph to plot filtered Ri considering consectutive/aligned wind events 
#    
    x = ax.contourf(dx, h, ri.transpose(), 100, locator=ticker.LogLocator())
    ax.plot(dx, thermo, linewidth=1, c='red', ls='-')
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')

    ax.set_xticklabels(dx, rotation=10)
    
    ax.set_ylim([z0,ma.trunc(np.max(cota))+1])
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    plt.colorbar(x, orientation='horizontal', label='Ri (-)')
    
    if z0 == 0:
        ax.set_ylabel('height above bottom (m)')   
    else:
        ax.set_ylabel('altitude (m)')
        ax.set_title( 'reference level for bottom: '+str(z0)+' m',loc='left',fontsize= 9)


    
def tempstructure(dx, time,h,thermo,temp,cota,z0,ax):
#
#   Plot Function: Temperature structure (time-serie)
#   Graph to plot temperature structure 
#     
    x = ax.contourf(dx,h,temp.transpose(), 300)
    
    ax.plot(dx,thermo,linewidth=1, c='w', ls='-')
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')
    
    ax.set_ylim([z0,ma.trunc(np.max(cota))+1])
    
    ax.set_xticklabels(dx, rotation=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
    plt.colorbar(x, orientation='horizontal', label='Temperature (°C)')
    
    if z0 == 0:
        ax.set_ylabel('height above bottom (m)')   
    else:
        ax.set_ylabel('altitude (m)')
        ax.set_title( 'reference level for bottom: '+str(z0)+' m',loc='left',fontsize= 9)


def tempstructure_paper(dx, iso1,iso2, time,h,thermo,temp,cota,z0,ax):                                         # retirar (depois)
#
#   Plot Function: Temperature structure (time-serie)
#   Graph to plot temperature structure 
#     
    x = ax.contourf(dx,h,temp.transpose(), 300)
    
    ax.plot(dx,iso1,linewidth=1, c='black', ls='-')
    ax.plot(dx,iso2,linewidth=1, c='red', ls='-')
    
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')
    
    ax.set_ylim([z0,ma.trunc(np.max(cota))+1])
    
    ax.set_xticklabels(dx, rotation=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
    plt.colorbar(x, orientation='horizontal', label='Temperature (°C)')
    
    if z0 == 0:
        ax.set_ylabel('height above bottom (m)')   
    else:
        ax.set_ylabel('altitude (m)')
        ax.set_title( 'reference level for bottom: '+str(z0)+' m',loc='left',fontsize= 9)

    
def tempstructure_zoom(dx, time,h,thermo,temp,cota, z0, ax):
#
#   Plot Function: Temperature structure near the thermocline (time-serie)
#   Graph to plot temperature structure highlighting the upper layer 
#      
    x = ax.contourf(dx,h,temp.transpose(), 300)

    ax.plot(dx,thermo,linewidth=1, c='w', ls='-')
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')
    
    ax.set_ylim([ma.trunc(np.min(thermo))-1,ma.trunc(np.max(cota))+1])
    
    ax.set_xticklabels(dx, rotation=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
      
    plt.colorbar(x, orientation='horizontal', label='Temperature (°C)')
    
    if z0 == 0:
        ax.set_ylabel('height above bottom (m)')   
    else:
        ax.set_ylabel('altitude (m)')
        ax.set_title( 'reference level for bottom: '+str(z0)+' m',loc='left',fontsize= 9)


def thorpe_scale(dx, time,h,thorpe,cota,z0,ax):
#
#   Plot Function: Thorpe scale (time-serie)
#   Graph to plot Thorpe-scale displacement 
#      
    x = ax.contourf(dx,h,thorpe.transpose(), 300)
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')
    
    ax.set_ylim([z0,ma.trunc(np.max(cota))+1])
    
    ax.set_xticklabels(dx, rotation=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
    plt.colorbar(x, orientation='horizontal', label='Thorpe dispalcement (m)')  

    if z0 == 0:
        ax.set_ylabel('height above bottom (m)')   
    else:
        ax.set_ylabel('altitude (m)')
        ax.set_title( 'reference level for bottom: '+str(z0)+' m',loc='left',fontsize= 9)


def movingaverage(interval, window_size):
#
#   Auxiliar Function: Convolution (moving average)
#   Function to calculate the moving average of a time-serie
#       
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')
   

def generation(cla,ax):
#
#   Plot Function: Lake mixing classification
#   Identify the temporal percentage of lake mixing regimes 
#   
    N       = len(cla)
    x       = range(N)
    width   = 1/1.5
    
    ax.bar(x,cla,width,color="navy")
    ax.set_ylim([0,100])
    ax.set_xticklabels(("-","Mixed","KH","IS/d","IS/s"))

def sensibility_metalimnion(delta,period,m,ax):
#
#   Plot Function: Sensibility analysis for layer thickness variation 
#   Graph to plot the BSIW period variation based on metalimnion thrshold (V2 mode) 
#  

    ax.plot(delta, period, linewidth=1, c='black', ls='-')      
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    ax.set_xlabel('metalimnion threshold (kg/m²)')  
    ax.set_ylabel('period (hours)')
    
    ax.set_title(r'$ \delta = $ '+str(m)+' kg/m²',loc='right', fontsize= 9)


def depth_sensitivity(xh,per,depth,labelx,ax):
#
#   Plot Function: Sensibility analysis for layer thickness variation 
#   Graph to plot the BSIW period variation based on layer thickness variation 
#       
    font = {'family': 'sans-serif',
    'color':  'black',
    'weight': 'normal',
    'size': 10,
    }
    
    ax.plot(xh, per, linewidth=1, c='navy', ls='-') 
    ax.text((np.nanmax(xh)+np.nanmin(xh))/2, (np.nanmax(per)+np.nanmin(per))/2, r'$h^* \equiv$'+str(depth), fontdict=font)
     
    ax.set_ylabel('period (h)')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    if(labelx == 1):
        ax.set_xlabel(r'$h^*(m)$') 
        
        
def densi_sensitivity(xp,per,rho,labelx,ax):
#
#   Plot Function: Sensibility analysis for water density variation 
#   Graph to plot the BSIW period variation based on water density variation 
#      
    font = {'family': 'sans-serif',
    'color':  'black',
    'weight': 'normal',
    'size': 10,
    }
    
    ax.plot(xp, per, linewidth=1, c='navy', ls='-')
    ax.text((np.nanmax(xp)+np.nanmin(xp))/2, (np.nanmax(per)+np.nanmin(per))/2, r'$\rho^* \equiv$'+str(rho), fontdict=font)
 
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    if(labelx==1):
        ax.set_xlabel(r'$\rho^* (kg/m³)$') 
    
    
def parabar_sensitivity(x,per,typ,limin1,lim2,limax1,ax):
#
#   Plot Function: Sensibility analysis considering general parameter
#   Graph to plot the BSIW period variation based on general parameter variation 
#        
    limin1 = limin1/60/60 # in hour
    limax1 = limax1/60/60
    
    ax.plot(x, per, linewidth=1, c='navy', ls='-')
    
    del1 = limax1 - limin1
    rec1 = mpatches.Rectangle((min(x), limin1), max(x), del1, color='green',alpha=0.5)
    ax.add_patch(rec1)
    
    z  = 2.58
    s = np.nanstd(lim2)
    nanlen = len(lim2) - np.isnan(np.array(lim2, dtype=np.float64)).sum()
    cl = (z * (s/np.sqrt(nanlen)))
    
    limin2 = np.nanmean(lim2) - cl
    del2 = 2*cl

    rec2 = mpatches.Rectangle((limin2, np.nanmin(per)), del2, np.nanmax(per), color='gray',alpha=0.5)
    ax.add_patch(rec2)

    ax.set_xlim([min(np.nanmin(x),limin2), max(np.nanmax(x),limin2 + del2)])
    ax.set_ylim([np.nanmin(per),np.nanmax(per)])
    
    if typ == 'dep':
        ax.set_ylabel('period (h)')
        ax.set_xlabel(r'$\bar{H}^{-1/2}$')
    elif typ == 'rho':
        ax.set_xlabel(r'$ (\rho_h/\Delta\rho)^{1/2}$')

def classification_genera(W1, W2, heh_found,ahe,tau, ax):
#
#   Plot Function: Basin-scale internal wave amplitude (theoretical result)
#   Graph to plot the amplitude of BSIW
#  
    try:
        label=['Spigel and Imberger (1980)', ' de Carvalho Bueno et al. (2020) '+r'$h_e/H=0.7$',' de Carvalho Bueno et al. (2020) '+r'$h_e/H=0.5$', ' de Carvalho Bueno et al. (2020) '+r'$h_e/H=0.2$',r'$W_{V1H1} $'+str(tau[0])+ '°C', r'$W_{wind} $'+str(tau[0])+ '°C']
    except IndexError:
        label=['Spigel and Imberger (1980)', ' de Carvalho Bueno et al. (2020) '+r'$h_e/H=0.7$',' de Carvalho Bueno et al. (2020) '+r'$h_e/H=0.5$', ' de Carvalho Bueno et al. (2020) '+r'$h_e/H=0.2$',r'$W_{V1H1} $', r'$W_{wind} $']
    
    wedderburn = np.linspace(0.1, 5,   num=1000, endpoint=True)
    wedd1      = np.linspace(0.1, 5, num=1000, endpoint=True)
    wedd2      = np.linspace(5, 1000, num=1000, endpoint=True)
    wedd       = np.concatenate((wedd1,wedd2))
    
    linear = 0.5/wedd
    curve1 = mod.bueno(10.5, 4.5, wedderburn,wedderburn)/10.5
    curve2 = mod.bueno(7.5, 7.5, wedderburn,wedderburn)/7.5
    curve3 = mod.bueno(3, 12, wedderburn,wedderburn)/3
    
    ax.plot(wedd[1:-1],linear[1:-1],linewidth=1, ls='--')
    ax.plot(wedderburn,curve1,linewidth=1, ls='-')
    ax.plot(wedderburn,curve2,linewidth=1, ls='-')
    ax.plot(wedderburn,curve3,linewidth=1, ls='-')
    
    ax.scatter(W1,ahe,s=20,marker='o',color='red')
    ax.scatter(W2,ahe,s=20,marker='x',color='navy')


    ax.legend(labels=label,loc='upper right',prop={'size': 8})
    
    ax.grid(axis='y',color='black',ls=":",lw=0.25)
    ax.grid(axis='x',color='black',ls=":",lw=0.25)

    ax.set_ylabel(r'$\zeta_o/h_e$'+' (-)', color='black')
    ax.set_xlabel(r'$Ri h_e/L (-)$', color='black')

    ax.set_xscale('log')
    ax.set_xlim([0.1,1000])
    ax.set_ylim([0,2])
    
def bueno_parameterization(W1, W2, norma, tau, ax):
#
#   Plot Function: Parameterization to predict BSIW amplitude (Bueno, 2020)
#   Graph to plot the amplitude of BSIW according to Bueno parameterization
#  
    try:
        label=['de Carvalho Bueno at al. (2020)', r'$W_{V1H1} $'+str(tau[0])+ '°C', r'$W_{wind} $'+str(tau[0])+ '°C']
    except IndexError:
        label=['de Carvalho Bueno at al. (2020)', r'$W_{V1H1} $', r'$W_{wind} $']

    wedderburn = np.linspace(0.1, 5, num=1000, endpoint=True)
    norm       = (wedderburn-6)**2
    
    ax.plot(wedderburn,norm,linewidth=1,color='black', ls=':')
    ax.scatter(W1,norma,s=20,marker='o',color='red')
    ax.scatter(W2,norma,s=20,marker='x',color='navy')

    ax.legend(labels=label,loc='lower left',prop={'size': 8})

    ax.grid(axis='y',color='black',ls=":",lw=0.25)
    ax.grid(axis='x',color='black',ls=":",lw=0.25)

    ax.set_ylabel(r'$2 f^2 log(\zeta_o/ (0.1~h_e))$'+' (-)', color='black')
    ax.set_xlabel(r'$Ri h_e/L (-)$', color='black')

    ax.set_xscale('log')
    ax.set_xlim([0.1,5])
