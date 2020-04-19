 # -*- coding: utf-8 -*-
"""
@ Created on Sun Jul 23 20:08:35 2017
@ Author: BUENO, R. (RAFAEL BUENO)

Modified:
"""

import math as ma
import numpy as np 
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches

from windrose import WindroseAxes
from scipy.interpolate import BSpline
from matplotlib import pyplot as plt
from internal_module import interpol
from matplotlib.ticker import FormatStrFormatter


# ---------------equations for degeneration of internal waves------------------


def equation_one (H, x, rhoe, rhoh, L, dh):
    
    nu = 1.1*10**(-6)  # water cinematic viscosity
    he = x*H
    
    glin = 9.81*(rhoh - rhoe)/rhoe
    cp   = np.sqrt(glin*(1-x)*he)
    
    aux1 = 2*nu*L/(3*H*cp) 
    aux2 = np.sqrt(np.pi*cp/(2*nu*L))
    aux3 = np.power(1-x,2)
    aux4 = 1/dh
    aux5 = x*(1-2*x)
    
    resp = aux1*(aux2*aux3 + aux4)/(aux5)
    
    return resp

def equation_two (H, x):
    
    
    aux1 = np.power(1-x,2)
    aux2 = np.power(x,3)
    aux3 = np.power(1-x,3)
    
    resp = np.sqrt(aux1/(aux2 + aux3))
    
    return resp

def equation_three (H, x, dh):
         
    aux1 = 1/x-1
    aux2 = dh/H
    
    resp = np.sqrt(2*(aux1*aux2))
    
    return resp

#-----------------------------------------------------------------------------

def model_plot(freq_mode,col,typ,vert,hori,maxima,ax):
    
    # freq_mode - define the weidth of the box-plot
    # psd1      - define the height of the box-plot 
    # col       - define the color of the box-plot
    # typ       - define the type model
    # vert      - mode of the verical mode
    # hori      - mode of the horizontal mode
    
    font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8,
        }
    
    maxima = maxima*10000
    delta  = freq_mode[0] - freq_mode[2]
    rec = mpatches.Rectangle((0.0, freq_mode[2]), maxima, delta, color=col,alpha=0.5)
    
    ax.add_patch(rec)
    
    ph, pm, eh, em = period_err(freq_mode[1],freq_mode[2])  
    ax.text(0.30*maxima, freq_mode[0],'V'+ str(vert) +'H'+ str(hori), fontdict=font)


def period_err(freq,fra_min):
    
    per = (1/freq)/(60*60)  # hours
    
    ph = ma.trunc(per)                              # hour
    pm = ma.trunc((per-ph)*60)                      # minute
    
    delta = (1/fra_min)/(60*60) - per
    
    eh = ma.trunc(delta)                             # hour
    em = ma.trunc((delta-eh)*60)                      # minute
    
    return ph,pm,eh,em

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def smooth_date(t,dt,y,limite=600):
    #
    # t  = timedate form 
    # dt = second (integer) form 
    #
    
    new_date = ["" for x in range(limite)]
    
    t_sm = np.array(dt)
    t_smooth = np.linspace(np.min(t_sm), np.max(t_sm), limite)
    index = value_nearest_vector(t_sm,t_smooth,limite)
    y_smooth = np.interp(t_smooth, dt, y)
    
    for x in range(limite):
        new_date[x] = t[int(index[x])]
        

    return new_date, t_smooth, y_smooth

def value_nearest_vector(array,value,lim):
    # find the nearest value of a 'value' from a 'array'
    # used to create a vector of index
    idx  = np.zeros(lim,float)
    for x in range(lim):
        idx[x] = (np.abs(array-value[x])).argmin()
        
    return idx

def correction_contourf(h,num,dj,lin,dci,limit):
    
    # interpolation function to plot a countourf graph 
    # trasnforms the elevation 2d array y into a 1d array y'
    #
    # input:
    #       dj   : significance of the elevation class (0,1 m is the stand val) 
    #       lin  : time lenght 
    #       z    : temperature/Ri (2d array with float information)
    #
    #       zl: temp/Ri/variable (2d array according to the elevation yl)
    #
    mx = np.amax(h)                  # maximum elevation in the y array
    mn = np.amin(h) 

    le = int((ma.ceil(mx)-ma.floor(mn))/dj)                 # yl vector lenght     
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
                nnum[t][z] = interpol(num[t][i],num[t][i+1],h[t][i],h[t][i+1],nh[z])  
            if i < (limit - 1): 
                if nh[z] < h[t][i+1]:
                    i=i+1

    return nnum, nh
    
    
def wavelet(dx,per, power, show, ax):
#--- Contour plot wavelet power spectrum

    levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16]
    x = ax.contourf(dx,per,np.log2(power),np.log2(levels), \
              extend='both', cmap=plt.get_cmap('plasma')) 
    

    ax.set_yscale('log')
    ax.set_ylabel('scale (hours)')
    ax.set_ylim([0,40])
    
    
    if show == 0: 
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
        plt.colorbar(x, orientation='vertical', \
                     label='wavelet (log2(m²/Hz))',aspect=40, ax=ax)

def wavelet_depth(dx,per, power, time, ax):
#--- Contour plot wavelet power spectrum

    levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16]
    x = ax.contourf(dx,per,np.log2(power),np.log2(levels), \
              extend='both', cmap=plt.get_cmap('plasma')) 
    
    
    ax.set_yscale('log')
    ax.set_ylabel('scale (hours)')
    ax.set_ylim([0,40])
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    plt.colorbar(x, orientation='vertical', \
                     label='wavelet (log2(°C²/Hz))',aspect=40, ax=ax)


def wavelet_resona(dx,per, power, v1, v2, ax):
#--- Contour plot wavelet power spectrum

    levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16]
    x = ax.contourf(dx,per,np.log2(power),np.log2(levels), \
              extend='both', cmap=plt.get_cmap('plasma')) 
    
    # 95% significance contour, levels at -99 (fake) and 1 (95% signif)
    #ax.contour(time,per,sig95,[1],color='k',linewidth=1)
    
    v1 = v1/60/60 # period in seconds to hours
    v2 = v2/60/60
    
    ax.plot(dx, v1, linewidth=1, c='black', ls='-', label='V1H1')
    ax.plot(dx, v2, linewidth=1, c='black', ls='--', label='V2H1')

    ax.legend(loc='upper right', prop={'size': 9})
    
    ax.set_yscale('log')
    ax.set_ylabel('scale (hours)')
    ax.set_ylim([0,40])
    
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    plt.colorbar(x, orientation='vertical', \
                     label='wavelet (log2((m/s)²/Hz))',aspect=40, ax=ax)


def graph_globalwavelet(f,wave,color,ax):
    ax.plot(wave, f, linewidth=1, c=color, ls='--')
    ax.set_yscale('log')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_xlabel('global wavelet spectrum')

def isotherm(date,t,y,color,iso,ax):

    
    ax.plot(date, y, linewidth=1, c=color, ls='-')
    
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    ax.set_xlim([date[0],date[-1]])
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_ylabel('isotherm '+str(iso)+ ' °C (m)')
    

def isodepth(date,t,y,color,depth,ax):
        
    ax.plot(date, y, linewidth=1, c=color, ls='-')
    
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    ax.set_xlim([date[0],date[-1]])
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_ylabel('temp. variation at '+str(round(np.nanmean(depth),2))+' m (°C)')



def multi_isotherms(date, y, iso, time,ax):
 
    color  = ['deepskyblue','dodgerblue','royalblue','darkblue']

    for i in range(4):
        if iso[i]!=-999:

            #t_smooth, t_numeral, y_smooth = smooth_date(date,time[i],y[i],len(y[i]))  
            #ax.plot(t_smooth, y_smooth, linewidth=1, c=color[i], ls='-', label = 'iso '+str(iso[i])+'°C')
            ax.plot(date, y[i], linewidth=1, c=color[i], ls='-', label = 'iso '+str(iso[i])+'°C')
           
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
    ax.set_xlim([date[0],date[-1]])
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.legend(loc='upper right', prop={'size': 8})
    ax.set_ylabel('depth (m)')   


def temp_bandpass(date, y, iso, ax):
    
    color= ['red','maroon','blue','navy']
    
    for i in range(4):
        if iso[i]!=-999:
            try:
                ax.plot(date, y[i], linewidth=1, c=color[i], ls='-', label='iso '+str(iso[i])+'°C')
            except ValueError:
                pass

    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    

    ax.set_xlim([date[0],date[-1]])
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.legend( loc='upper right', prop={'size': 9})
    ax.set_ylabel('displacement (m)')  

def nauticaldefinition(dw):

    dw = dw + 180
    dw = np.where(dw > 360, dw-360, dw)
    
    return dw


def wind_rose(ws,wd):
    
    wd = wd + 180
    wd = np.where(wd > 360, wd-360, wd) # nautica definition to wind direction
    
    
    ax = WindroseAxes.from_ax()
    rounda = int(max(ws))+1
    bins_range = np.linspace(0,rounda,7)
    ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white', bins=bins_range)
    ax.set_legend(fontsize=15)

def depth_bandpass(date, y, depth, time,s, ax):
    
   
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


def vector_time (a,j,ni):
    
    b =  np.zeros((ni),float)
    
    for i in range(ni):
        b[i] = a[i][j]
    return b    
    
def thermal_variation(date,depth, depth_point, temp, time, zoom, ax):
    # this does not account the surface elevation 
    points = len(time)
    
    #t_smooth, t_numeral, y_smooth = smooth_date(date,time,temp,points)  
    #ax.plot(t_smooth, y_smooth, linewidth=1, c='black', ls='-')
    ax.plot(date, temp, linewidth=1, c='black', ls='-')
    
    col = ['red','maroon','blue','navy']
    for i in range(4):  
        aux = depth[i]
        if aux[0] is not None:
            
            d = vector_time (temp,depth_point[i],points)

            #t_smooth, t_numeral, y_smooth = smooth_date(date,time,d,points)  
            #ax.plot(t_smooth, y_smooth, linewidth=1, c=col[i], ls='-', label=str(round(np.nanmean(depth[i]),1))+' m') 
            ax.plot(date, d, linewidth=1, c=col[i], ls='-', label=str(round(np.nanmean(depth[i]),1))+' m') 
      
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
    ax.set_xlim([date[0],date[-1]])
    ax.set_ylabel('water temperature (°C)')  
    
    ax.legend(loc='upper right',prop={'size': 9})
    
    if zoom == 1:
        ax.set_ylim([4,8])
    

def thermocline(date,t,y,color,thermo,ax):
    
    ax.plot(date, y, linewidth=1, c=color, ls='-')
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    ax.set_xlim([date[0],date[-1]])
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_ylabel(str(thermo)+ ' m (°C)')


def psd_thermocline(p,psd,color,depth,ax):

    ax.plot(psd,p, linewidth=1, c=color, ls='-')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_xlabel('PSD '+ str(depth) +' m (°C²/Hz)') 
    

    
def psd_isoz(tau, freq, psd, zv1h1,largelen,n,fo, ax):

    color = ['red','maroon','blue','navy']
    minimoy= 10**-6
    maximoy= 10**-3 

    maxima = 0

    for i in range(4):
        if(tau[i]!=-999):        
            if(largelen==1): 
                l0 = int(len(psd[i])/150)
                ax.plot(smooth(psd[i], l0), freq[i], linewidth=1, c=color[i], ls='-',label="iso "+str(tau[i])+"°C")
            else:
                ax.plot(psd[i], freq[i], linewidth=1, c=color[i], ls='-', label="iso "+str(tau[i])+"°C")
                
            max_x = valmax(psd[i],freq[i])
            if(max_x > maxima):
                maxima = max_x
        
    
    maxima=maxima+100
    ax.set_ylim([minimoy,maximoy])

    ax.set_yscale('log')
    ax.set_xscale('log')    
    
    ax.set_xlabel('PSD isotherms (m²/Hz)') 
    ax.set_ylabel('f(Hz)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    model_plot(zv1h1,'navy','z',1,1,maxima,ax)

    ax.set_xlim(right=10**ma.ceil(ma.log10(maxima)))
    
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[fo,fo],linewidth=1, c='black', ls='--')
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[n,n],linewidth=1, c='black', ls='--')
    
    ax.legend(loc='lower left')


def psd_multilayer(f,depth,ff,fi,n2,ax):
    

    ax.set_title('multi-layer spectrum (mean depth)', loc='left', fontsize=10)
    

    #cf = ax.contourf(f, depth, ff[:,::-1], 100, locator=ticker.LogLocator())
 
    
    ax.set_title('multi-layer spectrum (mean depth)', loc='left', fontsize=10)
    
    ff = np.array(ff)
    cf = ax.contourf(f, depth, np.log10(ff[::-1,:]), 100)
    plt.colorbar(cf,orientation='vertical',label='PSD (log(°C²/Hz))',ax=ax)
    
    ax.plot([fi,fi],[min(depth),max(depth)],linewidth=1, color='black', ls='--')
 
    if max(f) > 10**-2:
        ax.plot([n2,n2],[min(depth),max(depth)],linewidth=1, color='black', ls='--')
        
    ax.set_xscale('log')
    ax.set_xlim(left=10**-6)
    ax.set_ylabel('water level (m)')
    ax.set_xlabel('frequency (Hz)')



def psd_depthz(sen, seu, matrix_depth, f , psd, zv1h1, largelen, n,fo,ax):

    color = ['red','maroon','blue','navy']
    
    minimoy= 10**-6
    maximoy= 10**-3 
    
    maxima = 0
       
    l10 = int(len(psd[0])/150)
    

    for i in range(4):
        
        if(sen[i]==1):
            
            depth = np.nanmean(matrix_depth[i])
            if(largelen==1): 
                ax.plot(smooth(psd[seu[i]], l10), f[seu[i]], linewidth=1, c=color[i], ls='-', label=str(round(depth,2))+" m")
            else:
                ax.plot(psd[seu[i]], f[seu[i]], linewidth=1, c=color[i], ls='-', label=str(round(depth,2))+" m")
        
            max_x = valmax(psd[seu[i]],f[seu[i]])
            if(max_x > maxima):
                maxima = max_x
    
    maxima=maxima+50
    ax.set_ylim([minimoy,maximoy])

    ax.set_yscale('log')
    ax.set_xscale('log')
    
    ax.set_xlabel('PSD thermal fluctuations (°C²/Hz)') 
    ax.set_ylabel('f(Hz)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    model_plot(zv1h1,'navy','z',1,1,maxima,ax)
 
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[fo,fo],linewidth=1, c='black', ls='--')
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[n,n],linewidth=1, c='black', ls='--')    
    
    ax.set_xlim(right=10**ma.ceil(ma.log10(maxima)))

    ax.legend(loc='upper right',prop={'size': 9})

def valmax(psd,f):
    
    l    = len(psd)
    maxi = 0
    for i in range(l):
        
        if(f[i]>6*10**-6 and psd[i]>maxi):
            maxi = psd[i]
    
    return maxi
        
        

def psd_isox(tau, freq, psd, xv1h1, xv2h1, xv3h1,largelen,n,fo, ax):

    
    color = ['red','maroon','blue','navy']
    minimoy= 10**-6
    maximoy= 10**-3 

    maxima = 0

    for i in range(4):
        if(tau[i]!=-999):        
            if(largelen==1): 
                l0 = int(len(psd[i])/150)
                ax.plot(smooth(psd[i], l0), freq[i], linewidth=1, c=color[i], ls='-',label="iso "+str(tau[i])+"°C")
            else:
                ax.plot(psd[i], freq[i], linewidth=1, c=color[i], ls='-', label="iso "+str(tau[i])+"°C")
                
            max_x = valmax(psd[i],freq[i])
            if(max_x > maxima):
                maxima = max_x
        

    maxima = maxima + 50
    ax.set_ylim([minimoy,maximoy])
    
    ax.set_yscale('log')
    ax.set_xscale('log')

    
    ax.set_xlabel('PSD isotherms (m²/Hz)') 
    ax.set_ylabel('f(Hz)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    model_plot(xv1h1,'green','x',1,1,maxima,ax)
    model_plot(xv2h1,'red','x',2,1,maxima,ax)
    model_plot(xv3h1,'blue','x',3,1,maxima,ax)

    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[fo,fo],linewidth=1, c='black', ls='--')
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[n,n],linewidth=1, c='black', ls='--')

    ax.set_xlim(right=10**ma.ceil(ma.log10(maxima)))
    ax.legend(loc='lower left')

def psd_isox_depth(sen, matrix_depth, f, psd, xv1h1, xv2h1, xv3h1,largelen,n,fo, ax):


    color = ['red','maroon','blue','navy']
    
    minimoy= 10**-6
    maximoy= 10**-3 
    
    maxima = 0
    
    l0 = int(len(psd[0])/150)
    
    for i in range(4):
        
        if(sen[i]==1):
            
            depth = np.nanmean(matrix_depth[i])
            
            if(largelen==1): 
                ax.plot(smooth(psd[i], l0), f[i], linewidth=1, c=color[i], ls='-', label=str(round(depth,2))+" m")
            else:
                ax.plot(psd[i], f[i], linewidth=1, c=color[i], ls='-', label=str(round(depth,2))+" m")
        
            max_x = valmax(psd[i],f[i])
            if(max_x > maxima):
                maxima = max_x


    maxima = maxima + 50
    ax.set_ylim([minimoy,maximoy])
 
    ax.set_yscale('log')
    ax.set_xscale('log')    
    
    ax.set_xlabel('PSD isotherms (m²/Hz)') 
    ax.set_ylabel('f(Hz)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    model_plot(xv1h1,'green','x',1,1,maxima,ax)
    model_plot(xv2h1,'red','x',2,1,maxima,ax)
    model_plot(xv3h1,'blue','x',3,1,maxima,ax)
    
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[fo,fo],linewidth=1, c='black', ls='--')
    ax.plot([0,10**ma.ceil(ma.log10(maxima))],[n,n],linewidth=1, c='black', ls='--')

    ax.set_xlim(right=10**ma.ceil(ma.log10(maxima)))
    ax.legend(loc='lower left')    


def psd_wind(period,psd_wind,ax):

    ax.plot(psd_wind, period, linewidth=1, c='navy', ls='-')
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    ax.tick_params(axis='x', labelcolor='navy' )
    ax.set_xlabel('PSD wind ((m/s)²/Hz)', color='navy') 
    ax.set_ylabel('periodicity (hours)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)

    
def psd_sola(period,psd_sola,ax):

    ax.plot(psd_sola, period, linewidth=1, c='red', ls='-')
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    ax.tick_params(axis='x', labelcolor='red')
    
    ax.set_xlabel('PSD solar radiation ((W/m²)²/Hz)', color='red') 
    ax.set_ylabel('periodicity (hours)')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)


def coherence (c1w,c1r,f1w,f1r,tau,iso, largelen,ax):
    
    f1w = f1w/60/60 # hour to second
    f1r = f1r/60/60 # hour to second
    
    l10 = int(len(c1w)/100)
    
    if largelen==1:
        ax.plot(smooth(c1w,l10), f1w, linewidth=1, c='navy', ls='-', \
                label='cohe wind')
        ax.plot(smooth(c1r,l10), f1r, linewidth=1, c='blue', ls='-', \
                label='cohe rad')
    else:
        ax.plot(c1w, f1w, linewidth=1, c='navy', ls='-', label='cohe wind')
        ax.plot(c1r, f1r, linewidth=1, c='blue', ls='-', label='cohe rad')
    
    ax.set_xlim([0,1])
    ax.set_xlabel('coherence iso'+ str(tau[iso]) + '/' + 'met. data')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    ax.legend(loc='lower right')


def coh_depthwind (cow,c1r,faws,f1r,depth,largelen,s,ax):
        
    color = ['red','maroon','blue','navy']
      
    for i in range(4):
        if s[i] == 1:
            
            mean = np.nanmean(depth[i])
            faws[i] = faws[i]*(1/60)*(1/60) 
            f1r = f1r/60/60 # hour to second
            
            if largelen==1:
                l0 = int(len(cow[i])/100)
                ax.plot(smooth(cow[i],l0), faws[i], linewidth=1, c=color[i], ls='-', label='cohe wind ~ '+str(round(mean,1))+' m')
        
                if s[0] == 1:
                    ax.plot(smooth(c1r,l0), f1r, linewidth=1, c='black', ls='-', label='cohe rad. ~ '+str(round(depth[4],1))+' m')
            else:
                ax.plot(cow[i], faws[i], linewidth=1, c=color[i], ls='-', label='cohe wind ~ '+str(round(mean,1))+' m')

                if s[0] == 1:
                    
                    ax.plot(c1r, f1r, linewidth=1, c='black', ls='-', label='cohe rad. ~ '+str(round(mean,1))+' m')
    
    ax.set_xlim([0,1])
    ax.set_xlabel('coherence elev./met. data')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    ax.legend(loc='lower right',prop={'size': 9})


def code_index(r):
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
    
    c = ['red','dodgerblue','lightseagreen','plum','dimgray','gold']
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
                        ax.plot(1/fre[k],coh[k],linewidth=1, ls='--', color=c[k], label=str(t[i])+'/'+str(t[j])+'°C')


    ax.set_ylabel('coherence', color='gray')
    ax.set_xlabel('period (hours)')
    

def coherence_depth(matrix_depth,coh,fre,ax):

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


    c = ['darkred','navy','darkgreen','m','black','darkgoldenrod']
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
                        ax.scatter(1/f[conf[k]],abs(p[conf[k]]),marker='^', color=c[k], label=str(t[i])+'/'+str(t[j])+'°C')


    ax.set_ylabel('Phase (degrees)', color='black')
    ax.set_yticks([ 0, np.pi/2, np.pi])
    ax.set_yticklabels([ r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])


    ax.grid(axis='y',color='black',ls=":",lw=0.25)
    ax.legend(loc='upper right',prop={'size': 6})
    ax.set_ylabel('phase (radian)')



def phase_depth(matrix_depth,conf, ph,fr,ax):

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


def density(dx,t,pu,pd,pm,ax):

    ax.plot(dx, pu, linewidth=1, c='darkred', ls='-')
    ax.plot(dx, pm, linewidth=1, c='black', ls='-')
    ax.plot(dx, pd, linewidth=1, c='darkblue', ls='-')
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    ax.invert_yaxis()
    
    ax.set_ylabel('density (kg/m³)')
    ax.legend(['superficial', 'average','bottom'], loc='upper right',prop={'size': 8})

def radiation(date,t,y,ax):
    
    ax.plot(date, y, linewidth=1, c='darkgoldenrod', ls='-')
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    ax.set_ylabel('solar rad.(W/m²)')
    
def wind_direction(date,time,y,ys,ylow,ydi,wedd_up,wedd_low,ax):

    ax.plot(date, ydi, c='blue', linewidth=3, ls='-', label ='homogeneous wind events')
    
    ax.plot(date, y, linewidth=1, c='lightgray', ls='-')
    ax.plot(date, ys, linewidth=1, c='green', ls='-',label=str(round(wedd_up,1))+'< W < '+str(round(wedd_low,1)))
    ax.plot(date, ylow, linewidth=1, c='red', ls='-', label=str(round(min([wedd_up,1]),1))+'< W < 20')
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.set_ylabel('wind direction ($^{\circ}$)',color='gray') 
    
    ax.legend(loc='upper right', prop={'size': 8})
    ax.set_xlim([date[0],date[-1]])
    ax.set_ylim([0,360])
     
    
def wind(date,time,dire,velo,ax):

    l  = len(date)
    dt = time[1]-time[0]

    ax.plot(date, velo, linewidth=1, c='navy', ls='-')

    ax.set_xlim([date[0],date[-1]])
    ax.set_ylim([0,int(2*max(velo))])
    
    ax.tick_params(axis='y', labelcolor='navy' )
    ax.set_ylabel('wind velocity (m/s)',color='navy')  
    
    t_smooth, t_numeral, y_smooth = smooth_date(date,time,dire,int(l*dt))
    
    
    ax2 = ax.twinx() 
    ax2.plot(t_smooth, y_smooth, linewidth=1, c='gray', ls='-')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax2.tick_params(axis='y', labelcolor='gray' )
    ax2.set_ylabel('wind direction ($^{\circ}$)',color='gray') 
    ax2.set_ylim([-200,380])
    ax2.set_yticks([0, 180, 360])
    

def windstress (date,time, y, ax):

    ax.plot(date, y, linewidth=1, c='navy', ls='-')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2.00))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    
    ax.set_ylim(bottom=0)
    ax.set_ylabel('wind stress (N/m²)') 


def wedderburn (date,t,y,ax):
    
    ax.plot(date, y, linewidth=1, c='navy', ls='-')
    #ax.plot(date, smooth(y,len(y)/2), linewidth=1, c='navy', ls='-')
    
    ax.set_yscale('log')
    ax.set_ylabel('Wedd number') 
    
    ax.set_xlim([date[0],date[-1]])
    
    ax.set_xticklabels(date, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    

def richardson(dx,time,ri,dw,up,ax):
    
    
    ax.plot(dx, dw, linewidth=1, c='red', ls='-')
    ax.plot(dx, up, linewidth=1, c='red', ls='-')
    ax.plot(dx, ri, linewidth=1, c='navy', ls='-')
    
   
    ax.fill_between(dx, dw, up, where=up >= dw, facecolor='red', alpha=0.5, interpolate=True)
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)

    l= len(dx)
    ax.set_xlim([dx[0],dx[l-1]])

    ax.set_yscale('log')
    ax.set_ylabel(r'$Ri_{h_e}$')


def ri_compb(dx,time,ri,rid,dw,up,ax):
        
    ax.plot(dx, dw, linewidth=1, c='red', ls='-')
    ax.plot(dx, up, linewidth=1, c='red', ls='-')
    
    ax.plot(dx, ri, linewidth=1, c='navy', ls='-', label= 'standard Ri')
    ax.plot(dx, rid, linewidth=1, c='dodgerblue', ls='-', label='Ri duration + direction')
  
    ax.fill_between(dx, dw, up, where=up >= dw, facecolor='red', alpha=0.5, interpolate=True)
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.legend(loc='upper right', prop={'size': 9})

    l= len(dx)
    ax.set_xlim([dx[0],dx[l-1]])
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_yscale('log')
    ax.set_ylabel(r'$Ri_{h_e}$')


def richardson2d(dx,time,h,ri,thermo,cota,z0,ax):

    # lev_exp = np.arange(np.floor(np.log10(z.min())-1),np.ceil(np.log10(z.max())+1))

# cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())
    x = ax.contourf(dx, h, ri.transpose(), 100, locator=ticker.LogLocator())
    ax.plot(dx, thermo, linewidth=1, c='red', ls='-')
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')

    ax.set_xticklabels(dx, rotation=10)
    
    ax.set_ylim([z0,ma.trunc(np.max(cota))+1])
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    plt.colorbar(x, orientation='horizontal', label='Ri (-)')
    ax.set_ylabel('depth (m)')

    
def tempstructure(dx, time,h,thermo,temp,cota,z0,ax):
    

    x = ax.contourf(dx,h,temp.transpose(), 300)
    
    ax.plot(dx,thermo,linewidth=1, c='w', ls='-')
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')
    
    ax.set_ylim([z0,ma.trunc(np.max(cota))+1])
    
    ax.set_xticklabels(dx, rotation=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
        
    ax.set_ylabel('depth (m)')  
    plt.colorbar(x, orientation='horizontal', label='Temperature (°C)')
 
def thorpe_scale(dx, time,h,thorpe,cota,z0,ax):
    

    x = ax.contourf(dx,h,thorpe.transpose(), 300)
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')
    
    ax.set_ylim([z0,ma.trunc(np.max(cota))+1])
    
    ax.set_xticklabels(dx, rotation=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
        
    ax.set_ylabel('depth (m)')  
    plt.colorbar(x, orientation='horizontal', label='Thorpe Scale (°C)')
    
def tempstructure_zoom(dx, time,h,thermo,temp,cota,ax):
    

    x = ax.contourf(dx,h,temp.transpose(), 300)

    ax.plot(dx,thermo,linewidth=1, c='w', ls='-')
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')
    
    ax.set_ylim([ma.trunc(np.min(thermo))-1,ma.trunc(np.max(cota))+1])
    
    ax.set_xticklabels(dx, rotation=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
    
    ax.set_ylabel('depth (m)')   
    plt.colorbar(x, orientation='horizontal', label='Temperature (°C)')
    
   
def degeneration(p1,p2,p3,xa,xb,xc,ya,yb,yc,pe,ph,dh,H,ls,ax):
              # matalimnion thickness
    
    rhoe = pe
    rhoh = ph
    L    = ls
    
    Hdh  =  7.00
    dhe  = H/Hdh
    
    meta = dhe
    
    
    # boundaries lines
    xl = np.arange(0.01, 0.501, 0.001) 
    y1 = equation_one (H, xl, rhoe, rhoh, L, meta)
    y2 = equation_two (H, xl)
    y3 = equation_three (H, xl, meta)
    
    ax.plot(xl, y1, linewidth=1, c='blue', ls='-' )  
    ax.plot(xl, y2, linewidth=1, c='black', ls='-' )
    ax.plot(xl, y3, linewidth=1, c='red', ls='-' )

    cexa, ceya = np.mean(xa), np.mean(ya)  
    cexb, ceyb = np.mean(xb), np.mean(yb)
    cexc, ceyc = np.mean(xc), np.mean(yc)
    
    e1 = mpatches.Ellipse((cexa, ceya), 2*(xa[0]-cexa), 2*(ya[0]-ceya),\
                          color='red', label = "P1~"+str(p1[0].day)+'/'+\
                          str(p1[0].month)+'-'+str(p1[1].day)+'/'+\
                          str(p1[1].month) ,alpha=0.2)
    ax.add_artist(e1)
    
    e2 = mpatches.Ellipse((cexb, ceyb), 2*(xb[0]-cexb), 2*(yb[0]-ceyb),\
                          color='navy', label = "P2~"+str(p2[0].day)+'/'+\
                          str(p2[0].month)+'-'+str(p2[1].day)+'/'+\
                          str(p2[1].month), alpha=0.2)
    ax.add_artist(e2)
    
    e3 = mpatches.Ellipse((cexc, ceyc), 2*(xc[0]-cexc), 2*(yc[0]-ceyc),\
                          color='green', label= "P3~"+str(p3[0].day)+'/'+\
                          str(p3[0].month)+'-'+str(p3[1].day)+'/'+\
                          str(p3[1].month), alpha=0.2)
    ax.add_artist(e3)

    ax.legend(handles=[e1,e2,e3],loc='upper right',prop={'size': 9})
    ax.set_ylabel('1/W')
    ax.set_xlabel(r'$h_e/H$') 
    
    ax.set_xlim([0.1,0.5])
    ax.set_ylim([0.0,2.0])

def movingaverage(interval, window_size):
    
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')
   
def classification_genera(W, heh_found,ahe,tau, ax):


    heh = np.linspace(0.889, 210, num=500, endpoint=True)

    lin = 0.4996*(heh**(-1))

    label=['Spigel and Imberger (1980)', ' Bueno et al. (2020) '+r'$h_e/H=0.5$',' Bueno et al. (2020) '+r'$h_e/H=0.3$', "Analyzed period   "+ r'$h_e/H= $'+str(round(heh_found,2))]
    
    ax.plot(heh,lin,linewidth=1, ls='-')

    wedderburn = [100,	99,	98,	97,	96,	95,	94,	93,	92,	91,	90,	89,	88,	87,	86,	85,	84,	83,	82,	81,	80,	79,	78,	77,	76,	75,	74,	73,	72,	71,	70,	69,	68,	67,	66,	65,	64,	63,	62,	61,	60,	59,	58,	57,	56,	55,	54,	53,	52,	51,	50,	49,	48,	47,	46,	45,	44,	43,	42,	41,	40,	39,	38,	37,	36,	35,	34,	33,	32,	31,	30,	29,	28,	27,	26,	25,	24,	23,	22,	21,	20,	19,	18,	17,	16,	15,	14,	13,	12,	11,	10,	9.9,	9.8,	9.7,	9.6,	9.5,	9.4,	9.3,	9.2,	9.1,	9,	8.9,	8.8,	8.7,	8.6,	8.5,	8.4,	8.3,	8.2,	8.1,	8,	7.9,	7.8,	7.7,	7.6,	7.5,	7.4,	7.3,	7.2,	7.1,	7,	6.9,	6.8,	6.7,	6.6,	6.5,	6.4,	6.3,	6.2,	6.1,	6,	5.9,	5.8,	5.7,	5.6,	5.5,	5.4,	5.3,	5.2,	5.1,	5,	4.9,	4.8,	4.7,	4.6,	4.5,	4.4,	4.3,	4.2,	4.1,	4,	3.9,	3.8,	3.7,	3.6,	3.5,	3.4,	3.3,	3.2,	3.1,	3,	2.9,	2.8,	2.7,	2.6,	2.5,	2.4,	2.3,	2.2,	2.1,	2,	1.9,	1.8,	1.7,	1.6,	1.5,	1.4,	1.3,	1.2,	1.1,	1,	0.99,	0.98,	0.97,	0.96,	0.95,	0.94,	0.93,	0.92,	0.91,	0.9,	0.89,	0.88,	0.87,	0.86,	0.85,	0.84,	0.83,	0.82,	0.81,	0.8,	0.79,	0.78,	0.77,	0.76,	0.75,	0.74,	0.73,	0.72,	0.71,	0.7,	0.69,	0.68,	0.67,	0.66,	0.65,	0.64,	0.63,	0.62,	0.61,	0.6,	0.59,	0.58,	0.57,	0.56,	0.55,	0.54,	0.53,	0.52,	0.51,	0.5,	0.49,	0.48,	0.47,	0.46,	0.45,	0.44,	0.43,	0.42,	0.41,	0.4,	0.39,	0.38,	0.37,	0.36,	0.35,	0.34,	0.33,	0.32,	0.31,	0.3,	0.29,	0.28,	0.27,	0.26,	0.25,	0.24,	0.23,	0.22,	0.21,	0.2,	0.19,	0.18,	0.17,	0.16,	0.15,	0.14,	0.13,	0.12,	0.11,	0.1]
    curve1 = [0.00000000000132836,	0.00000000000172279,	0.00000000000223434,	0.00000000000289779,	0.00000000000375823,	0.00000000000487416,	0.00000000000632144,	0.00000000000819847,	0.0000000000106328,	0.0000000000137901,	0.0000000000178847,	0.0000000000231952,	0.0000000000300826,	0.000000000039015,	0.0000000000505998,	0.0000000000656244,	0.0000000000851102,	0.000000000110382,	0.000000000143158,	0.000000000185666,	0.000000000240795,	0.000000000312295,	0.000000000405024,	0.000000000525288,	0.000000000681262,	0.00000000088355,	0.0000000011459,	0.00000000148615,	0.00000000192744,	0.00000000249975,	0.00000000324201,	0.00000000420465,	0.00000000545314,	0.00000000707235,	0.00000000917234,	0.0000000118959,	0.0000000154281,	0.0000000200092,	0.0000000259505,	0.000000033656,	0.0000000436495,	0.0000000566104,	0.0000000734197,	0.0000000952202,	0.000000123494,	0.000000160163,	0.00000020772,	0.000000269398,	0.000000349391,	0.000000453135,	0.000000587684,	0.000000762185,	0.0000009885,	0.00000128201,	0.00000166268,	0.00000215638,	0.00000279666,	0.00000362706,	0.00000470403,	0.00000610076,	0.0000079122,	0.0000102615,	0.0000133083,	0.0000172596,	0.0000223841,	0.0000290299,	0.0000376485,	0.0000488253,	0.0000633195,	0.000082115,	0.000106487,	0.00013809,	0.000179065,	0.000232187,	0.000301051,	0.000390307,	0.000505976,	0.000655836,	0.000849937,	0.00110124,	0.001426438,	0.001846981,	0.002390363,	0.003091694,	0.003995603,	0.005158479,	0.006651005,	0.008560863,	0.010995341,	0.014083342,	0.017975989,	0.018415955,	0.01886585,	0.019325857,	0.019796158,	0.020276941,	0.020768391,	0.021270695,	0.021784041,	0.022308619,	0.022844618,	0.023392227,	0.023951638,	0.02452304,	0.025106623,	0.025702579,	0.026311096,	0.026932364,	0.027566571,	0.028213905,	0.028874551,	0.029548695,	0.030236519,	0.030938205,	0.03165393,	0.032383873,	0.033128206,	0.033887101,	0.034660724,	0.035449239,	0.036252807,	0.037071583,	0.037905718,	0.03875536,	0.039620648,	0.040501718,	0.041398702,	0.042311721,	0.043240893,	0.044186327,	0.045148128,	0.04612639,	0.047121199,	0.048132635,	0.049160768,	0.050205658,	0.051267356,	0.052345904,	0.053441333,	0.054553663,	0.055682904,	0.056829055,	0.057992101,	0.059172017,	0.060368765,	0.061582296,	0.062812545,	0.064059437,	0.065322882,	0.066602775,	0.067898998,	0.069211421,	0.070539897,	0.071884264,	0.073244348,	0.074619958,	0.076010887,	0.077416917,	0.07883781,	0.080273317,	0.08172317,	0.08318709,	0.084664779,	0.086155927,	0.087660208,	0.08917728,	0.090706788,	0.092248363,	0.093801621,	0.095366163,	0.096941581,	0.098527448,	0.100123329,	0.101728773,	0.103343321,	0.104966498,	0.106597822,	0.108236799,	0.109882926,	0.111535688,	0.113194564,	0.113360768,	0.113527027,	0.11369334,	0.113859708,	0.11402613,	0.114192605,	0.114359132,	0.114525712,	0.114692343,	0.114859025,	0.115025758,	0.115192541,	0.115359373,	0.115526253,	0.115693183,	0.11586016,	0.116027184,	0.116194255,	0.116361372,	0.116528534,	0.116695742,	0.116862994,	0.11703029,	0.117197629,	0.117365012,	0.117532436,	0.117699903,	0.11786741,	0.118034959,	0.118202547,	0.118370175,	0.118537842,	0.118705547,	0.11887329,	0.119041071,	0.119208888,	0.119376742,	0.119544631,	0.119712555,	0.119880514,	0.120048507,	0.120216533,	0.120384592,	0.120552684,	0.120720807,	0.120888962,	0.121057147,	0.121225362,	0.121393607,	0.12156188,	0.121730182,	0.121898512,	0.12206687,	0.122235253,	0.122403664,	0.122572099,	0.12274056,	0.122909045,	0.123077554,	0.123246086,	0.123414642,	0.123583219,	0.123751818,	0.123920438,	0.124089078,	0.124257739,	0.124426419,	0.124595117,	0.124763834,	0.124932569,	0.125101321,	0.125270089,	0.125438873,	0.125607673,	0.125776487,	0.125945316,	0.126114158,	0.126283013,	0.126451881,	0.126620761,	0.126789653,	0.126958555,	0.127127468,	0.12729639,	0.127465321,	0.127634261,	0.127803209,	0.127972164,	0.128141127,	0.128310095]
    curve2 = [4.84408E-21,	7.67339E-21,	1.21552E-20,	1.92548E-20,	3.0501E-20,	4.83158E-20,	7.65358E-20,	1.21238E-19,	1.9205E-19,	3.04222E-19,	4.8191E-19,	7.63382E-19,	1.20925E-18,	1.91555E-18,	3.03437E-18,	4.80666E-18,	7.61411E-18,	1.20613E-17,	1.9106E-17,	3.02653E-17,	4.79425E-17,	7.59445E-17,	1.20302E-16,	1.90567E-16,	3.01872E-16,	4.78187E-16,	7.57484E-16,	1.19991E-15,	1.90075E-15,	3.01093E-15,	4.76953E-15,	7.55529E-15,	1.19681E-14,	1.89584E-14,	3.00315E-14,	4.75722E-14,	7.53578E-14,	0.000000000000119372,	0.000000000000189095,	0.00000000000029954,	0.000000000000474493,	0.000000000000751633,	0.00000000000119064,	0.00000000000188606,	0.00000000000298767,	0.00000000000473268,	0.00000000000749692,	0.0000000000118757,	0.000000000018812,	0.0000000000297995,	0.0000000000472046,	0.0000000000747757,	0.00000000011845,	0.000000000187634,	0.000000000297226,	0.000000000470828,	0.000000000745826,	0.00000000118144,	0.00000000187149,	0.00000000296458,	0.00000000469612,	0.000000007439,	0.0000000117839,	0.0000000186666,	0.0000000295693,	0.00000004684,	0.000000074198,	0.000000117535,	0.000000186184,	0.00000029493,	0.00000046719,	0.000000740063,	0.00000117231,	0.00000185703,	0.00000294166,	0.0000046598,	0.00000738142,	0.0000116926,	0.0000185217,	0.000029339,	0.0000464734,	0.000073613,	0.000116598,	0.000184672,	0.000292465,	0.000463115,	0.000733177,	0.001160326,	0.001835338,	0.002900552,	0.004577829,	0.004791083,	0.005014162,	0.005247508,	0.005491582,	0.005746866,	0.006013859,	0.006293085,	0.006585088,	0.006890433,	0.007209713,	0.00754354,	0.007892555,	0.008257423,	0.008638836,	0.009037514,	0.009454206,	0.009889689,	0.010344771,	0.01082029,	0.011317117,	0.011836156,	0.012378342,	0.012944648,	0.013536078,	0.014153675,	0.014798517,	0.015471718,	0.01617443,	0.016907846,	0.017673193,	0.018471741,	0.019304795,	0.020173702,	0.021079848,	0.022024657,	0.023009594,	0.024036159,	0.025105892,	0.026220371,	0.027381208,	0.028590051,	0.029848581,	0.031158509,	0.032521578,	0.033939557,	0.03541424,	0.036947442,	0.038540996,	0.04019675,	0.041916562,	0.043702296,	0.045555816,	0.047478983,	0.049473646,	0.051541638,	0.053684767,	0.055904812,	0.058203513,	0.060582561,	0.063043595,	0.065588185,	0.068217832,	0.070933948,	0.073737853,	0.076630762,	0.079613774,	0.08268786,	0.085853855,	0.08911244,	0.09246414,	0.095909304,	0.099448101,	0.103080502,	0.106806277,	0.110624978,	0.114535936,	0.118538248,	0.12263077,	0.126812112,	0.131080631,	0.135434427,	0.139871338,	0.14438894,	0.148984548,	0.153655214,	0.158397731,	0.163208642,	0.168084238,	0.173020576,	0.178013479,	0.178515713,	0.179018465,	0.179521729,	0.180025502,	0.180529779,	0.181034555,	0.181539825,	0.182045585,	0.182551831,	0.183058557,	0.183565759,	0.184073433,	0.184581573,	0.185090175,	0.185599233,	0.186108744,	0.186618703,	0.187129104,	0.187639943,	0.188151215,	0.188662915,	0.189175038,	0.18968758,	0.190200535,	0.190713898,	0.191227665,	0.191741831,	0.19225639,	0.192771337,	0.193286669,	0.193802378,	0.194318462,	0.194834913,	0.195351728,	0.195868901,	0.196386428,	0.196904302,	0.19742252,	0.197941075,	0.198459962,	0.198979178,	0.199498715,	0.200018569,	0.200538735,	0.201059208,	0.201579982,	0.202101052,	0.202622413,	0.20314406,	0.203665987,	0.204188189,	0.20471066,	0.205233396,	0.205756391,	0.206279639,	0.206803136,	0.207326875,	0.207850853,	0.208375062,	0.208899498,	0.209424155,	0.209949029,	0.210474113,	0.210999402,	0.211524891,	0.212050573,	0.212576445,	0.2131025,	0.213628733,	0.214155138,	0.21468171,	0.215208443,	0.215735332,	0.216262372,	0.216789556,	0.217316879,	0.217844337,	0.218371922,	0.218899631,	0.219427456,	0.219955393,	0.220483436,	0.22101158,	0.221539819,	0.222068147,	0.222596559,	0.223125049,	0.223653611,	0.224182241,	0.224710933]

    ax.plot(wedderburn,curve1,linewidth=1, ls='-')
    ax.plot(wedderburn,curve2,linewidth=1, ls='-')

    ax.scatter(W,ahe,s=20,marker='o',color='red')


    ax.legend(labels=label,loc='upper right',prop={'size': 8})
    
    ax.grid(axis='y',color='black',ls=":",lw=0.25)
    ax.grid(axis='x',color='black',ls=":",lw=0.25)

    ax.set_ylabel(r'$\zeta_o/h_e$'+' (-)', color='black')
    ax.set_xlabel(r'$Ri h_e/L (-)$', color='black')

    ax.set_xscale('log')
    ax.set_xlim([0.01,100])
    ax.set_ylim([0,0.5])
    
def bueno_parameterization(W, norma, tau, ax):
    

    label=['bueno at al. (2020)', 'Analyzed period - '+str(tau[0])+ '°C']

    norm = [4.24835E-18,	6.3378E-18,	9.45489E-18,	1.4105E-17,	2.10422E-17,	3.13913E-17,	4.68304E-17,	6.98627E-17,	1.04223E-16,	1.55482E-16,	2.31952E-16,	3.46032E-16,	5.16219E-16,	7.70109E-16,	1.14887E-15,	1.71391E-15,	2.55685E-15,	3.81437E-15,	5.69038E-15,	8.48904E-15,	1.26642E-14,	1.88927E-14,	2.81846E-14,	4.20465E-14,	0.000000000000062726,	9.35762E-14,	0.000000000000139599,	0.000000000000208258,	0.000000000000310684,	0.000000000000463486,	0.00000000000069144,	0.00000000000103151,	0.00000000000153883,	0.00000000000229566,	0.00000000000342472,	0.00000000000510909,	0.00000000000762187,	0.0000000000113705,	0.0000000000169628,	0.0000000000253055,	0.0000000000377513,	0.0000000000563184,	0.0000000000840172,	0.000000000125339,	0.000000000186984,	0.000000000278947,	0.00000000041614,	0.000000000620808,	0.000000000926136,	0.00000000138163,	0.00000000206115,	0.00000000307488,	0.00000000458718,	0.00000000684327,	0.000000010209,	0.00000001523,	0.0000000227205,	0.0000000338949,	0.0000000505653,	0.0000000754346,	0.000000112535,	0.000000167883,	0.000000250452,	0.00000037363,	0.00000055739,	0.000000831528,	0.00000124049,	0.0000018506,	0.00000276076,	0.00000411857,	0.00000614417,	0.000009166,	0.000013674,	0.0000203991,	0.0000304316,	0.0000453979,	0.0000677241,	0.000101029,	0.00015071,	0.000224817,	0.00033535,	0.000500201,	0.000746029,	0.001112536,	0.001658801,	0.002472623,	0.00368424,	0.005486299,	0.008162571,	0.012128435,	0.01798621,	0.01870651,	0.019455085,	0.020232997,	0.021041347,	0.021881271,	0.022753943,	0.023660578,	0.024602428,	0.025580788,	0.026596994,	0.027652422,	0.028748496,	0.02988668,	0.031068484,	0.032295465,	0.033569223,	0.034891409,	0.036263716,	0.037687891,	0.039165723,	0.040699054,	0.042289772,	0.043939815,	0.045651171,	0.047425873,	0.049266006,	0.051173701,	0.053151136,	0.055200538,	0.057324176,	0.059524366,	0.061803466,	0.064163876,	0.066608036,	0.06913842,	0.071757542,	0.074467945,	0.077272202,	0.080172912,	0.083172696,	0.086274194,	0.089480059,	0.092792953,	0.096215542,	0.099750489,	0.103400451,	0.10716807,	0.111055967,	0.115066732,	0.119202922,	0.123467048,	0.127861566,	0.132388874,	0.137051293,	0.141851065,	0.14679034,	0.151871164,	0.157095469,	0.162465063,	0.167981615,	0.173646647,	0.179461519,	0.185427419,	0.191545349,	0.197816111,	0.204240302,	0.210818293,	0.217550224,	0.224435986,	0.231475217,	0.238667285,	0.246011284,	0.253506017,	0.261149994,	0.268941421,	0.276878195,	0.284957894,	0.293177779,	0.301534784,	0.310025519,	0.318646266,	0.327392983,	0.336261303,	0.345246539,	0.354343694,	0.36354746,	0.372852234,	0.382252125,	0.391740969,	0.40131234,	0.402273761,	0.403235934,	0.404198853,	0.405162509,	0.406126897,	0.40709201,	0.40805784,	0.40902438,	0.409991625,	0.410959566,	0.411928197,	0.41289751,	0.4138675,	0.414838158,	0.415809477,	0.416781451,	0.417754072,	0.418727333,	0.419701228,	0.420675748,	0.421650887,	0.422626637,	0.423602991,	0.424579942,	0.425557483,	0.426535606,	0.427514305,	0.428493571,	0.429473397,	0.430453776,	0.431434701,	0.432416164,	0.433398158,	0.434380675,	0.435363708,	0.43634725,	0.437331292,	0.438315828,	0.43930085,	0.440286351,	0.441272322,	0.442258757,	0.443245647,	0.444232986,	0.445220765,	0.446208977,	0.447197615,	0.44818667,	0.449176135,	0.450166003,	0.451156265,	0.452146914,	0.453137943,	0.454129343,	0.455121108,	0.456113228,	0.457105697,	0.458098506,	0.459091648,	0.460085115,	0.4610789,	0.462072994,	0.46306739,	0.464062079,	0.465057055,	0.466052309,	0.467047833,	0.468043619,	0.46903966,	0.470035948,	0.471032475,	0.472029233,	0.473026213,	0.474023409,	0.475020813,	0.476018415,	0.477016209,	0.478014186,	0.479012339,	0.48001066,	0.48100914,	0.482007772,	0.483006548,	0.484005459,	0.485004498,	0.486003658,	0.487002929,	0.488002303,	0.489001774,	0.490001333]
    wedderburn = [100,	99,	98,	97,	96,	95,	94,	93,	92,	91,	90,	89,	88,	87,	86,	85,	84,	83,	82,	81,	80,	79,	78,	77,	76,	75,	74,	73,	72,	71,	70,	69,	68,	67,	66,	65,	64,	63,	62,	61,	60,	59,	58,	57,	56,	55,	54,	53,	52,	51,	50,	49,	48,	47,	46,	45,	44,	43,	42,	41,	40,	39,	38,	37,	36,	35,	34,	33,	32,	31,	30,	29,	28,	27,	26,	25,	24,	23,	22,	21,	20,	19,	18,	17,	16,	15,	14,	13,	12,	11,	10,	9.9,	9.8,	9.7,	9.6,	9.5,	9.4,	9.3,	9.2,	9.1,	9,	8.9,	8.8,	8.7,	8.6,	8.5,	8.4,	8.3,	8.2,	8.1,	8,	7.9,	7.8,	7.7,	7.6,	7.5,	7.4,	7.3,	7.2,	7.1,	7,	6.9,	6.8,	6.7,	6.6,	6.5,	6.4,	6.3,	6.2,	6.1,	6,	5.9,	5.8,	5.7,	5.6,	5.5,	5.4,	5.3,	5.2,	5.1,	5,	4.9,	4.8,	4.7,	4.6,	4.5,	4.4,	4.3,	4.2,	4.1,	4,	3.9,	3.8,	3.7,	3.6,	3.5,	3.4,	3.3,	3.2,	3.1,	3,	2.9,	2.8,	2.7,	2.6,	2.5,	2.4,	2.3,	2.2,	2.1,	2,	1.9,	1.8,	1.7,	1.6,	1.5,	1.4,	1.3,	1.2,	1.1,	1,	0.99,	0.98,	0.97,	0.96,	0.95,	0.94,	0.93,	0.92,	0.91,	0.9,	0.89,	0.88,	0.87,	0.86,	0.85,	0.84,	0.83,	0.82,	0.81,	0.8,	0.79,	0.78,	0.77,	0.76,	0.75,	0.74,	0.73,	0.72,	0.71,	0.7,	0.69,	0.68,	0.67,	0.66,	0.65,	0.64,	0.63,	0.62,	0.61,	0.6,	0.59,	0.58,	0.57,	0.56,	0.55,	0.54,	0.53,	0.52,	0.51,	0.5,	0.49,	0.48,	0.47,	0.46,	0.45,	0.44,	0.43,	0.42,	0.41,	0.4,	0.39,	0.38,	0.37,	0.36,	0.35,	0.34,	0.33,	0.32,	0.31,	0.3,	0.29,	0.28,	0.27,	0.26,	0.25,	0.24,	0.23,	0.22,	0.21,	0.2,	0.19,	0.18,	0.17,	0.16,	0.15,	0.14,	0.13,	0.12,	0.11,	0.1]
 
    
    ax.plot(wedderburn,norm,linewidth=1,color='black', ls=':')
    ax.scatter(W,norma,s=20,marker='o',color='red')

    ax.legend(labels=label,loc='lower left',prop={'size': 8})

    ax.grid(axis='y',color='black',ls=":",lw=0.25)
    ax.grid(axis='x',color='black',ls=":",lw=0.25)

    ax.set_ylabel(r'$\zeta_o/(h_e f)$'+' (-)', color='black')
    ax.set_xlabel(r'$Ri h_e/L (-)$', color='black')

    ax.set_xscale('log')
    ax.set_xlim([0.01,100])
    ax.set_ylim([0,0.8])


def degeneration_evolution(p1,p2,p3,x,y,pe,ph,dh,H,ls,zoom,ax):
              # matalimnion thickness
    
    rhoe = pe
    rhoh = ph
    L    = ls
    
    Hdh  =  7.00
    dhe  = H/Hdh
    
    meta = dhe
        
    # boundaries lines
    xl = np.arange(0.01, 0.501, 0.001) 
    y1 = equation_one (H, xl, rhoe, rhoh, L, meta)
    y2 = equation_two (H, xl)
    y3 = equation_three (H, xl, meta)
    
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
    

    
    minix = min(np.nanmin(xa_moved),np.nanmin(xb_moved),np.nanmin(xc_moved))-0.05
    if(minix < 0.1):
        minix = 0.1
    maxix = max(np.nanmax(xa_moved),np.nanmax(xb_moved),np.nanmax(xc_moved))+0.05
    if(maxix > 0.5):
        maxix = 0.5
    miniy = min(np.nanmin(ya_moved),np.nanmin(yb_moved),np.nanmin(yc_moved))-0.05
    if(miniy <0):
        miniy=0
    maxiy = max(np.nanmax(ya_moved),np.nanmax(yb_moved),np.nanmax(yc_moved))+0.05
    if(maxiy >2):
        maxiy=2    
    
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
    
def generation(cla,ax):

    N       = len(cla)
    x       = range(N)
    width   = 1/1.5
    
    ax.bar(x,cla,width,color="navy")
    ax.set_ylim([0,100])
    ax.set_xticklabels(("-","Mixed","KH","IS/d","IS/s"))


def depth_sensitivity(xh,per,depth,labelx,ax):
    
    font = {'family': 'sans-serif',
    'color':  'black',
    'weight': 'normal',
    'size': 10,
    }
    
    ax.plot(xh, per, linewidth=1, c='navy', ls='-')
    
    ax.text((np.nanmax(xh)+np.nanmin(xh))/2, (np.nanmax(per)+np.nanmin(per))/2, r'$h^* \equiv$'+\
            str(depth), fontdict=font)
    
    
    ax.set_ylabel('period (h)')
    
    if(labelx == 1):
        ax.set_xlabel(r'$h^*(m)$') 
        
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
def densi_sensitivity(xp,per,rho,labelx,ax):
    
    font = {'family': 'sans-serif',
    'color':  'black',
    'weight': 'normal',
    'size': 10,
    }
    
    ax.plot(xp, per, linewidth=1, c='navy', ls='-')
    
    ax.text((np.nanmax(xp)+np.nanmin(xp))/2, (np.nanmax(per)+np.nanmin(per))/2, r'$\rho^* \equiv$'+\
            str(rho), fontdict=font)
    
    
    if(labelx==1):
        ax.set_xlabel(r'$\rho^* (kg/m³)$') 

        
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    
def parabar_sensitivity(x,per,typ,ax):
    
    ax.plot(x, per, linewidth=1, c='navy', ls='-')
    
    if typ == 'dep':
        ax.set_ylabel('period (h)')
        ax.set_xlabel(r'$\bar{H}^{-1/2}$')
    elif typ == 'rho':
        ax.set_xlabel(r'$ (\rho_h/\Delta\rho)^{1/2}$')

