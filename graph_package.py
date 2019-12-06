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
                     label='wavelet (log2(°C/Hz))',aspect=40, ax=ax)

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
                     label='wavelet (log2(°C/Hz))',aspect=40, ax=ax)


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
                     label='wavelet (log2((m/s)/Hz))',aspect=40, ax=ax)


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
    
    color=['deepskyblue','dodgerblue','royalblue','darkblue']
    
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
    ax.set_ylabel('depth (m)')  



def wind_rose(ws,wd):
    
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
    

    
def psd_isoz(tau, freq, psd, zv1h1,largelen, ax):

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
    
    
    ax.legend(loc='lower left')

def psd_depthz(sen, matrix_depth, f , psd, zv1h1, largelen, ax):

    color = ['red','maroon','blue','navy']
    
    minimoy= 10**-6
    maximoy= 10**-3 
    
    maxima = 0
    
    
    
    l10 = int(len(psd[0])/150)
    
    for i in range(4):
        
        if(sen[i]==1):
            
            depth = np.nanmean(matrix_depth[i])
            if(largelen==1): 
                ax.plot(smooth(psd[i], l10), f[i], linewidth=1, c=color[i], ls='-', label=str(round(depth,2))+" m")
            else:
                ax.plot(psd[i], f[i], linewidth=1, c=color[i], ls='-', label=str(round(depth,2))+" m")
        
            max_x = valmax(psd[i],f[i])
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
    
    ax.set_xlim(right=10**ma.ceil(ma.log10(maxima)))

    ax.legend(loc='upper right',prop={'size': 9})

def valmax(psd,f):
    
    l    = len(psd)
    maxi = 0
    for i in range(l):
        
        if(f[i]>6*10**-6 and psd[i]>maxi):
            maxi = psd[i]
    
    return maxi
        
        

def psd_isox(tau, freq, psd, xv1h1, xv2h1, xv3h1,largelen, ax):

    
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

    ax.set_xlim(right=10**ma.ceil(ma.log10(maxima)))
    ax.legend(loc='lower left')

def psd_isox_depth(sen, matrix_depth, f, psd, xv1h1, xv2h1, xv3h1,largelen, ax):


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
                        ax.plot(1/fre[k],coh[k],linewidth=1, color=c[k], label=str(t[i])+'/'+str(t[j])+'°C')


    ax.set_xlim([0,25])
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
    
def wind_direction(date,time,y,ys,ylow,ydi,ax):

    ax.scatter(date, ydi, color='blue', marker='x', label ='homo. dir.')
 
    ax.plot(date, y, linewidth=1, c='lightgray', ls='-')
    ax.plot(date, ys, linewidth=1, c='green', ls='-',label='W < 100')
    ax.plot(date, ylow, linewidth=1, c='red', ls='-', label='W < 20')
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
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


def richardson2d(dx,time,h,ri,thermo,cota,ax):

    # lev_exp = np.arange(np.floor(np.log10(z.min())-1),np.ceil(np.log10(z.max())+1))

# cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())
    x = ax.contourf(dx, h, ri.transpose(), 100, locator=ticker.LogLocator())
    ax.plot(dx, thermo, linewidth=1, c='red', ls='-')
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')

    ax.set_xticklabels(dx, rotation=10)
    
    ax.set_ylim([799,ma.trunc(np.max(cota))+1])
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    plt.colorbar(x, orientation='horizontal', label='Ri (-)')
    ax.set_ylabel('depth (m)')

    
def tempstructure(dx, time,h,thermo,temp,cota,ax):
    

    x = ax.contourf(dx,h,temp.transpose(), 300)
    
    ax.plot(dx,thermo,linewidth=1, c='w', ls='-')
    ax.plot(dx, cota, linewidth=1, c='blue', ls='-')
    
    ax.set_ylim([799,ma.trunc(np.max(cota))+1])
    
    ax.set_xticklabels(dx, rotation=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    
        
    ax.set_ylabel('depth (m)')  
    plt.colorbar(x, orientation='horizontal', label='Temperature (°C)')
    
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
    
    ax.plot([0.1,0.5],[1/12,1/12],   linewidth=1,   c='deepskyblue', ls='--' )
    ax.plot([0.1,0.5],[1/3,1/3],     linewidth=1,   c='deepskyblue', ls='--' )
    ax.plot([0.1,0.5],[1/1,1/1],     linewidth=1,   c='deepskyblue', ls='--' )
    ax.plot([0.1,0.5],[1/0.8,1/0.8], linewidth=1,   c='deepskyblue', ls='--' )
    
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
            ax.set_xlim([minix,maxix])
            ax.set_ylim([miniy,maxiy])        
  
    
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