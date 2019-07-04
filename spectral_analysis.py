import numpy as np
import graph_package as graph
import matplotlib.pyplot as plt
import internal_module as mod

def graph_wavelet(per,power,time,sig95,ax):
#--- Contour plot wavelet power spectrum

    levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16]
    ax.contourf(time,per,np.log2(power),np.log2(levels), \
              extend='both', cmap=plt.get_cmap('jet')) 
    
    # 95% significance contour, levels at -99 (fake) and 1 (95% signif)
    ax.contour(time,per,sig95,[1],color='k',linewidth=1)
    ax.set_yscale('log')
    ax.set_ylabel('scale (hours)')
    ax.set_xlabel('time (hours)')

def graph_wel_period(f,ff,color,ax):
    ax.plot(f, ff, linewidth=1, c=color, ls='--')
    
    ax.set_xlabel('f [Hz]')
    ax.set_xscale('log')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_ylabel('PSD')
 

def graph_globalwavelet(f,wave,color,ax):
    ax.plot(wave, f, linewidth=1, c=color, ls='--')
    #ax.plot(signif,f,'k--')
    ax.set_yscale('log')
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    ax.set_xlabel('Global Wavelet Spectrum')

  
dt = 0.01 # sample spacing (hours)

# sinal estacionário
ts = np.loadtxt('naoestacionario.txt', unpack=True)

# sinal não estacionário

# sinal estacionário 2d
#ts2= np.loadtxt('estacionario2d.txt', unpack=True)

f1, y1 = mod.fft(ts, dt)
freq, per, ff = mod.welch_method(ts, dt)
time_win, per_win, power_win, global_ws_win, sigwin = mod.wave_spectral(ts, dt, 'MORLET')



plt.figure(0)

ax1 = plt.subplot2grid((2,2),(0,0))
ax2 = plt.subplot2grid((2,2),(0,1))
ax3 = plt.subplot2grid((2,2),(1,0))
ax4 = plt.subplot2grid((2,2),(1,1))



graph.graph_psd_isotherm(f1,y1,'blue',10,ax1)
graph.graph_psd_isotherm(freq,ff,'blue',10,ax2)
graph.graph_wavelet(per_win,power_win,time_win,sigwin,ax3)
graph.graph_globalwavelet(per_win,power_win,'green',ax4)
