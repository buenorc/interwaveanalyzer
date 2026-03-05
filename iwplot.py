# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer - Plots module

Interwave Analyzer - Version 2 (2026) 
Plots module version: 2.260305

-------------------------------------------------------------------------------

de Carvalho Bueno, R; Bleninger, T. B.; Lorke, A. 
Internal wave analyzer for thermally stratified lakes
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

import os
import math as ma
import numpy as np 
import matplotlib.dates as mdates
import matplotlib.ticker as ticker

import matplotlib.patches as mpatches


import iwmod as mod

from matplotlib import cm
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle

import matplotlib.colors as colors



class MidpointNormalize(colors.Normalize):
    """
    Normalize velocity field for color mapping.
    
    Defines the relative position of zero velocity within
    the colormap scaling.
    """

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def _draw_colorbar_with_background(ax, contour, label, log=False, cbar_n_ticks=4):
    """
    Add a colorbar with a white background box between
    the main plot and the colorbar.
    """
    fig = ax.get_figure()
    cax = inset_axes(ax, width="28%", height="6%", loc="upper right", borderpad=0.6)

    # Draw background *in axes coordinates*, so it's positioned correctly
    # regardless of layout, aspect, or tight_layout.
    # We use ax.transAxes to keep the placement relative to the axes.
    rect = Rectangle(
        (0.68, 0.77),  # upper right region in axes coords (approx position of colorbar)
        width=0.40,
        height=0.24,
        transform=ax.transAxes,
        facecolor='white',
        edgecolor='none',
        zorder=2.5  # between plot (zorder~1) and colorbar (zorder~3)
    )
    ax.add_patch(rect)

    # Create colorbar
    cbar = fig.colorbar(contour, cax=cax, orientation='horizontal')

    if log:
        # use log ticks (10^x)
        cbar.locator = ticker.LogLocator(base=10.0, numticks=cbar_n_ticks)
        cbar.formatter = ticker.LogFormatterMathtext(base=10.0)
    else:
        # normal linear colorbar
        cbar.locator = ticker.MaxNLocator(nbins=cbar_n_ticks)
        cbar.formatter = ticker.FormatStrFormatter('%.2f')

    cbar.update_ticks()

    # make text visible
    cbar.ax.tick_params(color='black', labelcolor='black', labelsize=8)
    cbar.set_label(label, fontsize=9, color='black')

    # outline
    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(0.6)

    return cbar

def save_mode_timeseries(time_model,
                         hemod, mode2_nodes, mode3_nodes, mode4_nodes,
                         mode1Layer, mode2Layer, mode3Layer, mode4Layer,
                         text_dir, save_tasks):
    """
    Save interface depths AND layer thicknesses for modes 0–3.
    """

    # mode 0 (h1)
    try:
        h0 = np.array([np.asarray(x).flatten()[0] for x in hemod])
        arr0 = np.column_stack((time_model, h0))
    except Exception:
        arr0 = None

    save_tasks.append((
        os.path.join(text_dir, "mode1_interfaces.txt"),
        arr0,
        "time(hour)\th0(m)",
        "%0.8f %0.5f"
    ))

    try:
        # shape (N,2)
        arr0L = np.column_stack((time_model, np.array(mode1Layer)))
    except Exception:
        arr0L = None

    save_tasks.append((
        os.path.join(text_dir, "mode1_layers.txt"),
        arr0L,
        "time(hour)\tlayer1(m)\tlayer2(m)",
        "%0.8f %0.5f %0.5f"
    ))


    # mode 1 (h1, h2)
    try:
        arr1 = np.column_stack((
            time_model,
            np.array(mode2_nodes)[:, 0],
            np.array(mode2_nodes)[:, 1]
        ))
    except Exception:
        arr1 = None

    save_tasks.append((
        os.path.join(text_dir, "mode2_interfaces.txt"),
        arr1,
        "time(hour)\th1(m)\th2(m)",
        "%0.8f %0.5f %0.5f"
    ))

    try:
        arr1L = np.column_stack((time_model, np.array(mode2Layer)))
    except Exception:
        arr1L = None

    save_tasks.append((
        os.path.join(text_dir, "mode2_layers.txt"),
        arr1L,
        "time(hour)\tlayer1(m)\tlayer2(m)\tlayer3(m)",
        "%0.8f %0.5f %0.5f %0.5f"
    ))


    # mode 2 (h1,h2,h3)
    try:
        arr2 = np.column_stack((
            time_model,
            np.array(mode3_nodes)[:, 0],
            np.array(mode3_nodes)[:, 1],
            np.array(mode3_nodes)[:, 2]
        ))
    except Exception:
        arr2 = None

    save_tasks.append((
        os.path.join(text_dir, "mode3_interfaces.txt"),
        arr2,
        "time(hour)\th1(m)\th2(m)\th3(m)",
        "%0.8f %0.5f %0.5f %0.5f"
    ))

    try:
        arr2L = np.column_stack((time_model, np.array(mode3Layer)))
    except Exception:
        arr2L = None

    save_tasks.append((
        os.path.join(text_dir, "mode3_layers.txt"),
        arr2L,
        "time(hour)\tlayer1(m)\tlayer2(m)\tlayer3(m)\tlayer4(m)",
        "%0.8f %0.5f %0.5f %0.5f %0.5f"
    ))


    # mode 2 (h1,h2,h3,h4)
    try:
        arr3 = np.column_stack((
            time_model,
            np.array(mode4_nodes)[:, 0],
            np.array(mode4_nodes)[:, 1],
            np.array(mode4_nodes)[:, 2],
            np.array(mode4_nodes)[:, 3]
        ))
    except Exception:
        arr3 = None

    save_tasks.append((
        os.path.join(text_dir, "mode4_interfaces.txt"),
        arr3,
        "time(hour)\th1(m)\th2(m)\th3(m)\th4(m)",
        "%0.8f %0.5f %0.5f %0.5f %0.5f"
    ))

    # mode 3 (h1,h2,h3,h4)
    try:
        arr3L = np.column_stack((time_model, np.array(mode4Layer)))
    except Exception:
        arr3L = None

    save_tasks.append((
        os.path.join(text_dir, "mode4_layers.txt"),
        arr3L,
        "time(hour)\tlayer1(m)\tlayer2(m)\tlayer3(m)\tlayer4(m)\tlayer5(m)",
        "%0.8f %0.5f %0.5f %0.5f %0.5f %0.5f"
    ))

def equation_sd(H, x, rhoe, rhoh, L, dh):
    """
    Degeneration Regime 1-2 (vectorized).
    """
    x = np.asarray(x, dtype=float)
    nu = 1.1e-6
    he = x * H
    # reduced gravity (use absolute to keep sign safe)
    glin = 9.81 * (rhoh - rhoe) / np.maximum(rhoe, 1e-6)
    cp = np.sqrt(np.maximum(glin * (1.0 - x) * he, 1e-12))

    aux1 = 2.0 * nu * L / (3.0 * H * cp)                    # shape (x,)
    aux2 = np.sqrt(np.pi * cp / (2.0 * nu * L))
    aux3 = (1.0 - x) ** 2
    aux4 = 1.0 / np.maximum(dh, 1e-9)
    aux5 = x * (1.0 - 2.0 * x)

    # avoid divide-by-zero in aux5 — set large value where aux5 is ~0
    safe_aux5 = np.where(np.abs(aux5) < 1e-9, np.sign(aux5) * 1e-9 + 1e-9, aux5)

    y = aux1 * (aux2 * aux3 + aux4) / safe_aux5
    # clip to finite range for plotting
    y = np.where(np.isfinite(y), y, np.nan)
    return y


def equation_bo(H, x):
    """
    Degeneration Regime for supercritical flow (bores).
    """
    x = np.asarray(x, dtype=float)
    aux1 = (1.0 - x) ** 2
    aux2 = x ** 3
    aux3 = (1.0 - x) ** 3
    # avoid negative inside sqrt due to numerical issues
    denom = aux2 + aux3
    denom = np.where(denom <= 0, 1e-12, denom)
    out = np.sqrt(aux1 / denom)
    return out


def equation_kh(H, x, dh):
    """
    Degeneration Regime for Kelvin-Helmholtz billows.
    """
    x = np.asarray(x, dtype=float)
    aux1 = 1.0 / np.maximum(x, 1e-12) - 1.0
    aux2 = dh / np.maximum(H, 1e-12)
    val = 2.0 * np.sqrt(np.maximum(aux1 * aux2, 0.0))
    return val

 
def boundaries (typ, x, y, z,lim):
    """
    Find upper and lower boundary limits.
    """                  
    if (typ == 'lower'):
        value = min(np.nanmin(x),np.nanmin(y),np.nanmin(z))-0.05
        if(value < lim):
            value = lim
    else:
        value = max(np.nanmax(x),np.nanmax(y),np.nanmax(z))+0.05
        if(value >lim):
            value= lim
            
    return value            

def degeneration(p1, p2, p3, xa, xb, xc,
                 ya, yb, yc, pe, ph, meta, H, ls, ax):
    """
    Averaged degeneration of BSIW (95% c.l.)
    """

    # Theoretical curves (vectorized)
    xl = np.arange(0.01, 0.501, 0.001)
    y1 = equation_sd(H, xl, pe, ph, ls, meta)
    y2 = equation_bo(H, xl)
    y3 = equation_kh(H, xl, meta)

    ax.plot(xl, y1, linewidth=1.0, color='blue', label='sd (dissipative-steepening)')
    ax.plot(xl, y2, linewidth=1.0, color='black', label='bo (bore limit)')
    ax.plot(xl, y3, linewidth=1.0, color='red', label='kh (Kelvin-Helmholtz)')

    # Compute centers and ellipse sizes robustly
    def _ellipse_params(x_arr, y_arr):
        x_arr = np.asarray(x_arr, dtype=float)
        y_arr = np.asarray(y_arr, dtype=float)
        x_arr = x_arr[np.isfinite(x_arr)]
        y_arr = y_arr[np.isfinite(y_arr)]
        if x_arr.size == 0 or y_arr.size == 0:
            return np.nan, np.nan, 0.01, 0.01
        cex = np.mean(x_arr)
        cey = np.mean(y_arr)
        
        # Width/height: use peak-to-peak or small default
        w = max(np.ptp(x_arr), 0.01)
        h = max(np.ptp(y_arr), 0.05)
        
        # Scale widths to represent spread nicely (factor 2 to emulate previous code)
        return cex, cey, max(0.01, 2.0 * w), max(0.05, 2.0 * h)

    cexa, ceya, wa, ha = _ellipse_params(xa, ya)
    cexb, ceyb, wb, hb = _ellipse_params(xb, yb)
    cexc, ceyc, wc, hc = _ellipse_params(xc, yc)

    # Create ellipses (ensure positive widths/heights)
    e1 = mpatches.Ellipse((cexa, ceya), wa, ha, color='red', alpha=0.25)
    e2 = mpatches.Ellipse((cexb, ceyb), wb, hb, color='navy', alpha=0.25)
    e3 = mpatches.Ellipse((cexc, ceyc), wc, hc, color='green', alpha=0.25)

    ax.add_patch(e1)
    ax.add_patch(e2)
    ax.add_patch(e3)

    # Legend entries
    def _period_label(p):
        try:
            return f"{p[0].day}/{p[0].month} - {p[1].day}/{p[1].month}"
        except Exception:
            return str(p)

    label1 = "P1 ~ " + _period_label(p1)
    label2 = "P2 ~ " + _period_label(p2)
    label3 = "P3 ~ " + _period_label(p3)

    proxy1 = mpatches.Patch(color='red', alpha=0.25, label=label1)
    proxy2 = mpatches.Patch(color='navy', alpha=0.25, label=label2)
    proxy3 = mpatches.Patch(color='green', alpha=0.25, label=label3)

    ax.legend(handles=[proxy1, proxy2, proxy3], loc='upper right', prop={'size': 9})

    ax.set_ylabel('1/W')
    ax.set_xlabel(r'$h_e/H$')

    ax.set_xlim([0.1, 0.5])
    ax.set_ylim([0.0, 2.0])

    ax.grid(True, which='both', ls=':', lw=0.4)



def degeneration_evolution(p1, p2, p3, x, y, pe, ph, meta, H, ls, zoom, ax):
    """
    Degeneration of Basin-Scale Internal Waves (BSIW) in time.
    """

    # --- Theoretical curves ---
    xl = np.linspace(0.01, 0.5, 500)
    ax.plot(xl, equation_sd(H, xl, pe, ph, ls, meta), c='blue', lw=1.2, label='SD Regime')
    ax.plot(xl, equation_bo(H, xl), c='black', lw=1.0, label='BO Regime')
    ax.plot(xl, equation_kh(H, xl, meta), c='red', lw=1.0, label='KH Regime')

    # --- Smoothed experimental points ---
    window_size = max(1, len(x[0]) // 10)
    x_smooth, y_smooth = [], []
    colors = ['red', 'navy', 'green']
    labels = ['P1', 'P2', 'P3']

    for i in range(3):
        xi = movingaverage(x[i], window_size)
        yi = movingaverage(y[i], window_size)
        if xi is not None and yi is not None:
            x_smooth.append(xi)
            y_smooth.append(yi)
            ax.scatter(xi, yi, s=5, color=colors[i], label=labels[i], zorder=3)

    # Axis limits 
    def safe_bounds(arrays, func, default):
        try:
            vals = np.concatenate([np.asarray(a) for a in arrays if a is not None])
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                return default
            return func(vals)
        except Exception:
            return default

    minix = safe_bounds(x_smooth, np.min, 0.1)
    maxix = safe_bounds(x_smooth, np.max, 0.5)
    miniy = safe_bounds(y_smooth, np.min, 0.0)
    maxiy = safe_bounds(y_smooth, np.max, 2.0)

    if zoom == 'no':
        # Rectangle showing zoom region
        rect = mpatches.Rectangle((minix, miniy),
                                  maxix - minix, maxiy - miniy,
                                  lw=1, ls='--', fill=False, color='gray', alpha=0.7)
        ax.add_patch(rect)
        ax.set_xlim(0.1, 0.5)
        ax.set_ylim(0.0, 2.0)
        ax.set_ylabel('1/W')
        ax.legend(loc='upper right', fontsize=9)
    else:
        ax.set_xlim(minix, maxix)
        ax.set_ylim(miniy, maxiy)

    # Common axis formatting
    ax.set_xlabel(r'$h_e / H$')
    ax.grid(True, which='both', ls=':', lw=0.5, color='gray')

    
def modal_period(date, period_time, ax):
    """
    Theoretical modal periods (hours) for each vertical mode.
    
    """

    colors = ['black', 'dimgray', '#1f77b4', '#d62728', '#9467bd']

    nmodes = min(period_time.shape[1], 5)

    for i in range(nmodes):
        ax.plot(date, period_time[:, i], linewidth=1.2, color=colors[i],
            ls='-', label=f"Mode V{i+1}")

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel("Time")

    ax.set_xlim([date[0], date[-1]])

    ax.grid(True, which="both", color='black', ls=":", lw=0.25)

    ax.set_ylabel('Wave period (hours)')
    ax.legend(loc='upper right', fontsize=8)

    
def velocity_mode(dx, vel, h, per, i, ax, fig=None):
    """
    Plot arbitrary velocity mode using symmetric colormap and zero-centered normalization.
    """

    # Ensure arrays
    dx = np.asarray(dx)
    vel = np.asarray(vel)
    h = np.asarray(h)

    # Parent figure (for inset colorbar)
    if fig is None:
        fig = ax.get_figure()

    # Depth reversed (consistent with other depth-based plots)
    mab = h[:-1]

    ax.set_title(
        f"mode V{i+1}H1  -  $\\bar{{T}}$ = {per:.2f} hours",
        loc="left", fontsize=10)

    # Zero-centered normalization
    norm = MidpointNormalize(midpoint=0.)

    vel = vel.T
    pcm = ax.contourf(dx, mab, vel, levels=100, cmap="seismic", 
                      extend="both", norm=norm)

    ax.set_ylabel("Height above the bottom (m)")
    ax.set_xlabel("Time")

    # Format horizontal axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel('Time')

    cax = inset_axes( ax, width="30%", height="8%", loc="upper right",
        borderpad=1)

    cbar = fig.colorbar(pcm, cax=cax, orientation="horizontal")
    cbar.formatter = ticker.FormatStrFormatter('%.2f')
    cbar.locator = ticker.MaxNLocator(nbins=5)
    cbar.set_label("Arbitrary velocity (-)", fontsize=9)
    cbar.ax.tick_params(labelsize=7)
     
     

def maxi_significance(sig,maximum):
    """
    Maximum value for y-axis (analysis between PSD and significance curve)
    """    
    sig_max=maximum
    for i in range(4):
        if len(sig[i]) > 0:
            if sig_max < max(sig[i]):
                sig_max = max(sig[i])
    
    return sig_max
            

    
def model_plot(freq, col, maxi, ax, typ='log'):
    """
    Draw model-indicating rectangles inside a PSD plot.
    Marks frequency positions or bands for visual reference.
    """
    
    maxi = maxi * 1e3
    d = freq[0] - freq[2]  # frequency interval

    rec = mpatches.Rectangle((freq[2], 0.0), d, maxi, color=col, alpha=0.5)

        
    ax.add_patch(rec)
    
def significance(sig,wo,c,ax,typ='xax'):
    """
    Plot confidence interval indicator for PSD.
    """
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

def psd_iso(tau, freq, psd, v1h1, v2h1, v3h1, largelen, n, fo, wo, sig, ax):
    """
    Power Spectral Density of isotherms (m²/Hz)
    """
    color = ['black', 'dimgray', 'salmon', 'red']

    maxima = 0.0

    minimox = np.inf
    maximox = 1e-12

    # psd and freq are expected as lists/arrays per isotherm
    for i in range(4):
        if tau[i] != -999:
            ff = np.asarray(psd[i])
            f = np.asarray(freq[i])
            s = sig[i]
            w = wo[i]

            if largelen == 1:
                l0 = max(1, int(len(ff) / 150))
                plot_y = smooth(ff, l0)
                
            else:
                plot_y = ff

            ax.plot(f, plot_y, linewidth=1, c=color[i], ls='-', 
                    label=f"{tau[i]}°C")

            max_y = valmax(plot_y, f)
            maxima = max(maxima, max_y)

            # update freq limits (log-friendly)
            fmin = np.nanmin(f[np.isfinite(f)])
            fmax = np.nanmax(f[np.isfinite(f)])

            minimox = min(minimox, 10 ** (np.floor(np.log10(fmin))))
            maximox = max(maximox, 10 ** (np.ceil(np.log10(fmax))))

            # significance
            significance(w, s, color[i], ax)


    ax.set_xscale('log')
    ax.set_yscale('log')

    # Determine PSD limits based on data
    psd_min = np.nanmin([np.nanmin(np.asarray(p)) for p in psd if isinstance(p, (list, np.ndarray))])
    if not np.isfinite(psd_min) or psd_min <= 0:
        psd_min = 1e-8  # safe fallback
    
    ytop = 10 ** np.ceil(np.log10(maxima * 1.2))
    ybottom = 10 ** np.floor(np.log10(psd_min / 10))  # one decade below min PSD
    
    ax.set_ylim(ybottom, ytop)
    ax.set_xlim(minimox, maximox)
    
    ax.grid(True, which='both', color='black', ls=':', lw=0.25)


    try:
        model_plot(v1h1, 'green', ytop, ax)
    except Exception:
        pass
    try:
        model_plot(v2h1, 'blue', ytop, ax)
    except Exception:
        pass
    try:
        model_plot(v3h1, 'red', ytop, ax)
    except Exception:
        pass

    # aumentar right limit de acordo com significância
    maxima_sig = maxi_significance(sig, maxima)
    if np.isfinite(maxima_sig) and maxima_sig > 0:
        ax.set_ylim(top=10 ** np.ceil(np.log10(max(maxima_sig, maxima) * 1.2)))

    # linhas verticais fo e n até top do eixo
    ymax_current = ax.get_ylim()[1]
    ax.plot([fo, fo], [1e-12, ymax_current], linewidth=1, c='black', ls='--')
    ax.plot([n, n], [1e-12, ymax_current], linewidth=1, c='black', ls='--')

    ax.legend(loc='lower left')


def psd_variance(tau, freq, psd, n, fo, wo, sig, ax):
    """
    Plot Function: Variance-preserving PSD (m²).
    """
    color = ['black', 'dimgray', 'salmon', 'red']
    maxima = 0.0
    minimox = np.inf
    maximox = 1e-12

    for i in range(4):
        if tau[i] != -999:
            ff = np.asarray(psd[i])
            f = np.asarray(freq[i])
            s = sig[i]
            w = wo[i]

            yplot = ff * f  # variance-preserving (PSD * f)
            significance(s, w, color[i], ax, 'yax_var')

            ax.plot(f, yplot, linewidth=1, c=color[i], ls='-', label=f"{tau[i]}°C")

            max_y = max(yplot)
            maxima = max(maxima, max_y)

            # update freq limits (log-friendly)
            fmin = np.nanmin(f[np.isfinite(f)])
            fmax = np.nanmax(f[np.isfinite(f)])

            minimox = min(minimox, 10 ** (np.floor(np.log10(fmin))))
            maximox = max(maximox, 10 ** (np.ceil(np.log10(fmax))))

    ax.set_xscale('log')
    ymin, ymax = ax.get_ylim()

    ax.set_ylabel('Variance-preserving PSD (m²)')
    ax.set_xlabel('Frequency (Hz)')
        
    
    #maxima_sig = maxi_significance(sig, maxima)
    maxima_y = maxima * 1.2
    ax.set_ylim(0,maxima_y) 
    ax.set_xlim(minimox, maximox)

    ax.plot([fo, fo], [0, maxima_y], linewidth=1, c='black', ls='--')
    ax.plot([n, n], [0, maxima_y], linewidth=1, c='black', ls='--')

    ax.grid(True, which='both', color='black', ls=':', lw=0.25)
    ax.legend(loc='upper left')
    
    
def coherence(c1w, f1w, tau, iso, largelen, ax):
    """
    Coherence between isotherms and wind speed.
    """
    # smoothing if necessary
    if largelen == 1:
        l10 = max(1, int(len(c1w) / 100))
        yplot = smooth(c1w, l10)
    else:
        yplot = c1w

    ax.plot(f1w, yplot, linewidth=1, c='navy', ls='-', label=f'{tau[iso]} °C')
    ax.set_ylim(0, 1)
    
    ax.set_ylabel('Wind coherence (-)')
    ax.grid(True, which='both', color='black', ls=':', lw=0.25)
    ax.legend(loc='upper right')


def psd_depth(sen, matrix, freq, psd, v1h1, v2h1, v3h1,
              largelen, n, fo, ax):
    """
    Power Spectral Density of temperature variation (°C²/Hz)
    """

    colors = ['black', 'dimgray', 'salmon', 'red']

    # Limits to update
    psd_max = 0.0
    psd_min = np.inf
    freq_min = np.inf
    freq_max = 0.0

    # Loop sensors 
    for i in range(4):
        if sen[i] != 1:
            continue

        depth_i = np.nanmean(matrix[i])
        ff = np.asarray(psd[i], dtype=float)
        fr = np.asarray(freq[i], dtype=float)


        # Skip invalid or empty PSDs
        if ff.size == 0 or fr.size == 0:
            continue

        # Remove invalid/zero values for log-scale
        mask = np.isfinite(ff) & np.isfinite(fr) & (ff > 0) & (fr > 0)
        if np.sum(mask) == 0:
            continue

        ff = ff[mask]
        fr = fr[mask]

        # Smooth if requested
        if largelen == 1:
            l0 = max(1, int(len(ff) / 150))
            ff_plot = smooth(ff, l0)
        else:
            ff_plot = ff

        # Now compute limits from plotted PSD
        psd_max = max(psd_max, np.nanmax(ff_plot))
        psd_min = min(psd_min, np.nanmin(ff_plot))

        # Frequency bounds
        freq_min = min(freq_min, np.nanmin(fr))
        freq_max = max(freq_max, np.nanmax(fr))

        # PSD plot 
        ax.plot(
            fr, ff_plot,
            linewidth=1,
            color=colors[i],
            label=f"{round(depth_i,1)} m"
        )


    # Safe fallbacks
    if not np.isfinite(freq_min) or freq_min <= 0:
        freq_min = 1e-6
    if not np.isfinite(freq_max):
        freq_max = 1e-1
    if not np.isfinite(psd_min) or psd_min <= 0:
        psd_min = psd_max / 1e6
    if psd_max <= 0:
        psd_max = 1e-6

    # Log scaling limits
    left_lim   = 10 ** np.floor(np.log10(freq_min))
    right_lim  = 10 ** np.ceil(np.log10(freq_max))

    bottom_lim = 10 ** np.floor(np.log10(psd_min * 0.8))
    top_lim    = 10 ** np.ceil(np.log10(psd_max * 1.3))

    ax.set_yscale("log")   # PSD
    ax.set_xscale("log")   # Frequency

    ax.set_xlim(left_lim, right_lim)
    ax.set_ylim(bottom_lim, top_lim)

    # Model curves
    try: model_plot(v1h1, 'green', top_lim, ax)
    except: pass
    try: model_plot(v2h1, 'blue', top_lim, ax)
    except: pass
    try: model_plot(v3h1, 'red', top_lim, ax)
    except: pass

    # Vertical lines: Buoyancy frequency (n) and Coriolis (fo)
    ax.plot([fo, fo], [bottom_lim, top_lim], color='black', ls='--', lw=1)
    ax.plot([n,  n ], [bottom_lim, top_lim], color='black', ls='--', lw=1)

    ax.grid(True, which='both', color='black', ls=':', lw=0.25)

    ax.legend(loc='lower left', fontsize=8)

def smooth(y, box_pts):
    """
    Smooth time-series data.
    """
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def smooth_date(t,dt,y,limite=600):
    """
    Smooth time (date) signal.
    """  
    new_date = ["" for x in range(limite)]
    
    t_sm     = np.array(dt)
    t_smooth = np.linspace(np.min(t_sm), np.max(t_sm), limite)
    index    = value_nearest_vector(t_sm,t_smooth,limite)
    y_smooth = np.interp(t_smooth, dt, y)
    
    for x in range(limite):
        new_date[x] = t[int(index[x])]
        
    return new_date, t_smooth, y_smooth

def value_nearest_vector(array,value,lim):
    """
    Find nearest value in an array.
    """
    idx  = np.zeros(lim,float)
    
    for x in range(lim):
        idx[x] = (np.abs(array-value[x])).argmin()
        
    return idx

def interpolation(y1,y2,x1,x2,x):
    """
    Interpolate values from x and y arrays.
    Returns interpolated values based on input coordinates.
    """
    if(x1==x2):
        return abs(np.mean([y1,y2]))
    
    a = (y1-y2)/(x1-x2)
    b = y1 - a*x1
    
    y = a*x + b
    return y


def correction_contourf(h,num,dj,lin,dci,limit):
    """
    Transform 2D elevation array into 1D profile for plotting.
    Interpolates and reshapes elevation data for graph representation.
    """
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

def isotherm(date, t, y, color, iso, ax):
    """
    Isotherm variation over time.
    """
    # Extract scalar label safely (take first element if array)
    if isinstance(iso, (list, np.ndarray)):
        iso_val = float(np.ravel(iso)[0])
    else:
        iso_val = float(iso)

    # plot the isotherm time series (date vs y)
    ax.plot(date, y, linewidth=1.2, color=color, label=f'{iso_val:.1f} °C')
    ax.grid(True, which="both", color='gray', ls=":", lw=0.4)

    # y formatting
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax.set_ylabel('Isotherm depth (m)')

    # small legend
    ax.legend(fontsize=9, loc='upper right')
    
def wavelet_iso(date, per, power, ax, fig=None, mode='iso'):
    """
    Wavelet power spectrum of isotherm depth signal.
    """
    # ensure arrays are numpy
    date = np.asarray(date)
    per = np.asarray(per)
    power = np.asarray(power)

    # contour: expect power.shape == (len(per), len(date))
    pcMap = ax.contourf(date, per, power, levels=100, cmap='plasma', extend='both')

    ax.set_yscale('log')
    ax.set_ylim(0.5, 40)
    ax.set_ylabel('Scale (hours)')
    ax.set_xlabel('Time')
    
    # x formatting (show days, rotation)
    ax.set_xlim(date[0], date[-1])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)

    # get parent figure
    if fig is None:
        fig = ax.get_figure()
   
    if mode == 'iso':
        cax = inset_axes(ax, width="20%", height="8%", loc="upper right",
        borderpad=1)

        cbar = fig.colorbar(pcMap, cax=cax, orientation='horizontal')
        cbar.set_label('Wavelet Power (m²/Hz)', fontsize=9, color='white')
        cbar.ax.tick_params(labelsize=7, colors='white')
        cbar.locator = ticker.MaxNLocator(nbins=4)
    else:
        cax = inset_axes(ax, width="30%", height="4%", loc="upper right",
        borderpad=1)

        cbar = fig.colorbar(pcMap, cax=cax, orientation='horizontal')
        cbar.set_label('Wavelet Power (m²/Hz)', fontsize=8, color='white')
        cbar.ax.tick_params(labelsize=7, colors='white')
        cbar.locator = ticker.MaxNLocator(nbins=4)        
    
    cbar.outline.set_edgecolor('white')
    for spine in cbar.ax.spines.values():
        spine.set_edgecolor('white')


def thermocline(date, t, y, color, thermo, ax):
    """
    Thermocline temperature fluctuation.
    """
    # ensure arrays
    date = np.asarray(date)
    y = np.asarray(y, dtype=float)

    ax.plot(date, y, linewidth=1.2, color=color, label=f'Mean position: {thermo:.2f} m')
    ax.grid(True, which="both", color='gray', ls=":", lw=0.4)

    # x formatting: show dates but rotation/hide handled by caller
    ax.set_xlim(date[0], date[-1])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)

    # y formatting
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax.set_ylabel('Thermocline temperature (°C)')

    # small legend at top-right
    ax.legend(fontsize=9, loc='upper right')


def psd_thermocline(p, psd, color, depth, wo, sig, ax):
    """
    PSD of thermocline temperature.
    """

    ax.plot(psd, p, linewidth=1.2, color=color, ls='-')

    try:
        significance(sig, (1.0 / wo) / 3600.0, 'red', ax, typ='xax')
    except Exception:
        pass

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, which="both", color='gray', ls=":", lw=0.4)

    ax.set_xlabel(f'PSD {depth} m (°C²/Hz)')

    
def wind_direction(date, time, y, ys, ylow, ydi, wedd_up, wedd_low, ax):
    """
    Plot wind direction filtered by Wedderburn number.
    """

    # Base signal
    ax.plot(date, y, c='lightgray', lw=1, label='_nolegend_')

    # Homogeneous events
    ax.plot(date, ydi, c='blue', lw=3, label='Homogeneous wind events')

    # Filtered Wedderburn ranges
    ax.plot(date, ys, c='green', lw=1, label=f'20 < W < {wedd_low:.1f}')
    ax.plot(date, ylow, c='red', lw=1, label=f'{min(wedd_up, 1):.1f} < W < 20')

    # Formatting
    ax.set_ylabel(r'Wind direction $(^{\circ})$', color='gray')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.set_xlabel('Time')
    
    ax.set_xlim(date[0], date[-1])
    ax.set_ylim(0, 360)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, ls=':', lw=0.4, color='gray')
     
        
def wavelet_depth(date, per, power, time, ax, fig=None):
    """
    Wavelet power spectrum of temperature sensor signal.
    """

    # ensure numpy arrays
    date = np.asarray(date)
    per = np.asarray(per)
    power = np.asarray(power)

    # parent figure
    if fig is None:
        fig = ax.get_figure()

    # main contour
    pcm = ax.contourf(date, per, power, levels=100, cmap="plasma",
        extend="both")

    # Axes scaling + formatting
    ax.set_yscale("log")
    ax.set_ylim(0.5, 40)
    ax.set_ylabel("Scale (hours)")

    ax.set_xlim(date[0], date[-1])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.tick_params(axis="x", rotation=25)
    ax.set_xlabel("Time")

    cax = inset_axes(ax, width="30%", height="8%", loc="upper right", borderpad=1)

    cbar = fig.colorbar(pcm, cax=cax, orientation="horizontal")
    cbar.set_label("Wavelet Power (°C²/Hz)", fontsize=9, color='white')
    cbar.ax.tick_params(labelsize=7, colors='white')
    
    
    cbar.outline.set_edgecolor('white')
    for spine in cbar.ax.spines.values():
        spine.set_edgecolor('white')

def wavelet_resona(dx, per, power, v1, v2, ax, fig=None):
    """
    Wavelet resonance plot (wind vs modeled internal modes)

    """
    # Ensure arrays
    dx = np.asarray(dx)
    per = np.asarray(per)
    power = np.asarray(power)

    # Wavelet contour
    pcm = ax.contourf(dx, per, power, levels=100,
                      cmap='plasma', extend='both')

    # Model modes
    ax.plot(dx, v1, lw=1.0, c='black', ls='-',  label="V1H1")
    ax.plot(dx, v2, lw=1.0, c='red', ls='-', label="V2H1")
    ax.legend(loc='lower left', prop={'size': 9})

    # Axes formatting 
    ax.set_yscale('log')
    ax.set_ylim(0.5, 40)
    ax.set_ylabel("Scale (hours)")
    ax.set_xlim(dx[0], dx[-1])

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel('Time')

    # Figure handle 
    if fig is None:
        fig = ax.get_figure()

    cax = inset_axes(
        ax,
        width="30%",       # relative to axes
        height="8%",
        loc='upper right',
        borderpad=1.2
    )

    cbar = fig.colorbar(pcm, cax=cax, orientation='horizontal')
    cbar.set_label('Wavelet (m²/(s² Hz))', fontsize=9, color='white')
    cbar.ax.tick_params(labelsize=7, colors='white')
    
    
    cbar.outline.set_edgecolor('white')
    for spine in cbar.ax.spines.values():
        spine.set_edgecolor('white')


def temperature(date, t, y, color, depth, ax):
    """
    Plot Function: Temperature variation.
    """
    mean_depth = float(np.nanmean(depth))
    lab = f"At {mean_depth:.2f} m "
    ax.plot(date, y, linewidth=1.2, color=color, label=lab)
    ax.legend(loc='upper right', fontsize=9)
    
    # X axis formatting
    ax.set_xlim(date[0], date[-1])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis="x", rotation=25)

    # Y axis formatting
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax.set_ylabel("Temperature (°C)")

    # Grid
    ax.grid(True, which="both", color="gray", ls=":", lw=0.4)

def multi_isotherms(date, y, iso, time, z0, ax):
    """
    Iostherm variation
    Temporal variation of all analyzed isotherms.
    """

    colors = ['black', 'dimgray', 'salmon', 'red']

    # Plot valid isotherms only
    for i in range(4):
        if iso[i] != -999:
            ax.plot(date, y[i], lw=1, color=colors[i],
                    label=f"{iso[i]}°C")

    # Format horizontal axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel('Time')

    ax.set_xlim(date[0], date[-1])
    ax.grid(True, which="both", color='black', ls=":", lw=0.25)
    ax.legend(loc='upper right', fontsize=9)

    if z0 == 0:
        ax.set_ylabel('Height above bottom (m)')
    else:
        ax.set_ylabel('Altitude (m)')
        ax.set_title(f'Reference level for bottom: {z0} m',
                     loc='left', fontsize=9)

def temp_bandpass(date, y, iso, ax, fig=None):
    """
    Band-pass filtered isotherm displacement
    """

    color = ['black', 'dimgray', 'salmon', 'red']

    # Ensure parent figure
    if fig is None:
        fig = ax.get_figure()

    # Main plot 
    for i in range(4):
        if iso[i] != -999:
            try:
                ax.plot(
                    date,
                    y[i],
                    linewidth=1,
                    color=color[i],
                    label=f"{iso[i]}°C"
                )
            except Exception:
                pass

    # Axis formatting 
    ax.set_xlim(date[0], date[-1])
    ax.set_ylabel("Displacement (m)")

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel("Time")

    ax.grid(True, which="both", color='black', ls=":", lw=0.3)
    ax.legend(loc='upper right', fontsize=9)


def depth_bandpass(date, y, depth, time, s, ax, fig=None):
    """
    Bandpass filtered thermistor sensors
    Shows filtered thermal variation at each sensor depth.
    """

    # Ensure parent figure
    if fig is None:
        fig = ax.get_figure()

    colors = ['#d62728', '#8c1b13', '#1f77b4', '#0b1d51']  # red, maroon, blue, navy

    # Plot active sensors 
    for i in range(4):
        if s[i] != 1:
            continue

        if depth[i] is None:
            continue

        try:
            label_depth = round(float(np.nanmean(depth[i])), 0)
            ax.plot(date, y[i], linewidth=1.2, color=colors[i],
                label=f"{label_depth} m")
        except Exception:
            pass

    # Axis formatting
    ax.set_xlim(date[0], date[-1])
    ax.set_ylabel("Temperature variation (°C)")

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel('Time')

    ax.grid(True, which="both", color='black', ls=":", lw=0.25)
    ax.legend(loc='upper right', fontsize=9)

   
def averageTemperature(temp, low, lowsd, high, highsd, h, z0, water,ax):
    """
    Plot averaged temperature profile with confidence intervals and standard 
    deviations.
    """

    # Plot mean profile
    ax.plot(temp, h, color='black', linewidth=1, label='Mean temperature')
    ax.scatter(temp, h, color='red', label='Mean points')

    # Fill 95% confidence interval
    ax.fill_betweenx(h, low, high, where=high >= low,
                     facecolor='red', alpha=0.2, label='95% CI', 
                     interpolate=True)

    # Fill +/- 1 standard deviation
    ax.fill_betweenx(h, lowsd, highsd, where=highsd >= lowsd,
                     facecolor='blue', alpha=0.15, label='±1 SD', 
                     interpolate=True)

    # Grid, labels, and title
    ax.grid(True, which="both", color='black', ls=":", lw=0.25)
    ax.set_xlabel('Water temperature (°C)')
    
    xlims = ax.get_xlim()  # get current x-axis limits
    ax.plot(xlims, [water, water], lw=1, ls='-', color='blue', label='Surface')
    
    if z0 == 0:
        ax.set_ylabel('Height above bottom (m)', color='black')
    else:
        ax.set_ylabel('Altitude (m)',  color='black')
        ax.plot(xlims, [z0, z0], lw=1, ls='--', color='black', label='Bottom')
    
    ax.legend(loc='upper left', prop={'size': 9})
 
def averageBuoyancy(buoy, low, lowsd, high, highsd, h, z0, water, ax):
    """
    Plot averaged buoyancy profile with confidence intervals and standard 
    deviations.
    """

    # Mean profile
    ax.plot(buoy, h, color='black', linewidth=1, label='Mean buoyancy')
    ax.scatter(buoy, h, color='red', s=12, label='Mean points')

    # 1 standard deviation
    ax.fill_betweenx(
        h, lowsd, highsd, where=highsd >= lowsd,
        facecolor='blue', alpha=0.15, label='±1 SD', interpolate=True)

    # Mark region where SD goes below zero (dark red)
    below_zero = lowsd < 0
    if np.any(below_zero):
        ax.fill_betweenx(
            h, lowsd, 0, where=below_zero,
            facecolor='darkred', alpha=0.4,  interpolate=True)

    # 95% confidence interval (on top of SD)
    ax.fill_betweenx(
        h, low, high, where=high >= low,
        facecolor='red', alpha=0.25, label='95% CI', interpolate=True)

    ax.grid(True, which="both", color='black', ls=":", lw=0.25)
    ax.set_xlabel('Buoyancy frequency (Hz)')

    xlims = ax.get_xlim()  # get current x-axis limits
    ax.plot(xlims, [water, water], lw=1, ls='-', color='blue', label='Surface')

    if z0 != 0:
        ax.plot(xlims, [z0, z0], lw=1, ls='--', color='black', label='Bottom')



def thermal_variation(date, depth, depth_point, temp, time, ax):
    """
    Thermal fluctuation        
    Shows time-varying temperature at selected depth levels.

    """

    # Base temperature curve (average or reference)
    ax.plot(date, temp, linewidth=1.2, color="black")

    colors = ["#d62728", "#8c1b13", "#1f77b4", "#0b1d51"]  # red, maroon, blue, navy

    # Plot temperatures at selected depths
    for i in range(4):

        # skip invalid entries
        if depth[i] is None or depth_point[i] is None:
            continue

        idx = depth_point[i]
        if idx < 0 or idx >= temp.shape[1]:
            continue

        series = temp[:, idx]
        label_val = round(np.nanmean(depth[i]), 1)

        ax.plot(date, series, linewidth=1.2, color=colors[i],
            label=f"{label_val} m")

    ax.grid(True, which="both", color="black", ls=":", lw=0.25)

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.tick_params(axis="x", rotation=25)
    ax.set_xlabel("Time")

    ax.set_xlim(date[0], date[-1])
    ax.set_ylabel("Water temperature (°C)")

    ax.legend(loc="upper right", fontsize=9)

            

def psd_multilayer(f, depth, ff, fi, n2, v1, v2, z0, ax, fig=None):
    """
    Multi-layer PSD of temperature fluctuation.
    ff must be shape (nz, nf); depth is 1D (nz); f is 1D (nf)
    """

    # Ensure numpy arrays
    f = np.asarray(f)
    depth = np.asarray(depth)
    ff = np.asarray(ff)

    # Parent figure
    if fig is None:
        fig = ax.get_figure()

    # Avoid log10(-inf) warnings
    ff_safe = np.where(ff > 0, ff, np.nan)
    log_ff = np.log10(ff_safe[::-1, :])  # reverse depth for plotting

    # Main PSD contour
    pcm = ax.contourf(f, depth, log_ff, levels=120, cmap="plasma", extend="both")

    cax = inset_axes(ax, width="20%", height="4%", loc="upper center", borderpad=1)

    cbar = fig.colorbar(pcm, cax=cax, orientation="horizontal")
    cbar.set_label("PSD (log₁₀ °C²/Hz)", fontsize=8, color='white')
    cbar.ax.tick_params(labelsize=7, colors='white')
    cbar.locator = ticker.MaxNLocator(nbins=4)
    
    
    cbar.outline.set_edgecolor('white')
    for spine in cbar.ax.spines.values():
        spine.set_edgecolor('white')

    ax.plot([fi, fi], [depth.min(), depth.max()],
            color="black", ls="-", lw=1, label="Forcing frequency")

    ax.plot([n2, n2], [depth.min(), depth.max()],
            color="black", ls="--", lw=1, label="Buoyancy frequency")

    ax.plot([1/v1, 1/v1], [depth.min(), depth.max()],
            color="green", ls="--", lw=1, label="V1H1")
    
    ax.plot([1/v2, 1/v2], [depth.min(), depth.max()],
            color="blue", ls="--", lw=1, label="V2H1")
    
    ax.legend(loc='lower right',prop={'size': 9})



    ax.set_xscale("log")
    ax.set_xlabel("Frequency (Hz)")

    # Y-label depends on reference
    if z0 == 0:
        ax.set_ylabel("Height above bottom (m)")
    else:
        ax.set_ylabel("Altitude (m)")
        ax.set_title(f"Reference bottom level: {z0} m", loc="left", fontsize=9)

    ax.grid(True, which="both", ls=":", lw=0.3, color="gray")

def valmax(psd,f):
    """
    Find maximum spectral energy.
    Identifies the peak energy value in the analyzed spectrum.
    """
    maxi = 0    
    l    = len(psd)
    
    for i in range(l):    
        if(f[i]>6*10**-6 and psd[i]>maxi):
            maxi = psd[i]
            
    return maxi
        

def psd_wind(period, psd_wind, wo, sig, ax):
    """
    Power Spectral Density (PSD) of wind speed.
    """
    ax.plot(period[1:], psd_wind[1:], c='navy', lw=1.2, label='Wind')
    significance(sig, (1 / wo) / 3600, 'navy', ax, typ='yax')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(period.min(), period.max())
    ax.grid(True, which="both", color='gray', ls=":", lw=0.4)

    ax.tick_params(axis='y', labelcolor='navy')
    ax.set_ylabel('PSD wind ((m/s)²/Hz)', color='navy')
    ax.set_xlabel('Period (hours)')


    
def psd_sola(period, psd_sola, wo, sig, ax):
    """
    Power Spectral Density (PSD) of solar radiation.
    """
    ax.plot(period[1:], psd_sola[1:], c='red', lw=1.2, label='Solar')
    significance(sig, (1 / wo) / 3600, 'red', ax, typ='yax')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(period.min(), period.max())
    ax.grid(True, which="both", color='gray', ls=":", lw=0.4)

    ax.tick_params(axis='y', labelcolor='red')
    ax.set_ylabel('PSD solar radiation ((W/m²)²/Hz)', color='red')
    ax.set_xlabel('Periodicity (hours)')


def psd_level(f, ff, wo, sig, ax):
    """
    PSD of water elevation
    """
    ax.plot(f, ff, linewidth=1, c='black', ls='-')
    significance(sig, wo, 'red', ax, typ='yax')

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_ylabel('PSD water elevation (m²/Hz)') 
    ax.set_xlabel('Frequency (Hz)')
    ax.tick_params(axis='x', labelcolor='black')
    ax.tick_params(axis='y', labelcolor='black')
    ax.grid(True, which="both", color='black', ls=":", lw=0.25)



def code_index(r):
    """
    Define isotherm combination code.
    Maps isotherm pairs to a corresponding index value.
    """
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

def coherence_iso(t, coh, fre, ana1, ana2, ax):
    """
    Coherence between isotherms.
    """

    colors = ['black', 'dimgray', 'salmon', 'red', 'navy', 'blue']

    ax.grid(axis='y', color='black', ls=":", lw=0.25)
    ax.set_ylim(0, 1.01)
    ax.set_ylabel('coherence', color='gray')
    ax.set_xlabel('period (hours)')

    # Decode selected pairs only once (performance improvement)
    a1, b1 = int(ana1 // 10), int(ana1 % 10)
    a2, b2 = int(ana2 // 10), int(ana2 % 10)

    # Loop i<j only
    for i in range(4):
        if t[i] == -999:
            continue

        for j in range(i+1, 4):
            if t[j] == -999:
                continue

            # Check if this pair is requested
            if not (
                (a1 == i+1 and b1 == j+1) or
                (a2 == i+1 and b2 == j+1)
            ):
                continue
            

            k = code_index((i+1) * (j+1))
            period_h = (1.0 / fre[k]) / 3600.0

            ax.plot(period_h, coh[k], lw=1,
                    color=colors[k],
                    label=f"{t[i]}/{t[j]}°C")

    ax.legend(loc='upper right', fontsize=9)
    


def phase_iso(t, ph, fre, conf, ana1, ana2, ax):
    """
    Phase-shift between isotherms
    Graph to plot the phase-shift between isotherms.    """

    colors = ['black', 'dimgray', 'salmon', 'red', 'navy', 'blue']

    ax.grid(axis='y', color='black', ls=":", lw=0.25)
    ax.set_ylim([0, 1.01])

    # Decode the isotherm pairs selected in ana1 and ana2
    a1 = int(np.around(float(ana1) / 10.0, 0))
    b1 = ana1 - 10 * a1
    a2 = int(np.around(float(ana2) / 10.0, 0))
    b2 = ana2 - 10 * a2

    # Loop over all isotherm combinations
    for i in range(4):
        for j in range(i + 1, 4):

            # Only plot the pair selected
            if not (
                (a1 == i + 1 and b1 == j + 1) or
                (a2 == i + 1 and b2 == j + 1)
            ):
                continue

            # Skip invalid isotherms
            if t[i] == -999 or t[j] == -999:
                continue

            k = code_index((i + 1) * (j + 1))
            f = np.asarray(fre[k])
            p = np.asarray(ph[k])
            idxs = conf[k]

            # -ormalize idxs into a 1D integer index array 
            if isinstance(idxs, (tuple, list, np.ndarray)):
                idxs_arr = np.asarray(idxs)
                # If boolean mask, convert to integer indices
                if idxs_arr.dtype == bool:
                    idxs_arr = np.where(idxs_arr)[0]
                # If array of shape (n,1) or similar, ravel it
                idxs_arr = idxs_arr.ravel()
                if idxs_arr.size == 0:
                    continue
                # keep only integer-like entries
                try:
                    idxs_int = np.array([int(x) for x in idxs_arr])
                except Exception:
                    # if conversion fails, skip
                    continue
            else:
                # single scalar
                try:
                    idx_val = int(idxs)
                    idxs_int = np.array([idx_val])
                except Exception:
                    continue

            # Clip indices to valid range and remove negatives
            valid_mask = (idxs_int >= 0) & (idxs_int < len(f))
            idxs_int = idxs_int[valid_mask]
            if idxs_int.size == 0:
                continue

            # Compute periods (hours) and phases, then plot all significant points
            periods_h = (1.0 / f[idxs_int]) / 3600.0
            phases = np.abs(p[idxs_int])

            ax.scatter(periods_h, phases, marker='x', color=colors[k], label=f"{t[i]}/{t[j]}°C")

    ax.set_ylabel('Phase (radian)', color='black')
    ax.set_xlabel('Period (hours)')

    ax.set_yticks([0, np.pi / 2, np.pi])
    ax.set_yticklabels([r'$0$', r'$\frac{\pi}{2}$', r'$\pi$'])

    ax.grid(axis='y', color='black', ls=":", lw=0.25)
    ax.legend(loc='upper right', prop={'size': 9})



def density(dx,t,pu,pd,pe,ph,ax):
    """
    Water Density Structure Variation
    """

    ax.plot(dx, pu, color='darkred', linewidth=1.2, label='Surface')
    ax.plot(dx, pd, color='darkblue', linewidth=1.2, label='Bottom')
    ax.plot(dx, pe, color='red', linewidth=1.0, linestyle='--', label='Epilimnion')
    ax.plot(dx, ph, color='blue', linewidth=1.0, linestyle='--', label='Hypolimnion')

    ax.invert_yaxis()  
    ax.set_ylabel(r'Water density (kg m$^{-3}$)')
    ax.grid(True, linestyle=':', color='0.7', linewidth=0.5)

    ax.set_xticklabels(dx, rotation=25)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.set_xlabel('Time')

    ax.legend(loc='lower right', fontsize=9, ncol=2)

def radiation(date,t,y,ax):
    """
    Plot solar radiation signal.
    Generates a graph of the solar radiation time series.
    """
    ax.plot(date, y, linewidth=1, c='darkgoldenrod', ls='-')
    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    ax.set_ylabel('solar rad.(W/m²)')

    
def wind(date, time, dire, velo, ax):
    """
    Wind speed and direction (dual-axis)
    """
    # Basic checks
    if len(date) == 0 or len(velo) == 0:
        ax.text(0.5, 0.5, 'No wind data', ha='center', va='center', transform=ax.transAxes)
        return

    # Plot wind speed
    ax.plot(date, velo, color='navy', lw=1.2, label='Speed')
    ax.set_xlim(date[0], date[-1])
    ax.set_ylim(0, np.nanmax(velo) * 1.2)
    ax.set_ylabel('Wind speed (m/s)', color='navy', fontsize=9)
    ax.tick_params(axis='y', labelcolor='navy')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel('Time')

    # Wind direction on twin axis
    ax2 = ax.twinx()
    try:
        # Reuse your existing smoother if available
        t_smooth, t_numeral, y_smooth = smooth_date(date, time, dire, int(len(date) * (time[1]-time[0])))
        ax2.plot(t_smooth, y_smooth, color='black', lw=1.2, ls='-', label='Direction')
    except Exception:
        ax2.plot(date, dire, color='black', lw=1.0, alpha=0.6, ls='-', label='Direction')

    ax2.set_ylabel('Wind direction (°)', color='black', fontsize=9)
    ax2.tick_params(axis='y', labelcolor='black')
    ax2.set_ylim(-20, 380)
    ax2.set_yticks([0, 90, 180, 270, 360])

    # Combine legend
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=9)

    ax.grid(True, which='major', ls=':', lw=0.5, alpha=0.7)
    

def windstress(date, time, y, ax):
    """
    Wind stress signal over time.
    """
    ax.plot(date, y, c='navy', lw=1.2)

    ax.set_ylabel('Wind stress (N/m²)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(date[0], date[-1])

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)

    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.grid(True, which="both", color='gray', ls=":", lw=0.4)
    

def schmidt(date, time, y, ax):
    """
    Schmidt stability (J/m²).
    """
    ax.plot(date, y, c='navy', lw=1.2)

    ax.set_ylabel('Schmidt stability (J/m²)')
    ax.set_ylim(bottom=0)
    ax.set_xlim(date[0], date[-1])

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)

    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.grid(True, which="both", color='gray', ls=":", lw=0.4)
    ax.set_xlabel('Time')


def wedderburn(date, t, y, ln, ax):
    """
    Wedderburn number at thermocline depth.
    """
    ax.plot(date, y, c='navy', lw=1.2, label='Wedderburn number')

    ax.set_yscale('log')
    ax.set_ylabel('Wedderburn number',color='navy')
    ax.tick_params(axis='y', labelcolor='navy')
    ax.set_xlim(date[0], date[-1])

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel('Time')
    
    ax2 = ax.twinx()
    ax2.set_yscale('log')
    ax2.plot(date, ln, c='red', lw=1.2, label='Lake number')
    ax2.set_ylabel('Lake number',color='red', fontsize=9)
    ax2.tick_params(axis='y', labelcolor='red')
    
    ax.grid(True, which="both", color='gray', ls=":", lw=0.4)

def wedd_limit(dx, time, ri, dw, up, ax):
    """
    Wedderburn number at thermocline (internal wave regime)
    """
    ax.plot(dx, ri, color='navy', lw=1)
    ax.plot(dx, dw, color='red', lw=1, ls='-')
    ax.plot(dx, up, color='red', lw=1, ls='-')

    ax.fill_between(dx, dw, up, where=up >= dw, color='red', alpha=0.3, label='Internal seiche dominance')

    ax.set_yscale('log')
    ax.set_ylabel('Wedderbun number')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.grid(True, which='both', color='black', ls=':', lw=0.3)
    ax.legend(loc='upper right', fontsize=9)

    ax.set_xlim([dx[0], dx[-1]])

def wedd_compb(dx, time, ri, rid, dw, up, ax):
    """
    Plot filtered Richardson number (W) considering aligned/consecutive wind events.
    Uses same graphic style as updated Interwave plots.
    """

    # Main curves 
    ax.plot(dx, ri,  color='navy',       lw=1, label='Standard')
    ax.plot(dx, rid, color='dodgerblue', lw=1, label='Wind direction filter')

    # Threshold curves
    ax.plot(dx, dw, color='red', lw=1)
    ax.plot(dx, up, color='red', lw=1)

    # Filled region (wind-event detection zone)
    ax.fill_between(dx, dw, up, where=up >= dw,
                    facecolor='red', alpha=0.35)

    ax.set_xlim(dx[0], dx[-1])

    # Avoid issues in log scale (must be positive)
    ymin = np.nanmin([np.nanmin(ri), np.nanmin(rid), 1e-6])
    ymax = np.nanmax([np.nanmax(ri), np.nanmax(rid)])
    ax.set_ylim(max(ymin, 1e-6), ymax)

    ax.set_yscale('log')

    ax.set_ylabel('Wedderburn number (-)')
    ax.grid(True, which='both', color='gray', ls=":", lw=0.4)

    # Time axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel('Time')

    ax.legend(loc='upper right', fontsize=9, frameon=True,
              facecolor='white', edgecolor='gray')


def richardson2d(dx, time, h, ri, thermo, cota, z0, ax, fig=None, cbar_n_ticks=4):
    """
    2D Richardson number (filtered)
    """
    ri_valid = np.array(ri, dtype=float)
    ri_valid[~np.isfinite(ri_valid)] = np.nan
    ri_valid[ri_valid <= 0] = np.nan

    if np.all(np.isnan(ri_valid)):
        ax.text(0.5, 0.5, 'No valid Ri data', ha='center', va='center', 
                transform=ax.transAxes)
        return

    vmin = np.nanmin(ri_valid)
    vmax = np.nanmax(ri_valid)
    if not (np.isfinite(vmin) and np.isfinite(vmax) and vmin > 0):
        ax.text(0.5, 0.5, 'No valid Ri range', ha='center', va='center', 
                transform=ax.transAxes)
        return

    # Log normalization for color scale
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    levels = np.logspace(np.floor(np.log10(vmin)),np.ceil(np.log10(vmax)),100)

    contour = ax.contourf(dx, h, ri_valid.T, levels=levels, cmap='viridis',
                          norm=norm, extend='neither')

    # overlays
    ax.plot(dx, thermo, color='red', lw=1, label='Thermocline')
    ax.plot(dx, cota, color='blue', lw=1, label='Surface')

    ax.set_ylim([z0, np.trunc(np.nanmax(cota)) + 1])
    ax.set_xlim(dx[0], dx[-1])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel('Time')

    # Colorbar (log)
    _draw_colorbar_with_background(ax, contour, 'Richardson number (-)', 
                                   log=True, cbar_n_ticks=cbar_n_ticks)

    if z0 == 0:
        ax.set_ylabel('Height above bottom (m)', color='black')
    else:
        ax.set_ylabel('Altitude (m)',  color='black')
        ax.plot([dx[0],dx[-1]], [z0,z0], lw=1, ls='--', color='black', 
                label='Bottom')

    ax.legend(loc='lower left', fontsize=9, frameon=True, facecolor='white')

def buoyancy2d(dx, dx_mod, time, h, n, thermo, cota, z0, hmod, ax,
               fig=None, cbar_n_ticks=4):
    """
    2D Buoyancy Frequency (Hz)
    """
    contour = ax.contourf(dx, h, np.asarray(n).T, levels=100, cmap='viridis', 
                          extend='neither')
    
    try:
        ax.plot(dx_mod, hmod+z0, color='black', lw=1, 
                label='Thermocline (Decomposition)')
    except:
        pass

    ax.plot(dx, thermo, color='red', lw=1, label='Thermocline')
    ax.plot(dx, cota, color='blue', lw=1, label='Surface')

    ax.set_ylim([z0, np.trunc(np.nanmax(cota)) + 1])
    ax.set_xlim(dx[0], dx[-1])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=25)
    ax.set_xlabel('Time', color='black')

    # Colorbar with flexible white box behind it
    _draw_colorbar_with_background(ax, contour, 'Buoyancy frequency (Hz)', 
                                   log=False,  cbar_n_ticks=cbar_n_ticks)

    if z0 == 0:
        ax.set_ylabel('Height above bottom (m)', color='black')
    else:
        ax.set_ylabel('Altitude (m)',  color='black')
        ax.plot([dx[0],dx[-1]], [z0,z0], lw=1, ls='--', color='black', 
                label='Bottom')

    ax.tick_params(colors='black')
    ax.legend(loc='lower left', fontsize=9, frameon=True, facecolor='white')



def tempstructure2d(dx, dx_dec, time, h, temp, thermo, cota, z0, hemod, hvmod,
                    ax, mode='two-layer' ,fig=None, cbar_n_ticks=5):
    """
    2D Temperature structure (°C)
    """

    matlab_jet = cm.get_cmap('jet')

    temp2d = np.asarray(temp).T   

    contour = ax.contourf(dx, h, temp2d, levels=120, cmap=matlab_jet,
        extend='neither')

    # Limits 
    ax.set_ylim([z0, np.trunc(np.nanmax(cota)) + 1])
    ax.set_xlim(dx[0], dx[-1])

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=20)
    ax.set_xlabel('Time')

    # Colorbar (same style as buoyancy2d)
    _draw_colorbar_with_background(ax, contour, label='Temperature (°C)',
        log=False, cbar_n_ticks=cbar_n_ticks)

    if z0 == 0:
        ax.set_ylabel('Height above bottom (m)', color='black')
    else:
        ax.set_ylabel('Altitude (m)',  color='black')
        ax.plot([dx[0],dx[-1]], [z0,z0], lw=1, ls='--', color='black', 
                label='Bottom')
         
    hemod = np.array([ np.asarray(x).flatten()[0] for x in hemod ])
    h1mod = np.array(hvmod)[:,0]
    h2mod = np.array(hvmod)[:,1]
    
    ax.plot(dx, cota,   linewidth=1, color='blue', label='Surface')
    
    if mode == 'two-layer':
        ax.plot(dx_dec, h1mod+z0, linewidth=1, ls='-', color='gray')
        ax.plot(dx_dec, h2mod+z0, linewidth=1, ls='-', color='gray', 
                label='Metalimnion boundary')
    else:
        ax.plot(dx, thermo, linewidth=1, color='red', label='Thermocline')
    
    try:
        ax.plot(dx_dec, hemod+z0, linewidth=1, color='black', 
                label='Thermocline (Decomposition)')
    except:
        pass
    
    ax.legend(loc='lower left', fontsize=9, frameon=True, facecolor='white')

def tempstructure_paper(dx, iso1,iso2, time,h,thermo,temp,cota,z0,ax):                                         # retirar (depois)
    """
    Plot temperature structure time series.
    Generates a graph of the temporal evolution of the temperature profile.
    """  
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
    """
    Plot temperature structure near the thermocline.
    Generates a time-series graph highlighting the upper layer.
    """      
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


def thorpe_displacement2d(dx, h, thorpe, cota, z0, ax, cbar_n_ticks=4):
    """
    Plot Thorpe discplacement in a 2D pplot.
    """
    matlab_jet = cm.get_cmap('jet')

    thorpe2d = np.asarray(thorpe).T  # shape: (z, t)

    # Extract depth vector (1D)
    h = np.asarray(h)
    if h.ndim == 2:
        hvec = h[0, :]
    else:
        hvec = h

    contour = ax.contourf(dx, hvec, thorpe2d, levels=120, cmap=matlab_jet,
        extend='both')

    ax.plot(dx, cota, linewidth=1, color='blue')

    ax.set_ylim([z0, np.trunc(np.nanmax(cota)) + 1])
    ax.set_xlim(dx[0], dx[-1])

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.tick_params(axis='x', rotation=20)
    ax.set_xlabel('Time')

    _draw_colorbar_with_background(ax, contour, label='Thorpe displacement (m)',
        log=False, cbar_n_ticks=cbar_n_ticks)

    if z0 == 0:
        ax.set_ylabel('Height above bottom (m)', color='black')
    else:
        ax.set_ylabel('Altitude (m)',  color='black')
        ax.plot([dx[0],dx[-1]], [z0,z0], lw=1, ls='--', color='black', 
                label='Bottom')

    ax.legend(['Surface'], loc='lower left', fontsize=9,
              frameon=True, facecolor='white')


def thorpe_scale(dx, thorpe, ax):
    """
    1D Thorpe Scale time series
    """

    ax.plot(dx, thorpe, color='navy', lw=1.2)

    ax.set_ylabel('Thorpe scale (m)')
    ax.set_xlim(dx[0], dx[-1])

    ax.grid(True, which="both", color='gray', ls=":", lw=0.4)
 

def movingaverage(interval, window_size):
    """
    Compute moving average using convolution.
    """     
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')
   

def generation(cla, ax):
    """
    Lake mixing classification bar plot.
    """
    cla = np.asarray(cla, dtype=float)
    N = len(cla)
    x = np.arange(N)
    width = 0.6

    ax.bar(x, cla, width, color="navy")
    ax.set_ylim([0, 100])

 
    labels = ["Mixed", "KH", "IS/d", "IS/s"]

    ax.set_xticks(x)
    ax.set_xticklabels(labels[:N], rotation=0)

    ax.grid(axis='y', linestyle=':', linewidth=0.4)

def sensibility_metalimnion(delta,period,m,ax):
    """
    Plot sensitivity analysis for layer thickness variation.
    Shows BSIW period variation as a function of the metalimnion threshold (V2 mode).
    """
    ax.plot(delta, period, linewidth=1, c='black', ls='-')      
    
    ax.grid(True,which="both",color='black',ls=":",lw=0.25)
    
    ax.set_xlabel('metalimnion threshold (kg/m²)')  
    ax.set_ylabel('period (hours)')
    
    ax.set_title(r'$ \delta = $ '+str(m)+' kg/m²',loc='right', fontsize= 9)


def depth_sensitivity(xh, per, label, show_xlabel, ax):
    """
    Sensitivity analysis for layer thickness variation.
    Plots the variation of the BSIW period with layer thickness.
    """

    ax.plot(xh, per, lw=1.4, color='navy')

    # Center text annotation
    xmid = 0.5 * (np.nanmin(xh) + np.nanmax(xh))
    ymid = 0.5 * (np.nanmin(per) + np.nanmax(per))
    ax.text(xmid, ymid, label, fontsize=12, ha='center', va='center')

    # Axes formatting
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax.grid(True, linestyle=':', linewidth=0.4, color='gray')

    if show_xlabel:
        ax.set_xlabel(r'$h^* \, (m)$')
        
        
def densi_sensitivity(xrho, per, label, show_xlabel, ax):
    """
    Sensitivity analysis for density variation.
    """

    ax.plot(xrho, per, lw=1.4, color='darkred')

    xmid = 0.5 * (np.nanmin(xrho) + np.nanmax(xrho))
    ymid = 0.5 * (np.nanmin(per) + np.nanmax(per))
    ax.text(xmid, ymid, label, fontsize=12, ha='center', va='center')

    ax.set_ylabel('Period (h)')
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))  # 3–4 ticks only
    ax.grid(True, linestyle=':', linewidth=0.4, color='gray')

    if show_xlabel:
        ax.set_xlabel(r'$\rho^* \, (kg\, m^{-3})$')
    
    
def parabar_sensitivity(x, per, typ, limin1, lim2, limax1, ax):
    """
    Plot sensitivity of the BSIW period to a general parameter.
    """

    # Convert limits from seconds to hours
    lim1_min = limin1 / 3600
    lim1_max = limax1 / 3600

    # Main sensitivity curve
    ax.plot(x, per, linewidth=1, color="navy")

    # Rectangle region 1 (green)
    rect1 = mpatches.Rectangle(
        (x.min(), lim1_min),
        width=x.max() - x.min(),
        height=lim1_max - lim1_min,
        color="green",
        alpha=0.35
    )
    ax.add_patch(rect1)

    # Rectangle region 2 (gray): confidence interval
    z = 2.58  # 99% CI
    lim2 = np.asarray(lim2, float)
    s = np.nanstd(lim2)
    n = np.count_nonzero(~np.isnan(lim2))
    cl = z * s / np.sqrt(n)

    lim2_center = np.nanmean(lim2)
    lim2_min = lim2_center - cl
    lim2_width = 2 * cl

    rect2 = mpatches.Rectangle(
        (lim2_min, per.min()),
        width=lim2_width,
        height=per.max() - per.min(),
        color="gray",
        alpha=0.4
    )
    ax.add_patch(rect2)

    ax.set_xlim(
        min(x.min(), lim2_min),
        max(x.max(), lim2_min + lim2_width)
    )
    ax.set_ylim(per.min(), per.max())

    if typ == "dep":
        ax.set_ylabel("Wave period (h)")
        ax.set_xlabel(r"$\bar{H}^{-1/2}$")
    elif typ == "rho":
        ax.set_xlabel(r"$ (\rho_h/\Delta\rho)^{1/2}$")

    # Reduce ticks on y-axis for visibility
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

def classification_genera(W1, W2, heh_found, ahe, tau, ax):
    """
    Basin-scale internal wave amplitude (theoretical result)
    Graph to plot the amplitude of BSIW according to multiple models.
    """

    label = [
        'Spigel & Imberger (1980)',
        r'de Carvalho Bueno et al. (2020) $h_e/H=0.7$',
        r'de Carvalho Bueno et al. (2020) $h_e/H=0.5$',
        r'de Carvalho Bueno et al. (2020) $h_e/H=0.2$',
        rf'$W_{{V1H1}}$ {tau[0]}°C' if len(tau) > 0 else r'$W_{V1H1}$',
        rf'$W_{{wind}}$ {tau[0]}°C' if len(tau) > 0 else r'$W_{wind}$'
    ]

    # Domain setup
    wedderburn = np.linspace(0.1, 5, 1000)
    wedd2 = np.linspace(5, 1000, 1000)
    wedd = np.concatenate((wedderburn, wedd2))

    # Model curves
    linear = 0.5 / wedd
    curve1 = mod.bueno(10.5, 4.5, wedderburn, wedderburn) / 10.5
    curve2 = mod.bueno(7.5, 7.5, wedderburn, wedderburn) / 7.5
    curve3 = mod.bueno(3, 12, wedderburn, wedderburn) / 3

    # Plot models
    ax.plot(wedd, linear, '--', color='gray', linewidth=1)
    ax.plot(wedderburn, curve1, '-', color='black', linewidth=1)
    ax.plot(wedderburn, curve2, '-', color='dimgray', linewidth=1)
    ax.plot(wedderburn, curve3, '-', color='silver', linewidth=1)

    # Observations
    ax.scatter(W1, ahe, s=25, marker='o', color='red', label='Observed $W_{V1H1}$')
    ax.scatter(W2, ahe, s=25, marker='x', color='navy', label='Observed $W_{wind}$')

    ax.legend(labels=label, loc='upper right', prop={'size': 9})
    ax.grid(True, color='black', ls=":", lw=0.25)

    ax.set_xscale('log')
    ax.set_xlim(0.1, 1000)
    ax.set_ylim(0, 2)
    ax.set_xlabel(r'$Ri\, h_e/L$ (-)', color='black')
    ax.set_ylabel(r'$\zeta_o/h_e$ (-)', color='black')


def bueno_parameterization(W1, W2, norma, tau, ax):
    """
    de Carvalho Bueno (2020) parameterization
    Predicts the BSIW amplitude from Ri he/L using de Carvalho Bueno al. 2020.
    """

    label = [
        'de Carvalho Bueno et al. (2020)',
        rf'$W_{{V1H1}}$ {tau[0]}°C' if len(tau) > 0 else r'$W_{V1H1}$',
        rf'$W_{{wind}}$ {tau[0]}°C' if len(tau) > 0 else r'$W_{wind}$'
    ]

    wedderburn = np.linspace(0.1, 5, 1000)
    norm = (wedderburn - 6) ** 2

    # Model + points
    ax.plot(wedderburn, norm, color='black', ls=':', lw=1)
    ax.scatter(W1, norma, s=25, marker='o', color='red')
    ax.scatter(W2, norma, s=25, marker='x', color='navy')

    # Labels and grid
    ax.legend(labels=label, loc='lower left', prop={'size': 9})
    ax.grid(True, color='black', ls=":", lw=0.25)

    # Axis config
    ax.set_xscale('log')
    ax.set_xlim(0.1, 10**3)
    ax.set_xlabel(r'$Ri\, h_e/L$ (-)', color='black')
    ax.set_ylabel(r'$2 f^2 \log(\zeta_o / (0.1 h_e))$ (-)', color='black')