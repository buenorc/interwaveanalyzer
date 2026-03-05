# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer - Backend Main Module

Interwave Analyzer - Version 2 (2026) 
Backend Main module version: 2.260305

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

# --------------- Packages and modules ----------------------------------------
import os
import sys
import time
import random
import datetime
import platform
import subprocess
import numpy as np


import matplotlib
matplotlib.use('Agg')

from tkinter import Tk, Text, END
from matplotlib import gridspec
from matplotlib import pyplot as plt


# Interwave Analyzer packages
import iwplot as graph  # package of graphs
import iwmod  as mod    # package of functions
import iwmath as miw    # package of internal wave models
import iwload as load   # package of inputs 
import iwarn  as war    # package of warnings
import iwpars as pars   # package of addition parameters


# --------------------- Tk Text widget ----------------------------------------
class StdoutRedirector(object):
    def __init__(self, text_area):
        self.text_area = text_area

    def write(self, s):
        self.text_area.insert(END, s)
        self.text_area.see(END)

    def flush(self):
        pass

def main():
    print('Backed Main Module is running...')

    start_time = time.time()
    old_stdout = sys.stdout

    root = Tk()
    root.configure(background='white')
    root.title("Interwave Analyzer Running")
    root.geometry("550x530")

    if platform.system() == 'Windows':
        try:
            root.iconbitmap("iwcon.ico")
        except Exception:
            pass

    outputPanel = Text(root, wrap='word', height=30, width=100)
    outputPanel.grid(column=0, row=0, columnspan=2, sticky='NSWE', padx=5, pady=5)

    sys.stdout = StdoutRedirector(outputPanel)

    print("> Interwave Analyzer is starting the data processing... ")
    print("--------------------------------------------------------------------------------------")
    print("> Interwave Analyzer, version 2.00.0        January   2026")
    print("> ")
    print("--------------------------------------------------------------------------------------")
    print("> ")
    print("> Part I       Reading information from GUI... ")

    root.update()


    # Read temporary.txt file 
    with open('temporary.txt', 'r', encoding='utf-8') as fh:
        raw_lines = [ln.strip() for ln in fh.readlines()]

    it = iter(raw_lines)

    # File-paths and GUI basic values
    nam      = next(it)
    win      = next(it)
    sen_nam  = next(it)

    rw       = float(next(it))
    rad      = int(next(it))
    z0       = float(next(it))
    linang   = float(next(it))
    lat      = float(next(it))
    depi     = int(next(it))

    # Basin length/bathymetry type
    type_length = int(next(it))
    if type_length == 1:
        len_basin = float(next(it))
        fna       = None
    elif type_length == 2:
        fna       = next(it)
        len_basin = None
    else:
        len_basin = None
        fna       = None

    # Water level type
    ean_serie = int(next(it))
    if ean_serie == 1:
        ean_cota = float(next(it))
        nac      = -999
    elif ean_serie == 2:
        nac      = next(it)
        ean_cota = -999
    else:
        ean_cota = -999
        nac      = -999

    minval = float(next(it))

    # Filter
    filter_process = int(next(it))
    if filter_process == 2:
        low_per  = float(next(it))
        high_per = float(next(it))
    else:
        _ = next(it)
        _ = next(it)
        low_per, high_per = None, None

    # Window size
    autoff_wsuser = int(next(it))
    windsizeuser  = float(next(it)) * 24 * 60 * 60



    # Decomposition
    dt_decom = float(next(it))
    dt_decom = dt_decom / 60.0  # convert minutes to hours

    # Spectral analysis
    window = next(it)
    mother = next(it)
    window = window.strip().lower()
    mother = mother.strip()
    

    # Isotherms
    turn_iso = int(next(it))
    tau = np.zeros(4, float)
    tan = np.zeros(4, float)
    for i in range(4):
        numaux = float(next(it))
        if turn_iso == 0:
            tau[i] = -999
            tan[i] = 0
        else:
            tau[i] = numaux
            tan[i] = 1

    # Comparisons
    c1 = next(it)
    c2 = next(it)
    
    ana1 = mod.comparison_definition(c1)
    ana2 = mod.comparison_definition(c2)

    # Sensors
    turn_temp = int(next(it))
    sen = np.zeros(4, float)
    seu = np.zeros(4, int)
    for i in range(4):
        seu[i] = int(next(it))
        if turn_temp == 0 or seu[i] == -999:
            seu[i] = -999
            sen[i] = 0
        else:
            sen[i] = 1

    # Smoothing and output folder
    largelen    = int(next(it))
    output_path = next(it)
    if not output_path.endswith('/'):
        output_path += '/'

    # Shapefile (not used yet!)
    shape_checked = int(next(it))
    if shape_checked:
        shape_on = 1
        path_shape = next(it)
    else:
        shape_on = 0
        path_shape = None

    # Additional Parameters
    additional_params = {}
    for line in it:
        line = line.strip()
        if not line or line == '-999':
            continue
    
        if '=' in line:
            key, value = line.split('=', 1)
            additional_params[key.strip()] = value.strip()

    # Remove loading information file
    try:
        os.remove('temporary.txt')
    except Exception:
        pass

    root.update()

    # ---------------- Load data from file ------------------------------------
    os.makedirs(output_path, exist_ok=True)
    diag_fpath = os.path.join(output_path, 'diagnosis.txt')
    dig = open(diag_fpath, 'w', encoding='utf-8')
    dig.write("-------------------------------------------------------------------------------------\n")
    dig.write("Interwave Analyzer diagnosis\n")
    dig.write("-------------------------------------------------------------------------------------\n\n")

    print("> Part II      Loading information from files... ")
    root.update()

    lin, qt = load.file_len(nam, col='on')

    date, temp, serialt, dt = load.temper_read(nam, lin, qt)
    ean, h, tempa = load.serial_cota(serialt, nac, lin, qt, temp, sen_nam, ean_serie, ean_cota, z0)
    wind, dw, ra = load.serial_wind(serialt, win, rad, lin)
    
    seu = seu - 1
    
    type_length, fna, fna_trans, change_basin, warnTrans, missFile = \
        load.lengthType(type_length, fna, additional_params)

    if warnTrans:
        war.bathyMissing(dig, missFile)
    
    # Load geometry
    auxDepth = ean - z0
    longData, transData = load.loadData(type_length, fna, fna_trans, len_basin, change_basin, np.mean(auxDepth))
    
    # Apply orientation 
    if type_length == 3:
        longData = mod.basinOrientation(longData, change_basin)
        transData = mod.basinOrientation(transData, (change_basin + 90) % 360)
    
    elif type_length == 2:
        longData = mod.basinOrientation(longData, 270)
    
    elif type_length == 1:
        longData =mod.basinOrientation(longData, 270)
    

    print("> Part III     Defining constants and creating variables... ")
    root.update()
    
    # ------------- Organizing Additional Parameters --------------------------
    
    nameBasin = pars.extractName(additional_params)

    # ---- Defining parameters and constants + creating diagnosis file --------

    aux_iso = 1  # Iso number used to apply the corehence analysis with wind

    dj      = 0.1
    dci     = 1
    ddep    = 0.5
    drho    = 0.1
    lh_sen  = 100
    lp_sen  = 20

    if dt_decom == -999:
        dt_decom = dt
        war.decomp_default(dig, dt * 60)
    else:
        if dt_decom < dt:
            dt_decom = dt
            war.decomp_changed(dig, dt * 60)
        else:
            if ((dt_decom * 60) % (dt * 60) == 0):
                war.decomp_specified(dig, dt_decom * 60)
            else:
                dt_decom = int(dt*60*round(float(dt_decom * 60) / (dt * 60)))
                if dt_decom >= dt:
                    war.decomp_multiple(dig, dt_decom * 60)
                else:
                    dt_decom = dt
                    war.decomp_changed(dig, dt * 60)

    Ndeco = int(dt_decom / dt)

    # General parameters
    ls_fetch = np.empty(lin, float)
    
    glin = np.empty(lin, float)
    wast = np.empty(lin, float)
    strs = np.empty(lin, float)

    n          = np.empty(lin, float)
    schmidt    = np.empty(lin, float) 
    lakeNumber = np.empty(lin, float) 
    riw        = np.empty(lin, float)
    wedd       = np.empty(lin, float)
    wedd_inv   = np.empty(lin, float)

    thermo_temp = np.empty(lin, float)

    n_slope      = np.empty(lin, float)
    wave_slope   = np.empty(lin, float)
    longLeft     = np.empty(lin, float)
    longRight    = np.empty(lin, float)
    atLeft       = np.empty(lin, float)
    atRight      = np.empty(lin, float)
    transLeft    = np.empty(lin, float)
    transRight   = np.empty(lin, float)
    atcrossLeft  = np.empty(lin, float)
    atcrossRight = np.empty(lin, float)
    
    ht = np.empty(lin, float)
    he = np.empty(lin, float)
    hh = np.empty(lin, float)
    pe = np.empty(lin, float)
    ph = np.empty(lin, float)
    
    cond1 = np.empty(lin, float)
    cond2 = np.empty(lin, float)
    cond3 = np.empty(lin, float)

    Pbar = np.empty(lin, float)
    Hbar = np.empty(lin, float)

    v1mode = np.empty(lin, float)
    v2mode = np.empty(lin, float)

    h1 = np.empty(lin, float)
    h2 = np.empty(lin, float)
    h3 = np.empty(lin, float)
    p1 = np.empty(lin, float)
    p2 = np.empty(lin, float)
    p3 = np.empty(lin, float)

    pu = np.empty(lin, float)
    pd = np.empty(lin, float)
    hH = np.empty(lin, float)

    isoa = np.empty(lin, float)
    isob = np.empty(lin, float)
    isoc = np.empty(lin, float)
    isod = np.empty(lin, float)

    genera = np.empty(lin, float)
    iw_up = np.empty(lin, float)
    iw_dw = np.empty(lin, float)

    # Two-dim parameters: shape (lin, qt-1)
    riw2d = np.empty((lin, qt - 1), float)
    hzmid = np.empty((lin, qt - 1), float)

    dw_hom = np.empty(lin, float)
    dw_lit = np.empty(lin, float)
    dw_spi = np.empty(lin, float)

    wedu = np.empty(lin, float)
    wedi = np.empty(lin, float)

    # Initialize values that were previously assigned None in loop
    isoa.fill(np.nan)
    isob.fill(np.nan)
    isoc.fill(np.nan)
    isod.fill(np.nan)

    dw_hom.fill(np.nan)
    dw_lit.fill(np.nan)
    dw_spi.fill(np.nan)


    # ---------------------- Computing parameters -----------------------------
    print("> ")
    
    load_time = time.time()
    
    print("> Execution time for part I, II, and III: ")
    print("> " + str(round(load_time - start_time, 4)) + ' seconds')
    
    root.update()

    iw       = mod.velocityten(wind, rw)
    hmean    = np.mean(h, axis=0)

    print("--------------------------------------------------------------------------------------")
    print('> ')
    print("> Part IV      Starting computation... ")
    print('> ')
    root.update()

    # Temporal variables and counters
    consecutive_dire = 0
    consecutive_wind = 0
    ver_dire         = 0
    ver_wind         = 0

    error_thermo     = 0
    error_3layer     = 0
    imod             = 0
    period_time      = []
    vel_time         = []
    cpzin            = []
    time_model       = []
    date_model       = []
    h_deco           = []
    buoy             = []

    auxisa = auxisb   = auxisc   = auxisd    = None
    war3   = warmode1 = warmode2 = war3from2 = 0

    nanmax = np.nanmax


    # Progress reporting parameters
    progress_every = max(1, int(lin / 100))  # update roughly every 1%
    dt_sec = dt * 3600.0                     # seconds per dt
    
    he_t_lasts = {}
    
    hemod           = []
    mode2_nodes     = [] 
    mode3_nodes     = [] 
    mode4_nodes     = []
    mode1Layer      = []
    mode2Layer      = []
    mode3Layer      = []
    mode4Layer      = []
    refinedDepthArr = []
    
    
    # Define wind fetch for cases of 1D or 2D cases
    if type_length in [1, 2]:

        L_long = longData["dists"][0]
        ls_fetch = np.full(lin, L_long)
        
    longData = mod.computeBasinSlopes(longData)
    
    if transData is not None:
        transData = mod.computeBasinSlopes(transData)
    
    for t in range(lin):
        
        # Define wind fetch for cases of 3D case (transversal)
        if type_length == 3:
            ls_fetch[t] = mod.fetchVariable(longData, transData, dw[t])

        # Local references to arrays and row slices
        auxtem = tempa[t, :]
        auxh = h[t, :]
        auxiw = iw[t]
        auxean = ean[t]
        max_ean = nanmax(auxean)
        min_ean = z0

        auxtem_ordered = mod.sorting_1d(auxtem)
        rho_up    = mod.commission(auxtem_ordered[0])
        rho_bot   = mod.commission(auxtem_ordered[-1])
        
        
        Htotal     = auxDepth[t]

        # Decomposition: run when imod==0 or 1 
        if imod == 0 or imod == 1:
            auxh_deco = auxh - z0
            
            vel_aux, cpzin_aux, refined_depth, period_aux, cond = \
                miw.decomposition(ls_fetch[t], auxtem_ordered, \
                                  Htotal - auxh_deco, Htotal)
            
            if cond == 1:
                war.decomposition_overflow(dig)
            refined_depth = Htotal - refined_depth
            
            cross, crossThicker = mod.findLayers(vel_aux,refined_depth,Htotal)

            hemod.append(cross[0])        # Mode 0 - 1 zero crossing
            mode2_nodes.append(cross[1])
            mode3_nodes.append(cross[2])
            mode4_nodes.append(cross[3])
            
            mode1Layer.append(crossThicker[0])
            mode2Layer.append(crossThicker[1])
            mode3Layer.append(crossThicker[2])
            mode4Layer.append(crossThicker[3])

            
            refinedDepthArr.append(refined_depth)

            h_deco.append(np.concatenate((auxh_deco, [0])))
            date_model.append(date[t])
            time_model.append(t * dt)
            period_time.append(period_aux)
            vel_time.append(vel_aux)
            cpzin.append(cpzin_aux)
            imod = Ndeco
        else:
            imod = imod - 1

        n_slope[t] = np.sqrt(9.81*abs(rho_bot-rho_up)/(rho_up*Htotal))
        wave_slope[t] = mod.waveSlope(n_slope[t], period_aux[0])
        
        slopes_t = mod.thermoclineSlopes(type_length, longData, transData, auxean, cross[0][0], z0)

        longLeft[t]  = slopes_t["long"]["left"]
        longRight[t] = slopes_t["long"]["right"]

        atLeft[t]  = longLeft[t]/wave_slope[t]
        atRight[t] = longRight[t]/wave_slope[t]

        if type_length == 3:
            transLeft[t]  = slopes_t["trans"]["left"]
            transRight[t] = slopes_t["trans"]["right"]
        
            atcrossLeft[t]  = transLeft[t]/wave_slope[t]
            atcrossRight[t] = transRight[t]/wave_slope[t] 

        # Compute isotherms only if requested
        if (turn_iso == 1):
            if (tau[0] != -999):
                isoa[t] = mod.isotherms(tau[0], qt, auxh, auxtem_ordered, max_ean, min_ean, auxisa)
                auxisa = isoa[t]
            if (tau[1] != -999):
                isob[t] = mod.isotherms(tau[1], qt, auxh, auxtem_ordered, max_ean, min_ean, auxisb)
                auxisb = isob[t]
            if (tau[2] != -999):
                isoc[t] = mod.isotherms(tau[2], qt, auxh, auxtem_ordered, max_ean, min_ean, auxisc)
                auxisc = isoc[t]
            if (tau[3] != -999):
                isod[t] = mod.isotherms(tau[3], qt, auxh, auxtem_ordered, max_ean, min_ean, auxisd)
                auxisd = isod[t]

        # Two-layer structure
        he_t, hh_t, pe_t, ph_t, glin_t, n_t, pu_t, pd_t, error, he_t_lasts = mod.structure2layer(
        qt, auxh, auxtem_ordered, auxean, z0, he_t_lasts)
        
        error_thermo = error_thermo + error 

        ht[t] = ean[t] - he_t
        he[t] = he_t
        hh[t] = hh_t
        hH[t] = he_t / (he_t + hh_t) if (he_t + hh_t) != 0 else 0.0

        Hbar[t] = np.sqrt((he_t + hh_t) / (he_t * hh_t)) if (he_t * hh_t) != 0 else np.nan

        try:
            Pbar[t] = np.sqrt(ph_t / (ph_t - pe_t))
        except Exception:
            Pbar[t] = np.nan

        if hH[t] > 0.5:
            hH[t] = 1.0 - hH[t]

        strs[t], wast[t], riw[t] = mod.wind_parameters(auxiw, rw, pe_t, he_t, n_t, glin_t, auxean)

        # Store two-layer results
        pe[t] = pe_t
        ph[t] = ph_t
        glin[t] = glin_t
        n[t] = n_t
        pu[t] = pu_t
        pd[t] = pd_t

        # 2d parameters
        p, n2, hmid, gli2d = mod.thermal_stability(qt, auxh, auxean, auxtem)
        riw2 = mod.richardson(auxiw, rw, qt, auxh, pe_t, auxean, n2, hmid, p, gli2d)
        
        buoy.append(n2)
        

        riw2d[t, :] = riw2[:]
        hzmid[t, :] = hmid[:]

        # Wedderburn number
        wedd[t] = mod.wedderburn(glin_t, he_t, wast[t], ls_fetch[t])
        wedd_inv[t] = 1.0 / wedd[t] if wedd[t] != 0 else np.nan

        # Wind parametrization
        try:
            tbsiw = 2.0 * ls_fetch[t] / np.sqrt(glin_t * he_t * hh_t / (he_t + hh_t))
        except Exception:
            tbsiw = np.nan
            war.merian(dig)

        dw_hom[t] = np.nan

        max_spigel = ls_fetch[t] * (he_t + hh_t) / (4.0 * he_t * hh_t) if (he_t * hh_t) != 0 else np.nan
        min_spigel = 0.5 * np.sqrt((he_t + hh_t) / (he_t * hh_t)) if (he_t * hh_t) != 0 else np.nan

        wedd_aux = 1.0
        if (not np.isnan(min_spigel)) and (min_spigel < 1):
            wedd_aux = min_spigel

        # Compute winsize safely
        if np.isnan(tbsiw):
            if len(period_time) > 0:
                aux_tbsiw = np.array(period_time)
                
                try:
                    aux_tbsiw_val = aux_tbsiw[-1][0]
                except Exception:
                    aux_tbsiw_val = aux_tbsiw[-1]
            else:
                aux_tbsiw_val = (5 * 24 * 3600)  # fallback
            winsiz = int(0.25 / 2 * aux_tbsiw_val / dt_sec)
        else:
            winsiz = int(0.25 / 2 * tbsiw / dt_sec)

        points = t - winsiz
        polong = t - int(2 * winsiz)
        if points < 1:
            points = 0
        if polong < 1:
            polong = 0

        try:
            wedd_mean = np.nanmean(wedd[max(0, points):t + 1])
        except Exception:
            wedd_mean = wedd[t]

        if (wedd_mean < 20 and wedd_mean > wedd_aux):
            try:
                dw_lit[t] = mod.wind_average(dw[points:t + 1], iw[points:t + 1])
            except Exception:
                dw_lit[t] = dw[t]
        else:
            dw_lit[t] = np.nan

        if (wedd_mean < max_spigel) and (wedd_mean > min_spigel):
            try:
                dw_spi[t] = mod.wind_average(dw[points:t + 1], iw[points:t + 1])
            except Exception:
                dw_spi[t] = dw[t]

            if consecutive_wind > 0:
                try:
                    wmin, wmax = mod.wind_angle(np.mean(dw[polong:polong + winsiz + 1]), linang)
                except Exception:
                    wmin, wmax = 0, 360

                try:
                    dw10 = mod.wind_average(dw[points:t + 1], iw[points:t + 1])
                except Exception:
                    dw10 = dw[t]

                if wmin > wmax:
                    if (dw10 >= wmin or dw10 <= wmax):
                        consecutive_dire += 1
                        dw_hom[t] = 357
                    else:
                        if (consecutive_dire > ver_dire):
                            ver_dire = consecutive_dire
                        consecutive_dire = 0
                else:
                    if (dw10 >= wmin and dw10 <= wmax):
                        consecutive_dire += 1
                        dw_hom[t] = 357
                    else:
                        if (consecutive_dire > ver_dire):
                            ver_dire = consecutive_dire
                        consecutive_dire = 0

            consecutive_wind += 1
            if (consecutive_wind > ver_wind):
                ver_wind = consecutive_wind
            if (consecutive_dire > ver_dire):
                ver_dire = consecutive_dire
        else:
            consecutive_dire = 0
            consecutive_wind = 0
            dw_spi[t] = np.nan

        # 3-layer structure (suing strongst gradient methods)
        try:
            h1_t, h2_t, h3_t, p1_t, p2_t, p3_t, er3 = mod.structure3layer(qt, auxh, auxtem, minval, auxean, z0)
            error_3layer += er3
            h1[t], h2[t], h3[t] = h1_t, h2_t, h3_t
            p1[t], p2[t], p3[t] = p1_t, p2_t, p3_t
        except Exception:
            h1[t], h2[t], h3[t], p1[t], p2[t], p3[t] = mod.approx_layer(he_t, hh_t, pe_t, ph_t)
            war3from2 += 1

        if h1[t] == -999:
            war3 += 1
            h1[t], h2[t], h3[t], p1[t], p2[t], p3[t] = mod.approx_layer(he_t, hh_t, pe_t, ph_t)

        # Lake mixing classification
        genera[t],cond1[t], cond2[t], cond3[t] = mod.class_generation(riw[t], hh_t, he_t, ls_fetch[t])
        iw_dw[t], iw_up[t] = mod.iw_generation(wedd[t], hh_t, he_t, ls_fetch[t])

        # Schmidt Stability (J/m2)
        depthSensor = auxean-auxh
        schmidt[t] = mod.schmidtStability(auxtem_ordered,depthSensor , longData, transData)
        lakeNumber[t] = mod.lakeNumber(schmidt[t], wast[t], h1[t], h1[t]+h2[t], ph[t], longData, transData)    
    

        try:
            _, v1mode_t, _ = np.real(miw.disp_zmodel(pe_t, ph_t, he_t, hh_t, ls_fetch[t], 1))
            v1mode[t] = v1mode_t
        except Exception:
            v1mode[t] = np.nan

        try:
            _, v2mode_t, _ = np.real(miw.disp_xmodel3(p1[t], p2[t], p3[t], h1[t], h2[t], h3[t], ls_fetch[t], 2, 1))
            v2mode[t] = v2mode_t
        except Exception:
            v2mode[t] = np.nan

        if t > 0:
            try:
                if not np.isnan(v2mode[t]) and not np.isnan(v2mode[t - 1]) and abs(v2mode[t] - v2mode[t - 1]) > 5 * 60 * 60:
                    v2mode[t] = v2mode[t - 1]
                    warmode2 += 1
            except Exception:
                pass

            try:
                if not np.isnan(v1mode[t]) and not np.isnan(v1mode[t - 1]) and abs(v1mode[t] - v1mode[t - 1]) > 5 * 60 * 60:
                    v1mode[t] = v1mode[t - 1]
                    warmode1 += 1
            except Exception:
                pass

        # Progress update 
        if (t % progress_every) == 0:
            print(f"> Progress: {t}/{lin} ({(t / lin) * 100:.1f}%)")
            root.update()



    # Warnings computed inside of loop 
    
    if error_3layer > 0:
        war.three_layer(dig, 100.0 * error_3layer / lin)
    if warmode1 > 0:
        war.profile_structure(dig, 1, 100.0 * warmode1 / lin)
    if warmode2 > 0:
        war.profile_structure(dig, 2, 100.0 * warmode1 / lin)  # original used warmode1 in denominator
    if war3 > 0:
        war.metalimnion(dig, 100.0 * war3 / lin)
    if war3from2 > 0:
        war.threetotwo(dig, 100.0 * war3from2 / lin)
    if error_thermo > 0:
        war.thermocline(dig, 100.0 * error_thermo / lin)

    iso = [isoa, isob, isoc, isod]

    # Date formatting
    dx_mod = [datetime.datetime.strptime(d, '%Y/%m/%d/%H/%M') for d in date_model]
    dx = [datetime.datetime.strptime(d, '%Y/%m/%d/%H/%M') for d in date]


    # Convert lists to arrays
    time_model = np.array(time_model)
    period_time = np.array(period_time)
    vel_time = np.array(vel_time)
    cpzin = np.array(cpzin)
    buoy = np.array(buoy)
    
    # Interpolating decomp variables to dx size (original time-step)
    if dt_decom > dt:
        war.decomp_interpolation(dig)
        period_time = mod.resampleSignalMatrix(period_time, dx_mod, dx)
        mode2_nodes_interpolated = mod.resampleSignalMatrix(mode2_nodes, dx_mod, dx)
        hemod_inter = mod.resampleSignalMatrix(hemod, dx_mod, dx)
    else:
        mode2_nodes_interpolated = mode2_nodes
        hemod_inter = hemod


    # Internal seiche periods (time variation)
    T11  = mod.period_analysis(period_time[:,0]*3600)
    T21  = mod.period_analysis(period_time[:,1]*3600)
    T31  = mod.period_analysis(period_time[:,2]*3600)
    T41  = mod.period_analysis(period_time[:,3]*3600)
    

    # Splitting analyzed period into three groups
    group = [genera[i:i + int(lin / 3)] for i in range(0, len(genera), int(lin / 3))]
    hH_gp = [hH[i:i + int(lin / 3)] for i in range(0, len(hH), int(lin / 3))]
    wi_gp = [wedd_inv[i:i + int(lin / 3)] for i in range(0, len(wedd_inv), int(lin / 3))]

    h_lim1 = mod.ci(hH_gp[0])
    W_lim1 = mod.ci(wi_gp[0])
    h_lim2 = mod.ci(hH_gp[1])
    W_lim2 = mod.ci(wi_gp[1])
    h_lim3 = mod.ci(hH_gp[2])
    W_lim3 = mod.ci(wi_gp[2])


    dx_gp = [dx[i:i + int(lin / 3)] for i in range(0, len(dx), int(lin / 3))]

    P1 = [dx_gp[0][0], dx_gp[0][-1]]
    P2 = [dx_gp[1][0], dx_gp[1][-1]]
    P3 = [dx_gp[2][0], dx_gp[2][-1]]

    mean_temp, low_temp, high_temp, low_temp_sd, high_temp_sd = mod.stad_deviation(tempa)
    mean_buoy, low_buoy, high_buoy, low_buoy_sd, high_buoy_sd = mod.stad_deviation(buoy)

    mean_h  = np.mean(h, axis=0)
    mean_hm = np.mean(hzmid, axis=0)
    

    m_pe = np.nanmean(pe)
    m_ph = np.nanmean(ph)
    m_he = np.nanmean(he)
    m_hh = np.nanmean(hh)
    m_h1 = np.nanmean(h1)
    m_h2 = np.nanmean(h2)
    m_h3 = np.nanmean(h3)
    m_p1 = np.nanmean(p1)
    m_p2 = np.nanmean(p2)
    m_p3 = np.nanmean(p3)
    m_glin = np.nanmean(glin)
    m_wast = np.nanmean(wast)
    m_ean  = np.nanmean(ean)
    m_ls   = np.nanmean(ls_fetch)
    
    c_he = 1.96 * np.nanstd(he, ddof=1) / np.sqrt(np.sum(~np.isnan(he)))

    c_hh = 1.96 * np.nanstd(hh, ddof=1) / np.sqrt(np.sum(~np.isnan(hh)))
    c_h1 = 1.96 * np.nanstd(h1, ddof=1) / np.sqrt(np.sum(~np.isnan(h1)))
    c_h2 = 1.96 * np.nanstd(h2, ddof=1) / np.sqrt(np.sum(~np.isnan(h2)))
    c_h3 = 1.96 * np.nanstd(h3, ddof=1) / np.sqrt(np.sum(~np.isnan(h3)))

    c_pe = 1.96 * np.nanstd(pe, ddof=1) / np.sqrt(np.sum(~np.isnan(pe)))
    c_ph = 1.96 * np.nanstd(ph, ddof=1) / np.sqrt(np.sum(~np.isnan(ph)))
    c_p1 = 1.96 * np.nanstd(p1, ddof=1) / np.sqrt(np.sum(~np.isnan(p1)))
    c_p2 = 1.96 * np.nanstd(p2, ddof=1) / np.sqrt(np.sum(~np.isnan(p2)))
    c_p3 = 1.96 * np.nanstd(p3, ddof=1) / np.sqrt(np.sum(~np.isnan(p3)))


    m_n = np.nanmean(n) / (2.0 * np.pi)
    m_ht = np.nanmean(ht)

    wedd_lim_lower = m_ls * (m_he + m_hh) / (4.0 * m_he * m_hh)   if (m_he * m_hh) != 0 else np.nan
    wedd_lim_upper = 0.5 * np.sqrt((m_he + m_hh) / (m_he * m_hh)) if (m_he * m_hh) != 0 else np.nan

    try:
        m_dw_spi = np.nanmean(dw_spi)
    except Exception:
        war.homogeneous_condition(dig)
        m_dw_spi = -999

    print("> ")
    print(">         Parameters were defined")
    root.update()

    # Correction of parameters according to wind filtering
    wedu, wedi, fdura, fdire, dura, dire = mod.weddFilters(
    wedd=wedd,
    ver_wind=ver_wind,
    ver_dire=ver_dire,
    dt=dt,
    tbsiw=tbsiw,
    dw_hom=dw_hom,
    v1mode=period_time[:,0]*3600,  # Fundamental wave period (seconds)
    lin=lin
    )


    # Internal seiche modeling (V1H1)
    cp = np.sqrt(m_glin * m_he * (1.0 - m_he / (m_he + m_hh))) if (m_he + m_hh) != 0 else np.nan

    pv1h1 = T11
    pv2h1 = T21

    # coriolis & Bu
    omega_e = 7.2921e-5
    lat_rad = np.radians(lat)
    fc = 2.0 * omega_e * np.sin(lat_rad)
    fo = abs(fc / (2.0 * np.pi))

    if fc != 0:
        Bu = m_n ** 2 * (m_he + m_hh) ** 2 / (fc ** 2 * m_ls ** 2)
    else:
        Bu = 0

    f11, cf11 = 1.0 / T11, 1.0 / miw.coriolis_effect(T11, lat)
    f21, cf21 = 1.0 / T21, 1.0 / miw.coriolis_effect(T21, lat)
    f31, cf31 = 1.0 / T31, 1.0 / miw.coriolis_effect(T31, lat)
    

    xpe, v1h1_spe = miw.sensitivity_2layer(m_pe, drho, lp_sen, m_pe, m_ph, m_he, m_hh, ls_fetch, 1)
    xph, v1h1_sph = miw.sensitivity_2layer(m_ph, drho, lp_sen, m_pe, m_ph, m_he, m_hh, ls_fetch, 2)
    xhe, v1h1_she = miw.sensitivity_2layer(m_he, ddep, lh_sen, m_pe, m_ph, m_he, m_hh, ls_fetch, 3)
    xhh, v1h1_shh = miw.sensitivity_2layer(m_hh, ddep, lh_sen, m_pe, m_ph, m_he, m_hh, ls_fetch, 4)

    xp1, v2h1_sp1 = miw.sensitivity_3layer(m_p1, drho, lp_sen, m_p1, m_p2, m_p3, m_h1, m_h2, m_h3, ls_fetch, typ=1)
    xp2, v2h1_sp2 = miw.sensitivity_3layer(m_p2, drho, lp_sen, m_p1, m_p2, m_p3, m_h1, m_h2, m_h3, ls_fetch, typ=2)
    xh1, v2h1_sh1 = miw.sensitivity_3layer(m_h1, ddep, lp_sen, m_p1, m_p2, m_p3, m_h1, m_h2, m_h3, ls_fetch, typ=3)
    xh2, v2h1_sh2 = miw.sensitivity_3layer(m_h2, drho, lp_sen, m_p1, m_p2, m_p3, m_h1, m_h2, m_h3, ls_fetch, typ=4)

    rho_bar, dep_bar, Prho, Pdep = miw.sensitivity_dimension(m_ls, m_pe, m_ph, m_he, m_hh)

    print(">         Internal wave periods were estimated")
    root.update()

    # Theory of Generation of Basin-scale internal waves
    genera_0 = [np.count_nonzero(group[0] == 1), np.count_nonzero(group[0] == 2), np.count_nonzero(group[0] == 3), np.count_nonzero(group[0] == 4)]
    genera_1 = [np.count_nonzero(group[1] == 1), np.count_nonzero(group[1] == 2), np.count_nonzero(group[1] == 3), np.count_nonzero(group[1] == 4)]
    genera_2 = [np.count_nonzero(group[2] == 1), np.count_nonzero(group[2] == 2), np.count_nonzero(group[2] == 3), np.count_nonzero(group[2] == 4)]

    for i in range(4):
        genera_0[i] = genera_0[i] * 100.0 / len(group[0])
        genera_1[i] = genera_1[i] * 100.0 / len(group[1])
        genera_2[i] = genera_2[i] * 100.0 / len(group[2])

    print(">         Periods were classified according to lake mixing")
    root.update()

    # Spectral analysis of sensors and isotherms
    if int(autoff_wsuser) == 1:
        if pv1h1[1] != None:   
            windsizeuser = 10*pv1h1[1]
        else:
            windsizeuser = 5*24*60*60
            war.welch_average(dig,100*warmode1/lin)  


    # Prepare outputs lists 
    sensor_filtered = []
    timee = []
    per = []
    power_sensor = []
    freqe = []
    welch_sensor = []

    if filter_process == 1:
        low1, high1 = 1.0 / ((pv1h1[2] / 60.0) / 60.0), 1.0 / ((pv1h1[0] / 60.0) / 60.0)
    else:
        low1, high1 = 1.0 / high_per, 1.0 / low_per

    if turn_temp == 1:
        
        # PSD per sensor
        for index in range(qt):
            new = temp[:, index]
            fs = 1.0 / dt
            try:
                filtered_band = mod.butter_bandpass_filter(new, low1, high1, fs)
            except ValueError:
                filtered_band = None
                war.bandpass(dig, 'sensor')

            try:
                aux_time, aux_per, aux_power = mod.wave_spectral(new, dt, mother)
                aux_freq, _, aux_welch, _, _ = mod.welch_method(new, windsizeuser, window, dt)
            except Exception:
                war.spectral(dig, 'sensor')
                aux_time, aux_per, aux_power = 0, 0, 0
                aux_freq, aux_welch, aux_wr, aux_conf = 0, 0, 0, 0

            sensor_filtered.append(filtered_band)
            timee.append(aux_time)
            per.append(aux_per)
            power_sensor.append(aux_power)
            freqe.append(aux_freq)
            welch_sensor.append(aux_welch)


    # Spectral analysis for isotherms (if available)
    if filter_process == 1:
        lowcut, highcut = 1.0 / ((pv1h1[2] / 60.0) / 60.0), 1.0 / ((pv1h1[0] / 60.0) / 60.0)
    else:
        lowcut, highcut = 1.0 / high_per, 1.0 / low_per

    freq = []
    power = []
    time_temp = []
    per_temp = []
    band = []
    welch = []
    wr = []
    conf = []
    iso_corrected = np.zeros(4, float)

    if (turn_iso == 1):
        ean_norm = ean - np.mean(ean)
        for i in range(4):
            if np.isnan(iso[i]).all() and tau[i] != -999:
                war.isotherm_boundary(dig, str(tau[i]))
                tau[i] = -999

            if tau[i] != -999:
                try:
                    aux_time, aux_per, aux_power = mod.wave_spectral(iso[i], dt, mother)
                    aux_freq, _, aux_welch, aux_wr, aux_conf = mod.welch_method(iso[i], windsizeuser, window, dt)
                except Exception:
                    aux_time, aux_per, aux_power = 0, 0, 0
                    aux_freq, aux_welch, aux_wr, aux_conf = 0, 0, 0, 0
                    war.spectral(dig, 'isotherm')

                fs = 1.0 / dt
                try:
                    aux_band = mod.butter_bandpass_filter(iso[i], lowcut, highcut, fs)
                except ValueError:
                    aux_band = None
                    war.bandpass(dig, 'isotherm')

                iso_corrected[i] = np.mean(abs((iso[i] - np.mean(iso[i])) - ean_norm))
                time_temp.append(aux_time)
                per_temp.append(aux_per)
                power.append(aux_power)
                freq.append(aux_freq)
                welch.append(aux_welch)
                band.append(aux_band)
                wr.append(aux_wr)
                conf.append(aux_conf)
            else:
                time_temp.append(0)
                per_temp.append(0)
                power.append(0)
                freq.append(0)
                welch.append(0)
                band.append(0)
                wr.append(0)
                conf.append([0])

        banda = [np.nanmax(np.abs(band[i])) if band[i] is not None and not (isinstance(band[i], int) and band[i] == 0) else 0 for i in range(4)]
        amax = max(np.array(banda))
        aind = tau[np.where(banda == amax)]

    # Spectral analysis of solar radiation (if available)
    if rad == 1:
        solar_ws = 5 * 24 * 60 * 60
        try:
            _, wl_aper_sol, welch_sol, wr_sol, conf_sol = mod.welch_method(ra, solar_ws, window, dt)
        except Exception:
            wl_aper_sol, welch_sol, wr_sol, conf_sol = 0, 0, 0, 0
            war.spectral(dig, 'radiation')

    # Spectral analysis of wind intensity
    try:
        time_win, per_win, power_win = mod.wave_spectral(iw, dt, mother)
        _, wl_aper_win, welch_win, wr_win, conf_win = mod.welch_method(iw, windsizeuser, window, dt)
    except Exception:
        time_win, per_win, power_win = 0, 0, 0
        wl_aper_win, welch_win, wr_win, conf_win = 0, 0, 0, 0
        war.spectral(dig, 'wind')

    # Coherence (isotherms and meteorological data)
    if turn_temp == 1:
        
        # Compute depths only once for each selected sensor index
        d1 = mod.depths(seu[0], h, lin)
        d2 = mod.depths(seu[1], h, lin)
        d3 = mod.depths(seu[2], h, lin)
        d4 = mod.depths(seu[3], h, lin)
        depth = [d1, d2, d3, d4]

        s_filtered = []
        phws = []
        caws = []
        faws = []
        c95ws = []

        for i in range(4):
            if sen[i] == 1:
                new = temp[:, seu[i]]
                a1, a2, a3, a4 = mod.coherence_shift(new, iw, windsizeuser, dt)
                phws.append(a1)
                caws.append(a2)
                faws.append(a3)
                c95ws.append(a4)
                s_filtered.append(sensor_filtered[seu[i]])
            else:
                phws.append(None)
                caws.append(None)
                faws.append(None)
                c95ws.append(None)
                s_filtered.append(None)

        phij = []
        cohij = []
        fij = []
        c95ij = []
        for i in range(4):
            for j in range(4):
                if (i != j and j - i == 1):
                    if (sen[i] == 1 and sen[j] == 1):
                        ni = temp[:, seu[i]]
                        nj = temp[:, seu[j]]
                        ph_aux, coh, f, c95 = mod.coherence_shift(ni, nj, windsizeuser, dt)
                        phij.append(ph_aux)
                        cohij.append(coh)
                        fij.append(f)
                        c95ij.append(c95)
                    else:
                        phij.append(None)
                        cohij.append(None)
                        fij.append(None)
                        c95ij.append(None)

    if (turn_iso == 1):
        
        # Coherence of chosen isotherm with wind and radiation
        phw, caw, faw, c95aw = mod.coherence_shift(iso[aux_iso], iw, windsizeuser, dt)
        if rad == 1:
            phr, car, far, c95ar = mod.coherence_shift(iso[aux_iso], ra, windsizeuser, dt)

        phiso = []
        coiso = []
        friso = []
        cliso = []
        for i in range(4):
            for j in range(4):
                if (i != j and j > i):
                    if (tau[i] != -999) and (tau[j] != -999):
                        phaux, coaux, fiaux, c9aux = mod.coherence_shift(iso[i], iso[j], windsizeuser, dt)
                        phiso.append(phaux)
                        coiso.append(coaux)
                        friso.append(fiaux)
                        cliso.append(c9aux)
                    else:
                        phiso.append(None)
                        coiso.append(None)
                        friso.append(None)
                        cliso.append(None)

    # Contourf correction and mean temperature at thermocline depth
    templ, hl = graph.correction_contourf(h, temp, dj, lin, dci, qt)
    riwl, hlm = graph.correction_contourf(hzmid, riw2d, dj, lin, dci, qt - 1)
    buoy, hlm = graph.correction_contourf(hzmid, buoy, dj, lin, dci, qt - 1)

    tt, ithermo = mod.find_nearest(hl, m_ht)

    ri_wind = np.nanmin(mod.average(riw, (0.25 / 2) * (dura / (60 * 60) / dt)))
    we_wind = np.nanmin(mod.average(wedd, (0.25 / 2) * (dura / (60 * 60) / dt)))
 
    # Thorpe scales
    thermo_temp = templ[:, ithermo]
    tho, Lt = mod.thorpe_scale(temp, h)


    # Spectral analysis of temperature at thermocline depth 
    tthermo, pthermo, powerthermo = mod.wave_spectral(thermo_temp, dt, mother)
    freqthermo, wlthermo, welchthermo, wrthermo, confthermo = mod.welch_method(thermo_temp, windsizeuser, window, dt)
    freq_ean, wl_aper_ean, welch_ean, wr_ean, conf_ean = mod.welch_method(ean, windsizeuser, window, dt)
    
   

    print(">         Spectral analysis was computed")
    print("> ")
    root.update()

    anal_time = time.time()
    print("> Execution time for part IV: ")
    print("> " + str(round(anal_time - load_time, 4)) + ' seconds')
    root.update()

    print("--------------------------------------------------------------------------------------")
    print("> ")
    root.update()
    
    # ----------------- Generating results ------------------------------------
    
    print("> Part V       Generating text files... ")
    print("> ")
    root.update()
    
    # Create textfiles directory
    text_dir = os.path.join(output_path, 'textfiles')
    
    if os.path.exists(text_dir) and any(os.scandir(text_dir)):
        msg = "[WARN] Previous textfile outputs were overwritten.\n"
        dig.write(msg)
    else:
        os.makedirs(text_dir, exist_ok=True)
    

    def safe_savetxt(path, arr, header='', fmt='%0.8f', delimiter='\t', comments=''):

        try:
            if arr is None:
                raise ValueError("array is None")
            a = np.asarray(arr)
            
            if a.ndim == 1:
                a = a.reshape(-1, 1)
            
            if a.size == 0:
                raise ValueError("array is empty")
            np.savetxt(path, a, delimiter=delimiter, header=header, fmt=fmt, comments=comments)
            return True
        except Exception as e:
            msg = f"[WARN] Failed to save {os.path.basename(path)}: {e}\n"
            try:
                dig.write(msg)
            except Exception:
                pass

            return False
    
    # -------------------- Build list of save tasks ---------------------------
    save_tasks = []  # each entry: (path, array_like, header, fmt)
    
    # Isotherms
    if turn_iso == 1:
        for i in range(4):
            if tau[i] != -999:
                save_tasks.append((os.path.join(text_dir, f'spectral_isotherms{int(tau[i])}.txt'),
                                   np.column_stack((freq[i], welch[i])) if (isinstance(freq[i], (list, np.ndarray)) and isinstance(welch[i], (list, np.ndarray))) else None,
                                   'freq(Hz)\tPSD(m2/Hz)',
                                   '%0.8f %0.15f'))
                save_tasks.append((os.path.join(text_dir, f'isotherms{int(tau[i])}.txt'),
                                   np.column_stack((time_temp[i], iso[i])) if (isinstance(time_temp[i], (list, np.ndarray)) and isinstance(iso[i], (list, np.ndarray))) else None,
                                   'time(hour)\tisotherm(m)',
                                   '%0.8f %0.5f'))
                save_tasks.append((os.path.join(text_dir, f'bandpass_iso{int(tau[i])}.txt'),
                                   np.column_stack((time_temp[i], band[i])) if (isinstance(time_temp[i], (list, np.ndarray)) and isinstance(band[i], (list, np.ndarray))) else None,
                                   'time(hour)\tisotherm(m)',
                                   '%0.8f %0.5f'))
                save_tasks.append((os.path.join(text_dir, f'rednoise_isotherms{int(tau[i])}.txt'),
                                   np.column_stack((wr[i], conf[i])) if (isinstance(wr[i], (list, np.ndarray)) and isinstance(conf[i], (list, np.ndarray))) else None,
                                   'freq(Hz)\tPSD(m2/Hz)',
                                   '%0.8f %0.15f'))
    
    # Sensors, use seu indices and check presence
    if turn_temp == 1:
        for k in range(4):
            if sen[k] == 1 and (isinstance(seu[k], (int, np.integer)) and seu[k] >= 0):
                idx = int(seu[k])
                # spectral_sensor
                try:
                    arr_spec = np.column_stack((freqe[idx], welch_sensor[idx])) if (isinstance(freqe[idx], (list, np.ndarray)) and isinstance(welch_sensor[idx], (list, np.ndarray))) else None
                except Exception:
                    arr_spec = None
                save_tasks.append((os.path.join(text_dir, f'spectral_sensor{idx}.txt'),
                                   arr_spec,
                                   'freq(Hz)\tPSD(m2/Hz)',
                                   '%0.8f %0.15f'))
                # Sensor raw data (timee and temp series)
                try:
                    arr_sensor = np.column_stack((timee[idx], temp[:, idx])) if (isinstance(timee[idx], (list, np.ndarray)) and temp.shape[1] > idx) else None
                except Exception:
                    arr_sensor = None
                save_tasks.append((os.path.join(text_dir, f'sensor{idx}.txt'),
                                   arr_sensor,
                                   'time(hour)\ttemperature(oC)',
                                   '%0.8f %0.5f'))
                # Bandpass filtered sensor
                try:
                    arr_band = np.column_stack((timee[idx], sensor_filtered[idx])) if (isinstance(timee[idx], (list, np.ndarray)) and isinstance(sensor_filtered[idx], (list, np.ndarray))) else None
                except Exception:
                    arr_band = None
                save_tasks.append((os.path.join(text_dir, f'bandpass_sen{idx}.txt'),
                                   arr_band,
                                   'time(hour)\ttemperature(oC)',
                                   '%0.8f %0.5f'))

    
    # Isotherm pairs (coherence pairs)
    r = 0
    for i in range(4):
        for j in range(4):
            if (i != j and j > i) and (tau[i] != -999) and (tau[j] != -999):
                path = os.path.join(text_dir, f'iso_spectralpair_{int(tau[i])}_{int(tau[j])}.txt')
                try:
                    arr = np.column_stack((friso[r], phiso[r], coiso[r], cliso[r]))
                except Exception:
                    arr = None
                save_tasks.append((path, arr, 'frequency(hour)\t phase\t coherence\t coherenceg95', '%0.3f %0.8f'))
                r += 1
    

    time_win = np.asarray(time_win)
    riwl     = np.asarray(riwl)
    riwl_to_save = np.column_stack([time_win, riwl])
    
    save_tasks.append((
        os.path.join(text_dir, 'richardson.txt'),
        riwl_to_save,
        '000000\t' + '\t'.join(map(str, hlm)) if hlm is not None else '',
        '%0.15e'
    ))
    
    
    save_tasks.append((os.path.join(text_dir, 'wind.txt'),
                       np.column_stack((time_win, dw, iw, strs)) if (isinstance(time_win, (list, np.ndarray)) and isinstance(dw, (list, np.ndarray))) else None,
                       'time(hour)\tdirection(o)\tspeed(m/s)\tstress(N/m2)',
                       '%0.8f %0.1f %0.3f %0.6f'))
    save_tasks.append((os.path.join(text_dir, 'watercond_mode1.txt'),
                       np.column_stack((time_win, pe, ph, he, hh)) if (isinstance(time_win, (list, np.ndarray)) and isinstance(pe, (list, np.ndarray))) else None,
                       'time(hour)\tupper density (kg/m3)\tlower density (kg/m3)\tupper thickness (m)\tlower thickness (m)',
                       '%0.8f %0.6f %0.6f %0.2f %0.2f'))
    
    # spectral_wind optional
    try:
        save_tasks.append((os.path.join(text_dir, 'spectral_wind.txt'),
                           np.column_stack((wl_aper_win, welch_win)) if (isinstance(wl_aper_win, (list, np.ndarray)) and isinstance(welch_win, (list, np.ndarray))) else None,
                           'period(hour)\tPSD wind((m/s)2/Hz)',
                           '%0.3f %0.5f'))
    except Exception:
        pass
    

    
    save_tasks.extend([
        (os.path.join(text_dir, 'stability.txt'),
         np.column_stack((time_win, riw, wedd, iw_dw, iw_up, wedi, schmidt, lakeNumber)) if (isinstance(time_win, (list, np.ndarray))) else None,
         'time(hour)\t Ri(-)\tW (-)\tWmin (-)\tWmax (-)\tWfilt(-)\tSchmidt(J/m2)\tLakeNumber(-)',
         '%0.3f %0.8f %0.8f %0.8f %0.8f %0.8f %0.8f %0.8f'),
        (os.path.join(text_dir, 'buoyancy.txt'),
         np.column_stack((time_win, n, n_slope)) if (isinstance(time_win, (list, np.ndarray)) and isinstance(n, (list, np.ndarray))) else None,
         'time(hour)\t buoyancy frequency 1 (Hz)\t mean buoyancy frequency 2 (Hz)',
         '%0.3f %0.8f %0.8f'),
        (os.path.join(text_dir, 'thermocline.txt'),
         np.column_stack((time_win, he)) if (isinstance(time_win, (list, np.ndarray)) and isinstance(he, (list, np.ndarray))) else None,
         'time(hour)\tupper layer thickness(m)',
         '%0.3f %0.4f'),
        (os.path.join(text_dir, 'basin_length.txt'),
         np.column_stack((time_win, ls_fetch)) if (isinstance(time_win, (list, np.ndarray)) and isinstance(ls_fetch, (list, np.ndarray))) else None,
         'time(hour)\tbasin length at wind direction(m)',
         '%0.3f %0.4f'),
        (os.path.join(text_dir, 'thermocline_temperature.txt'),
         np.column_stack((time_win, thermo_temp)) if (isinstance(time_win, (list, np.ndarray)) and isinstance(thermo_temp, (list, np.ndarray))) else None,
         'time(hour)\ttemperature (oC)',
         '%0.3f %0.4f'),
        (os.path.join(text_dir, 'metalimnion_thickness.txt'),
         np.column_stack((time_win, h2)) if (isinstance(time_win, (list, np.ndarray)) and isinstance(h2, (list, np.ndarray))) else None,
         'time(hour)\tmetalimnion thickness(m)',
         '%0.3f %0.4f'),
        (os.path.join(text_dir, 'mean_profile.txt'),
         np.column_stack((mean_h, mean_temp, low_temp, high_temp, low_temp_sd, high_temp_sd)) if (isinstance(mean_h, (list, np.ndarray))) else None,
         'mab\ttemperature (dC)\tmin conf95 (dC)\tmax conf95 (dC)\tmin sd (dC)\tmax sd (dC)',
         '%0.3f %0.4f %0.4f %0.4f %0.4f %0.4f'),
        (os.path.join(text_dir, 'mean_buoyancy.txt'),
         np.column_stack((mean_hm, mean_buoy, low_buoy, high_buoy, low_buoy_sd, high_buoy_sd)) if isinstance(mean_hm, (list, np.ndarray)) else None,
         'mab\tbuoyancy frequency (Hz)\tmin conf95 (Hz)\tmax conf95 (Hz)\tmin sd (Hz)\tmax sd (Hz)',
         '%0.3f %0.6f %0.6f %0.6f %0.6f %0.6f'),
        (os.path.join(text_dir, 'wd_event.txt'),
         np.column_stack((time_win, dw_lit, dw_spi, dw_hom, fdire)) if (isinstance(time_win, (list, np.ndarray))) else None,
         'time(hour)\tliterature(o)\tSpigel(o)\thomogeneous(-)\tfhomog(-)',
         '%0.3f %0.2f %0.2f %0.2f %0.5f'),
        (os.path.join(text_dir, 'internalseiche_periods.txt'),
         np.column_stack((time_win, period_time[:,0], period_time[:,1], period_time[:,2], period_time[:,3], period_time[:,4])) if (len(time_win) > 0) else None,
         'time(hour)\tV1H1 (h)\tV2H1 (h)\tV3H1 (h)\tV4H1 (h)\tV5H1 (h)',
         '%0.3f %0.2f %0.2f %0.2f %0.2f %0.2f')
    ])

    graph.save_mode_timeseries(time_model, hemod, mode2_nodes, mode3_nodes,
    mode4_nodes, mode1Layer, mode2Layer, mode3Layer, mode4Layer, text_dir, 
    save_tasks)
    
    # Decomposition outputs (only if variables exist)
    try:
        save_tasks.append((os.path.join(text_dir, 'mab_decomp.txt'), np.column_stack((refined_depth.T)), 'water depth (m)\t', '%0.3f '))
        save_tasks.append((os.path.join(text_dir, 'time_decomp.txt'), np.column_stack((time_model.T)), 'time(hour)\t', '%0.3f '))
        save_tasks.append((os.path.join(text_dir, 'mab_decomp_oiginal.txt'), np.column_stack((h_deco)), None, '%0.3f '))
    except Exception:
        pass
    
    # velocity/cpzin arrays if present (use try to avoid crash)
    try:
        save_tasks.append((os.path.join(text_dir, 'uarbit_decomp_mode1.txt'), np.column_stack((vel_time[:,:,0])), None, '%0.15e '))
        save_tasks.append((os.path.join(text_dir, 'uarbit_decomp_mode2.txt'), np.column_stack((vel_time[:,:,1])), None, '%0.15e '))
        save_tasks.append((os.path.join(text_dir, 'uarbit_decomp_mode3.txt'), np.column_stack((vel_time[:,:,2])), None, '%0.15e '))
        save_tasks.append((os.path.join(text_dir, 'uarbit_decomp_mode4.txt'), np.column_stack((vel_time[:,:,3])), None, '%0.15e '))
        save_tasks.append((os.path.join(text_dir, 'uarbit_decomp_mode5.txt'), np.column_stack((vel_time[:,:,4])), None, '%0.15e '))
        save_tasks.append((os.path.join(text_dir, 'cpzinho_mode1.txt'), np.column_stack((cpzin[:,0,:])), None, '%0.15e '))
        save_tasks.append((os.path.join(text_dir, 'cpzinho_mode2.txt'), np.column_stack((cpzin[:,1,:])), None, '%0.15e '))
        save_tasks.append((os.path.join(text_dir, 'cpzinho_mode3.txt'), np.column_stack((cpzin[:,2,:])), None, '%0.15e '))
        save_tasks.append((os.path.join(text_dir, 'cpzinho_mode4.txt'), np.column_stack((cpzin[:,3,:])), None, '%0.15e '))
        save_tasks.append((os.path.join(text_dir, 'cpzinho_mode5.txt'), np.column_stack((cpzin[:,4,:])), None, '%0.15e '))
    except Exception:
        pass
    
    save_tasks.append((os.path.join(text_dir, 'thermo_wavelet.txt'),np.asarray(powerthermo),None,'%0.15e '))
    save_tasks.append((os.path.join(text_dir, 'thermo_wavelet_period.txt'), np.asarray(pthermo).reshape(-1, 1), 'period(hour)','%0.15e ' ))
    save_tasks.append((os.path.join(text_dir, 'thermo_psd.txt'), np.column_stack((wlthermo, welchthermo)),'Period(hour)\tPSD(oC2/Hz)','%0.15f '))
    save_tasks.append((os.path.join(text_dir, 'thermo_psd_siginificance.txt'),np.column_stack((wrthermo, confthermo)),'Period(Hz)\tConf(oC2/Hz)','%0.15f '))
        
    # ---------- Progress setup for GUI --------------------------------------
    # create a single progress line mark in the outputPanel for updates
    try:
        if "progress_mark" not in outputPanel.mark_names():
            outputPanel.insert(END, "\n")  
            outputPanel.mark_set("progress_mark", "end-1l linestart")  
    except Exception:
        pass
    
    total_tasks = len(save_tasks)
    progress_every = max(1, total_tasks // 100)  # update ~1% steps
    
    # Execute save tasks with progress
    for idx, (path, arr, header, fmt) in enumerate(save_tasks, start=1):
        
        prog_text = f"> Saving files: {idx}/{total_tasks} ({(idx/total_tasks)*100:.1f}%)"
        try:
            if "progress_mark" in outputPanel.mark_names():
                try:
                    outputPanel.delete("progress_mark linestart", "progress_mark lineend")
                except Exception:
                    pass
                outputPanel.insert("progress_mark linestart", prog_text + "\n")
                outputPanel.see("progress_mark")
            else:
                outputPanel.insert(END, prog_text + "\n")
                outputPanel.see(END)
        except Exception:
            print(prog_text, end="\r", flush=True)
    
        safe_savetxt(path, arr, header=header if header else '', fmt=fmt)
    
        if (idx % progress_every) == 0:
            try:
                root.update()
            except Exception:
                pass
    
    print("> ")
    try:
        if "progress_mark" in outputPanel.mark_names():

            try:
                outputPanel.delete("progress_mark linestart", "progress_mark lineend")
            except Exception:
                pass
            outputPanel.insert("progress_mark linestart", f"> Saved {total_tasks} files\n")
            outputPanel.see("progress_mark")
        else:
            outputPanel.insert(END, f"> Saved {total_tasks} files\n")
            outputPanel.see(END)
            
    except Exception:
        print(f"> Saved {total_tasks} files")
    
    print("> ")
    
    # Plots  
    print("> Part VI       Plotting graphs and making reports... ")
    root.update()
    
    if depi > 400:
        print(">               Attention: Figures of high quality (DPI > 400)")
        print(">               You can reduce DPI for a faster plotting")
        root.update()
    
    print(">               This may take few minutes")
    root.update()

    
    # Helper function to save figure safely
    def save_fig_safe(fname, dpi=depi):
        try:
            plt.tight_layout()
            plt.savefig(fname, dpi=dpi)
        except Exception as e:
            msg = f"[WARN] Failed saving figure {fname}: {e}\n"
            try: dig.write(msg)
            except: pass
        finally:
            plt.close()
    
    # ---------------- Isotherm and wavelet analysis ----------------
    if turn_iso == 1:
        for i in range(4):
            if tau[i] == -999:
                continue
    
            fig = plt.figure(figsize=(10, 5))
            gs = gridspec.GridSpec(2, 1)
    
            ax1 = plt.subplot(gs[0, 0]) 
            ax2 = plt.subplot(gs[1, 0], sharex=ax1)
    
            ax1.set_title('a)', loc='left')
            ax2.set_title('b)', loc='left')
    
            graph.isotherm(dx, time_temp[i], iso[i], 'navy', tau[i], ax1)
    
            # call wavelet_iso and let it grab fig if needed
            graph.wavelet_iso(dx, per_temp[i], power[i], ax2, fig=fig)
    
            # hide x tick labels on top plot (a)
            plt.setp(ax1.get_xticklabels(), visible=False)
            ax1.xaxis.set_tick_params(which='both', length=0)
    
            fig.tight_layout()
            save_fig_safe(os.path.join(output_path, f'iso{i}.png'), dpi=depi)

    # ---------------- Thermocline Analysis ----------------
    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(2, 3)
    
    ax1 = plt.subplot(gs[0, :2])
    ax2 = plt.subplot(gs[1, :2], sharex=ax1)
    ax3 = plt.subplot(gs[1, 2])
    
    ax1.set_title('a)', loc='left')
    ax2.set_title('b)', loc='left')
    ax3.set_title('c)', loc='left')
    
    graph.thermocline(dx, tthermo, thermo_temp, 'navy', round(m_ht, 2), ax1)
    graph.wavelet_iso(dx, pthermo, powerthermo, ax2, fig=fig)
    graph.psd_thermocline(wlthermo, welchthermo, 'navy', int(tt), wrthermo, confthermo, ax3)
    
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax1.xaxis.set_tick_params(which='both', length=0)
    
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'thermocline_analysis.png'), dpi=depi)
    plt.close(fig)   
    
    # ---------------- Temporal Analysis (Wind) ----------------
    fig = plt.figure(figsize=(10, 3))
    gs = plt.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0,0])
    
    graph.wind_direction(dx, time_win, dw, dw_spi, dw_lit, dw_hom,
        wedd_lim_upper, wedd_lim_lower, ax1)
    
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'temporal_analysis.png'), dpi=depi)  
    
    # ---------------- Wedderburn & Wind Stress ----------------
    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(2, 1)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    
    ax1.set_title('a)', loc='left')
    ax2.set_title('b)', loc='left')
    
    graph.windstress(dx, time_win, strs, ax1)
    graph.wedderburn(dx, time_win, wedd, lakeNumber, ax2)
    
    plt.setp(ax1.get_xticklabels(), visible=False)
    
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'stability.png'), dpi=depi)

    # ---------------- Schmidt stability --------------------
    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(1, 1)
    
    ax1 = fig.add_subplot(gs[0, 0])   

    graph.schmidt(dx, time_win, schmidt, ax1)
    
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'schmidt.png'), dpi=depi)
    
    # ---------------- Meteo Spectra ----------------
    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    
    graph.psd_wind(wl_aper_win, welch_win, wr_win, conf_win, ax1)
    
    if rad == 1:
        ax2 = ax1.twinx()
        graph.psd_sola(wl_aper_sol, welch_sol, wr_sol, conf_sol, ax2)
    
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'meteo_spectra.png'), dpi=depi)
    
    # ---------------- Level Spectra ----------------
    
    if ean_serie == 2:
        fig = plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 1, figure=fig)   
        ax1 = fig.add_subplot(gs[0, 0])
    
        graph.psd_level(freq_ean, welch_ean, wr_ean, conf_ean, ax1)
    
        fig.tight_layout()
        save_fig_safe(os.path.join(output_path, 'level_spectra.png'), dpi=depi)
    
    # ---------------- PSD & Coherence for Isotherms ----------------
    if turn_iso == 1:

        fig = plt.figure(figsize=(10, 6))
        gs = gridspec.GridSpec(2, 1, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
        
        ax1.set_title('a)', loc='left')
        ax2.set_title('b)', loc='left')
        graph.psd_iso(tau, freq, welch, f11, f21, f31, largelen, m_n, fo, wr, conf, ax1)
        graph.coherence(caw, faw, tau, aux_iso, largelen, ax2)
        
        ax1.set_ylabel('PSD isotherms (m²/Hz)')
        ax2.set_xlabel('Frequency (Hz)')
        plt.setp(ax1.get_xticklabels(), visible=False)
        
        fig.tight_layout()
        save_fig_safe(os.path.join(output_path, 'psd_nonhydro.png'), dpi=depi)
    

        fig = plt.figure(figsize=(10, 6))
        gs = gridspec.GridSpec(1, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
        
        ax1.set_title('a)', loc='left')
        ax2.set_title('b) Coriolis correction', loc='left')
        
        graph.psd_iso(tau, freq, welch, f11,  f21,  f31,  largelen, m_n, fo, wr, conf, ax1)
        graph.psd_iso(tau, freq, welch, cf11, cf21, cf31, largelen, m_n, fo, wr, conf, ax2)
        
        ax1.set_xlabel('Frequency (Hz)')
        ax2.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('PSD of isotherms (m²/Hz)')
        plt.setp(ax2.get_yticklabels(), visible=False)
        
        fig.tight_layout()
        save_fig_safe(os.path.join(output_path, 'psd_hydro_coriois.png'), dpi=depi)
    
        fig = plt.figure(figsize=(10, 6))
        gs = gridspec.GridSpec(2, 1, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
        
        ax1.set_title('a)', loc='left')
        ax2.set_title('b)', loc='left')
        
        graph.psd_iso(tau, freq, welch, f11,  f21,  f31,  largelen, m_n, fo, wr, conf, ax1)
        graph.psd_variance(tau, freq, welch, m_n, fo, wr, conf, ax2)
        
        ax1.set_ylabel('PSD isotherms (m²/Hz)')
        
        
        plt.setp(ax1.get_xticklabels(), visible=False)
        fig.tight_layout()
        save_fig_safe(os.path.join(output_path, 'psd.png'), dpi=depi)
    
    # ---------------- Degeneration, Classification, Evolution ----------------
    fig = plt.figure(figsize=(10, 7))
    gs = gridspec.GridSpec(3, 2, figure=fig, width_ratios=(2, 1))
    

    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])
    ax4 = fig.add_subplot(gs[2, 1])
    
    ax1.set_title('a)', loc='left')
    ax2.set_title('b) - P1', loc='left')
    ax3.set_title('c) - P2', loc='left')
    ax4.set_title('d) - P3', loc='left')
    

    graph.degeneration(P1, P2, P3, h_lim1, h_lim2, h_lim3,
                      W_lim1, W_lim2, W_lim3, m_pe, m_ph, m_h2,
                      m_he + m_hh, m_ls, ax1)
    
    graph.generation(genera_0, ax2)
    graph.generation(genera_1, ax3)
    graph.generation(genera_2, ax4)
    
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'degenera.png'), dpi=depi)
    
    # ---------------- Classification & Bueno Parameterization ----------------
    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(1, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    mode1 = pv1h1[1] / 3600 / 2
    average_wedda = np.nanmin(mod.average(wedd, mode1 * 0.25 / 2 / dt))
    parh = m_he / (m_he + m_hh)

    g = 12.156 * (parh ** 3) - 15.714 * (parh ** 2) + 2.8426 * parh + 2.0846
    f = g * np.exp((parh ** 2) / 0.25)
    paramet = 2 * f ** 2 * np.log(amax / (0.1 * m_he))

    graph.classification_genera(average_wedda, we_wind, parh, amax / m_he, aind, ax1)
    graph.bueno_parameterization(average_wedda, we_wind, paramet, aind, ax2)

    ax1.set_title("a)", loc='left', fontsize=10)
    ax2.set_title("b)", loc='left', fontsize=10)

    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'classification_evolution.png'), dpi=depi)
    
    # --------------- Degeneration evolution ----------------------------------
    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(1, 2, figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    ax1.set_title('a) Original Chart', loc='left')
    ax2.set_title('b) Zoomed-in View', loc='left')

    graph.degeneration_evolution(P1, P2, P3, hH_gp, wi_gp, m_pe, m_ph, m_h2, m_he + m_hh, m_ls, 'no', ax1)
    graph.degeneration_evolution(P1, P2, P3, hH_gp, wi_gp, m_pe, m_ph, m_h2, m_he + m_hh, m_ls, 'yes', ax2)

    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'degeneration_evolution.png'), dpi=depi)

    # ------------------- Buoyancy Frequency ----------------------------------
    fig = plt.figure(figsize=(10, 4))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    
    graph.buoyancy2d(dx, dx_mod, time_win, hlm, buoy, ht, ean, z0, hemod, ax1)

    plt.setp(ax1.get_xticklabels(), rotation=10, ha='right')
    ax1.grid(ls=':', color='gray', lw=0.3)

    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'buoyancy.png'), dpi=depi)   

    # --------------- Wedderburn and Richardson Numbers -----------------------
    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2], figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    
    ax1.set_title('a)', loc='left', fontsize=10)
    ax2.set_title('b)', loc='left', fontsize=10)
    
    graph.wedd_limit(dx, time_win, wedd, iw_dw, iw_up, ax1)
    graph.richardson2d(dx, time_win, hlm, riwl, ht, ean, z0, ax2)
    
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax1.grid(ls=':', color='gray', lw=0.3)
    ax2.grid(ls=':', color='gray', lw=0.3)

    fig.tight_layout()
    save_fig_safe(os.path.join(output_path,'richardson.png'), dpi=depi) 

    
    # ----------------- Time-series of water density --------------------------
    try:
  
        fig = plt.figure(figsize=(10, 3))
        gs = gridspec.GridSpec(1, 1, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])
        
        graph.density(dx, time_win, pu, pd, pe, ph, ax1)
    
        fig.tight_layout()
        fig.savefig(os.path.join(output_path, 'structure_thermo.png'), dpi=depi)
    
    except:
        war.plt_structure_thermo(dig)
    
    # ---- Time-serie of wind speed and wind direction ------------------------
    fig = plt.figure(figsize=(10, 3))
    gs = gridspec.GridSpec(1, 1, figure=fig)
    ax = fig.add_subplot(gs[0, 0]) 
    graph.wind(dx, time_win, dw, iw, ax)  
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path,'wind.png'), dpi=depi)
    
    # ---------------------- Temperature structure ----------------------------
    fig = plt.figure(figsize=(10, 4))
    gs = gridspec.GridSpec(1, 1, figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    graph.tempstructure2d(dx, dx_mod, time_win, hl, templ, ht, ean, z0, hemod, 
                          mode2_nodes, ax)
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'tempstructure.png'), dpi=depi)

    # ------------ Temperature structure multi-layer --------------------------
    fig = plt.figure(figsize=(10, 4))
    gs = gridspec.GridSpec(1, 1, figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    graph.tempstructure2d(dx, dx_mod, time_win, hl, templ, ht, ean, z0, hemod, 
                          mode2_nodes, ax, mode='multi-layer')
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'tempstructure-multi.png'), dpi=depi)
        
    # ------------------------- Thorpe scale ----------------------------------
    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2], figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ax1.set_title("a)", loc="left")
    ax2.set_title("b)", loc="left")
    
    graph.thorpe_scale(dx, Lt, ax1)
    graph.thorpe_displacement2d(dx, h, tho, ean, z0, ax2)
    plt.setp(ax1.get_xticklabels(), visible=False)
    fig.tight_layout()
    plt.savefig(os.path.join(output_path, 'thorpe_scale.png'), dpi=depi)
    
    # ------------ Sensitivity for wave period estimation (V1H1) --------------
    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(2, 2, figure=fig)
    
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]
    titles = ['a)', 'b)', 'c)', 'd)']  
    for ax, title in zip(axes, titles):
        ax.set_title(title, loc='left')
    
    ax1.set_ylabel('Wave period (h)')
    ax3.set_ylabel('Wave period (h)')
    
    ax1, ax2, ax3, ax4 = axes
    graph.depth_sensitivity(xhe, v1h1_she, r'$h_e$', 0,    ax1)
    graph.densi_sensitivity(xpe, v1h1_spe, r'$\rho_e$', 0, ax2)
    graph.depth_sensitivity(xhh, v1h1_shh, r'$h_h$', 1,    ax3)
    graph.densi_sensitivity(xph, v1h1_sph, r'$\rho_h$', 1, ax4)
    
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'sensitivity_fundamental.png'), dpi=depi)
    
    # ------------ Sensitivity internal wave parameters  ----------------------
    fig = plt.figure(figsize=(10, 4))
    gs  = gridspec.GridSpec(1, 2, wspace=0.30)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    
    ax1.set_title("a)", loc="left")
    ax2.set_title("b)", loc="left")
    
    graph.parabar_sensitivity(dep_bar, Pdep, "dep", pv1h1[0], Hbar, 
                              pv1h1[2], ax1)
    graph.parabar_sensitivity( rho_bar, Prho, "rho", pv1h1[0], Pbar, 
                              pv1h1[2], ax2)   
    
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'sensitivity_parameters.png'), dpi=depi)
    
    # ------------------ Mean temperature profile -----------------------------
    plt.figure(figsize=(10,4))
    ax1 = plt.subplot2grid((1,2),(0,0))
    ax2 = plt.subplot2grid((1,2),(0,1), sharey=ax1)
    ax2.tick_params(labelleft=False)
    ax1.set_title('a)',loc='left')
    ax2.set_title('b)',loc='left')
    graph.averageTemperature(mean_temp,low_temp,low_temp_sd,high_temp,high_temp_sd,mean_h,z0,m_ean,ax1)
    graph.averageBuoyancy(mean_buoy,low_buoy,low_buoy_sd,high_buoy,high_buoy_sd,mean_hm,z0,m_ean,ax2)
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path,'mean_temperature.png'), dpi=depi)
    
    # ------------ Wind-direction filtered Wedderbun number -------------------
    fig = plt.figure(figsize=(10, 4))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    graph.wedd_compb(dx, time_win, wedd, wedi, iw_dw, iw_up, ax1)
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path,'wedd_filtering.png'), dpi=depi)
    
    # ------- Isotherms and coherence/phase-shift between isotherms -----------
    if turn_iso == 1:

        fig = plt.figure(figsize=(10, 4))
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        graph.multi_isotherms(dx, iso, tau, time_temp, z0, ax1)
        fig.tight_layout()
        save_fig_safe(os.path.join(output_path,'isotherms.png'), dpi=depi)

        if ana1 != -999 or ana2 != -999:
    
            fig = plt.figure(figsize=(10, 4))
            gs = gridspec.GridSpec(1, 1)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = ax1.twinx()
            graph.coherence_iso(tau, coiso, friso, ana1, ana2, ax1)
            graph.phase_iso(tau, phiso, friso, cliso, ana1, ana2, ax2)
            fig.tight_layout()
            save_fig_safe(os.path.join(output_path,'coherence.png'), dpi=depi)
            
    
    # ----------------- Wind resonance on internal seiche ---------------------
    fig = plt.figure(figsize=(10, 4))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    graph.wavelet_resona(dx, per_win, power_win, period_time[:,0], period_time[:,1], ax1, fig=fig)
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, 'wind_resonance.png'), dpi=depi)
    
    # --------------- Band-pass filtered isotherms ----------------------------
    if turn_iso == 1:
        fig = plt.figure(figsize=(10, 4))
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
    
        if filter_process == 1:
            ax1.set_title(
                f'Band-pass filter – Bandwidth: '
                f'{round(pv1h1[2]/3600,1)} h to {round(pv1h1[0]/3600,1)} h',
                loc='left')
        elif filter_process == 2:
            ax1.set_title(
                f'Band-pass filter – Bandwidth: {round(low_per,1)} h to {round(high_per,1)} h',
                loc='left')
    
        graph.temp_bandpass(dx, band, tau, ax1, fig=fig)
        fig.tight_layout()
        save_fig_safe(os.path.join(output_path,'iso_bandpass.png'), dpi=depi)
    
    # --------- Temperature sensor and wavelet analysis -----------------------
    if turn_temp == 1:
        for jd in range(4):
            if seu[jd] >= 0:
    
                ids = seu[jd]
                dep = temp[:, ids]
    
                fig = plt.figure(figsize=(10, 5))
                gs = fig.add_gridspec(2, 1)
   
                ax1 = fig.add_subplot(gs[0, 0])
                ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    
                ax1.set_title("a)", loc="left")
                ax2.set_title("b)", loc="left")
    
                graph.temperature(dx, timee[ids], dep, "navy", depth[jd], ax1)
                graph.wavelet_depth(dx, per[ids], power_sensor[ids], timee[ids], ax2, fig=fig)
                
                fig.tight_layout()
                save_fig_safe(os.path.join(output_path, f"sensor{int(jd)}.png",), dpi=depi)
    
        # --------- PSD of time-series of temperature -------------------------
        fig = plt.figure(figsize=(10, 6))
        gs = gridspec.GridSpec(1, 2, figure=fig)
        
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
        
        ax1.set_title('a)', loc='left')
        ax2.set_title('b)', loc='left')
        
        graph.psd_depth(sen, depth, freqe, welch_sensor, f11, f21, f31,
            largelen, m_n, fo, ax1)
        
        graph.psd_depth(sen, depth, freqe, welch_sensor, cf11, cf21, cf31,
            largelen, m_n, fo, ax2)
        
        ax1.set_ylabel('PSD of temperature fluctuation (°C²/Hz)')
        ax1.set_xlabel('Frequency (Hz)')
        ax2.set_xlabel('Frequency (Hz)')
        plt.setp(ax2.get_yticklabels(), visible=False)
        
        fig.tight_layout()
        save_fig_safe(os.path.join(output_path, "psd_coriolis_depth.png"), dpi=depi)    
        
    
        # --------------- Band-pass filtered temperature ----------------------
        fig = plt.figure(figsize=(10, 4))
        gs = gridspec.GridSpec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        
        if filter_process == 1:
            ax1.set_title(
                f"Band-pass filter – Bandwidth: {round(pv1h1[2]/3600,1)} h to "
                f"{round(pv1h1[0]/3600,1)} h",
                loc='left')
        elif filter_process == 2:
            ax1.set_title(
                f"Band-pass filter – Bandwidth: {round(low_per,1)} h to "
                f"{round(high_per,1)} h",
                loc='left')
        
        graph.depth_bandpass(dx, s_filtered, depth, time_win, sen, ax1, fig=fig)
        
        fig.tight_layout()
        save_fig_safe(os.path.join(output_path, "depth_bandpass.png"), dpi=depi)

    # -------Thermal variation for each sensor of the thermistor chain --------
    fig = plt.figure(figsize=(10, 5))
    gs  = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    graph.thermal_variation(dx,depth,seu,temp,time_win,ax1)
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, "temperature_depth.png"), dpi=depi)

    
    # --------- PSD of temperature sensors along water depth ------------------
    fig = plt.figure(figsize=(10, 4))
    gs  = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0, 0])

    graph.psd_multilayer(freqe[0], hmean[::-1], welch_sensor, fo, m_n, 
                         pv1h1[1], pv2h1[1], z0,ax1) 
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, "psd_multi.png"), dpi=depi)
    
    # ---- Time-series of internal seiche periods from decomposition model ----
    fig = plt.figure(figsize=(10, 4))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])
    graph.modal_period(dx,period_time,ax)
    fig.tight_layout()
    save_fig_safe(os.path.join(output_path, "mode_period.png"),dpi=depi)

    
    # ----------- Decomposition velocity in arbitrary units -------------------
    try:
        for i in range(5):
            fig = plt.figure(figsize=(10, 4))
            gs = gridspec.GridSpec(1, 1)
            ax = fig.add_subplot(gs[0, 0])
    
            graph.velocity_mode(dx_mod, vel_time[:, :, i], refined_depth,
                np.mean(period_time[:, i]), i, ax, fig=fig)
    
            fig.tight_layout()
            save_fig_safe(os.path.join(output_path, f"velocity_arbitrary_mode{i}.png"),
                dpi=depi)
    except:
        war.decomposition_output(dig)
        
    # ---------- Dashboard ---------------------------------------------------
        
    print (">         Making the Interwave Analyzer's dashboard")
    root.update()
            
    # Configuring data to be analyzed via dashboard 
    
    durWave = round(100*(1-np.isnan(dw_spi).sum()/len(dx)),2)
    durDire = round(100*(1-np.isnan(dw_hom).sum()/len(dx)),2)
    
    if filter_process == 1:
        low_per  = round(pv1h1[2]/3600,1),
        high_per = round(pv1h1[0]/3600,1),            
    elif filter_process == 2:
        low_per  = round(low_per,2),
        high_per = round(high_per,2),

    # Dashboard data 

    np.savez(
        output_path+"dash_data.npz",
        
    
        # Time
        date=date,
        dx=dx,                     # shape (Nt,)
        dx_mod=dx_mod,
        nameBasin = nameBasin,
        
        numSensor = qt,
        windFreq  = dt*60,
        tempFreq  = dt*60,
              
        # Water level or references
        ean=ean,
        z0=z0,
        
        # Wind
        iw = iw,
        dw = dw,
        dura = round(dura/(60*60),2),
        dire = round(dire/(60*60),2), 
        
        m_dw_spi = round(m_dw_spi,0),
        
        fdura = round(fdura,3),
        fdire = round(np.nanmean(fdire),3),
        
        m_wast  = round(m_wast,5),
        wastmin = round(np.nanmin(wast),5),
        wastmax = round(np.nanmax(wast),5),

        # Mixing
        riw      = riw,
        ri_wind  = ri_wind,  
        wedd     = wedd,
        cond1    = cond1,        
        cond2    = cond2,
        cond3    = cond3,
        
        # Isotherms 
        iso=iso,                  # shape (Niso, Nt)
        tau=tau,                  # isotherm depths
        seu=seu,
        pe=pe,
        ph=ph,
    
        # Temperature and thermal structure
        temp=temp,                # shape (Nt, Nz)
        
        m_he = round(m_he,2),
        m_hh = round(m_hh,2),
        m_h1 = round(m_h1,2),
        m_h2 = round(m_h2,2),
        m_h3 = round(m_h3,2),

        m_pe = round(m_pe,2),
        m_ph = round(m_ph,2),
        m_p1 = round(m_p1,2),
        m_p2 = round(m_p2,2),
        m_p3 = round(m_p3,2),

        c_he = round(c_he,2),
        c_hh = round(c_hh,2),
        c_h1 = round(c_h1,2),
        c_h2 = round(c_h2,2),
        c_h3 = round(c_h3,2),

        c_pe = round(c_pe,2),
        c_ph = round(c_ph,2),
        c_p1 = round(c_p1,2),
        c_p2 = round(c_p2,2),
        c_p3 = round(c_p3,2),
        
        T11 = T11,
        T21 = T21,
        T31 = T31,
        T41 = T41,
        
        Bu = Bu, 
        cp = cp,     
        
        durWave = durWave,
        durDire = durDire,
        
        mean_h=mean_h,
        hemod=hemod,
        hemod_inter=hemod_inter,
        mode2_nodes=mode2_nodes,
        mode2_nodes_interpolated = mode2_nodes_interpolated,  
        
        
        # PSD
        fo = fo,
        m_n = m_n,
        
        low_per  = low_per,
        high_per = high_per,

        # Arrays that have none depending on user definitions
        freq  = np.array(mod.normalizeList(freq), dtype=object),
        welch = np.array(mod.normalizeList(welch), dtype=object),
        wr    = np.array(mod.normalizeList(wr), dtype=object),
        conf  = np.array(mod.normalizeList(conf), dtype=object),
        depth = np.array(mod.normalizeList(depth), dtype=object), # shape (Nz,)
        band  = np.array(mod.normalizeList(band), dtype=object),
        
        # Basin characteristics
        ls_fetch    = ls_fetch,
        type_length = type_length,
        longData    = longData,
        transData   = transData,
        
        schmidt     = schmidt,
        lakeNumber  = lakeNumber,
        
        n_slope      = n_slope,
        wave_slope   = wave_slope,
        longLeft     = longLeft,
        longRight    = longRight,
        transLeft    = transRight,
        transRight   = transRight,
        atLeft       = atLeft,
        atRight      = atRight,
        atcrossLeft  = atcrossLeft,
        atcrossRight = atcrossRight
         
    )


    port = random.randint(8050, 9000)
    if getattr(sys, 'frozen', False):
        # Running as EXE
        subprocess.Popen([
            sys.executable,
            os.path.abspath(output_path),
            str(port)
        ])
    else:
        # Running as normal Python script
        subprocess.Popen([
            sys.executable,
            "iwdash.py",
            os.path.abspath(output_path),
            str(port)
        ])
    
    # Final summary & cleanup 
    print("> ")
    root.update()
    plot_time = time.time()
    print("> Execution time for part V and VI: ")
    root.update()
    print("> " + str(round(plot_time - anal_time, 4)) + ' seconds')
    root.update()
    
    print("--------------------------------------------------------------------------------------")
    root.update()
    print('> ')
    root.update()
    print("> FINISHED            Interwave Analyzer ")
    root.update()
    print("> Check the following path for results:")
    root.update()
    print("> " + output_path)
    root.update()
    print("> For additional information:")
    root.update()
    print("> www.bit.ly/interwave_analyzer")
    root.update()



    
    # Close diagnosis file
    try:
        dig.close()
    except Exception:
        pass
    
    # Restore stdout and enter Tk mainloop
    sys.stdout = old_stdout
    try:
        root.mainloop()
    except Exception:
        pass
