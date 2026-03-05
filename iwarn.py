# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer - Warnings module

Interwave Analyzer - Version 2 (2026) 
Warnings module version: 2.260305

-------------------------------------------------------------------------------

de Carvalho Bueno, R; Bleninger, T. B.; Lorke, A. 
Internal wave analyzer for thermally stratified basins 
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

def bathyMissing(dig, filename):
    dig.write(f"[WARN] pathBathy file '{filename}' not found. Ignoring transverse basin axis\n")

def decomp_default(dig, dt):
    dig.write(f"[WARN] Decomposition model applied at each time step (dt={int(round(dt))} min)\n")

def decomp_interpolation(dig):
    dig.write("[WARN] Decomposition model outputs have been interpolated back to the original temporal grid.\n")

def decomp_changed(dig, dt):
    dig.write(f"[WARN] Decomposition resolution adjusted to match data (dt={int(round(dt))} min)\n")

def decomp_multiple(dig, dt):
    dig.write(f"[WARN] Decomposition dt adjusted to nearest multiple ({int(round(dt))} min)\n")

def decomp_specified(dig, dt):
    dig.write(f"[WARN] User-specified decomposition dt={int(round(dt))} min\n")

def coarse_fetch(dig):
    dig.write("[WARN] .fet data too coarse; fetch estimated from nearest mean direction\n")

def isotherm_boundary(dig, tau):
    dig.write(f"[WARN] {tau}°C isotherm deactivated (out of bounds entire period)\n")

def three_layer(dig, threrror):
    dig.write(f"[WARN] Metalimnion borders inefficient ({int(threrror)}%); thickness=1/3 depth\n")

def profile_structure(dig, mode, estimated):
    mode_name = "V1H1" if mode == 1 else "V2H1"
    dig.write(f"[WARN] {mode_name} period estimated from mean profiles ({int(estimated)}%)\n")

def merian(dig):
    dig.write("[WARN] Merian equation failed; Ri not filtered by wind duration\n")

def metalimnion(dig, estimated):
    dig.write(f"[WARN] Metalimnion borders undefined ({int(estimated)}%); thickness=5% depth\n")

def thermocline(dig, estimated):
    dig.write(f"[WARN] Thermocline by fast gradient approx ({int(estimated)}% affected)\n")

def threetotwo(dig, estimated):
    dig.write(f"[WARN] Three-layer replaced by two-layer model ({int(estimated)}%)\n")

def homogeneous_condition(dig):
    dig.write("[WARN] Homogeneous wind direction not identified\n")

def welch_average(dig, estimated):
    dig.write(f"[WARN] Welch averaging set to 5 days ({int(estimated)}% affected)\n")

def average(dig):
    dig.write("[WARN] Large wind event duration not found; mean value used\n")

def bandpass(dig, typ):
    if typ == 'sensor':
        dig.write("[WARN] Some sensors not band-pass filtered\n")
    if typ == 'isotherm':
        dig.write("[WARN] Some isotherms not band-pass filtered\n")

def spectral(dig, typ):
    if typ == 'sensor':
        dig.write("[WARN] Sensor spectral analysis failed (constant array)\n")
    if typ == 'isotherm':
        dig.write("[WARN] Isotherm spectral analysis failed (constant array)\n")
    if typ == 'radiation':
        dig.write("[WARN] Radiation constant; spectral analysis skipped\n")
    if typ == 'wind':
        dig.write("[WARN] Wind constant; spectral analysis skipped\n")

def plt_structure_thermo(dig):
    dig.write("[WARN] Figure 'structure_thermo.png' not generated\n")

def decomposition_overflow(dig):
    dig.write("[WARN] Decomposition unstable; results may be invalid\n")

def decomposition_output(dig):
    dig.write("[WARN] Decomposition unstable; output not generated\n")
    
