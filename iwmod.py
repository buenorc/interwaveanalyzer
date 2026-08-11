# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer - Functions Module

Interwave Analyzer - Version 2 (2026) 
Functions module version: 2.260305

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
# external modules

import math
import warnings
import numpy.matlib

import numpy as np
import scipy.special as sc

from scipy import signal, interpolate
from scipy.stats import t, sem
from scipy.stats.distributions import chi2

# internal modules
from   wavelib import wavelet

# -----------------------------------------------------------------------------
#  Numerical warning policy
# -----------------------------------------------------------------------------
#  Up to version 1.x this module executed ``warnings.simplefilter("error")``,
#  which promoted *every* warning - including the RuntimeWarnings that NumPy
#  emits for a division by zero or for the square root of a negative number -
#  into an exception.  Because those operations are unavoidable when a field
#  record contains gaps or a momentarily unstable profile, the program was
#  forced to wrap almost every formula in a try/except block, and a single
#  benign warning silently diverted the computation to a fallback value.
#
#  The policy is now the opposite: the offending operations are allowed to
#  produce NaN, and NaN is propagated through the NaN-aware reductions used
#  downstream.  A value that could not be computed therefore appears as a gap
#  in the results instead of being replaced by an undocumented substitute.
warnings.filterwarnings("ignore", category=RuntimeWarning)
np.seterr(divide='ignore', invalid='ignore', over='ignore', under='ignore')

# Gravitational acceleration (m/s2) and reference air density (kg/m3)
GRAVITY = 9.81
RHO_AIR = 1.225


def basinOrientation(data, orientation_deg):
    """
    Assigns orientation to dataset.

    Orientation represents the azimuth of the transect
    (used later for wind-fetch calculations).
    """

    data = data.copy()
    data["orientation"] = orientation_deg % 360

    return data

def fetchVariable(longData, transData, wind_dir):
    """
    Computes effective surface fetch based on wind direction.

    Inputs:
    longData  : Main basin data
    transData : Transverse basin data
    wind_dir  : Wind direction (degrees, nautical convention) [float]
        
    Outputs:
    ls_fetch : Effective surface basin length [float]
        
    """

    # Surface lengths (shallowest layer)
    L_long = longData["dists"][0]
    L_trans = transData["dists"][0]

    # Orientations
    theta_L = longData["orientation"]
    theta_T = transData["orientation"]

    # Convert to radians
    dtheta_L = np.deg2rad(wind_dir - theta_L)
    dtheta_T = np.deg2rad(wind_dir - theta_T)

    # Projection
    proj_L = abs(np.cos(dtheta_L))
    proj_T = abs(np.cos(dtheta_T))

    # Effective fetch
    ls_fetch = L_long * proj_L + L_trans * proj_T

    return ls_fetch



def computeBasinSlopes(data):
    """
    Compute segment slopes (degrees) between consecutive isolines.
    
    Inputs:
        data  : dict – dictionary containing slope arrays
                data["slope"]["left"]  : np.ndarray – left-side slopes (n-1)
                data["slope"]["right"] : np.ndarray – right-side slopes (n-1)
        depth : np.ndarray – depth array where depth[i] corresponds to depth[i+1]
    
    Outputs:
        slope_deg : np.ndarray – segment slope angles in degrees
    """

    depths = np.asarray(data["depths"])
    dists  = np.asarray(data["dists"])
    refs   = np.asarray(data["refs"])

    x_left  = refs
    x_right = refs + dists

    n = len(depths)

    slope_left  = np.zeros(n - 1)
    slope_right = np.zeros(n - 1)

    for i in range(n - 1):

        dz = depths[i + 1] - depths[i]

        dx_left = x_left[i + 1] - x_left[i]
        dx_right = x_right[i] - x_right[i + 1]

        # Left slope
        if dx_left == 0:
            slope_left[i] = 90.0
        else:
            slope_left[i] = np.degrees(np.arctan(dz / dx_left))

        # Right slope
        if dx_right == 0:
            slope_right[i] = 90.0
        else:
            slope_right[i] = np.degrees(np.arctan(dz / dx_right))

    data["slope"] = {"left": slope_left, "right": slope_right}

    return data

def slopeAtThermocline(data, z_therm):
    """
    Compute bathymetric slope (degrees) at thermocline depth
    for both basin sides.
    
    Inputs:
        data    : dict – dictionary containing bathymetric information
                  data["depths"]        : np.ndarray – depth array (m)
                  data["slope"]["left"] : np.ndarray – left-side slopes (degrees)
                  data["slope"]["right"]: np.ndarray – right-side slopes (degrees)
        z_therm : float – thermocline depth from surface (m)
    
    Outputs:
        slope_left  : float – bathymetric slope at z_therm (left side, degrees)
        slope_right : float – bathymetric slope at z_therm (right side, degrees)
    """

    depths = data["depths"]
    slopes_left = data["slope"]["left"]
    slopes_right = data["slope"]["right"]

    # If above first layer
    if z_therm <= depths[0]:
        return slopes_left[0], slopes_right[0]

    # If below deepest layer
    if z_therm >= depths[-1]:
        return slopes_left[-1], slopes_right[-1]

    # Find interval
    idx = np.searchsorted(depths, z_therm) - 1

    slope_left = slopes_left[idx]
    slope_right = slopes_right[idx]

    return slope_left, slope_right

def thermoclineSlopes(type_length,longData,transData,ean_t,hemod_t,z0):
    """
    Compute bathymetric slopes at thermocline depth
    for longitudinal and transverse transects.
    
    Inputs:
        type_length : int – basin configuration flag (1, 2, or 3)
    
        longData    : dict – longitudinal transect data
                      longData["depths"]        : np.ndarray – depth array (m)
                      longData["slope"]["left"] : np.ndarray – left-side slopes (degrees)
                      longData["slope"]["right"]: np.ndarray – right-side slopes (degrees)
    
        transData   : dict or None – transverse transect data (same structure as longData)
    
        ean_t       : float – surface elevation at time t (m)
        hemod_t     : float – epilimnion thickness at time t (m)
        z0          : float – reference vertical offset (m)
    
    Outputs:
        slope_long  : tuple(float, float) – (left, right) slopes at thermocline depth (longitudinal)
        slope_trans : tuple(float, float) or None – (left, right) slopes at thermocline depth (transverse)
    """

    # Thermocline depth from surface
    z_therm = ean_t - hemod_t

    # Longitudinal slopes
    long_left, long_right = slopeAtThermocline(longData, z_therm)

    slopes = {"long": {"left": long_left, "right": long_right}, "trans": None}

    # If type 3 = compute transverse
    if type_length == 3 and transData is not None:
        trans_left, trans_right = slopeAtThermocline(transData, z_therm)

        slopes["trans"] = {"left": trans_left, "right": trans_right}

    return slopes


def comparison_definition(c):
    """
    Map GUI label to isotherm pair code.
    
    Inputs:
        c : str – label from GUI (e.g., "Iso. 1-2")
    
    Outputs:
        code : int – corresponding isotherm pair code (e.g., 12, 13, ...)
               or -999 if label is not recognized
    """
    mapping = {
        "None": -999,
        "Iso. 1-2": 12,
        "Iso. 1-3": 13,
        "Iso. 1-4": 14,
        "Iso. 2-3": 23,
        "Iso. 2-4": 24,
        "Iso. 3-4": 34,
    }
    return mapping.get(c, c)

                 

def depths(profile: int, depth: np.ndarray, lin: int):
    """
    Get depth values for a specified sensor profile.

    Inputs:
        profile : int   – profile index (use <0 for none)
        depth   : 2D np.ndarray – depth matrix [rows, profiles]
        lin     : int   – number of rows to extract

    Outputs:
        list[np.ndarray or None] – selected depth array or None
    """
    if profile >= 0:
        aux = depth[:lin, profile].astype(float, copy=False)
    else:
        aux = None
    return [aux]

          

def commission(tau: np.ndarray):
    """
    Compute water density (kg/m³) from temperature (°C).

    Inputs:
        tau : float or np.ndarray – temperature in °C
    Outputs:
        np.ndarray – water density (kg/m³)
    """
    coeffs = [6.536336e-9, -1.120083e-6, 1.001685e-4,
              -9.09529e-3, 6.793952e-2, 999.842592]
    return np.polyval(coeffs, tau)


def interpolation(y1, y2, x1, x2, x, auxiso):
    """
    Linear interpolation between (x1, y1) and (x2, y2) for value x.

    Inputs:
        y1, y2 : float – dependent variable values
        x1, x2 : float – independent variable values
        x      : float – query point
        auxiso : float – fallback value if x is out of bounds
    Outputs:
        float – interpolated value
    """
    if y1 == y2 or x1 == x2:
        return (y1 + y2) / 2.0

    try:
        f = interpolate.interp1d([x1, x2], [y1, y2], bounds_error=True)
        return float(f(x))
    except ValueError:
        return auxiso


def isotherms(reqt, qt, h, tau, zmax, zmin, aux_iso):
    """
    Estimate depth of a target isotherm (reqt).

    Inputs:
        reqt : float – target temperature (°C)
        qt   : int   – number of points in profile
        h    : np.ndarray – depth vector
        tau  : np.ndarray – temperature profile
        zmax, zmin : float – max/min depth (unused but preserved)
        aux_iso : float – fallback interpolation value
    Outputs:
        float or None – estimated depth or None if undefined
    """
    tau_near, idx = find_nearest(tau, reqt)

    if tau_near > reqt:
        if idx >= qt - 1:
            return None
        return interpolation(h[idx], h[idx + 1], tau[idx], tau[idx + 1], reqt, aux_iso)
    else:
        if idx <= 0:
            return None
        return interpolation(h[idx - 1], h[idx], tau[idx - 1], tau[idx], reqt, aux_iso)


def thickness_decomposition(psi, zph, rhoph):
    """
    Decompose stratified profile into layer thickness and mean density.

    Inputs:
        psi   : 2D np.ndarray – modal structure (depth × mode)
        zph   : 1D np.ndarray – depth grid
        rhoph : 1D np.ndarray – density profile
    Outputs:
        (list[list[float]], list[list[float]]) – layer depths and densities
    """
    n_modes = psi.shape[1]
    n_depths = len(zph) - 1

    layer_depth, layer_rho = [], []

    for m in range(n_modes):
        laymode, layrho = [], []
        rho_acc, count = rhoph[0], 1

        for z in range(n_depths):
            if psi[z, m] * psi[z + 1, m] < 0:  # zero-crossing
                depth = interpolation(zph[z], zph[z + 1], psi[z, m], psi[z + 1, m], 0, 0)
                laymode.append(depth)
                layrho.append(rho_acc / count)
                rho_acc, count = 0.0, 0

            count += 1
            rho_acc += rhoph[z]

        layer_depth.append(laymode)
        layer_rho.append(layrho)

    return layer_depth, layer_rho

def deriva(y1, y2, x1, x2):
    """
    Compute discrete derivative (|dy| / |dx|).
    """
    dx = abs(x1 - x2)
    if dx == 0:
        dx = 1e-5
    return abs(y1 - y2) / dx


def consistency(ze, zh, h, z0):
    """
    Validate/repair metalimnion boundaries (ze upper, zh lower).

    Rules:
    - Ensure ze and zh are within [min(h), max(h)].
    - Ensure ze < zh; if inverted, try swapping, otherwise set sensible defaults.
    - Returns repaired (ze, zh, error) where error=0 means OK, 1 means repaired/adjusted.
    """
    error = 0
    hmin = float(np.nanmin(h))
    hmax = float(np.nanmax(h))

    # Guard against np.nan inputs
    if ze is None or zh is None:
        ze = (hmax + hmin) / 2.0
        zh = (hmax + hmin) / 2.0
        return ze, zh, 1

    # Clamp into range
    if not np.isfinite(ze) or ze < hmin or ze > hmax:
        ze = min(max(ze if np.isfinite(ze) else hmin, hmin), hmax)
        error = 1
    if not np.isfinite(zh) or zh < hmin or zh > hmax:
        zh = min(max(zh if np.isfinite(zh) else hmax, hmin), hmax)
        error = 1

    # Ensure ordering ze < zh
    if ze >= zh:
        # try swap if it fixes
        if zh < ze and (hmin <= zh <= hmax) and (hmin <= ze <= hmax):
            ze, zh = min(ze, zh), max(ze, zh)
            error = 1
        else:
            # fallback: place them at 1/3 and 2/3 of depth range above z0
            ze = z0 + (hmax - z0) * 2.0 / 3.0
            zh = z0 + (hmax - z0) * 1.0 / 3.0
            # ensure correct order: upper < lower
            if ze >= zh:
                ze, zh = zh, ze
            error = 1

    # final guard (very small gap)
    if zh - ze < 1e-6:
        # enforce minimal separation (1e-3 m)
        zh = ze + 1e-3
        error = 1

    return float(ze), float(zh), int(error)

def waveSlope(n_slope, Twave):
    """
    Determine wave slope.

    Inputs:
        n     : buoyancy frequency (Hz)
        Twave : wave period (hours)
    """
    Twave     = Twave*3600 # seconds
    sin_theta = 2*np.pi / (n_slope * Twave)   
    sin_theta = np.clip(sin_theta, -1, 1)
    wave_slope = np.degrees(np.arcsin(sin_theta))

    return wave_slope # degree 

def thermo_region(qt, h, tau):
    """
    Determine thermocline region (standard method).
    
    Inputs:
        qt  : int – number of sensors
        h   : np.ndarray – depth array (m)
        tau : np.ndarray – temperature profile (°C)
    
    Outputs:
        dpz : float – maximum vertical density gradient
        zt  : int – index of thermocline layer
        mid : float – mean thermocline depth (m)
        rho : np.ndarray – water density profile (kg/m³)
    """
    rho = commission(tau)
    dpz = deriva(rho[1], rho[0], h[1], h[0])
    zt = 0

    for z in range(qt - 1):
        new_dpz = deriva(rho[z + 1], rho[z], h[z + 1], h[z])
        if new_dpz > dpz:
            dpz, zt = new_dpz, z

    mid = abs(np.mean([h[zt], h[zt + 1]]))
    return dpz, zt, mid, rho



def thermocline_depth(qt, h, tau, last_values=None,
                      spike_limit=0.5,         # meters: absolute threshold to consider a candidate a spike
                      max_rate_per_step=1.0,   # meters: reject changes larger than this per single step (aggressive)
                      persist_thresh=2,        # how many consecutive detections required to accept new regime
                      max_history=9):
    """
    Determine thermocline depth (weighted) with aggressive persistence
    and rate-of-change filtering.
    
    - Rejects single-step spikes (downward or upward).
    - Accepts a new thermocline only after it persists
      persist_thresh consecutive steps.
    - Keeps small history for adaptive diagnostics.
    
    Inputs:
        qt                : int – number of sensors
        h                 : np.ndarray – depth array (m)
        tau               : np.ndarray – temperature profile (°C)
    
        last_values       : dict or None – memory structure between calls
                            (read and updated internally)
    
        spike_limit       : float – absolute spike threshold (m)
        max_rate_per_step : float – maximum allowed change per timestep (m)
        persist_thresh    : int – required consecutive detections to accept a new level
        max_history       : int – maximum number of accepted thermocline values stored
    
    Outputs:
        thermo        : float – accepted thermocline depth (m)
        rho           : np.ndarray – water density profile (kg/m³)
        error         : int – status/error flag
        current_cache : dict – updated memory structure for next call
    """
    error = 0
    dpz, z, mid, rho = thermo_region(qt, h, tau)

    # compute candidate thermo as before
    if z == 0 or z >= qt - 2:
        error = 1
        candidate = mid
    else:
        try:
            hplus = (h[z] - h[z + 2]) / 2.0
            hminu = (h[z - 1] - h[z + 1]) / 2.0

            drho = (rho[z + 1] - rho[z]) / (h[z] - h[z + 1])
            drho_plus = (rho[z + 2] - rho[z + 1]) / (h[z + 1] - h[z + 2])
            drho_minu = (rho[z] - rho[z - 1]) / (h[z - 1] - h[z])

            Dplus = hplus / (drho - drho_plus)
            Dminu = hminu / (drho - drho_minu)

            candidate = h[z + 1] * (Dplus / (Dminu + Dplus)) + h[z] * (Dminu / (Dminu + Dplus))

            if (not np.isfinite(candidate)) or (candidate < np.min(h)) or (candidate > np.max(h)):
                raise ValueError("Invalid candidate")
        except Exception:
            candidate = mid
            error = 1

    # initialize cache if needed
    cache = {} if last_values is None else dict(last_values)

    last_accepted = None
    if 'thermo' in cache and np.isfinite(cache['thermo']):
        last_accepted = float(cache['thermo'])

    # prepare persist state
    persist_state = cache.get('persist_state', None)
    # persist_state tracks a dict: {'candidate': float, 'count': int}
    if persist_state is None:
        persist_state = {'candidate': None, 'count': 0}

    # if no previous accepted value: accept candidate immediately (bootstrap)
    if last_accepted is None:
        accepted = float(candidate)
        persist_state = {'candidate': None, 'count': 0}
        error = max(error, 1)  # bootstrap: low confidence
    else:
        # compute difference from last accepted
        diff = candidate - last_accepted
        absdiff = abs(diff)

        # check hard rate-of-change guard
        if absdiff > max_rate_per_step:
            # treat as spike (too large to accept in one step)
            is_spike = True
        else:
            is_spike = absdiff > spike_limit

        if not is_spike:
            # Accept candidate (it is close enough)
            accepted = float(candidate)
            # reset persist_state
            persist_state = {'candidate': None, 'count': 0}
        else:
            # Candidate looks like a spike. Check persistence:
            # If same candidate as previous persist_state (within tiny tol), increment count
            tol = 1e-6
            prev_cand = persist_state.get('candidate', None)
            if prev_cand is not None and abs(prev_cand - candidate) < tol:
                persist_state['count'] += 1
            else:
                # new suspicious regime: start counting
                persist_state['candidate'] = float(candidate)
                persist_state['count'] = 1

            # Accept only if persisted enough times
            if persist_state['count'] >= max(persist_thresh, 1):
                accepted = float(candidate)
                # reset persist
                persist_state = {'candidate': None, 'count': 0}
            else:
                # keep last accepted value (reject spike)
                accepted = last_accepted
                error = 1  # flagged as corrected/rejected

    # update history (store accepted values, not rejected candidates)
    history = cache.get('history', [])
    history = [float(x) for x in history if np.isfinite(x)]
    history.append(accepted)
    if len(history) > max_history:
        history = history[-max_history:]

    # final guards: clamp inside depth range
    accepted = float(min(max(accepted, float(np.min(h))), float(np.max(h))))

    # prepare return cache
    current_cache = {
        'thermo': accepted,
        'rho': rho,
        'history': history,
        'persist_state': persist_state
    }

    return accepted, rho, int(error), current_cache

def metathick(qt, h, tau, minval, z0):
    """
    Define metalimnion boundaries (metalimnion thickness) robustly.
    
    Inputs:
        qt     : int – number of sensors
        h      : np.ndarray – depth array (m)
        tau    : np.ndarray – temperature profile (°C)
        minval : float – density gradient threshold (kg/m³/m)
        z0     : float – reference vertical level (m)
    
    Outputs:
        ze    : float – upper metalimnion limit (m)
        zh    : float – lower metalimnion limit (m)
        error : int – status flag (0 if consistent, 1 otherwise)
    """
    error = 0

    # Compute density and gradient (vectorised)
    tau = np.asarray(tau, dtype=float)[:qt]
    h = np.asarray(h, dtype=float)[:qt]

    rho = commission(tau)
    dh = np.diff(h)
    dh = np.where(dh == 0, 0.01, dh)
    drho_dz = np.abs(np.diff(rho) / dh)

    # Normalize gradient for shape recognition
    grad_norm = drho_dz / (np.max(drho_dz) + 1e-9)

    # Identify continuous region where gradient > threshold fraction 
    try:
        # Use relative threshold (e.g., 30% of max) to define metalimnion zone
        rel_thr = 0.3
        mask = grad_norm >= rel_thr

        if not np.any(mask):
            raise ValueError("No strong gradient zone found")

        idx = np.where(mask)[0]
        i_start, i_end = idx[0], idx[-1]

        ze = np.mean([h[i_start], h[i_start+1]])
        zh = np.mean([h[i_end], h[i_end+1]])

        # Small consistency correction
        if ze < zh:
            ze, zh = zh, ze

    except Exception:
        # Fallback: use maximum gradient position
        error = 1
        z = np.argmax(grad_norm)
        ze = np.mean([h[max(z-1, 0)], h[z]])
        zh = np.mean([h[min(z+1, qt-1)], h[z]])

    # Apply absolute minval threshold refinement 
    try:
        # Expand or contract based on true density gradient threshold
        valid = drho_dz > minval
        if np.any(valid):
            idv = np.where(valid)[0]
            ze = np.mean([h[idv[0]], h[min(idv[0]+1, qt-1)]])
            zh = np.mean([h[idv[-1]], h[min(idv[-1]+1, qt-1)]])
    except Exception:
        error = 1

    # Ensure within valid range 
    ze = np.clip(ze, z0, np.max(h))
    zh = np.clip(zh, z0, np.max(h))

    # Consistency check
    if (ze <= zh) or (ze - zh < 0.2):
        # fallback to approximate 3-layer equal division
        error = 1
        H = np.max(h)
        ze = z0 + (2/3)*(H - z0)
        zh = z0 + (1/3)*(H - z0)

    return ze, zh, error


def thermal_stability(qt, h, H, tau):
    """
    Compute thermal stratification parameters along depth.
    
    Inputs:
        qt  : int – number of depth layers
        h   : np.ndarray – depth array (m)
        H   : float – total depth (m) [unused, kept for compatibility]
        tau : np.ndarray – temperature (°C) or scalar property profile (size qt)
    
    Outputs:
        rho  : np.ndarray – water density profile (size qt, kg/m³)
        n2d  : np.ndarray – Brunt–Väisälä frequency N (size qt-1, s⁻¹)
        hmid : np.ndarray – mid-layer depths (size qt-1, m)
        glin : np.ndarray – reduced gravity profile (size qt-1, m/s²)

    Notes:
        ``n2d`` holds the buoyancy frequency N itself, in rad/s, not its
        square; the name is retained for backwards compatibility.

        The routine is fully vectorised.  In version 1.x the density was
        evaluated one sensor at a time and the layer properties in a second
        Python loop; since the function is called once per time step, that
        accounted for a substantial share of the total run time for long
        records.
    """
    tau = np.asarray(tau, dtype=float)[:qt]
    h = np.asarray(h, dtype=float)[:qt]

    rho = commission(tau)

    dh = np.abs(np.diff(h))
    dh = np.where(dh == 0, 0.01, dh)          # small safety offset

    hmid = np.abs(0.5 * (h[:-1] + h[1:]))
    glin = GRAVITY * np.abs(np.diff(rho)) / rho[1:]
    n2d = np.sqrt(glin / dh)

    return rho, n2d, hmid, glin


# =============================================================================
#  Basin hypsography
# =============================================================================

class BasinHypsography(object):
    """
    Horizontal area of the basin as a function of depth.

    The Interwave Analyzer transect files describe the basin along one or two
    orthogonal axes.  The horizontal area at a given depth is reconstructed
    from those lengths by assuming an elliptical cross-section,

        A(z) = pi/4 * L_long(z) * L_trans(z) ,

    which degenerates to a circle, A(z) = pi/4 * L_long(z)^2, when only the
    longitudinal transect is available.  This is an explicit modelling
    assumption: it reproduces the correct order of magnitude and, more
    importantly, the correct *shape* of the area-depth curve, which is what
    the Schmidt stability and the Lake Number depend on.  Where a measured
    hypsographic curve exists it should be supplied as the transect file, since
    A(z) is the only quantity used downstream.

    Inputs:
        longData : dict - longitudinal transect ("depths", "dists")
        transData: dict or None - transverse transect
        dz       : float - vertical resolution of the internal grid (m)

    Attributes:
        depth  : np.ndarray - internal depth grid (m below the surface)
        area   : np.ndarray - horizontal area at every grid depth (m2)
        A0     : float      - surface area (m2)
        Zcv    : float      - depth of the centre of volume (m)
        volume : float      - basin volume (m3)

    Notes:
        Version 1.x rebuilt this grid, interpolated both transects onto it and
        recomputed the centre of volume *inside the time loop*, once for the
        Schmidt stability and once again for the Lake Number, although the
        basin geometry does not change with time.  Precomputing it removes two
        interpolations and two reductions per time step.
    """

    def __init__(self, longData, transData=None, dz=0.1):
        self.dz = float(dz)

        depths = np.asarray(longData["depths"], dtype=float)
        dists = np.asarray(longData["dists"], dtype=float)

        order = np.argsort(depths)
        depths, dists = depths[order], dists[order]

        z_min = float(depths[0])
        z_max = float(depths[-1])
        if z_max <= z_min:
            z_max = z_min + self.dz

        self.depth = np.arange(z_min, z_max + self.dz, self.dz)

        L_long = np.interp(self.depth, depths, dists)

        if transData is None:
            self.area = 0.25 * np.pi * L_long ** 2
        else:
            t_depths = np.asarray(transData["depths"], dtype=float)
            t_dists = np.asarray(transData["dists"], dtype=float)
            t_order = np.argsort(t_depths)
            L_trans = np.interp(self.depth, t_depths[t_order],
                                t_dists[t_order])
            self.area = 0.25 * np.pi * L_long * L_trans

        self.A0 = float(self.area[0])
        area_sum = float(np.sum(self.area))
        self.Zcv = (float(np.sum(self.depth * self.area)) / area_sum
                    if area_sum > 0 else np.nan)
        self.volume = area_sum * self.dz

    def interp_temperature(self, temp, depths):
        """Interpolate a measured profile onto the internal depth grid."""
        temp = np.asarray(temp, dtype=float)
        depths = np.asarray(depths, dtype=float)

        finite = np.isfinite(temp) & np.isfinite(depths)
        if np.count_nonzero(finite) < 2:
            return np.full(self.depth.size, np.nan)

        order = np.argsort(depths[finite])
        return np.interp(self.depth, depths[finite][order],
                         temp[finite][order])

def schmidtStability(temp, depths, hypso):
    """
    Schmidt stability of the water column.

    Inputs:
        temp   : np.ndarray – water temperature profile (°C, size qt)
        depths : np.ndarray – depth of every sensor below the water surface
                 (m, size qt), increasing downwards
        hypso  : BasinHypsography – precomputed area-depth curve

    Outputs:
        St : float – Schmidt stability (J/m²)

    Method:
        Following Idso (1973),

            St = (g / A0) * integral_0^H (z - z_v) rho(z) A(z) dz ,

        with z_v the depth of the centre of volume and A0 the surface area.
        St is the work per unit surface area required to mix the water column
        to a uniform density without any change of heat content, and is
        therefore the natural measure of the strength of the stratification.
    """
    rho = commission(hypso.interp_temperature(temp, depths))

    if not np.any(np.isfinite(rho)) or not np.isfinite(hypso.Zcv) \
            or hypso.A0 <= 0:
        return np.nan

    integrand = (hypso.depth - hypso.Zcv) * rho * hypso.area * hypso.dz
    return float(GRAVITY / hypso.A0 * np.nansum(integrand))


def lakeNumber(St, uStar, metaT, metaB, rhoHyp, hypso):
    """
    Lake Number.

    Inputs:
        St     : float – Schmidt stability (J/m²)
        uStar  : float – friction velocity at the water surface (m/s)
        metaT  : float – depth of the top of the metalimnion (m below surface)
        metaB  : float – depth of the base of the metalimnion (m below surface)
        rhoHyp : float – mean hypolimnion density (kg/m³)
        hypso  : BasinHypsography – precomputed area-depth curve

    Outputs:
        Ln : float – Lake Number (dimensionless)

    Method:
        Following Imberger and Patterson (1990), the Lake Number compares the
        moment of the stabilising density distribution about the centre of
        volume with the moment of the destabilising wind stress,

            Ln = g * St_v * (z_t + z_b)
                 / ( 2 rho_h u*^2 A0^(3/2) z_v ) ,

        where St_v = St * A0 / g is the volumetric form of the stability.
        Ln >> 1 indicates a stratification strong enough to resist the wind;
        Ln of order 1 marks the onset of upwelling of hypolimnetic water at the
        upwind end of the basin.
    """
    if not np.isfinite(St) or not np.isfinite(uStar) or uStar <= 0:
        return np.nan
    if hypso.A0 <= 0 or not np.isfinite(hypso.Zcv) or hypso.Zcv <= 0:
        return np.nan
    if not np.isfinite(rhoHyp) or rhoHyp <= 0:
        return np.nan

    St_volume = St * hypso.A0 / GRAVITY

    return float(GRAVITY * St_volume * (metaT + metaB)
                 / (2.0 * rhoHyp * uStar ** 2 * hypso.A0 ** 1.5 * hypso.Zcv))


def heatContent(temp, depths, hypso, cp=4186.0):
    """
    Volumetric heat content of the basin, relative to 0 °C.

    Inputs:
        temp   : np.ndarray – temperature profile (°C)
        depths : np.ndarray – sensor depths below the surface (m)
        hypso  : BasinHypsography – area-depth curve
        cp     : float – specific heat capacity of water (J/kg/K)

    Outputs:
        (total, areal)
            total : float – heat content of the whole basin (J)
            areal : float – heat content per unit surface area (J/m²)

    Method:
            Q = cp * integral_0^H rho(z) T(z) A(z) dz .

        The areal form, Q/A0, is the quantity usually plotted, because it is
        independent of the size of the basin and can be compared directly with
        the surface heat flux.
    """
    profile = hypso.interp_temperature(temp, depths)
    rho = commission(profile)

    if not np.any(np.isfinite(profile)) or hypso.A0 <= 0:
        return np.nan, np.nan

    total = float(cp * np.nansum(rho * profile * hypso.area * hypso.dz))
    return total, total / hypso.A0


def basinEnergetics(temp, depths, hypso, cp=4186.0):
    """
    Time series of heat content and available potential energy.

    Inputs:
        temp   : np.ndarray – temperature field, shape (n_times, n_sensors)
        depths : np.ndarray – depth of every sensor below the water surface,
                 either shape (n_sensors,) when the water level is constant or
                 shape (n_times, n_sensors) when it varies (m)
        hypso  : BasinHypsography – area-depth curve
        cp     : float – specific heat capacity of water (J/kg/K)

    Outputs:
        dict with

            heat        : np.ndarray – basin heat content (J)
            heat_areal  : np.ndarray – heat content per surface area (J/m²)
            ape         : np.ndarray – available potential energy (J)
            ape_areal   : np.ndarray – the same per surface area (J/m²)
            rho_ref     : np.ndarray – reference density profile (kg/m³)
            depth_grid  : np.ndarray – depth grid of ``rho_ref`` (m)

    Method:
        The available potential energy of a basin-scale internal wave field is
        the energy stored in the displacement of the isopycnals away from their
        rest position.  It is therefore measured against the *time-mean*
        stratification of the analysed period, which is the natural rest state
        of the wave field,

            rho_ref(z) = < rho(z, t) >_t ,
            APE(t) = g * integral_0^H ( rho(z,t) - rho_ref(z) ) z A(z) dz .

        Defining the reference state instead by adiabatically sorting the
        instantaneous profile - the Lorenz construction - would measure only
        the static instability of that single profile, which for a stably
        stratified lake is zero almost everywhere and carries no information
        about the internal wave field.  The Thorpe scale already reports that
        instability separately.

        The heat content is the straightforward volume integral

            Q(t) = cp * integral_0^H rho(z,t) T(z,t) A(z) dz ,

        referred to 0 °C.  Its areal form Q/A0 is directly comparable with a
        surface heat flux in W/m².
    """
    temp = np.asarray(temp, dtype=float)
    depths = np.asarray(depths, dtype=float)

    n_times = temp.shape[0]
    constant_depth = depths.ndim == 1

    grid = hypso.depth
    volume = hypso.area * hypso.dz

    profiles = np.empty((n_times, grid.size), dtype=float)
    for t in range(n_times):
        profiles[t] = hypso.interp_temperature(
            temp[t], depths if constant_depth else depths[t])

    rho = commission(profiles)

    heat = cp * np.nansum(rho * profiles * volume[None, :], axis=1)

    rho_ref = np.nanmean(rho, axis=0)
    ape = GRAVITY * np.nansum((rho - rho_ref[None, :])
                              * grid[None, :] * volume[None, :], axis=1)

    area = hypso.A0 if hypso.A0 > 0 else np.nan

    return {'heat': heat, 'heat_areal': heat / area,
            'ape': np.abs(ape), 'ape_areal': np.abs(ape) / area,
            'rho_ref': rho_ref, 'depth_grid': grid}


def windWork(stress, speed, area, dt_seconds, efficiency=0.0016):
    """
    Cumulative wind energy delivered to the basin.

    Inputs:
        stress     : np.ndarray – surface wind stress (N/m²)
        speed      : np.ndarray – wind speed at 10 m (m/s)
        area       : float      – surface area of the basin (m²)
        dt_seconds : float      – analysis time step (s)
        efficiency : float      – fraction of the wind power that reaches the
                     stratified interior

    Outputs:
        (power, cumulative)
            power      : np.ndarray – instantaneous power available to the
                         stratification (W)
            cumulative : np.ndarray – time integral of that power (J)

    Method:
        The rate at which the wind does work on the water surface is

            P = A * tau * U10 ,

        of which only a small fraction reaches the interior; the remainder is
        dissipated in the surface layer.  The default efficiency of 0.0016
        follows the value commonly adopted for the mixing efficiency of wind
        energy in lakes.  The cumulative form is the quantity to compare with
        the available potential energy, so that the two together indicate
        whether the observed internal wave field can be sustained by the
        observed forcing.
    """
    stress = np.asarray(stress, dtype=float)
    speed = np.asarray(speed, dtype=float)

    power = efficiency * area * stress * np.abs(speed)
    cumulative = np.nancumsum(np.nan_to_num(power) * dt_seconds)

    return power, cumulative

def density_2layer(qt, h, tau, H, z0, last_cache=None):
    """
    Determine two-layer thermal structure using thermocline depth with memory.

    Inputs:
        qt         : int – number of sensors
        h          : np.ndarray – depths (m)
        tau        : np.ndarray – temperature profile (°C)
        H          : float – total water depth (m)
        z0         : float – reference level (m)
        last_cache : dict – previous results for memory persistence

    Outputs:
        he, hh : float – upper and lower layer thickness (m)
        pe, ph : float – average density of epilimnion/hypolimnion (kg/m³)
        pu, pd : float – surface and bottom density (kg/m³)
        error  : int – 0 if consistent, else 1
        cache  : dict – cached structure for next iteration
    """
    thermo, rho, error, cache = thermocline_depth(qt, h, tau, last_cache)

    hh = thermo - z0
    he = H - thermo

    if any(np.isnan([he, hh])) or he <= 0 or hh <= 0:
        if last_cache:
            thermo = last_cache.get('thermo', thermo)
            rho = last_cache.get('rho', rho)
            hh = thermo - z0
            he = H - thermo
            error = 1

    pe = np.nanmean(rho[h > thermo]) if np.any(h > thermo) else rho[0]
    ph = np.nanmean(rho[h <= thermo]) if np.any(h <= thermo) else rho[-1]
    pu, pd = rho[0], rho[-1]

    return he, hh, pe, ph, pu, pd, error, cache


def density_3layer(qt, h, tau, minval, H, z0):
    """
    Determine three-layer thermal structure (epilimnion, metalimnion, hypolimnion).

    Inputs:
        qt     : int – number of sensors
        h      : np.ndarray – depths (m)
        tau    : np.ndarray – temperature profile (°C)
        minval : float – metalimnion threshold (kg/m³/m)
        H      : float – maximum depth (m)
        z0     : float – reference level (m)

    Outputs:
        h1, h2, h3 : float – thickness of each layer (m)
        p1, p2, p3 : float – mean density in each layer (kg/m³)
        error      : int – 0 if consistent, else 1
    """
    ze, zh, error = metathick(qt, h, tau, minval, z0)
    rho = commission(np.asarray(tau, dtype=float))

    try:
        # Ensure monotonic depths
        h = np.sort(h)
        idx_epi = h >= ze
        idx_meta = (h < ze) & (h > zh)
        idx_hypo = h <= zh

        p1 = np.mean(rho[idx_epi]) if np.any(idx_epi) else rho[0]
        p2 = np.mean(rho[idx_meta]) if np.any(idx_meta) else np.mean([p1, rho[-1]])
        p3 = np.mean(rho[idx_hypo]) if np.any(idx_hypo) else rho[-1]

        h1 = H - ze
        h2 = ze - zh
        h3 = zh - z0

    except Exception:
        error = 1
        h1 = h2 = h3 = (H - z0) / 3
        n = len(rho)
        p1 = np.mean(rho[: n // 3])
        p2 = np.mean(rho[n // 3 : 2 * n // 3])
        p3 = np.mean(rho[2 * n // 3 :])

    return h1, h2, h3, p1, p2, p3, error


  
        
def structure2layer(qt, h, tau, H, z0, last_cache=None):
    """
    Compute two-layer structure and stability parameters with memory persistence.

    Inputs:
        qt         : int – number of sensors
        h          : np.ndarray – depths (m)
        tau        : np.ndarray – temperatures (°C)
        H          : float – total depth (m)
        z0         : float – reference level (m)
        last_cache : dict – previous structure (optional)

    Outputs:
        he, hh  : float – layer thicknesses (m)
        pe, ph  : float – layer densities (kg/m³)
        glin    : float – reduced gravity (m/s²)
        n       : float – Brunt–Väisälä frequency (s⁻¹)
        pu, pd  : float – surface and bottom densities (kg/m³)
        error   : int – 0 if consistent, else 1
        cache   : dict – cached structure for next iteration
    """
    he, hh, pe, ph, pu, pd, error, cache = density_2layer(qt, h, tau, H, z0, last_cache)

    try:
        glin = abs(9.81 * (ph - pe) / max(ph, 1e-6))
        n = math.sqrt(abs(glin / max(he, 1e-6)))
    except Exception:
        if last_cache:
            he = last_cache.get('he', he)
            glin = last_cache.get('glin', 1e-6)
            n = last_cache.get('n', 0)
            error = 1
        else:
            glin, n, error = 0, 0, 1

    current_cache = {
        'he': he,
        'hh': hh,
        'pe': pe,
        'ph': ph,
        'glin': glin,
        'n': n,
        'pu': pu,
        'pd': pd,
        'thermo': cache.get('thermo'),
        'rho': cache.get('rho'),
    }

    return he, hh, pe, ph, glin, n, pu, pd, error, current_cache

def structure3layer(qt, h, tau, minval, H, z0):
    """
    Compute thermal structure of a three-layer system.

    Inputs:
        qt     : int – number of sensors
        h      : np.ndarray – depths (m)
        tau    : np.ndarray – temperature profile (°C)
        minval : float – threshold for metalimnion detection
        H      : float – total depth (m)
        z0     : float – reference level (m)

    Outputs:
        h1, h2, h3 : float – layer thickness (m)
        p1, p2, p3 : float – mean density of layers (kg/m³)
        error      : int – 0 if consistent, else 1
    """
    return density_3layer(qt, h, tau, minval, H, z0)

def approx_layer(he, hh, pe, ph):
    """
    Estimate three-layer structure when derivative method fails.

    Inputs:
        he : float – upper (epilimnion) depth limit (m)
        hh : float – lower (hypolimnion) depth limit (m)
        pe : float – density at epilimnion (kg/m³)
        ph : float – density at hypolimnion (kg/m³)

    Outputs:
        h1, h2, h3 : float – approximate layer thicknesses (m)
        p1, p2, p3 : float – mean densities of layers (kg/m³)
    """

    h1 = he * 0.95
    h2 = 0.05 * (he + hh)
    h3 = hh * 0.95

    p1, p3 = pe, ph
    p2 = 0.5 * (pe + ph)

    return h1, h2, h3, p1, p2, p3

def wedderburn(glin, he, wast, ls):
    """
    Compute the Wedderburn Number, indicating wind-induced stability.

    Inputs:
        glin : float – reduced gravity (m/s²)
        he   : float – epilimnion thickness (m)
        wast : float – wind friction velocity (m/s)
        ls   : float – characteristic fetch length (m)

    Outputs:
        wedd : float – Wedderburn number (dimensionless), NaN when undefined

    Notes:
        Version 1.x returned ``None`` for a degenerate input.  Assigning
        ``None`` into the float array that collects the result raises a
        ``TypeError`` in NumPy, so the analysis aborted on the first calm or
        unstratified time step.  ``NaN`` is returned instead, which the
        NaN-aware statistics downstream handle correctly.
    """
    values = (glin, wast, he, ls)
    if not all(np.isfinite(v) for v in values) or any(v <= 0 for v in values):
        return np.nan

    return glin * he**2 / (ls * wast**2)


def chi2inv(p, nfft, nperseg, test=None):
    """
    Estimate the inverse cumulative chi-square distribution (percentile).

    Inputs:
        p       : float – probability level (0–1)
        nfft    : int   – FFT size
        nperseg : int   – segment length
        test    : optional flag for gamma-based test (default None)

    Output:
        float – chi-square inverse normalized by degrees of freedom
    """
    if test is None:
        nw2 = 2 * (2.5164 * (nfft / nperseg)) * 1.2
        return chi2.ppf(p, df=nw2) / nw2
    else:
        nw2 = nfft / nperseg
        return 2 * sc.gammaincinv(nw2, p) / nw2

def stad_deviation(data):
    """
    Compute mean, confidence interval bounds, and simple standard deviation.

    Inputs:
        data : np.ndarray – 2D array (time x variable)

    Outputs:
        m, lower, upper, sd_lower, sd_upper : np.ndarray – statistical metrics
    """
    ncols = data.shape[1]
    m = np.nanmean(data, axis=0)
    sd = np.array([ciout(data[:, i]) for i in range(ncols)])
    sd_simple = np.nanstd(data, axis=0)

    lower = m - sd
    upper = m + sd
    sd_lower = m - sd_simple
    sd_upper = m + sd_simple

    return m, lower, upper, sd_lower, sd_upper


def conflevel(Ax, npr, dt, rho, wr, nfft, nperseg):
    """
    Estimate chi-square confidence levels based on red-noise model.

    Inputs:
        Ax       : float – mean spectral density
        npr      : int   – number of frequencies
        dt       : float – sampling interval
        rho      : float – lag-1 autocorrelation coefficient
        wr       : np.ndarray – frequency array
        nfft     : int   – FFT size
        nperseg  : int   – segment length

    Output:
        tabtchi : np.ndarray – chi-square confidence levels
    """
    facchi95 = chi2inv(0.95, nfft, nperseg)
    fnyq = 1 / (2 * dt)  # Nyquist frequency

    cos_term = np.cos(np.pi * wr / fnyq)
    theored = (1 - rho**2) / (1 - 2 * rho * cos_term + rho**2)
    theored[0] = 0.0  # zero-frequency fix

    Art = np.nanmean(theored[1:])
    theored *= Ax / Art
    theored[0] = np.nanmean(theored)

    return theored * facchi95


def rhoAR1(datax):
    """
    Compute lag-1 autocorrelation coefficient for AR(1) process.

    Input:
        datax : np.ndarray – 1D data array

    Output:
        rho : float – lag-1 autocorrelation coefficient
    """
    datax = np.asarray(datax, dtype=float)
    datax -= np.nanmean(datax)
    num = np.nansum(datax[1:] * datax[:-1])
    den = np.nansum(datax[:-1] ** 2)
    return 0.0 if den == 0 else num / den


def RedConf(datax, dt, nsim, nperseg):
    """
    Compute red-noise confidence spectrum for time series.

    Inputs:
        datax   : np.ndarray – 1D data array
        dt      : float – sampling interval
        nsim    : int – number of red-noise simulations
        nperseg : int – segment length for Welch method

    Outputs:
        wr, tabtchi : np.ndarray – frequencies and confidence levels
    """
    rho = rhoAR1(datax)
    
    # w, Pxx = signal.welch(datax, fs=1/dt, nperseg=nperseg)
    # Ax = np.nanmean(Pxx)
    # npr = len(w)
    # tabtchi = conflevel(Ax, npr, dt, rho, w, nfft, nperseg)
    
    
    nt = len(datax)
    redtab = Rednoise(nt, rho, nsim)

    datan = datax - np.nanmean(datax)
    nfft = nperseg
    w, po = signal.welch(datan, fs=1 / dt, nperseg=nperseg)
    Ax = np.nanmean(po)

    wr, _ = signal.welch(redtab[:, 0] - np.nanmean(redtab[:, 0]), fs=1 / dt, nperseg=nperseg)
    npr = len(wr)
    tabtchi = conflevel(Ax, npr, dt, rho, wr, nfft, nperseg)
    
    return wr, tabtchi
    

def Rednoise(nt, rho, nsim):
    """
    Monte Carlo simulation of AR(1) red-noise series.

    Inputs:
        nt   : int – number of time steps
        rho  : float – lag-1 autocorrelation
        nsim : int – number of simulations

    Output:
        redtab : np.ndarray – red noise realizations (nt x nsim)
    """
    srho = np.sqrt(1 - rho**2)
    redtab = np.zeros((nt, nsim))
    white = np.random.randn(nt, nsim)
    for i in range(1, nt):
        redtab[i] = rho * redtab[i - 1] + srho * white[i]
    return redtab

    
def ciout(x):
    """
    Compute confidence interval (95%) width around mean.

    Input:
        x : np.ndarray – 1D data array

    Output:
        nivel : float – half-width of confidence interval
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return 0.0
    ci95 = t.interval(0.95, len(x) - 1, loc=np.mean(x), scale=sem(x))
    return ci95[1] - np.mean(x)

def average(arr, n):
    """
    Compute n-point averaged data.

    Inputs:
        arr : np.ndarray – input array
        n   : int – number of samples per average (0 = average all)

    Outputs:
        np.ndarray or float – averaged data
    """
    n = int(n)
    if n == 0:
        return np.nanmean(arr)
    end = (len(arr) // n) * n
    return np.nanmean(arr[:end].reshape(-1, n), axis=1)

def ci(x):
    """
    Compute 95% confidence interval of the mean.

    Inputs:
        x : np.ndarray – input data (may contain NaN)

    Outputs:
        lower, upper : float – lower and upper 95% confidence limits
    """
    x = np.asarray(x, dtype=float)
    valid = np.isfinite(x)
    size = np.count_nonzero(valid)
    if size == 0:
        return np.nan, np.nan

    sdev = np.nanstd(x)
    dev = 1.96 * sdev / np.sqrt(size)  # 95% confidence level
    mean = np.nanmean(x)
    return mean - dev, mean + dev



def velocityten(wz, z):
    """
    Convert wind velocity at z meters to equivalent velocity at 10 meters.

    Inputs:
        wz : np.ndarray – wind velocity at z meters (m/s)
        z  : float – measurement height (m)

    Outputs:
        w10 : np.ndarray – wind velocity scaled to 10 m (m/s)

    Method:
        Neutral logarithmic profile,

            U10 = Uz / ( 1 - sqrt(C_D)/kappa * ln(10/z) ) ,

        with kappa = 0.4.  The drag coefficient follows the law selected by
        :func:`drag_coefficient`, so the correction is continuous in wind
        speed; version 1.x used the discontinuous step law here as well.
    """
    if z == 10:
        return np.asarray(wz, dtype=float)

    wz = np.asarray(wz, dtype=float)
    kappa = 0.4
    offset = math.log(10.0 / z)

    Cd = drag_coefficient(wz)
    denominator = 1.0 - np.sqrt(Cd) / kappa * offset

    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(np.abs(denominator) > 1e-6, wz / denominator, wz)


def wind_average(wd, iw):
    """
    Compute intensity-weighted average wind direction.

    Inputs:
        wd : np.ndarray – wind direction (degrees)
        iw : np.ndarray – wind intensity (m/s)

    Outputs:
        mean : float – mean wind direction (degrees, 0–360)
    """
    wd_rad = np.deg2rad(wd)
    u_east = np.average(iw * np.sin(wd_rad))
    u_north = np.average(iw * np.cos(wd_rad))
    mean = np.degrees(np.arctan2(u_east, u_north)) % 360
    return mean


   

def windbar_fetch(dw_mean, linang, angle, dists):
    """
    Compute wind fetch contribution for given direction range.

    Inputs:
        dw_mean : float – mean wind direction (degrees)
        linang  : float – angular threshold for contribution (degrees)
        angle   : np.ndarray – direction array (degrees)
        dists   : np.ndarray – distance array (m)

    Outputs:
        Ldist : np.ndarray – distances where contribution occurs, NaN otherwise
    """
    dw_min, dw_max = wind_angle(dw_mean, linang)
    angle = np.asarray(angle, dtype=float)
    dists = np.asarray(dists, dtype=float)
    Ldist = np.full_like(angle, np.nan)

    if dw_max < dw_min:
        mask = (angle < dw_max) | (angle > dw_min)
    else:
        mask = (angle > dw_min) & (angle < dw_max)

    Ldist[mask] = dists[mask]
    return Ldist


def wind_angle(wind, linang):
    """
    Compute angular boundaries of wind direction contribution.

    Inputs:
        wind   : float – wind direction (degrees)
        linang : float – angular threshold (degrees)

    Outputs:
        dw_min, dw_max : float – lower and upper angular limits (degrees)
    """
    dw_max = (wind + linang) % 360
    dw_min = (wind - linang) % 360
    return dw_min, dw_max


#  Selected drag law; set once by the backend from the additional parameters
DRAG_LAW = 'smooth'


def drag_coefficient(w, law=None):
    """
    Neutral drag coefficient of the water surface at 10 m.

    Inputs:
        w   : float or np.ndarray – wind speed at 10 m (m/s)
        law : str or None – 'smooth' (default), 'legacy' or 'largepond'

    Outputs:
        float or np.ndarray – drag coefficient C_D10 (dimensionless)

    Method:
        'legacy'    reproduces version 1.x exactly: C_D = 1.0e-3 below 5 m/s
                    and 1.5e-3 above.  That step makes the surface stress jump
                    by 50 % as the wind crosses 5 m/s, and the jump propagates
                    into the friction velocity, the Richardson number, the
                    Wedderburn number and every classification derived from
                    them, producing steps in the time series that have no
                    physical counterpart.

        'smooth'    keeps the same two end values but joins them linearly
                    between 4 and 6 m/s,

                        C_D = 1.0e-3                       U <= 4
                        C_D = 1.0e-3 + 0.25e-3 (U - 4)     4 < U < 6
                        C_D = 1.5e-3                       U >= 6 ,

                    so the stress is continuous while the calibration of the
                    original formulation is preserved.  This is the default.

        'largepond' applies Large and Pond (1981),

                        C_D = 1.2e-3                       U <= 11
                        C_D = (0.49 + 0.065 U) 1e-3        U > 11 ,

                    which is the standard open-water parametrisation and is
                    preferable for large basins exposed to strong winds.
    """
    law = (law or DRAG_LAW).lower()
    w = np.abs(np.asarray(w, dtype=float))

    if law == 'legacy':
        return np.where(w > 5.0, 0.0015, 0.0010)

    if law == 'largepond':
        return np.where(w <= 11.0, 0.0012, (0.49 + 0.065 * w) * 1.0e-3)

    return np.clip(1.0e-3 + 0.25e-3 * (w - 4.0), 1.0e-3, 1.5e-3)


def wind_stress(w, law=None):
    """
    Surface wind stress.

    Inputs:
        w   : float or np.ndarray – wind speed at 10 m (m/s)
        law : str or None – drag law, see :func:`drag_coefficient`

    Outputs:
        float or np.ndarray – wind stress tau = rho_air C_D U^2 (N/m²)
    """
    w = np.maximum(np.abs(np.asarray(w, dtype=float)), 0.01)
    return drag_coefficient(w, law) * RHO_AIR * w ** 2

def wind_parameters(w, rw, pe, he, n, glin, H):
    """
    Compute 1D wind-related stability parameters.

    Inputs:
        w     : float – wind velocity (m/s)
        rw    : (unused) float – reserved for future use
        pe    : float – epilimnion density (kg/m³)
        he    : float – epilimnion thickness (m)
        n     : (unused) int – reserved for future use
        glin  : float – buoyancy gradient (s⁻²)
        H     : float – total depth (m)

    Outputs:
        stress : float – wind stress (N/m²)
        wast   : float – friction velocity in the water (m/s)
        riw    : float – bulk Richardson number of the surface layer
    """
    stress = float(wind_stress(w))

    if not np.isfinite(pe) or pe <= 0:
        return stress, np.nan, np.nan

    wast = math.sqrt(stress / pe)
    riw = glin * he / (wast ** 2) if wast > 0 else np.nan

    return stress, wast, riw

def richardson(w, rw, qt, h, pe, H, n2d, hmid, p, glin):
    """
    Compute depth-resolved Richardson number (2D).

    Inputs:
        w     : float – wind velocity (m/s)
        rw    : (unused) float – reserved for future use
        qt    : int – number of vertical layers
        h     : np.ndarray – depth array (m)
        pe    : float – air density (kg/m³)
        H     : float – total water depth (m)
        n2d   : np.ndarray – buoyancy frequency profile (s⁻²)
        hmid  : np.ndarray – mid-depths between layers (m)
        p     : np.ndarray – density profile (kg/m³)
        glin  : np.ndarray – buoyancy gradient profile (s⁻²)

    Outputs:
        riw2d : np.ndarray – Richardson number profile (along z-direction)

    Notes:
        The friction velocity is referenced to the mean density of each layer
        pair, and the length scale is the distance of the layer below the water
        surface, so that ``riw2d`` is the bulk Richardson number that would be
        obtained by treating every measurement level in turn as the base of the
        wind-mixed layer.
    """
    win_stress = float(wind_stress(w))

    rho_mean = np.abs(0.5 * (p[:-1] + p[1:]))
    wast = np.sqrt(win_stress / rho_mean)

    with np.errstate(divide='ignore', invalid='ignore'):
        riw2d = glin[:qt - 1] * (H - hmid[:qt - 1]) / (wast ** 2)

    return riw2d


def sorting_2d(data):
    """
    Sort each 1D slice of a 2D array in descending order
    while preserving NaN structure.
    
    Inputs:
        data : np.ndarray – 2D array (shape: [time, depth])
               May contain NaN values.
    
    Outputs:
        ordered : np.ndarray – array with same shape as data
                  Each row is sorted in descending order starting
                  from the first non-NaN index, with NaNs preserved
                  in their original positions.
    """
    ordered = np.full_like(data, np.nan)
    for i, row in enumerate(data):
        valid = ~np.isnan(row)
        if not np.any(valid):  # skip fully NaN rows
            continue
        first_idx = np.argmax(valid)
        sorted_vals = np.sort(row[valid])[::-1]  # descending
        ordered[i, first_idx:first_idx + len(sorted_vals)] = sorted_vals
    return ordered

def thorpe_scale(temp, h, tol=1e-12):
    """
    Compute Thorpe displacements and Thorpe scale.
    
    Inputs:
        temp : np.ndarray – temperature field (nt, nz)
        h    : np.ndarray – height above bottom (m, shape: nt x nz)
    
    Outputs:
        d  : np.ndarray – Thorpe displacements (m, shape: nt x nz)
        Lt : np.ndarray – Thorpe scale per time step             (RMS of overturning displacements, size nt, m)
    """
    temp = np.asarray(temp)
    h    = np.asarray(h)
    
    nt, nz = temp.shape
    d  = np.full((nt, nz), np.nan)
    Lt = np.full(nt, np.nan)
    
    for t in range(nt):
        prof = temp[t, :]
        hs   = h[t, :]
        
        # Skip if missing data
        if np.isnan(prof).any() or np.isnan(hs).any():
            continue
        
        # Check if column is already stable
        if np.all(np.diff(prof) <= tol):  
            
            # already monotonic decreasing with z (stable)
            d[t, :]  = 0.0
            Lt[t] = 0.0
            continue
        
        # Sort temperatures in stable order:
        order = np.argsort(prof)[::-1]   # descending temperature
        
        # inverse permutation
        inv = np.empty(nz, dtype=int)
        inv[order] = np.arange(nz)
        
        # Thorpe displacement:
        d_t = hs[inv] - hs
        d_t[np.abs(d_t) < tol] = 0.0
        
        d[t, :] = d_t
        
        mask = np.abs(d_t) > tol
        Lt[t] = np.sqrt(np.mean(d_t[mask]**2)) if np.any(mask) else 0.0

    return d, Lt

def sorting_1d(data):
    """
    Sort a 1D array in descending order while preserving NaN structure.
    
    Inputs:
        data : np.ndarray (1D)
            Array with possible NaN values.
            
    Outputs:
        auxiliary : np.ndarray (1D)
            Sorted array (descending) starting from the first non-NaN index.
    """
    valid = ~np.isnan(data)
    if not np.any(valid):
        return np.full_like(data, np.nan)
    first_idx = np.argmax(valid)
    sorted_vals = np.sort(data[valid])[::-1]
    auxiliary = np.full_like(data, np.nan)
    auxiliary[first_idx:first_idx + len(sorted_vals)] = sorted_vals
    return auxiliary



def find_nearest(array, value):
    """
    Find the nearest value and its index in an array, ignoring NaNs.
    
    Inputs:
        array : np.ndarray – 1D target array
        value : float – value to search for
    
    Outputs:
        nearest_value : float – closest element in array
                        (np.nan if all elements are NaN)
        idx           : int – index of nearest element
                        (-1 if all elements are NaN)
    """
    array = np.asarray(array)
    valid = ~np.isnan(array)
    if not np.any(valid):  # all NaN
        return np.nan, -1

    valid_array = array[valid]
    idx_valid = np.abs(valid_array - value).argmin()
    
    # Map back to original index
    idx = np.arange(len(array))[valid][idx_valid]
    return array[idx], idx
   
def class_generation(riw, hh, he, ls):
    """
    Classify lake mixing/internal wave regime.
    
    Inputs:
        riw : float – Wedderburn number or Richardson-like parameter
        hh  : float – hypolimnion thickness (m)
        he  : float – epilimnion thickness (m)
        ls  : float – basin length (m)
    
    Outputs:
        regime : int – mixing/internal wave regime classification
                 1 = stratification breakdown by mixing
                 2 = large internal displacement with billows/mixing
                 3 = dominant internal seiche
                 4 = small-amplitude internal seiche
    """
    aux1  = math.sqrt((hh + he) / he)
    aux2  = ls / (2 * he)
    
    cond1 = 1   
    cond2 = aux2 * aux1
    cond3 = aux1 * aux2**2
    
    if riw < 1:
        return 1, cond1, cond2, cond3
    
    if riw < aux2 * aux1:
        return 2, cond1, cond2, cond3
    if riw < aux1 * aux2**2:
        return 3, cond1, cond2, cond3
    return 4, cond1, cond2, cond3

def spigel(he, H, W1, W2):
    """
    Compute BSIW amplitude based on Spigel & Imberger (1980) theory.
    
    Inputs:
        he : float – epilimnion thickness (m)
        H  : float – total depth (m)
        W1 : float – basin width at the interface level (m)
        W2 : float – basin width at the surface (m)
    
    Outputs:
        zeta : float – estimated internal seiche amplitude (m)
    """
    zeta = 0.4996 * he / W1
    if zeta > H:
        zeta = 0.4996 * he / W2
    return zeta

def bueno(he, hh, W1, W2, typ=None):
    """
    Compute BSIW amplitude based on de Carvalho Bueno et al. (2020).
    
    Inputs:
        he  : float – epilimnion thickness (m)
        hh  : float – hypolimnion thickness (m)
        W1  : float or np.ndarray – basin width or scaling factor (interface level)
        W2  : float or np.ndarray – basin width or scaling factor (surface level)
        typ : str, optional – if 'amplitude', return single value;
                              otherwise return array result
    
    Outputs:
        amp : float or np.ndarray – estimated amplitude (m)
              following exponential model formulation
    """
    xi, k1, k2 = 0.1, 6, 0.25
    H = he / (he + hh)
    g = 12.156 * H**3 - 15.714 * H**2 + 2.8426 * H + 2.0846
    f = g * np.exp((H**2) / k2)

    if typ == 'amplitude':
        amp = xi * he * np.exp(((W1 - k1)**2) / (2 * f**2))
        amp_limit = 0.4996 * he / W1
        if amp_limit > he + hh:
            amp_limit = 0.4996 * he / W2
        return amp

    W1 = np.asarray(W1)
    return xi * he * np.exp(((W1 - k1)**2) / (2 * f**2))

def weddFilters(wedd, ver_wind, ver_dire, dt, tbsiw, dw_hom, v1mode, lin):
    """
    Apply duration and direction filters to WEDD index.
    
    Inputs:
        wedd     : np.ndarray – original W index (size lin)
        ver_wind : np.ndarray – counter of consecutive wind events
        ver_dire : np.ndarray – counter of consecutive direction-aligned wind events
        dt       : float – time step (hours)
        tbsiw    : float – internal-wave basin period (s)
        dw_hom   : np.ndarray – direction homogeneity flag array
        v1mode   : np.ndarray – mode-1 velocity scale array
        lin      : int – number of time steps
    
    Outputs:
        wedu  : np.ndarray – duration-filtered W index
        wedi  : np.ndarray – direction-filtered W index
        fdura : float – duration filter factor
        fdire : np.ndarray – direction filter factor (vector)
    """

    wedu = np.zeros(lin, dtype=float)
    wedi = np.zeros(lin, dtype=float)

    # Convert hours to seconds
    dura = ver_wind * dt * 3600.0
    dire = ver_dire * dt * 3600.0

    # Duration filter (fdura)
    if np.isnan(tbsiw) or tbsiw == 0:
        fdura = 1.0
    else:
        fdura = min(1.0, np.sqrt(4.0 * dura / tbsiw))
        if fdura <= 0 or np.isnan(fdura):
            fdura = 1.0

    # Direction filter (fdire)
    fdire = np.ones(lin, dtype=float)
    twind = dt * 3600.0  # accumulate wind duration

    for t in range(lin):

        if not np.isnan(dw_hom[t]) and dw_hom[t] == 357:

            if not np.isnan(v1mode[t]) and v1mode[t] != 0:
                fdire[t] = min(1.0, 4.0 * twind / v1mode[t])
            else:
                fdire[t] = 1.0

            twind += dt * 3600.0

        else:
            fdire[t] = 1.0
            twind = dt * 3600.0

    # Apply filters
    for t in range(lin):

        # Duration-filtered W
        if fdura != 0:
            try:
                wedu[t] = wedd[t] / fdura
            except Exception:
                wedu[t] = wedd[t]
        else:
            wedu[t] = wedd[t]

        # Direction-filtered W
        if fdire[t] != 0:
            wedi[t] = wedd[t] / fdire[t]
        else:
            wedi[t] = wedd[t]

    return wedu, wedi, fdura, fdire, dura, dire    

def iw_generation(wedd, hh, he, ls):
    """
    Identify regime 3 (internal seiche-dominant) generation limits.
    
    Inputs:
        wedd : float – Wedderburn number
        hh   : float – hypolimnion thickness (m)
        he   : float – epilimnion thickness (m)
        ls   : float – basin length (m)
    
    Outputs:
        lower_gene : float – lower Wedderburn limit for regime 3
        upper_gene : float – upper Wedderburn limit for regime 3
    """
    aux1 = math.sqrt((hh + he) / hh)
    lower_gene = 0.5 * aux1
    upper_gene = (ls / (4 * he)) * aux1**2
    return lower_gene, upper_gene


def period_analysis(periodSignal):
    """
    Analyzing time-series of wave periods for different vertical modes.

    Inputs:
        serie : np.ndarray – input time series

    """

    barT = np.nanmean(periodSignal)
    stdT = np.std(periodSignal)
    conT = 1.96 * np.nanstd(periodSignal, ddof=1) / np.sqrt(np.sum(~np.isnan(periodSignal)))
    
    return np.array([barT-stdT,barT,barT+stdT,conT])

def mask(serie):
    """
    Interpolate and replace NaN values in a 1D time series.

    Inputs:
        serie : np.ndarray – input time series (may contain NaNs)
    
    Outputs:
        xfiltered : np.ndarray – interpolated series (no NaNs)
    """
    serie = np.asarray(serie)
    m = np.isfinite(serie)
    if not np.any(m):  # all NaNs
        return np.zeros_like(serie)
    xi = np.arange(len(serie))
    return np.interp(xi, xi[m], serie[m])

def findLayers(vel_aux, refined_depth, H,modes=(0,1,2,3), eps=1e-4):
    """
    Identify zero-crossing depths for vertical modes of internal seiches.
        
    Outputs:
        crossings          : list – crossings[i] contains zero-crossing depths
                             (m above bed) for mode i
        crossingsThickness : list – crossingsThickness[i] contains layer
                             thicknesses for mode i
                             (length = expected[i] + 1)
    """

    # expected number of crossings per mode
    expected = {0: 1, 1: 2, 2: 3, 3: 4}

    crossings = {}
    crossingsThickness = {}

    for i in modes:
        v = vel_aux[:, i]

        # Treat very small values as zero
        v = np.where(np.abs(v) < eps, 0, v)

        # Locate sign changes
        sign_change = np.where(np.diff(np.sign(v)) != 0)[0]

        found = []

        for idx in sign_change:
            # Depths surrounding zero crossing
            z1, z2 = refined_depth[idx], refined_depth[idx+1]
            v1, v2 = v[idx], v[idx+1]

            # Linear interpolation for zero crossing
            if (v2 - v1) != 0:
                z_zero = z1 - v1 * (z2 - z1) / (v2 - v1)
                found.append(z_zero)

        # FORCE FIXED LENGTH

        n = expected[i]

        if len(found) < n:
            # Pad missing values with NaN
            found = found + [np.nan] * (n - len(found))

        elif len(found) > n:
            # Trim extra values
            found = found[:n]

        crossings[i] = found


        depths_from_surface = [
            H - x if not np.isnan(x) else np.nan
            for x in found
        ]

        # Sort so that smallest depth = shallowest layer boundary
        depths_from_surface_sorted = sorted(
            depths_from_surface,
            key=lambda x: (np.inf if np.isnan(x) else x)
        )

        # Build final boundaries from surface to bottom
        # Example Mode 1: [0, h_surface1, h_surface2, H]
        boundaries = [0] + depths_from_surface_sorted + [H]

        # Compute layer thickness from surface downward
        thickness = []
        for j in range(len(boundaries)-1):
            b1 = boundaries[j]
            b2 = boundaries[j+1]

            if np.isnan(b1) or np.isnan(b2):
                thickness.append(np.nan)
            else:
                thickness.append(b2 - b1)

        crossingsThickness[i] = thickness

    return crossings, crossingsThickness


def welch_method(serie, size, w, dt):
    """
    Compute power spectral density (PSD) using Welch’s method.

    Inputs:
        serie : np.ndarray – time-series data
        size  : float – window size (s)
        w     : str – window type (e.g., 'hann')
        dt    : float – time step (h)
    
    Outputs:
        freq : np.ndarray – frequency (Hz)
        per  : np.ndarray – period (h)
        ff   : np.ndarray – power spectral density (x²/Hz)
        wr   : np.ndarray – red-noise spectrum
        conf : np.ndarray – 95% confidence interval
    """
    serie = mask(serie)
    serie -= np.nanmean(serie)
    dt_sec = 3600 * dt
    n = int((size or 4 * 24 * 3600) / dt_sec)
    n = min(n, len(serie))

    freq, ff = signal.welch(serie, fs=1.0 / dt_sec, window=w, nperseg=n, detrend='linear')
    wr, conf = RedConf(serie, dt_sec, 10, n)  # Monte Carlo (10 sim)
    per = 1 / freq[1:] / 3600

    return freq[1:], per, ff[1:], wr, conf


def wave_spectral(si, dt, mother):
    """
    Compute wavelet transform of a time series.

    Inputs:
        si      : np.ndarray – time-series data
        dt      : float – time step (h)
        mother  : str – mother wavelet type
    
    Outputs:
        time  : np.ndarray – time axis (h)
        per   : np.ndarray – wavelet period (h)
        power : np.ndarray – wavelet power (x²)
    """
    si = mask(si)
    dt = float(dt)
    pad, dj, s0, j1 = 1, 0.25, 4 * dt, int(7 / 0.25)

    time = np.arange(len(si)) * dt
    data = (si - np.nanmean(si)) / np.nanstd(si)
    warnings.filterwarnings("ignore")

    wav, per, sca, coi = wavelet(data, dt, pad, dj, s0, j1, mother)
    return time, per, np.abs(wav) ** 2

def coherence_shift(s1, s2, window, dt):
    """
    Compute coherence and phase shift between two signals.

    Inputs:
        s1, s2 : np.ndarray – input signals
        window : float – window size (s)
        dt     : float – time step (h)
    
    Outputs:
        phase   : np.ndarray – phase shift (deg)
        cxy     : np.ndarray – magnitude-squared coherence
        fcoh    : np.ndarray – frequency (Hz)
        conf95  : np.ndarray – 95% confidence indices
    """
    dt_sec = 3600 * dt
    nperseg = int(window / dt_sec)

    fcoh, cxy = signal.coherence(s1, s2, fs=1 / dt_sec, nperseg=nperseg)
    _, pxy = signal.csd(s1, s2, fs=1 / dt_sec, nperseg=nperseg)

    phase = np.rad2deg(np.angle(pxy))
    phase[phase < -90] += 360

    edof = (len(s1) / (nperseg / 2)) * np.mean(cxy)
    gamma95 = 1 - (0.8) ** (1 / max(edof - 1, 1))
    conf95 = np.where(cxy > gamma95)

    return phase, cxy, fcoh, conf95

def butter_bandpass(lowcut, highcut, fs, order):
    """
    Compute coefficients of a Butterworth band-pass filter.

    Inputs:
        lowcut : float – lower cutoff frequency (Hz)
        highcut: float – upper cutoff frequency (Hz)
        fs     : float – sampling frequency (Hz)
        order  : int – filter order
    
    Outputs:
        sos : np.ndarray – second-order sections for filtering
    """
    nyq = 0.5 * fs
    sos = signal.butter(order, [lowcut / nyq, highcut / nyq],
                        btype='band', output='sos', analog=False)
    return sos

def butter_bandpass_filter(data, lowcut, highcut, fs):
    """
    Apply Butterworth band-pass filter to a signal.

    Inputs:
        data   : np.ndarray – signal to filter
        lowcut : float – lower cutoff frequency (Hz)
        highcut: float – upper cutoff frequency (Hz)
        fs     : float – sampling frequency (Hz)
    
    Outputs:
        f : np.ndarray – filtered signal
    """
    data = mask(data)
    data -= np.nanmean(data)
    sos = butter_bandpass(lowcut, highcut, fs, 1)
    return signal.sosfilt(sos, data)

def normalizeList(data):
    """
    Normalize a list of depth/freq/etc entries so that:
    - None = empty array
    - [None] = empty array
    - np.array([...]) = unchanged
    """

    normalized = []

    for item in data:
        # Case 1: item is None
        if item is None:
            normalized.append(np.array([]))
            continue

        # Case 2: item is [None]
        if isinstance(item, list) and len(item) == 1 and item[0] is None:
            normalized.append(np.array([]))
            continue

        # Case 3: item is list containing array(s)
        if isinstance(item, list):
            # assume first element is the array you want
            arr = item[0]
            if arr is None:
                normalized.append(np.array([]))
            else:
                normalized.append(np.asarray(arr))
            continue

        # Case 4: item is already an array
        if isinstance(item, np.ndarray):
            normalized.append(item)
            continue

        # Fallback
        normalized.append(np.asarray(item))

    return normalized

import numpy as np

def resampleSignalMatrix(values, dx_mod, dx):
    """
    Resample time-series array from dx_mod to dx.

    """

    values = np.asarray(values)

    # Convert datetimes to seconds
    t0 = dx_mod[0]
    t_old = np.array([(d - t0).total_seconds() for d in dx_mod])
    t_new = np.array([(d - t0).total_seconds() for d in dx])

    # 1D case 
    if values.ndim == 1:
        return np.interp(t_new, t_old, values)

    # Multi-dimensional case
    original_shape = values.shape
    values_flat = values.reshape(values.shape[0], -1)

    values_new_flat = np.array([
        np.interp(t_new, t_old, values_flat[:, i])
        for i in range(values_flat.shape[1])
    ]).T

    new_shape = (len(dx),) + original_shape[1:]
    return values_new_flat.reshape(new_shape)