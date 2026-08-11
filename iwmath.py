# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer - Modelling functions

Interwave Analyzer - Version 2 (2026) 
Modelling functions module version: 2.260305

-------------------------------------------------------------------------------

de Carvalho Bueno, R; Bleninger, T. B.; Lorke, A. 
Internal wave analyzer for thermally stratified lakess 
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
import numpy as np
import math as ma

from iwmod import thermal_stability


# Number of vertical modes returned by the decomposition model
N_MODES = 5

# Vertical resolution of the refined grid used by the decomposition model
N_LEVELS = 100

# Lower bound imposed on N^2 (s^-2).  A strictly positive value is required to
# keep the generalised eigenvalue problem well posed; 1e-8 s^-2 corresponds to
# a buoyancy period of about 17 hours, i.e. an essentially unstratified layer.
N2_FLOOR = 1.0e-8


# =============================================================================
#  Vertical mode decomposition
# =============================================================================

try:                                        # optional fast path
    from scipy.linalg import eigh_tridiagonal as _eigh_tridiagonal
except ImportError:                          # pragma: no cover
    _eigh_tridiagonal = None


def _smallest_eigenpairs(diagonal, offdiagonal, n_modes):
    """
    Smallest ``n_modes`` eigenpairs of a symmetric tridiagonal matrix.

    Inputs:
        diagonal    : np.ndarray - main diagonal
        offdiagonal : np.ndarray - first off-diagonal
        n_modes     : int        - number of eigenpairs required

    Outputs:
        (eigenvalues, eigenvectors) in ascending order, or ``(None, None)``
        when the decomposition fails.

    Notes:
        ``scipy.linalg.eigh_tridiagonal`` exploits the tridiagonal structure
        and can be restricted to the required subset of the spectrum, which is
        appreciably faster than forming the dense matrix.  The dense
        ``numpy.linalg.eigh`` is kept as a fallback so that the module remains
        usable without SciPy.
    """
    n = diagonal.size

    if _eigh_tridiagonal is not None:
        try:
            values, vectors = _eigh_tridiagonal(
                diagonal, offdiagonal, select='i',
                select_range=(0, min(n_modes, n) - 1))
            return values, vectors
        except Exception:
            pass

    matrix = np.zeros((n, n), dtype=float)
    idx = np.arange(n)
    matrix[idx, idx] = diagonal
    matrix[idx[:-1], idx[1:]] = offdiagonal
    matrix[idx[1:], idx[:-1]] = offdiagonal

    try:
        values, vectors = np.linalg.eigh(matrix)
    except np.linalg.LinAlgError:
        return None, None

    return values[:n_modes], vectors[:, :n_modes]


def buoyancy_profile(tau, h, H, n_levels=N_LEVELS, n2_floor=N2_FLOOR):
    """
    Build the squared buoyancy frequency on a uniform vertical grid.

    Inputs:
        tau      : np.ndarray - temperature profile (degC), ordered from the
                   surface downwards
        h        : np.ndarray - depth of every sensor below the water surface
                   (m), increasing
        H        : float      - total water depth (m)
        n_levels : int        - number of points of the refined grid
        n2_floor : float      - lower bound imposed on N^2 (s^-2)

    Outputs:
        (depth, n2)
            depth : np.ndarray - refined grid, 0 (surface) to H (bed), m
            n2    : np.ndarray - squared buoyancy frequency on that grid, s^-2

    Notes:
        The stratification is piecewise constant between measurement levels,

            N^2_j = (g / rho_{j+1}) * (rho_{j+1} - rho_j) / (h_{j+1} - h_j) ,

        with the shallowest value extended up to the surface and the deepest
        value extended down to the bed.  Version 1.x located each grid point by
        a Python ``while`` loop, and substituted the hard-coded values
        0.3 min^-1 or 0.08 min^-1 whenever the search fell outside the sensor
        range.  The mapping is now a single ``searchsorted`` call and the
        out-of-range substitution is replaced by an explicit extension of the
        outermost measured layers, so that no undocumented stratification is
        injected into the model.
    """
    tau = np.asarray(tau, dtype=float)
    h = np.asarray(h, dtype=float)

    finite = np.isfinite(tau) & np.isfinite(h)
    tau, h = tau[finite], h[finite]

    depth = np.linspace(0.0, H, n_levels)

    if tau.size < 2:
        return depth, np.full(n_levels, n2_floor)

    # Density from the equation of state, then layer stability
    rho = np.polyval([6.536336e-9, -1.120083e-6, 1.001685e-4,
                      -9.09529e-3, 6.793952e-2, 999.842592], tau)

    dz = np.diff(h)
    dz = np.where(np.abs(dz) < 1e-6, 1e-6, dz)
    n2_layer = 9.81 * np.diff(rho) / (rho[1:] * dz)

    # Unstable layers (denser water above lighter water) carry no restoring
    # force for the internal wave problem; they are clipped to the floor.
    n2_layer = np.clip(n2_layer, n2_floor, None)

    # Piecewise-constant mapping onto the refined grid
    index = np.clip(np.searchsorted(h, depth, side='right') - 1,
                    0, n2_layer.size - 1)
    n2 = n2_layer[index]

    return depth, n2


def decomposition(L, tau, h, H, n_modes=N_MODES, n_levels=N_LEVELS):
    """
    Vertical mode decomposition of a stratified basin.

    Inputs:
        L        : float      - basin length at the wind direction (m)
        tau      : np.ndarray - temperature profile (degC), surface to bed
        h        : np.ndarray - depth of every sensor below the surface (m)
        H        : float      - total water depth (m)
        n_modes  : int        - number of vertical modes to return
        n_levels : int        - vertical resolution of the refined grid

    Outputs:
        vel   : np.ndarray - horizontal velocity structure, shape
                (n_levels - 1, n_modes), arbitrary units
        conv  : list        - N^2(z) * W_m(z) for every mode, each of length
                ``n_levels``
        depth : np.ndarray - refined grid, 0 (surface) to H (bed), m
        peri  : list[float] - period of every mode (hours)
        cond  : int         - 1 when the solution is not trustworthy

    Method:
        For a linear internal wave of horizontal wavenumber k in a basin of
        constant depth, the vertical velocity structure W(z) obeys

            d2W/dz2 + k^2 ( N^2(z)/omega^2 - 1 ) W = 0 ,
            W(0) = W(H) = 0 ,

        where the term -k^2 W carries the non-hydrostatic correction.  Writing
        lambda = 1/omega^2 and discretising the second derivative on a uniform
        grid turns the problem into the symmetric generalised eigenvalue
        problem

            ( -D2 + k^2 I ) W = lambda ( k^2 diag(N^2) ) W ,

        i.e. A W = lambda B W with A symmetric positive definite and B
        diagonal positive.  Substituting W = B^(-1/2) y gives the standard
        symmetric problem

            ( B^(-1/2) A B^(-1/2) ) y = lambda y ,

        which is solved in one call of ``numpy.linalg.eigh``.  The period of
        each mode follows from

            T_m = 2 pi sqrt(lambda_m) ,

        so the fundamental vertical mode V1H1 is the *smallest* eigenvalue and
        the modes come out already ordered by increasing period, which is the
        physical ordering for internal seiches: the modal phase speed of a
        continuously stratified basin decreases with the mode number, hence
        V2H1 is slower - and therefore longer - than V1H1.

        Version 1.x solved the same equation by shooting: the recurrence was
        marched from the surface for a trial period, the period was incremented
        by one hour until the bottom boundary residual changed sign, and the
        root was then bracketed by bisection - separately for every mode and
        for every time step.  That procedure costs O(n_modes * n_iterations *
        n_levels) interpreted operations per profile, is sensitive to the
        arbitrary one-hour search increment (a mode whose period falls between
        two increments can be skipped entirely, and two modes closer than one
        hour can be returned twice), and diverges when the recurrence
        overflows.  The eigenvalue formulation returns all modes at once, is
        solved by LAPACK, cannot skip or duplicate a mode because the
        eigenvalues are ordered by construction, and has no tuning parameter.
    """
    tau = np.asarray(tau, dtype=float)
    h = np.asarray(h, dtype=float)
    cond = 0

    H = float(H)
    if not np.isfinite(H) or H <= 0:
        H = float(np.nanmax(h)) if np.any(np.isfinite(h)) else 1.0

    depth, n2 = buoyancy_profile(tau, h, H, n_levels=n_levels)
    dz = depth[1] - depth[0]

    if not np.isfinite(L) or L <= 0:
        L = 1.0
        cond = 1
    k = np.pi / L                       # first horizontal mode

    # ---- assemble the discrete operator on the interior points -------------
    # Interior unknowns are depth[1:-1]; W vanishes at both boundaries.
    n_in = n_levels - 2
    if n_in < n_modes + 1 or np.all(n2 <= N2_FLOOR * 1.001):
        # Either the grid is too coarse to resolve the requested modes, or the
        # water column carries no stratification at all.
        vel = np.full((n_levels - 1, n_modes), np.nan)
        conv = [np.full(n_levels, np.nan) for _ in range(n_modes)]
        return vel, conv, depth, [np.nan] * n_modes, 1

    n2_in = n2[1:-1]

    # A = -D2 + k^2 I  (symmetric tridiagonal)
    main = 2.0 / dz ** 2 + k ** 2
    off = -1.0 / dz ** 2

    # B = k^2 diag(N^2); the similarity transform uses s = sqrt(diag(B))
    s = k * np.sqrt(n2_in)
    inv_s = 1.0 / s

    # C = S^-1 A S^-1, still symmetric tridiagonal
    c_main = main * inv_s * inv_s
    c_off = off * inv_s[:-1] * inv_s[1:]

    eigenvalues, eigenvectors = _smallest_eigenpairs(c_main, c_off, n_modes)

    if eigenvalues is None:
        vel = np.full((n_levels - 1, n_modes), np.nan)
        conv = [np.full(n_levels, np.nan) for _ in range(n_modes)]
        return vel, conv, depth, [np.nan] * n_modes, 1

    peri = []
    w = np.zeros((n_modes, n_levels), dtype=float)

    for m in range(n_modes):
        lam = eigenvalues[m]
        if not np.isfinite(lam) or lam <= 0:
            peri.append(np.nan)
            w[m, :] = np.nan
            cond = 1
            continue

        peri.append(2.0 * np.pi * np.sqrt(lam) / 3600.0)     # hours

        profile = eigenvectors[:, m] * inv_s                 # back-substitute
        # Energy normalisation, consistent with version 1.x: sum(N^2 W^2) = 1
        norm = np.sqrt(np.sum(n2_in * profile ** 2))
        if norm > 0:
            profile = profile / norm

        # Deterministic sign: the shallowest non-negligible displacement is
        # taken positive, so successive time steps are directly comparable.
        significant = np.flatnonzero(np.abs(profile) > 1e-12)
        if significant.size and profile[significant[0]] < 0:
            profile = -profile

        w[m, 1:-1] = profile

    # Horizontal velocity structure is proportional to dW/dz
    hor = np.diff(w, axis=1)
    vel = np.transpose(hor)

    conv = [n2 * w[m, :] for m in range(n_modes)]

    if not np.all(np.isfinite(peri)):
        cond = 1

    return vel, conv, depth, peri, cond


def decomposition_shooting(L, tau, h, H):
    """
    Legacy shooting solver of the vertical structure equation.

    Kept for reproducibility of results published with Interwave Analyzer 1.x
    and for verification of :func:`decomposition`.  It is not used by the
    analysis pipeline; see :func:`decomposition` for the method actually
    applied and for the reasons behind the change.

    Inputs and outputs are identical to :func:`decomposition`.
    """
    cond = 0
    N = N_LEVELS
    increment = 3600
    k = np.pi / L

    a = len(tau)
    tau = np.concatenate((tau, [-9999]))
    h = np.concatenate((h, [H]))
    H = h[-1]

    refined_depth = np.linspace(0, H, N)
    dz = refined_depth[1] - refined_depth[0]

    bv = np.zeros((N), float)
    nsq = np.zeros((N), float)

    rho, buoy, hmid, glin = thermal_stability(a, h, H, tau)

    buoy = np.concatenate((buoy, [np.nan]))
    buoy = buoy * 60  # = [1/min]

    if h[0] < 50:
        for i in range(N):
            ii = 1
            while h[ii] < refined_depth[i]:
                ii = ii + 1

            if buoy[ii - 1] > -1:
                bv[i] = buoy[ii - 1]
            else:
                bv[i] = 0.3 if h[ii - 1] < 50 else 0.08

            nsq[i] = bv[i] ** 2 / 3600
    else:
        for i in range(N):
            ii = 1
            try:
                while h[ii] > refined_depth[i]:
                    ii = ii + 1
            except IndexError:
                ii = ii - 1

            if buoy[ii - 1] > -1:
                bv[i] = buoy[ii - 1]
            else:
                bv[i] = 0.3 if h[ii - 1] < 50 else 0.08

            nsq[i] = bv[i] ** 2 / 3600

    W = np.ones((N), float)
    W[0] = 0
    W[1] = 1

    finnew = 0
    p = 0
    pnew = p
    peri = []
    conv = []

    w = np.zeros((5, N), float)
    n = np.zeros((5, N), float)
    hor = np.zeros((5, N - 1), float)

    for m in range(5):
        e = 0.5
        while e >= 0:
            for i in range(2, N):
                f = 2 - k ** 2 * dz ** 2 * (nsq[i - 1] * (p ** 2)
                                            / (2 * np.pi) ** 2 - 1)
                W[i] = -W[i - 2] + f * W[i - 1]

            finold = finnew
            finnew = W[-1]
            e = finold * finnew

            pold = pnew
            pnew = p
            p = p + increment

        if finnew > 0:
            randoben, randunten = pnew, pold
        else:
            randoben, randunten = pold, pnew

        finold = finnew

        while abs(randunten - randoben) > 0.1:
            p = 0.5 * (randunten + randoben)

            for i in range(2, N):
                f = 2 - k ** 2 * dz ** 2 * (nsq[i - 1] * (p ** 2)
                                            / (2 * np.pi) ** 2 - 1)
                W[i] = -W[i - 2] + f * W[i - 1]

            finnew = W[-1]

            if finnew < 0:
                randunten = p
            else:
                randoben = p

        normw = np.sqrt(sum((W * nsq) * np.transpose(W)))
        if not np.isfinite(normw) or normw == 0:
            normw = 1

        for i in range(N):
            w[m, i] = W[i] / normw
            n[m, i] = np.sqrt(nsq[i]) * W[i] / normw

        for i in range(N - 1):
            hor[m, i] = w[m, i + 1] - w[m, i]

        finnew = finold

        peri.append(p / 3600)
        p = p + increment

        conv.append((nsq) * w[m, :])

    vel = np.transpose(hor)

    return vel, conv, refined_depth, peri, cond

# 1D non-hydrostatic analytical model for two-layer system
def disp_zmodel (pe,ph,he,hh,L,m):
    
    gamma = pe/ph
    
    Lmin = np.nanmin(L)
    Lave = np.nanmean(L)
    Lmax = np.nanmax(L)
    
    peri_min = biquadratic(Lmin,he,hh,gamma,m)
    peri_ave = biquadratic(Lave,he,hh,gamma,m)
    peri_max = biquadratic(Lmax,he,hh,gamma,m)
       
    return peri_min, peri_ave, peri_max

# 1D hydrostatic analytical model for three-layer system
def disp_xmodel3(p1,p2,p3,h1,h2,h3,L,vertical,m):
 
    gamma12 = p1/p2
    gamma13 = p1/p3
    gamma23 = p2/p3
       
    A = [[h1, h1,  h1], [h2*gamma12,  h2, h2], [h3*gamma13,  h3*gamma23, h3]]
    
    solv = np.linalg.eigvals(A)
    
    Lmin = np.nanmin(L)
    Lave = np.nanmean(L)
    Lmax = np.nanmax(L)
   
    pv1_min, pv2_min  = eigen3_values(Lmin,solv[0],solv[1],m)
    pv1_ave, pv2_ave  = eigen3_values(Lave,solv[0],solv[1],m)
    pv1_max, pv2_max  = eigen3_values(Lmax,solv[0],solv[1],m)
    
    if(vertical==1):
        return pv1_min,pv1_ave,pv1_max
    else:   
        return pv2_min,pv2_ave,pv2_max



def coriolis_effect(To, lat):
    """
    Estimate the internal seiche period including Coriolis effect.
    
    Parameters
    ----------
    To : float
        Wave period without Coriolis effect [seconds].
    lat : float
        Latitude [degrees].
    
    Returns
    -------
    Tcor : float
        Estimated period including Coriolis effect [seconds].
    """
    Omega = 7.292115e-5  # Earth's rotation rate [s^-1]
    phi = np.radians(lat)
    f = 2 * Omega * np.sin(phi)
    omega0 = 2 * np.pi / To
    omega = np.sqrt(omega0**2 + f**2)
    Tcor = 2 * np.pi / omega
    return Tcor

def biquadratic(L,he,hh,gamma,m):
#
#   Solver Function: Solver for 2-layer model
#
    g = 9.81
    k = m*np.pi/L
    
    th   = np.tanh(k*hh)
    te   = np.tanh(k*he)

    p     = (gamma*te*th + 1)/(k*th)
    q     = -g*(th+te)/th
    r     = -ma.pow(g,2)*(gamma-1)*k*te

    delta = ma.pow(q,2) - 4*p*r
    
    try:
        omega = np.sqrt((-q-np.sqrt(delta))/(2*p))
    except RuntimeWarning:
        return None
    
    peri  = 2*np.pi/omega
    
    return peri 


def eigen3_values(L,lambv1,lambv2,m):
#
#   Solver Function: Solver for 3-layer model
#    
    g     = 9.81   # m/s²

    try:
        peri_v1 =  2*L/(m*np.sqrt(g*lambv1))     # V1 interfacial period (sec) 
    except RuntimeWarning:
        peri_v1 = None
        
    try:
        peri_v2 =  2*L/(m*np.sqrt(g*lambv2))     # V2 interfacial period (sec)
    except RuntimeWarning:
        peri_v2 = None
        
    return peri_v1, peri_v2


def eigen4_values(L,lambv1,lambv2,lambv3,m):
#
#   Solver Function: Solver for 4-layer model
#   
    g     = 9.81   # m/s²

    try:
        peri_v1 =  2*L/(m*np.sqrt(g*lambv1))     # V1 interfacial period (sec) 
    except RuntimeWarning:
        peri_v1 = None
    try:
        peri_v2 =  2*L/(m*np.sqrt(g*lambv2))     # V2 interfacial period (sec)
    except RuntimeWarning:
        peri_v2 = None
    try:
        peri_v3 =  2*L/(m*np.sqrt(g*lambv3))     # V3 interfacial period (sec)
    except RuntimeWarning:
        peri_v3 = None
        
    return peri_v1, peri_v2, peri_v3



def sensitivity_2layer(mean,diff,N,pe,ph,he,hh,fetch,typ):
#
#   Sensitivity Function: Sensitivity analysis to check period variation 
#   based on layer thickness and water density changes (two-layer system)
#      
    x   = np.linspace(mean-diff, mean+diff, N)
    period = np.zeros((N),float)

    for i in range(N):
        
        if   typ == 1:
            _,per,_  = np.array(disp_zmodel(x[i], ph, he, hh,fetch,1))
        elif typ == 2:
            _,per,_  = np.array(disp_zmodel(pe, x[i], he, hh,fetch,1))
        elif typ == 3:
            _,per,_  = np.array(disp_zmodel(pe, ph, x[i], hh,fetch,1))
        elif typ == 4:
            _,per,_  = np.array(disp_zmodel(pe, ph, he, x[i],fetch,1))            
            
        period[i] = per     
        
    return x, period/60/60  

#  Sensitivity analysis to check period variation
def sensitivity_3layer(mean,diff,N,p1,p2,p3,h1,h2,h3,fetch,typ):
     
    aux = mean - diff
    x   = np.linspace(aux, mean+diff, N)
    period = np.zeros((N),float)

    for i in range(N):

        if   typ == 1:
            _,per,_  = np.array(disp_xmodel3(x[i],p2,p3,h1,h2,h3,fetch,2,1))
        elif typ == 2:
            _,per,_  = np.array(disp_xmodel3(p1,x[i],p3,h1,h2,h3,fetch,2,1))
        elif typ == 3:
            _,per,_  = np.array(disp_xmodel3(p1,p2,p3,x[i],h2,h3,fetch,2,1))
        elif typ == 4:
            _,per,_  = np.array(disp_xmodel3(p1,p2,p3,h1,x[i],h3,fetch,2,1))
                
        period[i] = per     
        
    return x, period/60/60  

# Compute sensitivity curves for BSIW period varying
def sensitivity_dimension(L, pe, ph, he, hh):

    N = 50
    g = 9.81

    mean_rho = np.sqrt(ph / (ph - pe))
    mean_dep = np.sqrt((he + hh) / (he * hh))

    drho = 4.0
    ddep = 0.10

    xrho = np.linspace(mean_rho - drho, mean_rho + drho, N)
    xdep = np.linspace(mean_dep - ddep, mean_dep + ddep, N)

    Crho = 2 * L / np.sqrt(g * he * hh / (he + hh))
    Cdep = 2 * L / np.sqrt(g * (ph - pe) / ph)

    yrho = Crho * xrho
    ydep = Cdep * xdep

    return xrho, xdep, yrho / 3600, ydep / 3600  # in hours