# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer - Dashboard module

Interwave Analyzer - Version 2 (2026) 
Dashboard module version: 2.260305

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
import sys
import dash
import webbrowser
import numpy as np
import plotly.graph_objects as go

# internalwave analyzer packages 
import iwmod as mod
import iwplot as graph 

from dash import Dash, html, dcc, Input, Output
from datetime import datetime

# -------- Load data ------------------------------------------

if len(sys.argv) > 1:
    output_path = sys.argv[1]
else:
    raise ValueError("No output path provided to dashboard.")

data_path = os.path.join(output_path, "dash_data.npz")

if not os.path.exists(data_path):
    raise FileNotFoundError(f"dash_data.npz not found in: {output_path}")

data = np.load(data_path, allow_pickle=True)

nameBasin   = data["nameBasin"]

numSensor = int(data["numSensor"])
windFreq  = int(data["windFreq"])
tempFreq  = int(data["tempFreq"])

date = data["date"]
start = datetime.strptime(date[0], "%Y/%m/%d/%H/%M")
end   = datetime.strptime(date[-1], "%Y/%m/%d/%H/%M")

dx     = data["dx"]
dx_mod = data["dx_mod"]

iso    = data["iso"]
tau    = data["tau"]

temp   = data["temp"]              # (Nt, Nz)
depth  = np.asarray(data["depth"]).squeeze()

seu    = data["seu"]
mean_h = data["mean_h"]

ean          = data["ean"]
z0           = data["z0"]
hemod        = data["hemod"]
hemod_inter  = data["hemod_inter"]
mode2_nodes  = data["mode2_nodes"]
mode2_nodes_interpolated  = data["mode2_nodes_interpolated"]

freq  = data["freq"]
welch = data["welch"]
wr    = data["wr"]
conf  = data["conf"]
fo    = data["fo"]
m_n   = data["m_n"]
Bu    = data["Bu"]
cp    = data["cp"]

band     = data["band"]
low_per  = data["low_per"]
high_per = data["high_per"]


iw    = data["iw"]
dw    = data["dw"]
wedd  = data["wedd"]
riw   = data["riw"]
cond1 = data["cond1"]
cond2 = data["cond2"]
cond3 = data["cond3"]

schmidt    = data["schmidt"]
lakeNumber = data["lakeNumber"]

dura     = data["dura"]
dire     = data["dire"]
m_dw_spi = data["m_dw_spi"]

fdura    = data["fdura"]
fdire    = data["fdire"]

m_wast   = data["m_wast"]
wastmin  = data["wastmin"]
wastmax  = data["wastmax"]

pe = data["pe"]
ph = data["ph"]

m_he = data["m_he"]
m_hh = data["m_hh"]
m_h1 = data["m_h1"]
m_h2 = data["m_h2"]
m_h3 = data["m_h3"]

m_pe = data["m_pe"]
m_ph = data["m_ph"]
m_p1 = data["m_p1"]
m_p2 = data["m_p2"]
m_p3 = data["m_p3"]

c_he = data["c_he"]
c_hh = data["c_hh"]
c_h1 = data["c_h1"]
c_h2 = data["c_h2"]
c_h3 = data["c_h3"]

c_pe = data["c_pe"]
c_ph = data["c_ph"]
c_p1 = data["c_p1"]
c_p2 = data["c_p2"]
c_p3 = data["c_p3"]

t11 = data['T11']
t21 = data['T21']
t31 = data['T31']
t41 = data['T41']

durDire = data["durDire"]
durWave = data["durWave"]

ls_fetch  = data["ls_fetch"]
type_length = int(data["type_length"])

longData = data["longData"].item()

transData_raw = data["transData"]
transData = None if transData_raw is None else transData_raw.item()

m_ean = np.nanmean(ean)

n_slope      = data["n_slope"]
longRight    = data["longLeft"]
longLeft     = data["longRight"]
transRight   = data["transRight"]
transLeft    = data["transRight"]
atRight      = data["atLeft"]
atLeft       = data["atRight"]
atcrossRight = data["atcrossLeft"]
atcrossLeft  = data["atcrossRight"]

analysis_period = (
    f"{start.day}/{start.month}/{start.year} at {start.hour}h "
    f"to {end.day}/{end.month}/{end.year} at {end.hour}h"
)

g = 9.81

# -------------------------- Dash App ----------------------------------------


if getattr(sys, 'frozen', False):
    base_path = os.path.dirname(sys.executable)
else:
    base_path = os.path.dirname(os.path.abspath(__file__))

assets_path = os.path.join(base_path, "assets")

app = dash.Dash(__name__, suppress_callback_exceptions=True)

# --------- Figures - Basin characteristcs ------------------------------------

def basin_2d_figure(type_length, longData, transData, wind_dir):

    fig = go.Figure()

    depths = longData["depths"]
    dists  = longData["dists"]
    refs   = longData["refs"]

    # Orientation 
    theta_long = np.deg2rad(90 - longData["orientation"])
    cosL = np.cos(theta_long)
    sinL = np.sin(theta_long)

    theta_trans = theta_long - np.pi/2
    cosT = np.cos(theta_trans)
    sinT = np.sin(theta_trans)

    basin_x = []
    basin_y = []


    # Draw isos
    for i in range(len(depths)):

        radius_x = dists[i] / 2.0
        center_long = refs[i] + radius_x

        if type_length in [1, 2]:
            radius_y = radius_x
            center_trans = 0.0
        else:
            trans_depths = transData["depths"]
            trans_dists  = transData["dists"]
            trans_refs   = transData["refs"]

            sort_idx = np.argsort(trans_depths)
            trans_depths = trans_depths[sort_idx]
            trans_dists  = trans_dists[sort_idx]
            trans_refs   = trans_refs[sort_idx]

            interp_dist = np.interp(depths[i], trans_depths, trans_dists)
            interp_ref  = np.interp(depths[i], trans_depths, trans_refs)

            radius_y = interp_dist / 2.0
            center_trans = interp_ref + radius_y

        theta = np.linspace(0, 2*np.pi, 200)

        x_local = center_long + radius_x * np.cos(theta)
        y_local = center_trans + radius_y * np.sin(theta)

        x_global = x_local * cosL + y_local * cosT
        y_global = x_local * sinL + y_local * sinT

        x_km = x_global / 1000
        y_km = y_global / 1000

        fig.add_trace(go.Scatter(
            x=x_km,
            y=y_km,
            mode="lines",
            line=dict(color="black" if i == 0 else "gray"),
            name=f"{depths[i]} m"
        ))

        basin_x.extend(x_km)
        basin_y.extend(y_km)

    # True center (surface layer)
    radius_x0 = dists[0] / 2.0
    center_long0 = refs[0] + radius_x0

    if type_length in [1, 2]:
        radius_y0 = radius_x0
        center_trans0 = 0.0
    else:
        interp_dist0 = np.interp(depths[0], trans_depths, trans_dists)
        interp_ref0  = np.interp(depths[0], trans_depths, trans_refs)
        radius_y0 = interp_dist0 / 2.0
        center_trans0 = interp_ref0 + radius_y0

    x0 = (center_long0 * cosL + center_trans0 * cosT) / 1000
    y0 = (center_long0 * sinL + center_trans0 * sinT) / 1000


    # lock axis range to basin only 
    xmin, xmax = min(basin_x), max(basin_x)
    ymin, ymax = min(basin_y), max(basin_y)

    fig.update_layout(
        xaxis=dict(
            title="X (km)",
            scaleanchor="y",
            range=[xmin, xmax]
        ),
        yaxis=dict(
            title="Y (km)",
            range=[ymin, ymax]
        )
    )


    # Transect 
    L = max(xmax - xmin, ymax - ymin) * 100  # extremely long

    fig.add_trace(go.Scatter(
        x=[x0 - L*cosL, x0 + L*cosL],
        y=[y0 - L*sinL, y0 + L*sinL],
        mode="lines",
        line=dict(color="red", width=3),
        name="Longitudinal transect",
        hoverinfo="skip"
    ))

    if type_length == 3:
        fig.add_trace(go.Scatter(
            x=[x0 - L*cosT, x0 + L*cosT],
            y=[y0 - L*sinT, y0 + L*sinT],
            mode="lines",
            line=dict(color="blue", width=3),
            name="Transversal transect",
            hoverinfo="skip"
        ))

    # Wind arrow 
    theta_w = np.deg2rad(90 - wind_dir)
    dx_arrow = np.cos(theta_w)
    dy_arrow = np.sin(theta_w)

    arrow_length = (xmax - xmin) * 0.15

    fig.add_annotation(
        x=x0 + arrow_length*dx_arrow,
        y=y0 + arrow_length*dy_arrow,
        ax=x0,
        ay=y0,
        xref="x",
        yref="y",
        axref="x",
        ayref="y",
        showarrow=True,
        arrowhead=3,
        arrowwidth=3,
        arrowcolor="green"
    )

    fig.update_layout(
        margin=dict(l=60, r=20, t=40, b=60),
        legend=dict(orientation="h", y=-0.2)
    )

    return fig


def basin_3d_figure(type_length, longData, transData):

    depths = longData["depths"]
    dists  = longData["dists"]
    refs   = longData["refs"]

    # Orientation 
    theta_long = np.deg2rad(90 - longData["orientation"])
    cosL = np.cos(theta_long)
    sinL = np.sin(theta_long)

    theta_trans = theta_long - np.pi/2
    cosT = np.cos(theta_trans)
    sinT = np.sin(theta_trans)

    # Local grid 
    x_local = np.linspace(
        np.min(refs),
        np.max(refs + dists),
        250
    )

    y_extent = np.max(dists)
    y_local = np.linspace(-y_extent, y_extent, 250)

    X_local, Y_local = np.meshgrid(x_local, y_local)
    Z = np.zeros_like(X_local)

    # Sort from outer to inner 
    order = np.argsort(dists)[::-1]

    depths_sorted = depths[order]
    dists_sorted  = dists[order]
    refs_sorted   = refs[order]

    # Build depth field 
    for i in range(len(depths_sorted) - 1):

        rxo = dists_sorted[i] / 2
        cx_o = refs_sorted[i] + rxo

        rxi = dists_sorted[i+1] / 2
        cx_i = refs_sorted[i+1] + rxi

        if type_length in [1, 2]:
            ryo = rxo
            ryi = rxi
            cy_o = 0
            cy_i = 0
        else:
            trans_depths = transData["depths"]
            trans_dists  = transData["dists"]
            trans_refs   = transData["refs"]

            sort_idx = np.argsort(trans_depths)
            trans_depths = trans_depths[sort_idx]
            trans_dists  = trans_dists[sort_idx]
            trans_refs   = trans_refs[sort_idx]

            interp_do = np.interp(depths_sorted[i],   trans_depths, trans_dists)
            interp_di = np.interp(depths_sorted[i+1], trans_depths, trans_dists)

            interp_ro = np.interp(depths_sorted[i],   trans_depths, trans_refs)
            interp_ri = np.interp(depths_sorted[i+1], trans_depths, trans_refs)

            ryo = interp_do / 2
            ryi = interp_di / 2

            cy_o = interp_ro + ryo
            cy_i = interp_ri + ryi

        E_outer = ((X_local - cx_o)/rxo)**2 + ((Y_local - cy_o)/ryo)**2
        E_inner = ((X_local - cx_i)/rxi)**2 + ((Y_local - cy_i)/ryi)**2

        mask = (E_outer <= 1) & (E_inner >= 1)

        Z[mask] = (
            depths_sorted[i] +
            (depths_sorted[i+1] - depths_sorted[i]) *
            (1 - E_outer[mask])
        )

    # Deepest core 
    r_last = dists_sorted[-1] / 2
    cx_last = refs_sorted[-1] + r_last

    if type_length in [1, 2]:
        ry_last = r_last
        cy_last = 0
    else:
        interp_d = np.interp(depths_sorted[-1], trans_depths, trans_dists)
        interp_r = np.interp(depths_sorted[-1], trans_depths, trans_refs)
        ry_last = interp_d / 2
        cy_last = interp_r + ry_last

    E_last = ((X_local - cx_last)/r_last)**2 + ((Y_local - cy_last)/ry_last)**2
    Z[E_last <= 1] = depths_sorted[-1]

    # Rotate 
    X_global = X_local * cosL + Y_local * cosT
    Y_global = X_local * sinL + Y_local * sinT

    X_plot = X_global / 1000
    Y_plot = Y_global / 1000


    # Force true x and y equal scale 
    x_min, x_max = np.min(X_plot), np.max(X_plot)
    y_min, y_max = np.min(Y_plot), np.max(Y_plot)

    x_span = x_max - x_min
    y_span = y_max - y_min
    max_span = max(x_span, y_span)

    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2

    fig = go.Figure(data=[go.Surface(
        x=X_plot,
        y=Y_plot,
        z=Z,
        colorscale="Blues",
        showscale=False
    )])

    fig.update_layout(
        scene=dict(
            xaxis=dict(
                title="X (km)",
                range=[x_mid - max_span/2, x_mid + max_span/2]
            ),
            yaxis=dict(
                title="Y (km)",
                range=[y_mid - max_span/2, y_mid + max_span/2]
            ),
            zaxis=dict(
                title="Depth (m)",
                autorange="reversed"
            ),
            aspectmode="manual",
            aspectratio=dict(x=1, y=1, z=0.6)
        ),
        margin=dict(l=0, r=0, t=40, b=0)
    )

    return fig


def fetch_figure(dx, ls_fetch, time_index):

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=dx,
        y=ls_fetch,
        mode="lines",
        name="Fetch (m)"
    ))

    # vertical cursor line
    fig.add_vline(
        x=dx[time_index],
        line_width=2,
        line_dash="dash",
        line_color="red"
    )

    fig.update_layout(
        xaxis_title="Time",
        yaxis_title="Basin length (m)",
        margin=dict(l=40, r=20, t=40, b=40)
    )

    return fig

def basin_slopes_boxes(longLeft, longRight,transLeft=None, transRight=None,type_length=2):


    longLeft_val  = float(np.nanmean(longLeft))
    longRight_val = float(np.nanmean(longRight))

    transLeft_val = float(np.nanmean(transLeft)) if (
        type_length == 3 and transLeft is not None
    ) else 0

    transRight_val = float(np.nanmean(transRight)) if (
        type_length == 3 and transRight is not None
    ) else 0


    # Define maximum scale 
    vmax = max(
        longLeft_val,
        longRight_val,
        transLeft_val,
        transRight_val
    )

    vmax = vmax * 1.2 if vmax > 0 else 1

    # Define rows 
    rows = [
        ("Mean longitudinal slope at thermolcine depth (Left)", longLeft_val, vmax, "(°)", "{:.4f}"),
        ("Mean longitudinal slope at thermolcine depth (Right)", longRight_val, vmax, "(°)", "{:.4f}")
    ]

    if type_length == 3:
        rows.extend([
            ("Mean transversal slope at thermolcine depth (Left)", transLeft_val, vmax, "(°)", "{:.4f}"),
            ("Mean transversal slope at thermolcine depth (Right)", transRight_val, vmax, "(°)", "{:.4f}")
        ])


    return html.Div(
        children=[
            html.Div(
                style={"display": "flex", "marginBottom": "10px"},
                children=[
                    graduated_bar(label, val, vmax, unit, fmt),
                    html.Div()
                ]
            )
            for label, val, vmax, unit, fmt in rows
        ]
    )


# ---------------- Figures - Thermal stratification  --------------------------

def isotherms_figure(time, iso, tau):
    fig = go.Figure()
    for i, z in enumerate(tau):
        fig.add_trace(go.Scatter(x=time, y=iso[i], mode="lines",
                                 name=f"Isotherm {z} oC"))
    fig.update_layout(
        xaxis=dict(title="Date", type="date"),
        yaxis=dict(title="Height above bed (m)"),
        hovermode="x unified"
    )
    return fig


def thermal_variation_figure(time, temp, depth_point, depth):


    fig = go.Figure()
    Nt, Nz = temp.shape

    # Background temperature profiles 
    for j in range(Nz):
        fig.add_trace(go.Scatter(
            x=time,
            y=temp[:, j],
            mode="lines",
            line=dict(color="black", width=0.6),
            opacity=0.15,
            showlegend=False,
            hoverinfo="skip"
        ))


    # Highlighted depths
    colors = ["#d62728", "#8c1b13", "#1f77b4", "#0b1d51"]

    for i in range(min(4, len(depth_point))):
        idx = depth_point[i]

        if idx is None:
            continue

        idx = int(idx)

        # Skip inactive or invalid indices
        if idx < 0 or idx >= Nz:
            continue

        # Skip empty depth arrays
        if depth[i] is None or len(depth[i]) == 0:
            label = "inactive"
        else:
            label = f"{np.nanmean(depth[i]):.1f} m"

        fig.add_trace(go.Scatter(
            x=time,
            y=temp[:, idx],
            mode="lines",
            line=dict(color=colors[i], width=2),
            name=label
        ))


    # Layout
    fig.update_layout(
        xaxis=dict(title="Date", type="date"),
        yaxis=dict(title="Water temperature (°C)"),
        hovermode="x unified",
        margin=dict(l=60, r=20, t=30, b=50)
    )

    return fig


def temperature_profile_with_N2(t_index):
    T = temp[t_index, :]
    z = mean_h

    rho = mod.commission(T)
    drhodz = np.gradient(rho, z)
    rho0 = np.nanmean(rho)
    arg = -(g / rho0) * drhodz
    N = np.sign(arg) * np.sqrt(np.abs(arg))

    fig = go.Figure()

    fig.add_trace(go.Scatter(x=T, y=z, mode="lines+markers",
                             name="Temperature"))
    fig.add_trace(go.Scatter(x=N, y=z, mode="lines",
                             name="N", xaxis="x2",
                             line=dict(dash="dot")))

    fig.update_layout(
        yaxis=dict(title="Height above bed (m)"),
        xaxis=dict(title="Temperature (°C)"),
        xaxis2=dict(
            title="Buoyancy frequency (Hz)",
            overlaying="x",
            side="top"
        ),
        margin=dict(l=70, r=40, t=40, b=60),
        legend=dict(orientation="h", y=-0.2)
    )
    return fig


def thermal_structure_lines_figure():
    h1 = mode2_nodes[:, 0] + z0
    h2 = mode2_nodes[:, 1] + z0
    he = np.array([np.asarray(x).flatten()[0] for x in hemod]) + z0

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=dx, y=ean, name="Surface"))
    fig.add_trace(go.Scatter(x=dx_mod, y=h1, name="Upper metalimnion"))
    fig.add_trace(go.Scatter(x=dx_mod, y=h2, name="Lower metalimnion"))
    fig.add_trace(go.Scatter(x=dx_mod, y=he, name="Thermocline"))

    fig.update_layout(
        xaxis=dict(title="Date",type="date"),
        yaxis=dict(title="Height above bed (m)"),
        hovermode="x unified"
    )
    return fig


def schmidt_figure(dx, schmidt):

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=dx,
        y=schmidt,
        mode="lines",
        name="Schmidt Stability"
    ))

    fig.update_layout(
        xaxis_title="Time",
        yaxis_title="Schmidt Stability (J/m²)",
        margin=dict(l=40, r=20, t=40, b=40)
    )

    return fig


def buoyancy_frequency_figure(dx, n_slope):

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=dx,
        y=n_slope,
        mode="lines",
        line=dict(width=2),
        name="Mean buoyancy frequency"
    ))

    fig.update_layout(
        xaxis_title="Time",
        yaxis_title="Mean buoyancy frequency (Hz)",
        margin=dict(l=60, r=20, t=40, b=60),
        hovermode="closest"
    )

    return fig

def thermal_metrics_boxes(
    m_he, m_hh, m_h1, m_h2, m_h3,
    m_pe, m_ph, m_p1, m_p2, m_p3, 
    c_he, c_hh, c_h1, c_h2, c_h3,
    c_pe, c_ph, c_p1, c_p2, c_p3, m_ean, z0, T):
    
    rho = mod.commission(T)
    rho_min = np.nanmin(rho)
    rho_max = np.nanmax(rho)
    
    rows = [
    
        # Thickness min = 0
        ("Epilimnion thickness",  m_he, c_he, 0, m_ean-z0, "m", "{:.2f}"),
        ("Hypolimnion thickness", m_hh, c_hh, 0, m_ean-z0, "m", "{:.2f}"),
        ("For three-layer system, epilimnion thickness",  m_h1, c_h1, 0, m_ean-z0, "m", "{:.2f}"),
        ("For three-layer system, metalimnion thickness", m_h2, c_h2, 0, m_ean-z0, "m", "{:.2f}"),
        ("For three-layer system, hypolimnion thickness", m_h3, c_h3, 0, m_ean-z0, "m", "{:.2f}"),
    
        # Density min = rho_min
        ("Mean epilimnion density",  m_pe, c_pe, rho_min, rho_max, "kg/m3", "{:.2f}"),
        ("Mean hypolimnion density", m_ph, c_ph, rho_min, rho_max, "kg/m3", "{:.2f}"),
        ("For three-layer system, mean epilimnion density",  m_p1, c_p1, rho_min, rho_max, "kg/m3", "{:.2f}"),
        ("For three-layer system, mean metalimnion density", m_p2, c_p2, rho_min, rho_max, "kg/m3", "{:.2f}"),
        ("For three-layer system, mean hypolimnion density", m_p3, c_p3, rho_min, rho_max, "kg/m3", "{:.2f}"),
    ]

    return html.Div(
        children=[
            html.Div(
                style={
                    "display": "flex",
                    "marginBottom": "10px"
                },
                children=[
                    graduated_bar_limited(label, val, conf, vmin, vmax, unit, fmt),
                    html.Div()
                ]
            )
            for label, val, conf, vmin, vmax, unit, fmt in rows
        ]
    )


# ---------------------- Figures - Thermal mixing -----------------------------

def wind_wedd_figure():
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=dx, y=iw, name="Wind (iw)"))
    fig.add_trace(go.Scatter(x=dx, y=wedd, name="Wedderburn",
                             yaxis="y2", line=dict(dash="dot")))
    fig.update_layout(
        xaxis=dict(type="date"),
        yaxis=dict(title="Wind (m/s)"),
        yaxis2=dict(title="Wedderburn", overlaying="y",
                    side="right", type="log"),
        hovermode="x"
    )
    return fig


def isotherms_zoom_figure(t, window=120):
    i0 = max(0, t - window)
    i1 = min(len(dx), t + window)

    fig = go.Figure()
       
    for i, z in enumerate(tau):
            fig.add_trace(go.Scatter(x=dx[i0:i1], y=iso[i, i0:i1], 
                        mode="lines", name=f"Isotherm {z} oC"))

    fig.update_layout(
        xaxis=dict(type="date"),
        yaxis=dict(title="Height above bed (m)"),
        hovermode="x unified"
    )
    return fig


def ri_signal_figure(ri, c1, c2, c3):

    fig = go.Figure()

    regimes = [
        (0,   c1, "Mixing is so rapid that billows have no time to form. The lake is for all practical purposes homogeneous.", "#2c7bb6"),
        (c1,  c2, "Interface shear, Kelvin-Helmholtz billows, and rapid deepening accompany large interface displacements. Internal seiche are not observed because wave period is larger than entrainment time", "#abd9e9"),
        (c2,  c3, "Internal seiching is the prominent feature of this regime", "#fdae61"),
        (c3,  max(c3*5, ri*1.3), "Boyancy dominates all processes, internal wave have short periods and small amplitudes", "#d7191c")
    ]

    for r0, r1, label, color in regimes:
        fig.add_trace(go.Bar(
            x=[r1 - r0],
            y=["Ri"],
            base=r0,
            orientation="h",
            marker=dict(color=color),
            hovertemplate=f"{label}<br>{r0:.2f} – {r1:.2f}<extra></extra>",
            showlegend=False
        ))

    # Ri marker
    fig.add_trace(go.Scatter(
        x=[ri],
        y=["Ri"],
        mode="markers+text",
        marker=dict(size=18, color="black"),
        text=[f"Ri = {ri:.2f}"],
        textposition="top center",
        showlegend=False
    ))

    fig.update_layout(
        xaxis=dict(
            title="Richardson number",
            type="log",
            range=[
                np.log10(max(1e-4, min(ri, c1)*0.5)),
                np.log10(max(c3*5, ri*1.5))
            ]
        ),
        yaxis=dict(showticklabels=False),
        height=220,
        margin=dict(l=60, r=40, t=60, b=40)
    )

    return fig

def lake_number_figure(dx, lakeNumber):

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=dx,
        y=lakeNumber,
        mode="lines",
        name="Lake Number",
    ))

    fig.update_layout(
        xaxis_title="Time",
        yaxis=dict(
            title="Lake Number",
            type="log"
        ),
        margin=dict(l=40, r=20, t=40, b=40)
    )

    return fig

def mixing_metrics_boxes(m_wast, wastmin, wastmax,dura, dire, fdura, fdire, t11):

    rows = [

        ("Mean friction velocity", m_wast, 0.01, "m/s", "{:.8f}"),
        ("Min friction velocity", wastmin, 0.01, "m/s", "{:.8f}"),
        ("Max friction velocity of the wind", wastmax, 0.01, "m/s", "{:.8f}"),
        ("Strongest wind duration", dura, t11[1]/60/60, "h", "{:.1f}"),
        ("Strongest wind duration considering homogeneous winds", dire, t11[1]/60/60, "h", "{:.1f}"),
        ("Wedderburn reduction (duration)", fdura, 1, "", "{:.2f}"),
        ("Wedderburn reduction (direction)", fdire, 1, "", "{:.2f}"),
    ]

    return html.Div(
        children=[
            html.Div(
                style={
                    "display": "flex",
                    "marginBottom": "10px"
                },
                children=[
                    graduated_bar(label, val, vmax, unit, fmt)
,
                    html.Div()  
                ]
            )
            for label, val, vmax, unit, fmt in rows

        ]
    )


# -------------- Figures - Internal waves analysis ----------------------------

def spectral_isotherm_figure(tau, freq, welch, wr, conf, m_n, fo):

    colors = ["black", "dimgray", "salmon", "red"]
    fig = go.Figure()

    for i in range(len(tau)):
        if tau[i] == -999:
            continue

        f = np.asarray(freq[i])
        p = np.asarray(welch[i])

        fig.add_trace(go.Scatter(
            x=f,
            y=p,
            mode="lines",
            line=dict(color=colors[i], width=1),
            name=f"{tau[i]} °C"
        ))

        # significance curve
        fig.add_trace(go.Scatter(
            x=wr[i],
            y=conf[i],
            mode="lines",
            line=dict(color=colors[i], dash="dash", width=1),
            showlegend=False
        ))

    # vertical reference lines
    fo = float(fo)
    m_n = float(m_n)
    
    shapes = [
        dict(
            type="line",
            x0=fo, x1=fo,
            y0=0, y1=1,
            xref="x",
            yref="paper",
            line=dict(color="black", dash="dash")
        ),
        dict(
            type="line",
            x0=m_n, x1=m_n,
            y0=0, y1=1,
            xref="x",
            yref="paper",
            line=dict(color="black", dash="dot")
        )
    ]

    fig.update_layout(
        xaxis=dict(
            title="Frequency (Hz)",
            type="log"
        ),
        yaxis=dict(
            title="Spectral energy (m²/Hz)",
            type="log"
        ),
        shapes=shapes,
        legend=dict(
            orientation="h",
            y=-0.2
        ),
        margin=dict(l=70, r=20, t=40, b=70),
        hovermode="closest"
    )

    return fig





def spectral_metrics_boxes(t11,t21,t31,t41,durDire, durWave):
    rows1 = [
        ("Theoretical fundamental internal seiche period", t11[1]/60/60, t11[3]/60/60, 0, t41[2]/60/60, "hours", "{:.2f}"),
        ("Theoretical internal seiche period of second vertical mode", t21[1]/60/60, 0, t21[3]/60/60, t41[2]/60/60, "hours", "{:.2f}"),
        ("Theoretical internal seiche period of thrid vertical mode",  t31[1]/60/60, 0, t31[3]/60/60, t41[2]/60/60, "hours", "{:.2f}"),
        ("Theoretical internal seiche period of fourth vertical mode", t41[1]/60/60, 0, t41[3]/60/60, t41[2]/60/60, "hours", "{:.2f}"),
    ]

    rows2 = [
        ("Theoretical phase speed of fundamental wave", cp, 1, "m/s", "{:.2f}"),
        ("Time of wind events favoring internal seiches", durWave, 100, "%", "{:.2f}"),
        ("Time of wind events favoring internal seiches considering only homogeneous winds direction", durDire, 100, "%", "{:.2f}"),
        ("Burger number", Bu, 2, "(-)", "{:.2f}"),
    ]
    return html.Div(
        children=[

            # ---- with confidence ----
            *[
                html.Div(
                    style={"marginBottom": "10px"},
                    children=graduated_bar_limited(
                        label, val, conf, vmin, vmax, unit, fmt
                    )
                )
                for label, val, conf, vmin, vmax, unit, fmt in rows1
            ],

            html.Hr(style={"margin": "15px 0"}),

            # ---- without confidence ----
            *[
                html.Div(
                    style={"marginBottom": "10px"},
                    children=graduated_bar(
                        label, val, vmax, unit, fmt
                    )
                )
                for label, val, vmax, unit, fmt in rows2
            ]
        ]
    )


def isotherm_bandpass_figure(dx, band, tau, low_per, high_per):


    fig = go.Figure()

    colors = ["black", "dimgray", "salmon", "red"]

    for i in range(4):
        if tau[i] != -999:
            fig.add_trace(go.Scatter(
                x=dx,
                y=band[i],
                mode="lines",
                line=dict(color=colors[i], width=1),
                name=f"{tau[i]}°C"
            ))

    fig.update_layout(
        xaxis=dict(
            title="Time",
            type="date",
            showgrid=True,
            gridcolor="rgba(0,0,0,0.15)"
        ),
        yaxis=dict(
            title="Displacement (m)",
            showgrid=True,
            gridcolor="rgba(0,0,0,0.15)"
        ),
        title=(
            f"Periods selected by the band-pass filter: {low_per[0]}–{high_per[0]} h"
        ),
        legend=dict(
            orientation="h",
            y=-0.25
        ),
        margin=dict(l=70, r=40, t=60, b=70),
        height=320
    )

    return fig

def degeneration_evolution_figure(he, H, W, pe, ph, meta, ls):


    fig = go.Figure()
    
    heH = he/H
    if heH > 0.5:
        heH = 1 - heH

    # x-axis domain
    xl = np.linspace(0.01, 0.5, 400)

    # Theoretical regimes
    fig.add_trace(go.Scatter(
        x=xl,
        y=graph.equation_sd(H, xl, pe, ph, ls, meta),
        mode="lines",
        name="SD regime",
        line=dict(color="blue")
    ))

    fig.add_trace(go.Scatter(
        x=xl,
        y=graph.equation_bo(H, xl),
        mode="lines",
        name="BO regime",
        line=dict(color="black")
    ))

    fig.add_trace(go.Scatter(
        x=xl,
        y=graph.equation_kh(H, xl, meta),
        mode="lines",
        name="KH regime",
        line=dict(color="red")
    ))

    # Current state 
    fig.add_trace(go.Scatter(
        x=[heH],
        y=[1.0 / W],
        mode="markers",
        marker=dict(size=14, color="red", line=dict(width=2, color="red")),
        name="Current state"
    ))

    fig.update_layout(
        xaxis_title=r"he/H",
        yaxis_title="1/W",
        hovermode="closest",
        yaxis=dict(range=[0, 2]),
        xaxis=dict(range=[0.1, 0.5]),
        margin=dict(l=60, r=20, t=30, b=50)
    )

    return fig

def slope_criticality_figure(dx, time_index, atLeft=None, atRight=None, atcrossLeft=None, atcrossRight=None, basin_type=2):


    # Theoretical curve 
    at = np.logspace(-3, 3, 600)

    def compute_y(val):
        return 0.9894553 + (0.1745722 - 0.9894553) / (
            1 + (val / 0.919315) ** 3.812312
        )

    y = compute_y(at)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=at,
        y=y,
        mode="lines",
        line=dict(color="black", width=2),
        name="Theoretical relation"
    ))

    # Longitudinal (blue) 
    if atLeft is not None:
        val = atLeft[time_index]
        if val > 0:
            fig.add_trace(go.Scatter(
                x=[val],
                y=[compute_y(val)],
                mode="markers",
                marker=dict(
                    size=11,
                    color="blue",
                    symbol="circle",
                    line=dict(color="black", width=1)
                ),
                name="Longitudinal – Left"
            ))

    if atRight is not None:
        val = atRight[time_index]
        if val > 0:
            fig.add_trace(go.Scatter(
                x=[val],
                y=[compute_y(val)],
                mode="markers",
                marker=dict(
                    size=11,
                    color="blue",
                    symbol="square",
                    line=dict(color="black", width=1)
                ),
                name="Longitudinal – Right"
            ))

    # Transversal (red) 
    if basin_type == 3:

        if atcrossLeft is not None:
            val = atcrossLeft[time_index]
            if val > 0:
                fig.add_trace(go.Scatter(
                    x=[val],
                    y=[compute_y(val)],
                    mode="markers",
                    marker=dict(
                        size=11,
                        color="red",
                        symbol="circle",
                        line=dict(color="black", width=1)
                    ),
                    name="Transversal – Left"
                ))

        if atcrossRight is not None:
            val = atcrossRight[time_index]
            if val > 0:
                fig.add_trace(go.Scatter(
                    x=[val],
                    y=[compute_y(val)],
                    mode="markers",
                    marker=dict(
                        size=11,
                        color="red",
                        symbol="square",
                        line=dict(color="black", width=1)
                    ),
                    name="Transversal – Right"
                ))

    # Layout 
    fig.update_layout(
        xaxis=dict(
            title="Slope criticality parameter (-)",
            type="log"
        ),
        yaxis=dict(
            title="Normalized wave amplitude (-)"
        ),
        legend=dict(
            orientation="h",
            y=-0.2
        ),
        margin=dict(l=60, r=20, t=40, b=60),
        hovermode="closest"
    )

    return fig

# ------------------ Graduated bar definitions --------------------------------

def graduated_bar(label, value, vmax, unit="", fmt="{:.3g}"):


    ratio = max(0, min(value / vmax, 1))

    # Color by magnitude
    if ratio < 0.33:
        color = "#4fa3ff"
    elif ratio < 0.66:
        color = "#4caf50"
    elif ratio < 0.85:
        color = "#ff9800"
    else:
        color = "#e53935"
    
    value_str = fmt.format(value)

    return html.Div(
        style={
            "border": "1px solid #d0d0d0",
            "borderRadius": "8px",
            "padding": "8px 14px",
            "backgroundColor": "#f8f9fa",
            "width": "100%"
        },
        children=[

            # Label
            html.Div(
                f"{label}: {value_str} {unit}",
                style={
                    "fontWeight": "600",
                    "marginBottom": "6px",
                    "fontSize": "13.5px",
                    "color": "#333"
                }
            ),

            # Bar
            html.Div(
                style={
                    "position": "relative",
                    "height": "14px",
                    "backgroundColor": "#e0e0e0",
                    "borderRadius": "7px",
                    "border": "1px solid #ccc"
                },
                children=[

                    # Filled part
                    html.Div(
                        style={
                            "width": f"{ratio*100:.1f}%",
                            "height": "100%",
                            "backgroundColor": color,
                            "borderRadius": "7px"
                        }
                    ),

                    # Current value marker
                    html.Div(
                        style={
                            "position": "absolute",
                            "left": f"{ratio*100:.1f}%",
                            "top": "-2px",
                            "width": "2px",
                            "height": "18px",
                            "backgroundColor": "#222"
                        }
                    )
                ]
            ),

            # Scale labels
            html.Div(
                style={
                    "display": "flex",
                    "justifyContent": "space-between",
                    "fontSize": "11px",
                    "color": "#555",
                    "marginTop": "2px"
                },
                children=[
                    html.Span("0"),
                    html.Span(f"{vmax:g} {unit}")
                ]
            )
        ]
    )

def graduated_bar_limited(label, value, conf, vmin, vmax, unit="", fmt="{:.3g}"):

    # Prevent division by zero
    if vmax == vmin:
        ratio = 0
    else:
        ratio = (value - vmin) / (vmax - vmin)

    ratio = max(0, min(ratio, 1))

    # Color by magnitude
    if ratio < 0.33:
        color = "#4fa3ff"
    elif ratio < 0.66:
        color = "#4caf50"
    elif ratio < 0.85:
        color = "#ff9800"
    else:
        color = "#e53935"
    
    value_str = fmt.format(value)
    limit_str = fmt.format(conf)

    return html.Div(
        style={
            "border": "1px solid #d0d0d0",
            "borderRadius": "8px",
            "padding": "8px 14px",
            "backgroundColor": "#f8f9fa",
            "width": "100%"
        },
        children=[

            html.Div(
                f"{label}: {value_str} {unit} ± {limit_str} {unit}",
                style={
                    "fontWeight": "600",
                    "marginBottom": "6px",
                    "fontSize": "13.5px",
                    "color": "#333"
                }
            ),

            html.Div(
                style={
                    "position": "relative",
                    "height": "14px",
                    "backgroundColor": "#e0e0e0",
                    "borderRadius": "7px",
                    "border": "1px solid #ccc"
                },
                children=[

                    html.Div(
                        style={
                            "width": f"{ratio*100:.1f}%",
                            "height": "100%",
                            "backgroundColor": color,
                            "borderRadius": "7px"
                        }
                    ),

                    html.Div(
                        style={
                            "position": "absolute",
                            "left": f"{ratio*100:.1f}%",
                            "top": "-2px",
                            "width": "2px",
                            "height": "18px",
                            "backgroundColor": "#222"
                        }
                    )
                ]
            ),

            html.Div(
                style={
                    "display": "flex",
                    "justifyContent": "space-between",
                    "fontSize": "11px",
                    "color": "#555",
                    "marginTop": "2px"
                },
                children=[
                    html.Span(f"{vmin:g}"),
                    html.Span(f"{vmax:g} {unit}")
                ]
            )
        ]
    )

# ---------------------------- Layout------------------------------------------

app.layout = html.Div(
    style={"width": "90%", "margin": "auto"},
    children=[

        # Header
        html.Div(
            style={
                "display": "flex",
                "alignItems": "center",
                "gap": "20px",
                "marginBottom": "10px"
            },
            children=[
                html.Img(src="/assets/iwlogo.png", style={"height": "60px"}),

                html.Div([
                    html.H1("Interwave Analyzer Dashboard", style={"margin": "0"}),
                    html.P(
                        "Analyzing hydrodynamics and mixing in lakes and reservoirs",
                        style={"margin": "0", "fontSize": "16px", "color": "#555"}
                    )
                ])
            ]
        ),

        html.Hr(),

        # Info box
        html.Div(
            style={"display": "flex", "gap": "20px", "marginBottom": "20px"},
            children=[

                html.Div(
                    style={
                        "backgroundColor": "#f8f9fa",
                        "border": "1px solid #ddd",
                        "borderRadius": "6px",
                        "padding": "12px 16px",
                        "flex": "1"
                    },
                    children=[
                        html.H4("General information", style={"marginTop": "0"}),
                        html.P([html.B("Basin name: "), nameBasin], style={"margin": "0"}),
                        html.P([html.B("Period analysis: "), analysis_period], style={"margin": "0"})
                    ]
                ),

                html.Div(
                    style={
                        "backgroundColor": "#f8f9fa",
                        "border": "1px solid #ddd",
                        "borderRadius": "6px",
                        "padding": "12px 16px",
                        "flex": "1"
                    },
                    children=[
                        html.H4("Sensor measurements", style={"marginTop": "0"}),
                        html.P([html.B("Number of temperature sensors: "), f"{numSensor}"], style={"margin": "0"}),
                        html.P([html.B("Time step of temperature data: "), f"{windFreq} min"], style={"margin": "0"}),
                    ]
                )
            ]
        ),

        dcc.Dropdown(
            id="page-selector",
            options=[
                {"label": "Basin characteristics", "value": "basin"},
                {"label": "Thermal stratification", "value": "thermal"},
                {"label": "Thermal mixing", "value": "mixing"},
                {"label": "Internal waves", "value": "spectral"}  
            ],
            value="basin",
            clearable=False,
            style={"width": "300px"}
        ),

        html.Div(id="page-content")
    ]
)

# ------------------------ Page Content ---------------------------------------

@app.callback(
    Output("page-content", "children"),
    Input("page-selector", "value")
)
def render_page(page):

    if page == "basin":
        return [
            html.H4("Basin characteristics"),
            
            html.Div(
                style={"display": "flex", "gap": "20px"},
                children=[
            
                    html.Div(
                        style={"flex": "1"},
                        children=[
                            dcc.Graph(id="basin-2d")
                        ]
                    ),
            
                    html.Div(
                        style={"flex": "1"},
                        children=[
                            dcc.Graph(
                                figure=basin_3d_figure(
                                    type_length,
                                    longData,
                                    transData
                                )
                            )
                        ]
                    )
                ]
            ),
            
            # Fetch below (full width)
            dcc.Graph(id="fetch-timeseries"),
            
            html.H4("Basin slopes"),

            basin_slopes_boxes(longLeft, longRight, transLeft, transRight,
                type_length)

        ]    

    if page == "thermal":
        return [

            html.H4("Isotherms"),
            dcc.Graph(figure=isotherms_figure(dx, iso, tau)),

            html.H4(id="temp-title",
                    children=f"Temperature variation: {dx[0]}"),

            html.Div(style={"display": "flex"}, children=[
                dcc.Graph(id="temp-timeseries",
                          figure=thermal_variation_figure(dx, temp, seu, depth),
                          style={"flex": "3"}),
                dcc.Graph(id="temp-profile", style={"flex": "1"})
            ]),

            html.H4("Thermal structure"),
            dcc.Graph(figure=thermal_structure_lines_figure()),
            
            dcc.Graph(figure=schmidt_figure(dx, schmidt)),
            
            dcc.Graph(figure=buoyancy_frequency_figure(dx, n_slope)),

            # Informative value
            html.H4("Thermal summary indicators"),

                thermal_metrics_boxes(
                    m_he, m_hh, m_h1, m_h2, m_h3,
                    m_pe, m_ph, m_p1, m_p2, m_p3, 
                    c_he, c_hh, c_h1, c_h2, c_h3,
                    c_pe, c_ph, c_p1, c_p2, c_p3, m_ean, z0, temp)
        ]

    if page == "mixing":
        return [
    
            html.H4("Wind forcing and mixing indicators"),
    
            # Wind and isotherms zoom 
            html.Div(
                style={"display": "flex", "gap": "20px"},
                children=[
    
                    # Wind and Wedderburn number
                    html.Div(
                        style={"flex": "3"},
                        children=[
                            dcc.Graph(
                                id="wind-graph",
                                figure=wind_wedd_figure()
                            )
                        ]
                    ),
    
                    # Isotherms zoom
                    html.Div(
                        style={"flex": "2"},
                        children=[
                            dcc.Graph(id="iso-zoom")
                        ]
                    )
                ]
            ),
    
            html.Br(),
    
            # Richardson (lake mixing regimes)
            html.Div(
                style={"width": "100%"},
                children=[
                    dcc.Graph(
                        id="ri-signal",
                        style={"height": "220px"}
                    )
                ]
            ),
            
            dcc.Graph(figure=lake_number_figure(dx, lakeNumber)),
            
            # Informative values
            html.H4("Mixing summary indicators"),

                mixing_metrics_boxes(
                    m_wast, wastmin, wastmax,
                    dura, dire, fdura, fdire, t11)
        ]
    
    if page == "spectral":
        return [

        html.Div(
            style={
                "display": "flex",
                "gap": "20px",
                "alignItems": "flex-start"
            },
            children=[

                # PSD
                html.Div(
                    style={"flex": "1"},
                    children=[
                        html.H4("Spectral energy of isotherms:"),
                        dcc.Graph(
                            figure=spectral_isotherm_figure(
                                tau=tau,
                                freq=freq,
                                welch=welch,
                                wr=wr,
                                conf=conf,
                                m_n=m_n,
                                fo=fo
                            ),
                            style={"height": "650px"}
                        )
                    ]
                ),

                # Summary bars
                html.Div(
                    style={"flex": "1"},
                    children=[
                        html.H4("Theoretical internal seiche periods"),
                        spectral_metrics_boxes(t11,t21,t31,t41, 
                                               durDire, durWave)
                    ]
                )
            ]
        ),
        html.Br(),
    
        # Band-pass filtered isotherms 
        dcc.Graph( figure=isotherm_bandpass_figure(dx=dx, band=band, tau=tau,
                low_per=low_per, high_per=high_per)
        ),
        
        html.Div([
            dcc.Graph(
                id="iso-timeseries",
                figure=isotherms_figure(dx, iso, tau)
                ),
    
            html.Div([
                dcc.Graph(
                    id="degeneration-graph",
                    style={"width": "50%", "display": "inline-block"}
                ),
                dcc.Graph(
                    id="secondary-graph",
                    style={"width": "50%", "display": "inline-block"}
                ),
            ])
        ])
    ]
    

# -------------------------------- Callbacks ----------------------------------

@app.callback(
    Output("basin-2d", "figure"),
    Output("fetch-timeseries", "figure"),
    Input("fetch-timeseries", "hoverData")
)
def update_basin(hoverData):

    # Default time index
    time_index = 0

    if hoverData is not None:
        time_index = hoverData["points"][0]["pointIndex"]

    wind_dir = dw[time_index]

    fig2d = basin_2d_figure(
        type_length,
        longData,
        transData,
        wind_dir
    )

    fig_fetch = fetch_figure(dx, ls_fetch, time_index)

    return fig2d, fig_fetch

@app.callback(
    Output("temp-title", "children"),
    Output("temp-profile", "figure"),
    Input("temp-timeseries", "hoverData")
)
def update_profile(hoverData):
    t = hoverData["points"][0]["pointIndex"] if hoverData else 0
    return (
        f"Temperature variation: {dx[t]}",
        temperature_profile_with_N2(t)
    )


@app.callback(
    Output("iso-zoom", "figure"),
    Output("ri-signal", "figure"),
    Input("wind-graph", "hoverData")
)
def update_mixing(hoverData):

    t = hoverData["points"][0]["pointIndex"] if hoverData else 0

    return (
        isotherms_zoom_figure(t),
        ri_signal_figure(
            riw[t],
            cond1[t],
            cond2[t],
            cond3[t]
        )
    )

@app.callback(
    Output("degeneration-graph", "figure"),
    Output("secondary-graph", "figure"),
    Input("iso-timeseries", "hoverData")
)
def update_degeneration(hoverData):

    if hoverData is None:
        raise dash.exceptions.PreventUpdate


    # Time index 
    t_idx = hoverData["points"][0]["pointIndex"]

    # Time-dependent values 
    ean_t = ean[t_idx]
    he_t  = np.asarray(hemod_inter[t_idx]).flatten()[0]
    H     = ean_t - z0
    W     = wedd[t_idx]
    pe_t  = pe[t_idx]
    ph_t  = ph[t_idx]
    meta  = mode2_nodes_interpolated[t_idx, 0] - mode2_nodes_interpolated[t_idx, 1]
    ls    = ls_fetch[t_idx]

    # Degeberation 
    fig_deg = degeneration_evolution_figure(
        he=he_t,
        H=H,
        W=W,
        pe=pe_t,
        ph=ph_t,
        meta=meta,
        ls=ls
    )

    # Slope critically parameter (At)
    fig_secondary = slope_criticality_figure(
        dx=dx,
        time_index=t_idx,
        atLeft=atLeft,
        atRight=atRight,
        atcrossLeft=atcrossLeft,
        atcrossRight=atcrossRight,
        basin_type=type_length
    )

    return fig_deg, fig_secondary



# ------------------------------- Run -----------------------------------------
def start_dashboard(output_path, port):

    global app

    webbrowser.open(f"http://127.0.0.1:{port}")
    app.run(debug=False, use_reloader=False, port=port)


if __name__ == "__main__":

    if len(sys.argv) > 2:
        port = int(sys.argv[2])
    else:
        port = 8050

    start_dashboard(output_path, port)
