# -*- coding: utf-8 -*-
"""
Created by Rafael de Carvalho Bueno 
Interwave Analyzer - Graphical User Interface module

Interwave Analyzer - Version 2 (2026) 
Graphical User Interface module version: 2.260305

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
import random
import platform
import subprocess

try:
    qt_backend = "PyQt5"
    qt_already_loaded = any("PyQt5" in m for m in sys.modules)

    from PyQt5 import QtCore, QtGui, QtWidgets
    
except:
    try:
        qt_backend = "PySide6"
        from PySide6 import QtCore, QtGui, QtWidgets
    except:
        print("[ERROR] No compatible Qt backend found.")

if getattr(sys, 'frozen', False) and len(sys.argv) > 2:
    import iwdash
    iwdash.start_dashboard(sys.argv[1], int(sys.argv[2]))
    sys.exit()

# Register aliases
sys.modules["QtCore"] = QtCore
sys.modules["QtGui"] = QtGui
sys.modules["QtWidgets"] = QtWidgets

# Adjust QAction location for compatibility
if qt_backend == "PyQt5":
    QtGui.QAction = QtWidgets.QAction

print(f"[INFO] Using Qt backend: {qt_backend}")


from QtCore import Qt
from QtGui import QIcon, QAction
from QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget, QVBoxLayout, QGridLayout,
    QLabel, QPushButton, QCheckBox, QLineEdit, QSpinBox, QFrame, QFileDialog,
    QMessageBox, QSizePolicy, QRadioButton, QButtonGroup, QComboBox, QTextEdit
)


import webbrowser

# ---------------------- Callback functions -----------------------------------
def open_button():
    global path_temp, path_meteo, path_senso, folder_path
    global path_fetch, path_level, path_bathy

    path_open, _ = QFileDialog.getOpenFileName(None, "Open Settings File", "", "SET Files (*.set)")

    if not path_open:
        return

    with open(path_open, 'r') as reader:
        lines = [line.strip() for line in reader]

    # Create a line iterator
    it = iter(lines)

    # Paths
    path_temp = next(it)
    path_meteo = next(it)
    path_senso = next(it)

    # Height of wind
    height_wind.setText(next(it))

    # Radiation checkbox
    radiation_on.setChecked(int(next(it)))

    # Reference level
    reference_level.setText(next(it))

    # Contribution wind variable
    contri_wind.setValue(int(float(next(it))))

    # Latitude
    latitude.setValue(int(float(next(it))))

    # DPI
    dpi_spin.setValue(int(next(it)))

    # Basin length selection
    typechoose_id = int(next(it))
    if typechoose_id == 1:
        rad1.setChecked(True)
        vari_len.setEnabled(True)
        file_len.setEnabled(False)
    elif typechoose_id == 2:
        rad2.setChecked(True)
        vari_len.setEnabled(False)
        file_len.setEnabled(True)

    if typechoose_id == 1:
        vari_len.setText(next(it))
    elif typechoose_id == 2:
        path_fetch = next(it)

    # Levels
    typelevel_id = int(next(it))
    if typelevel_id == 1:
        rad_level_uniform.setChecked(True)
        unif_level.setEnabled(True)
        path_level = ''
    elif typelevel_id == 2:
        rad_level_file.setChecked(True)
        unif_level.setEnabled(False)
        path_level = next(it)

    if typelevel_id == 1:
        unif_level.setText(next(it))

    # Meta
    meta.setText(next(it))

    # Filter
    filterchoose_id = int(next(it))
    filter_group.button(filterchoose_id).setChecked(True)
    if filterchoose_id == 2:
        lv1.setEnabled(True)
        hv1.setEnabled(True)
        lv1.setText(next(it))
        hv1.setText(next(it))
    else:
        next(it)  # Skip placeholders
        next(it)

    # Window size
    winsizechoose_id = int(next(it))
    winsize_group.button(winsizechoose_id).setChecked(True)
    if winsizechoose_id == 2:
        winsize.setEnabled(True)
        winsize.setText(next(it))
    else:
        next(it)

    # Decomposition
    deco_value = next(it)
    if deco_value != '-999':
        decomp.setChecked(True)
        deco_dt.setEnabled(True)
        deco_dt.setText(deco_value)
    else:
        decomp.setChecked(False)
        deco_dt.setEnabled(False)

    # Spectral analysis
    windows_text = next(it)
    index_win = welch_windows.index(windows_text) if windows_text in welch_windows else 0
    windows_combo.setCurrentIndex(index_win)

    mother_text = next(it)
    index_mother = wavelet_windows.index(mother_text) if mother_text in wavelet_windows else 0
    mother_combo.setCurrentIndex(index_mother)

    # Isotherms
    isotherm_button.setChecked(int(next(it)))
    for idx, (chk, le) in enumerate(zip([iso1, iso2, iso3, iso4], [isoa, isob, isoc, isod])):
        val = next(it)
        if val != '-999':
            chk.setChecked(True)
            le.setEnabled(True)
            le.setText(val)
        else:
            chk.setChecked(False)
            le.setEnabled(False)
            le.clear()

    # Comparisons
    compa1_val = next(it)
    compa2_val = next(it)
    if compa1_val in comparison_options:
        compa1_combo.setCurrentText(compa1_val)
    if compa2_val in comparison_options:
        compa2_combo.setCurrentText(compa2_val)

    # Sensors
    sensor_button.setChecked(int(next(it)))
    for chk, le in zip([sen1, sen2, sen3, sen4], [sena, senb, senc, send]):
        val = next(it)
        if val != '-999':
            chk.setChecked(True)
            le.setEnabled(True)
            le.setText(val)
        else:
            chk.setChecked(False)
            le.setEnabled(False)
            le.clear()

    # Smooth
    smooth_button.setChecked(int(next(it)))

    # Folder path
    folder_path = next(it)

    # Bathymetry
    try:
        bathy_checked = int(next(it))
        bat_on.setChecked(bool(bathy_checked))
        if bathy_checked:
            path_bathy = next(it).strip()
            bathy_file_btn.setEnabled(True)
        else:
            path_bathy = ''
            bathy_file_btn.setEnabled(False)
    except StopIteration:
        pass
    
    # Additional Parameters
    additional_lines = []
    for line in it:  
        line = line.strip()
        if line != "-999":
            additional_lines.append(line)
    additional_params.setPlainText('\n'.join(additional_lines))

def window_destroy():
    window.close()

def OpenUrl(url):
    webbrowser.open(url)

def AboutCallBack():
    QMessageBox.information(window, "About", " Interwave Analyzer - Version 2 (2026) \n Copyright (C) 2019 Rafael de Carvalho Bueno \n All rights reserved \n \n Developed by \n Rafael de Carvalho Bueno \n https://buenorc.github.io/ \n\n Improvements and betterments by \n Andreas Lorke \n Tobias Bleninger \n\n Report problems and improvements to email adresss below \n decarvalhobueno@gmail.com\n \n for mor information, see: \n https://buenorc.github.io/pages/interwave.html \n ")

def DashboardCallBack():

    file_path, _ = QFileDialog.getOpenFileName(
        window,
        "Select dash_data.npz",
        "",
        "NPZ Files (*.npz)"
    )

    if not file_path:
        return

    try:
        output_path = os.path.dirname(file_path)

        # Choose random free-ish port
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

    except Exception as e:
        QMessageBox.critical(
            window,
            "Dashboard Error",
            f"Could not start dashboard:\n{str(e)}"
        )

# ---------------------- Application setup ------------------------------------
vHeigh = 600

app = QApplication(sys.argv)
window = QMainWindow()
window.setGeometry(100, 100, 530, vHeigh)


def resource_path(relative_path):
    try:
        # PyInstaller stores files in _MEIPASS
        base_path = sys._MEIPASS
    except AttributeError:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

icon_path = resource_path("assets/iwcon.ico")
window.setWindowIcon(QIcon(icon_path))

window.setWindowTitle("Interwave Analyzer")

# ---------------------------- Menu bar ---------------------------------------
menubar = window.menuBar()

# File menu
filemenu = menubar.addMenu("File")

open_action = QAction("Open", window)
open_action.triggered.connect(open_button)
filemenu.addAction(open_action)

save_action = QAction("Save as...", window)
save_action.triggered.connect(lambda: save_settings(0))
filemenu.addAction(save_action)

filemenu.addSeparator()

exit_action = QAction("Exit", window)
exit_action.triggered.connect(window_destroy)
filemenu.addAction(exit_action)


# Help menu
helpmenu = menubar.addMenu("Help")

manual_action = QAction("Manual", window)
site = "https://buenorc.github.io/pages/interwave/user-manual.html"
manual_action.triggered.connect(lambda: OpenUrl(site))
helpmenu.addAction(manual_action)

dashboard_action = QAction("Dashboard", window)
dashboard_action.triggered.connect(DashboardCallBack)
helpmenu.addAction(dashboard_action)

about_action = QAction("About", window)
about_action.triggered.connect(AboutCallBack)
helpmenu.addAction(about_action)


central_widget = QWidget()
window.setCentralWidget(central_widget)

layout = QVBoxLayout()
central_widget.setLayout(layout)

tabControl = QTabWidget()
tabControl.setMaximumHeight(vHeigh) 
layout.addWidget(tabControl)


# ---------------------- Tab 1 - Input data -----------------------------------

tab1 = QWidget()
tab1_layout = QGridLayout()
tab1.setLayout(tab1_layout)
tabControl.addTab(tab1, "Input data")


# Input Files
openFile = "Open File"

def temperature_function():
    global path_temp
    path_temp, _ = QFileDialog.getOpenFileName(None, openFile, "", "TEM Files (*.tem)")
    print("Selected Temperature File:", path_temp)

def meteo_function():
    global path_meteo
    path_meteo, _ = QFileDialog.getOpenFileName(None, openFile, "", "MET Files (*.met)")
    print("Selected Meteo File:", path_meteo)


# Header
header_label = QLabel("<b>Input files</b>")
header_label.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
tab1_layout.addWidget(header_label, 0, 0, 1, 2, alignment=Qt.AlignLeft)

# Temperature
tab1_layout.addWidget(QLabel("Temperature data:"), 1, 0, Qt.AlignLeft)
temp_btn = QPushButton("Open File")
temp_btn.clicked.connect(temperature_function)
tab1_layout.addWidget(temp_btn, 1, 1, Qt.AlignLeft)

# Meteorological data
tab1_layout.addWidget(QLabel("Meteorological data:"), 2, 0, Qt.AlignLeft)
meteo_btn = QPushButton("Open File")
meteo_btn.clicked.connect(meteo_function)
tab1_layout.addWidget(meteo_btn, 2, 1, Qt.AlignLeft)

# Solar Radiation checkbox
radiation_on = QCheckBox("Solar Radiation")
radiation_on.setChecked(True)
tab1_layout.addWidget(radiation_on, 3, 0, Qt.AlignLeft)

# Height of wind measurements
tab1_layout.addWidget(QLabel("Height of the wind measurements (meters):"), 4, 0, Qt.AlignLeft)
height_wind = QLineEdit()
height_wind.setText("10")
tab1_layout.addWidget(height_wind, 4, 1)

# Wind direction contribution
tab1_layout.addWidget(QLabel("Wind direction contribution (°):"), 5, 0, Qt.AlignLeft)
contri_wind = QSpinBox()
contri_wind.setRange(0, 180)
contri_wind.setValue(20)
tab1_layout.addWidget(contri_wind, 5, 1)

# Latitude
tab1_layout.addWidget(QLabel("Latitude (°):"), 6, 0, Qt.AlignLeft)
latitude = QSpinBox()
latitude.setRange(-180, 180)
latitude.setValue(0)
tab1_layout.addWidget(latitude, 6, 1)

separator = QFrame()
separator.setFrameShape(QFrame.HLine)
separator.setFrameShadow(QFrame.Sunken)
tab1_layout.addWidget(separator, 7, 0, 1, 2)


# Basin length
header_label = QLabel("<b>Basin length and bathymetry:</b>")
header_label.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
tab1_layout.addWidget(header_label, 8, 0, 1, 2, alignment=Qt.AlignLeft)

# Radio buttons
typechoose_group = QButtonGroup(tab1)  # To group the radio buttons
rad1 = QRadioButton("Uniform depth")
rad2 = QRadioButton("File")
typechoose_group.addButton(rad1, 1)
typechoose_group.addButton(rad2, 2)

tab1_layout.addWidget(rad1, 10, 0, alignment=Qt.AlignLeft)
tab1_layout.addWidget(rad2, 11, 0, alignment=Qt.AlignLeft)

# Line edit for manual input
tab1_layout.addWidget(QLabel("Basin length (meters):"), 10, 0, alignment=Qt.AlignRight)
vari_len = QLineEdit()
vari_len.setEnabled(True)  # Initially enabled
tab1_layout.addWidget(vari_len, 10, 1)

# File selection button
file_len = QPushButton("Open File")
file_len.setEnabled(False)  # Initially disabled
tab1_layout.addWidget(file_len, 11, 1, alignment=Qt.AlignLeft)

# Functions 
def fetch_function():
    global path_fetch
    path_fetch, _ = QFileDialog.getOpenFileName(None, "Open File", "", "LEN Files (*.len)")
    print("Selected Basin File:", path_fetch)

file_len.clicked.connect(fetch_function)

def selected_len():
    id = typechoose_group.checkedId()
    if id == 2:  # File
        file_len.setEnabled(True)
        vari_len.setEnabled(False)
    elif id == 1:  # Uniform
        file_len.setEnabled(False)
        vari_len.setEnabled(True)

rad1.toggled.connect(selected_len)
rad2.toggled.connect(selected_len)

separator2 = QFrame()
separator2.setFrameShape(QFrame.HLine)
separator2.setFrameShadow(QFrame.Sunken)
tab1_layout.addWidget(separator2, 12, 0, 1, 2)  # row 12, spans 2 columns

# Sensor level
tab1_layout.addWidget(QLabel("Sensor level:"), 13, 0, alignment=Qt.AlignLeft)

senso_btn = QPushButton("Open File")
tab1_layout.addWidget(senso_btn, 13, 1, alignment=Qt.AlignLeft)

def senso_function():
    global path_senso
    path_senso, _ = QFileDialog.getOpenFileName(None, "Open File", "", "SEN Files (*.sen)")
    print("Selected Sensor File:", path_senso)

senso_btn.clicked.connect(senso_function)



#Horizontal separator
separator_level = QFrame()
separator_level.setFrameShape(QFrame.HLine)
separator_level.setFrameShadow(QFrame.Sunken)
tab1_layout.addWidget(separator_level, 15, 0, 1, 2)


# Water level data
header_level = QLabel("<b>Water level data:</b>")
header_level.setMaximumHeight(25)
tab1_layout.addWidget(header_level, 16, 0, 1, 2, alignment=Qt.AlignLeft)

# Radio buttons for input type 
level_group = QButtonGroup(tab1)  # group the radio buttons
rad_level_uniform = QRadioButton("Uniform")
rad_level_file = QRadioButton("File")
level_group.addButton(rad_level_uniform, 1)
level_group.addButton(rad_level_file, 2)

tab1_layout.addWidget(rad_level_uniform, 17, 0, alignment=Qt.AlignLeft)
tab1_layout.addWidget(rad_level_file, 18, 0, alignment=Qt.AlignLeft)

# Manual input line edit
unif_level = QLineEdit()
unif_level.setEnabled(True) 
unif_level.setMaximumHeight(25)
tab1_layout.addWidget(QLabel("Water level (meters):"), 17, 0, alignment=Qt.AlignRight)
tab1_layout.addWidget(unif_level, 17, 1)

# File selection button
file_level = QPushButton("Open File")
file_level.setEnabled(False)  # initially disabled
tab1_layout.addWidget(file_level, 18, 1, alignment=Qt.AlignLeft)

# Functions
def level_function():
    global path_level
    path_level, _ = QFileDialog.getOpenFileName(None, "Open File", "", "NIV Files (*.niv)")
    print("Selected Water Level File:", path_level)

file_level.clicked.connect(level_function)

def selected_level():
    id = level_group.checkedId()
    if id == 2:  # File
        file_level.setEnabled(True)
        unif_level.setEnabled(False)
    elif id == 1:  # Uniform
        file_level.setEnabled(False)
        unif_level.setEnabled(True)

rad_level_uniform.toggled.connect(selected_level)
rad_level_file.toggled.connect(selected_level)

# Default selection
rad_level_uniform.setChecked(True)
selected_level()

# Reference level 
tab1_layout.addWidget(QLabel("Level of reference (m)"), 14, 0, alignment=Qt.AlignLeft)

reference_level = QLineEdit()
reference_level.setText("0")
reference_level.setMaximumHeight(25)  # optional: keep row compact
tab1_layout.addWidget(reference_level, 14, 1)


# Horizontal separator
separator_bathy = QFrame()
separator_bathy.setFrameShape(QFrame.HLine)
separator_bathy.setFrameShadow(QFrame.Sunken)
tab1_layout.addWidget(separator_bathy, 23, 0, 1, 2)  # adjust row number


# Bathymetry checkbox
bat_on = QCheckBox("Include Shapefile")
bat_on.setChecked(False)
tab1_layout.addWidget(bat_on, 24, 0, alignment=Qt.AlignLeft)

# Bathymetry file selection button
bathy_file_btn = QPushButton("Open File")
bathy_file_btn.setEnabled(bat_on.isChecked())  # set initial state
tab1_layout.addWidget(bathy_file_btn, 24, 1, alignment=Qt.AlignLeft)

# Function
def bathy_file_function():
    global path_bathy
    path_bathy, _ = QFileDialog.getOpenFileName(None, "Open Bathymetry File", "", "Bathymetry Files (*.bth)")
    print("Selected Bathymetry File:", path_bathy)

bathy_file_btn.clicked.connect(bathy_file_function)

bat_on.toggled.connect(lambda checked: bathy_file_btn.setEnabled(checked))

# -------------- Tab 2 - Spectral analysis and definitions --------------------
tab2 = QWidget()
tab2_layout = QGridLayout()
tab2.setLayout(tab2_layout)
tabControl.addTab(tab2, "Spectral analysis and definitions")


# Density structure header
header_density = QLabel("<b>Density structure</b>")
tab2_layout.addWidget(header_density, 1, 0, 1, 2, alignment=Qt.AlignLeft)

# Metalimnion threshold
tab2_layout.addWidget(QLabel("Metalimnion threshold (kg/m³/m):"), 2, 0, alignment=Qt.AlignLeft)
meta = QLineEdit()
meta.setText("0.1")
tab2_layout.addWidget(meta, 2, 1)

# Separator
separator1 = QFrame()
separator1.setFrameShape(QFrame.HLine)
separator1.setFrameShadow(QFrame.Sunken)
tab2_layout.addWidget(separator1, 3, 0, 1, 2)

# Band pass filter header
header_filter = QLabel("<b>Band pass filter (Cutoff Frequency)</b>")
tab2_layout.addWidget(header_filter, 4, 0, 1, 2, alignment=Qt.AlignLeft)

# Radio buttons for filter choice 
filter_group = QButtonGroup(tab2)
rad_filter_model = QRadioButton("Defined by the Internal wave model")
rad_filter_manual = QRadioButton("Defined manually:")
filter_group.addButton(rad_filter_model, 1)
filter_group.addButton(rad_filter_manual, 2)

tab2_layout.addWidget(rad_filter_model, 5, 0, alignment=Qt.AlignLeft)
tab2_layout.addWidget(rad_filter_manual, 6, 0, alignment=Qt.AlignLeft)

# Default selection
rad_filter_model.setChecked(True)

# High-frequency period
tab2_layout.addWidget(QLabel("High-frequency band period (hour)"), 7, 0, alignment=Qt.AlignLeft)
lv1 = QLineEdit()
lv1.setEnabled(False)
tab2_layout.addWidget(lv1, 7, 1)

# Low-frequency period
tab2_layout.addWidget(QLabel("Low-frequency band period (hour)"), 8, 0, alignment=Qt.AlignLeft)
hv1 = QLineEdit()
hv1.setEnabled(False)
tab2_layout.addWidget(hv1, 8, 1)

# Function
def selected_filter():
    filter_type = filter_group.checkedId()
    if filter_type == 1:  # Internal wave model
        lv1.setEnabled(False)
        hv1.setEnabled(False)
    elif filter_type == 2:  # Manual
        lv1.setEnabled(True)
        hv1.setEnabled(True)

# Connect radio buttons to the function
rad_filter_model.toggled.connect(selected_filter)
rad_filter_manual.toggled.connect(selected_filter)

# Horizontal separator
separator2 = QFrame()
separator2.setFrameShape(QFrame.HLine)
separator2.setFrameShadow(QFrame.Sunken)
tab2_layout.addWidget(separator2, 11, 0, 1, 2)


# Spectral analysis parameter header
header_spectral = QLabel("<b>Spectral analysis parameter</b>")
tab2_layout.addWidget(header_spectral, 12, 0, 1, 2, alignment=Qt.AlignLeft)

# Fourier window function
tab2_layout.addWidget(QLabel("Fourier window function:"), 13, 0, alignment=Qt.AlignLeft)
welch_windows = ["Hamming", "Hann", "Blackman", "Flattop"]
windows_combo = QComboBox()
windows_combo.addItems(welch_windows)
tab2_layout.addWidget(windows_combo, 13, 1)

# Wavelet mother function
tab2_layout.addWidget(QLabel("Wavelet mother function:"), 14, 0, alignment=Qt.AlignLeft)
wavelet_windows = ["Morlet", "Paul", "DOG"]
mother_combo = QComboBox()
mother_combo.addItems(wavelet_windows)
tab2_layout.addWidget(mother_combo, 14, 1)

# Horizontal separator
separator3 = QFrame()
separator3.setFrameShape(QFrame.HLine)
separator3.setFrameShadow(QFrame.Sunken)
tab2_layout.addWidget(separator3, 15, 0, 1, 2)


# Window size (Welch-averaging)
header_window = QLabel("<b>Window Size (Welch-averaging)</b>")
tab2_layout.addWidget(header_window, 16, 0, 1, 2, alignment=Qt.AlignLeft)

# Radio buttons for window size
winsize_group = QButtonGroup(tab2)
rad_win_model = QRadioButton("Defined by the Internal wave model (10 times the V1H1 mode)")
rad_win_manual = QRadioButton("Defined manually (days):")
winsize_group.addButton(rad_win_model, 1)
winsize_group.addButton(rad_win_manual, 2)

tab2_layout.addWidget(rad_win_model, 17, 0, alignment=Qt.AlignLeft)
tab2_layout.addWidget(rad_win_manual, 18, 0, alignment=Qt.AlignLeft)

# Default selection
rad_win_model.setChecked(True)

# Entry for manual window size
winsize = QLineEdit()
winsize.setEnabled(False)
tab2_layout.addWidget(winsize, 18, 1)

# Function 
def selected_windowsize():
    if winsize_group.checkedId() == 1:
        winsize.setEnabled(False)
    else:
        winsize.setEnabled(True)

rad_win_model.toggled.connect(selected_windowsize)
rad_win_manual.toggled.connect(selected_windowsize)


# Horizontal separator
separator3 = QFrame()
separator3.setFrameShape(QFrame.HLine)
separator3.setFrameShadow(QFrame.Sunken)
tab2_layout.addWidget(separator3, 19, 0, 1, 2)


# Decomposition model
header_decomp = QLabel("<b>Decomposition model</b>")
tab2_layout.addWidget(header_decomp, 21, 0, 1, 2, alignment=Qt.AlignLeft)

# Temporal resolution checkbox
decomp = QCheckBox("Temporal resolution (min)")
tab2_layout.addWidget(decomp, 23, 0, alignment=Qt.AlignLeft)

# Entry for temporal resolution
deco_dt = QLineEdit()
deco_dt.setEnabled(False)
tab2_layout.addWidget(deco_dt, 23, 1)

# --- Function to enable/disable deco_dt based on checkbox ---
def selected_decomp(checked):
    deco_dt.setEnabled(checked)

decomp.toggled.connect(selected_decomp)


# --------------- Tab 3 - Isotherms and sensors definitios --------------------

tab3 = QWidget()
tab3_layout = QGridLayout()
tab3.setLayout(tab3_layout)
tabControl.addTab(tab3, "Isotherm Analysis")

# Isotherms analysis
header_iso = QLabel("<b>Isotherms analysis</b>")
tab3_layout.addWidget(header_iso, 0, 0, 1, 2, alignment=Qt.AlignLeft)

isotherm_button = QCheckBox("Enable isotherms analysis")
isotherm_button.setChecked(True)
tab3_layout.addWidget(isotherm_button, 2, 0, alignment=Qt.AlignLeft)

# Individual isotherms
isoa = QLineEdit()
isob = QLineEdit()
isoc = QLineEdit()
isod = QLineEdit()
isoa.setEnabled(False)
isob.setEnabled(False)
isoc.setEnabled(False)
isod.setEnabled(False)

iso1 = QCheckBox("Isotherm 1 (°C)")
iso2 = QCheckBox("Isotherm 2 (°C)")
iso3 = QCheckBox("Isotherm 3 (°C)")
iso4 = QCheckBox("Isotherm 4 (°C)")

tab3_layout.addWidget(iso1, 3, 0, alignment=Qt.AlignLeft)
tab3_layout.addWidget(isoa, 3, 1)
tab3_layout.addWidget(iso2, 4, 0, alignment=Qt.AlignLeft)
tab3_layout.addWidget(isob, 4, 1)
tab3_layout.addWidget(iso3, 5, 0, alignment=Qt.AlignLeft)
tab3_layout.addWidget(isoc, 5, 1)
tab3_layout.addWidget(iso4, 6, 0, alignment=Qt.AlignLeft)
tab3_layout.addWidget(isod, 6, 1)

# Enable/disable individual isotherms
def selected_isotherms(checked):
    iso1.setEnabled(checked)
    iso2.setEnabled(checked)
    iso3.setEnabled(checked)
    iso4.setEnabled(checked)
isotherm_button.toggled.connect(selected_isotherms)

iso1.toggled.connect(lambda checked: isoa.setEnabled(checked))
iso2.toggled.connect(lambda checked: isob.setEnabled(checked))
iso3.toggled.connect(lambda checked: isoc.setEnabled(checked))
iso4.toggled.connect(lambda checked: isod.setEnabled(checked))

# Comparison between isotherms
tab3_layout.addWidget(QLabel("Comparison between isotherms:"), 7, 0, alignment=Qt.AlignLeft)
comparison_options = ["None", "Iso. 1-2", "Iso. 1-3", "Iso. 1-4", "Iso. 2-3", "Iso. 2-4", "Iso. 3-4"]

compa1_combo = QComboBox()
compa1_combo.addItems(comparison_options)
tab3_layout.addWidget(QLabel("Comparison 1:"), 8, 0, alignment=Qt.AlignLeft)
tab3_layout.addWidget(compa1_combo, 8, 1)

compa2_combo = QComboBox()
compa2_combo.addItems(comparison_options)
tab3_layout.addWidget(QLabel("Comparison 2:"), 9, 0, alignment=Qt.AlignLeft)
tab3_layout.addWidget(compa2_combo, 9, 1)

separator_sen = QFrame()
separator_sen.setFrameShape(QFrame.HLine)
separator_sen.setFrameShadow(QFrame.Sunken)
tab3_layout.addWidget(separator_sen, 10, 0, 1, 2)


# Sensor analysis
header_sen = QLabel("<b>Sensor analysis</b>")
tab3_layout.addWidget(header_sen, 11, 0, 1, 2, alignment=Qt.AlignLeft)

sensor_button = QCheckBox("Enable sensor analysis")
sensor_button.setChecked(True)
tab3_layout.addWidget(sensor_button, 12, 0, alignment=Qt.AlignLeft)

# Individual sensors
sena = QLineEdit("1")
senb = QLineEdit("2")
senc = QLineEdit("3")
send = QLineEdit()
send.setEnabled(False)

sen1 = QCheckBox("Sensor 1 (-)")
sen2 = QCheckBox("Sensor 2 (-)")
sen3 = QCheckBox("Sensor 3 (-)")
sen4 = QCheckBox("Sensor 4 (-)")

sen1.setChecked(True)
sen2.setChecked(True)
sen3.setChecked(True)

tab3_layout.addWidget(sen1, 13, 0, alignment=Qt.AlignLeft)
tab3_layout.addWidget(sena, 13, 1)
tab3_layout.addWidget(sen2, 14, 0, alignment=Qt.AlignLeft)
tab3_layout.addWidget(senb, 14, 1)
tab3_layout.addWidget(sen3, 15, 0, alignment=Qt.AlignLeft)
tab3_layout.addWidget(senc, 15, 1)
tab3_layout.addWidget(sen4, 16, 0, alignment=Qt.AlignLeft)
tab3_layout.addWidget(send, 16, 1)

# Enable/disable sensors
def selected_sensor(checked):
    sen1.setEnabled(checked)
    sen2.setEnabled(checked)
    sen3.setEnabled(checked)
    sen4.setEnabled(checked)
sensor_button.toggled.connect(selected_sensor)

sen1.toggled.connect(lambda checked: sena.setEnabled(checked))
sen2.toggled.connect(lambda checked: senb.setEnabled(checked))
sen3.toggled.connect(lambda checked: senc.setEnabled(checked))
sen4.toggled.connect(lambda checked: send.setEnabled(checked))

separator_smooth = QFrame()
separator_smooth.setFrameShape(QFrame.HLine)
separator_smooth.setFrameShadow(QFrame.Sunken)
tab3_layout.addWidget(separator_smooth, 17, 0, 1, 2)


# Smoothing data
header_smooth = QLabel("<b>Smoothing data</b>")
tab3_layout.addWidget(header_smooth, 18, 0, 1, 2, alignment=Qt.AlignLeft)


smooth_button = QCheckBox("Smooth results (recommended for large data analysis)")
smooth_button.setChecked(False)
tab3_layout.addWidget(smooth_button, 19, 0, alignment=Qt.AlignLeft)


# ------------------ Tab 4 - Output and Run -----------------------------------
tab4 = QWidget()
tab4_layout = QGridLayout()
tab4.setLayout(tab4_layout)
tabControl.addTab(tab4, "Output")

# Output settings
header_output = QLabel("<b>Output settings</b>")
tab4_layout.addWidget(header_output, 0, 0, 1, 2, alignment=Qt.AlignLeft)

# Output folder
tab4_layout.addWidget(QLabel("Output folder:"), 2, 0, alignment=Qt.AlignLeft)

folder_path = ""
def output_folder():
    global folder_path
    folder_path = QFileDialog.getExistingDirectory()
    if folder_path:
        folder_line.setText(folder_path)

folder_line = QLineEdit()
folder_line.setReadOnly(True)
tab4_layout.addWidget(folder_line, 2, 1)

folder_btn = QPushButton("Select Folder")
folder_btn.clicked.connect(output_folder)
tab4_layout.addWidget(folder_btn, 2, 2, alignment=Qt.AlignLeft)

# Images quality
tab4_layout.addWidget(QLabel("Output images quality (dpi):"), 3, 0, alignment=Qt.AlignLeft)

dpi_spin = QSpinBox()
dpi_spin.setRange(50, 1000)
dpi_spin.setValue(800)
tab4_layout.addWidget(dpi_spin, 3, 1, alignment=Qt.AlignLeft)


# Horizontal separator
separator_output = QFrame()
separator_output.setFrameShape(QFrame.HLine)
separator_output.setFrameShadow(QFrame.Sunken)
tab4_layout.addWidget(separator_output, 4, 0, 1, 2)


header_output = QLabel("<b>Additional Parameters:</b>")
tab4_layout.addWidget(header_output, 6, 0, 1, 2, alignment=Qt.AlignLeft)

# Additional Parameters
additional_params = QTextEdit()
additional_params.setPlaceholderText("Enter additional parameters, one per line...\nCheck User Manual for correct codes")
tab4_layout.addWidget(additional_params, 7, 0, 1, 2)  # span 2 columns


# Run button
run_button = QPushButton("Run")
run_button.setFixedHeight(40)
run_button.setFixedWidth(150)
tab4_layout.addWidget(run_button, 9, 1, alignment=Qt.AlignLeft)


def save_settings(temporary):
    
    # try:
    #     from PySide6.QtWidgets import QFileDialog
    # except ImportError:
    #     from PyQt5.QtWidgets import QFileDialog

    global path_temp, path_meteo, path_senso, path_fetch, path_level, path_bathy, folder_path

    if temporary == 'temporary':
        data_file = 'temporary.txt'
        space = '\n'
    else:
        f = QFileDialog.getSaveFileName(None, "Save Settings", "", "SET Files (*.set)")
        data_file = str(f[0])
        if not data_file:
            return
        space = '\n'

    with open(data_file, 'w') as data:
        # Tab 0: Paths 
        data.write(str(path_temp)+space)
        data.write(str(path_meteo)+space)
        data.write(str(path_senso)+space)

        # Tab 1: Wind / Meteorology
        data.write(str(height_wind.text())+'\n')
        data.write(str(int(radiation_on.isChecked()))+'\n')
        data.write(str(reference_level.text())+'\n')
        data.write(str(contri_wind.value())+'\n')
        data.write(str(latitude.value())+'\n')
        data.write(str(dpi_spin.value())+'\n')

        # Tab 1: Basin length 
        typechoose_id = typechoose_group.checkedId()
        data.write(str(typechoose_id)+'\n')
        if typechoose_id == 1:
            data.write(str(vari_len.text())+'\n')
        elif typechoose_id == 2:
            data.write(str(path_fetch)+space)

        # Tab 1: Levels 
        typelevel_id = level_group.checkedId()
        data.write(str(typelevel_id)+'\n')
        if typelevel_id == 1:
            data.write(str(unif_level.text())+'\n')
        elif typelevel_id == 2:
            data.write(str(path_level)+space)

        # Tab 2: Density 
        data.write(str(meta.text())+'\n')

        # Tab 2: Band pass filter 
        filter_id = filter_group.checkedId()
        data.write(str(filter_id)+'\n')
        if filter_id == 2:
            data.write(str(lv1.text())+'\n')
            data.write(str(hv1.text())+'\n')
        else:
            data.write(str(-999)+'\n')
            data.write(str(-999)+'\n')

        # Tab 2: Window size 
        winsize_id = winsize_group.checkedId()
        data.write(str(winsize_id)+'\n')
        if winsize_id == 2:
            data.write(str(winsize.text())+'\n')
        else:
            data.write(str(-999)+'\n')

        # Tab 2: Decomposition
        if decomp.isChecked():
            data.write(str(deco_dt.text())+'\n')
        else:
            data.write(str(-999)+'\n')

        # Tab 2: Spectral analysis 
        data.write(str(windows_combo.currentText())+'\n')
        data.write(str(mother_combo.currentText())+'\n')

        #  Tab 3: Isotherms 
        data.write(str(int(isotherm_button.isChecked()))+'\n')
        data.write(str(isoa.text() if iso1.isChecked() else -999)+'\n')
        data.write(str(isob.text() if iso2.isChecked() else -999)+'\n')
        data.write(str(isoc.text() if iso3.isChecked() else -999)+'\n')
        data.write(str(isod.text() if iso4.isChecked() else -999)+'\n')
        data.write(str(compa1_combo.currentText())+'\n')
        data.write(str(compa2_combo.currentText())+'\n')

        # Tab 3: Sensors
        data.write(str(int(sensor_button.isChecked()))+'\n')
        data.write(str(sena.text() if sen1.isChecked() else -999)+'\n')
        data.write(str(senb.text() if sen2.isChecked() else -999)+'\n')
        data.write(str(senc.text() if sen3.isChecked() else -999)+'\n')
        data.write(str(send.text() if sen4.isChecked() else -999)+'\n')

        # Tab 3: Smoothing
        data.write(str(int(smooth_button.isChecked()))+'\n')

        # Tab 4: Output folder 
        data.write(str(folder_path)+space)

        # Tab 5: Bathymetry 
        data.write(str(int(bat_on.isChecked()))+'\n')
        if bat_on.isChecked() and 'path_bathy' in globals() and path_bathy:
            data.write(str(path_bathy)+space)
        else:
            data.write(str(-999)+space)
        
        # Additional Parameters 
        params_text = additional_params.toPlainText().strip()
        if params_text:

            for line in params_text.splitlines():
                data.write(line + '\n')
        else:
            data.write("-999\n")  

def export_data():
    import iwback as iwm
    print("Temporary file being created...")
    save_settings('temporary') 
    print("Temporary file saved. Running backend...")
    iwm.main()

run_button.clicked.connect(export_data)

if __name__ == "__main__":
    window.show()
    app.exec()