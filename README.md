## Interwave Analyzer 
### DOI: 10.5281/zenodo.4264575

The Interwave Analyzer is an open-software designed to investigate the occurrence of basin-scale internal waves in stratified lakes, providing a set of tools to estimate some physical indices to analyze the probability of internal waves activity and lake mixing. The program provides the classification of lake mixing, internal wave patterns, and degeneration of basin-scale internal waves based on well-established theories of internal waves. The Interwave Analyzer provides a powerful, easily accessible, and universal analysis of internal waves from instrumented buoy and meteorological data stations.

### Atual version: 1.01.2

### scientific paper (published with 1.00.3 version):
de Carvalho Bueno, Bleninger, and Lorke. **Internal wave analyzer for thermally stratified lakes**. *Environmental Modelling and Software*.

GuiHub (scientific paper repository): https://github.com/buenorc/espaper.git

### How to run:

To run the Interwave Analyzer’s scripts in a Python interpreter, the system must have:

* Python interpreter;

* Python packages: Numpy, Datetime, Reportlab, Scipy, Nitime, Matplotlib, and Tkinter;

The steps below describes the installation instruction for Anaconda’s users, which is the simplest way to run
the Interwave Analyzer’s codes. We recommend the use of Anaconda distribution since it will automatically
install almost all required additional Python packages. Interwave Analyzer also can be ran through other
Python interpreters, but additional packages installations is required.


To run Interwave Analyzer’s scripts directly in Anaconda interpreter, first download the Anaconda distribution
for Python 3.x.

1) Go to Anaconda website (https://www.anaconda.com/) and find the option for Anaconda distribution;
2) Choose the Python 3.x graphical installer version (note that there are three options for operating system:
Windows, macOS, and Linux);
3) Install the Anaconda interpreter;
4) After the installation, open the Anaconda Prompt (as administrator) and install the following packages
that are used by Interwave Analyzer and are not available in Anaconda:

* Nitime (0.7 or compatible): conda install -c conda-forge nitime
* Reportlab (3.5 or compatible): conda install conda-forge::reportlab

**Attention:** Some users reported that when installing ReportLab with the command conda install -c anaconda reportlab (for python 3.11 or higher), the installed version is outdated. This results in a DeprecationWarning that hinders the smooth operation of the Interwave Analyzer due to issues with the ReportLab library.

**Attention:** If you already have anaconda installed in your computer, make sure that the above packages
are installed and the anaconda version has Python 3.X. If you use another interpreter, make sure that the
following packages are installed in the your python interpreter: Numpy 1.16.3, Datetime 4.0.1, Reportlab
3.5.19, Scipy 1.2.1, Nitime 0.7, Matplotlib 3.1.0, and tk (tkinter) 8.6.8, or compatibles versions.

5) After the installation, go to https://sites.google.com/view/interwaveanalyzer/interwave-analyzer and click on
Code repository, or access directly our repository in https://github.com/buenorc/interwaveanalyzer;
6) Download all files .py available to download, including the raster-graphics file-format 0interwave.png
which is the logo used by the Interwave Analyzer’s report and the icon interwave icon, which is used
as Interwave Analyzer icon on the GUI;
7) Put everything in the same folder and run the script called GUIexecuter.py;
8) A graphical user interface (GUI) should be launched in seconds;
9) We recommend you to download the example files available at Interwave Analyzer’s website on Download
example. For a detailed tutorial using these files, see section 5.

**Attention:** For more information (user manual, how to obtain the pre-compiled version, FAQ, team, etc), please visit: https://sites.google.com/view/interwaveanalyzer/interwave-analyzer
