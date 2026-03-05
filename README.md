## Interwave Analyzer V2
### DOI: 10.5281/zenodo.4264575

The **Interwave Analyzer** is an open-source software designed to investigate the occurrence of basin-scale internal waves in stratified lakes and reservoirs. The program provides tools to estimate physical indices related to internal wave activity, lake mixing, and stratification dynamics based on established internal wave theories.

Using temperature time series from instrumented buoys together with meteorological data, the software identifies internal wave patterns, classifies mixing regimes, and evaluates the degeneration of basin-scale internal waves.

**Version 2** introduces a major refactoring of the codebase, improved computational performance, and the implementation of the **Additional Parameters Framework**, which expands the range of physical diagnostics that can be computed from the input data. The new version also incorporates **bathymetric information and basin asymmetry**, allowing improved estimation of key stability parameters such as the **Lake Number**, **Schmidt Stability**, and related indices of **stratification instability and mixing potential**. These developments also enable a more realistic representation of **sloping topography effects**, improving the physical interpretation of internal wave dynamics in natural lake basins.

### Current version
2.260303 (beta testing for Interwave Analyzer Version 2)

### Scientific paper (published with version 1.00.3)

de Carvalho Bueno, Bleninger, and Lorke.  
**Internal wave analyzer for thermally stratified lakes**. *Environmental Modelling and Software*.

Scientific paper repository: https://github.com/buenorc/espaper.git


## Installation
Interwave Analyzer can be executed in two different ways:
* **Python script version**
* **Executable version** (only for Windows users)

More informations about Executable version find on:
https://buenorc.github.io/pages/interwave.html


## Script Version (Python)

The script mode allows the Interwave Analyzer to run directly through a Python interpreter. No programming knowledge is required, executing the main script launches the graphical interface automatically.

### Required environment

To run Interwave Analyzer scripts, the system must have:

* Python 3.8 or higher
* Python packages: NumPy, SciPy, Matplotlib, Dash, Plotly, Tkinter, and Datetime

Most of these libraries are included in standard scientific Python distributions.

## Recommended Environment: Anaconda

We recommend using the **Anaconda Python distribution**, which simplifies installation and dependency management.

Download Anaconda:

https://www.anaconda.com/

### Installation steps

1) Go to the Anaconda website and download the **Anaconda Distribution for Python 3.x**.

2) Choose the graphical installer for your operating system (Windows, macOS, or Linux).

3) Install Anaconda using the recommended default settings.

4) After installation, open the **Anaconda Prompt**.

5) Install the additional packages required by Interwave Analyzer: **conda install -c conda-forge dash plotly**
   
**Attention:**  
If Anaconda is already installed on your computer, make sure that:

* Python 3.x is being used
* the required packages are installed
* the environment is activated before running the software


## Running the Scripts

1) Visit the Interwave Analyzer webpage:  
https://buenorc.github.io/pages/interwave.html

2) Click **Code Repository** to access the GitHub repository.

3) Download all `.py` source files.

4) Make sure the folder **assets** contains the following files: iwlogo.png, iwcon.ico, and style.css
   
5) Place all files in the same directory (with the assets folder inside this directory).

6) Open the Anaconda Prompt (or activate your Python environment), navigate to the project directory, and run: **python iwgui.py**

The graphical user interface should start automatically within a few seconds.

---

## Example Data

For first-time users, we recommend downloading the example datasets available on the Interwave Analyzer website under **Download Example**.

A detailed tutorial demonstrating the workflow using these datasets is provided in the **User Manual**.

---

## More Information

For documentation, downloads, FAQ, and updates, visit:

https://buenorc.github.io/pages/interwave.html


