# Introduction
[Section-properties](https://github.com/robbievanleeuwen/section-properties) examples.

# Installations

## [Winget](https://learn.microsoft.com/en-us/windows/package-manager/winget/)

## GIT

```
PS C:\> winget install -e --id Git.Git
```

## UV
```
PS C:\> winget install -e --id astral-sh.uv
```

## Python
```
C:\> uv python install 3.12
github\spe [main ≡]> uv venv
Using CPython 3.12.7
Creating virtual environment at: .venv
Activate with: .venv\Scripts\activate
github\spe [main ≡]> .venv\Scripts\activate
```

## [Spyder IDE](https://www.spyder-ide.org/)
Istalling using winget
```
PS C:\> winget install -e --id Spyder.Spyder
```
Tools/Preferences/Python interpreter github/spe/.venv/Scripts/python.exe

## Python packages

Dependencies can be installed using pip
 * spyder-kernels - needed for spyder integration
 * sectionproperties - needed for almost all examples, numba makes it faster after initial run
 * pygltflib - needed for gltf export
 * plyfile - needed for ply export
 * sympy - needed for  box_girder cell in paper_cells.py

Typical command needed after update of python is (uv in front if it is used)
```
pip install spyder-kernels==3.0.* pygltflib sectionproperties[numba] plyfile sympy
```

## Matlab

Recent version with following toolboxes is needed
 * Curve Fitting Toolbox
 * Statistics and Machine Learning Toolbox
 * Symbolic Math Toolbox

## GBTUL
[GBTUL](https://sites.fct.unl.pt/gbt/pages/gbtul)

## Using spe

 * Get code, e.g. by cloning or forking spe repository
```
PS C:\Users\simon> git clone https://github.com/simo-11/spe
```
 * Start spyder. Spyder can be started also from menus and paper_cells.m can be opened from menus
```
PS C:\Users\simon\spe> C:\ProgramData\spyder-6\envs\spyder-runtime\Scripts\spyder.exe paper_cells.py
```
 * Set current directory at top right corner to directory where code was loaded
 * Run first cell 'common config' and cells that you are interested in by pressing <img src="https://github.com/user-attachments/assets/2adfcb73-c4df-421a-88a4-7406afc74e39" width="20" alt="Ctrl+Return"/> or pressing Ctrl+Return
    * If you want to study detailed steps set breakpoints using left mouse button to set breakpoints <img src="https://github.com/user-attachments/assets/a955d6af-7800-42c0-9672-7b544b8a97d0" width="80"/>. Check shortcuts from Debug-menu
and run by pressing <img src="https://github.com/user-attachments/assets/19e3f974-932a-49b2-b7e9-d92f20eb0f6e" width="20"/>
    * Modify parameters and rerun
 * You can also use IPython Console to control execution. Use !h to get help on commands while debugging i.e. on IPdb promt.
```
In [1]: %runcell -n 'common config' C:/Users/simon/spe/paper_cells.py
In [2]: %debugcell -n rectangles C:/Users/simon/spe/paper_cells.py
IPdb [1]: c
```
### Notes on paper_cells.m

#### gbtul
Execution requires GBTUL to be available.
Location of GBT.exe can be configured using --gbt option.
Cell tries to locate it in parent directories.

 
