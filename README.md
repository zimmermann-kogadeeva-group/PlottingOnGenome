
# PlottingOnGenome

## Install

BLAST is required. You can follow the instructions to download and install
BLAST
[here](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html). 

Install with pip:
```
pip install git+https://git.embl.de/grp-zimmermann-kogadeeva/PlottingOnGenome.git
```
or with pipx (**recommended**):
```
pipx install git+https://git.embl.de/grp-zimmermann-kogadeeva/PlottingOnGenome.git
```
Whether using pip or pipx, python3 and git will need to be installed.

Building exe in Windows using PyInstaller package:
```
pyinstaller -F -w -n plotting_on_genome -p .\plotting_on_genome\ --hidden-import Bio.SearchIO.BlastIO  .\pyinstaller_entry.py
```

