
# PlottingOnGenome

Package to find and plot inserts alongside the genome or contigs.

BLAST is required. You can follow the instructions to download and install
BLAST
[here](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html). 

## Download 

Click
[here](https://oc.embl.de/index.php/s/Ahn7fKtJywiBDAX/download)
to download the Windows executable. 

If using this option, there is no need to follow any of the steps in the next
section.

## Install

Install with pip:
```
pip install git+https://git.embl.de/grp-zimmermann-kogadeeva/PlottingOnGenome.git
```
or with pipx (**recommended**):
```
pipx install git+https://git.embl.de/grp-zimmermann-kogadeeva/PlottingOnGenome.git
```
Whether using pip or pipx, python3 and git will need to be installed.

### PyInstaler

Building exe in Windows using PyInstaller package:
```
pyinstaller -F -w -n plotting_on_genome -p .\plotting_on_genome\ --hidden-import Bio.SearchIO.BlastIO  .\pyinstaller_entry.py
```

## Usage

This package can be used through GUI, CLI or through importing. 

All interfaces require the same set of inputs:
- Search term - the search_term is required for downloading the genome of your
  interest from ncbi.
- Email - your email address is required to download genomes from NCBI. 
- Forward suffix - suffix in your fasta files indicating forward sequences
  (default: `_F`)
- Reverse suffix - suffix in your fasta files indicating reverse sequences
  (default: `_R`)
- Filter threshold - threshold value for rejecting sequences with low coverage
  when determining whether a forward-reverse pair of sequences can be
  considered an insert
- Retmax - maximum number of records to be downloaded from NCBI
- Sequences file - location of your fasta file
- Output folder - location where to save all output files such as plots and the
  table of information on the inserts found

The package plots two types of plots: a circular graph of the whole genome or
linear graphs of individual inserts along with the genes found in the same
location on the genome.

### Example


