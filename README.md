
# PlottingOnGenome

Package to find and plot inserts alongside the genome or contigs.

BLAST is required. You can follow the instructions to download and install
BLAST
[here](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html). 

## Installation

```
pip install git+https://git.embl.org/grp-zimmermann-kogadeeva/PlottingOnGenome.git
```

## Running

### Docker

To run the streamlit app with docker:
```
docker run -it --rm -p 8501:8501 docker://registry.git.embl.org/grp-zimmermann-kogadeeva/plottingongenome:0.7.3
```

### Apptainer

To run the streamlit app with apptainer:
```
apptainer exec --no-home --writable-tmpfs --workdir /app docker://registry.git.embl.org/grp-zimmermann-kogadeeva/plottingongenome:0.7.3 streamlit run /app/app/app.py --server.headless true
```

