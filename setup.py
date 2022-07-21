#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PlottingOnGenome",
    version="0.0.1",
    author="Katarina Erbstein",
    author_email="katarina.erbstein@embl.de",
    description="Pipeline to plot sequences against genome",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://git.embl.de/grp-zimmermann-kogadeeva/PlottingOnGenome.git",
    python_requires=">=3.8",
    install_requires=[
        "matplotlib", "numpy", "biopython", "dna_features_viewer", "pandas", "tqdm"
    ],
    packages=setuptools.find_packages()
)
