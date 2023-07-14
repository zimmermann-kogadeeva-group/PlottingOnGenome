#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="plotting_on_genome",
    version="0.0.6",
    author="Bartosz Bartmanski",
    author_email="bartosz.bartmanski@embl.de",
    description="Package for making genome related plots",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://git.embl.de/grp-zimmermann-kogadeeva/PlottingOnGenome",
    python_requires='>=3.9',
    install_requires=[
        "matplotlib",
        "numpy",
        "biopython",
        "pandas",
        "tqdm",
        "dna_features_viewer",
        "PySimpleGUI"
    ],
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            "pog = plotting_on_genome.__cli__:main",
            "pog-gui = plotting_on_genome.gui:main"
        ],
    }
)
