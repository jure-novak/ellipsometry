# Ellipsometry - experimental data fitting

## Overview

The aim of this project is numerical analysis of a [data set from ellipsometric measurements](https://en.wikipedia.org/wiki/Ellipsometry)
 (psi and delta) in Python. 
 
 ## Features
Scripts find optimal fit of a set of parameters of an optical model of the sample using bounded non-linear fitting. Currently it is suitable for analysing measurements taken at one wavelength and multiple angles of incidence. It can be easily modified to perform analysis on data gathered at multiple wavelengths (so-called VASE ellipsometry).

Package *lmfit* used for fitting enables quick change of optimization algorithms, e.g. *leastsq*, *nelder*, *lbfgsb*.

## Setup & Installation
The scripts in this project use scientific computing libraries from open science platform Anaconda.
It can be downloaded [here](https://www.continuum.io/downloads).

**Package *[lmfit](https://lmfit.github.io/lmfit-py/intro.html)* is used for bounded non-linear fitting.**
It can be installed using command line: `conda install -c conda-forge lmfit`.
More info about *lmfit* package installation can be found [here](https://lmfit.github.io/lmfit-py/installation.html).

User should first run [*decaying sine wave* example](https://lmfit.github.io/lmfit-py/intro.html) from *lmfit* documentation to see if everything is working correctly.

## Usage
For data analysis use either plain Python notebook *Two layer model.py* or Jupyter notebook *Two layer model.ipynb*. The Jupyter notebook already contains some sample plots.


