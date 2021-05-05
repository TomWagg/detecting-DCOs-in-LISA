# Predictions for detecting BHNS and other double compact objects binaries with LISA

Hello there! This repository contains all of the code associated with with Wagg et al. (2021), both for the simulations so that you can reproduce the results, but also for the paper so that you can reproduce all of the plots. If you haven't read the paper yet, I recommend that you go and take a look before trying to see what's going on in the code.

## Table of contents
- [Required dependencies](#required-dependencies)
    - [Installing using Conda (recommended)](#installing-using-conda-recommended)
    - [Installing using pip](#installing-using-pip)
- [Required data](#required-data)
- [Using this repository](#using-this-repository)
    - [How to reproduce/adapt paper figures](how-to-reproduceadapt-paper-figures)
    - [How to reproduce results/run a new simulation](how-to-reproduce-resultsrun-a-new-simulation)
- [Repository map](#repository-map)

<hr style="border-width:5px">

## Required dependencies
This code makes heavy use of [LEGWORK](https://legwork.readthedocs.io/en/latest/), a Python package for determining the detectability of stellar-mass LISA sources that I wrote with Katie Breivik. It allows you evolve the orbits of binary sources, measure their strain and SNR as well as visualise the results.

We also use some more well known Python packages and list them below
- `numpy`
- `astropy`
- `scipy`
- `seaborn`
- `h5py`
- `matplotlib`
- `jupyter`
- `ipython`

### Installing using Conda (recommended)
You can install these packages using `conda` and create an environment for working with this code. Run the following to do so
```
conda create -n LISA_dcos numpy astropy scipy seaborn h5py matplotlib legwork jupyter ipython
conda activate LISA_dcos
```

### Installing using pip
If you don't have Anaconda then alternatively you can use pip to install everything by running
```
pip install numpy astropy scipy seaborn h5py matplotlib legwork jupyter ipython
```
## Required data
All data used in the paper is stored [here](https://zenodo.org/record/4699713) on Zenodo. You'll need to download all of this data for reproducing the figures. If you also want to run a new simulation then you'll need the data from Broekgaarden et al. (2021) which you can download [here](https://zenodo.org/record/4574727).

## Using this repository
This repository is set up so that it is very easy not only to reproduce every result and figure in the paper but also adapt the work for future studies. With this in mind, let's split into two sections depending on whether you just want to visualise the results differently or whether you want to produce new results entirely.

### How to reproduce/adapt paper figures

### How to reproduce results/run a new simulation

## Repository map
Not sure where to look? This map should help point you to the right folder!

```
detecting-DCOs-with-LISA
│
└── paper
│   │   This folder contains everything to do with the paper
│   │    
│   │   .xlsx files are excel files with the comparison of previous studies
│   │
│   └─── figures
│   │    │   A collection of every figure used in the paper
│   │    
│   └─── tex
│        │   All of the LaTeX files for the paper (split up by section)
│    
│   
└─── simulation
    │
    └─── data
    │    │   This is the folder where most of the code assumes the data from Zenodo is stored
    │    │
    └─── notebooks
    │    │   Jupyter notebooks that reproduce every figure in the paper
    │    │
    └─── output
    │    │   Where the simulation puts all of its output on completion
    │    │
    └─── slurm
    │    │   SLURM commands for running cluster jobs on the Harvard cluster
    │    │
    └─── src
         │   Main simulation code
```