# Predictions for detecting BHNS and other double compact objects binaries with LISA

Hello there! This repository contains all of the code associated with with Wagg et al. (2021), both for the simulations so that you can reproduce the results, but also for the paper so that you can reproduce all of the plots. If you haven't read the paper yet, I recommend that you go and take a look before trying to see what's going on in the code.

This code makes heavy use of [LEGWORK](https://legwork.readthedocs.io/en/latest/), a python package for determining the detectability of stellar-mass LISA sources that I wrote with Katie Breivik. You can [install LEGWORK](https://legwork.readthedocs.io/en/latest/install.html) using pip and will need it to run a lot of the code.

## How to use this repository

## Repository Map

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