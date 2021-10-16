<p align="center">
    <img width="100%", src="paper/figures/readme_header.png"/>
</p>

<p align="center">
    This repository is a companion to Wagg et al. (2021), a paper investigating the Galactic population of LISA detectable BH and NS binaries. We discuss the expected detection rate, variations when using different underlying binary physics assumptions, prospects for identifying the sources of signals and matching them to electromagnetic counterparts.
</p>
<p align="center">
    This repository contains all of the code associated with with Wagg et al. (2021), both for the simulations so that you can reproduce the results, but also for the paper so that you can reproduce all of the plots. If you've come this without reading the paper yet, I recommend that you go and take a look before trying to see what's going on in the code so it makes more sense.
</p>

## Table of contents
- [Required dependencies](#required-dependencies)
    - [Installing using Conda (recommended)](#installing-using-conda-recommended)
    - [Installing using pip](#installing-using-pip)
- [Required data](#required-data)
- [Using this repository](#using-this-repository)
    - [How to reproduce/adapt paper figures](#how-to-reproduceadapt-paper-figures)
    - [How to reproduce results/run a new simulation](#how-to-reproduce-resultsrun-a-new-simulation)
    - [Milky Way Model code](#milky-way-model-code)
    - [Colour Schemes](#colour_schemes)
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
All data used in the paper is stored [here](https://zenodo.org/record/4699713) on Zenodo. You'll need to download all of this data for reproducing the figures and be sure to follow the instructions on Zenodo for placing the data in the right place. If you also want to run a new simulation then you'll need the data from Broekgaarden et al. (2021) which you can download [here](https://zenodo.org/record/4574727).

## Using this repository
This repository is set up so that it is very easy not only to reproduce every result and figure in the paper but also adapt the work for future studies. With this in mind, let's split into two sections depending on whether you just want to visualise the results differently or whether you want to produce new results entirely.

### How to reproduce/adapt figures from the paper
Here's a quick guide on how you can reproduce any figure from the paper, as well as how to adapt those figures for your own use. First, for the majority of the figures in the paper you'll need to download the data by following the instructions on [Zenodo](https://zenodo.org/record/4699713). Once you've downloaded the data and put it in the right folder (`simulation/data`), you're ready to start making figures!

You're welcome to peruse the notebooks directly here on GitHub if you want to take a look at everything...but if you're looking for a particular figure then I recommend that you check the caption in the paper. There will be a small book symbol that links directly to the notebook that creates it. For example, let's consider [Figure 9](paper/figures/fig9_detection_rates.pdf) and say you wanted to only include the NSNS panel instead over DCO type. In this case, you need to open [this notebook](paper/figure_notebooks/detections.ipynb), which is linked in Figure 9's caption, and edit the code to remove the loop over DCO types so that you plot only the NSNSs instead.

In this same way you can pick any figure and follow the link in the caption to find a related notebook that will explain how to make the figure with commented code that you can edit.

Please remember to cite Wagg et al. (2021) should you use any figure/adapted figure in any public setting (papers/talks etc.). Thanks!
### How to reproduce results/run a new simulation
TODO

### Milky Way Model code
Want to use our Milky Way model in your work without coding it yourself? You're in luck, the `simulate_mw()` function in [simulation/src/galaxy.py](simulation/src/galaxy.py) will do exactly this. If you're interested in exploring the galaxy model in more detail I recommend you check out the [Galaxy Creation Station](paper/figure_notebooks/galaxy_creation_station.ipynb).

### Colour schemes
In case you're interested in using my colour schemes for the physics variations or formation channels:
- The colours used for each physics variation (as well as labels/descriptions) are contained in [simulation/src/variations.py](simulation/src/variations.py)
- The colours for the formation channels are defined in the formation channels notebook which can be found in [paper/figure_notebooks/formation_channels.ipynb](paper/figure_notebooks/formation_channels.ipynb)

## Repository map
Not sure where to look? This map should help point you to the right folder!

```
detecting-DCOs-with-LISA
│
└── paper
│   │   This folder contains everything to do with the paper
│   │    
│   │   .xlsx files are excel files with the comparison of previous studies│    
│   └─── figure_notebooks
│   |    |   Jupyter notebooks that reproduce every figure in the paper
│   │
│   └─── figures
│   │    │   A collection of every figure used in the paper
│   │    └─── extra_figures
│   │         │ Some figures that didn't make it into the paper but that are still interesting!
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
    └─── postprocessing_notebooks
    │    │   Jupyter notebooks that transform the direct simulation output into more usable formats
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
