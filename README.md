[![Build Status](https://travis-ci.org/ICB-DCM/LNAplusplus.svg?branch=master)](https://travis-ci.org/ICB-DCM/LNAplusplus)

# LNA++: a Fast C++ Implementation of the Linear Noise Approximation with first- and second-order sensitivities

The linear noise approximation (LNA) provides an approximate description of the statistical moments of stochastic chemical reaction networks, a commonly used modeling paradigm describing the probability distribution of systems of biochemical species in the intracellular environment. Unlike exact formulations, the LNA remains computationally feasible even for CRNs with many reaction events. The tractability of the LNA for approximation of CRNs makes it a common choice for inference of unknown chemical reaction parameters. However, this task is impeded by a lack of suitable inference tools for arbitrary CRN models. In particular, no available tool provides temporal cross-correlations, parameter sensitivities and efficient numerical integration.

LNA++ provides the means for fast derivation and simulation of the LNA including the computation of means, covariances, and temporal cross-covariances of species in the CRN. For efficient parameter estimation and uncertainty analysis, LNA++ implements 1st and 2nd order sensitivity equations. Interfaces are provided for easy integration with Matlab and Python.

More detailed information on background, requirements, installation and application is provided in [documentation/Documentation.pdf](documentation/Documentation.pdf).


## Requirements
LNA++ uses the standard C++ libraries, the CVODES Sundials library for numerical integration of ODEs and forward sensitivities, and the Blitz++ library for fast, templatized tensor operations.  

- C++ compiler
- Blitz++
- CVODES
- For python interface: Python 3.x with numpy, sympy, libsbml, pylab, ipdb
- For Matlab interface: Matlab (tested with R2017b)

## Getting started
After the successful installation of Blitz++ and CVODES, we suggest to have a look at the exmaples (in the folder ./examples). runBirthDeath provide an overview over the functionality of LNA++. runLinearChain and runWang2010 provide additional examples (with a longer run time).

## Contacts
- Justin Feigelman: JustinFeigelmanPhD@gmail.com
- Jan Hasenauer: jan.hasenauer@helmholtz-muenchen.de
- Carsten Marr: carsten.marr@helmholtz-muenchen.de
