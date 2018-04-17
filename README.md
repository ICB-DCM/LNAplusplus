# LNA++: a Fast C++ Implementation of the Linear Noise Approximation with first- and second-order sensitivities

The linear noise approximation (LNA) provides an approximate description of the statistical moments of stochastic chemical reaction networks, a commonly used modeling paradigm describing the probability distribution of systems of biochemical species in the intracellular environment. Unlike exact formulations, the LNA remains computationally feasible even for CRNs with many reaction events. The tractability of the LNA for approximation of CRNs makes it a common choice for inference of unknown chemical reaction parameters. However, this task is impeded by a lack of suitable inference tools for arbitrary CRN models. In particular, no available tool provides temporal cross-correlations, parameter sensitivities and efficient numerical integration.

LNA++ provides the means for fast derivation and simulation of the LNA including the computation of means, covariances, and temporal cross-covariances of species in the CRN. For efficient parameter estimation and uncertainty analysis, LNA++ implements 1st and 2nd order sensitivity equations. Interfaces are provided for easy integration with Matlab and Python.

More detailed information is available in [Tutorial.docx](Tutorial.docx) and [UserGuide.docx](UserGuide.docx).


## Installation

The contents of the LNA++ distribution zip file should be placed into a separate folder, such as in a folder within the home directory.  The archive contains all the necessary scripts and source code for building LNA++ Python modules or Matlab executable files.  No further building or installation is necessary until the creation of the model executables, as described below.

## Requirements
LNA++ uses the standard C++ libraries, the CVODES Sundials library for numerical integration of ODEs and forward sensitivities, and the Blitz++ library for fast, templatized tensor operations.  


- Compiler
  
  Compilation of C++ source code for use as a standalone executable, or with Matlab or Python requires the use of a suitable compiler such as gcc/g++ on Linux-based systems, or clang on mac OS.

- Blitz++
  
  Blitz++ version 0.9 or greater is required.  It can be obtained via the included source code in the ‘libraries’ folder, or from http://sourceforge.net.  It is also available via macports on mac OS, and the apt repository on Ubuntu.

- CVODES
  Numerical integration is computed using the CVODES Sundials package.  Source code is included in the ‘libraries’ folder, and is also available for download from http://computation.llnl.gov/casc/sundials/main.html.

  Note: CVODES should be compiled using the -fPIC flag to ensure cross-platform compatibility.  This can be enabled by specifying --with-pic to the configure script.

- Python
  
  In order to use LNA++ with Python, you must first install the necessary scientific computing modules, listed below.  Python libraries can be installed using package managers such as apt on Ubuntu, or macports/homebrew on Mac OS, or using the Pip3 Python package manager.

  LNA++ is currently only compatible with Python 3.0 or greater.  All libraries must thus be the Python3-compatible version.

  Note for Ubuntu Users:  the necessary Python packages can be installed using the command `sudo apt-get install libblitz0-dev libblitz0ldbl libblitz-doc python-sympy libsbml5-python`

  - Numpy

    Numpy is a widely used mathematics and numerical library. 

  - Sympy

    Sympy provides a framework for symbolic/algebraic computations, and is necessary for computing the linear noise approximation for arbitrary models.

- Libsbml

   Libsbml is required for translating SBML models for use with LNA++.

- (optional) pylab
  Pylab provides functions for plotting vectors and matrices, and might thus be useful.

- Matlab
  
  In order to use LNA++ with Matlab, the Matlab Symbolic Toolbox is required in order to compute the analytic components of the LNA and its sensitivities.  If the Symbolic Toolbox is not available, LNA++ can still be used in two ways:
    1. Python only.  Use the Python implementation (no Matlab) as described in the remainder of this document.
    2. Mixed Matlab/Python.  Both the Python and Matlab implementations of LNA++ generate C code for the analytical components of the LNA and its sensitivities.  It is therefore possible to generate the C code using the Python scripts, through the generateLNAComponents step, as described below.  The Matlab executable is then generated in the normal way, resuming with the compileLNA step.
