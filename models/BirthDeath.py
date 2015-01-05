# -*- coding: utf-8 -*-
"""
Created on Sat Aug  2 10:59:39 2014

@author: justinfeigelman
"""
# add python module to search path
from sys import path
path += ['../python']

import LNA
from LNA import *

# name of module
model = 'BirthDeath';

# define the stoichiometric matrix
S = Matrix([[1, -1, 0, 0], \
   [0, 0, 1, -1]])

# define symbolic variables   
# mRNA(m) and protein(p)
phi = symbols('m,p', real=True)

# reaction constants
Theta = symbols('k_m,k_p,g_m,g_p', real=True)

# reaction fluxes
f = lambda phi,t,Theta: \
   [Theta[0], Theta[2]*phi[0], Theta[1]*phi[0], Theta[3]*phi[1]]

# generate the LNA components
tups = generateLNAComponents(model, S, f, phi, Theta)

# generate the C code
npar = len(Theta) # number of parameters
compileLNA(model, S, tups, npar)

############### run some simulations and show results ############### 
path += ['modules']
import BirthDeathLNA
from numpy import arange


# define initial conditions
Y0 = [2, 200]
V0 = [2, 3, 150]

# model parameters
Theta = [20, 25, 10, 1]

# only observing protein (variable 2)
obsVar = 2
# no measurement error
merr = 0.0

# time span 0.0, 0.1, ... 10.0
# use numpy arange() function
tspan = arange(0,10,0.1).tolist()

# define the stoichiometric matrix
S = [[1,-1,0,0],
	[0,0,1,-1]]

# compute the mean and temporal auto-covariance matrix assuming steady-state mean and variance at t=0
MRE,Sigma = BirthDeathLNA.LNA(Theta,tspan, obsVar=obsVar)


# compute the mean and temporal auto-covariance matrix specifying the initial conditions for mean and variance
MRE,Sigma = BirthDeathLNA.LNA(Theta,tspan, Y0, V0)

# compute the mean and temporal auto-covariance matrix specifying the initial conditions for mean and variance
# compute also the first order sensitivities
MRE,Sigma,dMRE,dSigma  = BirthDeathLNA.LNA(Theta,tspan, Y0, V0, computeSens=True, obsVar=obsVar)

# compute the mean and temporal auto-covariance matrix specifying the initial conditions for mean and variance
# compute also the first and second order sensitivities
MRE,Sigma,dMRE,dSigma,d2MRE,d2Sigma = BirthDeathLNA.LNA(Theta,tspan, Y0, V0, obsVar=obsVar, computeSens=True, computeSens2=True)

## plots
#import matplotlib
import pylab
from numpy import squeeze

# MRE
pylab.plot(MRE.squeeze())
pylab.show()

# covariance matrix
pylab.matshow(squeeze(Sigma))
pylab.show()

# sensitivity w.r.t. first model parameter
pylab.matshow(squeeze(dSigma[0,0,:,:,1]))
pylab.show()

