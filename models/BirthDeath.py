# -*- coding: utf-8 -*-
"""
Created on Sat Aug  2 10:59:39 2014

@author: justinfeigelman
"""
# add python module to search path
from sys import path
path.append('../python')

from LNA import generateLNA

generateLNA('BirthDeath/BirthDeath.xml', 'BirthDeath', computeSS='BOTH')

############### run some simulations and show results ############### 
path += ['modules']
import BirthDeathLNA
from numpy import arange
import matplotlib.pyplot as plt
from numpy import squeeze


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

# compute the mean and temporal auto-covariance matrix assuming steady-state mean and variance at t=0
MRE,Sigma = BirthDeathLNA.LNA(Theta,tspan, obsVar=obsVar)

# MRE
plt.plot(MRE.squeeze())
plt.show()

# covariance matrix
plt.matshow(squeeze(Sigma))
plt.show()


# compute the mean and temporal auto-covariance matrix specifying the initial conditions for mean and variance
MRE,Sigma = BirthDeathLNA.LNA(Theta,tspan, Y0, V0, obsVar=obsVar)
# MRE of protein
plt.plot(MRE.squeeze())
plt.show()

# covariance matrix
plt.matshow(squeeze(Sigma))
plt.show()


# compute the mean and temporal auto-covariance matrix specifying the initial conditions for mean and variance
# compute also the first order sensitivities
MRE,Sigma,dMRE,dSigma  = BirthDeathLNA.LNA(Theta,tspan, Y0, V0, computeSens=True, obsVar=obsVar)

# sensitivity w.r.t. first model parameter
plt.matshow(squeeze(dSigma[0,0,:,:,1]))
plt.show()


# compute the mean and temporal auto-covariance matrix specifying the initial conditions for mean and variance
# compute also the first and second order sensitivities
MRE,Sigma,dMRE,dSigma,d2MRE,d2Sigma = BirthDeathLNA.LNA(Theta,tspan, Y0, V0, obsVar=obsVar, computeSens=True, computeSens2=True)

# MRE
plt.plot(MRE.squeeze())
plt.show()

# covariance matrix
plt.matshow(squeeze(Sigma))
plt.show()

# sensitivity w.r.t. first model parameter
plt.matshow(squeeze(dSigma[0,0,:,:,1]))
plt.show()

plt.matshow(squeeze(d2Sigma[0,0,3,3,:,:]))
plt.show()


