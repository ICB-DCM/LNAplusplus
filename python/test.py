from sys import path
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
MRE,Sigma,dMRE,dSigma,d2MRE,d2Sigma = BirthDeathLNA.LNA(Theta,tspan, Y0, V0, computeSens=True, computeSens2=True)

## plots
import matplotlib
import pylab
from numpy import squeeze

# MRE
pylab.plot(MRE[1,:].squeeze()) # just the protein
pylab.show()

# covariance matrix
pylab.matshow(squeeze(Sigma[1,1,:,:])) # just the protein
pylab.show()

# sensitivity w.r.t. first model parameter
pylab.matshow(squeeze(dSigma[0,0,:,:,1]))
pylab.show()
