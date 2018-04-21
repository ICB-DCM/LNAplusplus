#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is an LNA++ application example with a simple birth-death process 
involving a constitutively active DNA, and stochastic production and 
degradation of mRNA and protein:

Model of birth-death process:
==============================
         k_m
DNA      -->  DNA+mRNA
         g_m
mRNA     -->  0
         k_p
mRNA     -->  mRNA + protein
         g_p
protein  -->  0

A more detailed model description is provided in 
documentation/Documentation.pdf.

This example is used for performance testing of LNA++. We 
1) show how to perform exact simulation of the stochastic process and use LNA++
to compute the mean and temporal auto-covariance of the proteins, 
2) compute sensitivities of the mean protein copy number with respect to the model 
parameters and compare the computed results to a finite difference approximation,
3) perform multi-start optimization in order to compare parameter inference 
with and without the analytical sensitivities and 
4) use LNA++ to evaluate the likelihood landscape for this model.
"""

# Open plot windows or plot to file?
plotToFile = False
import sys
if len(sys.argv) > 1 and sys.argv[1] == '--headless':
    plotToFile = True
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

def plot(plotToFile, filename): 
    if plotToFile:
        plt.savefig(filename)
    else:
        plt.show()

# add LNA++ python module to search path
import os
scriptDir = os.path.dirname(os.path.realpath(__file__))
from sys import path
path.append(scriptDir + '/../python')


# Create the Python module for the birth / death system
from LNA import generateLNA
generateLNA(scriptDir + '/BirthDeath/BirthDeath.xml', 'BirthDeath', computeSS='BOTH')

# Import the generated module (adapt path first)
path.append(scriptDir + '/modules')
import BirthDeathLNA

import numpy as np

# model parameters
parameterNames =  ['$k_m$','$k_p$','$g_m$','$g_p$']
Theta = [20, 25, 10, 1]

# time span 0.0, 0.1, ... 10.0
tspan = np.arange(0, 10, 0.1).tolist()

# Simulate: IC = steady state
MRE,Var = BirthDeathLNA.LNA(Theta,tspan)

# plot results:
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
fig.suptitle('Simulation: IC = Steady state')
# mean of mRNA
ax1.plot(tspan, MRE[0])
ax1.set_xlabel('Time')
ax1.set_ylabel('mRNA')
ax1.set_title('Macroscopic mean of mRNA')
# auto-covariance of mRNA
ax2.matshow(Var[0, 0, :, :], extent=[np.min(tspan), np.max(tspan), np.max(tspan), np.min(tspan)])
ax2.set_xlabel('Time')
ax2.set_ylabel('Time')
ax2.set_title('Autocovariance of mRNA')
# mean of protein
ax3.plot(MRE[1])
ax3.set_xlabel('Time')
ax3.set_ylabel('Protein')
ax3.set_title('Macroscopic mean of protein')
# auto-covariance of protein
ax4.matshow(Var[1, 1, :, :], extent=[np.min(tspan), np.max(tspan), np.max(tspan), np.min(tspan)])
ax4.set_xlabel('Time')
ax4.set_ylabel('Time')
ax4.set_title('Autocovariance of protein')
plot(plotToFile, 'BirthDeath1.pdf')



# Simulate: IC = steady state; observable = mRNA and protein
# Both species are observed (= mRNA,protein)
ObsIndex = [1, 2]
# Variance of measurement noise
VarNoise = [10.0,50.0]

# solve LNA and compute measured distribution
MRE, Var = BirthDeathLNA.LNA(Theta, tspan, merr =VarNoise,obsVar=ObsIndex)

# plot results:
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
fig.suptitle('Simulation: IC = steady state; observable = mRNA and protein')
# mean of mRNA
ax1.plot(tspan, MRE[0])
ax1.set_xlabel('Time')
ax1.set_ylabel('mRNA')
ax1.set_title('Macroscopic mean of mRNA')
# auto-covariance of mRNA
ax2.matshow(Var[0, 0, :, :], extent=[np.min(tspan), np.max(tspan), np.max(tspan), np.min(tspan)])
ax2.set_xlabel('Time')
ax2.set_ylabel('Time')
ax2.set_title('Autocovariance of mRNA')
# mean of protein
ax3.plot(MRE[1])
ax3.set_xlabel('Time')
ax3.set_ylabel('Protein')
ax3.set_title('Macroscopic mean of protein')
# auto-covariance of protein
ax4.matshow(Var[1, 1, :, :], extent=[np.min(tspan), np.max(tspan), np.max(tspan), np.min(tspan)])
ax4.set_xlabel('Time')
ax4.set_ylabel('Time')
ax4.set_title('Autocovariance of protein')
plot(plotToFile, 'BirthDeath2.pdf')



# Simulate: IC = steady state; observable = protein
# observable is second species (= protein)
ObsIndex = [2]
# Variance of measurement noise
VarNoise = [50.0]

# solve LNA and compute measured distribution
MRE, Var = BirthDeathLNA.LNA(Theta, tspan, merr=VarNoise,obsVar=ObsIndex)

# plot results:
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Simulation: IC = steady state; observable = protein')
# mean of protein
ax1.plot(tspan, MRE[0])
ax1.set_xlabel('Time')
ax1.set_ylabel('Protein')
ax1.set_title('Macroscopic mean of protein')
# auto-covariance of protein
ax2.matshow(Var[0, 0, :, :], extent=[np.min(tspan), np.max(tspan), np.max(tspan), np.min(tspan)])
ax2.set_xlabel('Time')
ax2.set_ylabel('Time')
ax2.set_title('Autocovariance of protein')
plot(plotToFile, 'BirthDeath3.pdf')



# Simulate: IC = no steady state; observable = protein
# observable is second species (= protein)
ObsIndex = [2]
# Variance of measurement noise
VarNoise = [0.0]
# initial values (E(mRNA),E(Protein))
MRE0 = [2, 200];
# initial co-variances (cov[mRNA,mRNA],cov(mRNA,Protein),cov(Protein,Protein))
Var0 = [0, 0, 0];

# solve LNA and compute measured distribution
MRE, Var = BirthDeathLNA.LNA(Theta, tspan, Y0=MRE0, V0=Var0, merr=VarNoise, obsVar=ObsIndex)

# plot results:
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Simulation: IC = no steady state; observable = protein')
# mean of protein
ax1.plot(tspan, MRE[0])
ax1.set_xlabel('Time')
ax1.set_ylabel('Protein')
ax1.set_title('Macroscopic mean of protein')
# auto-covariance of protein
ax2.matshow(Var[0, 0, :, :], extent=[np.min(tspan), np.max(tspan), np.max(tspan), np.min(tspan)])
ax2.set_xlabel('Time')
ax2.set_ylabel('Time')
ax2.set_title('Autocovariance of protein')
plot(plotToFile, 'BirthDeath4.pdf')



# Simulate: IC = no steady state; observable = protein; sensitivity = 1st order

# solve LNA and compute measured distribution
MRE, Var, Sens_MRE, Sens_Var = BirthDeathLNA.LNA(Theta, tspan, Y0=MRE0, V0=Var0, merr=VarNoise, obsVar=ObsIndex, computeSens=True)

# plot results:
fig, axs = plt.subplots(2, len(Theta))
fig.suptitle('Simulate: IC = no steady state; observable = protein; sensitivity = 1st order')
for i in range(len(Theta)):
    (ax1, ax2) = axs[0][i], axs[1][i]
    # Sensitivity of mean of protein
    ax1.plot(tspan, Sens_MRE[0, i, :])
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Sensitivity of mean')
    ax1.set_title(parameterNames[i])
    # Sensitivity of auto-covariance of protein
    ax2.matshow(Sens_Var[0, 0, :, :, i], extent=[np.min(tspan), np.max(tspan), np.max(tspan), np.min(tspan)])
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Sensitivity of autocovariance')
    ax2.set_title(parameterNames[i])
plot(plotToFile, 'BirthDeath5.pdf')



# Simulate: IC = no steady state; observable = protein; sensitivity = 1st & 2nd order

# solve LNA and compute measured distribution
MRE, Var, Sens_MRE, Sens_Var,Sens2_MRE,Sens2_Var = BirthDeathLNA.LNA(Theta, tspan, Y0=MRE0, V0=Var0, merr=VarNoise, obsVar=ObsIndex, computeSens2=True)

# plot results:
fig, axs = plt.subplots(len(Theta), len(Theta))
fig.suptitle('Simulation of mean: IC = no steady state; observable = protein; sensitivity = 1st & 2nd order')
for i in range(len(Theta)):
    for j in range(len(Theta)):
        ax1 = axs[i][j]
        # Sensitivity of mean of protein
        ax1.plot(tspan, Sens2_MRE[0, i, j, :])
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Sensitivity of mean')
        ax1.set_title('(%s, %s)' % (parameterNames[i], parameterNames[j]))
plot(plotToFile, 'BirthDeath6.pdf')

fig, axs = plt.subplots(len(Theta), len(Theta))
fig.suptitle('Simulation of autocovariance: IC = no steady state; observable = protein; sensitivity = 1st & 2nd order')
for i in range(len(Theta)):
    for j in range(len(Theta)):
        ax1 = axs[i][j]
        # Sensitivity of mean of protein
        ax1.matshow(Sens2_Var[0, 0, i, j, :, :], extent=[np.min(tspan), np.max(tspan), np.max(tspan), np.min(tspan)])
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Time')
        ax1.set_title('(%s, %s)' % (parameterNames[i], parameterNames[j]))
plot(plotToFile, 'BirthDeath7.pdf')

