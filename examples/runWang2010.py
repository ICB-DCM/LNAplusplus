#!/usr/bin/env python3
"""
This example constructs the LNA for the model of DNA self-regulation
introduced by Wang et al. (2010) in "Parameter inference for discretely 
observed stochastic kinetic models using stochastic gradient descent."
BMC Syst. Biol., 4, 99.


Model:
======
R1: DNA + P2  -->  DNA.P2
R2: DNA.P2    -->  DNA + P2
R3: DNA       -->  DNA + RNA 
R4: RNA       -->  0
R5: P + P     -->  P2
R6: P2        -->  P + P
R7: RNA       -->  RNA + P 
R8: P         -->  0
"""

# add LNA++ Python module to search path
import os
import sys
scriptDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(scriptDir + '/../python')

import numpy as np
import time
import LNA
import sympy

# Create the Python module for the model
LNA.generateLNA(scriptDir + '/Wang2010.xml', 'Wang2010', 'NONE')
""" Note: For this model the symbolic calculation of the steady state fails as
the resulting equations are already to complex."""

# Extend path for LNA++-created modules and import
sys.path.append(LNA.lnaModulesDir)
import Wang2010LNA

## Define parameters, initial conditions and simulation time
Theta = [0.1, 0.7, 0.35, 0.3, 0.1, 0.9, 0.2, 0.1]
speciesNames = ['DNA', 'DNAP2', 'RNA', 'P', 'P2']
numSpecies = len(speciesNames)
MRE0  = [20.0, 1.0, 1.0, 1.0, 1.0]
Var0  = LNA.toLinear(np.zeros((numSpecies, numSpecies)))
tspan = np.linspace(0, 150, 100).tolist()

# Simulate model
[MRE,Var,sMRE,sVar,s2MRE,s2Var] = Wang2010LNA.LNA(Theta, tspan, merr=0.0, Y0=MRE0, V0=Var0, computeSens2=True)

# Plot results
import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 3)
fig.suptitle('Simulation')

for k in range(numSpecies):
    ax = fig.axes[k]
    # Assign mean and variance for the k-th species
    m = MRE[k, :]
    s = np.sqrt(Var[k, k, :, :].squeeze().diagonal())
    ax.fill_between(tspan + list(reversed(tspan)),
                np.append(m + s, (m-s)[::-1]), 
                facecolor='b', alpha=0.25, label='+/- std')
    ax.plot(tspan, m, 'b', label='mean')
    ax.set_xlabel('Time')
    ax.set_ylabel('Abundance')
    ax.set_title(speciesNames[k])
    ax.legend()
fig.axes[-1].axis('off')
plt.show()


# Test of cross-species sensitivities
i = 4
j = 4
eps_theta = 1e-3;
Theta_per = Theta 
Theta_per[i] *= 1 + eps_theta
[MRE_per,Var_per,sMRE_per,sVar_per,s2MRE_per,s2Var_per] = Wang2010LNA.LNA(Theta_per, tspan, merr=0.0, Y0=MRE0, V0=Var0, computeSens2=True)

k1 = 50-1;
k2 = 100-1;

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax1.imshow((Var_per[:, :, k1, k2] - Var[:, :, k1, k2]) / eps_theta)
ax1.set_title('finite differences')
ax2.imshow(sVar[:, :, k1, k2, i])
ax2.set_title('analytical sensitivities')
ax3.imshow((Var_per[:, :, k1, k2] - Var[:, :, k1, k2]) / eps_theta - sVar[:, :, k1, k2, i])
ax3.set_title('error')
plt.show()

# 2nd order test - diagonal
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax1.imshow(np.squeeze((sVar_per[1, 1, :, :, j] - sVar[1, 1, :, :, j]) / eps_theta))
ax1.set_title('finite differences')
ax2.imshow(np.squeeze(s2Var[1, 1, i, j, :, :]))
ax2.set_title('analytical sensitivities')
ax3.imshow(np.squeeze((sVar_per[1, 1, :, :, j] - sVar[1, 1, :, :, j]) / eps_theta) - np.squeeze(s2Var[1, 1, i, j, :, :]))
ax3.set_title('error')
plt.show()

# 2nd order test - off-diagonal
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax1.imshow((sVar_per[:, :, k1, k2, j] - sVar[:, :, k1, k2, j]) / eps_theta)
ax1.set_title('finite differences')
ax2.imshow(s2Var[:, :, i, j, k1, k2])
ax2.set_title('analytical sensitivities')
ax3.imshow((sVar_per[:, :, k1, k2, j] - sVar[:, :, k1, k2, j]) / eps_theta - s2Var[:, :, i, j, k1, k2])
ax3.set_title('error')
plt.show()
