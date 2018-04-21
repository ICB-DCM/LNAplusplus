#!/usr/bin/env python3
"""
This example constructs the LNA for the model of DNA self-regulation
introduced by Wang et al. 2010.

Model:
======
R1: DNA + P2  -->  DNA.P2
R2: DNA.P2    -->  DNA + P2
R3: DNA       -->  DNA + P 
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
LNA.generateLNA(scriptDir + '/Wang2010/Wang2010.xml', 'Wang2010', 'NONE')

# Extend path for LNA++-created modules and import
sys.path.append(LNA.lnaModulesDir)
import Wang2010LNA

## Define parameters, initial conditions and simulation time
Theta = [0.1, 0.7, 0.35, 0.3, 0.1, 0.9, 0.2, 0.1]
speciesNames = ['DNA', 'DNAP2', 'RNA', 'P', 'P2']
numSpecies = len(speciesNames)
MRE0  = [20.0, 0.0, 0.0, 0.0, 0.0]
Var0  = LNA.toLinear(np.zeros((numSpecies, numSpecies)))
tspan = np.linspace(0, 150, 100).tolist()

# Simulate model
[MRE,Var] = Wang2010LNA.LNA(Theta, tspan, merr=0.0, Y0=MRE0, V0=Var0)

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

