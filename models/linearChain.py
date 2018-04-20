#!/usr/bin/env python3
"""
This example constructs the LNA for a linear reaction chain of varying
length. This is used to assess the scalability of the symbolic
construction as well as the numerical simulation without sensitivities,
with 1st order sensitivities and with 2nd order sensitivities.

Model of linear reaction chain:
===============================
X_{i} -> X_{i+1} for i = 1 to N-1
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

def linearChainStoichiometricMatrix(N):
    """
    Generate the symbolic stoichiometric matrix for a linear reaction chain with N species: X_{i} <-> X_{i+1} for i = 1 to N-1
    """
    
    S = sympy.zeros(N, 2*N)
    
    # birth/death
    S[0, 0] = 1
    S[-1, -1] = -1

    # forward reactions
    for i in range(N-1):
        S[i:(i+2), i+1] = [-1, 1]

    # backward reactions
    for i in range(N-1):
        S[i:(i+2), N+i+1] = [1, -1]
    
    return S    

def linearChainFluxes(phi, t, Theta): 
    """
    Return the reaction fluxes for the linear chain as a function of phi, t and Theta:
    """
    res =  [None] * len(phi) * 2
    res[0] = Theta[0] # birth, 0-order
    for i in range(1, len(phi)): # birth
        res[i] = Theta[i] * phi[i-1]
    for i in range(len(phi)): # death
        res[len(phi) + i] = Theta[len(phi)+i] * phi[i]
    #print(res)
    return res
    
# Define evaluation parameters
R = 10 # number of simulation runs
N_max = 6 # maximum chain length (the minimum is 3)

# Initialize matrices for storing computation and run time
computeTimeMat = np.full([N_max-2, 1], np.nan) # time for model construction and compilation
runTimeMat0    = np.full([N_max-2, R], np.nan) # time for simulation w/o sensitivities
runTimeMat1    = np.full([N_max-2, R], np.nan) # time for simulation with 1st order sensitivities
runTimeMat2    = np.full([N_max-2, R], np.nan) # time for simulation with 2nd order sensitivities

# Extend path for LNA++-created modules
sys.path.append(LNA.lnaModulesDir)

tspan = np.arange(0, 2.5, 0.1).tolist()

# estimate compute and run time for models
for N in range(3, N_max + 1):
    phi   = sympy.symbols('phi:%d' % N)     # symbolic variables for macroscopic mean
    Theta = sympy.symbols('Theta:%d' % (2*N)) # symbolic variables for parameters
    sympy.assumptions.assume.global_assumptions.add(sympy.Q.positive(Theta))
    
    # Generate stochiometric matric and flux vector
    S = linearChainStoichiometricMatrix(N);
    
    # Generate model
    timeStart = time.time()
    modelName = 'chain%d' % N
    # Generate the symbolic C code:
    tups = LNA.generateLNAComponents(modelName, S, linearChainFluxes, phi, Theta, 'NONE')
    # Compile the code and generate the Python module:
    LNA.compileLNA(modelName, S, tups, len(Theta))
    computeTimeMat[N-3] = time.time() - timeStart;
    
    # load generated model
    import importlib
    curModule = importlib.import_module("chain%dLNA" % N)
    curLNA = getattr(curModule, 'LNA')
    
    # Run simulations
    print('Running simulations (N=%d)...' % N)
    Var0  = [0] * int(N*(N+1)/2)
    
    for i in range(R):
        # Generate random parameters and initial conditions
        Theta = np.random.uniform(0, 100, 2*N).tolist()
        MRE0  = np.random.uniform(0, 1000, N).tolist()
        
        # Simulate w/o sensitivities
        timeStart2 = time.time()
        # Y0=MRE0, V0=Var0, merr=VarNoise, obsVar=ObsIndex, computeSens=True
        [M,S] = curLNA(Theta, tspan, Y0=MRE0, V0=Var0)
        runTimeMat0[N-3, i] = time.time() - timeStart2        
        print('.')

        # Simulate with 1st order sensitivities
        timeStart2 = time.time()
        [M1,S1,dM,dS] = curLNA(Theta, tspan, Y0=MRE0, V0=Var0, computeSens=True)
        runTimeMat1[N-3, i] = time.time() - timeStart2
        print('.')
        
        # Simulate with 2nd order sensitivities
        timeStart2 = time.time()
        [M2,S2,dM,dS,d2M,d2S] = curLNA(Theta, tspan, Y0=MRE0, V0=Var0, computeSens2=True)
        runTimeMat2[N-3, i] = time.time() - timeStart2
        print('.')
    
# Plot results
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set_title('Computation time for model construction')
ax1.semilogy(range(3, N_max + 1), computeTimeMat, '.-')
ax1.set_xlabel('Network size')
ax1.set_ylabel('Model construction time (s)')

ax2.set_title('Computation time for model simulation');
print(range(3, N_max + 1), runTimeMat0)
ax2.semilogy(range(3, N_max + 1), np.median(runTimeMat0, axis=1), 'k.-', label='w/o sensitivities')
ax2.semilogy(range(3, N_max + 1), np.median(runTimeMat1, axis=1), 'b.-', label='1st order sensitivities')
ax2.semilogy(range(3, N_max + 1), np.median(runTimeMat2, axis=1), 'r.-', label='2nd order sensitivities')
ax2.set_xlabel('Network size')
ax2.set_ylabel('Median run time (s)')
ax2.legend()
plt.show()