# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 15:42:47 2016

@author: fjustin
"""

import re
from parse_sbml_stoichiometry import parse_file
from sympy import Matrix
import ipdb

def SBML2StoichProp(sbmlFile):
    try:
        f = open(sbmlFile)
    except IOError:
        print("Invalid file specified")
        return
    
    [species, reactions, parameters, reactionFormulas, S] = parse_file(f)
    S = Matrix(S)
    # convert reactionFormulas to function handle
    
    for k in range(len(reactionFormulas)):
        formula = reactionFormulas[k]
        for k2 in range(len(species)):
            formula = re.sub(r'\b%s\b'%species[k2], 'phi['+str(k2)+']', formula )
        for k2 in range(len(parameters)):
            formula = re.sub(r'\b%s\b'%parameters[k2], 'Theta['+str(k2)+']', formula)
        reactionFormulas[k] = formula

    lambdaStr = 'lambda phi, t, Theta: ['+','.join(reactionFormulas)+']'
    fHandle = eval(lambdaStr)
    
    return S, species, parameters, fHandle
    