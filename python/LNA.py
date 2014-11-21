# -*- coding: utf-8 -*-
#     LNA.py
#     (C) Justin Feigelman
#     justin.feigelman@helmholtz-muenchen.de
#     Helmholtz Zentrum Muenchen

# This module provides the functions necessary for generating the Python interface to the LNA++ package

from sympy import symbols, Matrix, solve, sqrt, diag, Symbol, printing
from sympy.matrices import *
import os

def jacobian(x,y):
	"Compute the jacobian of x with respect to y"

    f=x.transpose().reshape(len(x),1)
    return f.jacobian(y)
    
def triu(X):
	"Return the upper triangular portion of X as a list"
    if not isinstance(X, Matrix):
        X=Matrix(X)
    return sum([Matrix(X)[i,i:].tolist()[0] for i in range(Matrix(X).rows)],[])
    
def generateLNAComponents(modelName, S, reactionFlux, phi, Theta):
	"Generate all of the symbolic components necessary for the LNA method"
    t = symbols('t', real=True)
    Nvar        = S.rows

    phi0 = symbols('phi0:%d'%Nvar) # symbolic state variables
    Theta=Theta+phi0 # initial conditions are additional parameters for sensitivities
    
    F = Matrix(reactionFlux(phi, t, Theta))
    J           = F.jacobian(phi)
    dFdTheta    = F.jacobian(Theta)
    d2fdTheta2   = jacobian(dFdTheta,Theta)

    "sensitivities of A"
    A           = S*J
    dAdTheta    = jacobian(A,Theta)
    dAdPhi      = jacobian(A,phi)
    d2AdTheta2  = jacobian(dAdTheta, Theta)
    d2AdPhi2    = jacobian(dAdPhi, Theta)
    
    "sensitivities of E"
    E           = S*diag(*[sqrt(i[0]) for i in F.tolist()])
    dEdTheta    = jacobian(E,Theta)
    d2EdTheta2  = jacobian(dEdTheta, Theta)
    dEdPhi      = jacobian(E,phi)
    d2EdPhi2    = jacobian(dEdPhi, phi)
       
    "compute initial steady state and sensitivities"
    Y0 = solve( S*F, phi) 
    S0 = jacobian(Matrix(Y0.values()), Theta)
    S20 = jacobian(S0, Theta)    
    
    "compute steady state variance"
    V0, systemJacobian, MI = solveSS_var(A,E,F,S,phi,Theta)
    
    "SS variance sensitivities"
    SV0  = jacobian(V0, Theta)
    S2V0 = jacobian(SV0, Theta)

    "subsitute steady state values"
    subsPairs =  zip(Y0.keys(), Y0.values())
    V0 = V0.subs(subsPairs)
    SV0 = SV0.subs(subsPairs)
    S2V0 = S2V0.subs(subsPairs)

    
    "mixed derivatives"    
    d2AdPhidTheta = jacobian(dAdPhi, Theta);
    d2AdThetadPhi = jacobian(dAdTheta, phi);
    
    d2EdPhidTheta = jacobian(dEdPhi, Theta);
    d2EdThetadPhi = jacobian(dEdTheta, phi);
    
	# objects dictionary for autogenerating C code
    objs = {'reactionFlux':F,
            'J':J,
            'dFdTheta':dFdTheta,
            'd2fdTheta2':d2fdTheta2,
            'Afunc':A,
            'dAdTheta':dAdTheta,
            'dAdPhi':dAdPhi,
            'd2AdTheta2':d2AdTheta2,
            'd2AdPhi2':d2AdPhi2,
            'Efunc':E,
            'dEdTheta':dEdTheta,
            'd2EdTheta2':d2EdTheta2,
            'dEdPhi':dEdPhi,
            'd2EdPhi2':d2EdPhi2,
            'systemJacobian':systemJacobian,  
            'd2AdPhidTheta':d2AdPhidTheta,
            'd2AdThetadPhi':d2AdThetadPhi,
            'd2EdPhidTheta':d2EdPhidTheta,
            'd2EdThetadPhi':d2EdThetadPhi,
    }   
    
    "depends on phi, t, Theta"
    arg1 = [('phi',phi), ('t',t), ('Theta',Theta)]
    
    "depends on Theta"
    arg2 = [('Theta', Theta)]
    
    objs2 = {'S0':S0,
             'S20':S20,
            'SV0':SV0,
            'S2V0':S2V0,
            'V0':V0,
            'Y0':Matrix([Y0[p] for p in phi])}

    "depends on phi, t, Theta and gamma"
    gamma = symbols('gamma', real=True)
    arg3 = [('phi',phi), ('t',t), ('Theta',Theta), ('gamma',gamma)]
            
    objs3 = {'MI':MI}
    
    tups = [(objs, arg1), (objs2,arg2), (objs3,arg3)]
    return tups

def solveSS_var(A,E,F,S,phi,Theta): 
	'''solve for the steady state variance and components for the preconditioner'''

    N = len(phi)
    V = [[Symbol('V%d%d' % (j,i), real=True) for i in range(N)] for j in range(N)]
    dVdt = A*Matrix(V) + Matrix(V)*A.T + E*E.T
    
    V0 = solve(dVdt, sum(V,[]))
    if len(V0)==0:
        warning('Could not solve for steady state variance.  Setting initial variance to zero.')
        V0 = zeros(N)    
        
    VV = triu(V)
    V0 = Matrix([V0[x] for x in VV])
        
    "Solve for System Jacobian"
    Phi = [[Symbol('Phi%d%d' % (j,i), real=True) for i in range(N)] for j in range(N)]
    Phi2 = Matrix(Phi).reshape(1,N**2).tolist()[0]
    sysVar = list(phi)+VV+Phi2
    RHS = sum((S*F).tolist(),[]) + sum([dVdt[i,i:].tolist()[0] for i in range(dVdt.rows)],[]) + (A*Matrix(Phi)).reshape(1,N*N).tolist()[0]
    systemJacobian = jacobian(Matrix(RHS), sysVar)

    "preconditioner"
    gamma = Symbol('gamma',real=True)
    M = eye(systemJacobian.rows) - gamma*systemJacobian
    MI = M.inv()
    return (V0, systemJacobian, MI)    
    
def convertToC(model, S, tups, phi, t, Theta):
    '''generate the C code, for mex file executable in MATLAB'''

    # create output directory if necessary
    if not os.path.isdir('models/%s/C/' % model):
        os.mkdir('models/%s/C/' % model)
        
    for objs,args in tups:        
        for x in objs:               
            H,C=genCcode(objs[x], x, args)
            f=open("models/%s/C/%s.h" % (model,x), 'w')
            f.writelines(H)
            f.close()
            
            f=open("models/%s/C/%s.c" % (model,x), 'w')
            f.writelines(C)
            f.close()
        
	# generate MODEL_DEF.h
    f = open("models/%s/C/MODEL_DEF.h" % model, 'w')
    f.writelines('''#define STOICH_MAT %s
    #define NVAR %d
    #define NPAR %d
    #define NREACT %d''' % ( ','.join([str(item) for sublist in S.tolist() for item in sublist]), \
    S.rows, len(Theta), S.cols))
    f.close()
    
    generateModule(model, S)
        

def generateModule(model, S):
    '''Create the wrapper functions for the python module'''
    f = open("C++/pyModuleTemplate.cpp", 'r')
    src = f.readlines()
    f.close()
    
    src2 = [l.replace('myModule', model) for l in src]
    modFile = open("models/%s/%s_LNA.cpp" % (model, model), 'w')
    modFile.writelines(src2)
    modFile.close()
    
    '''Create the python setup script'''
    f = open("python/setupTemplate.py", 'r')
    src = f.readlines()
    f.close()
    
    src2 = [l.replace('myModule', model) for l in src]
    setupFile = open("models/%s/setup.py" % model, 'w')
    setupFile.writelines(src2)
    setupFile.close()
    
    from distutils.core import run_setup
    if not os.path.isdir('modules'):
        os.mkdir('modules')
    run_setup('models/BirthDeath/setup.py', 
        script_args=["build_ext", "install", "--install-lib=modules"])

    
    
def genCcode(f, fName, args ):
    '''generate a very simple C program from the matrix function passed in and the
    arguments necessary'''

    argsList = ', '.join([('const double *' if isinstance(q,tuple) else 'const double ') + p for p,q in args])
    N = len(f)    

    "header file"
    header = """
#ifndef %s_H
#define %s_H
#include "math.h"    
void %s(%s, double varOut[%d]);
#endif""" % (fName.upper(), fName.upper(), fName, argsList, N)
    
    #print header
    "source code"
    code = "#include \"%s.h\"\n" % fName
    code += "void %s(%s, double varOut[%d])\n" % (fName, argsList, N)
    #x = sum(f.reshape(len(f),1).tolist(),[])
    print(f)
    x = sum(f.transpose().reshape(len(f),1).tolist(),[])

    code += "{\n"
    for name, value in args:
        k = 0
        "unpack local variables"
        if isinstance(value,tuple):
            for v in value:            
                code += "double %s = %s[%d];\n" % (v, name, k)
                k += 1
            
    for i in range(N):
        print(x[i])
        #xOut = parseString(str(x[i]))                
        code += "varOut[%d] = %s;\n" % (i, printing.ccode(x[i]))
        
    code += "}"
    #print code
    
    return (header, code)
    
    
