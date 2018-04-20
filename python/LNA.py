# -*- coding: utf-8 -*-
#     LNA.py
#     (C) Justin Feigelman
#     justin.feigelman@helmholtz-muenchen.de
#     Helmholtz Zentrum Muenchen

# This module provides the functions necessary for generating the Python interface to the LNA++ package

from sympy import symbols, Matrix, solve, sqrt, diag, Symbol, printing, zeros
from sympy.matrices import *
from SBML2StoichProp import SBML2StoichProp
import numpy
import os
import ipdb

lnaRootDir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
lnaModelsDir = lnaRootDir + '/models'
lnaModulesDir = lnaModelsDir + '/modules'

def jacobian(x,y):
    "Compute the jacobian of x with respect to y"
    f=x.transpose().reshape(len(x),1)
    return f.jacobian(y)

def triu(X):
    "Return the upper triangular portion of X as a list"
    if not isinstance(X, Matrix):
        X=Matrix(X)
    return sum([Matrix(X)[i,i:].tolist()[0] for i in range(Matrix(X).rows)],[])

def toLinear(X):
	"Return the upper triangular portion of X in column-major ordering"
	if not isinstance(X, Matrix):
		X=Matrix(X)
	if X.shape[0] != X.shape[1]:
		raise(ValueError("Matrix must be square"))
	return sum([sum(X[0:i+1,i].tolist(),[]) for i in range(X.cols)],[])

	
def generateLNA(fileName, model, computeSS='NONE',include_dirs=[], lib_dirs=[]):
    """
    Create a LNA++ Python module for a model provided in an SBML file.
    
    Parameters:
    - sbmlPath: the path to the SBML model specification file
    - model: the desired name of the model 
    - computeSS (optional): one of the values:
      - 'Y0': LNA++ will attempt to calculate the steady state value of the MRE of the stochastic system specified. If `Y0` is not explicitly provided when invoking the model, the computed value will be automatically substituted. 
      - 'V0': LNA++ will attempt to calculate the steady state value of the variance of the stochastic system specified. If `V0` is not explicitly provided when invoking the model, the computed value will be automatically substituted.
      - 'BOTH': LNA++ will compute both `V0` and `Y0` if possible. If `Y0` or `V0` is not explicitly provided when invoking the model, the computed values will be automatically substituted.
      - 'NONE': LNA++ will not compute any steady state values. Invocation of the model without explicitly specifying the initial conditions for `Y0` and `V0` will result in an error message being generated.
    - include_dirs (optional): a Python list of directories to be searched for header files, specified as a list of strings.
    - lib_dirs (optional): a Python list of directories to be searched for libraries, specified as a list of strings.

    Using the created Python module:    
    
    After running generateLNA, a Python module is created and placed in the `modules/` subdirectory of LNA++. The module is named `modelLNA`, where `model' is replaced by the name of the model specified. The module can be imported using `import`, as long as the module is on the Python module path. To add this directory to the module path either modify the value of `sys.path` or add the modules directory to the `PYTHONPATH` environment variable. 
    
    Model simulations are run using the `LNA` function in the module created from the respective model.
    """
    try:
        S, species, parameters, fHandle = SBML2StoichProp(fileName)
    except IOError:
        print("Doing nothing!")
        return
        
    species = symbols(species)
    parameters = symbols(parameters)
        
    tups = generateLNAComponents(model, S, fHandle, species, parameters, computeSS)
    npar = len(parameters) # number of parameters
    compileLNA(model, S, tups, npar, include_dirs=include_dirs, lib_dirs=lib_dirs)


def generateLNAComponents(modelName, S, reactionFlux, phi, Theta, computeSS='NONE'):
    "Generate all of the symbolic components necessary for the LNA method"
    t = symbols('t', real=True)
    Nvar        = S.rows

    phi0 = symbols('phi0(:%d)'%Nvar) # symbolic state variables
    Theta=tuple(Theta)+phi0 # initial conditions are additional parameters for sensitivities

    phi=tuple(phi)
    Npar = len(Theta)

    F = Matrix(reactionFlux(phi, t, Theta))
    J           = F.jacobian(phi)
    dFdTheta    = F.jacobian(Theta)
    d2fdTheta2   = jacobian(dFdTheta,Theta)

    "sensitivities of A"
    A           = S*J
    dAdTheta    = jacobian(A,Theta)
    dAdPhi      = jacobian(A,phi)
    d2AdTheta2  = jacobian(dAdTheta, Theta)
    d2AdPhi2    = jacobian(dAdPhi, phi)

    "sensitivities of E"
    E           = S*diag(*[sqrt(i[0]) for i in F.tolist()])
    dEdTheta    = jacobian(E,Theta)
    d2EdTheta2  = jacobian(dEdTheta, Theta)
    dEdPhi      = jacobian(E,phi)
    d2EdPhi2    = jacobian(dEdPhi, phi)

    global COMPUTE_Y0, COMPUTE_V0
    COMPUTE_Y0 = True
    COMPUTE_V0 = True

    if computeSS is 'Y0' or computeSS is 'NONE':
        COMPUTE_V0 = False

    if computeSS is 'V0' or computeSS is 'NONE':
        COMPUTE_Y0 = False


    "compute initial steady state and sensitivities"
    if COMPUTE_Y0:
        Y0 = solve( S*F, phi)
        S0 = jacobian(Matrix(list(Y0.values())), Theta)
        # substitute values into the dictionary
        Y0 = Matrix([Y0[p] for p in phi])
        S20 = jacobian(S0, Theta)
    else:
        # need placeholders for Y0, S0, and S20
        Y0 = Matrix([0]*Nvar)
        S0 = zeros(Nvar,Npar)
        S20 = zeros(Nvar*Npar,Npar)

    "compute steady state variance"
    if COMPUTE_V0:
        V0 = solveSS_var(A,E,F,S,phi,Theta)
        V0 = V0.subs(list(zip(phi,phi0)))
        "SS variance sensitivities"
        SV0  = jacobian(V0, Theta)
        S2V0 = jacobian(SV0, Theta)
    else:
        # placeholders for V0, SV0, S2V0
        V0 = zeros(Nvar,Nvar)
        SV0 = zeros(Nvar*Nvar,Npar)
        S2V0 = zeros(Nvar*Nvar*Npar,Npar)

    if computeSS not in ['BOTH','V0','Y0','NONE']:
        raise ValueError('Invalid option for computeSS')

    "Solve for System Jacobian"
    N = len(phi)
    V = [[Symbol('V%d%d' % (j,i), real=True) for i in range(N)] for j in range(N)]
    VV = triu(V)

    Phi = [[Symbol('Phi%d%d' % (j,i), real=True) for i in range(N)] for j in range(N)]
    Phi2 = Matrix(Phi).reshape(1,N**2).tolist()[0]
    dVdt = A*Matrix(V) + Matrix(V)*A.T + E*E.T
    sysVar = list(phi)+VV+Phi2
    RHS = sum((S*F).tolist(),[]) + sum([dVdt[i,i:].tolist()[0] for i in range(dVdt.rows)],[]) + (A*Matrix(Phi)).reshape(1,N*N).tolist()[0]
    systemJacobian = jacobian(Matrix(RHS), sysVar)
    systemJacobian_diag = numpy.diag(systemJacobian)

#    "subsitute steady state values"
#    subsPairs =  [(key, Y0[key]) for key in Y0]
#    V0 = V0.subs(subsPairs)
#    SV0 = SV0.subs(subsPairs)
#    S2V0 = S2V0.subs(subsPairs)


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
            'systemJacobian_diag':systemJacobian_diag,
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
            'Y0':Y0}
    
    argTypes1 = ['VECTOR','SCALAR','VECTOR']
    argTypes2 = ['VECTOR']

    tups = [(objs,arg1,argTypes1), (objs2,arg2,argTypes2)]
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
    
    return V0

def compileLNA(model, S, tups, npar, include_dirs=[], lib_dirs=[]):
    '''generate the C code, for python module'''

    # create output directory if necessary
    modelDir = lnaModelsDir+ '/' + model
    if not os.path.isdir(modelDir):
        os.mkdir(modelDir)
    if not os.path.isdir('%s/C/' % modelDir):
        os.mkdir('%s/C/' % modelDir)

    for objs,args,argTypes in tups:
        for x in objs:
            H,C=genCcode(objs[x], x, args, argTypes)
            f=open("%s/C/%s.h" % (modelDir,x), 'w')
            f.writelines(H)
            f.close()

            f=open("%s/C/%s.c" % (modelDir,x), 'w')
            f.writelines(C)
            f.close()

    # generate MODEL_DEF.h
    f = open("%s/C/MODEL_DEF.h" % modelDir, 'w')
    f.writelines(['#define STOICH_MAT %s\n' % ','.join([str(item) for sublist in S.tolist() for item in sublist]),
	'#define NVAR %d\n' % S.rows, '#define NPAR %d\n' % npar, '#define NREACT %d\n' % S.cols])

    global COMPUTE_V0, COMPUTE_Y0
    if COMPUTE_Y0:
        f.writelines('#define COMPUTE_Y0\n')
    if COMPUTE_V0:
        f.writelines('#define COMPUTE_V0')

    f.close()

    generateModule(model, S, modelDir, include_dirs, lib_dirs)


def generateModule(model, S, modelDir = '.', include_dirs=[], lib_dirs=[]):
    '''Create the wrapper functions for the python module'''
    f = open(lnaRootDir + "/src/pyModuleTemplate.cpp", 'r')
    src = f.readlines()
    f.close()

    src2 = [l.replace('myModule', model) for l in src]
    modFile = open("%s/%s_LNA.cpp" % (modelDir, model), 'w')
    modFile.writelines(src2)
    modFile.close()

    '''Create the python setup script'''
    f = open(lnaRootDir + "/python/setupTemplate.py", 'r')
    src = f.readlines()
    f.close()

    src2 = [l.replace('myModule', model) for l in src]
    src2 = [l.replace('LNA_ROOT_DIR', lnaRootDir) for l in src2]

    if not isinstance(include_dirs,list):
    	include_dirs=[include_dirs]
    if not isinstance(lib_dirs,list):
    	lib_dirs=[lib_dirs]

    src2 = [l.replace('INCLUDE_DIRS',str(include_dirs)) for l in src2]
    src2 = [l.replace('LIB_DIRS',str(lib_dirs)) for l in src2]
    setupFile = open("%s/setup.py" % modelDir, 'w')
    setupFile.writelines(src2)
    setupFile.close()

    from distutils.core import run_setup
    if not os.path.isdir(lnaModulesDir):
        os.mkdir(lnaModulesDir)
    run_setup('%s/setup.py' % modelDir,
        script_args=["build_ext", "install", "--install-lib=%s" % lnaModulesDir])
    print('Finished')



def genCcode(f, fName, args, argTypes ):
    '''generate a very simple C program from the matrix function passed in and the
    arguments necessary'''

    argsList = ''
    for k in range(len(args)):
        if argTypes[k] is 'VECTOR':
            argsList +=  'const double * ' + args[k][0]
        else:
            argsList += 'const double ' + args[k][0]
        if k < len(args)-1:
            argsList += ','
            
#    argsList = ', '.join([('const double *' if isinstance(q,tuple) else 'const double ') + p for p,q in args])
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
    #print(f)
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
        #print(x[i])
        #xOut = parseString(str(x[i]))
        code += "varOut[%d] = %s;\n" % (i, printing.ccode(x[i]))

    code += "}"
    #print code

    return (header, code)


