/*
 * LNA_Onset.cpp
 *
 *  Created on: Oct 22, 2013
 *      Author: justinfeigelman
 *
*/

extern "C" {
	#include "mex.h"
}

#include "computeLinearNoise.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "cvodes/cvodes.h"
#include "nvector/nvector_serial.h"
#include "MODEL_DEF.h" // model specific

using namespace std;

static int nvar, npar;

inline bool isSingleton(const mxArray *arr) {
	return ( (mxGetNumberOfDimensions(arr)==2) && (std::min( mxGetN(arr), mxGetM(arr))==1) );
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {

	// check input/output  arguments

	/* input ***
	 * Theta
	 * time
	 * Y0
	 * V0
	 * merr
	 * observed variables
	*/

	static const char* usageMsg =
			"Usage:\n"
			"[MRE, Sigma, dMRE_dTheta, dSigma_dTheta, d2MRE_dTheta2, d2Sigma_dTheta2] = "
                        "LNA_Onset(Theta, timepoints, [Y0], [V0], [merr], [obsVar])\n"
			"Theta: vector of model parameters\n"
                        "timepoints: array of times for output\n"
                        "Y0: initial conditions for each variable\n"
                        "V0: initial variances (upper triangular portion)\n"
                        "merr: measurement variance for observed variables\n"
                        "obsVar: variables to output\n";

	if( (nrhs < 2) || (nrhs>6) || (nlhs>6)) {
		// require at least Theta and the time
		mexErrMsgIdAndTxt("LNA_Onset:nrhs", usageMsg);
	}

	/*
	 * process input arguments
	 */

	enum { ixTheta,
               ixTime,
               ixY0,
               ixV0,
               ixMerr,
               ixObsVar
	};

	// number of model parameters
	// from model def header
	static const int npar = NPAR+NVAR;
	static const int nvar = NVAR;
	static const int Nreact = NREACT;
	static double S_MAT[NREACT*NVAR] = {STOICH_MAT};
	static const MA2 S(S_MAT, shape(NVAR,NREACT), neverDeleteData);

	/* model parameters */
	if (!mxIsDouble(prhs[ixTheta]))
		mexErrMsgIdAndTxt("LNA_Onset:Theta", "Theta must be a numeric vector");
	if (!isSingleton(prhs[ixTheta]))
		mexErrMsgIdAndTxt("LNA_Onset:Theta", "Theta must be a numeric vector");
	double *Theta_ptr 			= (double*)mxGetPr(prhs[ixTheta]);
	double *Theta = new double[npar];
	for (int i=0; i<(npar-nvar); i++) {
		Theta[i] = Theta_ptr[i];
	}

	/* time span */
	Vector tspan;

	if (!mxIsDouble(prhs[ixTime]))
		mexErrMsgIdAndTxt("LNA_Onset:time", "time must be a numeric vector");
	if (!isSingleton(prhs[ixTheta]))
		mexErrMsgIdAndTxt("LNA_Onset:time", "time must be a numeric vector");
	double *tspan_ptr = (double*)mxGetPr(prhs[ixTime]);
	// length of tspan
	int N 		= (int)mxGetNumberOfElements(prhs[ixTime]);
	tspan.resize(N);
	for (int i=0; i<N; i++)
		tspan(i) = tspan_ptr[i];

	/* optional input arguments */
	double *Y0 = NULL, *V0 = NULL, *merr = NULL;
	Array<int,1> varObs;
	int nObsVar;
	LNA::SS_FLAG SS;

	int nargin;
	if (nrhs == 2)
		nargin = 2;
	else if (mxIsEmpty(prhs[2])) {
		nrhs = 2; // discard empty varargin
		nargin = 2;
	} else
		nargin = 2+mxGetNumberOfElements(prhs[2]);

	/* obsvar */
        if (nargin >= ixObsVar + 1 && !mxIsEmpty(mxGetCell(prhs[2], ixObsVar-2))) {

		if (!mxIsCell(prhs[2]))
			mexErrMsgIdAndTxt("LNA_Onset:obsVar", "invalid parameters. varargin must be a cell array");

		mxArray* obsVarArray = mxGetCell(prhs[2], ixObsVar-2);
		if (obsVarArray == NULL)
			mexErrMsgIdAndTxt("LNA_Onset:obsVar", "invalid parameters");

		if (!mxIsDouble(obsVarArray))
			mexErrMsgIdAndTxt("LNA_Onset:obsVar", "obsVar must be a numeric vector");
		if (!isSingleton(obsVarArray))
			mexErrMsgIdAndTxt("LNA_Onset:obsVar", "obsVar must be a numeric vector");

		nObsVar = mxGetNumberOfElements(obsVarArray);
		varObs.resize(nObsVar);

		// fill in observed variables
		double* obsVar_data = (double*)mxGetPr(obsVarArray);
		for (int i=0; i<nObsVar; i++) {
			varObs(i) = (int)obsVar_data[i] - 1;
			if ( (varObs(i) < 0) || (varObs(i) > (nvar-1)))
				mexErrMsgIdAndTxt("LNA_Onset:obsVar", "Invalid observable variable specified.");
		}
	} else {
		// default:  all variables
		varObs.resize(nvar);
		for (int i=0; i<nvar; i++) {
			varObs(i)=i;
		}
		nObsVar = nvar;
	}

	/* merr */
        if (nargin >= ixMerr + 1 && !mxIsEmpty(mxGetCell(prhs[2], ixMerr-2))) {
		const mxArray *merrArray = mxGetCell(prhs[2], ixMerr-2);
		if (!mxIsDouble(merrArray))
			mexErrMsgIdAndTxt("LNA_Onset:merr", "merr must be a numeric vector");
		if (!isSingleton(merrArray))
			mexErrMsgIdAndTxt("LNA_Onset:merr", "merr must be a numeric vector");
		if (mxGetNumberOfElements(merrArray) == 1) {
			double *merrVal = (double*)(mxGetPr(merrArray));
			if (*merrVal < 0.0)
				mexErrMsgIdAndTxt("LNA_Onset:merr", "Invalid measurement error specified.");

			// same error for all species
			merr = new double[nObsVar];
			for (int i=0; i<nObsVar; i++)
				merr[i] = *merrVal;
		} else if (mxGetNumberOfElements(merrArray) == nObsVar) {
			// specified for all species
			merr = new double[nObsVar];
			double *merrVal = (double*)(mxGetPr(merrArray));
			for (int i=0; i<nObsVar; i++) {
				merr[i] = merrVal[i];
				if (merr[i] < 0.0) {
					delete[] merr;
					mexErrMsgIdAndTxt("LNA_Onset:merr", "Invalid measurement error specified.");
				}
			}
		} else {
			mexErrMsgIdAndTxt("LNA_Onset:merr", "specified incorrect merr dimensions ");
		}
	} else {
		merr = new double[nObsVar];
		for (int i=0; i<nObsVar; i++)
			merr[i] = 0.0;
	}

	/* Y0 */
        if (nargin >= ixY0 + 1 && !mxIsEmpty(mxGetCell(prhs[2], ixY0-2))) {
		const mxArray *Y0_Array = mxGetCell(prhs[2], ixY0-2);
		if (!mxIsDouble(Y0_Array))
			mexErrMsgIdAndTxt("LNA_Onset:Y0", "Y0 must be a numeric vector");
		if (!isSingleton(Y0_Array))
			mexErrMsgIdAndTxt("LNA_Onset:Y0", "Y0 must be a numeric vector");
		if (mxGetNumberOfElements(Y0_Array) != nvar)
			mexErrMsgIdAndTxt("LNA_Onset:Y0", "Incorrect dimensions for Y0");
		Y0 = (double*)mxGetPr(Y0_Array);

	} else {
		#ifndef COMPUTE_Y0
		// if not Y0 was not computed when generating the model, then it has to be specified at run time
		mexErrMsgIdAndTxt("LNA_Onset:Y0", "Y0 must be explicitly specified");
		#endif
	}

	/* V0 */
        if (nargin >= ixV0 + 1 && !mxIsEmpty(mxGetCell(prhs[2], ixV0-2))) {
		const mxArray *V0_Array = mxGetCell(prhs[2], ixV0-2);
		if (!mxIsDouble(V0_Array))
			mexErrMsgIdAndTxt("LNA_Onset:V0", "V0 must be a numeric vector");
		if (!isSingleton(V0_Array))
			mexErrMsgIdAndTxt("LNA_Onset:V0", "V0 must be a numeric vector");
		if (mxGetNumberOfElements(V0_Array) != nvar*(nvar+1)/2)
			mexErrMsgIdAndTxt("LNA_Onset:V0", "Incorrect dimensions for V0");
		V0 = (double*)mxGetPr(V0_Array);
	} else {
		#ifndef COMPUTE_V0
		// if not V0 was not computed when generating the model, then it has to be specified at run time
		mexErrMsgIdAndTxt("LNA_Onset:V0", "V0 must be explicitly specified");
		#endif
	}


	/* check init. cond SS */
	if ((Y0 == NULL) && (V0 == NULL))
		SS = LNA::SS_BOTH;
	else if ((Y0 == NULL) && (V0 != NULL))
	   SS = LNA::SS_Y0;
	else if ((Y0 != NULL) && (V0 == NULL))
		SS = LNA::SS_V0;
	else
		SS = LNA::SS_NONE;

	/* compute sensitivities */
	bool computeSens = (nlhs > 2),
		computeSens2 = (nlhs > 4);


//	parameters pars;
//	pars.npar 	= npar;
//	pars.Theta 	= Theta;
//	pars.nvar 	= nvar;
//	pars.S 		= &MA_S;

	MA2 *Y 		= new MA2(nObsVar,N,ColumnMajorArray<2>());
	MA4 *Sigma 	= new MA4(nObsVar, nObsVar, N, N, ColumnMajorArray<4>());
	MA3 *traj_deriv = new MA3(nObsVar, npar, N,ColumnMajorArray<3>());
//	MA5 *dSigma 	= new MA5(nObsVar, nObsVar, N, N, npar,ColumnMajorArray<5>());
	MA5 *dSigma 	= new MA5(nObsVar, nObsVar, npar, N, N, ColumnMajorArray<5>());

	// second order sens.
	MA4 *Sens2_MRE = new MA4(nObsVar, npar, npar, N, ColumnMajorArray<4>());
	MA6 *Sens2_Var = new MA6(nObsVar, nObsVar, npar, npar, N, N, ColumnMajorArray<6>());

	*Y = 0; *Sigma=0; *traj_deriv=0; *dSigma=0; *Sens2_MRE=0, *Sens2_Var=0;

	outputStruct os;
	os.Y 			= Y;
	os.traj_deriv 	= traj_deriv;
	os.Sigma 		= Sigma;
	os.dSigma 		= dSigma;
	os.Sens2_MRE 	= Sens2_MRE;
	os.Sens2_Var  	= Sens2_Var;

	//static LNA lna(nvar, npar, S);
	LNA *lna = new LNA(nvar, npar, S);

	if (nlhs==0)
		return; // no computation necessary

	try {
		lna->computeLinearNoise(Y0, V0, Theta, computeSens, computeSens2,
						SS, os, tspan, varObs, merr);
		delete lna;
		}
	catch (exception &err) {
		mexPrintf("%s\n", err.what());
		mexErrMsgIdAndTxt("LNA_Onset:general", "Caught LNA++ Error");
	}

	/* MRE solution */
  if (nlhs > 0) {
    plhs[0] = mxCreateDoubleMatrix(nObsVar,N,mxREAL);
    if (plhs[0]==0)
    	mexErrMsgTxt("Error allocation output 0");
	void *Yout = mxGetData(plhs[0]);
	if (Y->size()!=nObsVar*N) {
		mexPrintf("Y.size %d, nvar*N %d\n", Y->size(), nvar*N);
		mexErrMsgIdAndTxt("LNA_Onset:output","Y.size()!=nvar*N");
	}
    memcpy( Yout, (void*)Y->data(), sizeof(double) * Y->size());
  }

  /* Temporal cross-covariance matrix */
    if (nlhs > 1) {
        const mwSize dims_Sigma[] = {static_cast<mwSize>(nObsVar),
                                     static_cast<mwSize>(nObsVar),
                                     static_cast<mwSize>(N),
                                     static_cast<mwSize>(N)};
		plhs[1] = mxCreateNumericArray( (mwSize)4, dims_Sigma, mxDOUBLE_CLASS, mxREAL);
		if (plhs[1]==0)
		    	mexErrMsgTxt("Error allocation output 1");
		void *SigmaData  = mxGetData(plhs[1]);
		// copy memory
		if(Sigma->size() != nObsVar*nObsVar*N*N) {
			mexErrMsgIdAndTxt("LNA_Onset:output","Sigma.size() != nObsVar*nObsVar*N*N");
		}
		memcpy(SigmaData, (void*)Sigma->data(), sizeof(double) * Sigma->size());

    }

    // Sensitivity of MRE
    if (nlhs > 2) {
        const mwSize dims_Sens[] = {static_cast<mwSize>(nObsVar),
                                    static_cast<mwSize>(npar),
                                    static_cast<mwSize>(N)};
    	plhs[2] = mxCreateNumericArray( (mwSize)3, dims_Sens, mxDOUBLE_CLASS, mxREAL);
    	if (plhs[2]==0)
    	    	mexErrMsgTxt("Error allocation output 2");
    	void *SensData = mxGetData(plhs[2]);
    	if( traj_deriv->size() != nObsVar*npar*N)
			mexErrMsgIdAndTxt("LNA_Onset:output","traj_deriv.size() != nObsVar*npar*N");

		memcpy(SensData, (void*)traj_deriv->data(), traj_deriv->size() * sizeof(double));
    }

    // Sensitivity of the covariance
    if (nlhs > 3) {
                const mwSize dims[] = {static_cast<mwSize>(nObsVar),
                                       static_cast<mwSize>(nObsVar),
									   static_cast<mwSize>(npar),
									   static_cast<mwSize>(N),
                                       static_cast<mwSize>(N)
                                       };
		plhs[3] = mxCreateNumericArray((mwSize)5, dims, mxDOUBLE_CLASS, mxREAL);
		if (plhs[3]==0)
		    	mexErrMsgTxt("Error allocation output 3");
		void *dSigmaData = mxGetData(plhs[3]);
		if( dSigma->size() != nObsVar*nObsVar*N*N*npar)
			mexErrMsgIdAndTxt("LNA_Onset:output","dSigma.size() != nObsVar*nObsVar*N*N*npar");

		memcpy(dSigmaData,(void*)dSigma->data(), sizeof(double) * dSigma->size());
    }

    // Second order sens. of MRE
    if (nlhs > 4) {
        const mwSize dims_Sens2[] = {static_cast<mwSize>(nObsVar),
                                     static_cast<mwSize>(npar),
                                     static_cast<mwSize>(npar),
                                     static_cast<mwSize>(N)};
    	plhs[4] = mxCreateNumericArray( (mwSize)4, dims_Sens2, mxDOUBLE_CLASS, mxREAL);
    	if (plhs[4]==0)
    	    	mexErrMsgTxt("Error allocation output 4");
    	void *Sens2Data = mxGetData(plhs[4]);
    	if( Sens2_MRE->size() != nObsVar*npar*npar*N)
			mexErrMsgIdAndTxt("LNA_Onset:output","Sens2_MRE.size() != nObsVar*npar*npar*N");
		memcpy(Sens2Data, (void*)Sens2_MRE->data(), Sens2_MRE->size() * sizeof(double));
    }

    // Second order sens. of Var
    if (nlhs > 5) {
        const mwSize dims_Sens2_Var[] = {static_cast<mwSize>(nObsVar),
                                         static_cast<mwSize>(nObsVar),
                                         static_cast<mwSize>(npar),
                                         static_cast<mwSize>(npar),
                                         static_cast<mwSize>(N),
                                         static_cast<mwSize>(N)};
    	plhs[5] = mxCreateNumericArray( (mwSize)6, dims_Sens2_Var, mxDOUBLE_CLASS, mxREAL);
    	if (plhs[5]==0)
    	    	mexErrMsgTxt("Error allocation output 5");
    	void *Sens2_Var_Data = mxGetData(plhs[5]);
    	if( Sens2_Var->size() != nObsVar*nObsVar*npar*npar*N*N)
			mexErrMsgIdAndTxt("LNA_Onset:output","Sens2_Var.size() != nObsVar*nObsVar*npar*npar*N*N");

		memcpy(Sens2_Var_Data, (void*)Sens2_Var->data(), Sens2_Var->size() * sizeof(double));
    }



    // clean up
    delete[] Theta;
    delete[] merr;
    delete Y;
    delete Sigma;
    delete traj_deriv;
    delete dSigma;
   	delete Sens2_MRE;
   	delete Sens2_Var;

}

