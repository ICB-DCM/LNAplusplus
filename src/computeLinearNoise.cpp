#ifndef __cplusplus
#define __cplusplus
#endif

/* Includes */
#include "computeLinearNoise.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <blitz/array.h>
#include <blitz/tv2fastiter.h>
using namespace blitz;


// MATLAB/Python generated code
extern "C" {
#include "V0.h"
#include "Y0.h"
#include "S0.h"
#include "S20.h"
#include "SV0.h"
#include "S2V0.h"
#include "reactionFlux.h"
#include "Afunc.h"
#include "dAdPhi.h"
#include "Efunc.h"
#include "dEdPhi.h"
#include "dFdTheta.h"
#include "d2fdTheta2.h"
#include "d2AdTheta2.h"
#include "d2AdPhi2.h"
#include "d2AdThetadPhi.h"
#include "d2AdPhidTheta.h"
#include "dAdTheta.h"
#include "dEdTheta.h"
#include "d2EdPhi2.h"
#include "d2EdTheta2.h"
#include "d2EdThetadPhi.h"
#include "d2EdPhidTheta.h"
//#include "systemJacobian.h"
//#include "MI.h"
#include "MODEL_DEF.h"
}

//#include "mex.h"

using namespace std;

/* CVODES parameters */
#define RTOL RCONST(1e-6)
#define ATOL RCONST(1e-6)
#define ATOLS RCONST(1e-6)

void **myUserData_ptr;
void *myUserData;

#include "cvodes/cvodes_impl.h"

int LNA::computeLinearNoise(const double* _y0, const double *_v0,
		const double *_Theta, const bool computeSens, const bool computeSens2,
		const SS_FLAG SS_flag, outputStruct &output, const Vector &tspan,
		const Array<int,1> &varObs, const double *merr)
{

	Range all 			= Range::all();

	const int N 		= tspan.size();
	const int NvarObs 	= varObs.size();
	double t0			= tspan(0);
	double *y0, *v0, *Theta;

	// sensitivities of the fundamental matrix
	MA5 dFund(nvar, nvar, N, N, npar);
	dFund = 0;

	// second order sensitivity of the fundamental matrix
	MA6 dFund2(nvar,nvar,N,N,npar,npar);

	// Block covariance matrix
	MA4 SigmaBlocks(nvar,nvar,N,N);
	SigmaBlocks = 0;

	// block covariance matrix (observed only)
	MA4 &Sigma = *(output.Sigma);

	// solution of the MREs
	MA2 Y(nvar, N);
	MA2 &Y_out  = *(output.Y);

	// sensitivites of the block covariance matrix (observed only)
	MA5 &dSigma 		= *(output.dSigma);

	// sensitivities of the MRE
	MA3 &traj_deriv 	= *(output.traj_deriv);

	MA4 &Sens2_MRE_out  = *(output.Sens2_MRE);
	MA6 &Sens2_Var_out 	= *(output.Sens2_Var);

	Theta = new double[npar];
	memcpy(Theta, _Theta, sizeof(double)*npar);

	/* set up initial conditions for the ODE */
	y0 = new double[nvar];
	v0 = new double[nvar*(nvar+1)/2];

	if (SS_flag == SS_Y0)
	{
		Y0(Theta, y0); // compute Y0
		memcpy(v0,_v0,sizeof(double)*nvar*(nvar+1)/2);
		// pad theta with the init conditions from y0
		for (int i=0; i<nvar; i++)
			Theta[npar-nvar+i] = y0[i];
	} else if (SS_flag == SS_V0) {
		memcpy(y0,_y0,sizeof(double)*nvar);
		for (int i=0; i<nvar; i++)
			Theta[npar-nvar+i] = y0[i];
		V0(Theta, v0); // compute V0: requires Y0 first
	} else if (SS_flag == SS_BOTH) {
		Y0(Theta, y0); // compute Y0
		for (int i=0; i<nvar; i++)
			Theta[npar-nvar+i] = y0[i];
		V0(Theta, v0); // compute V0: requires Y0 first
	} else {
		memcpy(y0,_y0,sizeof(double)*nvar);
		memcpy(v0,_v0,sizeof(double)*nvar*(nvar+1)/2);
		// pad theta with the init conditions from y0
		for (int i=0; i<nvar; i++)
			Theta[npar-nvar+i] = y0[i];
	}

	parameters pars = { nvar, npar, &S, Theta, computeSens, computeSens2 };

	// copy the initial conditions into the vector y for CVODES
	for (int i=0; i<nvar; i++) {
		NV_Ith_S(y,i) = y0[i];
	}

	// copy the initial variances into the vector y
	int k=0;
	for (int i=0; i<nvar; i++) {
		for (int j=0; j<=i; j++) {
			NV_Ith_S(y,nvar+k) = v0[k];
			k++;
		}
	}
	int ix = nvar+k;

	// fill in identity matrix for the fundamental matrix nvar*nvar
	k=0;
	for (int i=0; i<nvar; i++)
		for (int j=0; j<nvar; j++) {
//			NV_Ith_S(y,ix+i) = (i==j);
			NV_Ith_S(y,ix+k) = (i==j);
			k++;
		}

	// initial values for the output vector
	for (int i=0; i<nvar; i++)
		Y(i, 0) = y0[i];

	// initialize covariance matrix entries (only upper triangular region)
//	MA2 V(nvar*(nvar+1)/2, N, ColumnMajorArray<2>());
	MA2 V(nvar*(nvar+1)/2, N);
	V = 0;
	k=0;
	for (int i=0; i<nvar; i++) {
		for (int j=0; j<=i; j++) {
			V(k,0) = v0[k];
			k++;
		}
	}

	// the solution to the fundamental matrix ODE over all time points
	MA4 Fund(nvar,nvar,N,N);
	Fund=0;

	MA2 Identity(nvar,nvar);
	Identity = 0;
	for (int i=0; i<nvar; i++)
		Identity(i,i) = 1;

	// Fundamental matrix solution at t=0
	Fund(all,all,0,0) = Identity(all,all);

	// 1st order sensitivity of the macroscopic rate equation
	MA3 Sens_MRE(nvar,npar,N);
	Sens_MRE=0;
	// 1st order sensitivity of the variance
	MA4 Sens_Var(nvar,nvar,N,npar);
	Sens_Var=0;
	// 2nd order sensitivity of the MRE
	MA4 Sens2_MRE(nvar,npar,npar,N);
	Sens2_MRE=0;
	// 2nd order sensitivity of the variance
	MA5 Sens2_Var(nvar,nvar,N,npar,npar);
	Sens2_Var=0;

	/* setup the ODE system with the initial conditions, and activate the
	 * sensitivity computation if specified.
	 */
	setupCVODES(t0, pars);
	/* load the sensitivities into the data structure for CVODES */
	if (computeSens)
		setupSens(SS_flag, Sens_MRE, Sens_Var, pars);
	if (computeSens2)
		setupSens2(SS_flag, Sens2_MRE, Sens2_Var, pars);

	realtype tret;

	// outer loop for solution of ODE and sensitivities
	for (int i=0; i<N; i++)
	{
		// initial time
		realtype t1 = tspan(i);

		/* construct the initial conditions I for this integration step */

		// reinitialize the macroscopic rate equations and variance ODE from the cached values at the previous time point
		for (int species=0; species<nvar; species++)
			NV_Ith_S(I,species) = Y(species,i);
		ix=nvar;

		// reinitialize the variances from previous solution
		for (k=0; k<nvar*(nvar+1)/2; k++)
			NV_Ith_S(I,k+nvar) = V(k,i);
		ix=nvar+k;

		// reinitialize the Fundamental matrix from previous solution
		k=0;
		for (int u=0; u<nvar; u++)
			for (int v=0; v<nvar; v++)
			{
				NV_Ith_S(I,ix+k) = ((u==v) ? 1:0); // Identity matrix
				k++;
			}

		// diagonal is always the identity matrix
		Fund(all, all, i,i) = Identity(all,all);

		if (computeSens) {
			// copy cached sensitivity equation solution at time point i into the
			// initial condition Nvector array
			for (int lPar=0; lPar<npar; lPar++)
			{
				// MRE
				for (int lVar=0; lVar<nvar; lVar++)
					NV_Ith_S(yS0[lPar], lVar) = Sens_MRE(lVar,lPar,i);

				// Variance
				k=0; ix=nvar;
				for (int lVar1=0; lVar1<nvar; lVar1++) {
					for (int lVar2=0; lVar2<=lVar1; lVar2++) {
						NV_Ith_S(yS0[lPar], ix+k) = Sens_Var(lVar2,lVar1,i,lPar);
						k++;
					}
				}

				// Fundamental matrix.  sensitivities are zero at time 0.
				ix+=k; k=0;
				for (int lVar1=0; lVar1<nvar; lVar1++)
					for (int lVar2=0; lVar2<nvar; lVar2++) {
						NV_Ith_S(yS0[lPar], ix+k) = 0; //dFund(lVar1,lVar2,i,i,lPar);
						k++;
					}
			}
		}

		if (computeSens2) {
			// reinit from cached sensitivity equation solution at time point i
			for (int lPar1=0; lPar1<npar; lPar1++)
				for (int lPar2=0; lPar2<npar; lPar2++)
				{
					// MRE
					for (int lVar=0; lVar<nvar; lVar++)
						NV_Ith_S(yS0[lPar1*npar + lPar2 + npar], lVar) 	= Sens2_MRE(lVar,lPar1,lPar2,i);


					// Variance
					k=0; ix=nvar;
					for (int lVar1=0; lVar1<nvar; lVar1++) {
						for (int lVar2=0; lVar2<=lVar1; lVar2++) {
							NV_Ith_S(yS0[lPar1*npar + lPar2 + npar], ix+k) = Sens2_Var(lVar2,lVar1,i,lPar1,lPar2);
							k++;
						}
					}

					// Fundamental matrix
					ix+=k; k=0;
					for (int lVar1=0; lVar1<nvar; lVar1++)
						for (int lVar2=0; lVar2<nvar; lVar2++) {
							NV_Ith_S(yS0[lPar1*npar + lPar2 + npar], ix+k) = 0;
							k++;
						}
				}
		}
		// reinitialize with the new initial conditions vector
		CVodeReInit(cvode_mem, t1, I);
//		assert(((CVodeMem)cvode_mem)->cv_user_data == myUserData_ptr);

		if (computeSens || computeSens2) {
			// reinitialize with the new sensitivity vector
			flag = CVodeSensReInit(cvode_mem, CV_SIMULTANEOUS, yS0);
			if (check_flag(&flag, "CVodeSensReInit", 1))
				throw runtime_error("CVodeSensReInit");
		}

		for (int j=i+1; j<N; j++) {
			//outer loop for propagation of fluctuations

			// integrate to next time point
			realtype t2 = tspan(j);
			flag = CVode(cvode_mem, t2, x, &tret, CV_NORMAL);
			if (check_flag(&flag, "CVode", 1))
				throw runtime_error("CVode");

			if (computeSens) {
				// extract sensitivity solution
				flag = CVodeGetSens(cvode_mem, &tret, ySout);
				if (check_flag(&flag, "CVodeGetSens", 1))
					throw runtime_error("CVodeGetSens");
			}

			// results of integration
			double phi[nvar];
			static MA2 V_i_j(nvar,nvar), // covariance from t_i to t_j
					Fund_i_j(nvar,nvar); // Fundamental matrix from t_i to t_j

			// get RHS of ODE
			unpackYDot(x, phi, V_i_j, Fund_i_j);

			// results of sensitivity equation
			static MA2 Sens_MRE_j(nvar,npar);
			static MA3 Xi_i_j(nvar,nvar,npar), Sens_Var_j(nvar,nvar,npar);
			if (computeSens)
				unpackYSDot(ySout, Sens_MRE_j, Sens_Var_j, Xi_i_j ); // RHS of F.S.

			static MA3 Sens2_MRE_j(nvar,npar,npar);
			static MA4 Sens2_Var_j(nvar,nvar,npar,npar), Xi2_i_j(nvar,nvar,npar,npar);
			if (computeSens2) // RHS of higher order sensitivities
				unpackYS_HO(ySout, Sens2_MRE_j, Sens2_Var_j, Xi2_i_j );

			// store the solution of the integration (only from t_i to t_j)
			if (j==(i+1)) {
				// MRE
				for (int ix1=0; ix1<nvar; ix1++)
					Y(ix1,j) = phi[ix1];

				// Variance
				k=0;
				for (int lVar1=0; lVar1<nvar; lVar1++)
					for (int lVar2=0; lVar2<=lVar1; lVar2++) {
						V(k,j) = V_i_j(lVar2,lVar1);
						k++;
					}

				if (computeSens) {
					// Sensitivity of the MRE
					Sens_MRE(all, all, j) 		= Sens_MRE_j;
					// Sensitivity of the Variance
					Sens_Var(all, all, j, all) 	= Sens_Var_j;
				}
				if (computeSens2) {
					// second order sensitivities
					// MRE
					Sens2_MRE(all,all,all,j) 	= Sens2_MRE_j;
					// variance
					Sens2_Var(all,all,j,all,all) = Sens2_Var_j;
				}
			}

			// Sensitivity of fundamental matrix
			if (computeSens) {
				dFund(all, all, i, j, all) = Xi_i_j;
			}
			// second order sens. of fund. matr.
			if (computeSens2) {
				dFund2(all,all,i,j,all,all) = Xi2_i_j;
			}

			// Fundamental matrix solution
			Fund(all, all, i, j) = Fund_i_j;
		}
	}

	/* end of integration loop! */

	/* construct outputs */

	// construct tensor of variance solution at all time points
	MA3 Sigma_t(nvar,nvar,N);
	Sigma_t=0;

	k=0;
	for (int lVar1=0; lVar1<nvar; lVar1++)
		for (int lVar2=0; lVar2<=lVar1; lVar2++) {
			Sigma_t(lVar2, lVar1, all) = V(k,all);
			if (lVar2!=lVar1)
				Sigma_t(lVar1, lVar2, all) = V(k,all); //symmetric
			k++;
		}

#ifdef DEBUG
	cout << "Variance trajectory " << endl;
	cout << Sigma_t;
#endif

	/* construct the (upper triangular) block covariance matrix between all species at each time point */
	static MA2 SigmaBlocks_ii(nvar,nvar), Fund_i_ii(nvar,nvar),
			SigmaBlocks_i_ii(nvar,nvar);

	for (int i=0; i<N; i++) {
		// fill in the diagonal elements
		SigmaBlocks(all,all,i,i) 	= Sigma_t(all,all,i);
		SigmaBlocks_ii 				= Sigma_t(all,all,i);

		// off diagonal elements
		for (int ii=i+1; ii<N; ii++) {
			// get the fundamental matrix from time i to time ii
			Fund_i_ii = Fund(all,all,i,ii);

			firstIndex a; secondIndex b; thirdIndex c;
			SigmaBlocks_i_ii = sum(SigmaBlocks_ii(a,c)*Fund_i_ii(b,c),c); // transpose!
			SigmaBlocks(all,all,i,ii) = SigmaBlocks_i_ii;
		}
	}

	// tensor views
	static MA2 dSigmaBlocks_ii_r(nvar,nvar), dSigmaBlocks_i_ii_r(nvar,nvar),
		dFund_i_ii_r(nvar,nvar);


	MA5 dSigmaBlocks(nvar,nvar,N,N,npar);

	if (computeSens) {
		// sensitivity of the block covariance matrix
		dSigmaBlocks = 0;

		for (int lPar=0; lPar<npar; lPar++) {
			for (int i=0; i<N; i++) {
				// fill in diagonals
				dSigmaBlocks(all,all,i,i,lPar) 	= Sens_Var(all,all,i,lPar);
				dSigmaBlocks_ii_r 			 	= Sens_Var(all,all,i,lPar);

				SigmaBlocks_ii = SigmaBlocks(all,all,i,i);

				// off-diagonals
				for (int ii=i+1; ii<N; ii++) {
					// Fundamental matrix from i to ii
					Fund_i_ii = Fund(all,all,i,ii);

					// Derivative of Fundamental matrix from i to ii
					dFund_i_ii_r = dFund(all,all,i,ii,lPar);

					firstIndex a; secondIndex b; thirdIndex c;

					dSigmaBlocks_i_ii_r 	= sum(Fund_i_ii(a,c)*dSigmaBlocks_ii_r(c,b),c);
					dSigmaBlocks_i_ii_r  	+= sum( dFund_i_ii_r(a,c)*SigmaBlocks_ii(c,b),c);
					dSigmaBlocks_i_ii_r 	= dSigmaBlocks_i_ii_r.transpose(secondDim,firstDim);

					// copy back into the block covariance sensitivities matrix (upper triangular)
					dSigmaBlocks(all,all,i,ii,lPar) = dSigmaBlocks_i_ii_r;
				}
			}
		}
	}

	MA6 d2SigmaBlocks(nvar,nvar,N,N,npar,npar);
	static MA3 dSigmaBlocks_ii(nvar,nvar,npar);

	if (computeSens2) {
		d2SigmaBlocks = 0;

		for (int i=0; i<N; i++) {
			dSigmaBlocks_ii	= Sens_Var(all,all,i,all); // dVdTheta

			// tensor views
			static MA3 dFund_i_ii(nvar,nvar,npar), dSigmaBlocks_ii(nvar,nvar,npar);
			static MA4 d2Fund_i_ii(nvar,nvar,npar,npar), d2SigmaBlocks_ii(nvar,nvar,npar,npar),
					d2SigmaBlocks_i_ii(nvar,nvar,npar,npar);
			static MA5 vvvpp(nvar,nvar,nvar,npar,npar);

			// fill in diagonals
			d2SigmaBlocks(all,all,i,i,all,all) 	= Sens2_Var(all,all,i,all,all);
			d2SigmaBlocks_ii 			 		= Sens2_Var(all,all,i,all,all);

			SigmaBlocks_ii = SigmaBlocks(all,all,i,i);

			// off-diagonals

			for (int ii=i+1; ii<N; ii++) {
				// Fundamental matrix from i to ii
				Fund_i_ii = Fund(all,all,i,ii);

				// Derivative of Fundamental matrix from i to ii
				dFund_i_ii 	= dFund(all,all,i,ii,all);
				d2Fund_i_ii = dFund2(all,all,i,ii,all,all);

				firstIndex a; secondIndex b; thirdIndex c; fourthIndex d; fifthIndex e;
				vvvpp = d2SigmaBlocks_ii(a,c,d,e)*Fund_i_ii(b,c); // d2VdTheta_i_j*Phi'
				d2SigmaBlocks_i_ii 	= sum(vvvpp(a,b,e,c,d),e);
				vvvpp = dSigmaBlocks_ii(a,c,d)*dFund_i_ii(b,c,e); // dVdTheta_i*dPhi'dTheta_j
				d2SigmaBlocks_i_ii += sum(vvvpp(a,b,e,c,d),e);
//				vvvpp = dSigmaBlocks_ii(a,c,e)*dFund_i_ii(c,b,d); // dVdTheta_j*dPhidTheta_i
				vvvpp = dSigmaBlocks_ii(a,c,e)*dFund_i_ii(b,c,d); // dVdTheta_j*dPhi'dTheta_i
				d2SigmaBlocks_i_ii += sum(vvvpp(a,b,e,c,d),e);
//				vvvpp = SigmaBlocks_ii(a,c)*d2Fund_i_ii(c,b,d,e); // V*d2PhidTheta_i_j
				vvvpp = SigmaBlocks_ii(a,c)*d2Fund_i_ii(b,c,d,e); // V*d2Phi'dTheta_i_j
				d2SigmaBlocks_i_ii += sum(vvvpp(a,b,e,c,d),e);
				d2SigmaBlocks(all,all,i,ii,all,all) = d2SigmaBlocks_i_ii;
			}
		}
#ifdef DEBUG
		cout << "Sigma Blocks " << endl << SigmaBlocks << endl;
#endif
	}

	/* extract the observable variables */

	for (int lVar = 0; lVar<NvarObs; lVar++)
	{
		int var=varObs(lVar);
		Y_out(lVar,all) = Y(var,all);
	}

	for (int lVar1 = 0; lVar1<NvarObs; lVar1++) {
		int var1=varObs(lVar1);
		for (int lVar2 = 0; lVar2<NvarObs; lVar2++) {
			int var2=varObs(lVar2);

			for (int i=0; i<N; i++) {
				for (int ii=i; ii<N; ii++) {
					Sigma( lVar1, lVar2, i, ii) = SigmaBlocks(var1, var2, i, ii);

					if ((lVar1 == lVar2) && (i==ii))
						Sigma(lVar1,lVar2,i,ii) += merr[lVar1]; // add measurement error to variances
					if (ii>i)
						Sigma( lVar1, lVar2,ii,i ) = Sigma( lVar1, lVar2, i, ii ); // Symmetry in block cov. matrix
				}
			}
		}
	}

	// sensitivity of the block covariance matrix
	if (computeSens) {
		// loop over rows of small covariance
		for (int lVar1 = 0; lVar1<NvarObs; lVar1++)
			// loop over cols of small covariance
			for (int lVar2 = 0; lVar2<NvarObs; lVar2++) {
				int var1=varObs(lVar1), var2=varObs(lVar2);
				// loop over rows of big covariance
				for (int i=0; i<N; i++)
					// loop over cols of big covariance
					for (int ii=i; ii<N; ii++) {
						dSigma(lVar1, lVar2, i, ii, all) 	= dSigmaBlocks(var1, var2, i, ii, all);
						if (ii>i)
							dSigma(lVar1, lVar2, ii, i, all) = dSigmaBlocks(var1, var2, i, ii, all);  // Symmetry in block cov. matrix
					}
			}

		// derivatives of MRE for observables
		for (int lVar=0; lVar<NvarObs; lVar++)
		{
			const int var=varObs(lVar);
			traj_deriv(lVar, all, all) = Sens_MRE(var, all, all);
		}

#ifdef DEBUG
		cout << "sensitivity " << endl;
		cout << Sens_MRE;

		cout << "Fundamental matrix solution " << endl;
		cout << Fund;
#endif
	}

	if (computeSens2) {
		// MRE second order sens.
		for (int lVar=0; lVar<NvarObs; lVar++)
		{
			int var=varObs(lVar);
			Sens2_MRE_out(lVar, all, all, all) = Sens2_MRE(var, all, all, all);
		}

		// Variance second order sens.
		for (int lVar1 = 0; lVar1<NvarObs; lVar1++)
			// loop over cols of small covariance
			for (int lVar2 = 0; lVar2<NvarObs; lVar2++) {
				int var1=varObs(lVar1), var2=varObs(lVar2);
				// loop over rows of big covariance
				for (int i=0; i<N; i++)
					// loop over cols of big covariance
					for (int ii=i; ii<N; ii++) {
						Sens2_Var_out(lVar1, lVar2, all, all, i, ii) 		= d2SigmaBlocks(var1, var2, i, ii, all, all);
						if (ii>i)
							Sens2_Var_out(lVar1, lVar2, all, all, ii, i) 	= d2SigmaBlocks(var1, var2, i, ii, all, all);
					}
			}
	}

	delete[] y0;
	delete[] v0;
	delete[] Theta;

	return(0); // succesful return
}


int LNA::check_flag(void *flagvalue, const char *funcname, int opt)
{
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == 0) {
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1); }

	/* Check if flag < 0 */
	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			cerr << "\nSUNDIALS_ERROR: " << funcname << "failed with flag = " << *errflag << "\n\n";
			return(1); }}

	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == 0) {
		fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1); }

	return(0);
}


int LNA::fundRHS(realtype t, N_Vector yIn, N_Vector ydot, void *user_data) {
	parameters *par = (parameters*)user_data;
	static const int nvar 		= par->nvar;
	const double 	*Theta 	= par->Theta;
	// stoichiometric matrix
	static MA2   S	= par->S->copy();

	//  MRE solution
	double phi[nvar];
	// Fundamental Matrix
	static MA2 Phi(nvar,nvar);
	// current Variance
	static MA2 V(nvar,nvar);
	// get the current state out of y
	unpackYDot(yIn, phi, V, Phi);

	static const TinyVector<int, 2> Sdim = S.shape();
	static const int Nreact 	= Sdim(1);

	// get the reaction fluxes
	static double* flux = new double[Nreact];
	reactionFlux(phi, t,  Theta, flux);
	static Vector FluxVec(flux, shape(Nreact), deleteDataWhenDone);

	firstIndex a; secondIndex b;

	static Vector netChangeVec(nvar);
	// dot product!
	netChangeVec = sum(S(a,b)*FluxVec(b), b);

	// rhs for variances
	static double *A_mem = new double[nvar*nvar];
	Afunc(phi, t, Theta, A_mem);
	static MA2 A(A_mem, shape(nvar,nvar), neverDeleteData, ColumnMajorArray<2>());

#ifdef DEBUG
	cout << "A " << endl << A << endl;
#endif
	static double* E_mem = new double[nvar*Nreact];
	Efunc(phi, t, Theta, E_mem);
	static MA2 E(E_mem, shape(nvar,Nreact), deleteDataWhenDone, ColumnMajorArray<2>());

#ifdef DEBUG
	cout << "E " << endl << E << endl;
#endif

	static MA2 dVdt(nvar,nvar);
	thirdIndex c;

	dVdt  = sum(A(a,c)*V(c,b),c);
	dVdt += sum(V(a,c)*A(b,c),c); // A transpose
	dVdt += sum(E(a,c)*E(b,c),c); // E transpose

#ifdef DEBUG
	cout << "dVdt " << endl << dVdt << endl;
#endif

	// derivatives of the fundamental matrix
	static MA2 dPhidt(nvar,nvar);
	dPhidt = sum(A(a,c)*Phi(c,b),c);


#ifdef DEBUG
	cout << "Phi " << endl << Phi << endl;
	cout << "dPhidt = A*Phi" << endl;
	cout << dPhidt;
#endif

	packYDot(ydot, netChangeVec, dVdt, dPhidt);

#ifdef DEBUG
	cout << "Ydot: " ;
	for (int i=0; i<9; i++)
		cout << setw(15) << NV_Ith_S(ydot,i);
	cout << endl;

	cout << "Y (" << t << ") "  ;
	for (int i=0; i<9; i++)
		cout << setw(15) << NV_Ith_S(yIn,i);
	cout << endl;

	cout << "Ydot (" << t << ") "  ;
	for (int i=0; i<9; i++)
		cout << setw(15) << NV_Ith_S(ydot,i);
	cout << endl;

#endif

	return 0;

}

// system Jacobian df/dy
int LNA::Jac(long int N, realtype t,
		N_Vector y, N_Vector fy, DlsMat J, void *user_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

	/*
	parameters *par = (parameters*)user_data;
	int nvar 		= par->nvar;
	const double 	*Theta 	= par->Theta;


	double phi[nvar];
	static MA2 V(nvar,nvar), Phi(nvar,nvar);

	LNA::unpackYDot(y, phi, V, Phi);

	const int RHS_SIZE = (nvar*(nvar+3)/2 + nvar*nvar);

	double *jac_mem = new double[RHS_SIZE*RHS_SIZE];
	systemJacobian(phi,t,Theta,jac_mem);

	static MA2 Jacobian(jac_mem, shape(RHS_SIZE, RHS_SIZE), neverDeleteData);

	//	MA2 Jacobian(RHS_SIZE,RHS_SIZE,jac_mem);

	for (int i=0; i<RHS_SIZE; i++)
		for (int j=0; j<RHS_SIZE; j++)
			DENSE_ELEM(J, i, j) = Jacobian(i,j);

	delete[] jac_mem;
	*/
	return 0;
}


// preconditioner matrix
int LNA::Preconditioner(realtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
		realtype gamma, realtype delta, int lr, void *user_data, N_Vector tmp) {

	// compute the preconditioner matix P as M=I-gamma*J, and then solve for z by computing
	// inv(M)*r
	/*
	parameters *par = (parameters*)user_data;
	int nvar 		= par->nvar; //, npar=par->npar;
	const double 	*Theta 	= par->Theta;
	const int RHS_SIZE = (nvar*(nvar+3)/2 + nvar*nvar);

	static double *MI_mem = new double[RHS_SIZE*RHS_SIZE];

	double phi[nvar];
	static MA2 V(nvar,nvar), Phi(nvar,nvar);

	LNA::unpackYDot(y, phi, V, Phi);
	MI(phi,t,Theta,gamma, MI_mem);

	static MA2 myMI(MI_mem, shape(RHS_SIZE,RHS_SIZE), deleteDataWhenDone, ColumnMajorArray<2>());

	Vector myR(N_VGetArrayPointer_Serial(r), shape(RHS_SIZE), duplicateData);

	firstIndex a;
	secondIndex b;

	static Vector myZ(RHS_SIZE);
	myZ = sum( myMI(a,b)*myR(b), b);

	for (int i=0; i<RHS_SIZE; i++)
		NV_Ith_S(z,i) = myZ(i);
	*/

	return 0;


}

/* rhs to the sensitivity equations */

int LNA::sensRhs(int Ns, realtype t, N_Vector y, N_Vector ydot,
		N_Vector *yS, N_Vector *yS_dot,
		void *user_data, N_Vector tmp1, N_Vector tmp2) {

	parameters *par = (parameters*)user_data;
	static const int nvar 		= par->nvar, npar=par->npar;
	const double 	*Theta 	= par->Theta;
	// stoichiometric matrix
	static MA2   S	= par->S->copy();
	const bool computeSens2 = par->computeSens2;

	const TinyVector<int,2> Sdim = S.shape();

	static const int Nreact 	= Sdim(1); // number of reactions

	/* extract current state */
	//  MRE solution
	double phi[nvar];
	// Fundamental Matrix
	static MA2 Phi(nvar,nvar);
	// current Variance
	static MA2 V(nvar,nvar);
	// get the current state out of y
	unpackYDot(y, phi, V, Phi);

	static MA2 Sens_MRE(nvar,npar);
	static MA3 Sens_Var(nvar,nvar,npar), Xi(nvar,nvar,npar);
	Sens_MRE=0;
	Sens_Var=0;
	Xi=0;

	unpackYSDot(yS, Sens_MRE, Sens_Var, Xi);

	#ifdef DEBUG
	cout << "Sensitivity t=" << t << endl;
	cout << Sens_MRE;
	cout << "phi " << phi[0] << " " << phi[1] <<endl;

	cout << "Sens_Var: " << endl;
	cout << Sens_Var;

#endif
	static double *A_mem = new double[nvar*nvar];

	Afunc(phi, t, Theta, A_mem);
	static MA2 A(A_mem, shape(nvar, nvar), deleteDataWhenDone, ColumnMajorArray<2>());
//	static MA2 A(A_mem, shape(nvar, nvar), neverDeleteData, ColumnMajorArray<2>());

	// explicit derivative of reaction flux f on theta
	static double *dFdTheta_mem = new double[Nreact*npar];
	dFdTheta(phi, t, Theta, dFdTheta_mem);
	static MA2 mydfdtheta(dFdTheta_mem, shape(Nreact,npar), deleteDataWhenDone, ColumnMajorArray<2>());

	//	RHS of sensitivities for MRE
	static firstIndex i; secondIndex j; thirdIndex k;
	static MA2 dSdt(nvar,npar);
	dSdt  = sum( A(i,k)*Sens_MRE(k,j),k);
	dSdt += sum(S(i,k)*mydfdtheta(k,j),k);

#ifdef DEBUG
	cout << "dfdTheta" << endl;
	cout << mydfdtheta;

	cout << "dSdt at t=" << t<<endl;
	cout << dSdt;
#endif

	// RHS of variance sensitivities

	// dAdTheta
	static double *dAdTheta_mem = new double[nvar*nvar*npar];
	dAdTheta(phi, t, Theta, dAdTheta_mem);
	static MA3 mydAdTheta( dAdTheta_mem, shape(nvar,nvar,npar), deleteDataWhenDone, ColumnMajorArray<3>());

	// dEdTheta
	static double *dEdTheta_mem = new double[nvar*Nreact*npar];
	dEdTheta(phi,t,Theta,dEdTheta_mem);
	static MA3 mydEdTheta(dEdTheta_mem, shape(nvar, Nreact, npar), deleteDataWhenDone, ColumnMajorArray<3>());

	// E
	static double *E_mem = new double[nvar*Nreact];
	Efunc(phi,t,Theta,E_mem);
	static MA2 E(E_mem, shape(nvar,Nreact), deleteDataWhenDone, ColumnMajorArray<2>());

#ifdef DEBUG
	cout << "phi " << phi[0] << " " << phi[1] <<endl;
	cout << "t " << t << endl;
	cout << "Theta" <<endl;
	for (int l=0; l<npar; l++)
		cout << setw(10) << Theta[l] << endl;
	cout << "E" << endl;
	cout << E;
#endif

	static MA3 Sens_Var_dot(nvar,nvar,npar);

	static double *dAdPhi_mem = new double[nvar*nvar*nvar];
	dAdPhi(phi,t,Theta,dAdPhi_mem);
	static MA3 mydAdPhi(dAdPhi_mem, shape(nvar,nvar,nvar), deleteDataWhenDone, ColumnMajorArray<3>());

	static double *dEdPhi_mem = new double[nvar*Nreact*nvar];
	dEdPhi(phi,t,Theta,dEdPhi_mem);
	static MA3 mydEdPhi(dEdPhi_mem, shape(nvar,Nreact,nvar), deleteDataWhenDone, ColumnMajorArray<3>());

	/* compute the time derivative of the sensitivity of the covariance matrix */

	fourthIndex l;

	// Total derivative of A WRT Theta (3rd order Tensor)
	static MA3 dAdTheta_tot(nvar,nvar,npar);
	dAdTheta_tot = sum( mydAdPhi(i,j,l)*Sens_MRE(l,k),l) + mydAdTheta(i,j,k);

	// Total derivative of E WRT Theta
	static MA3 dEdTheta_tot(nvar,Nreact,npar);
	dEdTheta_tot = sum( mydEdPhi(i,j,l)*Sens_MRE(l,k),l) + mydEdTheta(i,j,k);

	// Tensor equation for time derivative of the sensitivities of the variance
	Sens_Var_dot  = sum(dAdTheta_tot(i,l,k)*V(l,j), l);
	Sens_Var_dot += sum(A(i,l)*Sens_Var(l,j,k),l);
	Sens_Var_dot += sum(Sens_Var(i,l,k)*A(j,l),l);
	Sens_Var_dot += sum(V(i,l)*dAdTheta_tot(j,l,k),l);
	Sens_Var_dot += sum(dEdTheta_tot(i,l,k)*E(j,l),l);
	Sens_Var_dot += sum(E(i,l)*dEdTheta_tot(j,l,k),l);
	//

	/* Sensitivities of the fundamental matrix */
	static MA3 Xi_dot(nvar,nvar,npar);

	static MA4 vvvp(nvar,nvar,nvar,npar);
	vvvp = dAdTheta_tot(i,k,l)*Phi(k,j);
	Xi_dot = sum(vvvp(i,j,l,k),l);
	vvvp = A(i,k)*Xi(k,j,l);
	Xi_dot += sum(vvvp(i,j,l,k),l);

	packYSDot(yS_dot, dSdt, Sens_Var_dot, Xi_dot);

	/* Higher order Sensitivities! */

	if (computeSens2) {
		static MA3 Sens2_MRE(nvar,npar,npar);
		static MA3 Sens2_MRE_dot(nvar,npar,npar);
		static MA4 Sens2_Var(nvar,nvar,npar,npar);
		static MA4 Sens2_Var_dot(nvar,nvar,npar,npar);
		static MA4 Xi2(nvar,nvar,npar,npar);
		static MA4 Xi2_dot(nvar,nvar,npar,npar);

		Sens2_MRE_dot 	= 0;
		Sens2_Var_dot 	= 0;
		Xi2_dot 		= 0;

		static double *d2fdTheta2_mem = new double[Nreact*npar*npar];
		d2fdTheta2(phi,t,Theta, d2fdTheta2_mem);
		static MA3 myd2fdTheta2(d2fdTheta2_mem, shape(Nreact,npar,npar), deleteDataWhenDone, ColumnMajorArray<3>());
		static MA3 myd2FdTheta2(nvar,npar,npar);

//		myd2FdTheta2=0;
//		static MA4 tmp1(nvar,npar,npar,Nreact);
//		tmp1 = ;
		firstIndex i; secondIndex j; thirdIndex k; fourthIndex l;
		myd2FdTheta2 = sum(S(i,l)*myd2fdTheta2(l,j,k),l);
//		cout << "t= " << t << endl;
//		cout << "d2fdTheta2 " << endl << myd2fdTheta2 << endl;
//		cout << "d2FdTheta2 " << endl << myd2FdTheta2 << endl;

		// TO DO:  need d2FdTheta2 = S*d2fdTheta2
		// get the current second order sensitivities
		unpackYS_HO(yS, Sens2_MRE, Sens2_Var, Xi2);

		// compute right hand side
		// MRE
		static MA4 vvpp(nvar,nvar,npar,npar);
		vvpp = dAdTheta_tot(i,j,l)*Sens_MRE(j,k); // dAdTheta_j*S_i

		Sens2_MRE_dot	=  sum(vvpp(i,l,j,k),l);
		vvpp = A(i,j)*Sens2_MRE(j,k,l);
		Sens2_MRE_dot 	+= sum(vvpp(i,l,j,k),l);
//
		vvpp = mydAdTheta(i,j,k)*Sens_MRE(j,l); // partial dAdTheta_i*S_j
//
		Sens2_MRE_dot 	+= sum(vvpp(i,l,j,k),l);
		Sens2_MRE_dot 	+= myd2FdTheta2(i,j,k);

		// second derivatives of A and E
		static double *d2AdTheta2_mem 		= new double[nvar*nvar*npar*npar];
		static double *d2AdPhi2_mem 		= new double[nvar*nvar*nvar*nvar];
		static double *d2AdThetadPhi_mem 	= new double[nvar*nvar*npar*nvar];
		static double *d2AdPhidTheta_mem 	= new double[nvar*nvar*nvar*npar];
		static double *d2EdTheta2_mem 		= new double[nvar*Nreact*npar*npar];
		static double *d2EdPhi2_mem 		= new double[nvar*Nreact*nvar*nvar];
		static double *d2EdThetadPhi_mem 	= new double[nvar*Nreact*npar*nvar];
		static double *d2EdPhidTheta_mem	= new double[nvar*Nreact*nvar*npar];

		// evaluate functions
		d2AdTheta2(phi,t,Theta,d2AdTheta2_mem);
		d2AdPhi2(phi,t,Theta,d2AdPhi2_mem);
		d2AdThetadPhi(phi,t,Theta,d2AdThetadPhi_mem);
		d2AdPhidTheta(phi,t,Theta,d2AdPhidTheta_mem);
		d2EdTheta2(phi,t,Theta,d2EdTheta2_mem);
		d2EdPhi2(phi,t,Theta,d2EdPhi2_mem);
		d2EdThetadPhi(phi,t,Theta,d2EdThetadPhi_mem);
		d2EdPhidTheta(phi,t,Theta,d2EdPhidTheta_mem);

		// Array objects
		static MA4 myd2AdTheta2(d2AdTheta2_mem, shape(nvar,nvar,npar,npar), deleteDataWhenDone, ColumnMajorArray<4>());
		static MA4 myd2AdPhi2(d2AdPhi2_mem, shape(nvar,nvar,nvar,nvar), deleteDataWhenDone, ColumnMajorArray<4>());
		static MA4 myd2AdThetadPhi(d2AdThetadPhi_mem, shape(nvar,nvar,npar,nvar), deleteDataWhenDone, ColumnMajorArray<4>());
		static MA4 myd2AdPhidTheta(d2AdPhidTheta_mem, shape(nvar,nvar,nvar,npar), deleteDataWhenDone, ColumnMajorArray<4>());
		static MA4 myd2EdTheta2(d2EdTheta2_mem, shape(nvar,Nreact,npar,npar), deleteDataWhenDone, ColumnMajorArray<4>());
		static MA4 myd2EdPhi2(d2EdPhi2_mem, shape(nvar,Nreact,nvar,nvar), deleteDataWhenDone, ColumnMajorArray<4>());
		static MA4 myd2EdThetadPhi(d2EdThetadPhi_mem, shape(nvar,Nreact,npar,nvar), deleteDataWhenDone, ColumnMajorArray<4>());
		static MA4 myd2EdPhidTheta(d2EdPhidTheta_mem, shape(nvar,Nreact,nvar,npar), deleteDataWhenDone, ColumnMajorArray<4>());

		// construct total second derivative of A and E WRT Theta
		static MA4 d2AdTheta2_tot(nvar,nvar,npar,npar);
		static MA4 d2EdTheta2_tot(nvar,nvar,npar,npar);

		static MA6 vvvvpp(nvar,nvar,nvar,nvar,npar,npar);
		static MA5 vvvpp(nvar,nvar,nvar,npar,npar);
		static MA4 vvvp(nvar,nvar,nvar,npar);
		static MA5 vvvvp(nvar,nvar,nvar,nvar,npar);

		fifthIndex	m; sixthIndex n;
		// total second derivative of A
		// compute dPhidTheta' * d2AdPhi2 * dPhidTheta
		vvvvp 			= Sens_MRE(k,m)*myd2AdPhi2(i,j,k,l);
		vvvp 			= sum(vvvvp(i,j,m,k,l),m);
		vvvpp 			= vvvp(i,j,k,l)*Sens_MRE(k,m);
		vvpp 			= sum(vvvpp(i,j,m,k,l),m);

		d2AdTheta2_tot 	= vvpp;
		vvvpp			= myd2AdPhidTheta(i,j,k,m)*Sens_MRE(k,l);

		d2AdTheta2_tot +=	sum(vvvpp(i,j,m,k,l), m);			// d2AdPhidTheta*dPhidTheta_i
		vvvpp 			=	mydAdPhi(i,j,k)*Sens2_MRE(k,l,m);

		d2AdTheta2_tot += 	sum(vvvpp(i,j,m,k,l), m);			// dAdPhi*d2PhidTheta_i_j

		vvvpp 			=	myd2AdThetadPhi(i,j,l,k)*Sens_MRE(k,m);
		d2AdTheta2_tot += 	sum(vvvpp(i,j,m,k,l),m);		// d2AdTheta_idPhi*dPhidTheta_j
		d2AdTheta2_tot +=	myd2AdTheta2(i,j,k,l);			// d2AdTheta_i_j

		// total second derivative of E
		// compute dPhidTheta' * d2EdPhi2 * dPhidTheta
		vvvvp 			= Sens_MRE(k,m)*myd2EdPhi2(i,j,k,l);
		vvvp 			= sum(vvvvp(i,j,m,k,l),m);
		vvvpp 			= vvvp(i,j,k,l)*Sens_MRE(k,m);
		vvpp 			= sum(vvvpp(i,j,m,k,l),m);

		d2EdTheta2_tot 	= vvpp;
		vvvpp			= myd2EdPhidTheta(i,j,k,m)*Sens_MRE(k,l);
		vvpp			= sum(vvvpp(i,j,m,k,l), m);

		d2EdTheta2_tot +=	vvpp;								// d2AdPhidTheta*dPhidTheta_i
		vvvpp 			=	mydEdPhi(i,j,k)*Sens2_MRE(k,l,m);

		d2EdTheta2_tot += 	sum(vvvpp(i,j,m,k,l), m);			// dAdPhi*d2PhidTheta_i_j

		vvvpp 			=	myd2EdThetadPhi(i,j,l,k)*Sens_MRE(k,m);
		d2EdTheta2_tot += 	sum(vvvpp(i,j,m,k,l),m);		// d2AdTheta_idPhi*dPhidTheta_j
		d2EdTheta2_tot +=	myd2EdTheta2(i,j,k,l);

		// Variance
//		Range all = Range::all();
		vvvpp = d2AdTheta2_tot(i,k,l,m)*V(k,j);
//		cout << "tmp5 j=0" << endl << tmp5(all,all,all,all,0) << endl;
		Sens2_Var_dot 	=	sum(vvvpp(i,j,m,k,l),m);			// d2AdTheta_l_m * V
//		cout << "d2AdTheta_l_m * V" << endl;
//		cout << "Sens2_Var_dot " << endl << Sens2_Var_dot << endl;

		vvvpp = dAdTheta_tot(i,k,l)*Sens_Var(k,j,m);
		Sens2_Var_dot	+= 	sum(vvvpp(i,j,m,k,l),m);			// dAdTheta_l * dVdTheta_m
//		cout << "dAdTheta_l * dVdTheta_m" << endl;
//		cout << "Sens2_Var_dot " << endl << Sens2_Var_dot << endl;

		vvvpp = dAdTheta_tot(i,k,m)*Sens_Var(k,j,l);
		Sens2_Var_dot 	+=	sum(vvvpp(i,j,m,k,l),m); 	// dAdTheta_m * dVdTheta_l
//		cout << "dAdTheta_m * dVdTheta_l" << endl;
//		cout << "Sens2_Var_dot " << endl << Sens2_Var_dot << endl;

//		cout << "A:" << endl << A << endl;
//		cout << "Sens2_Var" << endl << Sens2_Var << endl;
		vvvpp = A(i,k)*Sens2_Var(k,j,l,m);
		Sens2_Var_dot 	+=	sum(vvvpp(i,j,m,k,l),m);	// A*d2VdTheta_l_m (transpose necessary b/c of blitz++ partial reductions, see documentation)

//		cout << "A" << endl << A << endl;
//		cout << "Sens2_Var" << endl << Sens2_Var << endl;
//		cout << "A*d2VdTheta_l_m" << endl;
//		cout << "Sens2_Var_dot " << endl << Sens2_Var_dot << endl;
		vvvpp = Sens2_Var(i,k,l,m)*A(j,k);
		Sens2_Var_dot 	+= 	sum(vvvpp(i,j,m,k,l),m);	// d2VdTheta_l_m*A'

//		cout << "d2VdTheta_l_m*A\'" << endl;
//		cout << "Sens2_Var_dot " << endl << Sens2_Var_dot << endl;

		vvvpp = Sens_Var(i,k,l)*dAdTheta_tot(j,k,m);
		Sens2_Var_dot	+= 	sum(vvvpp(i,j,m,k,l),m);	// dVdTheta_i*dA'dTheta_j
//		cout << "dVdTheta_i*dA\'dTheta_j" << endl;
//		cout << "Sens2_Var_dot " << endl << Sens2_Var_dot << endl;

		vvvpp = Sens_Var(i,k,m)*dAdTheta_tot(j,k,l);
		Sens2_Var_dot	+= 	sum(vvvpp(i,j,m,k,l),m);	// dVdTheta_j*dA'dTheta_i
		vvvpp = V(i,k)*d2AdTheta2_tot(j,k,l,m);
		Sens2_Var_dot 	+=	sum(vvvpp(i,j,m,k,l),m);	// V*d2A'dTheta_i_j
		vvvpp = d2EdTheta2_tot(i,k,l,m)*E(j,k);
		Sens2_Var_dot 	+= 	sum(vvvpp(i,j,m,k,l),m); 	// d2EdTheta_i_j*E'
		vvvpp = dEdTheta_tot(i,k,l)*dEdTheta_tot(j,k,m);
		Sens2_Var_dot 	+= 	sum(vvvpp(i,j,m,k,l),m); 	// dEdTheta_i*dEdTheta_j'
		vvvpp = dEdTheta_tot(i,k,m)*dEdTheta_tot(j,k,l);
		Sens2_Var_dot 	+= 	sum(vvvpp(i,j,m,k,l),m); 	// dEdTheta_j*dEdTheta_i'
		vvvpp = E(i,k)*d2EdTheta2_tot(j,k,l,m);
		Sens2_Var_dot	+= 	sum(vvvpp(i,j,m,k,l),m);	// E*d2EdTheta_i_j'

//		cout << "Sens2_Var_dot " << endl << Sens2_Var_dot << endl;

		// Fundamental matrix
		vvvpp 			= d2AdTheta2_tot(i,k,l,m)*Phi(k,j);
		Xi2_dot 		= sum(vvvpp(i,j,m,k,l),m); 				// d2AdTheta_i_j*Phi
		vvvpp			= dAdTheta_tot(i,k,l)*Xi(k,j,m);
		Xi2_dot			+=	sum(vvvpp(i,j,m,k,l),m);			// dAdTheta_i*dPhidTheta_j
		vvvpp 			= dAdTheta_tot(i,k,m)*Xi(k,j,l);
		Xi2_dot 		+= 	sum(vvvpp(i,j,m,k,l),m);			// dAdTheta_j*dPhidTheta_i
		vvvpp 			= 	A(i,k)*Xi2(k,j,l,m);
		Xi2_dot 		+= 	sum(vvvpp(i,j,m,k,l),m);			// A*d2PhidTheta_i_j

		// pack the NVector with the computed derivatives

		packYSDot_HO(yS_dot, Sens2_MRE_dot, Sens2_Var_dot, Xi2_dot);
//		for (int lPar=0; lPar<npar; lPar++)
//		{
//			cout << lPar << ": ";
//			for (int lVar=0; lVar<nvar; lVar++)
//				cout << setw(20) << NV_Ith_S(yS_dot[lPar],lVar);
//			cout << endl;
//		}
//
//		for (int lPar1=0; lPar1<npar; lPar1++)
//			for (int lPar2=0; lPar2<npar; lPar2++)
//			{
//				cout << "(" << lPar1 << ", " << lPar2 << "): " ;
//				for (int lVar=0; lVar<nvar; lVar++)
//					cout << setw(20) << NV_Ith_S(yS_dot[lPar1*npar+lPar2+1],lVar);
//				cout << endl;
//			}
	}

	return 0;
}

void LNA::unpackYDot(N_Vector ydot, double *phi, MA2 &V, MA2 &Phi) {
	// get the parameters from the current solution vector
	// extract MRE

	//	const int *dims = V.getDims();
	const TinyVector<int,2> dims = V.shape();
	int nvar  = dims(0);

	for (int lVar1=0; lVar1<nvar; lVar1++) {
		phi[lVar1] = NV_Ith_S(ydot, lVar1);
	}

#ifdef DEBUG
	cout << "phi " << phi[0] << " " << phi[1] <<endl;
#endif

	int ix=nvar,k=0;

	// extract the variance solution
	for (int lVar1=0; lVar1<nvar; lVar1++) {
		for (int lVar2=0; lVar2<=lVar1; lVar2++) {
			V(lVar2,lVar1) = NV_Ith_S(ydot,ix+k);
			if (lVar2!=lVar1)
				V(lVar1,lVar2) = NV_Ith_S(ydot,ix+k);
			k++;
		}
	}
	ix += k;

#ifdef DEBUG
	cout << "current variance " << endl;
	//	V.print();
	cout << V;
#endif
	// extract the fundamental matrix solution
	k=0;
	for (int lVar1=0; lVar1<nvar; lVar1++)
		for (int lVar2=0; lVar2<nvar; lVar2++) {
			Phi(lVar2,lVar1) = NV_Ith_S(ydot, ix+k);
			k++;
		}
#ifdef DEBUG
	cout << "current fund matrix " << endl;
	//	Phi.print();
	cout << Phi;
#endif



}

//void packYDot(N_Vector ydot, const double *netChangeVec, const MA2 &dVdt, const MA2 &dPhidt)
void LNA::packYDot(N_Vector ydot, const Vector &netChangeVec, const MA2 &dVdt, const MA2 &dPhidt)
{
	const TinyVector<int,2> dims = dVdt.shape();

	int nvar  = dims(0);

	// MRE
	for (int lVar1=0; lVar1<nvar; lVar1++) {
		NV_Ith_S(ydot,lVar1) = netChangeVec(lVar1);
	}

	// load dVdt into the RHS vector
	int	k=0, ix;
	for (int lVar1=0; lVar1<nvar; lVar1++)
		for (int lVar2=0; lVar2<=lVar1; lVar2++) {
			NV_Ith_S(ydot, nvar+k) = dVdt(lVar2,lVar1);
			k++;
		}

	// derivative of fundamental matrix
	k=0;
	ix=nvar*(nvar+3)/2;
	for (int lVar1=0; lVar1<nvar; lVar1++)
		for (int lVar2=0; lVar2<nvar; lVar2++) {
			NV_Ith_S(ydot, ix+k) = dPhidt(lVar2,lVar1);
			k++;
		}

}

void LNA::unpackYSDot(const N_Vector *yS, MA2 &Sens, MA3 &dSigma, MA3 &Xi ) {
	// extract the current sensitivity
	// Sensitivities of the macroscopic rate equation
	const TinyVector<int,2> dims = Sens.shape();
	int nvar  = dims(0), npar = dims(1);

	for (int i=0; i<npar; i++) {
		for (int j=0; j<nvar; j++) {
			Sens(j,i) =  NV_Ith_S(yS[i], j);
		}
	}

	int ix,k;

	// extract the derivatives of the variances
	for (int lPar=0; lPar<npar; lPar++){
		k=0;
		for (int lVar1=0; lVar1<nvar; lVar1++) {
			for (int lVar2=0; lVar2<=lVar1; lVar2++) {
				dSigma(lVar2,lVar1,lPar) =  NV_Ith_S(yS[lPar], nvar+k);
#ifdef DEBUG
				cout << "NV_Ith_S(yS[" << lPar<< "]," <<  nvar+k<< "): " << NV_Ith_S(yS[lPar], nvar+k) << endl;
#endif
				if (lVar2!=lVar1)
					dSigma(lVar1,lVar2,lPar) =  NV_Ith_S(yS[lPar], nvar+k);
				k++;
			}
		}
	}

	// extract the derivatives of the fundamental matrix
	ix=nvar*(nvar+3)/2;
	for (int lPar=0; lPar<npar; lPar++) {
		k=0;
		for (int lVar1=0; lVar1<nvar; lVar1++)
			for (int lVar2=0; lVar2<nvar; lVar2++) {
				Xi(lVar2,lVar1,lPar) = NV_Ith_S(yS[lPar], ix+k);
				k++;
			}
	}

}

void LNA::packYSDot(N_Vector *ySdot, const MA2 &dSdt, const MA3 &d_dt_dSigma_dTheta, const MA3 &dXi_dt)
{
	//	const int *dims = d_dt_dSigma_dTheta.getDims();
	const TinyVector<int,3> dims = d_dt_dSigma_dTheta.shape();
	int nvar  = dims(0), npar = dims(2);

	// MRE
	for (int lPar=0; lPar<npar; lPar++)
		for (int lVar=0; lVar<nvar; lVar++) {
			NV_Ith_S(ySdot[lPar],lVar) = dSdt(lVar,lPar);
		}

	// Variances
	int ix,k;
	for (int lPar=0; lPar<npar; lPar++){
		k=0;
		for (int lVar1=0; lVar1<nvar; lVar1++) {
			for (int lVar2=0; lVar2<=lVar1; lVar2++) {
				NV_Ith_S(ySdot[lPar], nvar+k) = d_dt_dSigma_dTheta(lVar2,lVar1,lPar);
				k++;
			}
		}
	}

	// Fund
	// extract the derivatives of the fundamental matrix
	ix=nvar*(nvar+3)/2;
	for (int lPar=0; lPar<npar; lPar++) {
		k=0;
		for (int lVar1=0; lVar1<nvar; lVar1++)
			for (int lVar2=0; lVar2<nvar; lVar2++) {
				NV_Ith_S(ySdot[lPar], ix+k) = dXi_dt(lVar2,lVar1,lPar);
				k++;
			}
	}

}

void LNA::unpackYS_HO(N_Vector *yS, MA3 &Sens2_MRE, MA4 &Sens2_Var, MA4 &Xi2) {
	// extract the current sensitivity

	const TinyVector<int,3> dims = Sens2_MRE.shape();
	int nvar  = dims(0), npar = dims(1);

	// Sensitivities of the macroscopic rate equation

	for (int lPar1=0; lPar1<npar; lPar1++)
		for (int lPar2=0; lPar2<npar; lPar2++)
		{
			for (int lVar=0; lVar<nvar; lVar++)
			{
				Sens2_MRE(lVar,lPar1,lPar2) = NV_Ith_S(yS[lPar1*npar + lPar2 + npar], lVar);
			}
		}

	// extract the derivatives of the variances
	int k;

	for (int lPar1=0; lPar1<npar; lPar1++)
		for (int lPar2=0; lPar2<npar; lPar2++)
		{
			k=0;
			for (int lVar1=0; lVar1<nvar; lVar1++)
				for (int lVar2=0; lVar2<=lVar1; lVar2++) {
					Sens2_Var(lVar2,lVar1,lPar1,lPar2) 		= NV_Ith_S(yS[lPar1*npar + lPar2 + npar], k+nvar);
					if (lVar2!=lVar1)
						Sens2_Var(lVar1,lVar2,lPar1,lPar2) 	= NV_Ith_S(yS[lPar1*npar + lPar2 + npar], k+nvar); // symmetry
					k++;
				}
		}

	// extract the derivatives of the fundamental matrix
	int ix=nvar*(nvar+3)/2;
	for (int lPar1=0; lPar1<npar; lPar1++)
		for (int lPar2=0; lPar2<npar; lPar2++)
			{
			k=0;
			for (int lVar1=0; lVar1<nvar; lVar1++)
				for (int lVar2=0; lVar2<nvar; lVar2++) {
					Xi2(lVar2,lVar1,lPar1,lPar2) = NV_Ith_S(yS[lPar1*npar + lPar2 + npar], ix+k);
					k++;
			}
	}
}

void LNA::packYSDot_HO(N_Vector *ySdot, MA3 &Sens2_MRE_dot, MA4 &Sens2_Var_dot, MA4 &Xi2_dot)
{
	const TinyVector<int,3> dims = Sens2_MRE_dot.shape();
	int nvar  = dims(0), npar = dims(2);

	// MRE sensitivity
	for (int lPar1=0; lPar1<npar; lPar1++)
		for (int lPar2=0; lPar2<npar; lPar2++)
		{
			for (int lVar=0; lVar<nvar; lVar++)
			{
				NV_Ith_S(ySdot[lPar1*npar + lPar2 + npar], lVar) = Sens2_MRE_dot(lVar,lPar1,lPar2);
			}
		}

	// variances sensitivity
	int k;
	for (int lPar1=0; lPar1<npar; lPar1++)
		for (int lPar2=0; lPar2<npar; lPar2++)
		{
			k=0;
			for (int lVar1=0; lVar1<nvar; lVar1++)
				for (int lVar2=0; lVar2<=lVar1; lVar2++) {
					NV_Ith_S(ySdot[lPar1*npar + lPar2 + npar], k+nvar) = Sens2_Var_dot(lVar2,lVar1,lPar1,lPar2);
					k++;
				}
		}

	// Fund sensitivity
	int ix=nvar*(nvar+3)/2;
	for (int lPar1=0; lPar1<npar; lPar1++)
		for (int lPar2=0; lPar2<npar; lPar2++)
			{
			k=0;
			for (int lVar1=0; lVar1<nvar; lVar1++)
				for (int lVar2=0; lVar2<nvar; lVar2++) {
					NV_Ith_S(ySdot[lPar1*npar + lPar2 + npar], ix+k) = Xi2_dot(lVar2,lVar1,lPar1,lPar2);
					k++;
			}
	}
}

void LNA::initCVODES() {
	/* CVODES set up 												*/

	y = abstol = NULL;
	cvode_mem = 0;

	// initial conditions vector
	I = N_VNew_Serial(RHS_SIZE);
	// output of CVode at each iteration of solver
	x = N_VNew_Serial(RHS_SIZE);

	// allocate vectors for CVODES
	y = N_VNew_Serial(RHS_SIZE);
	if (check_flag((void *)y, "N_VNew_Serial", 0))
		throw runtime_error("N_VNew_Serial");

	abstol = N_VNew_Serial(RHS_SIZE);
	if (check_flag((void *)abstol, "N_VNew_Serial", 0))
		throw runtime_error("N_VNew_Serial");

	// Set the scalar relative tolerance
	reltol = RTOL;
	/* Set the vector absolute tolerance */
	for (int i=0; i<RHS_SIZE; i++) // ix = total number integration variables
		NV_Ith_S(abstol,i) = ATOL;

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0))
		throw runtime_error("CVodeCreate");

	 Call CVodeInit to initialize the integrator memory and specify the
		 * user's right hand side function in y'=f(t,y), the inital time t0, and
		 * the initial dependent variable vector y.
	realtype t0=0.0; // gets reinitialized when computeLinearNoise() is called.
	flag = CVodeInit(cvode_mem, &fundRHS, t0, y);
	if (check_flag(&flag, "CVodeInit", 1))
		throw runtime_error("CVodeInit");

		// max number of convergence failures?
	CVodeSetMaxConvFails(cvode_mem, 10);


	 Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances
	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1))
		throw runtime_error("CVodeSVtolerances");*/

	// SPECIFY THE LINEAR SOLVER

	/*	USE THE DENSE LINEAR SOLVER */



//	/* Set the preconditioner matrix */
//	flag = CVSpilsSetPreconditioner(cvode_mem, NULL, Preconditioner);
//	if (check_flag(&flag, "CVSpilsSetPreconditioner", 1))
//		throw runtime_error("CVSpilsSetPreconditioner");


	/* setup Sensitivity Analysis */

	/* construct vectors for the initial condition of the
	 * sensitivities.  one for first order system, another for
	 * second order system.
	 */
	yS0_1 = N_VCloneVectorArray_Serial(npar, y);
	yS0_2 = N_VCloneVectorArray_Serial((npar+1)*npar, y);
	yS0 = NULL; // set to the first order sens.

	// output vectors
	ySout_1 = N_VCloneVectorArray_Serial(npar, y);
	ySout_2 = N_VCloneVectorArray_Serial(npar*(npar+1), y);
	ySout = NULL;

	// set up absolute tolerances for sensitivity equations
	// first order
	abstol_vec1 = new realtype[npar];
	for (int i=0; i<npar; i++)
		abstol_vec1[i] = ATOLS;

	abstol_vec2 = new realtype[npar*(npar+1)];
	for (int i=0; i<(npar*(npar+1)); i++)
		abstol_vec2[i] = ATOLS;

}


void LNA::setupCVODES(const double t0, const parameters &pars) {

	/* Reinitialize the ODE solver with the correct initial conditions,
	 * parameters, and starting time*/

	if (cvode_mem != 0)
		CVodeFree(&cvode_mem);

	/* Call CVodeCreate to create the solver memory and specify the
		 * Backward Differentiation Formula and the use of a Newton iteration */
	// begin
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0))
		throw runtime_error("CVodeCreate");

	/* Call CVodeInit to initialize the integrator memory and specify the
		 * user's right hand side function in y'=f(t,y), the inital time t0, and
		 * the initial dependent variable vector y. */
//	realtype t0=0.0; // gets reinitialized when computeLinearNoise() is called.
	flag = CVodeInit(cvode_mem, &fundRHS, t0, y);
	if (check_flag(&flag, "CVodeInit", 1))
		throw runtime_error("CVodeInit");

		// max number of convergence failures?
	CVodeSetMaxConvFails(cvode_mem, 10);


	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1))
		throw runtime_error("CVodeSVtolerances");
	// end
//
//	flag = CVodeReInit(cvode_mem, t0, y);
//	if(check_flag(&flag, "CVodeReInit", 1))
//		throw runtime_error("CVodeReInit");

	/* Attach User Data (theta etc)*/
	flag = CVodeSetUserData(cvode_mem, (void*)&pars);
	if (check_flag(&flag, "CVodeSetUserData", 1))
		throw runtime_error("CVodeSetUserData");

	myUserData_ptr = &((CVodeMem)cvode_mem)->cv_user_data;
	myUserData = ((CVodeMem)cvode_mem)->cv_user_data;
//	void *tmp = &((CVodeMem)cvode_mem)->cv_user_data;

	// Use the Dense Linear Solver
	//flag = CVDense(cvode_mem, RHS_SIZE);
	//if (check_flag(&flag, "CVDense", 1))
	//	throw runtime_error("CVDense");

	//// Specify the Jacobian function
	//flag = CVDlsSetDenseJacFn(cvode_mem, &Jac);
	//if (check_flag(&flag, "CVDlsSetDenseJacFn", 1))
	//	throw runtime_error("CVDlsSetDenseJacFn");



	/* Use the Krylov subspace GMRES method */
	flag = CVSpgmr(cvode_mem, PREC_NONE, 0); // 0=default dimension 5  of subspace
	if (check_flag(&flag, "CVSpgmr", 1))
		throw runtime_error("CVSpgmr");
//
//	/* Set the preconditioner matrix */
//	flag = CVSpilsSetPreconditioner(cvode_mem, NULL, Preconditioner);
//	if (check_flag(&flag, "CVSpilsSetPreconditioner", 1))
//		throw runtime_error("CVSpilsSetPreconditioner");

	/* check if the state of the sensitivity computations changed
	 * since the last time the solver was called.
	 * if so, reinit the sens. solver, or disable it as necessary.
	 */

	/* initialize sensitivity solver */
	if (pars.computeSens2) {
		yS0 		= yS0_2;
		abstol_vec 	= abstol_vec2;
		ySout		= ySout_2;
		flag = CVodeSensInit(cvode_mem, npar*(npar+1), CV_SIMULTANEOUS, &sensRhs, yS0);
		if (check_flag(&flag, "CVodeSensInit", 1))
			throw runtime_error("CVodeSensInit");
		// set the tolerances
		flag = CVodeSensSStolerances(cvode_mem, reltol, abstol_vec);
		if(check_flag(&flag, "CVodeSensSStolerances", 1))
			throw runtime_error("Error in CVodeSensSStolerances");
	} else if (pars.computeSens) {
		yS0 		= yS0_1;
		abstol_vec 	= abstol_vec1;
		ySout		= ySout_1;

		flag = CVodeSensInit(cvode_mem, npar, CV_SIMULTANEOUS, &sensRhs, yS0);
		if (check_flag(&flag, "CVodeSensInit", 1))
			throw runtime_error("CVodeSensInit");

		flag = CVodeSensSStolerances(cvode_mem, reltol, abstol_vec);
		if(check_flag(&flag, "CVodeSensSStolerances", 1))
			throw runtime_error("Error in CVodeSensSStolerances");
	}

	/* update the state of the sensitivity computations for future calls */
	computeSens 	= pars.computeSens;
	computeSens2 	= pars.computeSens2;
}

void LNA::setupSens(SS_FLAG SS_flag, MA3 &Sens_MRE, MA4 &Sens_Var, const parameters &pars) {
	// initialize the first order sensitivity of the system
	// sets up Sens_MRE and Sens_Var
	const double *Theta 		= pars.Theta;

	int k;

	/* initial first order sensitivities of the MRE 				*/
	double *S_0 = new double[nvar*npar];
	if ( (SS_flag == SS_V0) || (SS_flag == SS_NONE) )
	{
		/* the initial conditions were specified, so the initial
		sensitivities of the MRE are zero */
//		for ()
		memset(S_0, 0, sizeof(double)*nvar*npar);
		for (int i=0; i<nvar; i++)
			for (int j=0; j<nvar; j++)
				Sens_MRE(j,npar-nvar+i,0) = (i==j); // sens of init cond wrt init cond is 1
	} else {
		// compute steady state for MRE, store in S_0
		S0(Theta, S_0);
	}

	k=0;
	// load the sensitivities of the MRE
	for (int i=0; i<npar-nvar; i++)
	{
		for (int j=0; j<nvar; j++) {
			Sens_MRE(j,i,0) = S_0[k];
			k++;
		}
	}
	delete[] S_0;

	/* initial first order sensitivities of the variance 			*/
	double *SV_0 = new double[ nvar*(nvar+1)/2 * npar];
	if ( (SS_flag==SS_Y0) || (SS_flag==SS_NONE) ) {
		// variance was specified explicitly, so sensitivity is zero at t=0
		memset(SV_0, 0, sizeof(double)*nvar*(nvar+1)/2 * npar);
	} else {
		// assume steady state variance
		SV0(Theta, SV_0);
	}

	k=0;
	for (int lPar=0; lPar<npar; lPar++)
		for (int lVar1=0; lVar1<nvar; lVar1++)
			for (int lVar2=0; lVar2<=lVar1; lVar2++)
			{
				Sens_Var(lVar2,lVar1,0,lPar) = SV_0[k];
				if (lVar1!=lVar2)
					Sens_Var(lVar1,lVar2,0,lPar) = SV_0[k]; // symmetric
				k++;
			}
	delete[] SV_0;

	iout = 0;  //tout = tspan[1];
}

void LNA::setupSens2(SS_FLAG SS_flag, MA4 &Sens2_MRE, MA5 &Sens2_Var, const parameters &pars) {
	// initialize the first order sensitivity of the system
	// sets up Sens2_MRE and Sens2_Var

	const double *Theta = pars.Theta;

	int k;
	// initialize higher order sensitivites
	//% initial sensitivities

	/* initial second order sensitivities of the MRE 				*/
	double *S2_0 = new double[nvar*npar*npar];
	if ( (SS_flag==SS_Y0) || (SS_flag==SS_BOTH) ) {
		// the system is initially in steady-state
		// compute the second order sensitivities
		S20(Theta, S2_0);
	}
	else{
		// Y0 was specified; the initial conditions don't depend on the params
		memset(S2_0, 0.0, sizeof(double)*nvar*npar*npar);
	}

	k=0;
	// load the second order sensitivities of the MRE
	for (int lPar1=0; lPar1<npar; lPar1++)
		for (int lPar2=0; lPar2<npar; lPar2++)
			for (int lVar=0; lVar<nvar; lVar++) {
				Sens2_MRE(lVar,lPar1,lPar2,0) = S2_0[k];
				k++;
			}
	delete[] S2_0;

	// initial second order sensitivities of the variance
	double *S2V_0 = new double[ nvar*(nvar+1)/2 *npar*npar];

	if ( (SS_flag == SS_V0) || (SS_flag == SS_BOTH)) {
		// assume steady state for the variances
		S2V0(Theta, S2V_0);
	}
	else {
		memset(S2V_0, 0.0, sizeof(double)*nvar*(nvar+1)/2 *npar*npar);
	}

	/* load the values into Sens2_Var */
	k=0;
	for (int lPar2=0; lPar2<npar; lPar2++)
		for (int lPar1=0; lPar1<npar; lPar1++)
			for (int lVar1=0; lVar1<nvar; lVar1++)
				for (int lVar2=0; lVar2<=lVar1; lVar2++)
				{
					Sens2_Var(lVar2,lVar1,0,lPar1,lPar2) 		= S2V_0[k];
					if (lVar1!=lVar2)
						Sens2_Var(lVar1,lVar2,0,lPar1,lPar2) 	= S2V_0[k]; // symmetric
					k++;
				}
	delete[] S2V_0;

	iout = 0;  //tout = tspan[1];

}
