///*
// * main.cpp
// *
// *  Created on: Oct 18, 2013
// *      Author: justinfeigelman
// */
//
//

//-fsanitize=address -fno-omit-frame-pointer

#include "computeLinearNoise.h"
#include <iostream>
#include <vector>

#include "cvodes/cvodes.h"
#include "nvector/nvector_serial.h"
#include "blitz/array.h"
#include "MODEL_DEF.h"

using namespace std;
using namespace blitz;

int main(int argc, char **argv ) {

	outputStruct os;

	int NvarObs = NVAR;
	Array<int,1> varObs(NVAR);
	varObs(0) = 0;
	varObs(1) = 0;
	varObs(2) = 0;
	varObs(3) = 0;
	varObs(4) = 0;


	double merr[] = {0.0, 0.0, 0.0, 0.0, 0.0}; // measurement error

	int N = 10;
	double tFinal = 1.0;

	Vector tspan(N);

	for (int i=0; i<N; i++)
		tspan(i) = (double)i*tFinal/(double)N;

	double Theta[] = {0.1, 0.7, 0.35, 0.3, 0.1, 0.9, 0.2, 0.1};

//	double y0[] = {20.0, 1.0, 1.0, 1.0, 1.0}; // initial conditions (WORKS)
	double y0[] = {20.0, 0.0, 0.0, 0.0, 0.0}; // initial conditions (THROWS EXCEPTION)
	double v0[15]; //  5*(5+1)/2 covariance terms

	for (int i=0; i<15; i++) {
		v0[i] = 0.0;
	}

	MA2 S(NVAR,NREACT);
	S = 	STOICH_MAT;

	int npar=NPAR, nvar=NVAR;

	os.Y 			= new MA2(NvarObs,N);
	os.traj_deriv	= new MA3(NvarObs, npar, N);
	os.Sigma 		= new MA4(NvarObs, NvarObs, N, N);
	os.dSigma		= new MA5(NvarObs, NvarObs, N, N, npar);

	// Second Order Sensitivities
	os.Sens2_MRE 	= new MA4(NvarObs, npar, npar, N);
	os.Sens2_Var 	= new MA6(NvarObs, NvarObs, npar, npar, N, N);

	// initialize all matrices to zero
	(*os.Y)=0;
	(*os.traj_deriv)=0;
	(*os.Sigma)=0;
	(*os.dSigma)=0;
	(*os.Sens2_MRE)=0;
	(*os.Sens2_Var)=0;

	LNA lna(nvar,npar,S);

//	LNA::SS_FLAG SS = LNA::SS_BOTH;
//	LNA::SS_FLAG SS = LNA::SS_V0;
	LNA::SS_FLAG SS = LNA::SS_NONE;



//	cout << "cs1: " << computeSens << " cs2: " << computeSens2 << endl;
	cout << "y0 ";
	for (int i=0; i<nvar; i++) {
		cout << setw(5) << y0[i];
	}

	cout << endl << "V0 ";
	for (int i=0; i<nvar; i++) {
		cout << setw(5) << v0[i];
	}

	cout << endl << "Theta ";
	for (int i=0; i<npar; i++) {
		cout << setw(5) << Theta[i];
	}
	cout << "SS_flag " << SS << endl;
	cout << "Time: ";
	for (int i=0; i<N; i++)
		cout << setw(5) << tspan(i) << endl;

	cout << endl << "var obs " ;
	for (int i=0; i<NvarObs; i++)
		cout << setw(5) << varObs(i);
	cout << endl;
//	for (int i=0; i<NvarObs; i++)
//			cout << setw(5) << merr[i];



	try {
//		lna.computeLinearNoise(y0, v0, Theta, true, true,
//				SS, os, tspan, varObs, &merr);
		lna.computeLinearNoise(y0, v0, Theta, true, true,
				SS, os, tspan, varObs, merr);



	} catch (int e)
	{
		cerr << "Caught exception " << e << endl;
	};

	cout << "MRE trajectory:" << endl;
	cout << *os.Y << endl;
//
	cout << "Variance:" << endl;
	cout << *os.Sigma << endl;
//
	cout << "Sensitivity MRE " << endl << *os.traj_deriv << endl;
//	Range all = Range::all();
//
//	for (int i=0; i < pars.npar; i++)
//		cout << "Sensitivity Var(" << i << ") " << endl << (*os.dSigma)(all,all,all,all,i) << endl;


// output whole thing as one messy array
//	cout << "second order sens. MRE " << endl << *os.Sens2_MRE << endl;

//	 output just time point 0
//	cout << "second order sens MRE (t=0) protein only:" << endl << (*os.Sens2_MRE)(Range::all(), Range::all(), Range::all(), 0) << endl;
	cout << "second order sens MRE (t=0) protein only:" << endl << (*os.Sens2_MRE)(Range::all(), 0, 1, Range::all()) << endl;

//	for (int i=0; i<N; i++)
//		for (int j=0; j<N; j++)
//			cout << "second order sens. Var(" << i << "," << j << ")" << endl <<
//				 (*os.Sens2_Var)(all,all,all,all,i,j) << endl;

	delete os.Y;
	delete os.traj_deriv;
	delete os.Sigma;
	delete os.dSigma;
	delete os.Sens2_MRE;
	delete os.Sens2_Var;

//	extern void *cvode_mem;
//	extern N_Vector y, abstol, x, I;
//	extern N_Vector *yS0, *ySout;
//


	return 0;
}
