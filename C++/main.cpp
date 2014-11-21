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

using namespace std;
using namespace blitz;

int main(int argc, char **argv ) {

	outputStruct os;

	int NvarObs = 2;
	Array<int,1> varObs(2);
	varObs(0) = 0;
	varObs(1) = 1;


	double merr[] = {0., 0.}; // measurement error

	int N = 101;
	Vector tspan(N);

	for (int i=0; i<N; i++)
		tspan(i) =i*0.1;

	double Theta[] = { 20., 25., 10., 1., 2., 200.};
	double y0[] = {20.,200.}; // initial conditions
	double v0[] = {2.,5.,150.};

	MA2 S(2,4);
	S = 	1, -1, 0, 0,
			0, 0, 1, -1;

	int npar=6, nvar=2;

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

	LNA lna(2,6,S);

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
