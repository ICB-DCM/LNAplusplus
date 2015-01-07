/*
 * datatypes.h
 *
 *  Created on: Oct 18, 2013
 *      Author: justinfeigelman
 */

#ifndef DATATYPES_H_
#define DATATYPES_H_

//#include "MultiArray.h"
#include <blitz/array.h>

using namespace blitz;

typedef Array<double,1> Vector;
typedef Array<double,2> MA2;
typedef Array<double,3> MA3;
typedef Array<double,4> MA4;
typedef Array<double,5> MA5;
typedef Array<double,6> MA6;


// will pass an object of type LNA into the usr_data for cvode
class LNA;

struct outputStruct {

	MA4 *Sigma;
	MA5 *dSigma;
	MA2 *Y;
//	MA3 *Sens;
	MA3 *traj_deriv;
	MA4 *Sens2_MRE;
	MA6 *Sens2_Var;

	outputStruct() {
		Sigma=0;
		dSigma=0;
		Y=0;
		traj_deriv=0;
		Sens2_MRE=0;
		Sens2_Var=0;
	}

};

// CVODES user data
typedef struct {
	const int nvar, npar;
	const MA2 	*S;
	const double 	*Theta;
	const bool computeSens, computeSens2;

	LNA *lna;

} parameters;

#endif /* DATATYPES_H_ */
