/*
 * computeLinearNoise.h
 *
 *  Created on: Oct 18, 2013
 *      Author: justinfeigelman
 */

#ifndef COMPUTELINEARNOISE_H_
#define COMPUTELINEARNOISE_H_

#include "datatypes.h"
#include <vector>

// CVODES
#include "cvodes/cvodes.h"
#include "nvector/nvector_serial.h"
#include "cvodes/cvodes_dense.h"
#include "cvodes/cvodes_spgmr.h"

//#include "gsl/gsl_matrix.h"

class LNA
{

public:

	LNA(const int nvar_, const int npar_, const MA2 &S_) :
	nvar(nvar_), npar(npar_), S(S_), RHS_SIZE(nvar_*(nvar_+3)/2 + nvar_*nvar_)
	{
		initCVODES();

		// default: no sensitivity computations
		computeSens=false;
		computeSens2=false;
	};

	~LNA() {
		N_VDestroy(I);
		N_VDestroy(x);
		N_VDestroy(y);
		N_VDestroy(abstol);
		N_VDestroyVectorArray(yS0_1, npar);
		N_VDestroyVectorArray(yS0_2, (npar+1)*npar);
		N_VDestroyVectorArray(ySout_1, npar);
		N_VDestroyVectorArray(ySout_2, (npar+1)*npar);
		CVodeFree(&cvode_mem);

		delete[] abstol_vec1;
		delete[] abstol_vec2;

		//gsl_matrix_free(myJ);
	}

	enum SS_FLAG {
		SS_Y0,
		SS_V0,
		SS_BOTH,
		SS_NONE
	};

	int computeLinearNoise(const double* y0, const double *V0,
			const double *Theta, const bool computeSens, const bool computeSens2,
			const SS_FLAG SS_flag, outputStruct &output, const Vector &tspan,
			const Array<int,1> &varObs, const double *merr);

	int check_flag(void *flagvalue, const char *funcname, int opt);

	// CVODEs function types must be static!
	static int fundRHS(realtype t, N_Vector yIn, N_Vector ydot, void *user_data);
	static int Jac(long int N, realtype t,
			N_Vector y, N_Vector fy, DlsMat J, void *user_data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
	static int Preconditioner(realtype t, N_Vector y, N_Vector fg, N_Vector r, N_Vector z,
			realtype gamma, realtype delta, int lr, void *user_data, N_Vector tmp);
	static int Preconditioner_diag(realtype t, N_Vector y, N_Vector fg, N_Vector r, N_Vector z,
			realtype gamma, realtype delta, int lr, void *user_data, N_Vector tmp);

	static int PreconditionerSetup(realtype t, N_Vector y, N_Vector fy, int jok, int *jcurPtr,
			realtype gamma, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

	static int sensRhs(int Ns, realtype t, N_Vector y, N_Vector ydot,
			N_Vector *yS, N_Vector *yS_dot,
			void *user_data, N_Vector tmp1, N_Vector tmp2);

	static void packYDot(N_Vector y, const Vector&, const MA2 &dVdt, const MA2 &dPhidt);
	static void unpackYDot(N_Vector y, double *phi, MA2 &V, MA2 &Phi);
	static void packYSDot(N_Vector *ySdot, const MA2 &dSdt, const MA3 &d_dt_dSigma_dTheta, const MA3 &dXi_dt);
	static void unpackYSDot(const N_Vector *yS, MA2 &Sens, MA3 &dSigma, MA3 &Xi );
	static void unpackYS_HO(N_Vector *yS, MA3 &Sens2_MRE, MA4 &Sens2_Var, MA4 &Xi2);
	static void packYSDot_HO(N_Vector *ySdot, MA3 &Sens2_MRE_dot, MA4 &Sens2_Var_dot, MA4 &Xi2_dot);

	void initCVODES();
	void setupCVODES(const double t0, const parameters &pars);
	void setupSens(SS_FLAG SS_flag, MA3 &Sens_MRE, MA4 &Sens_Var, const parameters &pars) ;
	void setupSens2(SS_FLAG SS_flag, MA4 &Sens2_MRE, MA5 &Sens2_Var, const parameters &pars);



private:

	const int nvar, npar, RHS_SIZE;
	const MA2 S;

	//static gsl_matrix *myJ; // system Jacobian

	bool computeSens, computeSens2;

	// CVODES globals
	 void 		*cvode_mem;
	 N_Vector 	y, abstol, x, I;
	 N_Vector 	*yS0, // initial conditions of sensitivity
	 	 *yS0_1, // size of first order sens.
	 	 *yS0_2, // size of second ordr sens.
	 	 *ySout, // output of sens.
	 	 *ySout_1, // first order
	 	 *ySout_2;

	 int flag, iout; //flagr,
//	 int rootsfound[2];
	 realtype reltol; //, tout;
	 realtype *abstol_vec, *abstol_vec1, *abstol_vec2;

};



#endif /* COMPUTELINEARNOISE_H_ */
