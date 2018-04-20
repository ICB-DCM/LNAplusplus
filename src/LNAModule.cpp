#include "Python.h"
#include "computeLinearNoise.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "cvodes/cvodes.h"
#include "nvector/nvector_serial.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

#include "MODEL_DEF.h"

//extern void *cvode_mem;
//extern N_Vector y, abstol, x, I;
//extern N_Vector *yS0, *ySout;

extern int main(int argc, char **argv );
static PyObject *LNAError;

static int assignDoubleFromPyList(PyObject *inList, double *&outList) {
	if (PyList_Check(inList)) {
		int n=PyList_Size(inList);
		outList = new double[n];

		for (int i=0; i<n; i++) {
			PyObject *list_item = PyList_GetItem(inList,i);
			if (PyLong_Check(list_item)) {
				// is an integer, convert to double
				outList[i] = PyLong_AsDouble(list_item);
			} else if (PyFloat_Check(list_item)) {
				// float type
				outList[i] = PyFloat_AsDouble(list_item);
			} else {
				// incorrect type, throw error
				PyErr_SetString(LNAError, "Must supply numeric value");
			}
		}
		return n;
	}
	else {
		// wasn't a list, just a single value
		if (PyLong_Check(inList)) {
			// is an integer, convert to double
			outList = new double;
			*outList = PyLong_AsDouble(inList);
		} else if (PyFloat_Check(inList)) {
			// float type
			outList = new double;
			*outList = PyFloat_AsDouble(inList);
		} else {
			// incorrect type, throw error
			PyErr_SetString(LNAError, "Must supply numeric value");
		}

		return 1;
	}
}

static void assignPyListFromDouble(double *inList, PyObject *&outList, int n) {
	outList = PyList_New(n);
	for (int i=0; i<n; i++) {
		PyList_SetItem(outList, i, PyFloat_FromDouble(inList[i]));
	}
}


static int assignIntFromPyList(PyObject *inList, int *&outList) {
	if (PyList_Check(inList)) {
		int n=PyList_Size(inList);
		outList = new int[n];
		for (int i=0; i<n; i++)
		{
			outList[i] = PyLong_AsLong(PyList_GetItem(inList,i));
		}
		return n;
	}
	else {
		outList = new int;
		*outList = PyLong_AsLong(inList);
		return 1;
	}
}

template <class T>
static void copyOut(PyObject *obj, T *data_ptr) {
    PyArrayObject* array = reinterpret_cast<PyArrayObject*>(obj);
    void *out_data = PyArray_DATA(array);
    memcpy( out_data, (void*)data_ptr->data(), sizeof(double) * data_ptr->size());
}


//extern int main(int argc, char **argv );
static PyObject *
LNA_LNA(PyObject *self, PyObject *args, PyObject *kwds)
{
	outputStruct os;
	double S_MAT[NREACT*NVAR] = {STOICH_MAT}; //{1,-1,0,0,0,0,1,-1};
	MA2 S(S_MAT, shape(NVAR,NREACT), neverDeleteData);

	const int Nreact = NREACT, nvar = NVAR;
	const int npar = NPAR+NVAR; // for the init. conditions

    int nObsVar;
    Array<int,1> varObs;

    Vector tspan;
    int N;

    PyObject *_y0=NULL, *_V0=NULL, *_Theta=NULL, *_merr=NULL, *_time=NULL, *_obsVar=NULL,
    		*_computeSens=0, *_computeSens2=0;

    bool computeSens=false, computeSens2=false;

	/* Expected input:
	 * 		Theta: numeric list
	 * 		time: numeric list
	 * 		y0: numeric list
	 * 		V0: numeric list
	 * 		merr: numeric
	 * 		obsVar: numeric list
	 * 		computeSens: boolean
	 * 		computeSens2: boolean
	*/
    static const char* argnames[] = {"Theta", "time", "Y0", "V0", "merr", "obsVar",
    		"computeSens", "computeSens2", NULL};

    // check input arguments format and parse
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|OOOO$OO", const_cast<char **>(argnames), &_Theta, &_time, &_y0, &_V0, &_merr,
    		 &_obsVar, &_computeSens, &_computeSens2))
    {
    	PyErr_SetString(LNAError, "Usage string.");
    	return NULL;
    }

    /* validate input arguments */
    double *y0=NULL, *V0=NULL;
    double *Theta=NULL;
    double *merr=NULL;

    /* check init. cond. ss. */

   LNA::SS_FLAG SS;

   // check that initial conditions were supplied, if they are required i.e. SS was not computed when the model was generated
   // MRE
#ifndef COMPUTE_Y0
   if (_y0 == NULL) {
	  PyErr_SetString(LNAError, "Initial conditions for MRE must be explicitly supplied");
	  return NULL;
   }
#endif

   // variance
#ifndef COMPUTE_V0
   if (_V0 == NULL) {
	  PyErr_SetString(LNAError, "Initial conditions for variance must be explicitly supplied");
	  return NULL;
   }
#endif

   // determine initial conditions based on arguments
   if ( (_y0 == NULL) && (_V0 == NULL)) {
	SS = LNA::SS_BOTH;
   } else if ((_y0 == NULL) && (_V0 != NULL)) {
	SS = LNA::SS_Y0;
	if (!(PyList_Check(_V0) || PyFloat_Check(_V0))) {
		PyErr_SetString(LNAError, "V0 has wrong type.");
		return NULL;
	}
	   assignDoubleFromPyList(_V0,V0);
   } else if ((_y0 != NULL) && (_V0 == NULL)) {
//	cout << "y0 is a list: " << PyList_Check(_y0) << endl;
	SS = LNA::SS_V0;
	if (!(PyList_Check(_y0) || PyFloat_Check(_y0))) {
		PyErr_SetString(LNAError, "y0 has wrong type.");
		return NULL;
	}
	   assignDoubleFromPyList(_y0,y0);
   } else {
	SS = LNA::SS_NONE;
	if (!(PyList_Check(_V0) || PyFloat_Check(_V0))) {
		PyErr_SetString(LNAError, "V0 has wrong type.");
		return NULL;
	}
	if (!( PyList_Check(_y0) || PyFloat_Check(_y0))) {
		PyErr_SetString(LNAError, "y0 has wrong type.");
		return NULL;
	}
	   assignDoubleFromPyList(_y0,y0);
	   assignDoubleFromPyList(_V0,V0);
   }

    // Theta
    if (!(PyList_Check(_Theta)) || PyFloat_Check(_Theta) || PyLong_Check(_Theta)) {
    	PyErr_SetString(LNAError, "Theta has wrong type.");
		return NULL;
    } else {
    	Theta = new double[npar];

    	if (PyList_Check(_Theta)) {
			for (int i=0; i<(npar-nvar); i++)
				if (PyFloat_Check(PyList_GetItem(_Theta,i)))
					Theta[i] = PyFloat_AsDouble(PyList_GetItem(_Theta,i));
				else Theta[i] = PyLong_AsDouble(PyList_GetItem(_Theta,i));
    	} else {
    		if (PyFloat_Check(_Theta))
    			Theta[0] = PyFloat_AsDouble(_Theta);
    		else
    			Theta[0] = PyLong_AsDouble(_Theta);
    	}
    }

    // time
    if (!PyList_Check(_time)) // no default value
    {
    	PyErr_SetString(LNAError, "Time must be a list");
		return NULL;
    } else {
        N = PyList_Size(_time);
        tspan.resize(N);
        for (int i=0; i<N; i++) {
        	tspan(i) = PyFloat_AsDouble(PyList_GetItem(_time,i));
        }
    }

    // observed variables
    if (_obsVar == NULL) {
    	// default all variables
    	varObs.resize(nvar);
    	for (int i=0; i<nvar; i++) {
    		varObs(i)=i;
    	}
    	nObsVar = nvar;
    }
    else {
		if (!((PyList_Check(_obsVar) || PyLong_Check(_obsVar) ))) {
			PyErr_SetString(LNAError, "obsVar has wrong type.");
			return NULL;
		}
		if (PyList_Check(_obsVar)) {
			nObsVar = PyList_Size(_obsVar);
			if (nObsVar > nvar)
			{
				PyErr_SetString(LNAError, "Too many observable variables specified.");
				return NULL;
			}
	    	varObs.resize(nObsVar);
	    	for (int i=0; i<nObsVar; i++) {
	    		PyObject *list_item = PyList_GetItem(_obsVar,i);
	    		if (!PyLong_Check(list_item)) {
	    			PyErr_SetString(LNAError, "obsVar must consist of integers only");
					return NULL;
	    		}
				varObs(i) = PyLong_AsLong(list_item)-1;
				if ( (varObs(i) < 0) || (varObs(i) > (nvar-1)))
				{
					PyErr_SetString(LNAError, "Invalid observable variable specified.");
					return NULL;
				}
	    	}
	    } else { // was a number
	    	nObsVar=1;
	    	varObs.resize(1);
	    	varObs(0) = PyLong_AsLong(_obsVar)-1;
	    }
	}


    // measurement error
	if (_merr == NULL) {
		// merr. default 0
		merr = new double[nObsVar];
    	for (int i=0; i<nObsVar; i++)
    		merr[i] = 0.0;

	} else {
		// parse measurement error
		if (!((PyList_Check(_merr) || PyFloat_Check(_merr) || PyLong_Check(_merr) ))) {
			PyErr_SetString(LNAError, "merr has wrong type.");
			return NULL;
		}
		int n;
		if (PyList_Check(_merr))
			n=PyList_Size(_merr);
		else
			n=1;

		if (n != nObsVar)
		{
			// was a list but of the wrong length
			PyErr_SetString(LNAError, "Number of measurement errors must be the same as the number of observed variables");
			return NULL;
		}
		else {
			assignDoubleFromPyList(_merr,merr);
		}

	}

    /* compute sensitivities ? */
    if (_computeSens != 0) {
		if (!PyBool_Check(_computeSens)) {
			PyErr_SetString(LNAError, "computeSens must be Boolean");
			return NULL;
		} else
			computeSens = PyObject_IsTrue(_computeSens);
    }

    if (_computeSens2 != 0) {
		if (!PyBool_Check(_computeSens2)) {
			PyErr_SetString(LNAError, "computeSens2 must be Boolean");
			return NULL;
		}
		else
			computeSens2 = PyObject_IsTrue(_computeSens2);
	}

	// always compute first order sens. if 2nd order is requested
	if (computeSens2 && !computeSens)
		computeSens = true;

	static LNA lna(nvar, npar, S);
//	cout << S << endl;
//	MA2 S_tmp(2,4);
//	S_tmp = 	1, -1, 0, 0,
//				0, 0, 1, -1;
//	cout << S_tmp << endl;
//	static LNA lna(2,6,S_tmp);

	/* temporary */
//	N = 101;
//	Vector tspan_tmp(N);
//
//	for (int i=0; i<N; i++)
//		tspan_tmp(i) =i*0.1;
	/* end */


	// output variables
	MA2 *Y 		= new MA2(nObsVar,N);
	MA4 *Sigma 	= new MA4(nObsVar, nObsVar, N, N);
	MA3 *traj_deriv	= new MA3(nObsVar, npar, N);
	MA5 *dSigma 	= new MA5(nObsVar, nObsVar, N, N, npar);

	// second order sens.
	MA4 *Sens2_MRE	= new MA4(nObsVar, npar, npar, N);
	MA6 *Sens2_Var 	= new MA6(nObsVar, nObsVar, npar, npar, N, N);

	// set to zero
    	*Y = 0; *Sigma=0; *traj_deriv=0; *dSigma=0; *Sens2_MRE=0; *Sens2_Var=0;

//	outputStruct os;
    	os.Y 			= Y;
    	os.traj_deriv 		= traj_deriv;
    	os.Sigma 		= Sigma;
	os.dSigma 		= dSigma;
	os.Sens2_MRE 		= Sens2_MRE;
	os.Sens2_Var  		= Sens2_Var;

//	/* temporary */
//	double Theta_tmp[] = { 20., 25., 10., 1., 2., 200.};
//	bool computeSens_tmp = true, computeSens2_tmp = true;
//	double y0_tmp[] = {20.,200.}; // initial conditions
//	double V0_tmp[] = {2.,5.,150.};
//	LNA::SS_FLAG SS_tmp = LNA::SS_BOTH;
//
//	double merr_tmp[] = {0.};
//	bool computeSens_tmp = false, computeSens2_tmp = false;
//	/* end */


	// debugging code
//	cout << "y0 ";
//	for (int i=0; i<nvar; i++) {
//		cout << setw(5) << y0[i];
//	}
//
//	cout << endl << "V0 ";
//	for (int i=0; i<nvar; i++) {
//		cout << setw(5) << V0[i];
//	}

//	cout << endl << "Theta ";
//	for (int i=0; i<npar; i++) {
//		cout << setw(5) << Theta[i];
//	}
//	cout << endl;
//
//	cout << "N " << N << endl;
//	cout << "cs1: " << computeSens << " cs2: " << computeSens2 << endl;
//	cout << "SS_flag " << SS << endl;
//	cout << "Time: ";
//	for (int i=0; i<N; i++)
//		cout << setw(5) << tspan(i) << endl;

//	cout << endl << "var obs " ;
//	for (int i=0; i<nObsVar; i++)
//		cout << setw(5) << varObs(i);
//	cout << endl;
//	for (int i=0; i<nObsVar; i++)
//			cout << setw(5) << merr[i];
//
//	cout << endl;

//	main(0,0);




	try {
//	    main(0,0);
//		Py_INCREF(Py_None);
//		return(Py_None);
//	 	Sigma

		lna.computeLinearNoise(y0, V0, Theta, computeSens, computeSens2,
				SS, os, tspan, varObs, merr);

//	     only tmp
//		lna.computeLinearNoise(y0_tmp, V0_tmp, Theta_tmp, computeSens_tmp, computeSens2_tmp,
//						SS_tmp, os, tspan_tmp, varObs_tmp, merr_tmp);

//		lna.computeLinearNoise(y0, V0, Theta_tmp, computeSens, computeSens2,
//								SS, os, tspan_tmp, varObs_tmp, merr_tmp);
//

	} catch(...) {
    	PyErr_SetString(LNAError, "Error encountered during computation.");
		return NULL;
	}



	/* assign output */
	// MRE
	long int dims_MRE[] = {nObsVar, N};
	PyObject *MRE_out = PyArray_SimpleNew(2, dims_MRE, NPY_DOUBLE);
	copyOut(MRE_out, Y);

//		Py_INCREF(Py_None);
//		return(Py_None);
	// Sigma

	long int dims_Sigma[] = {nObsVar, nObsVar, N, N};
	PyObject *Sigma_out = PyArray_SimpleNew(4, dims_Sigma, NPY_DOUBLE);
	copyOut(Sigma_out, Sigma);

	// MRE sensitivity
        PyObject *dMRE_out, *dSigma_out, *d2MRE_out, *d2Sigma_out, *ret;
        dMRE_out = dSigma_out = d2MRE_out = d2Sigma_out = ret = NULL;

	if (computeSens)
	{
		long int dims_dMRE[] = {nObsVar, npar, N};
		dMRE_out = PyArray_SimpleNew(3, dims_dMRE, NPY_DOUBLE);
		copyOut(dMRE_out, traj_deriv);

		long int dims_dSigma[] = {nObsVar, nObsVar, N, N, npar};
		dSigma_out = PyArray_SimpleNew(5, dims_dSigma, NPY_DOUBLE);
		copyOut(dSigma_out, dSigma);
	}

	if (computeSens2)
	{
		long int dims_d2MRE[] = {nObsVar, npar, npar, N};
		d2MRE_out = PyArray_SimpleNew(4, dims_d2MRE, NPY_DOUBLE);
		copyOut(d2MRE_out, Sens2_MRE);

		long int dims_d2Sigma[] = {nObsVar, nObsVar, npar, npar, N, N};
		d2Sigma_out = PyArray_SimpleNew(6, dims_d2Sigma, NPY_DOUBLE);
		copyOut(d2Sigma_out, Sens2_Var);
	}

	if (!computeSens) {
		ret = Py_BuildValue("OO", MRE_out, Sigma_out);
		Py_INCREF(MRE_out);
		Py_INCREF(Sigma_out);

	} else if (!computeSens2) {
		ret = Py_BuildValue("OOOO", MRE_out, Sigma_out, dMRE_out, dSigma_out);
		Py_INCREF(MRE_out);
		Py_INCREF(Sigma_out);
		Py_INCREF(dMRE_out);
		Py_INCREF(dSigma_out);
	} else {
		ret = Py_BuildValue("OOOOOO", MRE_out, Sigma_out, dMRE_out, dSigma_out, d2MRE_out, d2Sigma_out);
		Py_INCREF(MRE_out);
		Py_INCREF(Sigma_out);
		Py_INCREF(dMRE_out);
		Py_INCREF(dSigma_out);
		Py_INCREF(d2MRE_out);
		Py_INCREF(d2Sigma_out);
	}
	/* clean up */
	if (y0 != NULL)
		delete[] y0;
	if (V0 != NULL)
		delete[] V0;
	if (merr != NULL)
		delete[] merr;

	delete[] Theta; //TODO: restore
	Py_INCREF(ret);
	delete Y;
	delete Sigma;
	delete traj_deriv;
	delete dSigma;
	delete Sens2_MRE;
	delete Sens2_Var;
	
	return ret;
}

