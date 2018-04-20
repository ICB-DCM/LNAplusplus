
#include "LNAModule.cpp"

// Methods table
static PyMethodDef myModuleMethods[] = {
    {"LNA",  (PyCFunction)LNA_LNA,  METH_VARARGS | METH_KEYWORDS,
     "LNA(Theta, time, Y0 = None, V0 = None, merr = None, obsVar = None, computeSens = False, computeSens2 = False)"
     " -> MRE, Var [, dMRE, dSigma, [d2MRE, d2Sigma]]\n\n"
     "Compute the linear noise approximation.\n\n"
     "Parameters:\n"
     "-----------\n"
     "Theta (numeric list): model parameters\n"
     "time (numeric list): output times\n"
     "Y0 (numeric list): initial values for each species\n"
     "V0 (numeric list): initial values for each entry of the (upper triangular portion of the) covariance matrix in column-major ordering\n"
     "merr (numeric or numeric list): measurement error, single value for all observables or list with one value per observable\n"
     "obsVar (numeric list): list of observable state indices (1-based!!)\n"
     "computeSens (boolean): compute first-order sensitivities\n"
     "computeSens2 (boolean): compute second-order sensitivities\n"
     "\n"
     "Output:\n"
     "-------\n"
     "MRE: solution of the macroscopic rate equation for the observed variables specified\n"
     "Var: block cross-covariance matrix for the observed variables specified\n"
     "dMREdTheta: sensitivity tensor of the MRE for each observed variable with respect to each parameter, at each time point\n"
     "dVardTheta: sensitivity tensor of the block cross-covariance for each entry of the covariance matric w.r.t. all parameters\n"
     "d2MREdTheta2: second order sensitivity tensor of the MRE solution\n"
     "d2VardTheta2: second order sensitivity tensor of the block cross-covariance solution\n"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

// Module definition structure
static struct PyModuleDef myModule = {
   PyModuleDef_HEAD_INIT,
   "myModule",   /* name of module */
   "myModule  model created using LNA++", /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   myModuleMethods
};


// Module init function
PyMODINIT_FUNC
PyInit_myModuleLNA(void)
{
//    return PyModule_Create(&LNAmodule);
	PyObject *m = PyModule_Create(&myModule);
	if (m==NULL)
		return NULL;

	LNAError = PyErr_NewException("myModule.error", NULL, NULL);
    Py_INCREF(LNAError);
    PyModule_AddObject(m, "error", LNAError);

    import_array();
    return m;
}
