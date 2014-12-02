
#include "LNAModule.cpp"

// Methods table
static PyMethodDef BirthDeathMethods[] = {
    {"LNA",  (PyCFunction)LNA_LNA,  METH_VARARGS | METH_KEYWORDS,
     "Compute the linear noise approximation."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

// Module definition structure
static struct PyModuleDef BirthDeath = {
   PyModuleDef_HEAD_INIT,
   "BirthDeath",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   BirthDeathMethods
};


// Module init function
PyMODINIT_FUNC
PyInit_BirthDeathLNA(void)
{
//    return PyModule_Create(&LNAmodule);
	PyObject *m = PyModule_Create(&BirthDeath);
	if (m==NULL)
		return NULL;

	LNAError = PyErr_NewException("BirthDeath.error", NULL, NULL);
    Py_INCREF(LNAError);
    PyModule_AddObject(m, "error", LNAError);

    import_array();
    return m;
}
