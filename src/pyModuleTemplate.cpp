
#include "LNAModule.cpp"

// Methods table
static PyMethodDef myModuleMethods[] = {
    {"LNA",  (PyCFunction)LNA_LNA,  METH_VARARGS | METH_KEYWORDS,
     "Compute the linear noise approximation."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

// Module definition structure
static struct PyModuleDef myModule = {
   PyModuleDef_HEAD_INIT,
   "myModule",   /* name of module */
   NULL, /* module documentation, may be NULL */
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
