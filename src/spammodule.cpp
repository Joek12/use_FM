//
// Created by Joseph Kang on 2019-07-02.
//

#define PYSSIZE_T_CLEAN
#include "Python/Python.h"
#include "spammodule.h"

static PyObject *
spam_system(PyObject *self, PyObject *args){
    const char *command;
    int sts;

    if(!PyArg_ParseTuple(args, "s", &command))
        return NULL;

    sts = system(command);
    return PyLong_FromLong(sts);

}