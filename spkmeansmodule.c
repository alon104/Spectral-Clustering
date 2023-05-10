#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "spkmeans.h"


/*function that buildes a python list from each C vector (array)*/
PyObject* make_list(double* vector, int vector_dimension){
    int i;
    PyObject* py_coordiante;
    PyObject* py_vector = PyList_New(vector_dimension);
    for(i=0; i<vector_dimension; i++){
        py_coordiante = Py_BuildValue("d", vector[i]);
        PyList_SetItem(py_vector, i, py_coordiante);
    }
    return py_vector;
}

/*function that buildes a python matrix with from C nxm matrix*/
PyObject* make_matrix(double** matrix, int n, int m){
    PyObject* py_matrix = PyList_New(n);
    PyObject* py_vector;
    int i;
    for(i=0; i<n; i++){
        py_vector = make_list(matrix[i], m);
        PyList_SetItem(py_matrix, i, py_vector);
    }
    return py_matrix;
}

/*function that frees a C matrix's memory allocation*/
void free_matrix(double** matrix, int n){
    int i;
    for(i=0; i<n; i++){
        free(matrix[i]);
    }
    free(matrix);
}

/*function that buildes a C nxm matrix from pyhton matrix*/
double** C_matrix(PyObject* py_matrix, int n, int m){
    int i, j;
    double **matrix;
    double *row;
    PyObject *py_row;
    matrix = (double **) malloc(sizeof(double *)*n);
    if (matrix == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for(i=0; i<n; i++){
        row = (double *) malloc(sizeof(double)*m);
        if (row == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
        py_row = PyList_GetItem(py_matrix, i);
        for(j=0; j<m; j++){
            *(row+j) = PyFloat_AsDouble(PyList_GetItem(py_row, j));
        }
        *(matrix+i) = row;
    }
    return matrix;
}

/*jacobi interface for 'goal == jacobi/spk'*/
PyObject* jacobi_interface(PyObject* py_matrix, int n){
    double **matrix;
    double *** AV;
    PyObject *py_A;
    PyObject *py_V;
    PyObject *py_AV;

    /*building C matrix from python matrix*/
    matrix = C_matrix(py_matrix, n, n);

    /*run jacobi*/
    AV = fit_jacobi(matrix, n);

    /*make python ret list*/
    py_A = make_matrix(AV[0], n, n);
    py_V = make_matrix(AV[1], n, n);
    py_AV = PyList_New(2);
    PyList_SetItem(py_AV, 0, py_A);
    PyList_SetItem(py_AV, 1, py_V);
    
    /*free matrices returned from jacobi algorithm*/
    free_matrix(AV[0], n);
    free_matrix(AV[1], n);
    free(AV);

    return py_AV;
}

/*wam wraper*/
static PyObject* wam(PyObject *self, PyObject *args){   
    int n;
    int d;
    double **matrix;
    double** ret;
    PyObject *py_matrix;

    if (!PyArg_ParseTuple(args, "Oii", &py_matrix, &n, &d)){
        return NULL;
    }

    /*wam algorithm and converting back to python*/
    matrix = C_matrix(py_matrix, n, d);
    ret = fit_wam(matrix, n, d);
    py_matrix = make_matrix(ret, n, n);
    free_matrix(matrix, n);
    free_matrix(ret, n);

    return Py_BuildValue("O",py_matrix);
}

/*ddg wraper*/
static PyObject* ddg(PyObject *self, PyObject *args){
    int n;
    int d;
    double **matrix;
    double **ret;
    PyObject *py_matrix;

    if (!PyArg_ParseTuple(args, "Oii", &py_matrix, &n, &d)){
        return NULL;
    }

    /*ddg algorithm and printing*/
    matrix = C_matrix(py_matrix, n, d);
    ret = fit_ddg(matrix, n, d);
    py_matrix = make_matrix(ret, n, n);
    free_matrix(matrix, n);
    free_matrix(ret, n);

    return Py_BuildValue("O",py_matrix);
}

/*gl wraper*/
static PyObject* gl(PyObject *self, PyObject *args){
    int n;
    int d;
    double **matrix;
    double **ret;
    PyObject *py_matrix;

    if (!PyArg_ParseTuple(args, "Oii", &py_matrix, &n, &d)){
        return NULL;
    }

    /*gl algorithm and printing*/
    matrix = C_matrix(py_matrix, n, d);
    ret = fit_gl(matrix, n, d);
    py_matrix = make_matrix(ret, n, n);
    free_matrix(matrix, n);
    free_matrix(ret, n);

    return Py_BuildValue("O",py_matrix);
}

/*jacobi wraper*/
static PyObject* jacobi(PyObject *self, PyObject *args)
{   
    int n;
    PyObject *py_matrix;
    PyObject *py_AV;

    if (!PyArg_ParseTuple(args, "Oi", &py_matrix, &n)){
        return NULL;
    }

    py_AV = jacobi_interface(py_matrix, n);

    /*return the final list of the two matrices to python*/
    return Py_BuildValue("O",py_AV);
}

/*spk wraper*/
static PyObject* spk(PyObject *self, PyObject *args){
    int K;
    int maxiter;
    double EPS;
    int vector_dimension;
    int n;
    double **vectors_list;
    double **final_centroids;
    double **centroids_list;
    PyObject *Py_centroids;
    PyObject *Py_vectors_list;
    PyObject* ret_python_list;

    
    if (!PyArg_ParseTuple(args, "iidiiOO",&K, &maxiter, &EPS, &vector_dimension, &n, &Py_centroids, &Py_vectors_list)){
        return NULL;
    }

/*building array of centroids*/
    centroids_list = C_matrix(Py_centroids, K, K);

/*building array of all the vectors*/
    vectors_list = C_matrix(Py_vectors_list, n, K);

    final_centroids = (double **)fit_spk(K, maxiter, EPS, vector_dimension, n, centroids_list, vectors_list);

    /*building python list from C array*/
    ret_python_list = make_matrix(final_centroids, K, K);
    
    /*free centroids returned from kmeans algorithm*/
    free_matrix(final_centroids, K);
    free_matrix(vectors_list, n);

    /*return the final list of centroids to python*/
    return Py_BuildValue("O",ret_python_list);
}

static PyMethodDef kmeansspMethods[] = {
    {"wam",
    (PyCFunction) wam,
    METH_VARARGS,
    "The wam algorithm"
    },
    {"ddg",
    (PyCFunction) ddg,
    METH_VARARGS,
    "The ddg algorithm"
    },
    {"gl",
    (PyCFunction) gl,
    METH_VARARGS,
    "The gl algorithm"
    },
    {"jacobi",
    (PyCFunction) jacobi,
    METH_VARARGS,
    "The Jacobi algorithm"
    },
    {"spk",
    (PyCFunction) spk,
    METH_VARARGS,
    "The spk algorithm"
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef kmeansspmodule = 
{
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    kmeansspMethods
};

PyMODINIT_FUNC PyInit_mykmeanssp (void)
{
    PyObject *m;
    m = PyModule_Create(&kmeansspmodule);
    if (!m){
        return NULL;
    }
    return m;
};

