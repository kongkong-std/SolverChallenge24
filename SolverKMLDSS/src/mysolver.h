#ifndef MYSOLVER_H_
#define MYSOLVER_H_

// system file
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <complex.h>
#include "kml_solver.h"

// please add your code in this file
typedef struct my_solver
{
    // add struct array
    KmlSolverMatrixStore dss_store_a; // matrix original data
    KmlSolverMatrixStore dss_store_b; // rhs original data
    KmlSolverMatrixStore dss_store_x; // solution original data

    KmlSolverMatrixOption dss_option_a; // matrix option
    KmlSolverMatrixOption dss_option_b; // rhs option
    KmlSolverMatrixOption dss_option_x; // solution option

    KmlDssSolver *dss_solver;
    KmlSolverMatrix *dss_solver_a; // solver matrix
    KmlSolverMatrix *dss_solver_b; // solver rhs
    KmlSolverMatrix *dss_solver_x; // solver solution

    KmlDssInitOption dss_solver_init_option;        // solver initial option
    KmlDssAnalyzeOption dss_solver_analyze_option;  // solver analyze option
    KmlDssFactorizeOption dss_solver_factor_option; // solver factor option
    KmlDssSolveOption dss_solver_solve_option;      // solver solve option
    KmlDssInfo dss_solver_info;                     // solver information
} MySolver;

typedef struct my_solver_complex
{
    // add struct array
    KmlSolverMatrixStore dss_store_a; // matrix original data
    KmlSolverMatrixStore dss_store_b; // rhs original data
    KmlSolverMatrixStore dss_store_x; // solution original data

    KmlSolverMatrixOption dss_option_a; // matrix option
    KmlSolverMatrixOption dss_option_b; // rhs option
    KmlSolverMatrixOption dss_option_x; // solution option

    KmlDssSolver *dss_solver;
    KmlSolverMatrix *dss_solver_a; // solver matrix
    KmlSolverMatrix *dss_solver_b; // solver rhs
    KmlSolverMatrix *dss_solver_x; // solver solution

    KmlDssInitOption dss_solver_init_option;        // solver initial option
    KmlDssAnalyzeOption dss_solver_analyze_option;  // solver analyze option
    KmlDssFactorizeOption dss_solver_factor_option; // solver factor option
    KmlDssSolveOption dss_solver_solve_option;      // solver solve option
    KmlDssInfo dss_solver_info;                     // solver information
} MySolverComplex;

//! real system
void KMLRealSolverInitialize(MySolver * /*solver type*/, double * /*sol*/, int /*size of linear system*/,
                             const int * /*row ptr in csr*/, const int * /*col idx in csr*/,
                             const double * /*val in csr*/, const double * /*rhs*/);
void KMLRealSolverAnalyze(MySolver * /*solver type*/);
void KMLRealSolverFactor(MySolver * /*solver type*/);
void KMLRealSolverSolve(MySolver * /*solver type*/);
void KMLRealSolverQuery(MySolver * /*solver type*/);
void KMLRealSolverClean(MySolver * /*solver type*/);

//! complex system
void KMLComplexSolverInitialize(MySolverComplex * /*solver type*/, kml_complex_double * /*double complex sol*/,
                                int /*size of linear system*/, const int * /*row prt in csr*/, const int * /*col idx in csr*/,
                                const kml_complex_double * /*double complex val in csr*/,
                                const kml_complex_double * /*double complex rhs*/);
void KMLComplexSolverAnalyze(MySolverComplex * /*solver type*/);
void KMLComplexSolverFactor(MySolverComplex * /*solver type*/);
void KMLComplexSolverSolve(MySolverComplex * /*solver type*/);
void KMLComplexSolverQuery(MySolverComplex * /*solver type*/);
void KMLComplexSolverClean(MySolverComplex * /*solver type*/);
#endif