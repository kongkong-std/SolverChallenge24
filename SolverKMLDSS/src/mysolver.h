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

#define KML_DSS_IR_

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

KmlSolverMatrixType ParseMatrixType(const char * /*matrix type*/);
KmlDssRefineMethod ParseRefineMethod(const char * /*refine method*/);

//! real system
void KMLRealSolverMatrixCreate(MySolver * /*solver type*/, KmlSolverMatrixType /*matrix type*/,
                               int /*size of linear system*/,
                               int * /*row ptr in csr*/, int * /*col idx in csr*/,
                               double * /*val in csr*/);
void KMLRealSolverRHSCreate(MySolver * /*solver type*/, int /*size of linear system*/,
                            double * /*rhs*/);
void KMLRealSolverSOLCreate(MySolver * /*solver type*/, int /*size of linear system*/,
                            double * /*sol*/);
void KMLRealSolverInitialize(MySolver * /*solver type*/, int /*number of threads*/);
void KMLRealSolverAnalyze(MySolver * /*solver type*/, int /*number of threads of fill-in reduction*/);
void KMLRealSolverFactor(MySolver * /*solver type*/, double /*factor threshold*/);
void KMLRealSolverSolve(MySolver * /*solver type*/, KmlDssRefineMethod /*refine method*/, int /*refine maxit*/, double /*refine tol*/);
void KMLRealSolverQuery(MySolver * /*solver type*/);
void KMLRealSolverClean(MySolver * /*solver type*/);

//! complex system
void KMLComplexSolverMatrixCreate(MySolverComplex * /*solver type*/, KmlSolverMatrixType /*matrix type*/,
                                  int /*size of linear system*/,
                                  int * /*row prt in csr*/, int * /*col idx in csr*/,
                                  kml_complex_double * /*double complex val in csr*/);
void KMLComplexSolverRHSCreate(MySolverComplex * /*solver type*/, int /*size of linear system*/,
                               kml_complex_double * /*double complex rhs*/);
void KMLComplexSolverSOLCreate(MySolverComplex * /*solver type*/, int /*size of linear system*/,
                               kml_complex_double * /*double complex sol*/);
void KMLComplexSolverInitialize(MySolverComplex * /*solver type*/, int /*number of threads*/);
void KMLComplexSolverAnalyze(MySolverComplex * /*solver type*/, int /*number of threads of fill-in reduction*/);
void KMLComplexSolverFactor(MySolverComplex * /*solver type*/, double /*factor threshold*/);
void KMLComplexSolverSolve(MySolverComplex * /*solver type*/, KmlDssRefineMethod /*refine method*/, int /*refine maxit*/, double /*refine tol*/);
void KMLComplexSolverQuery(MySolverComplex * /*solver type*/);
void KMLComplexSolverClean(MySolverComplex * /*solver type*/);
#endif