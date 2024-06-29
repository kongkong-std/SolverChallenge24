#ifndef MYSOLVER_H_
#define MYSOLVER_H_

// #define DIRECT_SOLVER
#define ITERATIVE_SOLVER

#include <petscksp.h>
#include "file_process.h"

// please add your code in this file
typedef struct my_solver
{
    KSP ksp;
    PC pc;
    Mat solver_a;
    Vec solver_b, solver_x, solver_r; // rhs, solution, residual
} MySolver;

#ifdef DIRECT_SOLVER
//! real system
void preprocess(MySolver *solver, const int n, const int *row_ptr, const int *col_idx);

void direct_solver(MySolver *solver, const int n, const double *val);

void solve(MySolver *solver, const int n, const double *x, const double *b);

//! complex system
void preprocess_complex(MySolverComplex *solver, const int n, const int *row_ptr, const int *col_idx);

void direct_solver_complex(MySolverComplex *solver, const int n, const double *val, const double *val_im);

void solve_complex(MySolverComplex *solver, const int n, const double *x, const double *x_im, const double *b, const double *b_im);

#endif // DIRECT_SOLVER

#ifdef ITERATIVE_SOLVER
/*
 * assigning value to linear system
 */
void SolverPetscInitialize(int argc, char **argv, const Matrix *, const Vector *, MySolver *mysolver);

/*
 * solver preprocess
 */
void SolverPetscPreprocess(int argc, char **argv, MySolver *mysolver);

/*
 * solver solve
 */
void SolverPetscSolve(int argc, char **argv, MySolver *mysolver);

/*
 * solver || residual || _ {lp} check
 */
void SolverPetscResidualCheck(int argc, char **argv, MySolver *mysolver);

#if 0
/*
 * solver get linear system
 */
void SolverPetscGetLinearSystem(const MySolver *mysolver, int *m, int *n, int *nnz,
                                int **row_ptr, int **col_idx, double **val, double **x, double **b);
#endif // 0 solver get linear systemc

#endif // ITERATIVE_SOLVER

#endif // MYSOLVER_H_