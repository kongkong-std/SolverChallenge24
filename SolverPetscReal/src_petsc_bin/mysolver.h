#ifndef MYSOLVER_H_
#define MYSOLVER_H_

// #define DIRECT_SOLVER
#define ITERATIVE_SOLVER

#include <petscksp.h>

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
void SolverPetscInitialize(int argc, char **argv, MySolver *mysolver);

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

/*
 * solver solution file io
 */
void SolverPetscSolutionFileIO(MySolver *);

/*
 * solver destroy
 */
void SolverPetscDestroy(MySolver *);

#endif // ITERATIVE_SOLVER

#endif // MYSOLVER_H_