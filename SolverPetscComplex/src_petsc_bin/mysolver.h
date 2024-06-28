#ifndef MYSOLVER_H_
#define MYSOLVER_H_

// #define DIRECT_SOLVER
#define ITERATIVE_SOLVER

#define CHALLENGE_06 // preconditioner for 06

#include <petscksp.h>

// please add your code in this file
typedef struct my_solver
{
    KSP ksp;
    PC pc;
    Mat solver_a;
    Vec solver_b, solver_x, solver_r; // rhs, solution, residual
#ifdef CHALLENGE_06
    Mat solver_a_re, solver_a_im, solver_a_im_oppo;
    Mat solver_pc_re, solver_pc_im, solver_pc_im_oppo;
    Vec solver_b_re, solver_b_im;
    Vec solver_x_re, solver_x_im;
    Vec solver_r_re, solver_r_im;
    Mat solver_block_a, solver_block_pc;
    Vec solver_block_b, solver_block_x, solver_block_r;
#endif
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

#ifdef CHALLENGE_06
/*
 * real part matrix
 */
void SolverRealPartMatrix(const Mat * /*original complex petsc matrix*/, Mat * /*real part matrix*/);

/*
 * imaginary part matrix
 */
void SolverImaginaryPartMatrix(const Mat * /*original complex petsc matrix*/, Mat * /*imaginary part matrix*/);

/*
 * real part vector
 */
void SolverRealPartVector(const Vec * /*original complex petsc vector*/, Vec * /*real part vector*/);

/*
 * imaginary part vector
 */
void SolverImaginaryPartVector(const Vec * /*original complex petsc vector*/, Vec * /*imaginary part vector*/);
#endif // CHALLENGE_06

#if 0
/*
 * solver get linear system
 */
void SolverPetscGetLinearSystem(const MySolver *mysolver, int *m, int *n, int *nnz,
                                int **row_ptr, int **col_idx, double **val, double **x, double **b);
#endif // 0 solver get linear system

#endif // ITERATIVE_SOLVER

#endif // MYSOLVER_H_
