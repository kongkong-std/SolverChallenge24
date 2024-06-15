#ifndef MYSOLVER_H_
#define MYSOLVER_H_

// #define DIRECT_SOLVER
#define ITERATIVE_SOLVER

// please add your code in this file
struct MySolver
{
};

struct MySolver_complex
{
};

#ifdef DIRECT_SOLVER
//! real system
void preprocess(struct MySolver *solver, const int n, const int *row_ptr, const int *col_idx);

void direct_solver(struct MySolver *solver, const int n, const double *val);

void solve(struct MySolver *solver, const int n, const double *x, const double *b);

//! complex system
void preprocess_complex(struct MySolver_complex *solver, const int n, const int *row_ptr, const int *col_idx);

void direct_solver_complex(struct MySolver_complex *solver, const int n, const double *val, const double *val_im);

void solve_complex(struct MySolver_complex *solver, const int n, const double *x, const double *x_im, const double *b, const double *b_im);

#endif

#ifdef ITERATIVE_SOLVER
//! real system
void analyse(struct MySolver *solver, const int n, const int *row_ptr, const int *col_idx);

void preprocess(struct MySolver *solver, const int n, const double *val);

void iterative_solver(struct MySolver *solver, const int n, const double *x, const double *b);

//! complex system
void analyse_complex(struct MySolver_complex *solver, const int n, const int *row_ptr, const int *col_idx);

void preprocess_complex(struct MySolver_complex *solver, const int n, const double *val, const double *val_im);

void iterative_solver_complex(struct MySolver_complex *solver, const int n, const double *x, const double *x_im, const double *b, const double *b_im);

#endif

#endif