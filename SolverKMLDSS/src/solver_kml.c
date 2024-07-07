#include "mysolver.h"

void KMLRealSolverMatrixCreate(MySolver *solver, int n,
                               int *row_ptr, int *col_idx, double *val)
{
    int ierr = 0;

    // mat data process
    printf("\n>>>> kml-dss matrix begin assembling...\n");
    (solver->dss_store_a).indexType = KMLSS_INDEX_INT32;
    (solver->dss_store_a).valueType = KMLSS_VALUE_FP64;
    (solver->dss_store_a).nRow = n;
    (solver->dss_store_a).nCol = n;
    (solver->dss_store_a).format = KMLSS_MATRIX_STORE_CSR;
    (solver->dss_store_a).csr.rowOffset = row_ptr;
    (solver->dss_store_a).csr.colIndex = col_idx;
    (solver->dss_store_a).csr.value = val;
    (solver->dss_option_a).fieldMask = KMLSS_MATRIX_OPTION_TYPE;
    (solver->dss_option_a).type = KMLSS_MATRIX_GEN;
    ierr = KmlSolverMatrixCreate(&(solver->dss_solver_a), &(solver->dss_store_a), &(solver->dss_option_a));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR when create matrix A: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss matrix has been assembled!!!\n\n");
}

void KMLRealSolverRHSCreate(MySolver *solver, int n, double *rhs)
{
    int ierr = 0;

    // rhs data process
    printf(">>>> kml-dss rhs vector begin assembling...\n");
    (solver->dss_store_b).indexType = KMLSS_INDEX_INT32;
    (solver->dss_store_b).valueType = KMLSS_VALUE_FP64;
    (solver->dss_store_b).nRow = n;
    (solver->dss_store_b).nCol = 1;
    (solver->dss_store_b).format = KMLSS_MATRIX_STORE_DENSE_COL_MAJOR;
    (solver->dss_store_b).dense.value = rhs;
    (solver->dss_store_b).dense.ld = n;
    (solver->dss_option_b).fieldMask = KMLSS_MATRIX_OPTION_TYPE;
    (solver->dss_option_b).type = KMLSS_MATRIX_GEN;
    ierr = KmlSolverMatrixCreate(&(solver->dss_solver_b), &(solver->dss_store_b), &(solver->dss_option_b));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR when create rhs vector b: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss rhs vector has been assembled!!!\n\n");
}

void KMLRealSolverSOLCreate(MySolver *solver, int n, double *sol)
{
    int ierr = 0;

    // sol data process
    printf(">>>> kml-dss sol vector begin assembling...\n ");
    (solver->dss_store_x).indexType = KMLSS_INDEX_INT32;
    (solver->dss_store_x).valueType = KMLSS_VALUE_FP64;
    (solver->dss_store_x).nRow = n;
    (solver->dss_store_x).nCol = 1;
    (solver->dss_store_x).format = KMLSS_MATRIX_STORE_DENSE_COL_MAJOR;
    (solver->dss_store_x).dense.value = sol;
    (solver->dss_store_x).dense.ld = n;
    (solver->dss_option_x).fieldMask = KMLSS_MATRIX_OPTION_TYPE;
    (solver->dss_option_x).type = KMLSS_MATRIX_GEN;
    ierr = KmlSolverMatrixCreate(&(solver->dss_solver_x), &(solver->dss_store_x), &(solver->dss_option_x));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR when create solution x: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss sol vector has been assembled!!!\n\n");
}

void KMLRealSolverInitialize(MySolver *solver, int n_thread)
{
    int ierr = 0;

    // Init solver
    printf(">>>> kml-dss solver begin to initialize...\n");
    (solver->dss_solver_init_option).fieldMask = KMLDSS_INIT_OPTION_BWR_MODE | KMLDSS_INIT_OPTION_NTHREADS;
    (solver->dss_solver_init_option).bwrMode = KMLDSS_BWR_OFF;
    (solver->dss_solver_init_option).nThreads = n_thread;
    ierr = KmlDssInit(&(solver->dss_solver), &(solver->dss_solver_init_option));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssInit: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss solver has been initialized!!!\n\n");
}

void KMLRealSolverAnalyze(MySolver *solver, int n_thread_rdr)
{
    int ierr = 0;

    // Analyze
    printf(">>>> kml-dss solver begin to analyze...\n");
    (solver->dss_solver_analyze_option).fieldMask = KMLDSS_ANALYZE_OPTION_MATCHING_TYPE |
                                                    KMLDSS_ANALYZE_OPTION_RDR_TYPE |
                                                    KMLDSS_ANALYZE_OPTION_NTHREADS_RDR;
    (solver->dss_solver_analyze_option).matchingType = KMLDSS_MATCHING_OFF;
    (solver->dss_solver_analyze_option).rdrType = KMLDSS_RDR_KRDR;
    (solver->dss_solver_analyze_option).nThreadsRdr = n_thread_rdr;
    ierr = KmlDssAnalyze(solver->dss_solver, solver->dss_solver_a, &(solver->dss_solver_analyze_option));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssAnalyze: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss solver has been analyzed!!!\n\n");
}

void KMLRealSolverFactor(MySolver *solver, double factor_threshold)
{
    int ierr = 0;

    // Factorize
    printf(">>>> kml-dss solver begin to factorize...\n");
    (solver->dss_solver_factor_option).fieldMask = KMLDSS_FACTORIZE_OPTION_PERTURBATION_THRESHOLD;
    (solver->dss_solver_factor_option).perturbationThreshold = factor_threshold;
    ierr = KmlDssFactorize(solver->dss_solver, solver->dss_solver_a, &(solver->dss_solver_factor_option));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssFactorize: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss solver has been factorized!!!\n\n");
}

void KMLRealSolverSolve(MySolver *solver)
{
    int ierr = 0;

    // Solve
    printf(">>>> kml-dss solver begin to solver...\n");
    (solver->dss_solver_solve_option).fieldMask = KMLDSS_SOLVE_OPTION_SOLVE_STAGE | KMLDSS_SOLVE_OPTION_REFINE_METHOD;
    (solver->dss_solver_solve_option).stage = KMLDSS_SOLVE_ALL;
    (solver->dss_solver_solve_option).refineMethod = KMLDSS_REFINE_OFF;
    ierr = KmlDssSolve(solver->dss_solver, solver->dss_solver_b, solver->dss_solver_x, &(solver->dss_solver_solve_option));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssSolve: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss solver has been solved!!!\n\n");
}

void KMLRealSolverQuery(MySolver *solver)
{
#if 0
    // Output result x
    printf("Result of first factorize and solve:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%lf ", x[i]);
    }
    printf("\n");
#endif
    int ierr = 0;

    // Query
    printf(">>>> kml-dss solver begin to show info...\n");
    (solver->dss_solver_info).fieldMask = KMLDSS_INFO_PEAK_MEM;
    ierr = KmlDssQuery(solver->dss_solver, &(solver->dss_solver_info));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssQuery: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf("Peak memory is %ld Byte\n", (solver->dss_solver_info).peakMem);
    printf(">>>> kml-dss solver has shown info!!!\n\n");
}

void KMLRealSolverClean(MySolver *solver)
{
    int ierr = 0;

    // Destroy
    printf(">>>> kml-dss solver begin to destroy...\n");
    ierr = KmlDssClean(&(solver->dss_solver));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssClean: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss solver has been destroyed!!!\n\n");

    printf(">>>> kml-dss matrix begin to destroy...\n");
    ierr = KmlSolverMatrixDestroy(&(solver->dss_solver_a));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlSolverMatrixDestroy matrix: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss matrix has beed destroyed!!!\n\n");

    printf(">>>> kml-dss rhs begin to destroy...\n");
    ierr = KmlSolverMatrixDestroy(&(solver->dss_solver_b));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlSolverMatrixDestroy rhs: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss rhs has beed destroyed!!!\n\n");

    printf(">>>> kml-dss sol begin to destroy...\n");
    ierr = KmlSolverMatrixDestroy(&(solver->dss_solver_x));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlSolverMatrixDestroy solution: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss sol has beed destroyed!!!\n\n");
}

void KMLComplexSolverMatrixCreate(MySolverComplex *solver, int n,
                                  int *row_ptr, int *col_idx, kml_complex_double *val)
{
    int ierr = 0;

    // mat data process
    printf("\n>>>> kml-dss matrix begin assembling...\n");
    (solver->dss_store_a).indexType = KMLSS_INDEX_INT32;
    (solver->dss_store_a).valueType = KMLSS_VALUE_FP64C;
    (solver->dss_store_a).nRow = n;
    (solver->dss_store_a).nCol = n;
    (solver->dss_store_a).format = KMLSS_MATRIX_STORE_CSR;
    (solver->dss_store_a).csr.rowOffset = row_ptr;
    (solver->dss_store_a).csr.colIndex = col_idx;
    (solver->dss_store_a).csr.value = val;
    (solver->dss_option_a).fieldMask = KMLSS_MATRIX_OPTION_TYPE;
    (solver->dss_option_a).type = KMLSS_MATRIX_GEN;
    ierr = KmlSolverMatrixCreate(&(solver->dss_solver_a), &(solver->dss_store_a), &(solver->dss_option_a));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR when create matrix A: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss matrix has been assembled!!!\n\n");
}

void KMLComplexSolverRHSCreate(MySolverComplex *solver, int n,
                               kml_complex_double *rhs)
{
    int ierr = 0;

    // rhs data process
    printf(">>>> kml-dss rhs vector begin assembling...\n");
    (solver->dss_store_b).indexType = KMLSS_INDEX_INT32;
    (solver->dss_store_b).valueType = KMLSS_VALUE_FP64C;
    (solver->dss_store_b).nRow = n;
    (solver->dss_store_b).nCol = 1;
    (solver->dss_store_b).format = KMLSS_MATRIX_STORE_DENSE_COL_MAJOR;
    (solver->dss_store_b).dense.value = rhs;
    (solver->dss_store_b).dense.ld = n;
    (solver->dss_option_b).fieldMask = KMLSS_MATRIX_OPTION_TYPE;
    (solver->dss_option_b).type = KMLSS_MATRIX_GEN;
    ierr = KmlSolverMatrixCreate(&(solver->dss_solver_b), &(solver->dss_store_b), &(solver->dss_option_b));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR when create rhs vector b: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss rhs vector has been assembled!!!\n\n");
}

void KMLComplexSolverSOLCreate(MySolverComplex *solver, int n,
                               kml_complex_double *sol)
{
    int ierr = 0;

    // sol data process
    printf(">>>> kml-dss sol vector begin assembling...\n ");
    (solver->dss_store_x).indexType = KMLSS_INDEX_INT32;
    (solver->dss_store_x).valueType = KMLSS_VALUE_FP64C;
    (solver->dss_store_x).nRow = n;
    (solver->dss_store_x).nCol = 1;
    (solver->dss_store_x).format = KMLSS_MATRIX_STORE_DENSE_COL_MAJOR;
    (solver->dss_store_x).dense.value = sol;
    (solver->dss_store_x).dense.ld = n;
    (solver->dss_option_x).fieldMask = KMLSS_MATRIX_OPTION_TYPE;
    (solver->dss_option_x).type = KMLSS_MATRIX_GEN;
    ierr = KmlSolverMatrixCreate(&(solver->dss_solver_x), &(solver->dss_store_x), &(solver->dss_option_x));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR when create solution x: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss sol vector has been assembled!!!\n\n");
}

void KMLComplexSolverInitialize(MySolverComplex *solver, int n_thread)
{
    int ierr = 0;

    // Init solver
    printf(">>>> kml-dss solver begin to initialize...\n");
    (solver->dss_solver_init_option).fieldMask = KMLDSS_INIT_OPTION_BWR_MODE | KMLDSS_INIT_OPTION_NTHREADS;
    (solver->dss_solver_init_option).bwrMode = KMLDSS_BWR_OFF;
    (solver->dss_solver_init_option).nThreads = n_thread;
    ierr = KmlDssInit(&(solver->dss_solver), &(solver->dss_solver_init_option));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssInit: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss solver has been initialized!!!\n\n");
}

void KMLComplexSolverAnalyze(MySolverComplex *solver, int n_thread_rdr)
{
    int ierr = 0;

    // Analyze
    printf(">>>> kml-dss solver begin to analyze...\n");
    (solver->dss_solver_analyze_option).fieldMask = KMLDSS_ANALYZE_OPTION_MATCHING_TYPE |
                                                    KMLDSS_ANALYZE_OPTION_RDR_TYPE |
                                                    KMLDSS_ANALYZE_OPTION_NTHREADS_RDR;
    (solver->dss_solver_analyze_option).matchingType = KMLDSS_MATCHING_OFF;
    (solver->dss_solver_analyze_option).rdrType = KMLDSS_RDR_KRDR;
    (solver->dss_solver_analyze_option).nThreadsRdr = n_thread_rdr;
    ierr = KmlDssAnalyze(solver->dss_solver, solver->dss_solver_a, &(solver->dss_solver_analyze_option));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssAnalyze: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss solver has been analyzed!!!\n\n");
}

void KMLComplexSolverFactor(MySolverComplex *solver, double factor_threshold)
{
    int ierr = 0;

    // Factorize
    printf(">>>> kml-dss solver begin to factorize...\n");
    (solver->dss_solver_factor_option).fieldMask = KMLDSS_FACTORIZE_OPTION_PERTURBATION_THRESHOLD;
    (solver->dss_solver_factor_option).perturbationThreshold = factor_threshold;
    ierr = KmlDssFactorize(solver->dss_solver, solver->dss_solver_a, &(solver->dss_solver_factor_option));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssFactorize: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss solver has been factorized!!!\n\n");
}

void KMLComplexSolverSolve(MySolverComplex *solver)
{
    int ierr = 0;

    // Solve
    printf(">>>> kml-dss solver begin to solver...\n");
    (solver->dss_solver_solve_option).fieldMask = KMLDSS_SOLVE_OPTION_SOLVE_STAGE | KMLDSS_SOLVE_OPTION_REFINE_METHOD;
    (solver->dss_solver_solve_option).stage = KMLDSS_SOLVE_ALL;
    (solver->dss_solver_solve_option).refineMethod = KMLDSS_REFINE_OFF;
    ierr = KmlDssSolve(solver->dss_solver, solver->dss_solver_b, solver->dss_solver_x, &(solver->dss_solver_solve_option));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssSolve: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss solver has been solved!!!\n\n");
}

void KMLComplexSolverQuery(MySolverComplex *solver)
{
#if 0
    // Output result x
    printf("Result of first factorize and solve:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%lf ", x[i]);
    }
    printf("\n");
#endif
    int ierr = 0;

    // Query
    printf(">>>> kml-dss solver begin to show info...\n");
    (solver->dss_solver_info).fieldMask = KMLDSS_INFO_PEAK_MEM;
    ierr = KmlDssQuery(solver->dss_solver, &(solver->dss_solver_info));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssQuery: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf("Peak memory is %ld Byte\n", (solver->dss_solver_info).peakMem);
    printf(">>>> kml-dss solver has shown info!!!\n\n");
}

void KMLComplexSolverClean(MySolverComplex *solver)
{
    int ierr = 0;

    // Destroy
    printf(">>>> kml-dss solver begin to destroy...\n");
    ierr = KmlDssClean(&(solver->dss_solver));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlDssClean: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss solver has been destroyed!!!\n\n");

    printf(">>>> kml-dss matrix begin to destroy...\n");
    ierr = KmlSolverMatrixDestroy(&(solver->dss_solver_a));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlSolverMatrixDestroy matrix: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss matrix has beed destroyed!!!\n\n");

    printf(">>>> kml-dss rhs begin to destroy...\n");
    ierr = KmlSolverMatrixDestroy(&(solver->dss_solver_b));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlSolverMatrixDestroy rhs: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss rhs has beed destroyed!!!\n\n");

    printf(">>>> kml-dss sol begin to destroy...\n");
    ierr = KmlSolverMatrixDestroy(&(solver->dss_solver_x));
    if (ierr != KMLSS_NO_ERROR)
    {
        printf("ERROR in KmlSolverMatrixDestroy solution: %d\n", ierr);
        exit(EXIT_FAILURE);
    }
    printf(">>>> kml-dss sol has beed destroyed!!!\n\n");
}