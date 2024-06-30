/*******************************************************************************
 * Copyright (c) Huawei Technologies Co., Ltd. 2020-2021. All rights reserved.
 * Description: Part of KML library
 * Author: KML
 * Create: 2020
 ******************************************************************************/

#ifndef KML_DSS_DEFS_H_INCLUDED
#define KML_DSS_DEFS_H_INCLUDED

#include <stdint.h>
#include "kml_solver_defs.h"

#if defined(__cplusplus)
extern "C" {
#endif

typedef enum {
    KMLDSS_BWR_OFF,
    KMLDSS_BWR_FIXED_THREADS,
    KMLDSS_BWR_LAST,
} KmlDssBWRMode;

typedef enum {
    KMLDSS_MATCHING_OFF,
    KMLDSS_MATCHING_HUNGARIAN,
    KMLDSS_MATCHING_MC64,
    KMLDSS_MATCHING_SMC64,
    KMLDSS_MATCHING_LAST,
} KmlDssMatchingType;

typedef enum {
    KMLDSS_RDR_KRDR,
    KMLDSS_RDR_CUSTOM,
    KMLDSS_RDR_LAST,
} KmlDssRdrType;

typedef int (*KmlDssCustomRdrAlgo) (const KmlSolverMatrixStore *store, void *perm, void *iperm, void* arg);

typedef enum {
    KMLDSS_PIVOTING_OFF,
    KMLDSS_PIVOTING_LOCAL,
    KMLDSS_PIVOTING_LAST,
} KmlDssPivotingType;

typedef enum {
    KMLDSS_MATRIX_NO_TRANS,
    KMLDSS_MATRIX_TRANS,
    KMLDSS_MATRIX_LAST,
} KmlDssMatrixTransType;

typedef enum {
    KMLDSS_SOLVE_ALL,
    KMLDSS_SOLVE_FORWARD,
    KMLDSS_SOLVE_DIAGONAL,
    KMLDSS_SOLVE_BACKWARD,
    KMLDSS_SOLVE_REFINE,
    KMLDSS_SOLVE_LAST,
} KmlDssSolveStage;

typedef enum {
    KMLDSS_REFINE_OFF,
    KMLDSS_REFINE_CLASSICAL,
    KMLDSS_REFINE_LAST,
} KmlDssRefineMethod;

typedef enum {
    KMLDSS_INIT_OPTION_BWR_MODE = KMLSS_BIT(0),
    KMLDSS_INIT_OPTION_NTHREADS = KMLSS_BIT(1),
} KmlDssInitOptionField;

typedef struct {
    uint64_t fieldMask;
    KmlDssBWRMode bwrMode;
    int32_t nThreads;
} KmlDssInitOption;

typedef enum {
    KMLDSS_ANALYZE_OPTION_NTHREADS_ANALYZE = KMLSS_BIT(0),
    KMLDSS_ANALYZE_OPTION_NTHREADS_FACTORIZE = KMLSS_BIT(1),
    KMLDSS_ANALYZE_OPTION_NTHREADS_SOLVE = KMLSS_BIT(2),
    KMLDSS_ANALYZE_OPTION_MATCHING_TYPE = KMLSS_BIT(3),
    KMLDSS_ANALYZE_OPTION_NTHREADS_MATCHING = KMLSS_BIT(4),
    KMLDSS_ANALYZE_OPTION_RDR_TYPE = KMLSS_BIT(5),
    KMLDSS_ANALYZE_OPTION_NTHREADS_RDR = KMLSS_BIT(6),
    KMLDSS_ANALYZE_OPTION_CUSTOM_RDR_ALGO = KMLSS_BIT(7),
    KMLDSS_ANALYZE_OPTION_CUSTOM_RDR_ARGS = KMLSS_BIT(8),
    KMLDSS_ANALYZE_OPTION_RDR_PERM = KMLSS_BIT(9),
    KMLDSS_ANALYZE_OPTION_PIVOTING_TYPE = KMLSS_BIT(10),
    KMLDSS_ANALYZE_OPTION_PIVOTING_THRESHOLD = KMLSS_BIT(11),
    KMLDSS_ANALYZE_OPTION_BATCH_NRHS = KMLSS_BIT(12),
    KMLDSS_ANALYZE_OPTION_SCHUR_SIZE = KMLSS_BIT(13),
    KMLDSS_ANALYZE_OPTION_SCHUR_FORMAT = KMLSS_BIT(14),
    KMLDSS_ANALYZE_OPTION_MATIRX_TRANS_TYPE = KMLSS_BIT(15),
} KmlDssAnalyzeOptionField;

typedef struct {
    uint64_t fieldMask;
    int32_t nThreadsAnalyze;
    int32_t nThreadsFactorize;
    int32_t nThreadsSolve;
    KmlDssMatchingType matchingType;
    int32_t nThreadsMatching;
    KmlDssRdrType rdrType;
    int32_t nThreadsRdr;
    KmlDssCustomRdrAlgo customRdrAlgo;
    void *customRdrArgs;
    void *rdrPerm;
    KmlDssPivotingType pivotingType;
    double pivotingThreshold;
    int64_t batchNRHS;
    int64_t schurSize;
    KmlSolverMatrixStoreFormat schurFormat;
    KmlDssMatrixTransType matrixTransType;
} KmlDssAnalyzeOption;

typedef enum {
    KMLDSS_FACTORIZE_OPTION_PERTURBATION_THRESHOLD = KMLSS_BIT(0),
} KmlDssFactorizeOptionField;

typedef struct {
    uint64_t fieldMask;
    double perturbationThreshold;
} KmlDssFactorizeOption;

typedef enum {
    KMLDSS_SOLVE_OPTION_SOLVE_STAGE = KMLSS_BIT(0),
    KMLDSS_SOLVE_OPTION_REFINE_METHOD = KMLSS_BIT(1),
    KMLDSS_SOLVE_OPTION_REFINE_MAX_STEP = KMLSS_BIT(2),
    KMLDSS_SOLVE_OPTION_REFINE_TOLERANCE = KMLSS_BIT(3),
} KmlDssSolveOptionField;

typedef struct {
    uint64_t fieldMask;
    KmlDssSolveStage stage;
    KmlDssRefineMethod refineMethod;
    int32_t refineMaxSteps;
    double refineTolerance;
} KmlDssSolveOption;

typedef enum {
    KMLDSS_INFO_RDR_PERM = KMLSS_BIT(0),
    KMLDSS_INFO_SCHUR_NNZ = KMLSS_BIT(1),
    KMLDSS_INFO_SCHUR_MAT = KMLSS_BIT(2),
    KMLDSS_INFO_REFINE_STEPS = KMLSS_BIT(3),
    KMLDSS_INFO_PEAK_MEM = KMLSS_BIT(4),
    KMLDSS_INFO_FILL_IN = KMLSS_BIT(5),
    KMLDSS_INFO_FILL_IN_L = KMLSS_BIT(6),
    KMLDSS_INFO_FILL_IN_U = KMLSS_BIT(7),
} KmlDssInfoField;

typedef struct {
    uint64_t fieldMask;
    void *rdrPerm;
    int64_t schurNnz;
    KmlSolverMatrixStore *schurMat;
    int32_t refineSteps;
    int64_t peakMem;
    int64_t nFillIn;
    int64_t nFillInL;
    int64_t nFillInU;
} KmlDssInfo;

typedef struct KmlDssSolver KmlDssSolver;

#if defined(__cplusplus)
}
#endif

#endif // KML_DSS_DEFS_H_INCLUDED