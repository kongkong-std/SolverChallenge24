#ifndef SOLVER_PCSHELL_H_
#define SOLVER_PCSHELL_H_

#include <petscksp.h>

extern PetscErrorCode DefShellPCApply(PC, Vec, Vec);
extern PetscErrorCode DefShellPCSetup(PC);

#endif