/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_helix_controller_mex.h
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 29-Oct-2024 18:41:09
 */

#ifndef _CODER_HELIX_CONTROLLER_MEX_H
#define _CODER_HELIX_CONTROLLER_MEX_H

/* Include Files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
MEXFUNCTION_LINKAGE void mexFunction(int32_T nlhs, mxArray *plhs[],
                                     int32_T nrhs, const mxArray *prhs[]);

emlrtCTX mexFunctionCreateRootTLS(void);

void unsafe_helix_controller_mexFunction(int32_T nlhs, mxArray *plhs[7],
                                         int32_T nrhs, const mxArray *prhs[23]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for _coder_helix_controller_mex.h
 *
 * [EOF]
 */
