/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_helix_controller_api.h
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 24-Oct-2024 11:42:25
 */

#ifndef _CODER_HELIX_CONTROLLER_API_H
#define _CODER_HELIX_CONTROLLER_API_H

/* Include Files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <string.h>

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void helix_controller(real_T q[10], real_T dq[10], real_T qd[10],
                      real_T dqd[10], real_T ddqd[10], real_T d, real_T m,
                      real_T r, real_T kb, real_T ks, real_T bb, real_T bs,
                      real_T bm, real_T L0, real_T Kp[100], real_T KD[100],
                      real_T Kpx, real_T KDx, real_T xd[3], real_T dxd[3],
                      real_T dxr[3], real_T conv_pcc, real_T conv_motor,
                      real_T tau[10], real_T tau_r[10], real_T x[3],
                      real_T M[100], real_T C[10], real_T A[100], real_T cq[3]);

void helix_controller_api(const mxArray *const prhs[23], int32_T nlhs,
                          const mxArray *plhs[7]);

void helix_controller_atexit(void);

void helix_controller_initialize(void);

void helix_controller_terminate(void);

void helix_controller_xil_shutdown(void);

void helix_controller_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for _coder_helix_controller_api.h
 *
 * [EOF]
 */
