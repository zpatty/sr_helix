/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: PCC_jacobian.h
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 17:40:15
 */

#ifndef PCC_JACOBIAN_H
#define PCC_JACOBIAN_H

/* Include Files */
#include "helix_controller_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void PCC_jacobian(const emxArray_real_T *q, double d, double L0,
                  const emxArray_real_T *qd, double X[36], double J[18],
                  double T_q[16], double dJ[18]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for PCC_jacobian.h
 *
 * [EOF]
 */
